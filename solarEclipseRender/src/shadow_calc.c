// shadow_calc.c

// -------------------------------------------------
// Copyright 2019 Dominic Ford.

// This file is part of EclipseRender.

// EclipseRender is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// EclipseRender is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with EclipseRender.  If not, see <http://www.gnu.org/licenses/>.
// -------------------------------------------------

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "constants.h"
#include "settings.h"
#include "shadow_calc.h"

#include "mathsTools/julianDate.h"
#include "mathsTools/sphericalAst.h"

/**
 * siderealTime - Turns a Julian date into a sidereal time (in hours, at Greenwich)
 * @param JD [in] - Julian date
 * @return - Sidereal time in radians
 */
double siderealTime(double JD) {
    double T = (JD - 51545.0) / 36525.0; // See pages 87-88 of Astronomical Algorithms, by Jean Meeus
    return fmod(M_PI / 180 * (
            280.46061837 +
            360.98564736629 * (JD - 51545.0) +
            0.000387933 * T * T -
            T * T * T / 38710000.0
    ), 2 * M_PI);
}

/**
 * pos3D - Convert a latitude and longitude into a unit vector in Cartesian coordinates, where z-axis is orientated
 * in the direction of the north pole, and x-axis points towards zero point of longitude.
 * @param out [out] - A three-component Cartesian vector.
 * @param lat [in] - Latitude, degrees
 * @param lng [in] - Longitude, degrees
 */
void pos3D(double *out, double lat, double lng) {
    lat *= M_PI / 180;
    lng *= M_PI / 180;
    out[0] = cos(lng) * cos(lat);
    out[1] = sin(lng) * cos(lat);
    out[2] = sin(lat);
}

/**
 * earthTopocentricPositionICRF - Return the 3D position of a point on the Earth's surface, in ICRF coordinates,
 * relative to the solar system barycentre (the origin and coordinate system used by DE405).
 * @param out [out] - A three-component Cartesian vector.
 * @param lat [in] - Latitude, degrees
 * @param lng [in] - Longitude, degrees
 * @param radius_in_earth_radii [in] - The radial position, in Earth radii, of the location to query. Set to 1 for
 * Earth's surface.
 * @param pos_earth [in] - The 3D position of the centre of the Earth at the epoch, as quoted by DE405
 * @param epoch [in] - The Julian Day number when the calculation is to be performed
 * @param sidereal_time [in] - Sidereal time in degrees
 */
void earthTopocentricPositionICRF(double *out, double lat, double lng, double radius_in_earth_radii,
                                  const double *pos_earth, double epoch, double sidereal_time) {

    // In degrees, the declination and RA of the point on the sky directly above the requested point, for the
    // ecliptic of the epoch.
    const double lat_at_epoch = lat;
    const double lng_at_epoch = lng + sidereal_time;

    // Distance from Earth's centre, AU
    const double radius = sqrt(
            gsl_pow_2(RADIUS_EARTH_EQUATOR * cos(lat_at_epoch * M_PI / 180)) +
            gsl_pow_2(RADIUS_EARTH_POLE * sin(lat_at_epoch * M_PI / 180))
    ) / AU * radius_in_earth_radii;

    // Transform into J2000.0
    double lat_j2000, lng_j2000;
    ra_dec_to_j2000(lng_at_epoch * 12. / 180, lat_at_epoch, unix_from_jd(epoch), &lng_j2000, &lat_j2000);
    lng_j2000 *= 180. / 12;

    // Output position relative to the solar system barycentre, in ICRF coordinates
    pos3D(out, lat_j2000, lng_j2000);
    out[0] = out[0] * radius + pos_earth[0];
    out[1] = out[1] * radius + pos_earth[1];
    out[2] = out[2] * radius + pos_earth[2];
}

/**
 * _shadowFraction - Return the fraction of one circle covered by a smaller transiting circle. Formulae below taken from
 * Mandel & Agol (2002).
 * @param p [in] - The size ratio of the two circles (angular diameter moon / angular diameter sun)
 * @param z [in] - The angular separation of the centers of the two circles, in units of the angular radius of Sun.
 * @return - The fraction of the Sun covered by the Moon
 */
double _shadowFraction(double p, double z) {
    const double p2 = p * p;
    const double z2 = z * z;

    if ((1 + p) < z) {
        return 0;
    }

    if ((fabs(1 - p) < z) && (z <= (1 + p))) {
        const double kappa_0 = acos((p2 + z2 - 1) / (2 * p * z));
        const double kappa_1 = acos((1 - p2 + z2) / (2 * z));
        return 1 / M_PI * (p2 * kappa_0 + kappa_1 - sqrt((4 * z2 - gsl_pow_2(1 + z2 - p2)) / 4));
    }

    if (z <= (1 - p)) {
        return p2;
    }

    return 1;
}

/**
 * testIfAnnularEclipse - Test whether there is an annular eclipse at a particularly location on Earth at a particular
 * time.
 * @param lat - Latitude, degrees
 * @param lng - Longitude, degrees
 * @param JD - The epoch, as a Julian day number
 * @param radius - The radial position, in Earth radii, of the location to query. Set to 1 for Earth's surface.
 * @param pos_sun - The position of the Sun, as returned by ephemerisCompute, in AU, relative to solar system barycentre
 * @param pos_moon - The position of the Moon, as returned by ephemerisCompute
 * @param pos_earth - The position of the centre of the Earth, as returned by ephemerisCompute
 * @param sidereal_time - Sidereal time, radians
 * @return - Boolean flag indicating whether an annular eclipse is visible
 */
int testIfAnnularEclipse(double lat, double lng, double JD,
                         double radius, const double *pos_sun, const double *pos_moon, const double *pos_earth,
                         double sidereal_time) {
    double x0[3];
    earthTopocentricPositionICRF(x0, lat, lng, radius, pos_earth, JD, sidereal_time * 180 / M_PI);

    // Cosine rule on Sun - Earth - Moon triangle. Lengths in AU.
    const double a = hypot(hypot(x0[0] - pos_sun[0], x0[1] - pos_sun[1]), x0[2] - pos_sun[2]);
    const double b = hypot(hypot(x0[0] - pos_moon[0], x0[1] - pos_moon[1]), x0[2] - pos_moon[2]);
    const double c = hypot(hypot(pos_sun[0] - pos_moon[0], pos_sun[1] - pos_moon[1]), pos_sun[2] - pos_moon[2]);
    const double sun_moon_ang_sep = acos((a * a + b * b - c * c) / (2 * a * b)); // radians

    const double sun_ang_size = atan2(RADIUS_SUN, a * AU);
    const double moon_ang_size = atan2(RADIUS_MOON, b * AU);

    if (moon_ang_size >= sun_ang_size) return 0;

    if (sun_moon_ang_sep > sun_ang_size - moon_ang_size) return 0;

    return 1;
}

/**
 * getShadowFraction - Calculate the fraction of the Sun's disk which is covered by the Moon at a particularly location
 * on Earth at a particular time.
 * @param lat - Latitude, degrees
 * @param lng - Longitude, degrees
 * @param JD - The epoch, as a Julian day number
 * @param radius - The radial position, in Earth radii, of the location to query. Set to 1 for Earth's surface.
 * @param pos_sun - The position of the Sun, as returned by ephemerisCompute, in AU, relative to solar system barycentre
 * @param pos_moon - The position of the Moon, as returned by ephemerisCompute
 * @param pos_earth - The position of the centre of the Earth, as returned by ephemerisCompute
 * @param sidereal_time - Sidereal time, radians
 * @return - Fraction, in the range 0 to 1
 */
double getShadowFraction(double lat, double lng, double JD,
                         double radius, const double *pos_sun, const double *pos_moon, const double *pos_earth,
                         double sidereal_time) {
    double x0[3];
    earthTopocentricPositionICRF(x0, lat, lng, radius, pos_earth, JD, sidereal_time * 180 / M_PI);

    // Cosine rule on Sun - Earth - Moon triangle. Lengths in AU.
    const double a = hypot(hypot(x0[0] - pos_sun[0], x0[1] - pos_sun[1]), x0[2] - pos_sun[2]);
    const double b = hypot(hypot(x0[0] - pos_moon[0], x0[1] - pos_moon[1]), x0[2] - pos_moon[2]);
    const double c = hypot(hypot(pos_sun[0] - pos_moon[0], pos_sun[1] - pos_moon[1]), pos_sun[2] - pos_moon[2]);
    const double sun_moon_ang_sep = acos((a * a + b * b - c * c) / (2 * a * b)); // radians

    const double sun_ang_size = atan2(RADIUS_SUN, a * AU);
    const double moon_ang_size = atan2(RADIUS_MOON, b * AU);

    return _shadowFraction(moon_ang_size / sun_ang_size, sun_moon_ang_sep / sun_ang_size);
}

/**
 * allocate_shadow_map - Allocate an array of doubles to hold the fraction of the Sun's disk covered by the Moon at
 * every point in a 2D array. The array also stores the latitude and longitude of each pixel in the 2D array on the
 * Earth's surface, so that it can be used with any 2D projection.
 * @param x_size [in] - The horizontal size of the array
 * @param y_size [in] - The vertical size of the array
 * @return Uninitialised storage structure
 */
shadow_map *allocate_shadow_map(int x_size, int y_size) {
    shadow_map *output = (shadow_map *) malloc(sizeof(shadow_map));

    output->map = (double *) calloc(x_size * y_size, sizeof(double));
    output->lat = (double *) calloc(x_size * y_size, sizeof(double));
    output->lng = (double *) calloc(x_size * y_size, sizeof(double));

    return output;
}

/**
 * shadow_map_free - Free an array of doubles to hold the fraction of the Sun's disk covered by the Moon at
 * every point in a 2D array.
 * @param item [in] - The array to deallocate.
 */
void shadow_map_free(shadow_map *item) {
    free(item->map);
    free(item->lat);
    free(item->lng);
    free(item);
};

/**
 * calculate_where_sun_overhead - Calculate the latitude and longitude on Earth where the Sun is overhead at a
 * particular time.
 * @param lat_sun [out] - Latitude where Sun is overhead, radians
 * @param lng_sun [out] - Longitude where Sun is overhead, radians
 * @param sidereal_time [out] - Sidereal time at queried moment, radians
 * @param pos_sun [in] - The position of the Sun, as returned by ephemerisCompute, in AU, relative to solar system barycentre
 * @param pos_earth [in] - The position of the centre of the Earth, as returned by ephemerisCompute
 * @param jd [in] - Time point to query, specified as a Julian day number.
 */
void calculate_where_sun_overhead(double *lat_sun, double *lng_sun, double *sidereal_time,
                                  const double *pos_sun, const double *pos_earth, double jd) {

    double earth_sun_vector[3], earth_sun_unit_vector[3];

    earth_sun_vector[0] = pos_sun[0] - pos_earth[0];
    earth_sun_vector[1] = pos_sun[1] - pos_earth[1];
    earth_sun_vector[2] = pos_sun[2] - pos_earth[2];

    // Distance from Earth to Sun
    const double sun_dist = sqrt(
            gsl_pow_2(earth_sun_vector[0]) +
            gsl_pow_2(earth_sun_vector[1]) +
            gsl_pow_2(earth_sun_vector[2])
    );

    earth_sun_unit_vector[0] = earth_sun_vector[0] / sun_dist;
    earth_sun_unit_vector[1] = earth_sun_vector[1] / sun_dist;
    earth_sun_unit_vector[2] = earth_sun_vector[2] / sun_dist;

    const double dec_sun_j2000 = asin(earth_sun_unit_vector[2]) * 180 / M_PI;
    const double ra_sun_j2000 = atan2(earth_sun_unit_vector[1], earth_sun_unit_vector[0]) * 12 / M_PI;

    double dec_sun_at_epoch, ra_sun_at_epoch;
    ra_dec_from_j2000(ra_sun_j2000, dec_sun_j2000, unix_from_jd(jd), &ra_sun_at_epoch, &dec_sun_at_epoch);

    *sidereal_time = siderealTime(jd);  // radians
    *lat_sun = dec_sun_at_epoch * M_PI / 180;  // radians
    *lng_sun = ra_sun_at_epoch * M_PI / 12 - *sidereal_time;  // radians
}

/**
 * calculate_eclipse_map_2d - Calculate a map of the fraction of the Sun's disk covered by the Moon across a 2D flat
 * map of the world, at a given time instant.
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param jd [in] - The Julian day number of the current point in the simulation
 * @param pos_sun [in] - The position of the Sun, as returned by ephemerisCompute, in AU, relative to solar system barycentre
 * @param pos_earth [in] - The position of the centre of the Earth, as returned by ephemerisCompute
 * @param pos_moon [in] - The position of the Moon, as returned by ephemerisCompute
 * @param span_output [out] - Structure containing the start and end times of the partial and total/annular eclipses
 * @param greatest_shadow [out] - A map of the greatest degree of shadow seen in each pixel, over all the time points
 * sampled.
 * @return - Map of the Moon's shadow at a single time instant
 */
shadow_map *calculate_eclipse_map_2d(const settings *config,
                                     double jd, const double *pos_sun, const double *pos_earth, const double *pos_moon,
                                     time_span *span_output, shadow_map *greatest_shadow) {
    int x, y;

    // Calculate the latitude and longitude on Earth where the Sun is overhead
    double sidereal_time, lat_sun, lng_sun;
    calculate_where_sun_overhead(&lat_sun, &lng_sun, &sidereal_time, pos_sun, pos_earth, jd);

    // Calculate map of eclipse at time point
    shadow_map *shadow_map = allocate_shadow_map(config->x_size_2d, config->y_size_2d);

    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            // Work out the latitude and longitude represented by this point on the map
            double lat = 90 - (y * 180. / config->y_size_2d);
            double lng = (x * 360. / config->x_size_2d) + (-180) + 360;
            while (lng > 180) lng -= 360;

            // Convert point into 3D Cartesian space, with distances in Earth radii
            double pPoint[3];
            pos3D(pPoint, lat, lng);

            // Look up the point on the Earth's surface where the Sun is directly overhead
            double sun0[3];
            pos3D(sun0, lat_sun * 180 / M_PI, lng_sun * 180 / M_PI);

            // Look up the distance between the two points. If this is greater than sqrt(2), then Sun below horizon.
            const double d = gsl_pow_2(pPoint[0] - sun0[0]) +
                             gsl_pow_2(pPoint[1] - sun0[1]) +
                             gsl_pow_2(pPoint[2] - sun0[2]);
            const int night_time = (d > 2);

            // Superimpose shadow map over Earth
            double shadow;
            if (night_time) {
                // Do not show shadow on night side of Earth
                shadow = -1;
            } else {
                shadow = getShadowFraction(lat, lng, jd, 1, pos_sun, pos_moon, pos_earth, sidereal_time);
            }

            // Update shadow map array to show the eclipse magnitude at this location in the image
            const int offset = x + y * config->x_size_2d;
            shadow_map->map[offset] = shadow;
            shadow_map->lat[offset] = lat;  // degrees
            shadow_map->lng[offset] = lng;  // degrees

            // Update the greatest eclipse map to show the greatest eclipse magnitude reached in each pixel
            if (shadow > greatest_shadow->map[offset]) {
                greatest_shadow->map[offset] = shadow;

                if (shadow > span_output->greatest_eclipse_magnitude) {
                    span_output->greatest_eclipse_magnitude = shadow;  // In range 0-1
                    span_output->greatest_eclipse_latitude = lat;  // degrees
                    span_output->greatest_eclipse_longitude = lng;  // degrees
                }
            }

            // Keep track of the earliest and latest times when any partial eclipse is visible
            if (shadow > 0) {
                if ((span_output->partial_start <= 0) || (span_output->partial_start > jd)) {
                    span_output->partial_start = jd;
                }
                if ((span_output->partial_end <= 0) || (span_output->partial_end < jd)) {
                    span_output->partial_end = jd;
                }
            }

            // Keep track of the earliest and latest times when any total (or annular) eclipse is visible
            if (shadow > 0.98) {
                if ((span_output->total_start <= 0) || (span_output->total_start > jd)) {
                    span_output->total_start = jd;
                }
                if ((span_output->total_end <= 0) || (span_output->total_end < jd)) {
                    span_output->total_end = jd;
                }
            }
        }

    // Return output
    return shadow_map;
}

/**
 * calculate_eclipse_map_3d - Calculate a map of the fraction of the Sun's disk covered by the Moon across a 3D image
 * of the world, projected as seen from the Sun.
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param jd [in] - The Julian day number of the current point in the simulation
 * @param pos_sun [in] - The position of the Sun, as returned by ephemerisCompute, in AU, relative to solar system barycentre
 * @param pos_earth [in] - The position of the centre of the Earth, as returned by ephemerisCompute
 * @param pos_moon [in] - The position of the Moon, as returned by ephemerisCompute
 * @return - Map of the Moon's shadow at a single time instant
 */
shadow_map *calculate_eclipse_map_3d(const settings *config, double jd, const double *pos_sun, const double *pos_earth,
                                     const double *pos_moon) {
    int x, y;

    // Calculate the latitude and longitude on Earth where the Sun is overhead
    double sidereal_time, lat_sun, lng_sun;
    calculate_where_sun_overhead(&lat_sun, &lng_sun, &sidereal_time, pos_sun, pos_earth, jd);

    // Draw spherical map of Earth
    shadow_map *shadow_map = allocate_shadow_map(config->x_size_3d, config->y_size_3d);

    for (y = 0; y < config->y_size_3d; y++)
        for (x = 0; x < config->x_size_3d; x++) {
            double zn, p_radius;

            // Project point in 3D space where we measure eclipse fraction, measured in Earth radii with z-axis
            // pointing out of the screen

            // x-axis starts out point horizontally across the screen
            double xn = -(x - config->x_size_3d / 2.) / config->earth_pixel_radius;

            // y-axis starts out pointing vertically up the screen
            double yn = -(y - config->y_size_3d / 2.) / config->earth_pixel_radius;
            double n = gsl_pow_2(xn) + gsl_pow_2(yn);

            // z-axis points out of the screen
            if (n >= 1) {
                // This pixel is not on the Earth's surface, so project point above Earth
                zn = 0;
                p_radius = sqrt(n);

                // Renormalise (xn,yn,zn) is have unit length
                xn /= p_radius;
                yn /= p_radius;
            } else {
                // This pixel is on the Earth's surface
                zn = -sqrt(1 - n);
                p_radius = 1;
            }

            // RA - hour angle
            const double rotang1 = -M_PI / 2 + lng_sun;
            const double rotang2 = -M_PI / 2 - lat_sun;

            // Tip (y,z) by declination
            const double x2 = xn;
            const double y2 = yn * cos(rotang2) + zn * sin(rotang2);
            const double z2 = -yn * sin(rotang2) + zn * cos(rotang2);

            // Tip (x,y) to get the right central longitude
            const double xf = x2 * cos(rotang1) - y2 * sin(rotang1);
            const double yf = x2 * sin(rotang1) + y2 * cos(rotang1);
            const double zf = z2;

            // Coordinates now orientated with z-axis pointing out of north pole
            // longitude
            double lng = atan2(yf, xf) / M_PI * 180;
            double lat = asin(zf) / M_PI * 180;

            double shadow = getShadowFraction(lat, lng, jd, p_radius, pos_sun, pos_moon, pos_earth, sidereal_time);

            if (p_radius > 1) {
                lng = lat = GSL_NAN;
            }

            const int offset = x + y * config->x_size_3d;
            shadow_map->map[offset] = shadow;
            shadow_map->lat[offset] = lat;
            shadow_map->lng[offset] = lng;
        }
    return shadow_map;
}