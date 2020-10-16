// shadow_calc.c

// -------------------------------------------------
// Copyright 2019-2020 Dominic Ford.

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
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "constants.h"
#include "settings.h"
#include "shadow_calc.h"

#include "mathsTools/julianDate.h"
#include "mathsTools/sphericalAst.h"
#include "projection.h"

/**
 * Return the quantity delta_T at epoch JD (seconds)
 *
 * delta_T = TT - UT
 * TT = UT + delta_T
 * UT = TT - delta_T
 *
 * @param JD [in] - Julian day at which to evaluate delta_T
 * @return delta_T, measured in seconds
 */
double delta_t(double JD) {
    // See https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    const double y = (JD - 1721059.5) / 365.2425;

    if (y < -500) {
        const double u = (y - 1820) / 100;
        return -20 + 32 * pow(u, 2);
    }

    if (y < 500) {
        const double u = y / 100;
        return 10583.6 - 1014.41 * u + 33.78311 * pow(u, 2) - 5.952053 * pow(u, 3) - 0.1798452 * pow(u, 4) +
               0.022174192 * pow(u, 5) + 0.0090316521 * pow(u, 6);
    }

    if (y < 1600) {
        const double u = (y - 1000) / 100;
        return 1574.2 - 556.01 * u + 71.23472 * pow(u, 2) + 0.319781 * pow(u, 3)
               - 0.8503463 * pow(u, 4) - 0.005050998 * pow(u, 5) + 0.0083572073 * pow(u, 6);
    }

    if (y < 1700) {
        const double t = y - 1600;
        return 120 - 0.9808 * t - 0.01532 * pow(t, 2) + pow(t, 3) / 7129;
    }

    if (y < 1800) {
        const double t = y - 1700;
        return 8.83 + 0.1603 * t - 0.0059285 * pow(t, 2) + 0.00013336 * pow(t, 3) - pow(t, 4) / 1174000;
    }

    if (y < 1860) {
        const double t = y - 1800;
        return 13.72 - 0.332447 * t + 0.0068612 * pow(t, 2) + 0.0041116 * pow(t, 3) - 0.00037436 * pow(t, 4)
               + 0.0000121272 * pow(t, 5) - 0.0000001699 * pow(t, 6) + 0.000000000875 * pow(t, 7);
    }

    if (y < 1900) {
        const double t = y - 1860;
        return 7.62 + 0.5737 * t - 0.251754 * pow(t, 2) + 0.01680668 * pow(t, 3)
               - 0.0004473624 * pow(t, 4) + pow(t, 5) / 233174;
    }

    if (y < 1920) {
        const double t = y - 1900;
        return -2.79 + 1.494119 * t - 0.0598939 * pow(t, 2) + 0.0061966 * pow(t, 3) - 0.000197 * pow(t, 4);
    }

    if (y < 1941) {
        const double t = y - 1920;
        return 21.20 + 0.84493 * t - 0.076100 * pow(t, 2) + 0.0020936 * pow(t, 3);
    }

    if (y < 1961) {
        const double t = y - 1950;
        return 29.07 + 0.407 * t - pow(t, 2) / 233 + pow(t, 3) / 2547;
    }

    if (y < 1986) {
        const double t = y - 1975;
        return 45.45 + 1.067 * t - pow(t, 2) / 260 - pow(t, 3) / 718;
    }

    if (y < 2005) {
        const double t = y - 2000;
        return 63.86 + 0.3345 * t - 0.060374 * pow(t, 2) + 0.0017275 * pow(t, 3) + 0.000651814 * pow(t, 4)
               + 0.00002373599 * pow(t, 5);
    }

    if (y < 2050) {
        const double t = y - 2000;
        return 62.92 + 0.32217 * t + 0.005589 * pow(t, 2);
    }

    if (y < 2150) {
        return -20 + 32 * pow(((y - 1820) / 100), 2) - 0.5628 * (2150 - y);
    }

    const double u = (y - 1820) / 100;
    return -20 + 32 * pow(u, 2);
}

/**
 * siderealTime - Turns a Julian date into a sidereal time (in hours, at Greenwich)
 * @param JD [in] - Julian date (TT)
 * @return - Sidereal time in radians
 */
double siderealTime(double JD) {
    const double UT = JD - delta_t(JD) / 86400.;
    // printf("Delta T = %.1f\n", delta_t(JD));
    const double T = (UT - 2451545.0) / 36525.0; // See pages 87-88 of Astronomical Algorithms, by Jean Meeus
    return fmod(M_PI / 180 * (
            280.46061837 +
            360.98564736629 * (UT - 2451545.0) +
            0.000387933 * T * T -
            T * T * T / 38710000.0
    ), 2 * M_PI);
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

    // In radians, the geodetic coordinates of the requested point, rotated to place RA=0 at longitude 0
    const double lat_geodetic = lat * M_PI / 180;
    const double lng_geodetic = (lng + sidereal_time) * M_PI / 180;

    // Position in WGS84 coordinate system
    // See <https://en.wikipedia.org/wiki/Reference_ellipsoid>
    const double altitude = 0;  // metres

    const double n = gsl_pow_2(RADIUS_EARTH_EQUATOR) / sqrt(gsl_pow_2(RADIUS_EARTH_EQUATOR * cos(lat_geodetic)) +
                                                            gsl_pow_2(RADIUS_EARTH_POLE * sin(lat_geodetic)));
    const double x = (n + altitude) * cos(lng_geodetic) * cos(lat_geodetic);
    const double y = (n + altitude) * sin(lng_geodetic) * cos(lat_geodetic);
    const double z = (gsl_pow_2(RADIUS_EARTH_POLE / RADIUS_EARTH_EQUATOR) * n + altitude) * sin(lat_geodetic);

    // Work out RA and Dec of star above this point, for ecliptic of epoch
    const double radius_geoid = sqrt(gsl_pow_2(x) + gsl_pow_2(y) + gsl_pow_2(z));  // metres

    const double lat_at_epoch = asin(z / radius_geoid);  // planetocentric coordinates; radians
    const double lng_at_epoch = atan2(y, x);  // planetocentric coordinates; radians

    // Transform into J2000.0
    double lat_j2000, lng_j2000; // radians
    ra_dec_to_j2000(lng_at_epoch * 12. / M_PI, lat_at_epoch * 180 / M_PI, unix_from_jd(epoch), &lng_j2000, &lat_j2000);
    lng_j2000 *= M_PI / 12;
    lat_j2000 *= M_PI / 180;

    // Output position relative to the solar system barycentre, in ICRF coordinates
    const double radius_requested = radius_geoid * radius_in_earth_radii / AU;  // AU
    out[0] = cos(lat_j2000) * cos(lng_j2000) * radius_requested + pos_earth[0];
    out[1] = cos(lat_j2000) * sin(lng_j2000) * radius_requested + pos_earth[1];
    out[2] = sin(lat_j2000) * radius_requested + pos_earth[2];
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

    output->jd_partial_start = (double *) calloc(x_size * y_size, sizeof(double));
    output->jd_total_start = (double *) calloc(x_size * y_size, sizeof(double));
    output->jd_total_end = (double *) calloc(x_size * y_size, sizeof(double));
    output->jd_partial_end = (double *) calloc(x_size * y_size, sizeof(double));

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
    free(item->jd_partial_start);
    free(item->jd_total_start);
    free(item->jd_total_end);
    free(item->jd_partial_end);
    free(item);
}

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

            // Work out whether the Sun is above the horizon at this point
            const double ang_dist_sun = angDist_RADec(lng * M_PI / 180, lat * M_PI / 180, lng_sun, lat_sun);

            const int night_time = (ang_dist_sun > M_PI / 2);

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

                // Update the output structure which details the magnitude and place where greatest eclipse occurs
                if (shadow > span_output->greatest_eclipse_magnitude) {
                    span_output->greatest_eclipse_magnitude = shadow;  // In range 0-1
                    span_output->greatest_eclipse_latitude = lat;  // degrees
                    span_output->greatest_eclipse_longitude = lng;  // degrees
                }
            }

            // Keep track of the earliest and latest times when any partial eclipse is visible
            if (shadow > 0) {

                // Update time span of partial eclipse in this pixel
                if (greatest_shadow->jd_partial_start[offset] == 0) {
                    greatest_shadow->jd_partial_start[offset] = jd;
                }

                // Update time span of partial eclipse in this pixel
                greatest_shadow->jd_partial_end[offset] = jd;

                // Update global record of time span of partial eclipse
                if ((span_output->partial_start <= 0) || (span_output->partial_start > jd)) {
                    span_output->partial_start = jd;
                }
                if ((span_output->partial_end <= 0) || (span_output->partial_end < jd)) {
                    span_output->partial_end = jd;
                }
            }

            // Keep track of the earliest and latest times when any total (or annular) eclipse is visible
            if (shadow > 0.999) {
                // Update time span of total eclipse in this pixel
                if (greatest_shadow->jd_total_start[offset] == 0) {
                    greatest_shadow->jd_total_start[offset] = jd;
                }

                // Update time span of total eclipse in this pixel
                greatest_shadow->jd_total_end[offset] = jd;

                // Update global record of time span of total eclipse
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
            double lng, lat, p_radius;
            project_3d(config, x, y, lng_sun, lat_sun, &lng, &lat, &p_radius);

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
