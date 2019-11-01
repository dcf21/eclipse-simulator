// map_greatest_eclipse.c

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
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "coreUtils/errorReport.h"
#include "mathsTools/julianDate.h"
#include "mathsTools/sphericalAst.h"
#include "constants.h"
#include "duration.h"
#include "ephemeris.h"
#include "map_greatest_eclipse.h"
#include "settings.h"
#include "shadow_calc.h"

/**
 * parameters - structure for passing variables to the function my_f, which evaluates the perpendicular distance of
 * a test point from the Sun-Moon line.
 */
typedef struct parameters {
    const double *pos_sun;
    const double *pos_earth;
    const double *pos_moon;
    double JD; // TT
    double sidereal_time;
    double lat_sun, lng_sun;
    double prev_sun_ang_dist;
} parameters;

/**
 * Subtract two three-component vectors.
 * @param out [out] - Output a-b
 * @param a [in] - Vector a
 * @param b [in] - Vector b
 * @return Pointer to out
 */
double *subtract(double *out, const double *a, const double *b) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
    return out;
}

/**
 * Calculate the magnitude of a three-component vector.
 * @param a [in] - Vector a
 * @return Length of vector a
 */
double magnitude(const double *a) {
    return sqrt(gsl_pow_2(a[0]) + gsl_pow_2(a[1]) + gsl_pow_2(a[2]));
}

double point_line_distance(const double earth_surface[3], const double pos_sun[3], const double pos_moon[3]) {
    const double sun_moon_vector[3] = {pos_moon[0] - pos_sun[0],
                                       pos_moon[1] - pos_sun[1],
                                       pos_moon[2] - pos_sun[2]
    };

    const double dist_sun_moon = magnitude(sun_moon_vector);
    const double sun_moon_unit_vector[3] = {
            (pos_moon[0] - pos_sun[0]) / dist_sun_moon,
            (pos_moon[1] - pos_sun[1]) / dist_sun_moon,
            (pos_moon[2] - pos_sun[2]) / dist_sun_moon
    };

    const double vector_moon_earth[3] = {
            earth_surface[0] - pos_moon[0],
            earth_surface[1] - pos_moon[1],
            earth_surface[2] - pos_moon[2],
    };

    const double dot_product = (vector_moon_earth[0] * sun_moon_unit_vector[0] +
                                vector_moon_earth[1] * sun_moon_unit_vector[1] +
                                vector_moon_earth[2] * sun_moon_unit_vector[2]);

    const double closest_point_on_line[3] = {
            pos_moon[0] + dot_product * sun_moon_unit_vector[0],
            pos_moon[1] + dot_product * sun_moon_unit_vector[1],
            pos_moon[2] + dot_product * sun_moon_unit_vector[2],
    };

    const double perpendicular_distance = sqrt(
            gsl_pow_2(earth_surface[0] - closest_point_on_line[0]) +
            gsl_pow_2(earth_surface[1] - closest_point_on_line[1]) +
            gsl_pow_2(earth_surface[2] - closest_point_on_line[2])
    );

    return perpendicular_distance;
}

/**
 * Objective function to minimise. This evaluates the perpendicular distance of a test point from the Sun-Moon line.
 * It returns NaN if the test-point is round the back of the Earth as seen from the Sun, where there is a second point
 * which falls on the Sun-Moon line, but which doesn't experience a solar eclipse owing to it being night time.
 * @param v [in] - Vector a = (lng, lat) in radians
 * @param params [in] - Null pointer to a structure of type <parameters>
 * @return Perpendicular distance of a test point from the Sun-Moon line
 */
// Objective function to minimise
double my_f(const gsl_vector *v, void *params) {
    // Extract parameters from null pointer
    const parameters *p = (parameters *) params;

    // Extract longitude and latitude of test point, in radians
    const double lng = gsl_vector_get(v, 0);
    const double lat = gsl_vector_get(v, 1);

    // Check latitude is within legitimate range
    if (lat < -M_PI / 2) return GSL_NAN;
    if (lat > M_PI / 2) return GSL_NAN;

    // Check angular distance of point from where Sun is overhead
    const double ang_dist_sun = angDist_RADec(lng, lat, p->lng_sun, p->lat_sun);

    // If point is on night side of the Earth, it is not allowed
    if (ang_dist_sun > M_PI / 2) {
        return GSL_NAN;
    }

    // If point is much further from where Sun is overhead than previous point, we're in danger of flipping round to
    // the far side of the Earth...
    if (gsl_finite(p->prev_sun_ang_dist) && (ang_dist_sun > p->prev_sun_ang_dist + 5 * M_PI / 180)) {
        return GSL_NAN;
    }

    // Project point into Cartesian coordinates
    double pos[3];
    earthTopocentricPositionICRF(pos, lat * 180 / M_PI, lng * 180 / M_PI, 1,
                                 p->pos_earth, p->JD, p->sidereal_time * 180 / M_PI);

    // Calculate perpendicular distance of point from the Sun-Moon line
    const double perpendicular_dist = point_line_distance(pos, p->pos_sun, p->pos_moon);

    return perpendicular_dist;
}

/**
 * eclipse_duration_from_path - Look up the maximum duration of the eclipse at the central point, at a given time.
 * @param paths [in] - The path of the eclipse, as returned by <map_greatest_eclipse>, which includes duration data
 * @param jd [in] - The time point to query, specified as a Julian day number (TT)
 * @param is_total [out] - Flag indicating whether the eclipse is total or annular
 * @return The duration of the eclipse at this time point, in seconds
 */
double eclipse_duration_from_path(const eclipse_path_list *paths, double jd, int *is_total) {
    for (int i = 0; i < paths->path_count; i++) {
        const eclipse_path *path = &paths->paths[i];
        for (int j = 1; j < path->point_count; j++) {
            const path_point *p0 = &path->path[j - 1];
            const path_point *p1 = &path->path[j];

            if ((jd >= p0->jd) && (jd <= p1->jd)) {
                // Weighted linear interpolation
                const double w0 = p1->jd - jd;
                const double w1 = jd - p0->jd;
                const double duration = (w0 * p0->duration + w1 * p1->duration) / (w0 + w1);

                *is_total = path->is_total;
                return duration;
            }
        }
    }

    *is_total = 0;
    return 0;
}

/**
 * Loop over an ephemeris, producing a text file listing the position of greatest eclipse on the Earth's surface over
 * time. We use a GSL numerical minimiser to find the position of the Earth's surface (oblate spheroid) which is closest
 * to the Sun-Moon line.
 * @param config [in] - settings
 * @param ephemeris [in] - an ephemeris computed using <ephemerisCompute>
 * @return A list of the path segments of total and/or annular periods in this eclipse
 */
eclipse_path_list *map_greatest_eclipse(const settings *config, const ephemeris *ephemeris) {
    // Open file to write path of greatest eclipse to
    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/maximumEclipsePath.json", config->output_dir);
    FILE *f = fopen(fname, "wb");

    // Create data structure for storing the path of greatest eclipse
    eclipse_path_list *paths = (eclipse_path_list *) malloc(sizeof(eclipse_path_list));
    paths->path_count = 0;

    // Work out a critical (angular) distance across the Earth's surface that maximum eclipse needs to move before
    // we write out a new data point. We use this to reduce the number of data points that we write to disk.
    const double earth_circumference = 2 * M_PI * RADIUS_EARTH_EQUATOR;
    const double critical_distance = 20e3; // metres
    const double critical_ang_dist = (critical_distance / earth_circumference) * 2 * M_PI;

    // The last point along the path that we wrote to disk. Use this to work out when the path has moved far enough
    // that we need to write out a new point
    const path_point null_point = {GSL_NAN, GSL_NAN, GSL_NAN};
    path_point previous_point_written = null_point;

    // The last point that we calculated. Use this to write out one final data point to disk when the path ends.
    path_point previous_point_calculated = null_point;
    int previous_point_total = 0;
    double prev_sun_ang_dist = GSL_NAN;

    // How many samples are there in the input ephemeris between those we analyse
    const int eclipse_path_search_stride = (int) round(config->eclipse_path_search_time_resolution /
                                                       config->time_resolution);

    // Loop over video frames
    int frame_counter = -1;
    for (int j = 0; j < ephemeris->point_count; j += eclipse_path_search_stride) {
        frame_counter++;

        // Create data structure for passing to the function we are going to minimise
        parameters p;

        // Read ephemeris
        p.JD = ephemeris->jd_start + ephemeris->jd_step * j;
        p.pos_sun = ephemeris->data[j].sun_pos;
        p.pos_earth = ephemeris->data[j].earth_pos;
        p.pos_moon = ephemeris->data[j].moon_pos;
        p.prev_sun_ang_dist = prev_sun_ang_dist;

        // Calculate the latitude and longitude on Earth where the Sun is overhead
        calculate_where_sun_overhead(&p.lat_sun, &p.lng_sun, &p.sidereal_time, p.pos_sun, p.pos_earth, p.JD);

        // Create a GSL function minimiser to find the point on the surface of the Earth closest to the Sun-Moon line
        const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer *s = NULL;
        gsl_vector *ss, *x;
        gsl_multimin_function minex_func;

        size_t iter = 0;
        int status;
        double size;

        // Starting point
        x = gsl_vector_alloc(2);
        if (gsl_finite(previous_point_calculated.jd) && (prev_sun_ang_dist * 180. / M_PI < 85)) {
            gsl_vector_set(x, 0, previous_point_calculated.longitude);
            gsl_vector_set(x, 1, previous_point_calculated.latitude);
        } else {
            gsl_vector_set(x, 0, p.lng_sun);
            gsl_vector_set(x, 1, p.lat_sun);
        }

        // Set initial step sizes to 0.001 radians in latitude and longitude
        ss = gsl_vector_alloc(2);
        gsl_vector_set_all(ss, 0.001);

        // Initialize method and iterate
        minex_func.n = 2;
        minex_func.f = my_f;
        minex_func.params = &p;

        s = gsl_multimin_fminimizer_alloc(T, 2);
        gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

        do {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);

            if (status)
                break;

            size = gsl_multimin_fminimizer_size(s);
            status = gsl_multimin_test_size(size, 1e-8);
        } while (status == GSL_CONTINUE && iter < 5000);

        // Look up latitude and longitude of minimisation result
        double lng = gsl_vector_get(s->x, 0);
        double lat = gsl_vector_get(s->x, 1);

        while (lng <= -M_PI) lng += 2 * M_PI;
        while (lng > M_PI) lng -= 2 * M_PI;

        // Store position of midpoint of eclipse - this is used as the starting point when we are mapping contours
        if (frame_counter == ephemeris->point_count / eclipse_path_search_stride / 2) {
            paths->latitude_midpoint = lat;
            paths->longitude_midpoint = lng;
        }

        // Test whether this point lies on Earth-Moon line
        const double perpendicular_distance = my_f(s->x, &p); // distance in AU

        // Free up resources
        gsl_vector_free(x);
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free(s);

        // Update <prev_sun_ang_dist>
        const double sun_ang_dist = angDist_RADec(lng, lat, p.lng_sun, p.lat_sun);
        prev_sun_ang_dist = sun_ang_dist;

        // If point not on Earth-Moon line, then total eclipse has moved off the Earth's surface
        if ((!gsl_finite(perpendicular_distance)) || (perpendicular_distance > 0.001 * RADIUS_EARTH_EQUATOR / AU)) {
            // Write one final point at the end of the line
            if (gsl_finite(previous_point_calculated.jd) &&
                (previous_point_calculated.jd != previous_point_written.jd)) {
                eclipse_path *path = &paths->paths[paths->path_count - 1];
                path->jd_end = previous_point_calculated.jd;
                path->path[path->point_count] = previous_point_calculated;
                path->point_count++;
                if (path->point_count >= MAX_PATH_LENGTH) ephem_fatal(__FILE__, __LINE__, "Path too long");
            }

            // Break line
            previous_point_written = null_point;
            previous_point_calculated = null_point;
            // fprintf(f,
            //         "## perpendicular_distance = %10.4f (iter %5ld lng %10.5f  lat %10.5f sun_dist %10.5f lng_sun %10.5f  lat_sun %10.5f)\n",
            //         perpendicular_distance / RADIUS_EARTH_EQUATOR * AU,
            //         iter,
            //         lng * 180 / M_PI, lat * 180 / M_PI,
            //         sun_ang_dist * 180 / M_PI,
            //         p.lng_sun * 180 / M_PI, p.lat_sun * 180 / M_PI);
            continue;
        }

        // Project point into Cartesian coordinates
        double pos[3], p0[3], p1[3];
        earthTopocentricPositionICRF(pos, lat * 180 / M_PI, lng * 180 / M_PI, 1,
                                     p.pos_earth, p.JD, p.sidereal_time * 180 / M_PI);

        // Test whether greatest eclipse is a total eclipse or an annular eclipse
        const double moon_sun_dist_ratio = (magnitude(subtract(p0, p.pos_moon, pos)) /
                                            magnitude(subtract(p1, p.pos_sun, pos)));
        const double moon_sun_size_ratio = RADIUS_MOON / RADIUS_SUN;
        const double moon_sun_ang_size_ratio = moon_sun_size_ratio / moon_sun_dist_ratio;

        const int total_eclipse = (moon_sun_ang_size_ratio >= 1);

        // If we have moved from total eclipse to annular eclipse, break the path of greatest eclipse
        if (total_eclipse != previous_point_total) {
            // Write one final point at the end of the line
            if (gsl_finite(previous_point_calculated.jd) &&
                (previous_point_calculated.jd != previous_point_written.jd)) {
                eclipse_path *path = &paths->paths[paths->path_count - 1];
                path->jd_end = previous_point_calculated.jd;
                path->path[path->point_count] = previous_point_calculated;
                path->point_count++;
                if (path->point_count >= MAX_PATH_LENGTH) ephem_fatal(__FILE__, __LINE__, "Path too long");
            }

            previous_point_written = null_point;
            previous_point_calculated = null_point;
            previous_point_total = total_eclipse;
        }

        // Test if we should write this point
        int write_point = 0;
        if (!gsl_finite(previous_point_written.jd)) {
            // If previous point was not on Earth's surface, start new line
            paths->path_count++;
            if (paths->path_count >= MAX_PATH_ITEMS) ephem_fatal(__FILE__, __LINE__, "Path count overflow");
            eclipse_path *path = &paths->paths[paths->path_count - 1];
            path->jd_start = p.JD;
            path->is_total = total_eclipse;
            path->point_count = 0;
            // fprintf(f, "\n");
            write_point = 1;
        } else {
            // If we've moved more than critical distance since last point was written, it's time to write a new point
            const double ang_distance = angDist_RADec(previous_point_written.longitude, previous_point_written.latitude,
                                                      lng, lat);
            write_point = (ang_distance >= critical_ang_dist);
        }

        // Update <previous_point_calculated>
        previous_point_calculated.latitude = lat;
        previous_point_calculated.longitude = lng;
        previous_point_calculated.jd = p.JD;

        // Write this point, if required
        if (write_point) {
            // Write this point to <paths>
            eclipse_path *path = &paths->paths[paths->path_count - 1];
            path->jd_end = previous_point_calculated.jd;
            path->path[path->point_count] = previous_point_calculated;
            path->point_count++;
            if (path->point_count >= MAX_PATH_LENGTH) ephem_fatal(__FILE__, __LINE__, "Path too long");

            // Update <previous_point_written>
            previous_point_written.latitude = lat;
            previous_point_written.longitude = lng;
            previous_point_written.jd = p.JD;

            // Write this path to text file
            // fprintf(f,
            //         "%16.1f  %10.5f  %10.5f  %d",
            //         unix_from_jd(p.JD), lng * 180 / M_PI, lat * 180 / M_PI, total_eclipse
            // );
            // fprintf(f,
            //         " # iter %5ld sun_dist %10.5f",
            //         iter, sun_ang_dist * 180 / M_PI
            // );
            // fprintf(f, "\n");
        }
    }

    // Write path out to file
    double maximum_duration = 0;
    path_point *previous_point = NULL;

    fprintf(f, "{\"paths\":[");

    for (int i = 0; i < paths->path_count; i++) {
        eclipse_path *path = &paths->paths[i];

        if (i > 0) fprintf(f, ",");
        fprintf(f, "[%d,[", path->is_total);
        int first_point = 1;

        if (previous_point != NULL) {
            if (!first_point) fprintf(f, ",");
            // Unix timestamp (UT; not TT), longitude/deg, latitude/deg, duration/sec
            fprintf(f,
                    "[%.1f,%.5f,%.5f,%.1f]",
                    unix_from_jd(previous_point->jd) - delta_t(previous_point->jd),
                    previous_point->longitude * 180 / M_PI, previous_point->latitude * 180 / M_PI,
                    0.
            );
            first_point = 0;
        }

        for (int j = 0; j < path->point_count; j++) {
            path_point *p = &path->path[j];
            const double duration = eclipse_duration(config, ephemeris, p->longitude, p->latitude, p->jd,
                                                     path->is_total);
            p->duration = duration;

            if (duration > maximum_duration) maximum_duration = duration;

            if (!first_point) fprintf(f, ",");
            // Unix timestamp (UT; not TT), longitude/deg, latitude/deg, duration/sec
            fprintf(f,
                    "[%.1f,%.5f,%.5f,%.1f]",
                    unix_from_jd(p->jd) - delta_t(p->jd),
                    p->longitude * 180 / M_PI, p->latitude * 180 / M_PI,
                    duration
            );

            first_point = 0;
            previous_point = p;
        }
        fprintf(f, "]]");
    }

    // Close output file
    fprintf(f, "],\n\"path_segments\":%d,\"duration\":%.1f}\n", paths->path_count, maximum_duration);
    fclose(f);

    // Return paths we computed
    return paths;
};
