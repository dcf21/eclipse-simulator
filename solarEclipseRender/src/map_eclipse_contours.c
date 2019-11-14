// map_eclipse_contours.c

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

#include "coreUtils/errorReport.h"
#include "mathsTools/sphericalAst.h"

#include "constants.h"
#include "ephemeris.h"
#include "map_greatest_eclipse.h"
#include "map_eclipse_contours.h"
#include "settings.h"
#include "shadow_calc.h"

/**
 * objective_function - Define the objective function that we are plotting contours of, and the various inputs that
 * it requires
 */
typedef struct objective_function {
    int (*func)(const struct objective_function *, int, int);

    const settings *config;
    const ephemeris *ephemeris;
    const time_span *eclipse_time_span; // The time span of the partial and total phases of the eclipse.
    const double contour_tracing_grid_resolution; // degrees. The steps in lat and lng that we take to trace contours.
    const double contour_tracing_critical_distance; // metres. Output points along each contour at this specing.
    const double contour_tracing_time_resolution; // seconds. Check eclipse magnitude with this time resolution.
    const double contour_level; // fraction of Sun's disk covered; 0-1
} objective_function;

/**
 * trace_contour_by_eclipse_level - Trace a contour of the eclipse magnitude based on the fraction of the Sun's disk
 * covered by the Moon.
 *
 * @param f [in] - Settings for the contour tracing process
 * @param x [in] - The x pixel in the grid that we are sampling. Perhaps in units of hundredths of a degree.
 * @param y [in] - The x pixel in the grid that we are sampling. Perhaps in units of hundredths of a degree.
 * @return - Boolean indicating whether this point is inside contour
 */
int trace_contour_by_eclipse_level(const objective_function *f, int x, int y) {
    const ephemeris *ephemeris = f->ephemeris;

    // Calculate the first and last point numbers in the ephemeris structure that we need to analyse
    int point_start = (int) floor((f->eclipse_time_span->partial_start - f->ephemeris->jd_start) /
                                  f->ephemeris->jd_step);
    int point_end = (int) ceil((f->eclipse_time_span->partial_end - f->ephemeris->jd_start) /
                               f->ephemeris->jd_step);

    // Check that limits are within range
    if (point_start < 0) point_start = 0;
    if (point_end > ephemeris->point_count) point_end = ephemeris->point_count;

    // Convert pixel number into a longitude and latitude in degrees
    const double longitude_deg = x * f->contour_tracing_grid_resolution;
    const double latitude_deg = y * f->contour_tracing_grid_resolution;

    if ((latitude_deg <= -90) || (latitude_deg >= 90)) return 0;

    const int j_step = (int) (f->contour_tracing_time_resolution / f->config->time_resolution);

    // loop over the ephemeris points that we are to analyse
    double maximum_shadow_fraction = 0;
    for (int j = point_start; j < point_end; j += j_step) {
        const double jd = ephemeris->jd_start + ephemeris->jd_step * j;
        const ephemeris_point p = ephemeris->data[j];

        // Calculate the latitude and longitude on Earth where the Sun is overhead
        double sidereal_time, lat_sun, lng_sun; // all in radians
        calculate_where_sun_overhead(&lat_sun, &lng_sun, &sidereal_time, p.sun_pos, p.earth_pos, jd);

        // Shadow fraction is zero if the Sun is below the horizon
        const double sun_ang_dist = angDist_RADec(longitude_deg * M_PI / 180, latitude_deg * M_PI / 180, lng_sun,
                                                  lat_sun);
        if (sun_ang_dist >= M_PI / 2) continue;

        // Work out the fraction of the Sun covered by the Moon
        const double shadow_fraction = getShadowFraction(latitude_deg, longitude_deg, jd, 1,
                                                         p.sun_pos, p.moon_pos, p.earth_pos,
                                                         sidereal_time);

        // Keep track of maximum eclipse
        if (shadow_fraction > maximum_shadow_fraction) maximum_shadow_fraction = shadow_fraction;
    }

    // printf("Testing (%10.5f,%10.5f)...   %.4f\n", longitude_deg, latitude_deg, maximum_shadow_fraction);
    return maximum_shadow_fraction >= f->contour_level;
}

/**
 * trace_contour_by_annularity - Trace a contour of the eclipse magnitude based on whether an annular eclipse is
 * visible.
 *
 * @param f [in] - Settings for the contour tracing process
 * @param x [in] - The x pixel in the grid that we are sampling. Perhaps in units of hundredths of a degree.
 * @param y [in] - The x pixel in the grid that we are sampling. Perhaps in units of hundredths of a degree.
 * @return - Boolean indicating whether this point is inside contour
 */
int trace_contour_by_annularity(const objective_function *f, int x, int y) {
    const ephemeris *ephemeris = f->ephemeris;

    // Calculate the first and last point numbers in the ephemeris structure that we need to analyse
    int point_start = (int) floor((f->eclipse_time_span->partial_start - f->ephemeris->jd_start) /
                                  f->ephemeris->jd_step);
    int point_end = (int) ceil((f->eclipse_time_span->partial_end - f->ephemeris->jd_start) /
                               f->ephemeris->jd_step);

    // Check that limits are within range
    if (point_start < 0) point_start = 0;
    if (point_end > ephemeris->point_count) point_end = ephemeris->point_count;

    // Convert pixel number into a longitude and latitude in degrees
    const double longitude_deg = x * f->contour_tracing_grid_resolution;
    const double latitude_deg = y * f->contour_tracing_grid_resolution;

    if ((latitude_deg <= -90) || (latitude_deg >= 90)) return 0;

    const int j_step = (int) (f->contour_tracing_time_resolution / f->config->time_resolution);

    // loop over the ephemeris points that we are to analyse
    for (int j = point_start; j < point_end; j += j_step) {
        const double jd = ephemeris->jd_start + ephemeris->jd_step * j;
        const ephemeris_point p = ephemeris->data[j];

        // Calculate the latitude and longitude on Earth where the Sun is overhead
        double sidereal_time, lat_sun, lng_sun; // all in radians
        calculate_where_sun_overhead(&lat_sun, &lng_sun, &sidereal_time, p.sun_pos, p.earth_pos, jd);

        // Shadow fraction is zero if the Sun is below the horizon
        const double sun_ang_dist = angDist_RADec(longitude_deg * M_PI / 180, latitude_deg * M_PI / 180, lng_sun,
                                                  lat_sun);
        if (sun_ang_dist >= M_PI / 2) continue;

        // For annular eclipses, look at whether the Moon is contained within the Sun's disk
        const int is_annular = testIfAnnularEclipse(latitude_deg, longitude_deg, jd, 1,
                                                    p.sun_pos, p.moon_pos, p.earth_pos,
                                                    sidereal_time);

        if (is_annular) {
            // printf("Testing (%10.5f,%10.5f)...   IN\n", longitude_deg, latitude_deg);
            return 1;
        }
    }
    // printf("Testing (%10.5f,%10.5f)...   OUT\n", longitude_deg, latitude_deg);
    return 0; // no annular eclipse found
}

/**
 * test_pixel_in_direction - From the latest point along the contour, look one square to the side to see whether that
 * neighbouring square is inside or outside the contour.
 *
 * @param f [in] - Settings for the contour tracing process
 * @param x_pos [in] - The x pixel in the grid of the latest point along the contour. Perhaps in units of hundredths
 * of a degree.
 * @param y_pos [in] - The y pixel in the grid of the latest point along the contour. Perhaps in units of hundredths
 * of a degree.
 * @param direction [in] - The direction we are currently facing. 0=right; 1=up; 2=left; 3=down
 * @return - Boolean indicating whether the pixel being looked at is inside the contour
 */
int test_pixel_in_direction(const objective_function *f, int x_pos, int y_pos, int direction) {
    // looking right
    if (direction == 0) {
        if (f->func(f, x_pos + 1, y_pos) > 0) return 1;
        return 0;
    }

    // looking up
    if (direction == 1) {
        if (f->func(f, x_pos, y_pos - 1) > 0) return 1;
        return 0;
    }

    // looking left
    if (direction == 2) {
        if (f->func(f, x_pos - 1, y_pos) > 0) return 1;
        return 0;
    }

    // looking down
    if (direction == 3) {
        if (f->func(f, x_pos, y_pos + 1) > 0) return 1;
        return 0;
    }

    // illegal direction
    logging_fatal(__FILE__, __LINE__, "Illegal direction");
    exit(1);
}

/**
 * line_status - Structure used when tracing a contour to keep track of the last point which we wrote to disk.
 */
typedef struct line_status {
    FILE *output;
    double lat_prev, lng_prev;
} line_status;

/**
 * register_point - Add a new point along the current contour. We may or may not write the point out to disk, depending
 * whether it is more than a critical distance away from the last point that we wrote out
 *
 * @param f [in] - Settings for the contour tracing process
 * @param s - The current line_status for this contour
 * @param x [in] - The position in longitude of the new point along the contour
 * @param y [in] - The position in latitude of the new point along the contour
 */
void register_point(const objective_function *f, line_status *s, int x, int y) {
    // Work out a critical (angular) distance across the Earth's surface that maximum eclipse needs to move before
    // we write out a new data point. We use this to reduce the number of data points that we write to disk.
    const double earth_circumference = 2 * M_PI * RADIUS_EARTH_EQUATOR;
    const double critical_distance = f->contour_tracing_critical_distance; // metres
    const double critical_ang_dist = (critical_distance / earth_circumference) * 2 * M_PI;

    // Convert pixel in the grid we are exploring into a latitude and longitude (in radians)
    double longitude = x * f->contour_tracing_grid_resolution * M_PI / 180;
    double latitude = y * f->contour_tracing_grid_resolution * M_PI / 180;

    // Test if we should write this point
    int write_point = 0;
    if (!gsl_finite(s->lng_prev)) {
        write_point = 1;
    } else {
        // If we've moved more than critical distance since last point was written, it's time to write a new point
        const double ang_distance = angDist_RADec(s->lng_prev, s->lat_prev, longitude, latitude);
        write_point = (ang_distance >= critical_ang_dist);
    }

    if (write_point) {
        if (gsl_finite(s->lng_prev)) fprintf(s->output, ",");
        s->lng_prev = longitude;
        s->lat_prev = latitude;
        fprintf(s->output, "[%.2f,%.2f]", longitude * 180 / M_PI, latitude * 180 / M_PI);
        fflush(s->output);
    }
}

/**
 * trace_contour - Main entry point for tracing a contour of eclipse magnitude.
 * @param output - The file we are to write the contour to
 * @param f [in] - Settings for the contour tracing process
 * @param latitude_start [in] - The latitude, in degrees, of any point which lies inside the contour
 * @param longitude_start [in] - The longitude, in degrees, of any point which lies inside the contour
 */
void trace_contour(FILE *output, const objective_function *f, double latitude_start, double longitude_start) {
    int y_scan;

    // Pixel starting position
    const int centre_x = (int) (longitude_start / f->contour_tracing_grid_resolution);
    const int centre_y_initial = (int) (latitude_start / f->contour_tracing_grid_resolution);
    // printf("Contour centre (%6d, %6d)\n", centre_x, centre_y_initial);

    // Make sure that central starting point is actually inside contour
    int centre_y = (int) (-100 / f->contour_tracing_grid_resolution);
    for (int i = 0; i < 500; i++) {
        if (f->func(f, centre_x, centre_y_initial + i)) {
            centre_y = centre_y_initial + i;
            break;
        }
        if (f->func(f, centre_x, centre_y_initial - i)) {
            centre_y = centre_y_initial - i;
            break;
        }
    }
    // printf("Contour adjusted centre (%6d, %6d)\n", centre_x, centre_y);

    // If we couldn't find a starting point inside the contour, return an empty contour
    if (centre_y < -90 / f->contour_tracing_grid_resolution) {
        fprintf(output, "[]");
        return;
    }

    // Starting position is in the middle of the eclipse, so track up or down in longitude to find edge
    const int direction = (centre_y >= 0) ? -1 : 1;
    y_scan = centre_y;
    while (f->func(f, centre_x, y_scan + direction)) y_scan += direction;

    const int start_x = centre_x;
    const int start_y = y_scan;
    // printf("Contour start point (%6d, %6d) direction %d\n", start_x, start_y, direction);

    // Start tracing the output of a new contour
    fprintf(output, "[");
    line_status s = {output, GSL_NAN, GSL_NAN};

    // Start tracing clockwise around this land mass

    // Direction 0 is right
    // Direction 1 is up
    // Direction 2 is left
    // Direction 3 is down

    int direction_facing = (direction > 0) ? 2 : 0; // start off heading clockwise around contour
    int direction_looking = (direction_facing + 1) % 4; // start off looking up
    int x = start_x;
    int y = start_y;
    const int start_direction = direction_facing;

    // Iterate until we come back to our starting pixel
    int iteration = 0;
    while ((iteration == 0) || (x != start_x) || (y != start_y) ||
           (direction_facing != start_direction)) {
        // Log message
        // printf("Tracing... (%6d,%6d) %2d\n", x, y, direction_facing);

        // Make sure we don't iterate unreasonably long
        iteration++;
        if (iteration > 1e8) {
            logging_fatal(__FILE__, __LINE__, "Contour tracing in infinite loop");
        }

        // If we're looking at a pixel containing land, we change direction to explore it
        const int value_looking_at = test_pixel_in_direction(f, x, y, direction_looking);

        if (value_looking_at) {
            switch (direction_facing) {
                case 0:
                    register_point(f, &s, x, y);
                    break;
                case 1:
                    register_point(f, &s, x, y + 1);
                    break;
                case 2:
                    register_point(f, &s, x + 1, y);
                    break;
                case 3:
                    register_point(f, &s, x + 1, y);
                    break;

            }
            direction_facing = direction_looking;
            direction_looking = (direction_facing + 1) % 4;
        }

        // If we're about to move into a square which is not land
        const int value_moving_into = test_pixel_in_direction(f, x, y, direction_facing);

        if (!value_moving_into) {
            switch (direction_facing) {
                case 1:
                    register_point(f, &s, x, y);
                    break;
                case 2:
                    register_point(f, &s, x, y + 1);
                    break;
                case 3:
                    register_point(f, &s, x + 1, y);
                    break;
                case 0:
                    register_point(f, &s, x + 1, y);
                    break;

            }
            direction_facing = (direction_facing + 3) % 4;
            direction_looking = (direction_facing + 1) % 4;
            continue;
        }

        // If we're looking outside contour, and moving within contour, we move forwards one pixel
        switch (direction_facing) {
            case 0:
                x++;
                break;
            case 1:
                y--;
                break;
            case 2:
                x--;
                break;
            case 3:
                y++;
                break;

        }

        // Check that longitude remains within range
        if (x <= -180 / f->contour_tracing_grid_resolution) x += (int) (360 / f->contour_tracing_grid_resolution);
        if (x > 180 / f->contour_tracing_grid_resolution) x -= (int) (360 / f->contour_tracing_grid_resolution);
    }

    // Finish tracing the outline of this contour
    fprintf(output, "]");

}

/**
 * map_eclipse_contours - Main entry point for tracing contours of an eclipse.
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param ephemeris [in] - an ephemeris computed using <ephemerisCompute>
 * @param paths [in] - The path of this eclipse, as calculated by <map_greatest_eclipse>
 * @param eclipse_time_span [in] - The time span of the partial and total phases of this eclipse, as calculated by
 * <calculate_eclipse_map_2d>
 */
void map_eclipse_contours(const settings *config, const ephemeris *ephemeris,
                          const eclipse_path_list *paths, const time_span *eclipse_time_span) {
    // Open file to write path of greatest eclipse to
    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/maximumEclipseContours.json", config->output_dir);
    FILE *output = fopen(fname, "wb");
    fprintf(output, "[");
    int first = 1;

    // Write contours for total and annular eclipses
    for (int i = 0; i < paths->path_count; i++) {
        const eclipse_path path_segment = paths->paths[i];
        objective_function f = {NULL, config, ephemeris, eclipse_time_span,
                                0.01, 15e3, 0.5, 1};

        if (path_segment.is_total) {
            f.func = trace_contour_by_eclipse_level;
        } else {
            f.func = trace_contour_by_annularity;
        }

        // Calculate the midpoint of the contour, in pixel coordinates
        const int midpoint = path_segment.point_count / 2;
        const double midpoint_latitude = path_segment.path[midpoint].latitude * 180 / M_PI;
        const double midpoint_longitude = path_segment.path[midpoint].longitude * 180 / M_PI;

        // Write output
        fprintf(output, "%s\n", first ? "" : ",");
        first = 0;
        fprintf(output, "[\"%s eclipse\",", path_segment.is_total ? "Total" : "Annular");
        trace_contour(output, &f, midpoint_latitude, midpoint_longitude);
        fprintf(output, "]\n");
    }

    // Write contours at various eclipse magnitudes
    for (double level = 1e-8; level < 0.9; level += 0.2) {
        objective_function f = {trace_contour_by_eclipse_level, config, ephemeris, eclipse_time_span,
                                0.025, 40e3, 4, level};

        // Calculate the midpoint of the contour, in pixel coordinates
        const double midpoint_latitude = eclipse_time_span->greatest_eclipse_latitude;
        const double midpoint_longitude = eclipse_time_span->greatest_eclipse_longitude;

        // Write output
        char contour_title[FNAME_LENGTH];
        if (level > 0.05) {
            sprintf(contour_title, ">%d%% eclipse", (int) (level * 100));
        } else {
            sprintf(contour_title, "Partial eclipse");
        }
        fprintf(output, "%s\n", first ? "" : ",");
        first = 0;
        fprintf(output, "[\"%s\",", contour_title);
        trace_contour(output, &f, midpoint_latitude, midpoint_longitude);
        fprintf(output, "]\n");
    }

    // Close output file
    fprintf(output, "]\n");
    fclose(output);
};
