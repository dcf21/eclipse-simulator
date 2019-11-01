// main.c

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "argparse/argparse.h"

#include "coreUtils/asciiDouble.h"
#include "coreUtils/strConstants.h"
#include "coreUtils/errorReport.h"

#include "country_lookup.h"
#include "ephemeris.h"
#include "jpeg_in.h"
#include "make_binary_map.h"
#include "map_greatest_eclipse.h"
#include "map_eclipse_contours.h"
#include "render_2d.h"
#include "render_3d.h"
#include "settings.h"
#include "shadow_calc.h"

/**
 * usage - Command-line usage instructions for when user supplies the --help switch
 */
static const char *const usage[] = {
        "ephem.bin [options] [[--] args]",
        "ephem.bin [options]",
        NULL,
};

int main(int argc, const char **argv) {
    int j;
    jpeg_ptr earthDay, earthNight;

    // Initialise sub-modules
    if (DEBUG) ephem_log("Initialising eclipse renderer.");

    // Turn off GSL's automatic error handler
    gsl_set_error_handler_off();

    // Default options
    settings config = default_settings();

    // Scan commandline options for any switches
    struct argparse_option options[] = {
            OPT_HELP(),
            OPT_GROUP("Basic options"),
            OPT_FLOAT('a', "jd_min", &config.jd_min,
                      "The Julian day number at which the ephemeris should begin; TT"),
            OPT_FLOAT('b', "jd_max", &config.jd_max,
                      "The Julian day number at which the ephemeris should end; TT"),
            OPT_STRING('t', "title", &config.title,
                       "The title of this solar eclipse event."),
            OPT_STRING('o', "output", &config.output_dir,
                       "The directory to store output in."),

            OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse,
                      "\nEclipse Renderer",
                      "\n");
    argc = argparse_parse(&argparse, argc, argv);

    if (argc != 0) {
        int i;
        for (i = 0; i < argc; i++) {
            printf("Error: unparsed argument <%s>\n", *(argv + i));
        }
        ephem_fatal(__FILE__, __LINE__, "Unparsed arguments");
    }

    // Create output directory
    char cmd[FNAME_LENGTH];
    sprintf(cmd, "mkdir -p %s", config.output_dir);
    system(cmd);

    // Create ephemeris for the Sun. Earth and Moon
    ephemeris *ephemeris;
    fetch_ephemeris(&config, &ephemeris);

    // Make data structure to store time lapse maps of eclipse magnitude across the world
    const int binary_snapshot_count = ceil((config.jd_max - config.jd_min) /
                                           (config.binary_map_time_resolution / 86400));
    unsigned char *binary_eclipse_maps = calloc(
            (int) (360 * config.binary_map_angular_resolution) *
            (int) (180 * config.binary_map_angular_resolution) *
            binary_snapshot_count, 1);

    // Initialise country lookup
    country_lookup_handle *cl = country_lookup_init();

    // Read Earth images
    if (DEBUG) { ephem_log("Reading Earth images."); }
    earthDay = jpeg_get(SRCDIR "/../earth_day.jpg");
    earthNight = jpeg_get(SRCDIR "/../earth_night.jpg");
    if ((earthDay.xsize <= 0) || (earthNight.xsize <= 0)) {
        ephem_fatal(__FILE__, __LINE__, "Could not open Earth images.");
        exit(1);
    }

    // Array for storing the greatest fraction of shadow at every point on flat map
    shadow_map *greatest_shadow = allocate_shadow_map(config.x_size_2d, config.y_size_2d);
    time_span timeSpan = {0, 0, 0, 0, 0, 0, 0};

    // Work out how many ephemeris steps between each video frame
    const int video_ephemeris_stride = (int) round(config.video_time_resolution / config.time_resolution);

    // Select the middle frame of the video to use as a poster image
    config.poster_image_frame = ephemeris->point_count / video_ephemeris_stride / 2;

    // Work out the path of greatest eclipse across the world
    eclipse_path_list *paths = map_greatest_eclipse(&config, ephemeris);
    // exit(0);

    // Loop over video frames
    config.frame_counter = -1;
    for (j = 0; j < ephemeris->point_count; j += video_ephemeris_stride) {
        config.frame_counter++;

        // Read ephemeris
        const double JD = ephemeris->jd_start + ephemeris->jd_step * j;
        const double *pos_sun = ephemeris->data[j].sun_pos;
        const double *pos_earth = ephemeris->data[j].earth_pos;
        const double *pos_moon = ephemeris->data[j].moon_pos;

        shadow_map *shadow_map_2d = calculate_eclipse_map_2d(&config, JD, pos_sun, pos_earth, pos_moon, &timeSpan,
                                                             greatest_shadow);
        render_2d_eclipse_map(&config, JD, earthDay, earthNight, shadow_map_2d, paths);
        update_binary_map(&config, binary_eclipse_maps, shadow_map_2d);
        shadow_map_free(shadow_map_2d);

        shadow_map *shadow_map_3d = calculate_eclipse_map_3d(&config, JD, pos_sun, pos_earth, pos_moon);
        render_3d_eclipse_map(&config, JD, earthDay, shadow_map_3d, paths);
        shadow_map_free(shadow_map_3d);
    }

    // Make an image mapping the greatest magnitude of the eclipse across a 2D map of the world
    render_2d_maximum_extent(&config, earthDay, greatest_shadow);

    // Make a small teaser image showing where the eclipse is visible
    render_2d_eclipse_icon(cl, &config, greatest_shadow);

    // Compute contours of constant eclipse magnitude across the world
    map_eclipse_contours(&config, ephemeris, paths, &timeSpan);

    // Work out the maximum extent of the eclipse by country
    country_lookup_max_eclipse(cl, &config, greatest_shadow);

    // Output a binary file containing a series of snapshots of the eclipse's magnitude across the world
    output_binary_map(&config, ephemeris, binary_eclipse_maps);

    // Return times when the eclipse begins and ends, in partial and total phases
    // All times written as Julian day numbers, in UT
    const double delta_t_days = delta_t(timeSpan.partial_start) / 86400;
    printf("%.10f %.10f %.10f %.10f\n",
           timeSpan.partial_start - delta_t_days,
           timeSpan.partial_end - delta_t_days,
           timeSpan.total_start - delta_t_days,
           timeSpan.total_end - delta_t_days
           );

    // Clean up and exit
    if (DEBUG) ephem_log("Freeing data structures.");
    shadow_map_free(greatest_shadow);
    country_lookup_free(cl);
    jpeg_dealloc(&earthDay);
    jpeg_dealloc(&earthNight);
    if (DEBUG) ephem_log("Terminating normally.");
    return 0;
}
