// ephemeris.c

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

#include <stdio.h>
#include <stdlib.h>

#include "coreUtils/strConstants.h"
#include "coreUtils/errorReport.h"

#include "ephemeris.h"
#include "settings.h"

void fetch_ephemeris(const settings *config, ephemeris **output) {

    // Path to ephemeris calculator
    char ephemeris_compute_path[FNAME_LENGTH];

    sprintf(ephemeris_compute_path,
            "%s/../../ephemeris-compute-de430/bin/serial/ephem.bin",
            SRCDIR);

    // Compute data points at this interval (measured in days)
    const double jd_step = config->time_resolution / 86400.;

    // Data structure to store ephemeris for Sun, Earth and Moon
    ephemeris *x = (ephemeris *) malloc(sizeof(ephemeris));

    // Generous estimate of how many lines we expect ephemerisCompute to return
    const int point_count_max = (int) (2 + (config->jd_max - config->jd_min) / jd_step);

    // Allocate data to hold the ephemeris
    x->data = (ephemeris_point *) malloc(point_count_max * sizeof(ephemeris_point));

    // Use ephemerisCompute to track the path of this object
    char ephemeris_compute_command[FNAME_LENGTH];

    snprintf(ephemeris_compute_command, FNAME_LENGTH, "%.2048s "
                                                      "--jd_min %.14f "
                                                      "--jd_max %.14f "
                                                      "--jd_step %.14f "
                                                      "--output_format 0 "
                                                      "--output_constellations 0 "
                                                      "--output_binary 1 "
                                                      "--objects \"sun,moon,earth\" ",
             ephemeris_compute_path,
             config->jd_min, config->jd_max, jd_step);

    FILE *ephemeris_data = popen(ephemeris_compute_command, "r");
    // printf("%s\n", ephemeris_compute_command);

    // Loop over the lines returned by ephemerisCompute
    int line_counter = 0;
    while ((!feof(ephemeris_data)) && (!ferror(ephemeris_data))) {
        int items_fetched;

        // Read columns of data
        items_fetched = (int)fread(&x->data[line_counter].sun_pos, sizeof(double), 3, ephemeris_data);
        if (items_fetched == 0) break;
        if (items_fetched < 3) logging_fatal(__FILE__, __LINE__, "Too few items");
        items_fetched = (int)fread(&x->data[line_counter].moon_pos, sizeof(double), 3, ephemeris_data);
        if (items_fetched < 3) logging_fatal(__FILE__, __LINE__, "Too few items");
        items_fetched = (int)fread(&x->data[line_counter].earth_pos, sizeof(double), 3, ephemeris_data);
        if (items_fetched < 3) logging_fatal(__FILE__, __LINE__, "Too few items");

        // Increment data point counter
        line_counter++;
    }

    // Throw an error if we got no data
    if (line_counter == 0) {
        logging_fatal(__FILE__, __LINE__, "ephemerisCompute returned no data");
        exit(1);
    }

    // Record how many lines of data were returned
    x->point_count = line_counter;
    x->jd_start = config->jd_min;
    x->jd_end = config->jd_max;
    x->jd_step = jd_step;

    if (DEBUG) {
        sprintf(temp_err_string, "Read list of %d ephemeris points from ephemeris.", line_counter);
        logging_log(temp_err_string);
    }

    // Return output
    *output = x;
}

