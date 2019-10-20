// make_binary_map.c

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

#include "coreUtils/strConstants.h"

#include "ephemeris.h"
#include "make_binary_map.h"
#include "settings.h"
#include "shadow_calc.h"

// We make a binary file containing time-lapse maps of the eclipse magnitude across the world. This is used to show
// visitors the progress of the eclipse at a series of time points at their location. These binary maps are at a
// lower time and angular resolution than the videos.

/**
 * update_binary_map - This is called for each video frame we generate, to see if we should store a binary snapshot
 * of the eclipse's magnitude across the world at this point in time.
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param eclipse_maps [out] - The data structure in which we store the binary snapshots
 * @param shadow_map [in] - A binary map of the eclipse magnitude across the world
 */
void update_binary_map(const settings *config, unsigned char *eclipse_maps, const shadow_map *shadow_map) {
    int x, y;

    // Calculate how many video frames lie between subsequent binary snapshots
    const int binary_map_video_frame_stride = (int) round(config->binary_map_time_resolution /
                                                          config->video_time_resolution);

    // Do not proceed if this frame is not to be made into a binary snapshot
    if ((config->frame_counter % binary_map_video_frame_stride) != 0) return;

    // Update eclipse profile data structure

    // Number of pixels along the longitude direction of each binary snapshot
    const int stride0 = (int) (360 * config->binary_map_angular_resolution);

    // Number of pixels along the latitude direction of each binary snapshot
    const int stride1 = (int) (180 * config->binary_map_angular_resolution);

    // The frame number of this binary snapshot
    const int k = config->frame_counter / binary_map_video_frame_stride;

    // Loop over all the pixels in this binary snapshot
    for (y = 0; y < stride1; y++)
        for (x = 0; x < stride0; x++) {
            // Work out offset of this pixel in the output structure
            const int p0 = (int) (((double) x) / stride0 * config->x_size_2d);
            const int p1 = (int) (((double) y) / stride1 * config->y_size_2d);

            // Set pixel value
            const int offset = p0 + p1 * config->x_size_2d;
            const double shadow_fraction = gsl_max(0, shadow_map->map[offset]);  // may be -1 when in Earth shadow
            eclipse_maps[x + stride0 * (y + stride1 * k)] = (int) (100 * shadow_fraction);
        }
}

/**
 * output_binary_map - Having finished compiling an array of binary snapshots, now write it to disk.
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param ephemeris [in] - The ephemeris used to construct this eclipse simulation
 * @param eclipse_maps [in] - The data structure in which we store the binary snapshots
 */
void output_binary_map(const settings *config, const ephemeris *ephemeris, unsigned char *eclipse_maps) {
    int x, y, j;

    // Output eclipse profile binary file
    // Binary file commences with array of [lat][lng] with file position pointers for where data starts for each position on the globe
    // At this file position is a null-terminated structure of the form [int] [char] [char] [char] ...
    // double is the unix time at which data begins. chars are eclipse extent (percent) at 10 minute intervals
    // when extent reaches zero, the track stops

    // Number of pixels along the longitude direction of each binary snapshot
    const int stride0 = (int) (360 * config->binary_map_angular_resolution);

    // Number of pixels along the latitude direction of each binary snapshot
    const int stride1 = (int) (180 * config->binary_map_angular_resolution);

    // How many binary snapshots have we created?
    const int binary_snapshot_count = ceil((config->jd_max - config->jd_min) /
                                           (config->binary_map_time_resolution / 86400));

    // Reprocess the data we have collected into a compacted data structure
    unsigned char *binary_eclipse_maps_2 = calloc(stride0 * stride1 * (binary_snapshot_count + 4), 1);

    // How many samples are there in the input ephemeris between successive binary snapshots
    const int binary_map_ephemeris_stride = (int) round(config->binary_map_time_resolution /
                                                        config->time_resolution);

    const int jmax = ephemeris->point_count / binary_map_ephemeris_stride;
    int *ptrArray = (int *) binary_eclipse_maps_2;
    int ptr = stride0 * stride1 * sizeof(int);
    for (y = 0; y < stride1; y++)
        for (x = 0; x < stride0; x++) {
            int state = 0;
            for (j = 0; j < jmax; j++) {
                int c = x + stride0 * (y + stride1 * j);
                unsigned char f = eclipse_maps[c];
                if ((!state) && f) {
                    ptrArray[x + y * stride0] = ptr;

                    const int offset = j * binary_map_ephemeris_stride;
                    const double JD = ephemeris->jd_start + ephemeris->jd_step * offset;
                    const double utc = 86400.0 * (JD - 2440587.5);

                    *((double *) (binary_eclipse_maps_2 + ptr)) = utc;
                    ptr += sizeof(double);
                    state = 1;
                } else if (!state) continue;
                binary_eclipse_maps_2[ptr++] = f;
                if (state && (!f)) { break; } // End of eclipse has been reached
            }
        }

    // Write out binary snapshot data
    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/maximumEclipse.dat", config->output_dir);
    FILE *f = fopen(fname, "wb");
    fwrite((void *) binary_eclipse_maps_2, 1, ptr, f);
    fclose(f);
}