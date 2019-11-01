// settings.h

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

#ifndef SETTINGS_H
#define SETTINGS_H 1

#include "coreUtils/strConstants.h"

typedef struct settings {

    /**
     * output - The directory to store output in
     */
    char *output_dir;

    /**
     * title - The title of this eclipse event, e.g. 'Annular Solar Eclipse'. Useful for debugging...
     */

    char *title;

    /**
     * Dimensions of the video showing the eclipse on a flat map of the world
     */
    int x_size_2d, y_size_2d;

    /**
     * Dimensions of the video showing the eclipse on a 3D globe of the world
     */
    int x_size_3d, y_size_3d;

    /**
     * Radius of the Earth in the 3D video
     */
    int earth_pixel_radius;

    /**
     * Dimensions of eclipse teaser image, used on the page <https://in-the-sky.org/eclipses.php>
     */
    int x_size_teaser, y_size_teaser;

    /**
     * Shadow of Moon lets 50% of Earth light remain in videos
     */
    double moon_shadow_fade_fraction;


    /**
     * Color to use to shade the Moon's shadow on the Earth
     */
    int shadow_col_r;
    int shadow_col_g;
    int shadow_col_b;

    /**
     * Color to use to shade the small dot where a total eclipse is visible
     */
    int totality_col_r;
    int totality_col_g;
    int totality_col_b;

    /**
     * Time span for eclipse simulation, specified in TT
     */
    double jd_min, jd_max;

    /**
     * Time resolution for searching for maximum eclipse at any given location
     */
    double time_resolution;

    /**
     * Time resolution for searching path a maximum eclipse
     */
    double eclipse_path_search_time_resolution;

    /**
     * Time resolution for video frames
     */
    double video_time_resolution;

    /**
     * Time resolution for binary time-lapse snapshots of eclipse
     */

    double binary_map_time_resolution;

    /**
     * Within the binary time-lapse snapshots, sample this number of points per degree of latitude and longitude
     */
    double binary_map_angular_resolution;

    /**
     * month_names - Names of the months, as they should be displayed in the date/time display in the bottom corner
     */
    const char **month_names;

    /**
     * font_name - The font to use to label the eclipse diagrams
     */
    char font_name[FNAME_LENGTH];

    /**
     * poster_image_frame - The frame number to use as a poster image for the eclipse videos
     */
    int poster_image_frame;

    /**
     * frame_counter - The frame currently being worked on
     */
    int frame_counter;
} settings;

/**
 * time_span - Structure for storing the time span of a solar eclipse.
 * All times stored as Julian day numbers, in TT.
 */
typedef struct time_span {
    double total_start;
    double total_end;
    double partial_start;
    double partial_end;

    double greatest_eclipse_magnitude;  // In range 0-1
    double greatest_eclipse_longitude;  // degrees
    double greatest_eclipse_latitude;  // degrees
} time_span;

settings default_settings();

#endif