// settings.c

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

#include <string.h>

#include "settings.h"

settings default_settings() {
    settings output;

    output.output_dir = "/tmp/eclipse_demo/";
    output.title = "Undefined";

    output.x_size_2d = 1600;
    output.y_size_2d = 800;

    output.x_size_3d = 800;
    output.y_size_3d = 800;

    output.earth_pixel_radius = 320;

    output.x_size_teaser = 400;
    output.y_size_teaser = 200;

    output.moon_shadow_fade_fraction = 0.5;

    output.shadow_col_r = 16;
    output.shadow_col_g = 16;
    output.shadow_col_b = 16;

    output.jd_min = 2449483.217363426;
    output.jd_max = 2449483.217365426;

    output.time_resolution = 0.1;

    output.eclipse_path_search_time_resolution = 10;
    output.video_time_resolution = 100;
    output.binary_map_time_resolution = 600;
    output.binary_map_angular_resolution = 4;

    static const char *monthNames[] = {"x", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov",
                                       "Dec"};
    output.month_names = monthNames;

    return output;
}
