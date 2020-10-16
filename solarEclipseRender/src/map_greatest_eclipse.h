// map_greatest_eclipse.h

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

#ifndef MAP_ECLIPSE_CONTOURS_H
#define MAP_ECLIPSE_CONTOURS_H 1

#include "ephemeris.h"
#include "settings.h"

// Data structure for storing the path of greatest eclipse

#define MAX_PATH_LENGTH 2500
#define MAX_PATH_ITEMS 20

typedef struct path_point {
    double latitude, longitude; // radians
    double jd; // TT
    double duration; // seconds
} path_point;

typedef struct eclipse_path {
    int is_total;
    double jd_start, jd_end;
    int point_count;
    path_point path[MAX_PATH_LENGTH];
} eclipse_path;

typedef struct eclipse_path_list {
    eclipse_path paths[MAX_PATH_ITEMS];
    int path_count;
    double latitude_midpoint, longitude_midpoint;
    double maximum_duration;
} eclipse_path_list;

double eclipse_duration_from_path(const eclipse_path_list *paths, double jd, int *is_total);

eclipse_path_list *map_greatest_eclipse(const settings *config, const ephemeris *ephemeris, time_span *span_output);

void eclipse_position_from_path(const eclipse_path_list *paths, double jd, double *lng_out, double *lat_out);

#endif
