// map_eclipse_contours.h

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

#ifndef MAP_GREATEST_ECLIPSE_H
#define MAP_GREATEST_ECLIPSE_H 1

#include "ephemeris.h"
#include "map_greatest_eclipse.h"
#include "settings.h"

typedef struct contour_line {
    double longitude[50000]; // radians
    double latitude[50000]; // radians
    int point_count;
} contour_line;

typedef struct contour_line_list {
    contour_line line[10];
    int contour_count;
} contour_line_list;

contour_line_list *map_eclipse_contours(const settings *config, const ephemeris *ephemeris,
                                        const eclipse_path_list *paths, const time_span *eclipse_time_span);

#endif
