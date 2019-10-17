// ephemeris.c

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

#ifndef EPHEMERIS_H
#define EPHEMERIS_H 1

#include "settings.h"

//! A structure defining a point in an ephemeris
typedef struct ephemeris_point {
    double sun_pos[3];
    double earth_pos[3];
    double moon_pos[3];
} ephemeris_point;

//! A structure defining an ephemeris of a solar system object
typedef struct ephemeris {
    double jd_start, jd_end, jd_step; // Julian day numbers
    int point_count;
    ephemeris_point *data;
} ephemeris;


void fetch_ephemeris(const settings *config, ephemeris **output);

#endif

