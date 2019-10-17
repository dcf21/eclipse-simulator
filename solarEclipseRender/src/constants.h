// constants.h

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

#ifndef CONSTANTS_H
#define CONSTANTS_H 1

/**
 * Physical constants
 */

#define MOON_FUDGE_FACTOR 1.0    /* facility here to increase the size of the Moon to account for atmospheric effects */
#define RADIUS_MOON   1737.10e3 * MOON_FUDGE_FACTOR /* metres */
#define RADIUS_SUN    695700e3    /* metres */
#define AU            149597871e3 /* metres */

#define RADIUS_EARTH_EQUATOR 6378.1e3 /* metres */
#define RADIUS_EARTH_POLE    6356.8e3 /* metres */
#define RADIUS_EARTH_MEAN    6371.0e3 /* metres */

#endif
