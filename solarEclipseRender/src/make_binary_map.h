// make_binary_map.h

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

#ifndef MAKE_BINARY_MAP_H
#define MAKE_BINARY_MAP_H 1

#include "ephemeris.h"
#include "settings.h"
#include "shadow_calc.h"

void update_binary_map(const settings *config, unsigned char *eclipse_maps, const shadow_map *shadow_map);

void output_binary_map(const settings *config, const ephemeris *ephemeris, unsigned char *eclipse_maps);

#endif
