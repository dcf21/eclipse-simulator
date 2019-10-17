// render_3d.c

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

#ifndef RENDER_3D_H
#define RENDER_3D_H 1

#include "map_greatest_eclipse.h"
#include "settings.h"
#include "shadow_calc.h"

void render_3d_eclipse_map(settings *config, double jd, jpeg_ptr earthDay,
                           const shadow_map *shadow_map, const eclipse_path_list *eclipse_path);

#endif
