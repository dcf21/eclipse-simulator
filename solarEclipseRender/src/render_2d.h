// render_2d.h

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

#ifndef RENDER_2D_H
#define RENDER_2D_H 1

#include "country_lookup.h"
#include "map_greatest_eclipse.h"
#include "settings.h"
#include "shadow_calc.h"
#include "map_eclipse_contours.h"

void render_2d_eclipse_map(settings *config, double jd, jpeg_ptr earthDay, jpeg_ptr earthNight,
                           const shadow_map *shadow_map, const eclipse_path_list *eclipse_path);

void render_2d_maximum_extent(const country_lookup_handle *cl, const settings *config,
                              const contour_line_list *contours,
                              const shadow_map *greatest_shadow,
                              const eclipse_path_list *eclipse_path, const char *format);

void render_2d_eclipse_icon(const country_lookup_handle *cl, const settings *config,
                            const shadow_map *greatest_shadow);

#endif