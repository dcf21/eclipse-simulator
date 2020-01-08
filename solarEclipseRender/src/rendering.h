// rendering.h

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

#ifndef RENDERING_H
#define RENDERING_H 1

#include <cairo/cairo.h>
#include "shadow_calc.h"

//! Define an RGB colour to use to draw a particular item on a chart
typedef struct colour {
    double red, grn, blu;
} colour;

void shadowContoursLabelPositions(const double *contourList, const shadow_map *shadow,
                                  int x_offset, int x_size, int y_size,
                                  int **label_position_x, int **label_position_y,
                                  int *previous_label_position_x, int *previous_label_position_y);

void drawShadowContours(unsigned char *frame, const double *contourList, const shadow_map *shadow,
                        int *label_position_x, int *label_position_y,
                        int x_offset, int x_wrap, int stride, int x_size, int y_size);

void chart_label(cairo_t *cairo_draw, colour colour, const char *label,
                 double x_canvas, double y_canvas, int h_align, int v_align,
                 double font_size, int font_bold, int font_italic);

#endif