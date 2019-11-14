// rendering.c

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

#include <stdlib.h>
#include <stdint.h>

#include <gsl/gsl_math.h>

#include <cairo/cairo.h>

#include "rendering.h"

#define X_FIX(x) (x + x_size) % x_size

uint32_t color_blend(const colour *new_colour, const uint32_t *old_colour, double fraction) {
    // Cairo's ARGB32 pixel format stores pixels as 32-bit ints, with alpha in most significant byte.
    unsigned int red_1 = (unsigned) new_colour->red;
    unsigned int green_1 = (unsigned) new_colour->grn;
    unsigned int blue_1 = (unsigned) new_colour->blu;

    unsigned int red_2 = (*old_colour >> (unsigned) 16) & (unsigned) 255;
    unsigned int green_2 = (*old_colour >> (unsigned) 8) & (unsigned) 255;
    unsigned int blue_2 = (*old_colour) & (unsigned) 255;

    unsigned int red = (unsigned) (red_1 * fraction + red_2 * (1 - fraction));
    unsigned int green = (unsigned) (green_1 * fraction + green_2 * (1 - fraction));
    unsigned int blue = (unsigned) (blue_1 * fraction + blue_2 * (1 - fraction));

    return ((uint32_t) blue +  // blue
            ((uint32_t) green << (unsigned) 8) +  // green
            ((uint32_t) red << (unsigned) 16) + // red
            ((uint32_t) 255 << (unsigned) 24)  // alpha
    );

}

void set_pixel(unsigned char *frame, int x_size, int stride, int x, int y, const colour *colour, double fraction) {
    const int x_pos = X_FIX(x);
    const int y_pos = y;
    const int offset = x_pos * 4 + y_pos * stride;
    uint32_t blended_color = color_blend(colour, (const uint32_t *) &frame[offset], fraction);
    *(uint32_t *) &frame[offset] = blended_color;
}

int test_pixel(const shadow_map *shadow, int x, int y, int x_size, double level) {
    return ((shadow->map[y * x_size + X_FIX(x)] > level) &&
            ((shadow->map[(y - 1) * x_size + X_FIX(x)] <= level) ||
             (shadow->map[(y - 1) * x_size + X_FIX(x - 1)] <= level) ||
             (shadow->map[(y) * x_size + X_FIX(x - 1)] <= level) ||
             (shadow->map[(y + 1) * x_size + X_FIX(x - 1)] <= level) ||
             (shadow->map[(y + 1) * x_size + X_FIX(x)] <= level) ||
             (shadow->map[(y + 1) * x_size + X_FIX(x + 1)] <= level) ||
             (shadow->map[(y) * x_size + X_FIX(x + 1)] <= level) ||
             (shadow->map[(y - 1) * x_size + X_FIX(x + 1)] <= level)
            )
    );
}

void shadowContoursLabelPositions(const double *contourList, const shadow_map *shadow,
                                  int x_offset, int x_size, int y_size,
                                  int **label_position_x, int **label_position_y,
                                  int *previous_label_position_x, int *previous_label_position_y) {
    int y, x_output, i, j;

    *label_position_x = malloc(256 * sizeof(int));
    *label_position_y = malloc(256 * sizeof(int));

    for (i = 0; i < 256; i++) {
        (*label_position_x)[i] = (*label_position_y)[i] = -10000;
    }

    for (i = 0; contourList[i] >= 0; i++) {
        const double level = contourList[i] / 100.;
        double best_radius = 1e6;

        // Do not place contour labels right at the edges of the image
        for (y = 25; y < y_size - 25; y++)
            for (x_output = 25; x_output < x_size - 25; x_output++) {
                const int x = X_FIX(x_output + x_offset);
                if (test_pixel(shadow, x, y, x_size, level)) {

                    double radius = sqrt(gsl_pow_2(x_output - x_size / 2.) + gsl_pow_2(y - y_size / 2.));

                    if (previous_label_position_x != NULL) {
                        for (j = 0; contourList[j] >= 0; j++) {
                            if (previous_label_position_x[j] > 0) {
                                radius += sqrt(gsl_pow_2(x_output - previous_label_position_x[j]) +
                                               gsl_pow_2(y - previous_label_position_y[j]));
                            }
                        }
                    }

                    if (radius < best_radius) {
                        best_radius = radius;
                        (*label_position_x)[i] = x_output;
                        (*label_position_y)[i] = y;
                    }
                }
            }
    }
}

/**
 * drawShadowContours - Draw contours onto a GD image surface, wherever the values in the array <shadow> pass any
 * of the thresholds in the array <contourList>.
 * @param frame - The image surface onto which to trace the contours.
 * @param shadow - The array of shadow fractions within each pixel, which we are to draw contours from
 * @param x_offset - Shift the output diagram horizontally by some number of pixels to place a longitude other than
 *                   zero at the centre.
 * @param stride - The number of bytes separating consecutive rows in <frame>.
 * @param bold - Boolean flag indicating whether to make contours bold
 * @param x_size - The horizontal pixel size of the output image
 * @param y_size - The vertical pixel size of the output image
 */
void drawShadowContours(unsigned char *frame, const double *contourList, const shadow_map *shadow,
                        int *label_position_x, int *label_position_y,
                        int x_offset, int stride, int x_size, int y_size) {
    int x, y, i, j;

    const colour yellow = {255, 255, 0};
    const colour red = {255, 0, 0};

    for (i = 0; contourList[i] >= 0; i++) {
        double level = contourList[i] / 100.;

        const colour color = (level > 0.01) ? yellow : red;
        for (y = 1; y < y_size - 1; y++)
            for (x = 0; x < x_size; x++) {
                if (test_pixel(shadow, x, y, x_size, level)) {
                    const int x_output = X_FIX(x - x_offset);

                    int mask_pixel = 0;
                    for (j = 0; contourList[j] >= 0; j++) {
                        const double radius = (gsl_pow_2(x_output - label_position_x[j]) +
                                               gsl_pow_2(y - label_position_y[j]));
                        const double critical_radius = 16;
                        if (radius < critical_radius * critical_radius) {
                            mask_pixel = 1;
                        }
                    }
                    if (mask_pixel) continue;

                    set_pixel(frame, x_size, stride, x_output, y, &color, 1);

                    set_pixel(frame, x_size, stride, x_output - 1, y, &color, 0.25);
                    set_pixel(frame, x_size, stride, x_output + 1, y, &color, 0.25);

                    set_pixel(frame, x_size, stride, x_output, y - 1, &color, 0.25);
                    set_pixel(frame, x_size, stride, x_output, y + 1, &color, 0.25);

                    set_pixel(frame, x_size, stride, x_output - 1, y - 1, &color, 0.1);
                    set_pixel(frame, x_size, stride, x_output + 1, y - 1, &color, 0.1);
                    set_pixel(frame, x_size, stride, x_output - 1, y + 1, &color, 0.1);
                    set_pixel(frame, x_size, stride, x_output + 1, y + 1, &color, 0.1);
                }
            }
    }
}

//! chart_label - Write a text label onto a cairo page immediately.
//! \param cairo_draw - The cairo drawing context
//! \param colour - The colour to use to write the text
//! \param label - The string of the text label
//! \param x_canvas - The horizontal position of the label
//! \param y_canvas - The vertical position of the label
//! \param h_align - The horizontal alignment of the label
//! \param v_align - The vertical alignment of the label
//! \param font_size - Font size of label (pixels)
//! \param font_bold - Boolean indicating whether to render label in bold
//! \param font_italic - Boolean indicating whether to render label in italic
//! \return - None

void chart_label(cairo_t *cairo_draw, colour colour, const char *label,
                 double x_canvas, double y_canvas, int h_align, int v_align,
                 double font_size, int font_bold, int font_italic) {

    // Select font
    cairo_text_extents_t extents;
    cairo_set_font_size(cairo_draw, font_size);
    cairo_select_font_face(cairo_draw, "Sans",
                           font_italic ? CAIRO_FONT_SLANT_ITALIC : CAIRO_FONT_SLANT_NORMAL,
                           font_bold ? CAIRO_FONT_WEIGHT_BOLD : CAIRO_FONT_WEIGHT_NORMAL);

    // Measure text bounding box
    cairo_text_extents(cairo_draw, label, &extents);

    switch (h_align) {
        case 1:
            x_canvas -= extents.width + extents.x_bearing; // right align
            break;
        case 0:
            x_canvas -= extents.width / 2 + extents.x_bearing; // centre align
            break;
        case -1:
            x_canvas -= extents.x_bearing; // left align
            break;
        default:
            break;
    }

    switch (v_align) {
        case 1:
            y_canvas -= extents.height + extents.y_bearing; // top align
            break;
        case 0:
            y_canvas -= extents.height / 2 + extents.y_bearing; // centre align
            break;
        case -1:
            y_canvas -= extents.y_bearing; // bottom align
            break;
        default:
            break;
    }

    // Render the text label itself
    cairo_set_source_rgb(cairo_draw, colour.red, colour.grn, colour.blu);
    cairo_move_to(cairo_draw, x_canvas, y_canvas);
    cairo_show_text(cairo_draw, label);
}
