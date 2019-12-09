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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <cairo/cairo.h>

#include <gsl/gsl_math.h>

#include "jpeg/jpeg.h"

#include "coreUtils/errorReport.h"
#include "coreUtils/strConstants.h"

#include "mathsTools/julianDate.h"

#include "map_greatest_eclipse.h"
#include "projection.h"
#include "rendering.h"
#include "render_3d.h"
#include "settings.h"
#include "shadow_calc.h"

/**
 * render_3d_eclipse_map - Render a 3D snapshot of a globe of the world as viewed from the direction of the Sun, with
 * the magnitude of a solar eclipse overlaid as contours. This is ideal for turning into an animation of the eclipse's
 * process if snapshots are turned into a video.
 *
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param jd [in] - The Julian day number of the current point in the simulation (TT)
 * @param earthDay [in] - A JPEG image of the world in daylight
 * @param pos_sun [in] - The position of the Sun, as returned by ephemerisCompute, in AU, relative to solar system barycentre
 * @param pos_earth [in] - The position of the centre of the Earth, as returned by ephemerisCompute
 * @param shadow_map [in] - A binary map of the eclipse magnitude across the world
 * @param eclipse_path [in] - The path of greatest eclipse, with duration at each point
 */
void render_3d_eclipse_map(settings *config, double jd, jpeg_ptr earthDay,
                           const double *pos_sun, const double *pos_earth,
                           const shadow_map *shadow_map, const eclipse_path_list *eclipse_path) {
    int x, y;

    // Generate image RGB data
    const int stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, config->x_size_3d);

    // Allocate buffer for image data
    unsigned char *pixel_data = malloc(stride * config->y_size_3d);

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_3d; y++)
        for (x = 0; x < config->x_size_3d; x++) {
            // Pixel's offset inside shadow_map
            const int offset = x + y * config->x_size_3d;

            jpeg_ptr *srcimg = &earthDay;

            // Look up eclipse magnitude in this pixel
            double shadow = shadow_map->map[offset];

            colour colour_this = {0, 0, 0};

            if (!gsl_isnan(shadow_map->lng[offset])) {
                // Lat & lng of this pixel are finite, which means pixel lies on surface of Earth
                // Otherwise this pixel lies off the side of the globe of the world
                int p0 = (int) ((shadow_map->lng[offset] + 180.) * srcimg->xsize / 360.);
                int p1 = (int) ((90. - shadow_map->lat[offset]) * srcimg->ysize / 180.);
                if (p0 >= srcimg->xsize) p0 -= srcimg->xsize;
                if (p1 >= srcimg->ysize) p1 = srcimg->ysize - 1;

                const colour colour_earth = {
                        ((int) srcimg->data_red[p0 + p1 * srcimg->xsize]),
                        ((int) srcimg->data_grn[p0 + p1 * srcimg->xsize]),
                        ((int) srcimg->data_blu[p0 + p1 * srcimg->xsize])
                };

                colour_this = colour_earth;
            }

            // Superimpose shadow map over Earth
            if (shadow > 0) {
                // If this pixel experiences a partial eclipse, shade it accordingly
                colour_this.red = (int) (config->moon_shadow_fade_fraction * colour_this.red +
                                         (1 - config->moon_shadow_fade_fraction) * config->shadow_col_r);
                colour_this.grn = (int) (config->moon_shadow_fade_fraction * colour_this.grn +
                                         (1 - config->moon_shadow_fade_fraction) * config->shadow_col_g);
                colour_this.blu = (int) (config->moon_shadow_fade_fraction * colour_this.blu +
                                         (1 - config->moon_shadow_fade_fraction) * config->shadow_col_b);
            }

            // Set pixel color
            const int output_offset = x * 4 + y * stride;

            *(uint32_t *) &pixel_data[output_offset] = ((uint32_t) colour_this.blu +  // blue
                                                        ((uint32_t) colour_this.grn << (unsigned) 8) +  // green
                                                        ((uint32_t) colour_this.red << (unsigned) 16) + // red
                                                        ((uint32_t) 255 << (unsigned) 24)  // alpha
            );
        }

    // Overlay contours of eclipse magnitude on top of the map
    const double contourList[] = {80, 60, 40, 20, 1e-6, -1};
    int *label_position_x, *label_position_y;
    static int previous_label_position_x[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    static int previous_label_position_y[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    shadowContoursLabelPositions(contourList, shadow_map, 0,
                                 config->x_size_3d, config->y_size_3d,
                                 &label_position_x, &label_position_y,
                                 previous_label_position_x, previous_label_position_y);
    drawShadowContours(pixel_data, contourList, shadow_map, label_position_x, label_position_y,
                       0, 0, stride, config->x_size_3d, config->y_size_3d);

    // Turn bitmap data into a Cairo surface
    cairo_surface_t *surface = cairo_image_surface_create_for_data(pixel_data, CAIRO_FORMAT_ARGB32,
                                                                   config->x_size_3d, config->y_size_3d,
                                                                   stride);

    cairo_t *cairo_draw = cairo_create(surface);

    // Label contours
    for (int i = 0; contourList[i] >= 0; i++)
        if (label_position_x[i] >= 0) {
            const double level = contourList[i];
            char text[8];
            sprintf(text, "%.0f%%", level);

            const colour yellow = {255, 255, 0};
            const colour red = {255, 0, 0};
            const colour colour_this = (level > 0.01) ? yellow : red;

            chart_label(cairo_draw, colour_this, text, label_position_x[i], label_position_y[i],
                        0, 0, 13, 1, 0);
        }
    memcpy(previous_label_position_x, label_position_x, sizeof(previous_label_position_x));
    memcpy(previous_label_position_y, label_position_y, sizeof(previous_label_position_y));
    free(label_position_x);
    free(label_position_y);

    // Mark the position of central eclipse
    double lng_central, lat_central;
    eclipse_position_from_path(eclipse_path, jd, &lng_central, &lat_central);

    if (gsl_finite(lng_central)) {
        const int cross_size = 4;

        // Calculate the latitude and longitude on Earth where the Sun is overhead
        double sidereal_time, lat_sun, lng_sun;
        calculate_where_sun_overhead(&lat_sun, &lng_sun, &sidereal_time, pos_sun, pos_earth, jd);

        int x_centre, y_centre;
        inv_project_3d(config, &x_centre, &y_centre, lng_sun, lat_sun, lng_central, lat_central);

        cairo_set_source_rgb(cairo_draw, 0, 1, 0);
        cairo_new_path(cairo_draw);
        cairo_move_to(cairo_draw, x_centre - cross_size, y_centre - cross_size);
        cairo_line_to(cairo_draw, x_centre + cross_size, y_centre + cross_size);
        cairo_move_to(cairo_draw, x_centre - cross_size, y_centre + cross_size);
        cairo_line_to(cairo_draw, x_centre + cross_size, y_centre - cross_size);
        cairo_stroke(cairo_draw);
    }

    // Look up duration of the eclipse
    int is_total;
    const double duration = eclipse_duration_from_path(eclipse_path, jd, &is_total);

    // Get date components (in UT; not TT)
    int year, month, day, hour, min, status = 0;
    double sec;
    inv_julian_day(jd - delta_t(jd) / 86400.,
                   &year, &month, &day, &hour, &min, &sec, &status, temp_err_string);

    // Write the time and date in bottom left corner of the image
    char text[FNAME_LENGTH];
    colour yellow = {255, 255, 0};

    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.6);
    cairo_rectangle(cairo_draw,
                    0, config->y_size_3d - 70,
                    190, 70);
    cairo_fill(cairo_draw);

    sprintf(text, "%d %s %d UTC", day, config->month_names[month], year);
    chart_label(cairo_draw, yellow, text, 95, config->y_size_3d - 17, 0, 0, 17, 1, 0);

    sprintf(text, "%02d:%02d", hour, min);
    chart_label(cairo_draw, yellow, text, 95, config->y_size_3d - 49, 0, 0, 28, 1, 0);

    // Write copyright text in bottom right corner of the image
    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.6);
    cairo_rectangle(cairo_draw,
                    config->x_size_3d - 240, config->y_size_3d - 66,
                    240, 66);
    cairo_fill(cairo_draw);

    sprintf(text, "\u00A9 Dominic Ford 2012\u20132019");
    chart_label(cairo_draw, yellow, text, config->x_size_3d - 12, config->y_size_3d - 44, 1, 0, 14, 1, 0);

    sprintf(text, "https://in-the-sky.org/");
    chart_label(cairo_draw, yellow, text, config->x_size_3d - 12, config->y_size_3d - 18, 1, 0, 14, 1, 0);

    // Write duration in the top left corner of the image
    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.6);
    cairo_rectangle(cairo_draw,
                    0, 0,
                    166, 70);
    cairo_fill(cairo_draw);

    sprintf(text, "Central duration");
    chart_label(cairo_draw, yellow, text, 83, 20, 0, 0, 15, 1, 0);

    if (duration > 0) {
        sprintf(text, is_total ? "Total" : "Annular");
        chart_label(cairo_draw, yellow, text, 45, 50, 0, 0, 15, 1, 0);

        sprintf(text, "%dm%02ds", (int) (duration / 60), (int) duration % 60);
        chart_label(cairo_draw, yellow, text, 125, 50, 0, 0, 15, 1, 0);
    } else {
        sprintf(text, "\u2014");
        chart_label(cairo_draw, yellow, text, 83, 50, 0, 0, 16, 1, 0);
    }

    // Write output image
    char output_filename[FNAME_LENGTH];
    sprintf(output_filename, "%s/frameB%06d.png", config->output_dir, config->frame_counter);
    cairo_surface_write_to_png(surface, output_filename);

    // If this frame is at the midpoint of the eclipse we output a special "poster" frame which acts as a teaser image
    // for the animation.
    if (config->frame_counter == config->poster_image_frame) {
        sprintf(text, "Please wait \u2013 loading...");
        chart_label(cairo_draw, yellow, text, (int) (config->x_size_3d / 2), (int) (config->y_size_3d / 2),
                    0, 0, 16, 1, 0);

        // Write out poster image
        sprintf(output_filename, "%s/solarEclipseB.png", config->output_dir);
        cairo_surface_write_to_png(surface, output_filename);
    }

    cairo_destroy(cairo_draw);
    cairo_surface_finish(surface);
    free(pixel_data);
}
