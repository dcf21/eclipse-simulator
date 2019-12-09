// render_2d.c

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
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <cairo/cairo-pdf.h>

#include <gsl/gsl_math.h>

#include "jpeg/jpeg.h"

#include "coreUtils/errorReport.h"
#include "coreUtils/strConstants.h"

#include "mathsTools/julianDate.h"

#include "country_lookup.h"
#include "map_greatest_eclipse.h"
#include "rendering.h"
#include "render_2d.h"
#include "settings.h"
#include "shadow_calc.h"
#include "map_eclipse_contours.h"

/**
 * render_2d_eclipse_map - Render a 2D flat snapshot map of the magnitude of a solar eclipse across the world. This is
 * ideal for turning into an animation of the eclipse's process if snapshots are turned into a video.
 *
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param jd [in] - The Julian day number of the current point in the simulation (TT)
 * @param earthDay [in] - A JPEG image of the world in daylight
 * @param earthNight [in] - A JPEG image of the world by night
 * @param shadow_map [in] - A binary map of the eclipse magnitude across the world
 * @param eclipse_path [in] - The path of greatest eclipse, with duration at each point
 */
void render_2d_eclipse_map(settings *config, double jd, jpeg_ptr earthDay, jpeg_ptr earthNight,
                           const shadow_map *shadow_map, const eclipse_path_list *eclipse_path) {
    int x, y;

    // Generate image RGB data
    const int stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, config->x_size_2d);

    // Allocate buffer for image data
    unsigned char *pixel_data = malloc(stride * config->y_size_2d);

    // Horizontally shift map of the world to put the point of greatest eclipse at the centre
    int x_offset = (int) (config->solar_longitude_at_midpoint * config->x_size_2d / 360);
    if (x_offset > config->x_size_2d / 2) x_offset -= config->x_size_2d;
    if (x_offset < -config->x_size_2d / 2) x_offset += config->x_size_2d;

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            const int offset = x + y * config->x_size_2d;

            // Work out latitude / longitude of this pixel
            const double latitude = shadow_map->lat[offset];
            double longitude = shadow_map->lng[offset];

            jpeg_ptr *srcimg_day = &earthDay;
            jpeg_ptr *srcimg_night = &earthNight;

            // Fetch colour of day time pixel
            int p0_day = (int) ((longitude + 180.) * srcimg_day->xsize / 360.);
            int p1_day = (int) ((90. - latitude) * srcimg_day->ysize / 180.);
            if (p0_day >= srcimg_day->xsize) p0_day -= srcimg_day->xsize;
            if (p1_day >= srcimg_day->ysize) p1_day = srcimg_day->ysize - 1;

            const colour colour_day = {
                    ((int) srcimg_day->data_red[p0_day + p1_day * srcimg_day->xsize]),
                    ((int) srcimg_day->data_grn[p0_day + p1_day * srcimg_day->xsize]),
                    ((int) srcimg_day->data_blu[p0_day + p1_day * srcimg_day->xsize])
            };

            // Fetch colour of night time pixel
            int p0_night = (int) ((longitude + 180.) * srcimg_night->xsize / 360.);
            int p1_night = (int) ((90. - latitude) * srcimg_night->ysize / 180.);
            if (p0_night >= srcimg_night->xsize) p0_night -= srcimg_night->xsize;
            if (p1_night >= srcimg_night->ysize) p1_night = srcimg_night->ysize - 1;

            const colour colour_night = {
                    ((int) srcimg_night->data_red[p0_night + p1_night * srcimg_night->xsize]),
                    ((int) srcimg_night->data_grn[p0_night + p1_night * srcimg_night->xsize]),
                    ((int) srcimg_night->data_blu[p0_night + p1_night * srcimg_night->xsize])
            };
            // Look up eclipse magnitude in the current pixel
            const double shadow = shadow_map->map[x + y * config->x_size_2d];

            // Test whether this pixel is on the day side of the Earth, or the night side
            const int night_time = (shadow < 0);

            colour colour_this;
            if (night_time) {
                // Night time
                const colour colour_this_night = {
                        (colour_night.red * 0.8 + colour_day.red * 0.2),
                        (colour_night.grn * 0.8 + colour_day.grn * 0.2),
                        (colour_night.blu * 0.8 + colour_day.blu * 0.2)
                };
                colour_this = colour_this_night;
            } else {
                // Day time
                const colour colour_this_day = {
                        (colour_night.red * 0.1 + colour_day.red * 0.9),
                        (colour_night.grn * 0.1 + colour_day.grn * 0.9),
                        (colour_night.blu * 0.1 + colour_day.blu * 0.9)
                };
                colour_this = colour_this_day;
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
            const int output_offset = ((x - x_offset + config->x_size_2d) % config->x_size_2d) * 4 + y * stride;

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
    shadowContoursLabelPositions(contourList, shadow_map, x_offset,
                                 config->x_size_2d, config->y_size_2d,
                                 &label_position_x, &label_position_y,
                                 previous_label_position_x, previous_label_position_y);
    drawShadowContours(pixel_data, contourList, shadow_map, label_position_x, label_position_y,
                       x_offset, 1, stride, config->x_size_2d, config->y_size_2d);

    // Turn bitmap data into a Cairo surface
    cairo_surface_t *surface = cairo_image_surface_create_for_data(pixel_data, CAIRO_FORMAT_ARGB32,
                                                                   config->x_size_2d, config->y_size_2d,
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
        const int x_centre = ((int) ((lng_central + 180) / 360. * config->x_size_2d - x_offset + config->x_size_2d) %
                              config->x_size_2d);
        const int y_centre = (int) ((90 - lat_central) / 180 * config->y_size_2d);

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

    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.5);
    cairo_rectangle(cairo_draw,
                    0, config->y_size_2d - 70,
                    190, 70);
    cairo_fill(cairo_draw);

    sprintf(text, "%d %s %d UTC", day, config->month_names[month], year);
    chart_label(cairo_draw, yellow, text, 95, config->y_size_2d - 17, 0, 0, 17, 1, 0);

    sprintf(text, "%02d:%02d", hour, min);
    chart_label(cairo_draw, yellow, text, 95, config->y_size_2d - 49, 0, 0, 28, 1, 0);

    // Write copyright text in bottom right corner of the image
    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.5);
    cairo_rectangle(cairo_draw,
                    config->x_size_2d - 240, config->y_size_2d - 66,
                    240, 66);
    cairo_fill(cairo_draw);

    sprintf(text, "\u00A9 Dominic Ford 2012\u20132019");
    chart_label(cairo_draw, yellow, text, config->x_size_2d - 12, config->y_size_2d - 44, 1, 0, 14, 1, 0);

    sprintf(text, "https://in-the-sky.org/");
    chart_label(cairo_draw, yellow, text, config->x_size_2d - 12, config->y_size_2d - 18, 1, 0, 14, 1, 0);

    // Write duration in the top left corner of the image
    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.5);
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
    sprintf(output_filename, "%s/frameA%06d.png", config->output_dir, config->frame_counter);
    cairo_surface_write_to_png(surface, output_filename);

    // If this frame is at the midpoint of the eclipse we output a special "poster" frame which acts as a teaser image
    // for the animation.
    if (config->frame_counter == config->poster_image_frame) {
        sprintf(text, "Please wait \u2013 loading...");
        chart_label(cairo_draw, yellow, text, (int) (config->x_size_2d / 2), (int) (config->y_size_2d / 2),
                    0, 0, 16, 1, 0);

        // Write out poster image
        sprintf(output_filename, "%s/solarEclipseA.png", config->output_dir);
        cairo_surface_write_to_png(surface, output_filename);
    }

    cairo_destroy(cairo_draw);
    cairo_surface_finish(surface);
    free(pixel_data);
}

void trace_eclipse_contour(cairo_t *cairo_draw, const settings *config, double x_offset,
                           const contour_line *contour) {


    double previous_x = -1, previous_y = -1;
    int pen_up = 1;
    for (int j = 0; j < contour->point_count; j++) {
        // Look up the coordinates of this point
        const double longitude = contour->longitude[j] * 180 / M_PI;
        const double latitude = contour->latitude[j] * 180 / M_PI;

        // Project this point onto the 2D canvas
        const double y_point = (90 - latitude) / 180. * config->y_size_2d;
        double x_point = (longitude + 180.) / 360 * config->x_size_2d - x_offset;

        while (x_point < 0) x_point += config->x_size_2d;
        while (x_point >= config->x_size_2d) x_point -= config->x_size_2d;

        // If we have flipped from one side of the screen onto the other, break the line
        if ((previous_x >= 0) && (fabs(x_point - previous_x) > 0.5 * config->x_size_2d)) {
            cairo_line_to(cairo_draw,
                          (previous_x > 0.5 * config->x_size_2d) ? config->x_size_2d : 0,
                          previous_y);
            cairo_line_to(cairo_draw,
                          (previous_x > 0.5 * config->x_size_2d) ? config->x_size_2d : 0,
                          (previous_y > 0.5 * config->y_size_2d) ? config->y_size_2d : 0);
            cairo_line_to(cairo_draw,
                          (previous_x > 0.5 * config->x_size_2d) ? 0 : config->x_size_2d,
                          (previous_y > 0.5 * config->y_size_2d) ? config->y_size_2d : 0);
            cairo_line_to(cairo_draw,
                          (previous_x > 0.5 * config->x_size_2d) ? 0 : config->x_size_2d,
                          y_point);
        }

        // Add this point to the line
        if (pen_up) {
            cairo_move_to(cairo_draw, x_point, y_point);
            pen_up = 0;
        } else {
            cairo_line_to(cairo_draw, x_point, y_point);
        }
        previous_x = x_point;
        previous_y = y_point;
    }
}

/**
 * render_2d_maximum_extent - Render a flat 2D world map of the maximum extent of the eclipse
 * @param cl [in] - Data structure used to look up which country any given (lat,lng) is within
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param contours [in] - A list of contours of eclipse magnitude
 * @param greatest_shadow [in] - A binary map of the eclipse magnitude across the world
 * @param eclipse_path [in] - The path of greatest eclipse, with duration at each point
 * @param format [in] - The graphics output format to write. Options are png, pdf or svg.
 */
void render_2d_maximum_extent(const country_lookup_handle *cl, const settings *config,
                              const contour_line_list *contours,
                              const shadow_map *greatest_shadow,
                              const eclipse_path_list *eclipse_path, const char *format) {
    int x, y;

    // Generate image RGB data
    const int stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, config->x_size_2d);

    // Allocate buffer for image data
    unsigned char *pixel_data = malloc(stride * config->y_size_2d);

    // Horizontally shift map of the world to put the point of greatest eclipse at the centre
    int x_offset = (int) (config->solar_longitude_at_midpoint * config->x_size_2d / 360);
    if (x_offset > config->x_size_2d / 2) x_offset -= config->x_size_2d;
    if (x_offset < -config->x_size_2d / 2) x_offset += config->x_size_2d;

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            // Work out the latitude and longitude of this pixel on the Earth
            double latitude = 90 - (y * 180. / config->y_size_2d);
            double longitude = (x * 360. / config->x_size_2d) - 180;

            while (longitude < -180) longitude += 360;
            while (longitude >= 180) longitude -= 360;

            // Look up whether this pixel is land or sea
            int country = test_if_land_or_sea(cl, longitude, latitude);

            // Work out color of this pixel
            colour colour_this = {0, 0, 0};

            if (country > 0) {

                // Land
                const colour colour_land = {188, 205, 173};
                colour_this = colour_land;

            } else if (country == 0) {

                // Sea
                const colour colour_sea = {139, 189, 224};
                colour_this = colour_sea;

            }

            // Set pixel color
            const int output_offset = ((x - x_offset + config->x_size_2d) % config->x_size_2d) * 4 + y * stride;

            *(uint32_t *) &pixel_data[output_offset] = ((uint32_t) colour_this.blu +  // blue
                                                        ((uint32_t) colour_this.grn << (unsigned) 8) +  // green
                                                        ((uint32_t) colour_this.red << (unsigned) 16) + // red
                                                        ((uint32_t) 255 << (unsigned) 24)  // alpha
            );
        }

    // Overlay contours of eclipse magnitude on top of the map
    const double contourList[] = {80, 60, 40, 20, 1e-6, -1};
    int *label_position_x, *label_position_y;
    shadowContoursLabelPositions(contourList, greatest_shadow, x_offset,
                                 config->x_size_2d, config->y_size_2d,
                                 &label_position_x, &label_position_y,
                                 NULL, NULL);

    // Turn bitmap data into a Cairo surface
    cairo_surface_t *surface = cairo_image_surface_create_for_data(pixel_data, CAIRO_FORMAT_ARGB32,
                                                                   config->x_size_2d, config->y_size_2d,
                                                                   stride);

    // Create output surface
    char output_filename[FNAME_LENGTH];
    sprintf(output_filename, "%s/maximumEclipse.%s", config->output_dir, format);
    const double dots_per_inch = 200.;
    const double points_per_dot = 72 / dots_per_inch;

    // Create either PNG, SVG or PDF surface
    cairo_surface_t *output_surface = NULL;
    if (strcmp(format, "png") == 0) {
        output_surface = surface;
        surface = NULL;
    } else if (strcmp(format, "pdf") == 0) {
        output_surface = cairo_pdf_surface_create(output_filename, config->x_size_2d * points_per_dot,
                                                  config->y_size_2d * points_per_dot);
    } else if (strcmp(format, "svg") == 0) {
        output_surface = cairo_svg_surface_create(output_filename, config->x_size_2d * points_per_dot,
                                                  config->y_size_2d * points_per_dot);
    }

    // Create drawing context
    cairo_t *cairo_draw = cairo_create(output_surface);

    // If necessary, paste bitmap image onto vector graphics surface
    if (surface != NULL) {
        // We need to do some coordinate transformations in order to display the image at the right scale...

        // Scale our coordinate system so the vector graphics axes match the bitmap graphics axes
        cairo_scale(cairo_draw, points_per_dot, points_per_dot);

        // Paint PNG image onto destination canvas, where it fills up a space of dimensions width x height
        cairo_set_source_surface(cairo_draw, surface, 0, 0);
        cairo_rectangle(cairo_draw, 0, 0, config->x_size_2d, config->y_size_2d);
        cairo_fill(cairo_draw);

        cairo_set_source_rgb(cairo_draw, 0, 0, 0);
    }

    const colour black = {0, 0, 0};

    // Label contours
    for (int i = 0; contourList[i] >= 0; i++)
        if (label_position_x[i] >= 0) {
            const double level = contourList[i];
            char text[8];
            sprintf(text, "%.0f%%", level);

            chart_label(cairo_draw, black, text, label_position_x[i], label_position_y[i],
                        0, 0, 14, 1, 0);
        }

    // Draw path of central eclipse
    {
        cairo_set_source_rgb(cairo_draw, 0, 0, 0);
        cairo_set_line_width(cairo_draw, 1.4);
        cairo_new_path(cairo_draw);

        double previous_x = -1;
        int pen_up = 1;
        for (int i = 0; i < eclipse_path->path_count; i++) {
            for (int j = 0; j < eclipse_path->paths[i].point_count; j++) {
                // Look up the coordinates of this point
                const double longitude = eclipse_path->paths[i].path[j].longitude * 180 / M_PI;
                const double latitude = eclipse_path->paths[i].path[j].latitude * 180 / M_PI;

                // Project this point onto the 2D canvas
                const double y_point = (90 - latitude) / 180. * config->y_size_2d;
                double x_point = (longitude + 180.) / 360 * config->x_size_2d - x_offset;

                while (x_point < 0) x_point += config->x_size_2d;
                while (x_point >= config->x_size_2d) x_point -= config->x_size_2d;

                // If we have flipped from one side of the screen onto the other, break the line
                if ((previous_x >= 0) && (fabs(x_point - previous_x) > 0.5 * config->x_size_2d)) {
                    pen_up = 1;
                }

                // Add this point to the line
                if (pen_up) {
                    cairo_move_to(cairo_draw, x_point, y_point);
                    pen_up = 0;
                } else {
                    cairo_line_to(cairo_draw, x_point, y_point);
                }
                previous_x = x_point;
            }
        }
        cairo_stroke(cairo_draw);
    }

    // Shade eclipse contours
    cairo_save(cairo_draw);
    cairo_set_fill_rule(cairo_draw, CAIRO_FILL_RULE_EVEN_ODD);
    cairo_set_source_rgba(cairo_draw, 0, 0, 0, 0.15);
    cairo_new_path(cairo_draw);
    for (int index = 0; index < contours->contour_count; index++) {
        cairo_new_sub_path(cairo_draw);
        trace_eclipse_contour(cairo_draw, config, x_offset, &contours->line[index]);
        cairo_close_path(cairo_draw);
    }
    cairo_fill(cairo_draw);

    // Make clipping region around contour labels before stroking the lines
    cairo_new_path(cairo_draw);
    for (int i = 0; contourList[i] >= 0; i++)
        if (label_position_x[i] >= 0) {
            cairo_new_sub_path(cairo_draw);
            cairo_arc(cairo_draw, label_position_x[i], label_position_y[i], 18, 0, 2 * M_PI);
        }
    // We want to clip outside the circles where the labels are, so add an extra path around the outside of plot
    cairo_new_sub_path(cairo_draw);
    cairo_rectangle(cairo_draw, 0, 0, config->x_size_2d, config->y_size_2d);
    cairo_clip(cairo_draw);

    // Draw eclipse contours
    for (int index = 0; index < contours->contour_count; index++) {
        const contour_line *contour = &contours->line[index];

        // Set colour and line width for this contour
        cairo_set_source_rgb(cairo_draw, 0, 0, 0);
        cairo_set_line_width(cairo_draw, (contour->eclipse_magnitude > 0.01) ? 1 : 2);
        cairo_new_path(cairo_draw);

        trace_eclipse_contour(cairo_draw, config, x_offset, contour);

        // Close path
        cairo_close_path(cairo_draw);

        // Stroke path
        cairo_stroke(cairo_draw);
    }
    cairo_restore(cairo_draw);

    // Get date components (in UT; not TT)
    int year, month, day, hour, min, status = 0;
    double sec;
    double jd = (config->jd_min + config->jd_max) / 2;
    inv_julian_day(jd - delta_t(jd) / 86400.,
                   &year, &month, &day, &hour, &min, &sec, &status, temp_err_string);

    // Write the time and date in bottom left corner of the image
    char text[FNAME_LENGTH];

    cairo_set_source_rgba(cairo_draw, 1, 1, 1, 0.5);
    cairo_rectangle(cairo_draw,
                    0, config->y_size_2d - 95,
                    190, 95);
    cairo_fill_preserve(cairo_draw);
    cairo_set_line_width(cairo_draw, 1);
    cairo_set_source_rgb(cairo_draw, 0.25, 0.25, 0.25);
    cairo_stroke(cairo_draw);

    sprintf(text, "Greatest Eclipse");
    chart_label(cairo_draw, black, text, 95, config->y_size_2d - 80, 0, 0, 15, 1, 0);

    sprintf(text, "%d %s %d UTC", day, config->month_names[month], year);
    chart_label(cairo_draw, black, text, 95, config->y_size_2d - 17, 0, 0, 17, 1, 0);

    sprintf(text, "%02d:%02d", hour, min);
    chart_label(cairo_draw, black, text, 95, config->y_size_2d - 49, 0, 0, 28, 1, 0);

    // Write copyright text in bottom right corner of the image
    cairo_set_source_rgba(cairo_draw, 1, 1, 1, 0.5);
    cairo_rectangle(cairo_draw,
                    config->x_size_2d - 240, config->y_size_2d - 66,
                    240, 66);
    cairo_fill_preserve(cairo_draw);
    cairo_set_line_width(cairo_draw, 1);
    cairo_set_source_rgb(cairo_draw, 0.25, 0.25, 0.25);
    cairo_stroke(cairo_draw);

    sprintf(text, "\u00A9 Dominic Ford 2012\u20132019");
    chart_label(cairo_draw, black, text, config->x_size_2d - 12, config->y_size_2d - 44, 1, 0, 14, 1, 0);

    sprintf(text, "https://in-the-sky.org/");
    chart_label(cairo_draw, black, text, config->x_size_2d - 12, config->y_size_2d - 18, 1, 0, 14, 1, 0);

    // Write duration in the top left corner of the image
    cairo_set_source_rgba(cairo_draw, 1, 1, 1, 0.5);
    cairo_rectangle(cairo_draw,
                    0, 0,
                    166, (eclipse_path->maximum_duration > 0) ? 70 : 35);
    cairo_fill_preserve(cairo_draw);
    cairo_set_line_width(cairo_draw, 1);
    cairo_set_source_rgb(cairo_draw, 0.25, 0.25, 0.25);
    cairo_stroke(cairo_draw);

    sprintf(text, eclipse_path->maximum_duration > 0 ? "Greatest duration" : "Partial eclipse");
    chart_label(cairo_draw, black, text, 83, 20, 0, 0, 15, 1, 0);

    if (eclipse_path->maximum_duration > 0) {
        sprintf(text, "%dm%02ds", (int) (eclipse_path->maximum_duration / 60),
                (int) eclipse_path->maximum_duration % 60);
        chart_label(cairo_draw, black, text, 125, 50, 0, 0, 15, 1, 0);

        cairo_set_source_rgba(cairo_draw, 1, 1, 1, 0.5);
        cairo_rectangle(cairo_draw,
                        config->x_size_2d - 240, 0,
                        240, 35);
        cairo_fill_preserve(cairo_draw);
        cairo_set_line_width(cairo_draw, 1);
        cairo_set_source_rgb(cairo_draw, 0.25, 0.25, 0.25);
        cairo_stroke(cairo_draw);

        chart_label(cairo_draw, black, config->title, config->x_size_2d - 120, 20, 0, 0, 15, 1, 0);
    }

    // Write output image
    if (strcmp(format, "png") == 0) {
        cairo_surface_write_to_png(output_surface, output_filename);
    }

    cairo_destroy(cairo_draw);
    cairo_surface_finish(output_surface);
    if (surface != NULL) cairo_surface_finish(surface);
    free(label_position_x);
    free(label_position_y);
    free(pixel_data);
}

/**
 * render_2d_eclipse_icon - Render a small icon of the area where this eclipse is visible, as is seen on the page
 * <https://in-the-sky.org/eclipses.php>.
 * @param cl [in] - Data structure used to look up which country any given (lat,lng) is within
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param greatest_shadow [in] - A binary map of the eclipse magnitude across the world
 */
void render_2d_eclipse_icon(const country_lookup_handle *cl, const settings *config,
                            const shadow_map *greatest_shadow) {
    int x, y;

    // Generate image RGB data
    const int stride = cairo_format_stride_for_width(CAIRO_FORMAT_ARGB32, config->x_size_teaser);

    // Allocate buffer for image data
    unsigned char *pixel_data = malloc(stride * config->y_size_teaser);

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_teaser; y++)
        for (x = 0; x < config->x_size_teaser; x++) {
            // Work out the latitude and longitude of this pixel on the Earth
            double latitude = 90 - (y * 180. / config->y_size_teaser);
            double longitude = (x * 360. / config->x_size_teaser) - 180;

            while (longitude < -180) longitude += 360;
            while (longitude >= 180) longitude -= 360;

            // Look up whether this pixel is land or sea
            int country = test_if_land_or_sea(cl, longitude, latitude);

            // Work out the offset of this point within the shadow map, which has dimensions (x_size_2d by y_size_2d)
            const int offset = (int) (floor(x * (double) config->x_size_2d / config->x_size_teaser) +
                                      floor(y * (double) config->y_size_2d / config->y_size_teaser) *
                                      config->x_size_2d);
            double shadow = greatest_shadow->map[offset];

            const int point_within_eclipse = (shadow > 0.001);

            // Work out color of this pixel
            colour colour_this = {0, 0, 0};

            if (point_within_eclipse && (country > 0)) {

                // Land, with eclipse
                const colour colour_eclipsed_land = {150, 230, 150};
                colour_this = colour_eclipsed_land;

            } else if ((!point_within_eclipse) && (country > 0)) {

                // Land, without eclipse
                const colour colour_uneclipsed_land = {113, 173, 113};
                colour_this = colour_uneclipsed_land;

            } else if (point_within_eclipse && (country == 0)) {

                // Sea, with eclipse
                const colour colour_eclipsed_sea = {160, 160, 255};
                colour_this = colour_eclipsed_sea;

            } else if ((!point_within_eclipse) && (country == 0)) {

                // Sea, without eclipse
                const colour colour_uneclipsed_sea = {120, 120, 192};
                colour_this = colour_uneclipsed_sea;

            }

            // Set pixel color
            const int output_offset = x * 4 + y * stride;

            *(uint32_t *) &pixel_data[output_offset] = ((uint32_t) colour_this.blu +  // blue
                                                        ((uint32_t) colour_this.grn << (unsigned) 8) +  // green
                                                        ((uint32_t) colour_this.red << (unsigned) 16) + // red
                                                        ((uint32_t) 255 << (unsigned) 24)  // alpha
            );
        }

    // Turn bitmap data into a Cairo surface
    cairo_surface_t *surface = cairo_image_surface_create_for_data(pixel_data, CAIRO_FORMAT_ARGB32,
                                                                   config->x_size_teaser, config->y_size_teaser,
                                                                   stride);

    cairo_t *cairo_draw = cairo_create(surface);

    // Write output image
    char output_filename[FNAME_LENGTH];
    sprintf(output_filename, "%s/solarEclipseC.png", config->output_dir);
    cairo_surface_write_to_png(surface, output_filename);

    cairo_destroy(cairo_draw);
    cairo_surface_finish(surface);
    free(pixel_data);
}
