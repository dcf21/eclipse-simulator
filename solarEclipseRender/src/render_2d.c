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

#include <math.h>
#include <gd.h>
#include "jpeg_in.h"

#include "coreUtils/errorReport.h"
#include "coreUtils/strConstants.h"

#include "mathsTools/julianDate.h"

#include "country_lookup.h"
#include "map_greatest_eclipse.h"
#include "render_2d.h"
#include "settings.h"
#include "shadow_calc.h"

/**
 * drawShadowContours - Draw contours onto a GD image surface, wherever the values in the array <shadow> pass any
 * of the thresholds in the array <contourList>.
 * @param frame - The GD image surface onto which to trace the contours.
 * @param shadow - The array of shadow fractions within each pixel, which we are to draw contours from
 * @param xs - The horizontal pixel size of the GD image
 * @param ys - The vertical pixel size of the GD image
 */
void drawShadowContours(gdImagePtr frame, const shadow_map *shadow, int xs, int ys) {
    int x, y, i;
    double contourList[] = {80, 60, 40, 20, 0.0001, -1};

    for (i = 0; contourList[i] >= 0; i++) {
        double level = contourList[i] / 100.;
        int color = (level > 0.01) ? gdTrueColor(200, 200, 0) : gdTrueColor(220, 0, 0);
        for (y = 1; y < ys - 2; y++)
            for (x = 1; x < xs - 2; x++) {
                if (shadow->map[y * xs + x] > level)
                    if ((shadow->map[(y - 1) * xs + (x)] <= level) || (shadow->map[(y - 1) * xs + (x - 1)] <= level) ||
                        (shadow->map[(y) * xs + (x - 1)] <= level) || (shadow->map[(y + 1) * xs + (x - 1)] <= level) ||
                        (shadow->map[(y + 1) * xs + (x)] <= level) || (shadow->map[(y + 1) * xs + (x + 1)] <= level) ||
                        (shadow->map[(y) * xs + (x + 1)] <= level) || (shadow->map[(y - 1) * xs + (x + 1)] <= level)) {
                        gdImageSetPixel(frame, x, y, color);
                        gdImageSetPixel(frame, x + 1, y, color);
                        gdImageSetPixel(frame, x, y + 1, color);
                        gdImageSetPixel(frame, x + 1, y + 1, color);
                    }
            }
    }
}

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

    // Draw flat map of Earth onto this canvas
    gdImagePtr frame = gdImageCreateTrueColor(config->x_size_2d, config->y_size_2d);

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            const int offset = x + y * config->x_size_2d;

            jpeg_ptr *srcimg_day = &earthDay;
            jpeg_ptr *srcimg_night = &earthNight;

            // Fetch colour of day time pixel
            int p0_day = (int) ((shadow_map->lng[offset] + 180.) * srcimg_day->xsize / 360.);
            int p1_day = (int) ((90. - shadow_map->lat[offset]) * srcimg_day->ysize / 180.);
            if (p0_day >= srcimg_day->xsize) p0_day -= srcimg_day->xsize;
            if (p1_day >= srcimg_day->ysize) p1_day = srcimg_day->ysize - 1;
            int c0_day = ((int) srcimg_day->data_red[p0_day + p1_day * srcimg_day->xsize]);
            int c1_day = ((int) srcimg_day->data_grn[p0_day + p1_day * srcimg_day->xsize]);
            int c2_day = ((int) srcimg_day->data_blu[p0_day + p1_day * srcimg_day->xsize]);

            // Fetch colour of night time pixel
            int p0_night = (int) ((shadow_map->lng[offset] + 180.) * srcimg_night->xsize / 360.);
            int p1_night = (int) ((90. - shadow_map->lat[offset]) * srcimg_night->ysize / 180.);
            if (p0_night >= srcimg_night->xsize) p0_night -= srcimg_night->xsize;
            if (p1_night >= srcimg_night->ysize) p1_night = srcimg_night->ysize - 1;
            int c0_night = ((int) srcimg_night->data_red[p0_night + p1_night * srcimg_night->xsize]);
            int c1_night = ((int) srcimg_night->data_grn[p0_night + p1_night * srcimg_night->xsize]);
            int c2_night = ((int) srcimg_night->data_blu[p0_night + p1_night * srcimg_night->xsize]);

            // Look up eclipse magnitude in the current pixel
            const double shadow = shadow_map->map[x + y * config->x_size_2d];

            // Test whether this pixel is on the day side of the Earth, or the night side
            const int night_time = (shadow < 0);

            int c0, c1, c2;
            if (night_time) {
                c0 = c0_night * 0.8 + c0_day * 0.2;
                c1 = c1_night * 0.8 + c1_day * 0.2;
                c2 = c2_night * 0.8 + c2_day * 0.2;
            } else { // Day time
                c0 = c0_night * 0.1 + c0_day * 0.9;
                c1 = c1_night * 0.1 + c1_day * 0.9;
                c2 = c2_night * 0.1 + c2_day * 0.9;
            }

            // Superimpose shadow map over Earth
            if (shadow > 0.98) {
                // If this pixel experiences a total eclipse, shade it accordingly
                c0 = (int) (config->moon_shadow_fade_fraction * c0 +
                            (1 - config->moon_shadow_fade_fraction) * config->totality_col_r);
                c1 = (int) (config->moon_shadow_fade_fraction * c1 +
                            (1 - config->moon_shadow_fade_fraction) * config->totality_col_g);
                c2 = (int) (config->moon_shadow_fade_fraction * c2 +
                            (1 - config->moon_shadow_fade_fraction) * config->totality_col_b);
            } else if (shadow > 0) {
                // If this pixel experiences a partial eclipse, shade it accordingly
                c0 = (int) (config->moon_shadow_fade_fraction * c0 +
                            (1 - config->moon_shadow_fade_fraction) * config->shadow_col_r);
                c1 = (int) (config->moon_shadow_fade_fraction * c1 +
                            (1 - config->moon_shadow_fade_fraction) * config->shadow_col_g);
                c2 = (int) (config->moon_shadow_fade_fraction * c2 +
                            (1 - config->moon_shadow_fade_fraction) * config->shadow_col_b);
            }

            // Set pixel color
            int color = gdTrueColor(c0, c1, c2);
            gdImageSetPixel(frame, x, y, color);
        }

    // Overlay contours of eclipse magnitude on top of the map
    drawShadowContours(frame, shadow_map, config->x_size_2d, config->y_size_2d);

    // Look up duration of the eclipse
    int is_total;
    const double duration = eclipse_duration_from_path(eclipse_path, jd, &is_total);

    // Get date components (in UT; not TT)
    int year, month, day, hour, min, status = 0;
    double sec;
    inv_julian_day(jd - delta_t(jd) / 86400.,
                   &year, &month, &day, &hour, &min, &sec, &status, temp_err_string);

    // Write the time and date in bottom left corner of the image
    int brect[8];
    char scratch[FNAME_LENGTH];
    gdImageFilledRectangle(frame, 0, config->y_size_2d - 70, 175, config->y_size_2d,
                           gdTrueColorAlpha(0, 0, 0, (int) (gdAlphaTransparent * 0.5)));
    sprintf(scratch, "%d %s %d UTC", day, config->month_names[month], year);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 15, 0, 10, config->y_size_2d - 11,
                    scratch);
    sprintf(scratch, "%02d:%02d", hour, min);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 24, 0, 28, config->y_size_2d - 38,
                    scratch);

    // Write copyright text in bottom right corner of the image
    sprintf(scratch, "&copy; Dominic Ford 2012-2019");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, 0, 0, scratch);
    gdImageFilledRectangle(frame, config->x_size_2d - brect[4] - 16, config->y_size_2d - 70, config->x_size_2d,
                           config->y_size_2d,
                           gdTrueColorAlpha(0, 0, 0, (int) (gdAlphaTransparent * 0.5)));
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, config->x_size_2d - brect[4] - 8,
                    config->y_size_2d - 38,
                    scratch);

    sprintf(scratch, "https://in-the-sky.org/");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, 0, 0, scratch);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, config->x_size_2d - brect[4] - 8,
                    config->y_size_2d - 11,
                    scratch);

    // Write duration in the top left corner of the image
    sprintf(scratch, "Central duration");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, 0, 0, scratch);
    gdImageFilledRectangle(frame, 0, 0, brect[4] + 16, 70,
                           gdTrueColorAlpha(0, 0, 0, (int) (gdAlphaTransparent * 0.5)));
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, 8, 20, scratch);

    if (duration > 0) {
        sprintf(scratch, is_total ? "Total" : "Annular");
        gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, 0, 0, scratch);
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, 41 - brect[4] / 2, 52,
                        scratch);

        sprintf(scratch, "%dm%02ds", (int) (duration / 60), (int) duration % 60);
        gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 15, 0, 0, 0, scratch);
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 15, 0, 146 - brect[4], 52,
                        scratch);
    } else {
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, brect[4] / 2, 52,
                        "&ndash;");
    }

    // Write output image
    char outfname[FNAME_LENGTH];
    sprintf(outfname, "%s/frameA%06d.jpg", config->output_dir, config->frame_counter);
    FILE *f = fopen(outfname, "wb");
    gdImageJpeg(frame, f, 95);
    fclose(f);

    // If this frame is at the midpoint of the eclipse we output a special "poster" frame which acts as a teaser image
    // for the animation.
    if (config->frame_counter == config->poster_image_frame) {
        int brect[8];
        char scratch[FNAME_LENGTH];
        sprintf(scratch, "Please wait &ndash; loading...");
        gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 16, 0, 0, 0, scratch);
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 16, 0,
                        (config->x_size_2d - brect[4]) / 2,
                        (config->y_size_2d - brect[7]) / 2, scratch);

        // Write out poster image
        char outfname[FNAME_LENGTH];
        sprintf(outfname, "%s/solarEclipseA.jpg", config->output_dir);
        FILE *f = fopen(outfname, "wb");
        gdImageJpeg(frame, f, 95);
        fclose(f);
    }

    gdImageDestroy(frame);
}

/**
 * render_2d_maximum_extent - Render a flat 2D world map of the maximum extent of the eclipse
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param earthDay [in] - A JPEG image of the world in daylight
 * @param greatest_shadow [in] - A binary map of the eclipse magnitude across the world
 */
void render_2d_maximum_extent(settings *config, jpeg_ptr earthDay, const shadow_map *greatest_shadow) {
    int x, y;

    // Make map of maximum eclipse extent
    gdImagePtr frame = gdImageCreateTrueColor(config->x_size_2d, config->y_size_2d);

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            // Work out the latitude and longitude of this pixel on the Earth
            double lat = 90 - (y * 180. / config->y_size_2d);
            double lng = (x * 360. / config->x_size_2d) + (-180) + 360;
            while (lng > 180) lng -= 360;

            // Convert this into a pixel position on the image <earthDay>
            int p0 = (int) ((lng + 180.) * earthDay.xsize / 360.);
            int p1 = (int) ((90. - lat) * earthDay.ysize / 180.);

            // Make sure longitude wraps sensibly
            if (p0 >= earthDay.xsize) p0 -= earthDay.xsize;
            if (p1 >= earthDay.ysize) p1 = earthDay.ysize - 1;

            // Look up the color of this pixel in <earthDay>
            int c0 = ((int) earthDay.data_red[p0 + p1 * earthDay.xsize]);
            int c1 = ((int) earthDay.data_grn[p0 + p1 * earthDay.xsize]);
            int c2 = ((int) earthDay.data_blu[p0 + p1 * earthDay.xsize]);

            // Look up the maximum eclipse fraction in this pixel
            double shadow = greatest_shadow->map[x + y * config->x_size_2d];

            if (shadow > 0.98) {
                // If this pixel experiences a total eclipse, shade it accordingly
                c0 = config->moon_shadow_fade_fraction * c0 +
                     (1 - config->moon_shadow_fade_fraction) * config->totality_col_r;
                c1 = config->moon_shadow_fade_fraction * c1 +
                     (1 - config->moon_shadow_fade_fraction) * config->totality_col_g;
                c2 = config->moon_shadow_fade_fraction * c2 +
                     (1 - config->moon_shadow_fade_fraction) * config->totality_col_b;
            } else if (shadow > 0.001) {
                // If this pixel experiences a partial eclipse, shade it accordingly
                c0 = config->moon_shadow_fade_fraction * c0 +
                     (1 - config->moon_shadow_fade_fraction) * config->shadow_col_r;
                c1 = config->moon_shadow_fade_fraction * c1 +
                     (1 - config->moon_shadow_fade_fraction) * config->shadow_col_g;
                c2 = config->moon_shadow_fade_fraction * c2 +
                     (1 - config->moon_shadow_fade_fraction) * config->shadow_col_b;
            }

            // Set pixel color
            int color = gdTrueColor(c0, c1, c2);
            gdImageSetPixel(frame, x, y, color);
        }

    // Overlay contours of eclipse magnitude on top of the map
    drawShadowContours(frame, greatest_shadow, config->x_size_2d, config->y_size_2d);

    // Write the time and date in bottom corners of the image
    int brect[8];
    char scratch[FNAME_LENGTH];

    sprintf(scratch, "&copy; Dominic Ford 2012-2019");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, 0, 0, scratch);
    gdImageFilledRectangle(frame, config->x_size_2d - brect[4] - 16, config->y_size_2d - 70, config->x_size_2d,
                           config->y_size_2d,
                           gdTrueColorAlpha(0, 0, 0, (int) (gdAlphaTransparent * 0.5)));
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, config->x_size_2d - brect[4] - 8,
                    config->y_size_2d - 35, scratch);

    sprintf(scratch, "https://in-the-sky.org/");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, 0, 0, scratch);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, config->x_size_2d - brect[4] - 8,
                    config->y_size_2d - 11, scratch);

    // Write output image
    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/maximumEclipse.jpg", config->output_dir);
    FILE *f = fopen(fname, "wb");
    gdImageJpeg(frame, f, 95);
    fclose(f);
    gdImageDestroy(frame);
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

    // Make an icon map of maximum eclipse extent
    gdImagePtr frame = gdImageCreateTrueColor(config->x_size_teaser, config->y_size_teaser);

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_teaser; y++)
        for (x = 0; x < config->x_size_teaser; x++) {
            // Work out the latitude and longitude of this pixel on the Earth
            double lat = 90 - (y * 180. / config->y_size_teaser);
            double lng = (x * 360. / config->x_size_teaser) + (-180) + 360;

            // Look up whether this pixel is land or sea
            int country = test_if_land_or_sea(cl, lng, lat);

            // Work out the offset of this point within the shadow map, which has dimensions (x_size_2d by y_size_2d)
            const int offset = (int) (floor(x * (double) config->x_size_2d / config->x_size_teaser) +
                                      floor(y * (double) config->y_size_2d / config->y_size_teaser) *
                                      config->x_size_2d);
            double shadow = greatest_shadow->map[offset];

            const int point_within_eclipse = (shadow > 0.001);

            // Work out color of this pixel
            int c0 = 0, c1 = 0, c2 = 0;

            if (point_within_eclipse && (country > 0)) {
                // Land, with eclipse
                c0 = 150;
                c1 = 230;
                c2 = 150;
            } else if ((!point_within_eclipse) && (country > 0)) {
                // Land, without eclipse
                c0 = 90;
                c1 = 138;
                c2 = 90;
            } else if (point_within_eclipse && (country == 0)) {
                // Sea, with eclipse
                c0 = 160;
                c1 = 160;
                c2 = 255;
            } else if ((!point_within_eclipse) && (country == 0)) {
                // Sea, without eclipse
                c0 = 96;
                c1 = 96;
                c2 = 153;
            }

            // Set pixel color
            int color = gdTrueColor(c0, c1, c2);
            gdImageSetPixel(frame, x, y, color);
        }

    // Write output image
    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/solarEclipseC.png", config->output_dir);
    FILE *f = fopen(fname, "wb");
    gdImagePng(frame, f);
    fclose(f);
    gdImageDestroy(frame);
}
