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

#include <gd.h>
#include <gsl/gsl_math.h>

#include "jpeg_in.h"

#include "coreUtils/errorReport.h"
#include "coreUtils/strConstants.h"

#include "mathsTools/julianDate.h"

#include "map_greatest_eclipse.h"
#include "render_2d.h"
#include "render_3d.h"
#include "settings.h"
#include "shadow_calc.h"

/**
 * render_3d_eclipse_map - Render a 3D snapshot of a globe of the world as viewed from the direction of the Sun, with
 * the magnitude of a solar eclipse overlaid as contours. This is ideal for turning into an animation of the eclipse's
 * process if snapshots are turned into a video.
 *
 * @param config [in] - The settings for this eclipse simulation, including, e.g. the output image size
 * @param jd [in] - The Julian day number of the current point in the simulation
 * @param earthDay [in] - A JPEG image of the world in daylight
 * @param shadow_map [in] - A binary map of the eclipse magnitude across the world
 * @param eclipse_path [in] - The path of greatest eclipse, with duration at each point
 */
void render_3d_eclipse_map(settings *config, double jd, jpeg_ptr earthDay,
                           const shadow_map *shadow_map, const eclipse_path_list *eclipse_path) {
    int x, y;

    // Draw spherical map of Earth onto this canvas
    gdImagePtr frame = gdImageCreateTrueColor(config->x_size_3d, config->y_size_3d);

    // Loop over all the pixels in the image we are to produce
    for (y = 0; y < config->y_size_3d; y++)
        for (x = 0; x < config->x_size_3d; x++) {
            // Pixel's offset inside shadow_map
            const int offset = x + y * config->x_size_3d;

            jpeg_ptr *srcimg = &earthDay;

            // Look up eclipse magnitude in this pixel
            double shadow = shadow_map->map[offset];

            unsigned int c0 = 0, c1 = 0, c2 = 0;

            if (!gsl_isnan(shadow_map->lng[offset])) {
                // Lat & lng of this pixel are finite, which means pixel lies on surface of Earth
                // Otherwise this pixel lies off the side of the globe of the world
                int p0 = (int) ((shadow_map->lng[offset] + 180.) * srcimg->xsize / 360.);
                int p1 = (int) ((90. - shadow_map->lat[offset]) * srcimg->ysize / 180.);
                if (p0 >= srcimg->xsize) p0 -= srcimg->xsize;
                if (p1 >= srcimg->ysize) p1 = srcimg->ysize - 1;

                c0 = ((int) srcimg->data_red[p0 + p1 * srcimg->xsize]);
                c1 = ((int) srcimg->data_grn[p0 + p1 * srcimg->xsize]);
                c2 = ((int) srcimg->data_blu[p0 + p1 * srcimg->xsize]);
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
    drawShadowContours(frame, shadow_map, config->x_size_3d, config->y_size_3d);

    // Look up duration of the eclipse
    int is_total;
    const double duration = eclipse_duration_from_path(eclipse_path, jd, &is_total);

    // Get date components
    int year, month, day, hour, min, status = 0;
    double sec;
    inv_julian_day(jd, &year, &month, &day, &hour, &min, &sec, &status, temp_err_string);

    // Write the time and date in bottom left corner of the image
    int brect[8];
    char scratch[FNAME_LENGTH];
    //gdImageFilledRectangle(frame,0,y_size_3d-70,190,y_size_3d,gdTrueColorAlpha(0,0,0, (int)(gdAlphaTransparent*0.5)));
    sprintf(scratch, "%d %s %d UTC", day, config->month_names[month], year);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 15, 0, 10, config->y_size_3d - 11,
                    scratch);
    sprintf(scratch, "%02d:%02d", hour, min);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 24, 0, 28, config->y_size_3d - 38,
                    scratch);

    // Write copyright text in bottom right corner of the image
    sprintf(scratch, "https://in-the-sky.org/");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, 0, 0, scratch);
    //gdImageFilledRectangle(frame,x_size_2d-brect[4]-16,y_size_2d-70,x_size_2d,y_size_2d,gdTrueColorAlpha(0,0,0, (int)(gdAlphaTransparent*0.5)));
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, config->x_size_3d - brect[4] - 8,
                    config->y_size_3d - 11,
                    scratch);
    sprintf(scratch, "&copy; Dominic Ford 2012-2019");
    gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, 0, 0, scratch);
    gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 12, 0, config->x_size_3d - brect[4] - 8,
                    config->y_size_3d - 38,
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
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, 33 - brect[4] / 2, 52,
                        scratch);

        sprintf(scratch, "%dm%02ds", (int)(duration / 60), (int)duration % 60);
        gdImageStringFT(NULL, brect, gdTrueColor(255, 255, 0), config->font_name, 15, 0, 0, 0, scratch);
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 15, 0, 142 - brect[4], 52,
                        scratch);
    } else {
        gdImageStringFT(frame, brect, gdTrueColor(255, 255, 0), config->font_name, 13, 0, brect[4] / 2, 52,
                        "&ndash;");
    }

    // Write output image
    char outfname[FNAME_LENGTH];
    sprintf(outfname, "%s/frameB%06d.jpg", config->output_dir, config->frame_counter);
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
                        (config->x_size_3d - brect[4]) / 2,
                        (config->y_size_3d - brect[7]) / 2, scratch);

        // Write out poster image
        char outfname[FNAME_LENGTH];
        sprintf(outfname, "%s/solarEclipseB.jpg", config->output_dir);
        FILE *f = fopen(outfname, "wb");
        gdImageJpeg(frame, f, 95);
        fclose(f);
    }

    gdImageDestroy(frame);
}
