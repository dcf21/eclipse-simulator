// country_lookup.c

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "coreUtils/asciiDouble.h"
#include "coreUtils/strConstants.h"
#include "coreUtils/errorReport.h"

#include "png/image.h"
#include "mathsTools/julianDate.h"

#include "country_lookup.h"
#include "settings.h"
#include "shadow_calc.h"
#include "map_greatest_eclipse.h"

//! country_lookup_init - Load externally-created PNG file where the colours of pixels indicate the countries which
//! lie at any given latitude and longitude.
//!
//! \return - A country lookup handle object

country_lookup_handle *country_lookup_init() {
    int j, x, y;
    FILE *f;

    // Create country lookup handle
    country_lookup_handle *cl = (country_lookup_handle *) malloc(sizeof(country_lookup_handle));

    // Read world map image
    image_ptr worldMapImg = image_get(SRCDIR "../worldMap/worldMap.png");

    // Make sure that handle object is empty to begin with
    for (j = 0; j < COUNTRYLISTLEN; j++) {
        cl->countryList[j].red = cl->countryList[j].grn = cl->countryList[j].blu = 0;

        cl->countryList[j].max_eclipse =
        cl->countryList[j].jd_partial_start =
        cl->countryList[j].jd_total_start =
        cl->countryList[j].jd_total_end =
        cl->countryList[j].jd_partial_end = 0;
    }

    // This text file lists the colours used to shade each country in the PNG file
    f = fopen(SRCDIR "../worldMap/countryList.dat", "r");
    if (f == NULL) {
        logging_fatal(__FILE__, __LINE__, "Could not open list of country info.");
        exit(1);
    }

    // Loop over the contents of the text file, line by line
    while ((!feof(f)) && (!ferror(f))) {
        char *lineptr, line[FNAME_LENGTH];
        file_readline(f, line);
        str_strip(line, line);
        if ((line[0] == '\0') || (line[0] == '#')) continue;
        lineptr = line;
        int uid = (int) get_float(lineptr, NULL);
        if ((uid <= 0) || (uid >= COUNTRYLISTLEN)) {
            logging_fatal(__FILE__, __LINE__, "Illegal country UID.");
            exit(1);
        }

        // Read the colour used to shade each country
        lineptr = next_word(lineptr);
        cl->countryList[uid].red = (int) get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].grn = (int) get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].blu = (int) get_float(lineptr, NULL);

        // Read the latitude and longitude of the capital city of each country
        lineptr = next_word(lineptr);
        cl->countryList[uid].lat = get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].lng = get_float(lineptr, NULL);

        // Read the name of the country
        lineptr = next_word(lineptr);
        strcpy(cl->countryList[uid].name, lineptr);
    }
    fclose(f);

    // Convert world map PNG image to an array of ints
    cl->worldMapWidth = worldMapImg.xsize;
    cl->worldMapHeight = worldMapImg.ysize;
    cl->worldMapArray = malloc(cl->worldMapWidth * cl->worldMapHeight * sizeof(int));

    for (y = 0; y < cl->worldMapHeight; y++)
        for (x = 0; x < cl->worldMapWidth; x++) {
            const int offset = x + y * worldMapImg.xsize;
            const unsigned char red = (unsigned char) worldMapImg.data_red[offset];
            const unsigned char grn = (unsigned char) worldMapImg.data_grn[offset];
            const unsigned char blu = (unsigned char) worldMapImg.data_blu[offset];
            cl->worldMapArray[y * cl->worldMapWidth + x] = 0;
            for (j = 1000; j < COUNTRYLISTLEN; j++)
                if (red || grn || blu)
                    if ((cl->countryList[j].red == red) && (cl->countryList[j].grn == grn) &&
                        (cl->countryList[j].blu == blu)) {
                        cl->worldMapArray[y * cl->worldMapWidth + x] = j;
                        break;
                    }
        }
    image_dealloc(&worldMapImg);

    // Return the country lookup handle object we have just populated
    return cl;
}

//! country_lookup_free - Free the storage associated with a country lookup object
//!
//! @param cl [in] - The country lookup object to free

void country_lookup_free(country_lookup_handle *cl) {
    free(cl->worldMapArray);
    free(cl);
}

//! test_if_land_or_sea - Test whether a particular latitude and longitude on Earth lies on land or sea
//!
//! \param cl [in] - The country lookup object to use
//! \param longitude [in] - The longitude to query, in degrees
//! \param latitude [in] - The latitude to query, in degrees
//! \return Boolean flag, indicating 1 for land, or 0 for sea.

int test_if_land_or_sea(const country_lookup_handle *cl, double longitude, double latitude) {
    while (longitude > 180) longitude -= 360;
    const int p0 = (int) ((longitude + 180.) * cl->worldMapWidth / 360.);
    const int p1 = (int) ((90. - latitude) * cl->worldMapHeight / 180.);
    const int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];

    return country > 0;
}

//! country_lookup_max_eclipse - Compile a JSON file of the maximum extent of the eclipse in each country
//!
//! \param cl [in] - The country lookup object to use
//! \param config [in] - The settings for this eclipse simulation
//! \param greatest_shadow [in] - A binary map of the eclipse magnitude across the world
//! \param eclipse_path [in] - The path of greatest eclipse, with duration at each point

void country_lookup_max_eclipse(country_lookup_handle *cl, const settings *config, const shadow_map *greatest_shadow,
                                const eclipse_path_list *eclipse_path) {
    int x, y;
    FILE *f;

    // Work out the maximum extent of the eclipse by country

    // First look up the extent of the eclipse in the capital city of each country
    for (x = 0; x < COUNTRYLISTLEN; x++)
        if ((cl->countryList[x].red) || (cl->countryList[x].grn) || (cl->countryList[x].blu)) {
            const double lat = cl->countryList[x].lat;
            const double lng = cl->countryList[x].lng;
            int p0 = (int) ((lng + 180.) * config->x_size_2d / 360.);
            int p1 = (int) ((90. - lat) * config->y_size_2d / 180.);
            if (p0 >= config->x_size_2d) p0 -= config->x_size_2d;
            if (p0 < 0) p0 += config->x_size_2d;
            if (p1 >= config->y_size_2d) p1 = config->y_size_2d - 1;

            const int offset = p0 + p1 * config->x_size_2d;
            double shadow = greatest_shadow->map[offset];
            cl->countryList[x].max_eclipse = shadow;
            cl->countryList[x].jd_partial_start = greatest_shadow->jd_partial_start[offset];
            cl->countryList[x].jd_total_start = greatest_shadow->jd_total_start[offset];
            cl->countryList[x].jd_total_end = greatest_shadow->jd_total_end[offset];
            cl->countryList[x].jd_partial_end = greatest_shadow->jd_partial_end[offset];
        }

    // Loop over the pixels in the 2D bitmap map of the greatest extent of the eclipse, noting the greatest value
    // within each country
    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            const int offset = x + y * config->x_size_2d;
            double shadow = greatest_shadow->map[offset];
            if (shadow > 0) {
                const int p0 = (int) (((double) x) / config->x_size_2d * cl->worldMapWidth);
                const int p1 = (int) (((double) y) / config->y_size_2d * cl->worldMapHeight);
                const int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];

                // Update maximum magnitude of eclipse in this country
                if (country > 0) {
                    if (cl->countryList[country].max_eclipse < shadow) {
                        cl->countryList[country].max_eclipse = shadow;
                    }

                    // Update earliest time of partial eclipse in this country
                    if ((greatest_shadow->jd_partial_start[offset] > 0) &&
                        ((cl->countryList[country].jd_partial_start == 0) ||
                         (greatest_shadow->jd_partial_start[offset] < cl->countryList[country].jd_partial_start))) {
                        cl->countryList[country].jd_partial_start = greatest_shadow->jd_partial_start[offset];
                    }

                    // Update earliest time of total eclipse in this country
                    if ((greatest_shadow->jd_total_start[offset] > 0) &&
                        ((cl->countryList[country].jd_total_start == 0) ||
                         (greatest_shadow->jd_total_start[offset] < cl->countryList[country].jd_total_start))) {
                        cl->countryList[country].jd_total_start = greatest_shadow->jd_total_start[offset];
                    }

                    // Update latest time of total eclipse in this country
                    if ((greatest_shadow->jd_total_end[offset] > 0) &&
                        ((cl->countryList[country].jd_total_end == 0) ||
                         (greatest_shadow->jd_total_end[offset] > cl->countryList[country].jd_total_end))) {
                        cl->countryList[country].jd_total_end = greatest_shadow->jd_total_end[offset];
                    }

                    // Update latest time of partial eclipse in this country
                    if ((greatest_shadow->jd_partial_end[offset] > 0) &&
                        ((cl->countryList[country].jd_partial_end == 0) ||
                         (greatest_shadow->jd_partial_end[offset] > cl->countryList[country].jd_partial_end))) {
                        cl->countryList[country].jd_partial_end = greatest_shadow->jd_partial_end[offset];
                    }
                }
            }
        }

    // Now trace along the path of the central eclipse
    // Shadow fraction of 101% means an annular eclipse
    // Shadow fraction of 102% means a total eclipse
    for (int i = 0; i < eclipse_path->path_count; i++) {
        const int is_total = eclipse_path->paths[i].is_total;
        for (int j = 0; j < eclipse_path->paths[i].point_count; j++) {
            // Look up the coordinates of this point
            const path_point *point = &eclipse_path->paths[i].path[j];
            const double longitude = point->longitude * 180 / M_PI;
            const double latitude = point->latitude * 180 / M_PI;

            int p0 = (int) ((longitude + 180.) * cl->worldMapWidth / 360.);
            int p1 = (int) ((90. - latitude) * cl->worldMapHeight / 180.);

            if (p0 >= cl->worldMapWidth) p0 -= cl->worldMapWidth;
            if (p0 < 0) p0 += cl->worldMapWidth;
            if (p1 >= cl->worldMapHeight) p1 = cl->worldMapHeight - 1;

            const int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];

            if (country > 0) {
                const double shadow_new = is_total ? 1.02 : 1.01;

                // Update maximum magnitude of eclipse in this country
                if (shadow_new > cl->countryList[country].max_eclipse) {
                    cl->countryList[country].max_eclipse = shadow_new;
                }

                // Update earliest time of total eclipse in this country
                if ((cl->countryList[country].jd_total_start == 0) ||
                    (point->jd < cl->countryList[country].jd_total_start)) {
                    cl->countryList[country].jd_total_start = point->jd;
                }

                // Update latest time of total eclipse in this country
                if ((cl->countryList[country].jd_total_end == 0) ||
                    (point->jd > cl->countryList[country].jd_total_end)) {
                    cl->countryList[country].jd_total_end = point->jd;
                }
            }
        }
    }

    // Start writing JSON output
    {
        char output_filename[FNAME_LENGTH];
        sprintf(output_filename, "%s/maximumEclipse.json", config->output_dir);
        f = fopen(output_filename, "w");
        fprintf(f, "{\"byuid\":{");
        int first = 1;
        for (x = 1000; x < 2000; x++)
            if ((cl->countryList[x].red) || (cl->countryList[x].grn) || (cl->countryList[x].blu))
                if (cl->countryList[x].max_eclipse > 0.01) {
                    if (!first) fprintf(f, ",");

                    // Time limits of eclipse are output as unix times, in UT (not TT)
                    fprintf(f, "\"%d\":[%d,%.0f,%.0f,%.0f,%.0f]",
                            x, (int) (cl->countryList[x].max_eclipse * 100),

                            ((cl->countryList[x].jd_partial_start == 0) ? 0 :
                             (unix_from_jd(cl->countryList[x].jd_partial_start) -
                              delta_t(cl->countryList[x].jd_partial_start))),

                            ((cl->countryList[x].jd_total_start == 0) ? 0 :
                             (unix_from_jd(cl->countryList[x].jd_total_start) -
                              delta_t(cl->countryList[x].jd_total_start))),

                            ((cl->countryList[x].jd_total_end == 0) ? 0 :
                             (unix_from_jd(cl->countryList[x].jd_total_end) -
                              delta_t(cl->countryList[x].jd_total_end))),

                            ((cl->countryList[x].jd_partial_end == 0) ? 0 :
                             (unix_from_jd(cl->countryList[x].jd_partial_end) -
                              delta_t(cl->countryList[x].jd_partial_end)))
                    );
                    first = 0;
                }
        fprintf(f, "},\"extras\":{");
        first = 1;
        for (x = 2000; x < COUNTRYLISTLEN; x++)
            if (cl->countryList[x].max_eclipse > 0.01) {
                if (!first) fprintf(f, ",");

                // Time limits of eclipse are output as unix times, in UT (not TT)
                fprintf(f, "\"%s\":[%d,%.0f,%.0f,%.0f,%.0f]",
                        cl->countryList[x].name, (int) (cl->countryList[x].max_eclipse * 100),

                        ((cl->countryList[x].jd_partial_start == 0) ? 0 :
                         (unix_from_jd(cl->countryList[x].jd_partial_start) -
                          delta_t(cl->countryList[x].jd_partial_start))),

                        ((cl->countryList[x].jd_total_start == 0) ? 0 :
                         (unix_from_jd(cl->countryList[x].jd_total_start) -
                          delta_t(cl->countryList[x].jd_total_start))),

                        ((cl->countryList[x].jd_total_end == 0) ? 0 :
                         (unix_from_jd(cl->countryList[x].jd_total_end) -
                          delta_t(cl->countryList[x].jd_total_end))),

                        ((cl->countryList[x].jd_partial_end == 0) ? 0 :
                         (unix_from_jd(cl->countryList[x].jd_partial_end) -
                          delta_t(cl->countryList[x].jd_partial_end)))
                );
                first = 0;
            }
        fprintf(f, "}}");
        fclose(f);
    }
}

