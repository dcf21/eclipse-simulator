// country_lookup.c

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "coreUtils/asciiDouble.h"
#include "coreUtils/strConstants.h"
#include "coreUtils/errorReport.h"

#include "png/image.h"

#include "country_lookup.h"
#include "settings.h"
#include "shadow_calc.h"
#include "map_greatest_eclipse.h"

country_lookup_handle *country_lookup_init() {
    int j, x, y;
    FILE *f;

    // Create country lookup handle
    country_lookup_handle *cl = (country_lookup_handle *) malloc(sizeof(country_lookup_handle));

    // Read world map image
    image_ptr worldMapImg = image_get(SRCDIR "../worldMap/worldMap.png");

    // Read information about countries
    for (j = 0; j < COUNTRYLISTLEN; j++) {
        cl->countryList[j].red = cl->countryList[j].grn = cl->countryList[j].blu = 0;
        cl->countryList[j].max_eclipse = 0;
    }
    f = fopen(SRCDIR "../worldMap/countryList.dat", "r");
    if (f == NULL) {
        logging_fatal(__FILE__, __LINE__, "Could not open list of country info.");
        exit(1);
    }
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
        lineptr = next_word(lineptr);
        cl->countryList[uid].red = (int) get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].grn = (int) get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].blu = (int) get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].lat = get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        cl->countryList[uid].lng = get_float(lineptr, NULL);
        lineptr = next_word(lineptr);
        strcpy(cl->countryList[uid].name, lineptr);
    }
    fclose(f);

    // Convert world map image to an array of ints
    cl->worldMapWidth = worldMapImg.xsize;
    cl->worldMapHeight = worldMapImg.ysize;
    cl->worldMapArray = malloc(cl->worldMapWidth * cl->worldMapHeight * sizeof(int));

    for (y = 0; y < cl->worldMapHeight; y++)
        for (x = 0; x < cl->worldMapWidth; x++) {
            const int offset = x + y * worldMapImg.xsize;
            const unsigned char red = (unsigned char)worldMapImg.data_red[offset];
            const unsigned char grn = (unsigned char)worldMapImg.data_grn[offset];
            const unsigned char blu = (unsigned char)worldMapImg.data_blu[offset];
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

    return cl;
}

void country_lookup_free(country_lookup_handle *cl) {
    free(cl->worldMapArray);
    free(cl);
}

int test_if_land_or_sea(const country_lookup_handle *cl, double lng, double lat) {
    while (lng > 180) lng -= 360;
    const int p0 = (int) ((lng + 180.) * cl->worldMapWidth / 360.);
    const int p1 = (int) ((90. - lat) * cl->worldMapHeight / 180.);
    const int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];

    return country > 0;
}

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
            double shadow = greatest_shadow->map[p0 + p1 * config->x_size_2d];
            cl->countryList[x].max_eclipse = shadow;
        }

    // Loop over the pixels in the 2D bitmap map of the greatest extent of the eclipse, noting the greatest value
    // within each country
    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            double shadow = greatest_shadow->map[x + y * config->x_size_2d];
            if (shadow > 0) {
                const int p0 = (int) (((double) x) / config->x_size_2d * cl->worldMapWidth);
                const int p1 = (int) (((double) y) / config->y_size_2d * cl->worldMapHeight);
                const int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];
                if ((country > 0) && (cl->countryList[country].max_eclipse < shadow))
                    cl->countryList[country].max_eclipse = shadow;
            }
        }

    // Now trace along the path of the central eclipse
    // Shadow fraction of 101% means an annular eclipse
    // Shadow fraction of 102% means a total eclipse
    for (int i = 0; i < eclipse_path->path_count; i++) {
        const int is_total = eclipse_path->paths[i].is_total;
        for (int j = 0; j < eclipse_path->paths[i].point_count; j++) {
            // Look up the coordinates of this point
            const double longitude = eclipse_path->paths[i].path[j].longitude * 180 / M_PI;
            const double latitude = eclipse_path->paths[i].path[j].latitude * 180 / M_PI;

            int p0 = (int) ((longitude + 180.) * cl->worldMapWidth / 360.);
            int p1 = (int) ((90. - latitude) * cl->worldMapHeight / 180.);

            if (p0 >= cl->worldMapWidth) p0 -= cl->worldMapWidth;
            if (p0 < 0) p0 += cl->worldMapWidth;
            if (p1 >= cl->worldMapHeight) p1 = cl->worldMapHeight - 1;

            const int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];

            double shadow_new = is_total ? 1.02 : 1.01;
            if (shadow_new > cl->countryList[country].max_eclipse) {
                cl->countryList[country].max_eclipse = shadow_new;
            }
        }
    }

    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/maximumEclipse.json", config->output_dir);
    f = fopen(fname, "w");
    fprintf(f, "{\"byuid\":{");
    y = 0;
    for (x = 1000; x < 2000; x++)
        if ((cl->countryList[x].red) || (cl->countryList[x].grn) || (cl->countryList[x].blu))
            if (cl->countryList[x].max_eclipse > 0.01) {
                if (y) fprintf(f, ",");
                fprintf(f, "\"%d\":%d", x, (int) (cl->countryList[x].max_eclipse * 100));
                y = 1;
            }
    fprintf(f, "},\"extras\":{");
    y = 0;
    for (x = 2000; x < COUNTRYLISTLEN; x++)
        if (cl->countryList[x].max_eclipse > 0.01) {
            if (y) fprintf(f, ",");
            fprintf(f, "\"%s\":%d", cl->countryList[x].name, (int) (cl->countryList[x].max_eclipse * 100));
            y = 1;
        }
    fprintf(f, "}}");
    fclose(f);
}

