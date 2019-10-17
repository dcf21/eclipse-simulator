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
#include <unistd.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include <gd.h>

#include "coreUtils/asciiDouble.h"
#include "coreUtils/strConstants.h"
#include "coreUtils/errorReport.h"

#include "mathsTools/julianDate.h"

#include "country_lookup.h"
#include "jpeg_in.h"
#include "settings.h"
#include "shadow_calc.h"

country_lookup_handle *country_lookup_init() {
    int j, x, y;
    FILE *f;

    // Create country lookup handle
    country_lookup_handle *cl = (country_lookup_handle *) malloc(sizeof(country_lookup_handle));

    // Read world map image
    f = fopen("/home/dcf21/projects/website4/auto/tmp/worldMap/worldMap.png", "rb");
    gdImagePtr worldMapImg = gdImageCreateFromPng(f);
    fclose(f);

    // Read information about countries
    for (j = 0; j < COUNTRYLISTLEN; j++) {
        cl->countryList[j].red = cl->countryList[j].grn = cl->countryList[j].blu = 0;
        cl->countryList[j].maxeclipse = 0;
    }
    f = fopen("/home/dcf21/projects/website4/auto/tmp/worldMap/countryList.dat", "r");
    if (f == NULL) {
        ephem_fatal(__FILE__, __LINE__, "Could not open list of country info.");
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
            ephem_fatal(__FILE__, __LINE__, "Illegal country UID.");
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
    cl->worldMapWidth = worldMapImg->sx;
    cl->worldMapHeight = worldMapImg->sy;
    cl->worldMapArray = malloc(cl->worldMapWidth * cl->worldMapHeight * sizeof(int));

    for (y = 0; y < cl->worldMapHeight; y++)
        for (x = 0; x < cl->worldMapWidth; x++) {
            int color = gdImageGetPixel(worldMapImg, x, y);
            unsigned char red = gdImageRed(worldMapImg, color);
            unsigned char grn = gdImageGreen(worldMapImg, color);
            unsigned char blu = gdImageBlue(worldMapImg, color);
            cl->worldMapArray[y * cl->worldMapWidth + x] = 0;
            for (j = 1000; j < COUNTRYLISTLEN; j++)
                if (red || grn || blu)
                    if ((cl->countryList[j].red == red) && (cl->countryList[j].grn == grn) &&
                        (cl->countryList[j].blu == blu)) {
                        cl->worldMapArray[y * cl->worldMapWidth + x] = j;
                        break;
                    }
        }
    gdImageDestroy(worldMapImg);

    return cl;
}

void country_lookup_free(country_lookup_handle *cl) {
    free(cl->worldMapArray);
    free(cl);
}

int test_if_land_or_sea(const country_lookup_handle *cl, double lng, double lat) {

    while (lng > 180) lng -= 360;
    int p0 = (int) ((lng + 180.) * cl->worldMapWidth / 360.);
    int p1 = (int) ((90. - lat) * cl->worldMapHeight / 180.);
    int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];

    return country > 0;
}

void country_lookup_max_eclipse(country_lookup_handle *cl, const settings *config, const shadow_map *greatest_shadow) {
    int x, y;
    FILE *f;

    // Work out the maximum extent of the eclipse by country
    for (x = 0; x < COUNTRYLISTLEN; x++)
        if ((cl->countryList[x].red) || (cl->countryList[x].grn) || (cl->countryList[x].blu)) {
            double lat = cl->countryList[x].lat;
            double lng = cl->countryList[x].lng;
            int p0 = (int) ((lng + 180.) * config->x_size_2d / 360.);
            int p1 = (int) ((90. - lat) * config->y_size_2d / 180.);
            if (p0 >= config->x_size_2d) p0 -= config->x_size_2d;
            if (p0 < 0) p0 += config->x_size_2d;
            if (p1 >= config->y_size_2d) p1 = config->y_size_2d - 1;
            double shadow = greatest_shadow->map[p0 + p1 * config->x_size_2d];
            if (shadow > 1) {
                if (DEBUG) {
                    sprintf(temp_err_string,
                            "Error. Home position of country <%d> at (%.1f,%.1f) pixel (%d,%d) with shadow %.1f.", x,
                            lat, lng, p0, p1, shadow);
                    ephem_log(temp_err_string);
                }
            }
            cl->countryList[x].maxeclipse = shadow;
        }

    for (y = 0; y < config->y_size_2d; y++)
        for (x = 0; x < config->x_size_2d; x++) {
            double shadow = greatest_shadow->map[x + y * config->x_size_2d];
            if (shadow > 0) {
                int p0 = ((double) x) / config->x_size_2d * cl->worldMapWidth;
                int p1 = ((double) y) / config->y_size_2d * cl->worldMapHeight;
                int country = cl->worldMapArray[p1 * cl->worldMapWidth + p0];
                if ((country > 0) && (cl->countryList[country].maxeclipse < shadow))
                    cl->countryList[country].maxeclipse = shadow;
            }
        }

    char fname[FNAME_LENGTH];
    sprintf(fname, "%s/maximumEclipse.json", config->output_dir);
    f = fopen(fname, "w");
    fprintf(f, "{\"byuid\":{");
    y = 0;
    for (x = 1000; x < 2000; x++)
        if ((cl->countryList[x].red) || (cl->countryList[x].grn) || (cl->countryList[x].blu))
            if (cl->countryList[x].maxeclipse > 0.01) {
                if (y) fprintf(f, ",");
                fprintf(f, "\"%d\":%d", x, (int) (cl->countryList[x].maxeclipse * 100));
                y = 1;
            }
    fprintf(f, "},\"extras\":{");
    y = 0;
    for (x = 2000; x < COUNTRYLISTLEN; x++)
        if (cl->countryList[x].maxeclipse > 0.01) {
            if (y) fprintf(f, ",");
            fprintf(f, "\"%s\":%d", cl->countryList[x].name, (int) (cl->countryList[x].maxeclipse * 100));
            y = 1;
        }
    fprintf(f, "}}");
    fclose(f);
}

