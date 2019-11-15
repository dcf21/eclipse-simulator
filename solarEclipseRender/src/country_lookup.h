// country_lookup.h

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

#ifndef COUNTRY_LOOKUP_H
#define COUNTRY_LOOKUP_H 1

#include "settings.h"
#include "shadow_calc.h"
#include "map_greatest_eclipse.h"

typedef struct countryInfo {
    unsigned char red, grn, blu;
    double lat, lng, max_eclipse;
    char name[64];
} countryInfo;

#define COUNTRYLISTLEN 2100

typedef struct country_lookup_handle {
    countryInfo countryList[COUNTRYLISTLEN];
    int worldMapWidth;
    int worldMapHeight;
    int *worldMapArray;
} country_lookup_handle;

country_lookup_handle *country_lookup_init();

void country_lookup_free(country_lookup_handle *cl);

int test_if_land_or_sea(const country_lookup_handle *cl, double longitude, double latitude);

void country_lookup_max_eclipse(country_lookup_handle *cl, const settings *config, const shadow_map *greatest_shadow,
                                const eclipse_path_list *eclipse_path);

#endif
