// shadow_calc.h

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

#ifndef SHADOW_CALC_H
#define SHADOW_CALC_H 1

typedef struct shadow_map {
    double *map;
    double *lat, *lng;  // degrees
} shadow_map;

double siderealTime(double JD);

void earthTopocentricPositionICRF(double *out, double lat, double lng, double radius_in_earth_radii,
                                  const double *pos_earth, double epoch, double sidereal_time);

int testIfAnnularEclipse(double lat, double lng, double JD,
                         double radius, const double *pos_sun, const double *pos_moon, const double *pos_earth,
                         double sidereal_time);

double getShadowFraction(double lat, double lng, double JD,
                         double radius, const double *pos_sun, const double *pos_moon, const double *pos_earth,
                         double sidereal_time);

shadow_map *allocate_shadow_map(int x_size, int y_size);

void shadow_map_free(shadow_map *item);

void calculate_where_sun_overhead(double *lat_sun, double *lng_sun, double *sidereal_time,
                                  const double *pos_sun, const double *pos_earth, double jd);

shadow_map *calculate_eclipse_map_2d(const settings *config,
                                     double jd, const double *pos_sun, const double *pos_earth,
                                     const double *pos_moon,
                                     time_span *span_output, shadow_map *greatest_shadow);

shadow_map *calculate_eclipse_map_3d(const settings *config,
                                     double jd, const double *pos_sun, const double *pos_earth,
                                     const double *pos_moon);

#endif
