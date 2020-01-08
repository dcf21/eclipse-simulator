// projection.c

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

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "projection.h"
#include "settings.h"


void project_3d(const settings *config, int x, int y, double lng_sun, double lat_sun,
                double *lng_out, double *lat_out, double *p_radius_out) {
    double zn, p_radius;

    // Project point in 3D space where we measure eclipse fraction, measured in Earth radii with z-axis
    // pointing out of the screen

    // x-axis starts out point horizontally across the screen
    double xn = -(x - config->x_size_3d / 2.) / config->earth_pixel_radius;

    // y-axis starts out pointing vertically up the screen
    double yn = -(y - config->y_size_3d / 2.) / config->earth_pixel_radius;
    double n = gsl_pow_2(xn) + gsl_pow_2(yn);

    // z-axis points out of the screen
    if (n >= 1) {
        // This pixel is not on the Earth's surface, so project point above Earth
        zn = 0;
        p_radius = sqrt(n);

        // Re-normalise (xn,yn,zn) is have unit length
        xn /= p_radius;
        yn /= p_radius;
    } else {
        // This pixel is on the Earth's surface
        zn = -sqrt(1 - n);
        p_radius = 1;
    }

    // RA - hour angle
    const double rotang1 = -M_PI / 2 + lng_sun;
    const double rotang2 = -M_PI / 2 - lat_sun;

    // Tip (y,z) by declination
    const double x2 = xn;
    const double y2 = yn * cos(rotang2) + zn * sin(rotang2);
    const double z2 = -yn * sin(rotang2) + zn * cos(rotang2);

    // Tip (x,y) to get the right central longitude
    const double xf = x2 * cos(rotang1) - y2 * sin(rotang1);
    const double yf = x2 * sin(rotang1) + y2 * cos(rotang1);
    const double zf = z2;

    // Coordinates now orientated with z-axis pointing out of north pole
    // longitude
    *lng_out = atan2(yf, xf) / M_PI * 180;
    *lat_out = asin(zf) / M_PI * 180;
    *p_radius_out = p_radius;
}

void inv_project_3d(const settings *config, int *x_out, int *y_out, double lng_sun, double lat_sun,
                    double lng, double lat) {

    lng *= M_PI / 180;
    lat *= M_PI / 180;

    // RA - hour angle
    const double rotang1 = -M_PI / 2 + lng_sun;
    const double rotang2 = -M_PI / 2 - lat_sun;

    const double xf = cos(lng) * cos(lat);
    const double yf = sin(lng) * cos(lat);
    const double zf = sin(lat);

    // Tip (x,y) to get the right central longitude
    const double x2 = xf * cos(rotang1) + yf * sin(rotang1);
    const double y2 = -xf * sin(rotang1) + yf * cos(rotang1);
    const double z2 = zf;

    // Tip (y,z) by declination
    const double xn = x2;
    const double yn = y2 * cos(rotang2) - z2 * sin(rotang2);
    //const double zn = y2 * sin(rotang2) + z2 * cos(rotang2);

    *x_out = (int) (config->x_size_3d / 2. - xn * config->earth_pixel_radius);
    *y_out = (int) (config->y_size_3d / 2. - yn * config->earth_pixel_radius);
}