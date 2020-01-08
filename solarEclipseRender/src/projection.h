// projection.h

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

#ifndef PROJECTION_H
#define PROJECTION_H 1

#include "settings.h"

void project_3d(const settings *config, int x, int y, double lng_sun, double lat_sun,
                double *lng_out, double *lat_out, double *p_radius_out);

void inv_project_3d(const settings *config, int *x_out, int *y_out, double lng_sun, double lat_sun,
                    double lng, double lat);

#endif
