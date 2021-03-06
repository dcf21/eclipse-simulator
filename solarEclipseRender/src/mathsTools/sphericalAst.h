// sphericalAst.h
// 
// -------------------------------------------------
// Copyright 2015-2020 Dominic Ford
//
// This file is part of EclipseRender.
//
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

#ifndef SPHERICALAST_H
#define SPHERICALAST_H 1

double angDist_ABC(double xa, double ya, double za, double xb, double yb, double zb, double xc, double yc, double zc);

double angDist_RADec(double ra0, double dec0, double ra1, double dec1);

void ra_dec_from_j2000(double ra0, double dec0, double utc_new, double *ra_out, double *dec_out);

void ra_dec_to_j2000(double ra1, double dec1, double utc_old, double *ra_out, double *dec_out);

void find_mean_position(double lng0, double lat0, double weight0,
                        double lng1, double lat1, double weight1,
                        double *lng_mean, double *lat_mean);

#endif

