// duration.c

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
#include <math.h>

#include <gsl/gsl_math.h>

#include "mathsTools/sphericalAst.h"

#include "duration.h"
#include "ephemeris.h"
#include "settings.h"
#include "shadow_calc.h"

/// eclipse_duration - Calculate the duration of an eclipse, as viewed from a particular place on Earth
///
/// \param config [in] - settings
/// \param ephemeris [in] - an ephemeris computed using <ephemerisCompute>
/// \param longitude [in] - The latitude at which to calculate the eclipse's duration (radians)
/// \param latitude [in] - The longitude at which to calculate the eclipse's duration (radians)
/// \param jd_midpoint [in] - The approximate midpoint of the eclipse as viewed from that location (Julian day)
/// \param is_total [in] - Boolean flag indicating whether to calculate the duration of totality or annularity
/// \return Duration, seconds
double eclipse_duration(const settings *config, const ephemeris *ephemeris,
                        double longitude, double latitude, double jd_midpoint, int is_total) {
    // search for 10 minutes before / after maximum eclipse
    const double search_span = 600; // seconds

    // start and end points of total / annular eclipse at this location
    double jd_eclipse_start = jd_midpoint;
    double jd_eclipse_end = jd_midpoint;

    // start and end points for our search process
    const double jd_search_start = jd_midpoint - (search_span / 86400);
    const double jd_search_end = jd_midpoint + (search_span / 86400);

    // first and last point numbers in the ephemeris structure that we need to analyse
    const int point_start = (int) ((jd_search_start - ephemeris->jd_start) / ephemeris->jd_step);
    const int point_end = (int) ((jd_search_end - ephemeris->jd_start) / ephemeris->jd_step);

    // loop over the ephemeris points that we are to analyse
    for (int j = point_start; j < point_end; j++) {
        const double jd = ephemeris->jd_start + ephemeris->jd_step * j;
        const ephemeris_point p = ephemeris->data[j];

        // Calculate the latitude and longitude on Earth where the Sun is overhead
        double sidereal_time, lat_sun, lng_sun; // all in radians
        calculate_where_sun_overhead(&lat_sun, &lng_sun, &sidereal_time, p.sun_pos, p.earth_pos, jd);

        // Shadow fraction is zero if the Sun is below the horizon
        const double sun_ang_dist = angDist_RADec(longitude, latitude, lng_sun, lat_sun);
        if (sun_ang_dist >= M_PI / 2) continue;

        if (is_total) {
            // For total eclipses, the duration of the eclipse is defined by when the shadow fraction is 100%
            const double shadow_fraction = getShadowFraction(latitude * 180 / M_PI, longitude * 180 / M_PI, jd, 1,
                                                             p.sun_pos, p.moon_pos, p.earth_pos,
                                                             sidereal_time);

            if (shadow_fraction == 1) {
                if (jd < jd_eclipse_start) jd_eclipse_start = jd;
                if (jd > jd_eclipse_end) jd_eclipse_end = jd;
            }
        } else {
            // For annular eclipses, look at whether the Moon is contained within the Sun's disk
            const int is_annular = testIfAnnularEclipse(latitude * 180 / M_PI, longitude * 180 / M_PI, jd, 1,
                                                        p.sun_pos, p.moon_pos, p.earth_pos,
                                                        sidereal_time);

            if (is_annular) {
                if (jd < jd_eclipse_start) jd_eclipse_start = jd;
                if (jd > jd_eclipse_end) jd_eclipse_end = jd;
            }
        }
    }
    return (jd_eclipse_end - jd_eclipse_start) * 86400; // seconds
}
