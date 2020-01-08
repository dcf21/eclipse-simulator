# Makefile for solar eclipse renderer
# -------------------------------------------------
# Copyright 2019-2020 Dominic Ford.

# This file is part of EclipseRender.

# EclipseRender is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EclipseRender is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with EclipseRender.  If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------

cmake_minimum_required(VERSION 3.6)
project(solarEclipseRender)

include_directories(src)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -std=c99 -DSRCDIR='\"${CMAKE_SOURCE_DIR}/src\"' -DDCFVERSION='\"x\"' -DDEBUG=0 -D MEMDEBUG1=0 -D MEMDEBUG2=0 -D NORISESET=0")

set(SOURCE_FILES
        src/argparse/argparse.c
        src/argparse/argparse.h
        src/coreUtils/asciiDouble.c
        src/coreUtils/asciiDouble.h
        src/coreUtils/errorReport.c
        src/coreUtils/errorReport.h
        src/coreUtils/strConstants.h
        src/mathsTools/julianDate.c
        src/mathsTools/julianDate.h
        src/mathsTools/sphericalAst.c
        src/mathsTools/sphericalAst.h
        src/constants.h
        src/country_lookup.c
        src/country_lookup.h
        src/duration.c
        src/duration.h
        src/ephemeris.c
        src/ephemeris.h
        src/jpeg/jpeg.h
        src/jpeg/jpeg_in.c
        src/make_binary_map.c
        src/make_binary_map.h
        src/map_eclipse_contours.c
        src/map_eclipse_contours.h
        src/map_greatest_eclipse.c
        src/map_greatest_eclipse.h
        src/png/image.h
        src/png/image_in.c
        src/projection.c
        src/projection.h
        src/rendering.c
        src/rendering.h
        src/render_2d.c
        src/render_2d.h
        src/render_3d.c
        src/render_3d.h
        src/settings.c
        src/settings.c
        src/shadow_calc.c
        src/shadow_calc.h
        src/main.c)

add_executable(solarEclipseRender ${SOURCE_FILES})

target_link_libraries(solarEclipseRender mysqlclient gsl gslcblas cairo png z jpeg m)
