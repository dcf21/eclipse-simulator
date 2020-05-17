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

CWD=$(shell pwd)

VERSION = 2.0
DATE    = 15/09/2019
PATHLINK= /

WARNINGS= -Wall -Wno-format-truncation -Wno-unused-result
COMPILE = $(CC) $(WARNINGS) -g -c -I $(CWD)/src `mysql_config --cflags`
LIBS    = -lcairo -lgsl -lgslcblas -lpng -lz -ljpeg -lm
LINK    = $(CC) $(WARNINGS) -g

OPTIMISATION = -O3 -march=native

DEBUG   = -D DEBUG=1 -D MEMDEBUG1=1 -D MEMDEBUG2=0 -D NORISESET=0
NODEBUG = -D DEBUG=0 -D MEMDEBUG1=0 -D MEMDEBUG2=0 -D NORISESET=0

LOCAL_SRCDIR = src
LOCAL_OBJDIR = obj
LOCAL_BINDIR = bin

CORE_FILES = argparse/argparse.c coreUtils/asciiDouble.c coreUtils/errorReport.c mathsTools/julianDate.c mathsTools/sphericalAst.c country_lookup.c duration.c ephemeris.c jpeg/jpeg_in.c make_binary_map.c map_eclipse_contours.c map_greatest_eclipse.c png/image_in.c projection.c rendering.c render_2d.c render_3d.c settings.c shadow_calc.c main.c

CORE_HEADERS = argparse/argparse.h coreUtils/asciiDouble.h coreUtils/errorReport.h coreUtils/strConstants.h mathsTools/julianDate.h mathsTools/sphericalAst.h constants.h country_lookup.h duration.h ephemeris.h jpeg/jpeg.h make_binary_map.h map_eclipse_contours.h map_greatest_eclipse.h png/image.h projection.h rendering.h render_2d.h render_3d.h settings.h shadow_calc.h

CORE_SOURCES          = $(CORE_FILES:%.c=$(LOCAL_SRCDIR)/%.c)
CORE_OBJECTS          = $(CORE_FILES:%.c=$(LOCAL_OBJDIR)/%.o)
CORE_OBJECTS_DEBUG    = $(CORE_OBJECTS:%.o=%.debug.o)
CORE_HFILES           = $(CORE_HEADERS:%.h=$(LOCAL_SRCDIR)/%.h) Makefile

ALL_HFILES = $(CORE_HFILES)

SWITCHES = -D DCFVERSION=\"$(VERSION)\"  -D DATE=\"$(DATE)\"  -D PATHLINK=\"$(PATHLINK)\"  -D SRCDIR=\"$(CWD)/$(LOCAL_SRCDIR)/\"

all: $(LOCAL_BINDIR)/eclipseRender.bin $(LOCAL_BINDIR)/debug/eclipseRender.bin

#
# General macros for the compile steps
#

$(LOCAL_OBJDIR)/%.o:         $(LOCAL_SRCDIR)/%.c $(ALL_HFILES)
	mkdir -p $(LOCAL_OBJDIR) $(LOCAL_OBJDIR)/argparse $(LOCAL_OBJDIR)/coreUtils $(LOCAL_OBJDIR)/jpeg $(LOCAL_OBJDIR)/mathsTools $(LOCAL_OBJDIR)/png
	$(COMPILE) $(OPTIMISATION) $(NODEBUG) $(SWITCHES) $< -o $@

$(LOCAL_OBJDIR)/%.debug.o:   $(LOCAL_SRCDIR)/%.c $(ALL_HFILES)
	mkdir -p $(LOCAL_OBJDIR) $(LOCAL_OBJDIR)/argparse $(LOCAL_OBJDIR)/coreUtils $(LOCAL_OBJDIR)/jpeg $(LOCAL_OBJDIR)/mathsTools $(LOCAL_OBJDIR)/png
	$(COMPILE) $(OPTIMISATION) $(DEBUG)   $(SWITCHES) $< -o $@

#
# Make the binaries
#

$(LOCAL_BINDIR)/eclipseRender.bin: $(CORE_OBJECTS)
	mkdir -p $(LOCAL_BINDIR)
	$(LINK) $(OPTIMISATION) $(CORE_OBJECTS) $(LIBS) -o $(LOCAL_BINDIR)/eclipseRender.bin

$(LOCAL_BINDIR)/debug/eclipseRender.bin: $(CORE_OBJECTS_DEBUG)
	mkdir -p $(LOCAL_BINDIR)/debug
	echo "The files in this directory are binaries with debugging options enabled: they produce activity logs called 'ephem.log'. It should be noted that these binaries can up to ten times slower than non-debugging versions." > $(LOCAL_BINDIR)/debug/README
	$(LINK) $(OPTIMISATION) $(CORE_OBJECTS_DEBUG) $(LIBS) -o $(LOCAL_BINDIR)/debug/eclipseRender.bin

#
# Clean macros
#

clean:
	rm -vfR $(LOCAL_OBJDIR) $(LOCAL_BINDIR)

afresh: clean all

