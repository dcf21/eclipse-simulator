#!/usr/bin/python3
# -*- coding: utf-8 -*-
# solarEclipses.py
#
# The python script in this file makes eclipse simulations.
#
# Copyright (C) 2012-2019 Dominic Ford <dcf21-www@dcford.org.uk>
#
# This code is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# You should have received a copy of the GNU General Public License along with
# this file; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA  02110-1301, USA

# ----------------------------------------------------------------------------

"""
Create videos of where solar eclipses are visible, and also JSON data files needed to simulate them.
"""

import argparse
import json
import logging
import os
import re
import sys
import time
from math import floor

import dask
import solarEclipses_makeKml
from dcf_ast import month_name, julian_day, unix_from_jd

pid = os.getpid()

src_path = os.getcwd()
out_path = os.path.join(src_path, "output")


@dask.delayed
def make_ephemeris(event_key, event_title, event_duration, j_day, year, month, day):
    """
    Make a diagram of the visibility of a particular eclipse.

    :param event_key:
        The eventKey for this eclipse in the table <inthesky_events>
    :param event_title:
        The title for this eclipse event in the table <inthesky_events>
    :param event_duration:
        The duration of this event (seconds) as taken from Fred Espenak
    :param j_day:
        The Julian day of the eclipse
    :type j_day:
        float 
    :param year:
        The year of the eclipse
    :type year:
        int 
    :param month:
        The month number of the eclipse (1-12)
    :type month:
        int 
    :param day:
        The day of the month at the central point of the eclipse (1-31)
    :type day:
        int 
    :return: 
        A dictionary of the start and end times of the eclipse
    """

    id1 = "{}_{}".format(pid, j_day)
    id2 = "{:04d}{:02d}{:02d}".format(year, month, day)

    # Status update
    logging.info("Working on {:02d}/{:02d}/{:04d}".format(day, month, year))

    # Create temporary working directory
    tmp = "/tmp/eclipse_{}".format(id1)
    os.system("mkdir -p {}".format(tmp))
    os.system("rm -Rf {}/*".format(tmp))

    # Construct command-line arguments to pass to eclipse-rendering tool
    eclipse_renderer = os.path.join(src_path, "solarEclipseRender/bin/eclipseRender.bin")
    jd_min = j_day - 3.2 / 24
    jd_max = j_day + 3.2 / 24

    cmd = "{exe} --jd_min {min:.3f} --jd_max {max:.3f} --title \"{title}\" --output {tmp} > {tmp}/stdout". \
        format(exe=eclipse_renderer, min=jd_min, max=jd_max, title=event_title, tmp=tmp)
    logging.info("Running command: {}".format(cmd))
    time_start = time.time()
    os.system(cmd)
    time_end = time.time()

    # Read time that eclipse starts / end
    times = open("{tmp}/stdout".format(tmp=tmp)).read().split()

    if len(times) < 7:
        logging.error("FAILURE in command: {}".format(cmd))
        raise RuntimeError("Eclipse simulator failed")

    output = {
        "event_key": event_key,
        "partial_start": unix_from_jd(float(times[0])),
        "partial_end": unix_from_jd(float(times[1])),
        "total_start": unix_from_jd(float(times[2])),
        "total_end": unix_from_jd(float(times[3])),
        "greatest_eclipse_magnitude": float(times[4]),
        "greatest_eclipse_latitude": float(times[5]),
        "greatest_eclipse_longitude": float(times[6])
    }

    # Convert frames into MP4 video files
    os.system("cd {} ; "
              "ffmpeg -nostats -loglevel panic -r 10 -i frameA%06d.png -codec:v libx264 -crf 28 {}/solar_{}_A.mp4".
              format(tmp, out_path, id2))

    os.system("cd {} ; "
              "ffmpeg -nostats -loglevel panic -r 10 -i frameB%06d.png -codec:v libx264 -crf 28 {}/solar_{}_B.mp4".
              format(tmp, out_path, id2))

    os.system("cd {} ; "
              "ffmpeg -nostats -loglevel panic -r 10 -i frameA%06d.png -acodec vorbis -vcodec libtheora -q:v 5 "
              "{}/solar_{}_A.ogg".format(tmp, out_path, id2))

    os.system("cd {} ; "
              "ffmpeg -nostats -loglevel panic -r 10 -i frameB%06d.png -acodec vorbis -vcodec libtheora -q:v 5 "
              "{}/solar_{}_B.ogg".format(tmp, out_path, id2))

    # Move images and JSON files to final destination
    os.system("cd {} ; mv solarEclipseA.png {}/solar_{}_A.png".format(tmp, out_path, id2))
    os.system("cd {} ; mv solarEclipseB.png {}/solar_{}_B.png".format(tmp, out_path, id2))
    os.system("cd {} ; mv solarEclipseC.png {}/solar_{}_C.png".format(tmp, out_path, id2))
    os.system("cd {} ; mv maximumEclipse.png {}/solar_{}.png".format(tmp, out_path, id2))
    os.system("cd {} ; mv maximumEclipse.pdf {}/solar_{}.pdf".format(tmp, out_path, id2))
    os.system("cd {} ; mv maximumEclipse.svg {}/solar_{}.svg".format(tmp, out_path, id2))
    os.system("cd {} ; mv maximumEclipse.json {}/solar_{}.json".format(tmp, out_path, id2))
    os.system("cd {} ; mv maximumEclipse.dat {}/solar_{}.dat".format(tmp, out_path, id2))

    # Merge JSON files
    path_json = json.loads(open("{}/maximumEclipsePath.json".format(tmp)).read())
    contours_json = json.loads(open("{}/maximumEclipseContours.json".format(tmp)).read())

    with open("{}/solar_{}_path.json".format(out_path, id2), "w") as f:
        f.write(json.dumps({
            'path': path_json['paths'],
            'contours': contours_json,
            'compute_time': int(time_start),
            'compute_duration': int(time_end - time_start),
            'event_title': event_title,
            'event_month': "{d} {m} {y}".format(d=day, m=month_name[month - 1], y=year),
            'event_duration_espenak': event_duration,  # seconds
            'event_duration_ford': path_json['duration'],  # seconds
            'path_segments': path_json['path_segments'],
            'path_midpoint_latitude': path_json['midpoint_latitude'],
            'path_midpoint_longitude': path_json['midpoint_longitude'],
            "event_key": event_key,
            "partial_start": int(output['partial_start']),
            "partial_end": int(output['partial_end']),
            "total_start": int(output['total_start']),
            "total_end": int(output['total_end']),
            "greatest_eclipse_magnitude": output['greatest_eclipse_magnitude'],
            "greatest_eclipse_latitude": output['greatest_eclipse_latitude'],
            "greatest_eclipse_longitude": output['greatest_eclipse_longitude'],
            'host': os.uname()[1]
        }, separators=(',', ':')))

    # Clean up temporary files
    os.system("rm -Rf {}".format(tmp))

    return output


def solar_eclipses(year_min, year_max):
    """
    Create videos of where solar eclipses are visible, and also JSON data files needed to simulate them.

    :return:
        None
    """

    # Create target directory
    os.system("mkdir -p {}".format(out_path))
    os.system("rm -f {}/solar*".format(out_path))

    # Yuck, but this has to be done
    os.environ['TZ'] = 'UTC'
    time.tzset()

    # Compile renderer
    os.system("cd {} ; ./prettymake".format(os.path.join(src_path, "solarEclipseRender")))

    # Look through list of solar eclipses, and pick out events within specified range of years. Then make animations
    # of all these eclipses
    data = "solarEclipses.dat"

    # Loop over all events
    eclipse_times_delayed = []
    for line in open(data):
        line = line.strip()
        if (not line) or (line[0] == "#"):
            continue
        words = line.split()
        ev_year = float(words[2])
        if (ev_year < year_min) or (ev_year > year_max):
            continue
        ev_mc = words[3]
        ev_mon = month_name.index(ev_mc)
        ev_day = float(words[4])
        ev_hms = words[5]
        [hours_str, mins_str, secs_str] = ev_hms.split(":")
        ev_hours = float(hours_str)
        ev_mins = float(mins_str)
        # ev_secs = float(secs_str)

        # Round to start on a round multiple of 10-minutes
        ev_mins = floor(ev_mins / 10) * 10
        ev_secs = 0

        jd0 = julian_day(year=int(ev_year), month=int(ev_mon + 1), day=int(ev_day),
                         hour=int(ev_hours), minute=int(ev_mins), sec=float(ev_secs))

        # Make string describing eclipse type
        eclipse_types = {
            'A': 'Annular solar eclipse',
            'H': 'Hybrid solar eclipse',
            'P': 'Partial solar eclipse',
            'T': 'Total solar eclipse'
        }

        eclipse_type = eclipse_types[words[9][0]]

        # Look up eclipse duration
        eclipse_duration = 0
        test = re.match(r"(\d\d)m(\d\d)s", words[-1])
        if test:
            eclipse_duration = int(test.group(1)) + int(test.group(2)) / 60  # minutes

        # Make diagrams of eclipse
        eclipse_times_delayed.append(make_ephemeris(event_key="eclipse",
                                                    event_title=eclipse_type,
                                                    event_duration=eclipse_duration * 60,
                                                    j_day=jd0,
                                                    year=int(ev_year),
                                                    month=int(ev_mon + 1),
                                                    day=int(ev_day)))

    # Run eclipse simulations in parallel
    dask.compute(*eclipse_times_delayed)


# Do it right away if we're run as a script
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        stream=sys.stdout,
                        format='[%(asctime)s] %(levelname)s:%(filename)s:%(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.info(__doc__.strip())

    # Read command-line arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--year-min', dest='year_min', type=int, default=2019,
                        help="The earliest year for which to compute eclipse simulations.")
    parser.add_argument('--year-max', dest='year_max', type=int, default=2019,
                        help="The latest year for which to compute eclipse simulations.")
    args = parser.parse_args()

    # We require the ephemerisCompute code. Make sure we have this now.
    # We run a dummy ephemeris, in order to make sure that binaries versions of text files are generated now.
    ephemeris_compute_path = os.path.join(src_path, "ephemeris-compute")
    if not os.path.exists(ephemeris_compute_path):
        logging.info("Installing dependant code from <https://github.com/dcf21/ephemeris-compute.git>")
        os.system("git clone https://github.com/dcf21/ephemeris-compute.git")
        os.system("cd ephemeris-compute ; ./setup.sh ; cd bin ; ./ephem.bin")

    # Make eclipse simulations
    solar_eclipses(year_min=args.year_min, year_max=args.year_max)
    solarEclipses_makeKml.make_kml()
