#!/usr/bin/python3
# -*- coding: utf-8 -*-
# solarEclipses_makeKml.py
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
Make KML files describing the paths of solar eclipses, from JSON equivalents
"""

import glob
import json
import logging
import os
import sys

src_path = os.getcwd()
out_path = os.path.join(src_path, "output")


class KMLTemplate:
    """
    Template to use for producing KML files describing the paths of eclipses
    """

    kml_document = """\
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <Style id="eclipseOuter">
      <PolyStyle>
        <color>30000000</color>
        <fill>1</fill>
        <outline>1</outline>
      </PolyStyle>
    </Style>
    <Style id="eclipseInner">
      <PolyStyle>
        <color>18000000</color>
        <fill>1</fill>
        <outline>1</outline>
      </PolyStyle>
    </Style>
    <Style id="eclipseTotal">
      <LineStyle>
        <color>ff0000ff</color>
        <width>2</width>
      </LineStyle>
    </Style>
    <Style id="eclipseAnnular">
      <LineStyle>
        <color>ff0ff0000</color>
        <width>2</width>
      </LineStyle>
    </Style>
    {body}
  </Document>
</kml>
"""

    kml_path = """
    <Placemark>
      <name>{path_title}</name>
      <styleUrl>#{path_style}</styleUrl>
      <LineString>
        <coordinates>
          {point_string}
        </coordinates>
      </LineString>
    </Placemark>
"""

    kml_polygon = """
    <Placemark>
      <name>{polygon_title}</name>
      <styleUrl>#{polygon_style}</styleUrl>
      <Polygon>
        <outerBoundaryIs>
          <LinearRing>
            <coordinates>
              {point_string}
            </coordinates>
          </LinearRing>
        </outerBoundaryIs>
      </Polygon>
    </Placemark>
"""


def approx(flt, n=2):
    """
    Round a floating point number to a set number of decimal places in order to make more compact JSON output.

    :param flt:
        Input floating-point value
    :param n:
        Number of decimal places to round to
    :return:
        Rounded floating point value
    """
    divisor = float(pow(10, n))
    if (type(flt) != int) and (type(flt) != float):
        print("Warning: approx() passed argument of type <{}>".format(type(flt)))
        return 0
    return int(flt * divisor) / divisor


def generate_point_string(point_list):
    """
    Turn a list of (longitude, latitude) points, each in degrees, into a text string for inclusion in a KML file.

    KML files don't like paths switching from longitude +180 to -180 and vice versa, so we smooth over these transitions
    by instead going out of the range +/-180 degrees.

    :param point_list:
        List of (longitude, latitude) points, each in degrees
    :return:
        Textual representation
    """
    previous_longitude = None
    point_string = ""

    for point in point_list:
        longitude, latitude = point

        # Check whether longitude has flipped across the discontinuity at +/- 180 degrees
        if previous_longitude is not None:
            while longitude <= previous_longitude - 180:
                longitude += 360
            while longitude > previous_longitude + 180:
                longitude -= 360
        previous_longitude = longitude

        # Add point to KML file
        point_string += "{0:.4f},{1:.4f}\n".format(longitude, latitude)
    return point_string


def make_kml():
    """
    Make KML files describing the paths of solar eclipses, from JSON equivalents.

    :return:
        None
    """

    # Read JSON files produced by <solarEclipses.py> describing the path of each solar eclipse
    json_in_path = os.path.join(out_path, "*_path.json")
    filenames = glob.glob(json_in_path)

    json_data = []
    for filename in filenames:
        item = json.loads(open(filename).read())
        item['filename'] = filename
        json_data.append(item)

    # Sort items in order of computational time, and name the worst offenders
    json_data.sort(key=lambda x: x['compute_duration'])
    json_data.reverse()
    item_count = 10
    for number, item in enumerate(json_data[:item_count]):
        logging.info("Longest compute times {n}/{c}: {k} - {t} ({d:.1f} min)".
                     format(n=number + 1, c=item_count, k=item['event_key'], t=item['event_title'],
                            d=item['compute_duration'] / 60.))

    # Sort items in order of duration error, and name the worst offenders
    json_data.sort(key=lambda x: abs(x['event_duration_ford'] - x['event_duration_espenak']))
    json_data.reverse()
    item_count = 25
    for number, item in enumerate(json_data[:item_count]):
        logging.info("Worst durations {n}/{c}: {k} - {t} ({d1:5.1f} vs {d2:5.1f})".
                     format(n=number + 1, c=item_count, k=item['event_key'], t=item['event_title'],
                            d1=item['event_duration_ford'], d2=item['event_duration_espenak']))

    # Sort items into chronological order
    json_data.sort(key=lambda x: x['filename'])

    # Compile list of paths of all total / annular eclipses
    paths_all = []
    for item in json_data:
        # Ignore partial eclipses, with no path
        if len(item['path']) == 0:
            continue

        # Calculate midpoint lat/lng
        flattened_path = [point for path in item['path'] for point in path[1]]
        midpoint = [approx(x) for x in flattened_path[int(len(flattened_path) / 2)][1:3]]

        # Shorten paths
        path_list = []
        for x in item['path']:
            points = [[approx(p[1]), approx(p[2])] for p in x[1][::5]]
            path_list.append([x[0], points])

        # Append item to list
        paths_all.append({
            'paths': path_list,
            'event_key': item['event_key'],
            'event_title': item['event_title'],
            'month': item['event_month'],
            'year': int(item['event_month'][-4:]),
            'max_duration': item['event_duration_espenak'],
            'mean_position': midpoint
        })

    # Write list of paths
    filename = os.path.join(out_path, "paths_all.json")
    with open(filename, "w") as f:
        f.write(json.dumps(paths_all, separators=(',', ':')))

    # Produce KML files
    kml_template = KMLTemplate()
    for item in json_data:
        kml_body = ""
        # Add paths where eclipse is total or annular
        for path in item['path']:
            point_string = generate_point_string(point_list=[point[1:3] for point in path[1]])
            kml_body += kml_template.kml_path.format(path_title="Central path, {}".format(item['event_month']),
                                                     path_style="eclipseTotal" if path[0] else "eclipseAnnular",
                                                     point_string=point_string)

        # Add contours
        for contour in item['contours']:
            if len(contour[1]) > 0:
                point_string = generate_point_string(point_list=[point[0:2] for point in contour[1]])
                kml_body += kml_template.kml_polygon.format(polygon_title=contour[0],
                                                            polygon_style=(
                                                                "eclipseOuter" if contour[0] == "Partial eclipse"
                                                                else "eclipseInner"),
                                                            point_string=point_string)

        # Write KML output
        filename = os.path.join(out_path, "{}.kml".format(item['filename'][:-5]))
        with open(filename, "w") as f:
            f.write(kml_template.kml_document.format(body=kml_body))


# Do it right away if we're run as a script
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        stream=sys.stdout,
                        format='[%(asctime)s] %(levelname)s:%(filename)s:%(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.info(__doc__.strip())

    make_kml()
