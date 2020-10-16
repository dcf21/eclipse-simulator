# eclipse-simulator

Eclipse simulator is a command-line tool for producing animations of solar
eclipses. It was written to produce all of the diagrams, KML path files, and
video animations that appear on <https://in-the-sky.org/>. An example of the
kinds of media it produces can be seen on this page:
<https://in-the-sky.org/news.php?id=20201214_09_100>

The simulations are based on the DE430 planetary ephemeris computed by the Jet
Propulsion Laboratory (JPL). The position of the Sun, Earth and Moon are
extracted from the DE430 files using EphemerisCompute
(<https://github.com/dcf21/ephemeris-compute-de430>), written by the author and
also freely available for download.

They assume that the Earth and Moon are both ellipsoids with fixed polar and
equatorial radii, and do not take into account the irregular topography of
either body. All eclipse predictions are made at sea level.  In practice, this
means that the predictions presented here are inaccurate by at most of few
seconds.

### Supported operating systems

`eclipse-simulator` is written in C and runs in Linux, MacOS, and other
Unix-like operating systems. A front-end user interface is written in python3.
The simulator does not run under Windows.

### License

This code is distributed under the Gnu General Public License. It is (C)
Dominic Ford 2012 - 2020.

### Set up

Before you start, the simulator needs to download various data from the
internet, including the tool `ephemerisCompute` from the author's GitHub
repository, and the DE430 ephemeris files.

This can be done with the Python script `solarEclipses.py`. The total download
size will be around 500 MB. These files are only downloaded the first time this
script is run.

The script will then proceed to run some demo simulations; by default it
simulates the eclipses of 2020 and 2021. Each simulation takes around 45
minutes on a single core of a modern computer (e.g. Intel i7 8700), or around
120 minutes on a machine dating from around 2010. If your computer has multiple
cores, multiple simulations will run on parallel on the various cores.

### Simulating particular eclipses

`eclipse-simulator` may be used to simulate any solar eclipse within the time
span of the DE430 planetary ephemeris - i.e. 1600 to 2200 AD. To simulate a
particular eclipse, pass the command-line switches `--year-min` and
`--year-max` to the Python script `solarEclipses.py`. For example:

```
./solarEclipses.py --year-min 2025 --year-max 2030
```

