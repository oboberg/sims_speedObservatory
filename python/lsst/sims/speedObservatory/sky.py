from __future__ import division

import numpy as np
import ephem
from datetime import datetime
from datetime import timedelta
import numbers
import time

import utils
from telescope import Telescope

from matplotlib import pyplot as plt
import palpy

__all__ = ["nightLength", "nightStart", "nightEnd", "nightNum", 
           "raOfMeridian", "phaseOfMoon", "getExpTime", "unix2lst", 
           "radec2altaz", "altaz2radec"]

# set up an observer to calculate sun rising/setting times
tel = ephem.Observer()

# this uses the default longitude and latitude values for the telescope.
# The alternative is to have methods take a Telescope object, which would
# allow lon/lat to be configurable, but the fewer arguments the better,
# and it seems unlikely that we'll need to try out different lon/lats
tel.lat = Telescope.latitude
tel.lon = Telescope.longitude
# start observing when sun is 12 deg below horizon
tel.horizon = "-12:00:00"
sun = ephem.Sun()

# cache night starts and ends since ephem.next_rising()/setting() are slow
nightStarts = {}
nightEnds = {}


def nightLength(surveyStartTime, nightNum):
    """Finds the length of the `nightNum`th night after `surveyStartTime`

    Parameters
    ----------
    surveyStartTime : float
        The time that the survey starts as a unix timestamp.
    nightNum : int
        The number (0 indexed) of nights after the start of the survey
        to get the night length of.

    Returns
    -------
    The length of the night (in seconds) that starts `nightNum` days after
    the day during which surveyStartTime falls.
    """
    return nightEnd(surveyStartTime, nightNum) - nightStart(surveyStartTime, nightNum)

def nightStart(surveyStartTime, nightNum):
    """Gives the start time of the `nightNum`th night

    Parameters
    ----------
    surveyStartTime : float
        The time that the survey starts as a unix timestamp.
    nightNum : int
        The number (0 indexed) of nights after the start of the survey
        whose start time is desired. Note that the 0th night starts in the
        evening of whatever day `surveyStartTime` falls on.

    Returns
    -------
    The start time (as a unix timestamp) of the `nightNum`th night
    """

    # use this (and the corresponding line in nightEnd) to return
    # constant-length nights (good for debugging)
    #return surveyStartTime + nightNum * 3600*24

    # check if we already know when this night starts
    if (surveyStartTime, nightNum) in nightStarts:
        return nightStarts[(surveyStartTime, nightNum)]

    # get the next setting after midnight on nightNum days after surveyStartTime
    startDate = datetime.fromtimestamp(surveyStartTime)
    startDate = startDate.replace(hour=12, minute=0, second=0)
    curDate = startDate + timedelta(nightNum)
    dublinJD = tel.next_setting(sun, start=curDate)

    # get the time as a unix timestamp
    unix = utils.mjd2unix(utils.djd2mjd(dublinJD))

    # cache the result
    nightStarts[(surveyStartTime, nightNum)] = unix
    return unix

def nightEnd(surveyStartTime, nightNum):
    """Gives the end time of the `nightNum`th night

    Parameters
    ----------
    surveyStartTime : float
        The time that the survey starts as a unix timestamp.
    nightNum : int
        The number (0 indexed) of nights after the start of the survey
        whose start time is desired. Note that the 0th night starts in the
        evening of whatever day `surveyStartTime` falls on.

    Returns
    -------
    The end time (as a unix timestamp) of the `nightNum`th night
    """

    # use this (and the corresponding line in nightStart) to return
    # constant-length nights (good for debugging)
    #return surveyStartTime + nightNum * 3600*24 + 12 * 3600

    # check if cached
    if (surveyStartTime, nightNum) in nightEnds:
        return nightEnds[(surveyStartTime, nightNum)]

    # get the next setting after midnight on nightNum days after surveyStartTime
    startDate = datetime.fromtimestamp(surveyStartTime)
    startDate = startDate.replace(hour=12, minute=0, second=0)
    curDate = startDate + timedelta(nightNum)
    dublinJD = tel.next_rising(sun, start=curDate)

    # get the time as a unix timestamp
    unix = utils.mjd2unix(utils.djd2mjd(dublinJD))

    # cache the result
    nightEnds[(surveyStartTime, nightNum)] = unix
    return unix

def nightNum(surveyStartTime, time):
    """Gives the integer night number corresponding to `time`

    Parameters
    ----------
    surveyStartTime : float
        The time that the survey starts as a unix timestamp.
    time : float
        The time whose night number we wish to calculate (as a unix timestamp)

    Returns
    -------
    The number of nights after `surveyStartTime` that `time` falls in

    Notes
    -----
    This function returns a night number even if `time` falls during the day.
    """
    return int((time - surveyStartTime) / 3600 / 24)

def raOfMeridian(time):
    """Gives the RA of the meridian at `time`

    Parameters
    ----------
    time : float
        The time of interest as a unix timestamp.

    Returns
    -------
    The RA (in radians) of the meridian at `time` at the Telescope site.
    """
    # get the RA of zenith
    ra, dec = altaz2radec(np.pi/2, 0., time)
    return ra

def radecOfSun(time):
    """Gives the RA and declination of the sun at `time`

    Parameters
    ----------
    time : float
        The time of interest as a unix timestamp.

    Returns
    -------
    A tuple (RA, dec) representing the location of the sun at `time`.
    """
    sun = ephem.Sun(utils.mjd2djd(utils.unix2mjd(time)))
    return (sun.ra, sun.dec)

def radecOfMoon(time):
    """Gives the RA and declination of the moon at `time`

    Parameters
    ----------
    time : float
        The time of interest as a unix timestamp.

    Returns
    -------
    A tuple (RA, dec) representing the location of the moon at `time`.
    """
    # ephem needs djds
    moon = ephem.Moon(utils.mjd2djd(utils.unix2mjd(time)))
    return (moon.ra, moon.dec)

def phaseOfMoon(time):
    """Gives the phase of the moon at `time`

    Parameters
    ----------
    time : float
        The time of interest

    Returns
    -------
    The phase of the moon (as a float between 0 and 1)
    """
    # ephem needs djds
    moon = ephem.Moon(utils.mjd2djd(utils.unix2mjd(time)))
    return moon.moon_phase

def getExpTime(ra, dec, otherstuff = None):
    """Gives the exposure time to use

    TODO I put this method here because I was originally assuming that if
    the exposure time were to vary, it would depend only on RA/dec and
    therefore would be the responsibility of this module to calculate.
    However there may be other reasons to vary exposure time, and if so,
    it may be better if this module had (a) function(s) to return more general
    information (like what region of the sky, or airmass, or seeing, etc)
    that could be used by the scheduler to decide on exposure time.

    Parameters
    ----------
    ra : float
        The RA of the observation
    dec : float
        The declination of the observation
    otherstuff : object
        This is a placeholder

    Returns
    -------
    The number of seconds that the shutter should stay open when pointing
    at the supplied RA and dec
    """
    return 30


def unix2lst(longitude, time):
    """Calculates the local sidereal time of the supplied time

    Parameters
    ----------
    longitude : float
        The longitude of the observatory in radians.
    time : float
        The time of interest as a unix timestamp.

    Returns
    -------
    The local sidereal time (LST) at the given `longitude` and `time`.
    """
    mjd = utils.unix2mjd(time)
    lst = palpy.gmst(mjd) + longitude
    lst %= 2*np.pi
    return lst

def _checkCoordInput(coord1, coord2):
    # helper for use in radec2altaz and altaz2radec
    # inputs (either ra/dec or alt/az)
    # must be ndarrays of equal size or both floats
    if isinstance(coord1, np.ndarray):
        isNumpy = True
        assert(isinstance(coord2, np.ndarray))
        assert(coord1.shape == coord2.shape)
    else:
        isNumpy = False
        assert(isinstance(coord1, numbers.Number))
        assert(isinstance(coord2, numbers.Number))
    return isNumpy


def radec2altaz(ra, dec, time):
    """Converts from RA/dec to altitude/azimuth

    Patameters
    ----------
    ra : float
        The RA of interest in radians.
    dec : float
        The declination of interest in radians.
    time : float
        The time of interest as a unix timestamp.

    Returns
    -------
    The altitude and azimuth corresponding to the provided `ra` and `dec`.
    The shapes of the returned values equal the shapes of the provided `ra`
    and `dec`.

    Notes
    -----
    As of this writing (June 18 2017), this method is a bottleneck in
    performance of the scheduler that I (Daniel Rothchild) am writing.
    Currently, the method takes about 16us to convert one ra/dec to alt/az.
    """

    # inputs must be ndarrays of equal size or both floats
    # this helper returns whether the inputs are numpy arrays or not
    isNumpy = _checkCoordInput(ra, dec)

    # code adapted from lsst-ts/ts_astrosky_model and lsst-ts/ts_dateloc
    lst = unix2lst(Telescope.longitude, time)
    ha = lst - ra
    if isNumpy:
        az, alt = palpy.de2hVector(ha.flatten().astype(float),
                                   dec.flatten().astype(float), Telescope.latitude)
        alt = alt.reshape(ra.shape)
        az = az.reshape(ra.shape)
    else:
        az, alt = palpy.de2h(float(ha), float(dec), Telescope.latitude)
    return alt, az

def altaz2radec(alt, az, time):
    """Converts from altitude/azimuth to RA/dec

    Patameters
    ----------
    alt : float
        The altitude of interest in radians.
    az : float
        The azimuth of interest in radians.
    time : float
        The time of interest as a unix timestamp.

    Returns
    -------
    The RA and dec corresponding to the provided `alt` and `az`.
    The shapes of the returned values equal the shapes of the provided `alt`
    and `az`.
    """

    # this helper returns whether the inputs are numpy arrays or not
    isNumpy = _checkCoordInput(alt, az)
    
    lst = unix2lst(Telescope.longitude, time)

    if isNumpy:
        ha, dec = palpy.dh2eVector(az.flatten().astype(float),
                                   alt.flatten().astype(float), Telescope.latitude)
        ha = ha.reshape(az.shape)
        dec = dec.reshape(az.shape)
    else:
        ha, dec = palpy.dh2e(float(az), float(alt), Telescope.latitude)

    ra = lst - ha
    return ra, dec
