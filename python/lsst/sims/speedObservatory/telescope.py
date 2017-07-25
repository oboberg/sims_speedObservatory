from __future__ import division
import numpy as np

__all__ = ["Telescope"]

class Telescope:
    """This class provides a kinematic model of the telescope

    The functionality provided is currently (June 18 2017) limited to
    calculating slew times. Putting the code into a class allows a
    simulator to consider many `Telescope`s at once, each with
    different parameters.

    Attributes
    ----------
    latitude : float
        The latitude of the telescope in radians. This attribute should not
        be changed on Telescope instances, since lsst.sims.speedObservatory.sky
        assumes the default value.
    longitude : float
        The longitude of the telescope in radians. This attribute should not
        be changed on Telescope instances, since lsst.sims.speedObservatory.sky
    filters : list of str
        The filters available on this telescope. This attribute should not
        be modified on instances.
    filterId : dict (keys are str, values are int)
        An integer ID assigned to each filter. This attribute should not
        be modified on instances.
    fovWidth : float
        The width of the field of view in radians.
    domSlitDiam : float
        The width of the dome slit in radians. TODO this assumes a circular
        dome slit, but I think it may actually be rectangular.
    raftWidth : float
        The width of a single raft in radians.
    minRotation : float
        The minimum angle of the rotator in radians.
    maxRotation : float
        The maximum angle of the rotator in radians.
    minAlt : float
        The minimum altitude that the telescope can point at in radians.
    maxAlt : float
        The maximum altitude that the telescope can point at in radians.
    domAltMaxSpeed : float
        The maximum rate of change in dome altitude in radians/s.
    domAltAccel : float
        The maximum acceleration in dome altitude in radians/s/s.
    domAltDecel : float
        The maximum deceleration in dome altitude in radians/s/s.
    domAzMaxSpeed : float
        The maximum rate of change in dome azimuth in radians/s.
    domAzAccel : float
        The maximum acceleration in dome azimuth in radians/s/s.
    domAzDecel : float
        The maximum deceleration in dome azimuth in radians/s/s.
    telAltMaxSpeed : float
        The maximum rate of change in telescope altitude in radians/s.
    telAltAccel : float
        The maximum acceleration in telescope altitude in radians/s/s.
    telAltDecel : float
        The maximum deceleration in telescope altitude in radians/s/s.
    telAzMaxSpeed : float
        The maximum rate of change in telescope azimuth in radians/s.
    telAzAccel : float
        The maximum acceleration in telescope azimuth in radians/s/s.
    telAzDecel : float
        The maximum deceleration in telescope azimuth in radians/s/s.
    rotMaxSpeed : float
        The maximum rate of change in rotator in radians/s.
    rotAccel : float
        The maximum acceleration of the rotator in radians/s/s.
    rotDecel : float
        The maximum deceleration of the rotator in radians/s/s.
    settleTime : float
        The number of seconds it takes for the instrument to settle.
    readoutTime : float
        The number of seconds it takes to read out an exposure from the camera.
    filterChangeTime : float
        The number of seconds it takes to change filters.
    """
    # non-configurable (AstronomicalSky uses these and I don't
    # want to slow it down by requiring it create a Telescope object)
    # (don't modify these on an instance)
    latitude = np.radians(-(30 + 14 / 60 + 40.7 / 3600))
    longitude = np.radians(-(70 + 44 / 60 + 57.9 / 3600))
    filters = ["u", "g", "r", "i", "z", "y"]
    filterId = {}
    for i, filter in enumerate(filters):
        filterId[filter] = i

    def __init__(self):
        self.fovWidth = np.radians(3.5)
        self.domSlitDiam = 2 * self.fovWidth
        self.raftWidth = self.fovWidth / 5
        self.minRotation = -np.pi/2
        self.maxRotation = np.pi/2

        self.minAlt = np.radians(30)
        # values from http://ops2.lsst.org/docs/current/system.html
        self.maxAlt = np.radians(86.5)

        # Kinematic and delay parameters for slew time computation

        # speed in rads/sec
        # acceleration in rads/second**2
        self.domAltMaxSpeed = np.radians(1.75)
        self.domAltAccel = np.radians(0.875)
        self.domAltDecel = np.radians(0.875)

        self.domAzMaxSpeed = np.radians(1.5)
        self.domAzAccel = np.radians(0.75)
        self.domAzDecel = np.radians(0.75)

        self.telAltMaxSpeed = np.radians(3.5)
        self.telAltAccel = np.radians(3.5)
        self.telAltDecel = np.radians(3.5)
        # assume accel == decel for calculations below
        # (easy to change but they are the same anyway and I'm  lazy)
        assert(self.telAltAccel == self.telAltDecel)

        self.telAzMaxSpeed = np.radians(7.0)
        self.telAzAccel = np.radians(7.0)
        self.telAzDecel = np.radians(7.0)
        assert(self.telAzAccel == self.telAzDecel)

        # not used in slew calculation
        self.rotMaxSpeed = np.radians(3.5)
        self.rotAccel = np.radians(1.0)
        self.rotDecel = np.radians(1.0)

        self.settleTime = 3
        self.readoutTime = 2
        self.filterChangeTime = 120

    def calcSlewTime(self, alt1, az1, filter1, alt2, az2, filter2, laxDome=False):
        """Calculates ``slew'' time

        Calculates the ``slew'' time necessary to get from alt1/az1/filter1
        to alt2/az2/filter2. The time returned is actually the time between
        the end of an exposure at alt1/az1 and the beginning of an exposure
        at alt2/az2, since it includes readout time in the ``slew'' time.

        Parameters
        ----------
        alt1 : float
            The altitude of the initial pointing.
        az1 : float
            The azimuth of the initial pointing.
        filter1 : str
            The filter used in the initial observation.
        alt2 : float
            The altitude of the destination pointing.
        az2 : float
            The azimuth of the destination pointing.
        filter2 : str
            The filter to be used in the destination observation.
        laxDome : boolean
            If True, allow the dome to creep, model a dome slit, and don't
            require the dome to settle in azimuth. If False, adhere to the way
            SOCS calculates slew times (as of June 21 2017).

        Returns
        -------
        The number of seconds between the two specified exposures.

        Notes
        -----
        This method should really be called `calcInterExposureTime`, but to be
        consistent with other code/documentation, I've called it `calcSlewTime`.
        """

        # FYI this takes on the order of 10us for 1 slew calculation

        # TODO also assumes we never max out the cable wrap-around constraint
        deltaAlt = np.abs(alt2 - alt1)
        deltaAz  = np.abs(az2 - az1)

        deltaAz = min(deltaAz, np.abs(deltaAz - 2*np.pi))

        def uamSlewTime(d, vmax, a):
            # if you accelerate uniformely to telAltMaxSpeed
            # and then slow down uniformely to zero, you'll travel
            # a distance v_max^2 / a
            if d < vmax**2 / a:
                # to travel a distance d/2 while accelerating at a rate a,
                # it takes time sqrt(2(d/2)/a)
                slewTime = 2 * np.sqrt(d / a)
            else:
                # the time to accelerate/decelerate to/from v_max is 2v_max / a
                # and the distance covered in those two steps is v_max^2 / a
                # so the total time is the accel/decel time plus the remaining
                # distance over v_max
                slewTime = 2 * vmax / a + (d - vmax**2 / a) / vmax
            return slewTime

        telAltSlewTime = uamSlewTime(deltaAlt, self.telAltMaxSpeed, self.telAltAccel)
        telAzSlewTime  = uamSlewTime(deltaAz,  self.telAzMaxSpeed,  self.telAzAccel)
        totTelTime = max(telAltSlewTime, telAzSlewTime)

        # open loop optics correction
        olTime = deltaAlt / np.radians(3.5)
        totTelTime += olTime

        if totTelTime > 0:
            # OL and settle can happen at the same time
            totTelTime += max(0, self.settleTime - olTime)
        # readout puts a floor on tel time
        totTelTime = max(self.readoutTime, totTelTime)

        # now compute dome slew time
        if laxDome:
            # model dome creep, dome slit, and no azimuth settle

            # if we can fit both exposures in the dome slit, do so
            if deltaAlt**2 + deltaAz**2 < self.fovWidth**2:
                totDomTime = 0
            else:
                # else, we take the minimum time from two options:
                # 1. assume we line up alt in the center of the dome slit so we
                #    minimize distance we have to travel in azimuth.
                # 2. line up az in the center of the slit
                # also assume:
                # * that we start out going maxspeed for both alt and az
                # * that we only just barely have to get the new field in the
                #   dome slit in one direction, but that we have to center the
                #   field in the other (which depends which of the two options used)
                # * that we don't have to slow down until after the shutter
                #   starts opening
                domDeltaAlt = deltaAlt

                # on each side, we can start out with the dome shifted away from
                # the center of the field by an amount domSlitRadius - fovRadius
                domDeltaAz = deltaAz - 2 * (self.domSlitDiam/2 - self.fovWidth/2)

                domAltSlewTime = domDeltaAlt / self.domAltMaxSpeed
                domAzSlewTime  = domDeltaAz  / self.domAzMaxSpeed

                totDomTime1 = max(domAltSlewTime, domAzSlewTime)

                domDeltaAlt = deltaAlt - 2 * (self.domSlitDiam/2 - self.fovWidth/2)
                domDeltaAz  = deltaAz
                domAltSlewTime = domDeltaAlt / self.domAltMaxSpeed
                domAzSlewTime  = domDeltaAz  / self.domAzMaxSpeed
                totDomTime2 = max(domAltSlewTime, domAzSlewTime)

                totDomTime = min(totDomTime1, totDomTime2)

        else:
            # the above models a dome slit and dome creep. However, it appears that
            # SOCS requires the dome to slew exactly to each field and settle in az
            domAltSlewTime = uamSlewTime(deltaAlt, self.domAltMaxSpeed, self.domAltAccel)
            domAzSlewTime  = uamSlewTime(deltaAz,  self.domAzMaxSpeed,  self.domAzAccel)
            if domAzSlewTime > 0:
                domAzSlewTime += 1 # dome takes 1 second to settle in az
            totDomTime = max(domAltSlewTime, domAzSlewTime)


        slewTime = max(totTelTime, totDomTime)

        # include filter change time if necessary
        if filter1 != filter2:
            slewTime = max(slewTime, self.filterChangeTime)

        # closed loop optics correction
        if deltaAlt >= np.radians(9):
            if slewTime == domAzSlewTime:
                # XXX more SOCS-similitude code
                # I thought we can do closed loop w/o waiting for dome az settle
                # so this would be += 19, but apparently not?
                slewTime += 20
            else:
                slewTime += 20

        return slewTime
