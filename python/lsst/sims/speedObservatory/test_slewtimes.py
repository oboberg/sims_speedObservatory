import numpy as np
import healpy as hp
import copy
from lsst.sims.ocs.observatory import MainObservatory
from lsst.sims.ocs.configuration import ObservingSite, Observatory
from telescope import Telescope
# Generate a few matrices to make it easy to generate slewtime maps

from matplotlib import pyplot as plt
import time

def opsimSlewTime(model, alt1=None, alt2=None, az1=None, az2=None, rot1=None, rot2=None):
    # park the observatory, then go from alt1,az1 to alt2,az2
    if alt1 is None:
        alt1 = model.current_state.alt
    if alt2 is None:
        alt2 = model.current_state.alt
    if az1 is None:
        az1 = model.current_state.az
    if az2 is None:
        az2 = model.current_state.az
    if rot1 is None:
        rot1 = model.current_state.rot
    if rot2 is None:
        rot2 = model.current_state.rot

    model.park()
    model.slew_altaz(0., np.radians(alt1), np.radians(az1), np.radians(rot1), model.current_state.filter)
    initial_slew_state = copy.deepcopy(model.current_state)

    model.slew_altaz(model.current_state.time, np.radians(alt2), np.radians(az2), np.radians(rot2),
                     model.current_state.filter)

    final_slew_state = copy.deepcopy(model.current_state)
    slew_time = final_slew_state.time - initial_slew_state.time

    return slew_time

# initialize opsim observatory
# Easy to just make an observatory and pull out the configured model
mo = MainObservatory(ObservingSite())
mo.configure(Observatory())
model = mo.model

# initialize alt-sched observatory
tel = Telescope()

numRand = 10000
randMax = np.radians(30)
deltaAlts = np.random.rand(numRand) * np.sqrt(randMax)
deltaAzes = np.random.rand(numRand) * np.sqrt(randMax)
opsim_slew_times = []
fast_slew_times = []
for deltaAlt, deltaAz in zip(deltaAlts, deltaAzes):
    # make altitude start at 45 degs so opsim doesn't complain
    # if we ask it to point at the ground
    opsim_slew_time = opsimSlewTime(model, alt1=45,
                                    alt2=np.degrees(deltaAlt)+45,
                                    az1=0,
                                    az2=np.degrees(deltaAz))
    fast_slew_time = tel.calcSlewTime(0, 0, 'u',
                                      deltaAlt, deltaAz, 'u')
    opsim_slew_times.append(opsim_slew_time)
    fast_slew_times.append(fast_slew_time)
opsim_slew_times = np.array(opsim_slew_times)
fast_slew_times = np.array(fast_slew_times)
plt.scatter(fast_slew_times, fast_slew_times - opsim_slew_times)
plt.title("Error for random slews")
plt.xlabel("Slew time (seconds)")
plt.ylabel("Error (fast - opsim) (seconds)")

# calculate slew times using both methods over alt and az grids

# we'll consider slews between az=azMin...azMax; and alt=altMin..altMax
# at a resolution of delta
azMin  = np.radians(0)
azMax  = np.radians(90) # set to 90 instead of 360 to save time
altMin = np.radians(20)
altMax = np.radians(80)
delta  = np.radians(0.5)

azimuths = np.arange(azMin, azMax+delta, delta)
altitudes = np.arange(altMin, altMax+delta, delta)

# opsim_alt_array, for example, will hold the time, according to opsimSlewTime,
# that it takes to slew between every altitudes.size^2 combindation of altitudes
# while keeping the azimuth constant. The other 3 arrays are similar
opsim_alt_array = np.zeros((altitudes.size, altitudes.size), dtype=float)
opsim_az_array  = np.zeros((azimuths.size,  azimuths.size),  dtype=float)

fast_alt_array = np.zeros((altitudes.size, altitudes.size), dtype=float)
fast_az_array  = np.zeros((azimuths.size,  azimuths.size),  dtype=float)

# populate the opsim and fast azimuth arrays
for i in range(azimuths.size):
    for j in range(azimuths.size):
        opsim_slew_time = opsimSlewTime(model, az1=np.degrees(azimuths[i]),
                                               az2=np.degrees(azimuths[j]))
        opsim_az_array[i, j] = opsim_slew_time

        fast_slew_time = tel.calcSlewTime(0, azimuths[i], 'u',
                                          0, azimuths[j], 'u')
        fast_az_array[i, j] = fast_slew_time

# and populate the opsim and fast altitude arrays
for i in range(altitudes.size):
    for j in range(altitudes.size):
        opsim_slew_time = opsimSlewTime(model, alt1=np.degrees(altitudes[i]),
                                               alt2=np.degrees(altitudes[j]))
        opsim_alt_array[i, j] = opsim_slew_time

        fast_slew_time = tel.calcSlewTime(altitudes[i], 0, 'u',
                                          altitudes[j], 0, 'u')
        fast_alt_array[i, j] = fast_slew_time

# convenience function for plotting the arrays
def plot(arr, minRad, maxRad, title="", xlabel="", ylabel="", clim=None):
    plt.figure(title)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.imshow(arr[::-1,:], interpolation='none',
               extent=np.degrees([minRad,maxRad,minRad,maxRad]),
               clim=clim)
    plt.colorbar()

# plot opsim, fast, and the difference for alt and az arrays
plot(opsim_az_array, azMin, azMax, "SOCS Azimuth Slew Times",
     "Final azimuth (deg)", "Initial azimuth (deg)")
plot(fast_az_array, azMin, azMax, "speedObservatory Azimuth Slew Times",
     "Final azimuth (deg)", "Initial azimuth (deg)")
plot(fast_az_array - opsim_az_array, azMin, azMax,
     "Azimuth Slew Time Error (speedObservatory - SOCS)",
     "Final azimuth (deg)", "Initial azimuth (deg)", clim=(-0.01, 0.01))

plot(opsim_alt_array, altMin, altMax, "SOCS Altitude Slew Times",
     "Final altitude (deg)", "Initial altitude (deg)")
plot(fast_alt_array, altMin, altMax, "speedObservatory Altitude Slew Times",
     "Final altitude (deg)", "Initial altitude (deg)")
plot(fast_alt_array - opsim_alt_array, altMin, altMax,
     "Altitude Slew Time Error (speedObservatory - SOCS)",
     "Final altitude (deg)", "Initial altitude (deg)", clim=(-0.01, 0.01))

plt.show()
