import numpy as np
import healpy as hp
import copy
from lsst.sims.ocs.observatory import MainObservatory
from lsst.sims.ocs.configuration import ObservingSite, Observatory
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
    slew_time = (final_slew_state.time - initial_slew_state.time, "seconds")

    return slew_time

# initialize opsim observatory
# Easy to just make an observatory and pull out the configured model
mo = MainObservatory(ObservingSite())
mo.configure(Observatory())
model = mo.model

# initialize alt-sched observatory
speedObs = SpeedObservatory()


# calculate slew times using both methods over alt and az grids
delta = 5

azimuths = np.arange(0, 360+delta, delta)
altitudes = np.arange(20., 90.+delta, delta)

opsim_alt_array = np.zeros((altitudes.size, altitudes.size), dtype=float)
opsim_az_array  = np.zeros((azimuths.size,  azimuths.size), dtype=float)

fast_alt_array = np.zeros((altitudes.size, altitudes.size), dtype=float)
fast_az_array  = np.zeros((azimuths.size,  azimuths.size), dtype=float)

for i in range(azimuths.size):
    for j in range(azimuths.size):
        opsim_slew_time = opsimSlewTime(model, az1=azimuths[i], az2=azimuths[j])
        opsim_az_array[i, j] += opsim_slew_time[0]

        fast_slew_time = speedObs.calcSlewTime(np.array([0,az1]), 'u', 
                                               np.array([0,az2]), 'u')
        fast_az_array[i, j] += fast_slew_time

plt.imshow(opsim_az_array, interpolation='none', extent=[0,360,0,360])
plt.imshow(alt_az_array,   interpolation='none', extent=[0,360,0,360])
plt.show()

for i in range(altitudes.size):
    for j in range(altitudes.size):
        st = slewtime(model, alt1=altitudes[i], alt2=altitudes[j])
        alt_array[i, j] += st[0]
