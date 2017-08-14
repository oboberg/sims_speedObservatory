import numpy as np
import unittest
import lsst.sims.speedObservatory as speedo
import lsst.utils.tests
from lsst.sims.featureScheduler import empty_observation


class TestSpeedObs(unittest.TestCase):

    def teststubb(self):
        so = speedo.Speed_observatory()
        # Check that we can get a status
        status = so.return_status()

        # Check that we can get an observation
        obs = empty_observation()
        obs['dec'] = np.radians(-30.)
        obs['filter'] = 'r'
        obs['exptime'] = 30.
        obs['nexp'] = 2.
        result = so.observe(obs)

        status2 = so.return_status()

        assert(status2['mjd'] > status['mjd'])

        obs['dec'] = np.radians(-35.)
        result = so.observe(obs)

        assert(result['airmass'] >= 1.)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
