import numpy as np
import unittest
import lsst.sims.speedObservatory as speedo
import lsst.utils.tests


class TestSpeedObs(unittest.TestCase):

    def teststubb(self):
        so = speedo.Speed_observatory()
        pass


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
