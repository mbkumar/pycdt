import unittest

from pymatgen.matproj.rest import MPRester
from pycdt.utils.get_vbm import get_vbm

class GetVBMTest(unittest.TestCase):
    def setUp(self):
        mpid = 'mp-23070'   
        mpr = MPRester()
        bs = mpr.get_bandstructure_by_material_id(mpid)
        self.vbm = bs.get_vbm()['energy']

    def test_vbm(self):
        self.assertAlmostEqual(self.vbm, 0.0, places=1)


if __name__ == "__main__":
    unittest.main()
