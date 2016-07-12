# coding: utf-8

from __future__ import division

__author__ = "Danny Broberg, Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "June 23, 2016"

import unittest
import os

import numpy as np

from pymatgen.io.vasp.outputs import Locpot
from pycdt.corrections.utils import *

bl_path = os.path.join('..', '..', '..', 'test_files', 'bLOCPOT.gz')
dl_path = os.path.join('..', '..', '..', 'test_files', 'dLOCPOT.gz')

class StructureFunctionsTest(unittest.TestCase):
    def setUp(self):
        self.bl = Locpot.from_file(bl_path)
        self.dl = Locpot.from_file(dl_path)
        self.bs = self.bl.structure
        self.ds = self.dl.structure
        self.a = self.bs.lattice.matrix[0]
        self.b = self.bs.lattice.matrix[1]
        self.c = self.bs.lattice.matrix[2]

    def test_cleanlat(self):
        self.assertAlmostEqual(
            cleanlat(self.bs.lattice.matrix),
            [5.750183, 5.750183, 5.750183])

    def test_genrecip(self):
        brecip = map(np.array, [[-1.0926931033637688, 0.0, 0.0],
                                [0.0, -1.0926931033637688, 0.0],
                                [0.0, 0.0, -1.0926931033637688],
                                [0.0, 0.0, 1.0926931033637688],
                                [0.0, 1.0926931033637688, 0.0],
                                [1.0926931033637688, 0.0, 0.0]])
        recip =  genrecip(self.a, self.b, self.c, 1.3)
        for i in range(len(brecip)):
            self.assertTrue(np.isclose(recip[i], brecip[i]).all())

    def test_generate_reciprocal_vectors_squared(self):
        brecip = [1.1939782181387439 for i in range(6)]
        self.assertAlmostEqual(
            generate_reciprocal_vectors_squared(
                self.a, self.b, self.c, 1.3),
            brecip)

    def test_closestsites(self):
        pos = [0.0000, 2.8751, 2.8751]
        bsite, dsite = closestsites(self.bs, self.ds, pos)
        self.assertAlmostEqual(list(bsite[0].coords),
                               list(self.bs.sites[1].coords))
        self.assertAlmostEqual(list(dsite[0].coords),
                               list(self.ds.sites[0].coords))

    def test_find_defect_pos(self):
        """
        TODO: might be good to include more defect poscars to test all
        types of defect finding.
        """
        bdefect, ddefect = find_defect_pos(self.bs, self.ds)
        self.assertAlmostEqual(list(bdefect), [0, 0, 0])
        self.assertIsNone(ddefect)

if __name__ == '__main__':
    unittest.main()


