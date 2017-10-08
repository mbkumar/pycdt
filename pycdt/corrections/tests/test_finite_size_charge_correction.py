# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "Jan 26, 2016"

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp.outputs import Locpot
from pymatgen.entries.computed_entries import ComputedStructureEntry

from pycdt.corrections.finite_size_charge_correction import *
from pycdt.core.defects_analyzer import ComputedDefect

#Paths to locpots we are testing on
bl_path = os.path.join('..', '..', '..', 'test_files', 'bLOCPOT.gz')
dl_path = os.path.join('..', '..', '..', 'test_files', 'dLOCPOT.gz')

class FiniteSizeChargeCorrectionTest(PymatgenTest):
    """
    Test functions for getting freysoldt and kumagai as well as ChargeCorrection class itself
    if one breaks the other will likely break as well...
    """
    def setUp(self):
        self.bl=Locpot.from_file(bl_path)
        self.dl=Locpot.from_file(dl_path)
        self.bs=self.bl.structure
        self.ds=self.dl.structure
        self.bulk_entry = ComputedStructureEntry(
                self.bs, 100, data={'locpot_path': bl_path, 'encut': 520})
        self.defect_entry = ComputedStructureEntry(
                self.ds, 100, data={'locpot_path': dl_path, 'encut': 520,
                                    'charge': -3})
        self.computed_defect = ComputedDefect(
                self.defect_entry, self.bs.sites[0], charge=-3,
                name='vac_1_Ga')
        self.kbi = KumagaiBulkInit(self.bs, self.bl.dim, 15,
                                   optgamma=3.49423226983)

    def test_get_correction_freysoldt(self):
        freyout = get_correction_freysoldt(self.computed_defect,
                                           self.bulk_entry, 15)
        self.assertAlmostEqual(freyout[0], 3.99126)
        self.assertIsInstance(freyout[1], Locpot)

    def test_get_correction_kumagai(self):
        kumagaiout = get_correction_kumagai(self.computed_defect, './',
                                            self.kbi, bulk_locpot=self.bl)
        self.assertAlmostEqual(kumagaiout, 4.24073)

    def test_chargecorrectionclass(self):
        # Could make this a seperate unit test class,
        # but kept it here to speed things up...
        cc = ChargeCorrection(15, bl_path, dl_path, -3, pure_locpot=self.bl,
                              defect_locpot=self.dl, optgamma=3.49423226983,
                              KumagaiBulk=self.kbi)
        self.assertIsInstance(cc, ChargeCorrection)
        freyout = cc.freysoldt()
        self.assertAlmostEqual(freyout, 3.99126)
        #Kumagai IS BROKEN
        # kumagaiout = cc.kumagai()
        # self.assertAlmostEqual(kumagaiout, 4.24073)


import unittest
if __name__ == '__main__':
    unittest.main()
