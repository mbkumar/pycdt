# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "May 6, 2015"

import os
import unittest
import tarfile
from shutil import copyfile

import numpy as np
from monty.serialization import loadfn, dumpfn
from monty.json import MontyDecoder, MontyEncoder
from monty.tempfile import ScratchDir

from pymatgen import __file__ as initfilep
from pymatgen.core import Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.lattice import Lattice
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.vasp import Locpot
from pymatgen.util.testing import PymatgenTest

from pycdt.core.defects_analyzer import freysoldt_correction_from_paths, \
        kumagai_correction_from_paths

PMG_TEST_DIR = os.path.join(
        os.path.split(os.path.split(initfilep)[0])[0], 'test_files')
TEST_DIR = os.path.abspath(os.path.join(
        __file__, '..', '..', '..', '..', 'test_files'))

class FilePathCorrectionsTest(PymatgenTest):
    def test_freysoldt_and_kumagai(self):
        with ScratchDir("."):
            # setup with fake Locpot object copied over
            copyfile(os.path.join(TEST_DIR, "test_path_files.tar.gz"),
                     "test_path_files.tar.gz")
            tar = tarfile.open("test_path_files.tar.gz")
            tar.extractall()
            tar.close()
            blocpot = Locpot.from_file(os.path.join(TEST_DIR, "bLOCPOT.gz"))
            blocpot.write_file(os.path.join("test_path_files", "bulk", "LOCPOT"))
            dlocpot = Locpot.from_file(os.path.join(TEST_DIR, "dLOCPOT.gz"))
            dlocpot.write_file(os.path.join("test_path_files", "sub_1_Sb_on_Ga",
                                            "charge_2", "LOCPOT"))

            fcc = freysoldt_correction_from_paths(
                    os.path.join("test_path_files", "sub_1_Sb_on_Ga", "charge_2"),
                    os.path.join("test_path_files", "bulk"), 18.12, 2, plot=True)
            self.assertAlmostEqual(fcc, -1.2435280589593547 )
            self.assertTrue(os.path.exists(
                os.path.join("test_path_files", "sub_1_Sb_on_Ga", "charge_2",
                             "Sub_Sb_on_Ga_mult32_chg_2_axis1_freysoldtplot.pdf")))

            kcc = kumagai_correction_from_paths(
                    os.path.join("test_path_files", "sub_1_Sb_on_Ga", "charge_2"),
                    os.path.join("test_path_files", "bulk"),
                    18.12, 2, plot=True)
            self.assertAlmostEqual(kcc, 0.6387768530616106)
            self.assertTrue(os.path.exists(
                os.path.join("test_path_files", "sub_1_Sb_on_Ga", "charge_2",
                             "Sub_Sb_on_Ga_mult32_chg_2_kumagaiplot.pdf")))


if __name__ == '__main__':
    unittest.main()
