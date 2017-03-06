# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "Jan 14, 2017"

import unittest
import os
import glob

from monty.json import MontyDecoder
from monty.tempfile import ScratchDir
from monty.serialization import loadfn
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pycdt.core.chemical_potentials import ChemPotAnalyzer

file_loc = os.path.join('..', '..', '..', 'test_files')

class ChemPotAnalyzerTest(unittest.TestCase):
    def setUp(self):
        self.bulk_comp = Composition("Cr2O3")

    def test_mp_entries_stable_mpid(self):
        """
        Test the case where all the entries are from MP and the 
        mpid corresponds to entry on convex hull
        """
        mpid = "mp-19399"
        cpa = ChemPotAnalyzer(self.bulk_comp)
        chem_pot = cpa.analyze_GGA_chempots(mpid=mpid)
        print (chem_pot)

    def test_mp_entries_unstable_mpid(self):
        """
        Test the case where all the entries are from MP and the 
        mpid corresponds to entry not on convex hull
        """
        mpid = "mp-776999"
        cpa = ChemPotAnalyzer(self.bulk_comp)
        chem_pot = cpa.analyze_GGA_chempots(mpid=mpid)
        print (chem_pot)

    def test_mp_entries_structure(self):
        """
        Test the case where all the entries are from MP and the 
        entry corresponds to local calculation (with no mpid)
        """
        pass

    def test_local(self):
        """
        Test the case where all the entries are from local 
        calculation.
        """
        pass


if __name__ == '__main__':
    unittest.main()
