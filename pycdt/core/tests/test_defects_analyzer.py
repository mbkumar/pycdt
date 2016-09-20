# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "May 6, 2015"

import unittest
import os

from monty.serialization import loadfn, dumpfn
from monty.json import MontyDecoder, MontyEncoder
from monty.tempfile import ScratchDir
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.lattice import Lattice
from pycdt.core.defects_analyzer import ComputedDefect, DefectsAnalyzer

file_loc = os.path.join('..', '..', '..', 'test_files')

class ComputedDefectTest(unittest.TestCase):
    def setUp(self):
        entry_file = os.path.join(file_loc, 'vac_Cr2o3_struct_entry.json')
        entry = loadfn(entry_file, cls=MontyDecoder)
        lattice = Lattice([[9.995004137201189, -2.1469568e-08, 0.0],
                           [-4.997501105922451, 8.655927903729987, 0.0],
                           [0.0, 0.0, 13.67956098598296]])
        coords = [0.1666665000000016, 0.3333334999999984, 0.014505185117094302]
        site_in_bulk = PeriodicSite('Cr', coords, lattice)
        multiplicity = 12
        supercell_size = [2, 2, 1]
        q = -2
        q_corr =  0.98716
        o_corr = 0.39139821874799996
        name = "vac_1_Cr"
        self.com_def = ComputedDefect(
            entry, site_in_bulk, multiplicity, supercell_size, q, q_corr,
            o_corr, name)

    def test_as_from_dict(self):
        d = self.com_def.as_dict()
        comp_def = ComputedDefect.from_dict(d)
        self.assertIsInstance(comp_def, ComputedDefect)
        with ScratchDir('.'):
            dumpfn(self.com_def, 'tmp.json', cls=MontyEncoder)
            comp_def = loadfn('tmp.json', cls=MontyDecoder)
            self.assertIsInstance(comp_def, ComputedDefect)


class DefectsAnalyzerTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_from_dict(self):
        pass

    def test_as_dict(self):
        pass

    def test_add_parsed_defect(self):
        pass

    def test_change_charge_correction(self):
        pass

    def test_change_other_correction(self):
        pass

    def test_correct_bg_simple(self):
        pass

    def test_get_transition_levels(self):
        pass

    def test_correct_bg(self):
        pass

    def test_get_defect_occupancies(self):
        pass

    def test_get_formation_energies(self):
        pass

    def test_get_defects_concentration(self):
        pass


if __name__ == '__main__':
    unittest.main()


