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

from pymatgen.core.structure import Structure
from pycdt.core.defectsmaker import *

class GetOptimizedScScaleTest(unittest.TestCase):
    def setUp(self):
        self.gaas_prim_struct = Structure.from_file('POSCAR_GaAs')

    def test_biggercell_wanted(self):
        pass

    def test_smallercell_wanted(self):
        pass


class DefectChargerSemiconductorTest(unittest.TestCase):
    def setUp(self):
        self.gaas_struct = Structure.from_file('POSCAR_GaAs')
        self.def_charger = DefectChargerSemiconductor(self.gaas_struct)

    def test_vacancy_charges(self):
        """
        Reference: PRB 71, 125207 (2005)
        """
        ga_vac_qs = self.def_charger.get_charges('vacancy', 'Ga')
        as_vac_qs = self.def_charger.get_charges('vacancy', 'As')
        self.assertIn(0, ga_vac_qs)
        self.assertIn(-3, ga_vac_qs)
        self.assertIn(1, as_vac_qs)
        self.assertIn(-3, as_vac_qs)

    def test_antisite_charges(self):
        ga_on_as_qs = self.def_charger.get_charges('antisite', 'Ga', 'As')
        as_on_ga_qs = self.def_charger.get_charges('antisite', 'Ga', 'As')
        pass

    def test_substitution_charges(self):
        s_impurity_qs = self.def_charger.get_charges('substitution', 'S')
        se_impurity_qs = self.def_charger.get_charges('substitution', 'Se')
        mg_impurity_qs = self.def_charger.get_charges('substitution', 'Mg')
        zn_impurity_qs = self.def_charger.get_charges('substitution', 'Zn')

    def test_interstitial_charges(self):
        """
        References:
        N interstitial: +1 to -3 [J. Phys.: Condens. Matter 20 (2008) 235231]
        Mn intersititial:
        Self interstital: arxiv.org/pdf/1101.1413.pdf
        """
        mn_qs = self.def_charger.get_charges('interstitial', 'Mn')
        n_qs = self.def_charger.get_charges('interstitial', 'N')
        self.assertIn(1, n_qs)
        self.assertIn(-3, n_qs)
        self_qs = self.def_charger.get_charges('interstitial')
        self.assertIn(-1, self_qs)
        self.assertIn(2, self_qs)
        pass


class DefectChargerInsulatorTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_vacancy_charges(self):
        pass

    def test_antisite_charges(self):
        pass

    def test_substitution_charges(self):
        pass

    def test_interstitial_charges(self):
        pass


class ChargedDefectsStructuresTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_vacancies(self):
        pass

    def test_antisites(self):
        pass

    def test_substitutions(self):
        pass

    def test_interstitials(self):
        pass


    def test_anion_on_cation_substitution(self):
        """
        Check the charge states for cases like As substitution on Ga site
        """
        pass


if __name__ == '__main__':
    unittest.main()
