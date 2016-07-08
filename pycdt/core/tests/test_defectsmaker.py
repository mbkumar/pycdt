# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "June 6, 2016"

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
        as_on_ga_qs = self.def_charger.get_charges('antisite', 'As', 'Ga')
        self.assertIn(0, ga_on_as_qs)
        self.assertIn(2, ga_on_as_qs)
        self.assertIn(0, as_on_ga_qs)
        self.assertIn(-3, as_on_ga_qs)


    def test_substitution_charges(self):
        """
        TODO: This needs test
        """
        s_impurity_qs = self.def_charger.get_charges('substitution', 'As', 'S')
        se_impurity_qs = self.def_charger.get_charges('substitution', 'As', 'Se')
        mg_impurity_qs = self.def_charger.get_charges('substitution', 'Ga', 'Mg')
        zn_impurity_qs = self.def_charger.get_charges('substitution', 'Ga', 'Zn')

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
        self_qs = self.def_charger.get_charges('interstitial', 'Ga')
        self.assertIn(-1, self_qs)
        self.assertIn(2, self_qs)


class DefectChargerInsulatorTest(unittest.TestCase):
    def setUp(self):
        cr2o3_struct = Structure.from_file('POSCAR_Cr2O3')
        self.def_charger = DefectChargerInsulator(cr2o3_struct)

    def test_vacancy_charges(self):
        """
        For insulators, the range of defect charges expected is
        -A to 0 or 0 to B, where A is cation charge in its native
        oxidation state in the compound, and B is the corresponding
        anion charge.
        """
        cr_vac_qs = self.def_charger.get_charges('vacancy', 'Cr')
        o_vac_qs = self.def_charger.get_charges('vacancy', 'O')
        self.assertIn(0, cr_vac_qs)
        self.assertIn(-3, cr_vac_qs)
        self.assertNotIn(-4, cr_vac_qs)
        self.assertNotIn(1, cr_vac_qs)
        self.assertIn(0, o_vac_qs)
        self.assertIn(2, o_vac_qs)
        self.assertNotIn(3, o_vac_qs)
        self.assertNotIn(-1, o_vac_qs)

    def test_antisite_charges(self):
        """
        Anitisites are not expected for insulators.
        Skipping this.
        """
        pass

    def test_substitution_charges(self):
        ti_on_cr_qs = self.def_charger.get_charges('substitution', 'Cr', 'Ti')
        self.assertIn(0, ti_on_cr_qs)
        self.assertIn(1, ti_on_cr_qs)
        self.assertNotIn(-1, ti_on_cr_qs)
        self.assertNotIn(2, ti_on_cr_qs)
        mg_on_cr_qs = self.def_charger.get_charges('substitution', 'Cr', 'Mg')
        self.assertIn(-1, mg_on_cr_qs)
        self.assertNotIn(0, mg_on_cr_qs)
        self.assertNotIn(-2, mg_on_cr_qs)
        self.assertNotIn(1, mg_on_cr_qs)

    def test_interstitial_charges(self):
        ti_inter_qs = self.def_charger.get_charges('interstitial', 'Ti')
        self.assertIn(0, ti_inter_qs)
        self.assertIn(4, ti_inter_qs)
        self.assertNotIn(-1, ti_inter_qs)
        self.assertNotIn(5, ti_inter_qs)



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
