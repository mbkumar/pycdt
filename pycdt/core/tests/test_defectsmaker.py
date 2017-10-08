# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "June 6, 2016"

import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.core import PeriodicSite
from pycdt.core.defectsmaker import *

file_loc = os.path.join('..', '..', '..', 'test_files')

class GetOptimizedScScaleTest(PymatgenTest):
    def setUp(self):
        self.gaas_prim_struct = Structure.from_file(
                os.path.join(file_loc, 'POSCAR_GaAs'))

    def test_biggercell_wanted(self):
        lattchange = get_optimized_sc_scale(self.gaas_prim_struct, 300)
        self.assertEqual([5, 5, 5], lattchange)
        lattchange = get_optimized_sc_scale(self.gaas_prim_struct, 100)
        self.assertEqual([3, 3, 3], lattchange)


class DefectChargerSemiconductorTest(PymatgenTest):
    def setUp(self):
        self.gaas_struct = Structure.from_file(
            os.path.join(file_loc, 'POSCAR_GaAs'))
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
        As Self interstital: arxiv.org/pdf/1101.1413.pdf
        """
        n_qs = self.def_charger.get_charges('interstitial', 'N')
        self.assertIn(1, n_qs)
        self.assertIn(-3, n_qs)
        self_qs = self.def_charger.get_charges('interstitial', 'As')
        self.assertIn(-1, self_qs)
        self.assertIn(1, self_qs)


class DefectChargerInsulatorTest(PymatgenTest):
    def setUp(self):
        cr2o3_struct = Structure.from_file(
                os.path.join(file_loc, 'POSCAR_Cr2O3'))
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


class ChargedDefectsStructuresTest(PymatgenTest):
    def setUp(self):
        self.gaas_struct = Structure.from_file(
            os.path.join(file_loc, 'POSCAR_GaAs'))
        self.ga_site = self.gaas_struct.sites[0]
        self.as_site = self.gaas_struct.sites[1]

    def test_simple_initialization(self):
        #test simple attributes
        CDS = ChargedDefectsStructures(self.gaas_struct)
        self.assertIsInstance(CDS.struct, Structure)
        self.assertIsInstance(CDS.defect_charger, DefectChargerSemiconductor)
        self.assertFalse(CDS.substitutions)
        superstruct = CDS.defects['bulk']['supercell']['structure']
        self.assertIsInstance(superstruct, Structure)
        cellsize = CDS.defects['bulk']['supercell']['size']
        self.assertEqual([4, 4, 4], cellsize)

        #test native (intrinsic) defect generation
        self.assertEqual('vac_1_Ga', CDS.defects['vacancies'][0]['name'])
        self.assertEqual(self.ga_site, CDS.defects['vacancies'][0]['unique_site'])
        self.assertEqual('vac_2_As', CDS.defects['vacancies'][1]['name'])
        self.assertEqual(self.as_site, CDS.defects['vacancies'][1]['unique_site'])

        self.assertEqual('as_1_As_on_Ga', CDS.defects['substitutions'][0]['name'])
        self.assertEqual(self.ga_site, CDS.defects['substitutions'][0]['unique_site'])
        self.assertEqual('as_2_Ga_on_As', CDS.defects['substitutions'][1]['name'])
        self.assertEqual(self.as_site, CDS.defects['substitutions'][1]['unique_site'])


    def test_extra_initialization(self):
        CDS = ChargedDefectsStructures(self.gaas_struct, cellmax = 513,
                                            struct_type='insulator',
                                            antisites_flag=False)
        self.assertIsInstance(CDS.defect_charger, DefectChargerInsulator)
        cellsize = CDS.defects['bulk']['supercell']['size']
        self.assertEqual([5, 5, 5], cellsize)
        self.assertFalse(len(CDS.defects['substitutions'])) #testing antisite flag


    def test_subs_and_interstits(self):
        # latt = self.gaas_struct.lattice
        # pregen_intersite = [PeriodicSite('As', [0.4750, 0.4750, 0.5250], latt),
        #                     PeriodicSite('As', [0.7250, 0.7250, 0.7750], latt)]
        # CDS = ChargedDefectsStructures(self.gaas_struct,
        #                                substitutions={'Ga':['Si','In'], 'As':['Sb']},
        #                                intersites= pregen_intersite,
        #                                include_interstitials=True)
        CDS = ChargedDefectsStructures(self.gaas_struct,
                                       substitutions={'Ga':['Si','In'], 'As':['Sb']})
        self.assertEqual('sub_1_Si_on_Ga', CDS.defects['substitutions'][0]['name'])
        self.assertEqual(self.ga_site, CDS.defects['substitutions'][0]['unique_site'])
        self.assertEqual('sub_1_In_on_Ga', CDS.defects['substitutions'][1]['name'])
        self.assertEqual(self.ga_site, CDS.defects['substitutions'][1]['unique_site'])
        self.assertEqual('sub_2_Sb_on_As', CDS.defects['substitutions'][2]['name'])
        self.assertEqual(self.as_site, CDS.defects['substitutions'][2]['unique_site'])

        #NEXT do intersites


if __name__ == '__main__':
    import unittest
    unittest.main()
