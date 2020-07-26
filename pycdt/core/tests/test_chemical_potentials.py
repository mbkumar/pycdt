# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "Jan 14, 2017"

import os
import copy
import inspect
import unittest
from shutil import copyfile

from monty.serialization import loadfn
from monty.tempfile import ScratchDir

from pycdt.core.chemical_potentials import ChemPotAnalyzer, MPChemPotAnalyzer, \
    UserChemPotAnalyzer, UserChemPotInputGenerator, get_mp_chempots_from_dpd

from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.defects.thermodynamics import DefectPhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.ext.matproj import MPRester
from pymatgen.util.testing import PymatgenTest

TEST_DIR = os.path.abspath(os.path.join(
    __file__, '..', '..', '..', '..', 'test_files'))

# class DpdAnalyzerTest(PymatgenTest):
#     def setUp(self):
#         dpd_path = os.path.split(inspect.getfile(DefectPhaseDiagram))[0]
#         pymatgen_test_files = os.path.abspath(os.path.join(dpd_path, 'tests'))
#         vbm_val = 2.6682
#         gap = 1.5
#         entries = list(loadfn(os.path.join(pymatgen_test_files, "GaAs_test_defentries.json")).values())
#         for entry in entries:
#             entry.parameters.update({'vbm': vbm_val})
#
#         self.dpd = DefectPhaseDiagram(entries, vbm_val, gap)
#
#     def test_get_mp_chempots_from_dpd(self):
#         cps = get_mp_chempots_from_dpd(self.dpd)
#         self.assertEqual(set([u'As-GaAs', u'Ga-GaAs']), set(cps.keys()))
#         self.assertEqual([ -4.66, -3.73],
#                                [ round(cps['As-GaAs'][Element('As')],2),
#                                  round(cps['As-GaAs'][Element('Ga')],2)] )
#         self.assertEqual([-5.36, -3.03],
#                                [ round(cps['Ga-GaAs'][Element('As')],2),
#                                  round(cps['Ga-GaAs'][Element('Ga')],2)] )
#
#
#
# class ChemPotAnalyzerTest(PymatgenTest):
#     def setUp(self):
#         self.CPA = ChemPotAnalyzer()
#
#     def test_get_chempots_from_pda(self):
#         with MPRester() as mp:
#             bulk_ce = mp.get_entry_by_material_id('mp-2534')
#             bulk_species_symbol = [s.symbol for s in bulk_ce.composition.elements]
#             entries = mp.get_entries_in_chemsys(bulk_species_symbol)
#
#         #intrinsic chempot test
#         pd = PhaseDiagram(entries)
#         self.CPA.bulk_ce = bulk_ce
#         gaas_cp = self.CPA.get_chempots_from_pd(pd)
#         self.assertEqual(set([u'As-GaAs', u'Ga-GaAs']), set(gaas_cp.keys()))
#         self.assertEqual([ -4.66, -3.73],
#                                [ round(gaas_cp['As-GaAs'][Element('As')],2),
#                                  round(gaas_cp['As-GaAs'][Element('Ga')],2)] )
#         self.assertEqual([-5.36, -3.03],
#                                [ round(gaas_cp['Ga-GaAs'][Element('As')],2),
#                                  round(gaas_cp['Ga-GaAs'][Element('Ga')],2)] )
#
#     def test_diff_bulk_sub_phases(self):
#         fl = ['GaAs', 'Sb', 'GaSb', 'Ga']
#         blk, blknom, subnom = self.CPA.diff_bulk_sub_phases(fl, 'Sb')
#         self.assertEqual(['Ga', 'GaAs'], blk)
#         self.assertEqual('Ga-GaAs', blknom)
#         self.assertEqual('GaSb-Sb', subnom)
#
#
# class MPChemPotAnalyzerTest(PymatgenTest):
#     def setUp(self):
#         self.MPCPA = MPChemPotAnalyzer()
#
#     def test_analyze_GGA_chempots(self):
#         """
#         Test-cases to cover:
#             (0) (stable + mpid given = tested in the CPA.test_get_chempots_from_pda unit test)
#             (i) user's computed entry is stable w.r.t MP phase diagram; full_sub_approach = False
#             (ii) user's computed entry is stable w.r.t MP phase diagram (and not currently in phase diagram);
#                 full_sub_approach = False
#             (iii) not stable, composition exists in list of stable comps of PD; full_sub_approach = False
#             (iv) not stable, composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
#             (v) one example of full_sub_approach = True... (just do for case with user's entry being stable)
#         """
#         with MPRester() as mp:
#             bulk_ce = mp.get_entry_by_material_id('mp-2534') #simulates an entry to play with for this unit test
#
#         #test (i) user's computed entry is stable w.r.t MP phase diagram; full_sub_approach = False
#         bce = copy.copy(bulk_ce)
#         self.MPCPA = MPChemPotAnalyzer(bulk_ce=bce, sub_species=set(['Sb', 'In']))
#         cp_fsaf = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
#         self.assertEqual(set([u'Ga-GaAs-GaSb-In',u'As-GaAs-InAs-SbAs']), set(cp_fsaf.keys()))
#         true_answer = [round(cp_fsaf['Ga-GaAs-GaSb-In'][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.44, -5.36, -3.03, -2.75],
#                                true_answer)
#         true_answer = [round(cp_fsaf['As-GaAs-InAs-SbAs'][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.16, -4.66, -3.73, -3.21],
#                                true_answer)
#
#         #test (v) one example of full_sub_approach = True... (just doing for case with user's entry being stable)
#         cp_fsat = self.MPCPA.analyze_GGA_chempots(full_sub_approach=True)
#         self.assertEqual(set(['Ga-GaAs-GaSb-In', 'GaAs-InAs-InSb-Sb', 'GaAs-In-InAs-InSb',
#                               'As-GaAs-InAs-SbAs', 'GaAs-GaSb-In-InSb', 'GaAs-InAs-Sb-SbAs',
#                               'GaAs-GaSb-InSb-Sb']), set(cp_fsat.keys()))
#         true_answer = [round(cp_fsat["Ga-GaAs-GaSb-In"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.44, -5.36, -3.03, -2.75],
#                          true_answer)
#         true_answer = [round(cp_fsat["GaAs-InAs-InSb-Sb"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.13, -4.91, -3.48, -2.96],
#                          true_answer)
#         true_answer = [round(cp_fsat["GaAs-In-InAs-InSb"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.33, -5.12, -3.27, -2.75],
#                          true_answer)
#         true_answer = [round(cp_fsat["As-GaAs-InAs-SbAs"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.16, -4.66, -3.73, -3.21],
#                          true_answer)
#         true_answer = [round(cp_fsat["GaAs-GaSb-In-InSb"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.33, -5.25, -3.14, -2.75],
#                          true_answer)
#         true_answer = [round(cp_fsat["GaAs-InAs-Sb-SbAs"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.13, -4.69, -3.7, -3.18],
#                          true_answer)
#         true_answer = [round(cp_fsat["GaAs-GaSb-InSb-Sb"][Element(elt)],2) for elt in ['Sb', 'As', 'Ga', 'In']]
#         self.assertEqual([-4.13, -5.05, -3.34, -2.96],
#                          true_answer)
#
#         #test (iii) not stable, composition exists in list of stable comps of PD; full_sub_approach = False
#         us_bce=ComputedEntry(Composition({'Ga':1, 'As':1}), -8)
#         self.MPCPA = MPChemPotAnalyzer(bulk_ce=us_bce)
#         cp_fsaf_us_ce = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
#         self.assertEqual(set(['As-GaAs', 'Ga-GaAs']), set(cp_fsaf_us_ce.keys()))
#         self.assertEqual([-4.66, -3.73],
#                          [round(cp_fsaf_us_ce['As-GaAs'][Element(elt)],2) for elt in ['As', 'Ga']])
#         self.assertEqual([-5.36, -3.03],
#                          [round(cp_fsaf_us_ce['Ga-GaAs'][Element(elt)],2) for elt in ['As', 'Ga']])
#
#         #test (ii) user's computed entry is stable w.r.t MP phase diagram
#         #       AND  composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
#         s_cdne_bce = ComputedEntry(Composition({'Ga':2, 'As':3}), -21.5)
#         self.MPCPA = MPChemPotAnalyzer(bulk_ce=s_cdne_bce)
#         cp_fsaf_s_cdne = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
#         self.assertEqual(set(['Ga2As3-GaAs', 'As-Ga2As3']), set(cp_fsaf_s_cdne.keys()))
#         for tested, answer in zip([round(cp_fsaf_s_cdne['Ga2As3-GaAs'][Element(elt)],2) for elt in ['As', 'Ga']],
#                                   [-4.73, -3.66]):
#             self.assertAlmostEqual(answer,tested)
#         self.assertEqual([-4.66, -3.76],
#                          [round(cp_fsaf_s_cdne['As-Ga2As3'][Element(elt)],2) for elt in ['As', 'Ga']])
#
#         #test (iv) not stable, composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
#         #       case a) simple 2D phase diagram
#         us_cdne_bce_a = ComputedEntry(Composition({'Ga':2, 'As':3}), -20.)
#         self.MPCPA = MPChemPotAnalyzer(bulk_ce=us_cdne_bce_a)
#         cp_fsaf_us_cdne_a = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
#         self.assertEqual(set(['As-GaAs']), set(cp_fsaf_us_cdne_a.keys()))
#         self.assertEqual([-4.66, -3.73],
#                          [round(cp_fsaf_us_cdne_a['As-GaAs'][Element(elt)],2) for elt in ['As', 'Ga']])
#
#         #       case b) larger phase diagram
#         us_cdne_bce_b = ComputedEntry(Composition({'Ga':2, 'As':3, 'Sb':2}), -20.)
#         self.MPCPA = MPChemPotAnalyzer(bulk_ce=us_cdne_bce_b)
#         cp_fsaf_us_cdne_b = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
#         self.assertEqual(set(['GaAs-Sb-SbAs']), set(cp_fsaf_us_cdne_b.keys()))
#         self.assertEqual([-4.13, -4.69, -3.70],
#                          [round(cp_fsaf_us_cdne_b['GaAs-Sb-SbAs'][Element(elt)],2) for elt in ['Sb', 'As', 'Ga']])
#
#
#     def test_get_chempots_from_composition(self):
#         bulk_comp = Composition("Cr2O3")
#         self.MPCPA = MPChemPotAnalyzer() #reinitalize
#         cro_cp = self.MPCPA.get_chempots_from_composition(bulk_comp)
#         self.assertEqual(set(['Cr2O3-CrO2', 'Cr-Cr2O3']), set(cro_cp.keys()))
#         self.assertAlmostEqual(-14.9, round(cro_cp['Cr2O3-CrO2'][Element('Cr')],2))
#         self.assertAlmostEqual(-5.42, round(cro_cp['Cr2O3-CrO2'][Element('O')],2))
#         self.assertAlmostEqual(-9.65, round(cro_cp['Cr-Cr2O3'][Element('Cr')],2))
#         self.assertAlmostEqual(-8.92, round(cro_cp['Cr-Cr2O3'][Element('O')],2))
#
#     def test_get_mp_entries(self):
#         #name mp-id of GaAs system and get mp-entries...
#         self.MPCPA = MPChemPotAnalyzer(mpid= 'mp-2534')
#         self.MPCPA.get_mp_entries()
#         ents = self.MPCPA.entries
#         self.assertEqual(set(['bulk_derived', 'subs_set']), set(ents.keys()))
#         self.assertTrue(len(ents['bulk_derived']))


class UserChemPotAnalyzerTest(PymatgenTest):
    def setUp(self):
        with MPRester() as mp:
            self.bulk_ce = mp.get_entry_by_material_id('mp-2534')
        self.UCPA = UserChemPotAnalyzer(bulk_ce = self.bulk_ce)
        self.UCPA_sub = UserChemPotAnalyzer(bulk_ce = self.bulk_ce, sub_species=["In"])

    def test_read_phase_diagram_and_chempots(self):
        #set up a local phase diagram object...
        # test non mp case,
        with ScratchDir('.'):
            #os.mkdir('PhaseDiagram')
            os.makedirs(os.path.join("PhaseDiagram", "Ga"))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_Ga'),
                     os.path.join("PhaseDiagram", "Ga", "vasprun.xml"))
            os.mkdir(os.path.join('PhaseDiagram', 'As'))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_As'),
                     os.path.join("PhaseDiagram", "As", "vasprun.xml"))
            os.mkdir(os.path.join("PhaseDiagram", "GaAs"))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_GaAs'),
                     os.path.join("PhaseDiagram", "GaAs", "vasprun.xml"))
            cp = self.UCPA.read_phase_diagram_and_chempots(
                    full_sub_approach=False, include_mp_entries=False)
            self.assertEqual(set(['As-GaAs', 'Ga-GaAs']), set(cp.keys()))
            self.assertEqual(
                    [-5.36, -4.29],
                    [round(cp['As-GaAs'][Element(elt)], 2) for elt in ['As', 'Ga']])
            self.assertEqual(
                    [-6.04, -3.61],
                    [round(cp['Ga-GaAs'][Element(elt)], 2) for elt in ['As', 'Ga']])

        # followed by an case where MP needs to supplement...
        with ScratchDir('.'):
            os.mkdir('PhaseDiagram')
            #NO Ga entry included this time
            os.mkdir(os.path.join("PhaseDiagram", "As"))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_As'), 
                     os.path.join("PhaseDiagram", "As", "vasprun.xml"))
            os.mkdir(os.path.join("PhaseDiagram", "GaAs"))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_GaAs'),
                     os.path.join("PhaseDiagram", "GaAs", "vasprun.xml"))
            cp = self.UCPA.read_phase_diagram_and_chempots(
                    full_sub_approach=False, include_mp_entries=True)
            self.assertEqual(set(['As-GaAs', 'Ga-GaAs']), set(cp.keys()))
            self.assertEqual(
                    [-5.36, -4.29],
                    [round(cp['As-GaAs'][Element(elt)], 2) for elt in ['As', 'Ga']])
            self.assertEqual(
                    [-6.62, -3.03],
                    [round(cp['Ga-GaAs'][Element(elt)], 2) for elt in ['As', 'Ga']])

        #quick and dirty test for finding extrinsic defects...
        with ScratchDir('.'):
            os.mkdir('PhaseDiagram')
            #NO Ga entry or In entry this time
            os.mkdir(os.path.join("PhaseDiagram", "As"))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_As'),
                     os.path.join("PhaseDiagram", "As", "vasprun.xml"))
            os.mkdir(os.path.join("PhaseDiagram", "GaAs"))
            copyfile(os.path.join(TEST_DIR, 'vasprun.xml_GaAs'),
                     os.path.join("PhaseDiagram", "GaAs", "vasprun.xml"))
            cp = self.UCPA_sub.read_phase_diagram_and_chempots(
                    full_sub_approach=False, include_mp_entries=True)
            self.assertEqual(set(['As-GaAs-In', 'Ga-GaAs-In']), set(cp.keys()))

        #TODO: add a vasprun.xml with extrinsic specie to verify user-supplied extrinsic defect


class UserChemPotInputGeneratorTest(PymatgenTest):
    def setUp(self):
        self.UCPIGT = UserChemPotInputGenerator(Composition({'Ga': 1, 'As': 1}))

    def test_setup_phase_diagram_calculations(self):
        # test that files get craeted for some file of interest...
        with ScratchDir('.'):
            self.UCPIGT.setup_phase_diagram_calculations(full_phase_diagram=True)
            # self.assertTrue(os.path.exists('PhaseDiagram')) # redundant
            # self.assertTrue(os.path.exists('PhaseDiagram/mp-11_As')) # redundant
            self.assertTrue(os.path.exists(
                    os.path.join("PhaseDiagram", "mp-11_As", "POSCAR")))
            # self.assertTrue(os.path.exists('PhaseDiagram/mp-142_Ga')) # redundant
            self.assertTrue(os.path.exists(
                    os.path.join("PhaseDiagram", "mp-142_Ga", "POSCAR")))
            # self.assertTrue(os.path.exists('PhaseDiagram/mp-2534_GaAs/')) # redundant
            self.assertTrue(os.path.exists(
                    os.path.join("PhaseDiagram", "mp-2534_GaAs", "POSCAR")))


if __name__ == '__main__':
    unittest.main()
