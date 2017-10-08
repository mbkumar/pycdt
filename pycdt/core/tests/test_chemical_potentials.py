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
from shutil import copyfile

from monty.tempfile import ScratchDir
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.composition import Composition
from pycdt.core.chemical_potentials import *

from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
#from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.entries.computed_entries import ComputedEntry

file_loc = os.path.abspath(os.path.join('..', '..', '..', 'test_files'))

class ChemPotAnalyzerTest(PymatgenTest):
    def setUp(self):
        self.CPA = ChemPotAnalyzer()

    def test_get_chempots_from_pda(self):
        with MPRester() as mp:
            bulk_ce = mp.get_entry_by_material_id('mp-2534')
            bulk_species_symbol = [s.symbol for s in bulk_ce.composition.elements]
            entries = mp.get_entries_in_chemsys(bulk_species_symbol)

        #intrinsic chempot test
        pd = PhaseDiagram(entries)
        pda = pd # PDAnalyzer(pd)
        self.CPA.bulk_ce = bulk_ce
        gaas_cp = self.CPA.get_chempots_from_pda(pda)
        self.assertEqual([u'As_rich', u'Ga_rich'], gaas_cp.keys())
        self.assertEqual({u'As': -4.6580705550000001, u'Ga': -3.7317319750000006}, gaas_cp['As_rich'])
        self.assertEqual({u'As': -5.3589730550000008, u'Ga': -3.030829475}, gaas_cp['Ga_rich'])

    def test_diff_bulk_sub_phases(self):
        fl = ['GaAs', 'Sb', 'GaSb', 'Ga']
        blk, blknom, subnom = self.CPA.diff_bulk_sub_phases(fl, 'Sb')
        self.assertEqual( ['Ga', 'GaAs'], blk)
        self.assertEqual( 'Ga-GaAs', blknom)
        self.assertEqual( 'GaSb-Sb', subnom)


class MPChemPotAnalyzerTest(PymatgenTest):
    def setUp(self):
        self.MPCPA = MPChemPotAnalyzer()

    def test_analyze_GGA_chempots(self):
        """
        Test-cases to cover:
            (0) (stable + mpid given = tested in the CPA.test_get_chempots_from_pda unit test)
            (i) user's computed entry is stable w.r.t MP phase diagram; full_sub_approach = False
            (ii) user's computed entry is stable w.r.t MP phase diagram (and not currently in phase diagram);
                full_sub_approach = False
            (iii) not stable, composition exists in list of stable comps of PD; full_sub_approach = False
            (iv) not stable, composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
            (v) one example of full_sub_approach = True... (just do for case with user's entry being stable)
        """
        with MPRester() as mp:
            bulk_ce = mp.get_entry_by_material_id('mp-2534') #simulates an entry to play with for this unit test

        #test (i) user's computed entry is stable w.r.t MP phase diagram; full_sub_approach = False
        bce = copy.copy(bulk_ce)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=bce, sub_species=set(['Sb', 'In']))
        cp_fsaf = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set([u'GaSb-Ga-In', u'InAs-GaSb-In', u'InAs-GaSb-Sb', u'InAs-As-Sb']), set(cp_fsaf.keys()))
        self.assertAlmostEqual({u'Sb': -4.4526680150000004, u'As': -5.3589730550000008, u'Ga': -3.030829475,
                                u'In': -2.7215626899999998}, cp_fsaf['GaSb-Ga-In'])
        self.assertAlmostEqual({u'Sb': -4.2476110999999994, u'As': -5.1539161399999998, u'Ga': -3.235886390000001,
                                u'In': -2.7215626899999998},  cp_fsaf['InAs-GaSb-In'])
        self.assertAlmostEqual({u'Sb': -4.1243295800000004, u'As': -5.0306346200000007, u'Ga': -3.35916791,
                                u'In': -2.8448442099999989}, cp_fsaf['InAs-GaSb-Sb'])
        self.assertAlmostEqual({u'Sb': -4.1243295800000004, u'As': -4.6580705550000001, u'Ga': -3.7317319750000006,
                                u'In': -3.2174082749999995}, cp_fsaf['InAs-As-Sb'])

        #test (v) one example of full_sub_approach = True... (just doing for case with user's entry being stable)
        cp_fsat = self.MPCPA.analyze_GGA_chempots(full_sub_approach=True)
        self.assertEqual(set(['GaSb-Ga-In','InAs-As-Sb', 'InSb-InAs-In', 'InSb-GaSb-In',
                              'InSb-InAs-Sb', 'InSb-GaSb-Sb']), set(cp_fsat.keys()))
        self.assertAlmostEqual({u'Sb': -4.4526680150000004, u'As': -5.3589730550000008, u'Ga': -3.030829475,
                                u'In': -2.7215626899999998}, cp_fsat['GaSb-Ga-In'])
        self.assertAlmostEqual({u'Sb': -4.1243295800000004, u'As': -4.6580705550000001, u'Ga': -3.7317319750000006,
                                u'In': -3.2174082749999995}, cp_fsat['InAs-As-Sb'])
        self.assertAlmostEqual({u'Sb': -4.3857206700000004, u'As': -5.1539161399999998, u'Ga': -3.235886390000001,
                                u'In': -2.7215626899999998}, cp_fsat['InSb-InAs-In'])
        self.assertAlmostEqual({u'Sb': -4.3857206700000004, u'As': -5.2920257100000008, u'Ga': -3.09777682,
                                u'In': -2.7215626899999998}, cp_fsat['InSb-GaSb-In'])
        self.assertAlmostEqual({u'Sb': -4.1243295800000004, u'As': -4.8925250499999997, u'Ga': -3.497277480000001,
                                u'In': -2.9829537799999999}, cp_fsat['InSb-InAs-Sb'])
        self.assertAlmostEqual({u'Sb': -4.1243295800000004, u'As': -5.0306346200000007, u'Ga': -3.35916791,
                                u'In': -2.9829537799999999}, cp_fsat['InSb-GaSb-Sb'])

        #test (iii) not stable, composition exists in list of stable comps of PD; full_sub_approach = False
        us_bce=ComputedEntry(Composition({'Ga':1, 'As':1}), -8)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=us_bce)
        cp_fsaf_us_ce = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set(['As_rich', 'Ga_rich']), set(cp_fsaf_us_ce.keys()))
        self.assertAlmostEqual({u'As': -4.6580705550000001, u'Ga': -3.7317319750000006}, cp_fsaf_us_ce['As_rich'])
        self.assertAlmostEqual({u'As': -5.3589730550000008, u'Ga': -3.030829475}, cp_fsaf_us_ce['Ga_rich'])

        #test (ii) user's computed entry is stable w.r.t MP phase diagram
        #       AND  composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
        s_cdne_bce = ComputedEntry(Composition({'Ga':2, 'As':3}), -21.5)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=s_cdne_bce)
        cp_fsaf_s_cdne = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set(['GaAs_rich', 'As_rich']), set(cp_fsaf_s_cdne.keys()))
        self.assertAlmostEqual(-4.7203949399999985, cp_fsaf_s_cdne['GaAs_rich']['As'])
        self.assertAlmostEqual(-3.6694075900000023, cp_fsaf_s_cdne['GaAs_rich']['Ga'])
        self.assertAlmostEqual(-4.6580705550000001, cp_fsaf_s_cdne['As_rich']['As'])
        self.assertAlmostEqual(-3.7628941674999994, cp_fsaf_s_cdne['As_rich']['Ga'])
        #self.assertAlmostEqual({u'As': -4.7203949399999985, u'Ga': -3.6694075900000023}, cp_fsaf_s_cdne['GaAs_rich'])
        #self.assertAlmostEqual({u'As': -4.6580705550000001, u'Ga': -3.7628941674999994}, cp_fsaf_s_cdne['As_rich'])

        #test (iv) not stable, composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
        #       case a) simple 2D phase diagram
        us_cdne_bce_a = ComputedEntry(Composition({'Ga':2, 'As':3}), -20.)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=us_cdne_bce_a)
        cp_fsaf_us_cdne_a = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertAlmostEqual(set(['Ga-GaAs']), set(cp_fsaf_us_cdne_a.keys()))
        self.assertAlmostEqual({u'As': -5.3589730550000008, u'Ga': -3.030829475}, cp_fsaf_us_cdne_a['Ga-GaAs'])

        #       case b) larger phase diagram
        us_cdne_bce_b = ComputedEntry(Composition({'Ga':2, 'As':3, 'Sb':1}), -20.)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=us_cdne_bce_b)
        cp_fsaf_us_cdne_b = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertAlmostEqual(set(['As-GaAs-Sb']), set(cp_fsaf_us_cdne_b.keys()))
        self.assertAlmostEqual({u'Sb': -4.1243295800000004, u'As': -4.6580705550000001,
                                u'Ga': -3.7317319750000006}, cp_fsaf_us_cdne_b['As-GaAs-Sb'])


    def test_get_chempots_from_composition(self):
        bulk_comp = Composition("Cr2O3")
        self.MPCPA = MPChemPotAnalyzer() #reinitalize
        cro_cp = self.MPCPA.get_chempots_from_composition(bulk_comp)
        self.assertEqual([u'CrO2_rich', u'Cr_rich'], cro_cp.keys())
        self.assertAlmostEqual(-14.638073584999995, cro_cp['CrO2_rich']['Cr'])
        self.assertAlmostEqual(-5.5171868950000045, cro_cp['CrO2_rich']['O'])
        self.assertAlmostEqual(-9.6385552600000004, cro_cp['Cr_rich']['Cr'])
        self.assertAlmostEqual(-8.850199111666667, cro_cp['Cr_rich']['O'])

    def test_get_mp_entries(self):
        #name mp-id of GaAs system and get mp-entries...
        # NOTE this test will start breaking
        #   if more mp entries added to this phase diagram in future??
        self.MPCPA = MPChemPotAnalyzer(mpid= 'mp-2534')
        self.MPCPA.get_mp_entries()
        ents = self.MPCPA.entries
        self.assertEqual(['bulk_derived', 'subs_set'], ents.keys())
        self.assertEqual(15, len(ents['bulk_derived']))


class UserChemPotAnalyzerTest(PymatgenTest):
    def setUp(self):
        with MPRester() as mp:
            self.bulk_ce = mp.get_entry_by_material_id('mp-2534')
        self.UCPA = UserChemPotAnalyzer(bulk_ce = self.bulk_ce)

    def test_read_phase_diagram_and_chempots(self):
        #set up a local phase diagram object...
        # test non mp case,
        with ScratchDir('.'):
            os.mkdir('PhaseDiagram')
            os.mkdir('PhaseDiagram/Ga')
            copyfile( os.path.join(file_loc, 'vasprun.xml_Ga'), 'PhaseDiagram/Ga/vasprun.xml')
            os.mkdir('PhaseDiagram/As')
            copyfile( os.path.join(file_loc, 'vasprun.xml_As'), 'PhaseDiagram/As/vasprun.xml')
            os.mkdir('PhaseDiagram/GaAs')
            copyfile( os.path.join(file_loc, 'vasprun.xml_GaAs'), 'PhaseDiagram/GaAs/vasprun.xml')
            cp = self.UCPA.read_phase_diagram_and_chempots(full_sub_approach=False, include_mp_entries=False)
            self.assertEqual(set(['As_rich', 'Ga_rich']), set(cp.keys()))
            self.assertAlmostEqual({u'As': -5.3585191833333328, u'Ga': -4.2880321141666675}, cp['As_rich'])
            self.assertAlmostEqual({u'As': -6.0386246900000007, u'Ga': -3.6079266075}, cp['Ga_rich'])

        # followed by an case where MP needs to supplement...
        with ScratchDir('.'):
            os.mkdir('PhaseDiagram')
            #NO Ga entry included this time
            os.mkdir('PhaseDiagram/As')
            copyfile(os.path.join(file_loc, 'vasprun.xml_As'), 'PhaseDiagram/As/vasprun.xml')
            os.mkdir('PhaseDiagram/GaAs')
            copyfile(os.path.join(file_loc, 'vasprun.xml_GaAs'), 'PhaseDiagram/GaAs/vasprun.xml')
            cp = self.UCPA.read_phase_diagram_and_chempots(full_sub_approach=False, include_mp_entries=True)
            self.assertEqual(set(['As_rich', 'Ga_rich']), set(cp.keys()))
            self.assertAlmostEqual({u'As': -5.3585191833333328, u'Ga': -4.2880321141666675}, cp['As_rich'])
            self.assertAlmostEqual({u'As': -6.6157218225000003, u'Ga': -3.030829475}, cp['Ga_rich'])


class UserChemPotInputGeneratorTest(PymatgenTest):
    def setUp(self):
        self.UCPIGT = UserChemPotInputGenerator(Composition({'Ga':1, 'As':1}))

    def test_setup_phase_diagram_calculations(self):
        # test that files get craeted for some file of interest...
        with ScratchDir('.'):
            self.UCPIGT.setup_phase_diagram_calculations(full_phase_diagram = True)
            self.assertTrue( os.path.exists('PhaseDiagram/'))
            self.assertTrue( os.path.exists('PhaseDiagram/mp-11_As'))
            self.assertTrue( os.path.exists('PhaseDiagram/mp-11_As/POSCAR'))
            self.assertTrue( os.path.exists('PhaseDiagram/mp-142_Ga'))
            self.assertTrue( os.path.exists('PhaseDiagram/mp-142_Ga/POSCAR'))
            self.assertTrue( os.path.exists('PhaseDiagram/mp-2534_GaAs/'))
            self.assertTrue( os.path.exists('PhaseDiagram/mp-2534_GaAs/POSCAR'))


if __name__ == '__main__':
    import unittest
    unittest.main()
