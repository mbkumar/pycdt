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
import unittest
from shutil import copyfile

from monty.tempfile import ScratchDir
from pymatgen.core import Composition, Element
from pycdt.core.chemical_potentials import ChemPotAnalyzer, MPChemPotAnalyzer, \
    UserChemPotAnalyzer, UserChemPotInputGenerator

from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry

file_loc = os.path.abspath(os.path.join(
    __file__, '..', '..', '..', '..', 'test_files'))

class ChemPotAnalyzerTest(unittest.TestCase):
    def setUp(self):
        self.CPA = ChemPotAnalyzer()

    def test_get_chempots_from_pda(self):
        with MPRester() as mp:
            bulk_ce = mp.get_entry_by_material_id('mp-2534')
            bulk_species_symbol = [s.symbol for s in bulk_ce.composition.elements]
            entries = mp.get_entries_in_chemsys(bulk_species_symbol)

        #intrinsic chempot test
        pd = PhaseDiagram(entries)
        self.CPA.bulk_ce = bulk_ce
        gaas_cp = self.CPA.get_chempots_from_pd(pd)
        self.assertEqual(set([u'GaAs-As', u'GaAs-Ga']), set(gaas_cp.keys()))
        self.assertEqual([ -4.6580705550000001, -3.7317319750000006],
                               [ gaas_cp['GaAs-As'][Element('As')],
                                 gaas_cp['GaAs-As'][Element('Ga')]] )
        self.assertEqual([-5.352569090000001, -3.03723344],
                               [ gaas_cp['GaAs-Ga'][Element('As')],
                                 gaas_cp['GaAs-Ga'][Element('Ga')]] )

    def test_diff_bulk_sub_phases(self):
        fl = ['GaAs', 'Sb', 'GaSb', 'Ga']
        blk, blknom, subnom = self.CPA.diff_bulk_sub_phases(fl, 'Sb')
        self.assertEqual( ['Ga', 'GaAs'], blk)
        self.assertEqual( 'Ga-GaAs', blknom)
        self.assertEqual( 'GaSb-Sb', subnom)


class MPChemPotAnalyzerTest(unittest.TestCase):
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
        self.assertEqual(set([u'GaSb-Ga-In-GaAs', u'InAs-GaSb-In-GaAs',
                              u'InAs-GaSb-Sb-GaAs', u'InAs-As-Sb-GaAs']), set(cp_fsaf.keys()))
        true_answer = [cp_fsaf['GaSb-Ga-In-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.44626405, -5.352569090000001, -3.03723344, -2.72488125],
                               true_answer)
        true_answer = [cp_fsaf['InAs-GaSb-In-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.243837989999999, -5.15014303, -3.239659500000001, -2.72488125],
                               true_answer)
        true_answer = [cp_fsaf['InAs-GaSb-Sb-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.127761275, -5.0340663150000005, -3.3557362150000003, -2.8409579649999994],
                               true_answer)
        true_answer = [cp_fsaf['InAs-As-Sb-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.127761275, -4.658070555, -3.7317319750000006, -3.2169537249999998],
                               true_answer)

        #test (v) one example of full_sub_approach = True... (just doing for case with user's entry being stable)
        cp_fsat = self.MPCPA.analyze_GGA_chempots(full_sub_approach=True)
        self.assertEqual(set(['GaSb-Ga-In-GaAs','InAs-As-Sb-GaAs', 'InSb-InAs-In-GaAs', 'InSb-GaSb-In-GaAs',
                              'InSb-InAs-Sb-GaAs', 'InSb-GaSb-Sb-GaAs']), set(cp_fsat.keys()))
        true_answer = [cp_fsat['GaSb-Ga-In-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.44626405, -5.352569090000001, -3.03723344, -2.72488125],
                               true_answer)
        true_answer = [cp_fsat['InAs-As-Sb-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.127761275, -4.658070555, -3.7317319750000006, -3.2169537249999998],
                         true_answer)
        true_answer = [cp_fsat['InSb-InAs-In-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.3819529699999995, -5.15014303, -3.239659500000001, -2.72488125],
                         true_answer)
        true_answer = [cp_fsat['InSb-GaSb-In-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.3819529699999995, -5.28825801, -3.101544520000001, -2.72488125],
                         true_answer)
        true_answer = [cp_fsat['InSb-InAs-Sb-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.127761275, -4.895951335, -3.4938511950000004, -2.9790729449999995],
                         true_answer)
        true_answer = [cp_fsat['InSb-GaSb-Sb-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga', 'In']]
        self.assertEqual([-4.127761275, -5.0340663150000005, -3.3557362150000003, -2.9790729449999995],
                         true_answer)

        #test (iii) not stable, composition exists in list of stable comps of PD; full_sub_approach = False
        us_bce=ComputedEntry(Composition({'Ga':1, 'As':1}), -8)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=us_bce)
        cp_fsaf_us_ce = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set(['As-GaAs', 'Ga-GaAs']), set(cp_fsaf_us_ce.keys()))
        self.assertEqual([-4.6580705550000001, -3.7317319750000006],
                         [cp_fsaf_us_ce['As-GaAs'][Element(elt)] for elt in ['As', 'Ga']])
        self.assertEqual([-5.352569090000001, -3.03723344],
                         [cp_fsaf_us_ce['Ga-GaAs'][Element(elt)] for elt in ['As', 'Ga']])

        #test (ii) user's computed entry is stable w.r.t MP phase diagram
        #       AND  composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
        s_cdne_bce = ComputedEntry(Composition({'Ga':2, 'As':3}), -21.5)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=s_cdne_bce)
        cp_fsaf_s_cdne = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set(['Ga2As3-GaAs', 'As-Ga2As3']), set(cp_fsaf_s_cdne.keys()))
        for tested, answer in zip([cp_fsaf_s_cdne['Ga2As3-GaAs'][Element(elt)] for elt in ['As', 'Ga']],
                                  [-4.720394939999998, -3.669407590000002]):
            self.assertAlmostEqual( answer,tested)
        self.assertEqual([-4.6580705550000001, -3.7628941674999994],
                         [cp_fsaf_s_cdne['As-Ga2As3'][Element(elt)] for elt in ['As', 'Ga']])

        #test (iv) not stable, composition DOESNT exists in list of stable comps of PD; full_sub_approach = False
        #       case a) simple 2D phase diagram
        us_cdne_bce_a = ComputedEntry(Composition({'Ga':2, 'As':3}), -20.)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=us_cdne_bce_a)
        cp_fsaf_us_cdne_a = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set(['GaAs-As']), set(cp_fsaf_us_cdne_a.keys()))
        self.assertEqual([-4.658070555, -3.7317319750000006],
                         [cp_fsaf_us_cdne_a['GaAs-As'][Element(elt)] for elt in ['As', 'Ga']])

        #       case b) larger phase diagram
        us_cdne_bce_b = ComputedEntry(Composition({'Ga':2, 'As':3, 'Sb':1}), -20.)
        self.MPCPA = MPChemPotAnalyzer( bulk_ce=us_cdne_bce_b)
        cp_fsaf_us_cdne_b = self.MPCPA.analyze_GGA_chempots(full_sub_approach=False)
        self.assertEqual(set(['As-Sb-GaAs']), set(cp_fsaf_us_cdne_b.keys()))
        self.assertEqual([-4.127761275, -4.658070555, -3.7317319750000006],
                         [cp_fsaf_us_cdne_b['As-Sb-GaAs'][Element(elt)] for elt in ['Sb', 'As', 'Ga']])


    def test_get_chempots_from_composition(self):
        bulk_comp = Composition("Cr2O3")
        self.MPCPA = MPChemPotAnalyzer() #reinitalize
        cro_cp = self.MPCPA.get_chempots_from_composition(bulk_comp)
        self.assertEqual(set(['CrO2-Cr2O3', 'Cr2O3-Cr']), set(cro_cp.keys()))
        self.assertAlmostEqual(-14.635303979999982, cro_cp['CrO2-Cr2O3'][Element('Cr')])
        self.assertAlmostEqual(-5.51908629500001, cro_cp['CrO2-Cr2O3'][Element('O')])
        self.assertAlmostEqual(-9.63670386, cro_cp['Cr2O3-Cr'][Element('Cr')])
        self.assertAlmostEqual(-8.851486374999999, cro_cp['Cr2O3-Cr'][Element('O')])

    def test_get_mp_entries(self):
        #name mp-id of GaAs system and get mp-entries...
        # NOTE this test will start breaking
        #   if more mp entries added to this phase diagram in future.
        self.MPCPA = MPChemPotAnalyzer(mpid= 'mp-2534')
        self.MPCPA.get_mp_entries()
        ents = self.MPCPA.entries
        self.assertEqual(set(['bulk_derived', 'subs_set']), set(ents.keys()))
        self.assertEqual(18, len(ents['bulk_derived']))


class UserChemPotAnalyzerTest(unittest.TestCase):
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
            self.assertEqual(set(['As-GaAs', 'Ga-GaAs']), set(cp.keys()))
            self.assertEqual([ -5.3585191833333328, -4.2880321141666675],
                             [cp['As-GaAs'][Element(elt)] for elt in ['As', 'Ga']])
            self.assertEqual([-6.0386246900000007, -3.6079266075],
                             [cp['Ga-GaAs'][Element(elt)] for elt in ['As', 'Ga']])

        # followed by an case where MP needs to supplement...
        with ScratchDir('.'):
            os.mkdir('PhaseDiagram')
            #NO Ga entry included this time
            os.mkdir('PhaseDiagram/As')
            copyfile(os.path.join(file_loc, 'vasprun.xml_As'), 'PhaseDiagram/As/vasprun.xml')
            os.mkdir('PhaseDiagram/GaAs')
            copyfile(os.path.join(file_loc, 'vasprun.xml_GaAs'), 'PhaseDiagram/GaAs/vasprun.xml')
            cp = self.UCPA.read_phase_diagram_and_chempots(full_sub_approach=False, include_mp_entries=True)
            self.assertEqual(set(['As-GaAs', 'Ga-GaAs']), set(cp.keys()))
            self.assertEqual([-5.3585191833333328, -4.2880321141666675],
                             [cp['As-GaAs'][Element(elt)] for elt in ['As', 'Ga']])
            self.assertEqual([-6.609317857500001, -3.03723344],
                             [cp['Ga-GaAs'][Element(elt)] for elt in ['As', 'Ga']])


class UserChemPotInputGeneratorTest(unittest.TestCase):
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
