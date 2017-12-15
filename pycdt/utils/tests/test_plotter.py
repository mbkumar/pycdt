# coding: utf-8

from __future__ import division

__author__ = "Nils E. R. Zimmermann"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Nils. E. R. Zimmermann"
__email__ = "nils.e.r.zimmermann@gmail.com"
__status__ = "Development"
__date__ = "July 19, 2017"

import os

from pymatgen.core.structure import PeriodicSite, Structure, Lattice
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pycdt.utils.plotter import DefectPlotter
from pycdt.core.defects_analyzer import ComputedDefect, DefectsAnalyzer

class DefectPlotterTest(unittest.TestCase):
    def setUp(self):
        l = Lattice([[3.52,0.0,2.033], [1.174,3.32,2.033], \
                [0.0,0.0,4.066]])
        s_bulk = Structure(l, ['Ga', 'As'], \
                [[0.0000, 0.0000, 0.0000], \
                [0.2500, 0.2500, 0.2500]])
        e_bulk = 0.5
        s_vacAs = Structure(l, ['Ga'], \
                [[0.0000, 0.0000, 0.0000]])
        e_vacAs = 2.0
        entry_bulk = ComputedStructureEntry(s_bulk, e_bulk)
        entry_defect = ComputedStructureEntry(s_vacAs, e_vacAs)
        self.da = DefectsAnalyzer(entry_bulk, 0.0, {'As': 0,'Ga': 0}, 1.0)
        d_vacAs = ComputedDefect(entry_defect, s_bulk[1], multiplicity=1.0, charge=0.0, name='vac_1_As')
        self.da.add_computed_defect(d_vacAs)
        self.dp = DefectPlotter(self.da)

    def test_get_plot_form_energy(self):
        self.dp.get_plot_form_energy().savefig('test.pdf')
        self.assertTrue(os.path.exists('test.pdf'))
        os.system('rm test.pdf')

    def test_plot_conc_temp(self):
        self.dp.plot_conc_temp().savefig('test.pdf')
        self.assertTrue(os.path.exists('test.pdf'))
        os.system('rm test.pdf')

    def test_plot_carriers_ef(self):
        self.dp.plot_carriers_ef().savefig('test.pdf')
        self.assertTrue(os.path.exists('test.pdf'))
        os.system('rm test.pdf')

    def tearDown(self):
        self.da


import unittest
if __name__ == '__main__':
    unittest.main()
