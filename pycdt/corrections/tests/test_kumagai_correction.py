# coding: utf-8

from __future__ import division

__author__ = "Danny Broberg, Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "June 23, 2016"

import unittest
import os
import numpy as np

from pymatgen.io.vasp.outputs import Locpot
from pycdt.corrections.kumagai_correction import *

# Paths to files we are testing on
# bl_path = 'bLOCPOT.gz'
# dl_path = 'dLOCPOT.gz'
# kad_path = 'testKumagaiData.json'
bl_path = os.path.join('..', '..', '..', 'test_files', 'bLOCPOT.gz')
dl_path = os.path.join('..', '..', '..', 'test_files', 'dLOCPOT.gz')
kad_path = os.path.join('..', '..', '..', 'test_files', 'testKumagaiData.json')

class KumagaiBulkInitANDCorrectionTest(unittest.TestCase):
    #TODO: also might want to test outcar Kumagai method...
    def setUp(self):
        self.bl = Locpot.from_file(bl_path)
        self.dl = Locpot.from_file(dl_path)
        self.bs = self.bl.structure
        self.ds = self.dl.structure
        self.kbi = KumagaiBulkInit(self.bs, self.bl.dim, 15,
                                   optgamma=3.49423226983)
        self.kc = KumagaiCorrection(15, -3, 3.49423226983, self.kbi.g_sum,
                                    self.bs, self.ds, bulk_locpot=self.bl,
                                    defect_locpot=self.dl)

    def test_find_optimal_gamma(self):
        self.assertEqual(self.kbi.find_optimal_gamma(), 3.4942322698305639)

    def test_reciprocal_sum(self):
        # This is initialized with KBI.
        # Not sure how else to test besides size of list of vectors
        self.assertEqual(self.kbi.g_sum.size, 884736)
        self.assertAlmostEqual(self.kbi.g_sum[0][0][0], 0.050661706751775192)

    def test_pc(self):
        self.assertAlmostEqual(self.kc.pc(), 2.1315841582145407)

    def test_potalign(self):
        self.assertAlmostEqual(self.kc.potalign(), 2.1091426308966001)

    def test_correction(self):
        self.assertAlmostEqual(self.kc.correction(), 4.24073)

    def test_plot(self):
        tmpforplot = {'C': {'r': [1, 2, 3], 'Vqb': [0.1, 0.2, 0.3],
                            'Vpc': [-0.05, -0.15, -0.25]},
                      'EXTRA': {'wsrad': 1, 'potalign': 0.05,
                                'lengths': (3, 3, 3)}}
        KumagaiCorrection.plot(tmpforplot, 'TMP')
        self.assertTrue(os.path.exists('TMP_kumagaisiteavgPlot.pdf'))
        os.system('rm TMP_kumagaisiteavgPlot.pdf')

    def test_plot_from_datfile(self):
        KumagaiCorrection.plot_from_datfile(name=kad_path, title='TMP')
        self.assertTrue(os.path.exists('TMP_kumagaisiteavgPlot.pdf'))
        os.system('rm TMP_kumagaisiteavgPlot.pdf')

    #NOTE there are here because g_sum is pre computed for it
    def test_get_sum_at_r(self):
        val = get_g_sum_at_r(self.kbi.g_sum, self.bs, self.bl.dim,
                             [0.1, 0.1, 0.1])
        self.assertAlmostEqual(val, 0.04795055159361078)

    def test_anisotropic_madelung_potential(self):
        val = anisotropic_madelung_potential(
                self.bs, self.bl.dim, self.kbi.g_sum, [0.1, 0.1, 0.1],
                [[15, 0.1, -0.1], [0.1, 13, 0], [-0.1, 0, 20]], -3,
                self.kbi.gamma, self.kbi.tolerance)
        self.assertAlmostEqual(val, -4.2923511216202419)

    def test_anisotropic_pc_energy(self):
        val = anisotropic_pc_energy(
                self.bs, self.kbi.g_sum,
                [[15, 0.1, -0.1], [0.1, 13, 0], [-0.1, 0, 20]], -3,
                self.kbi.gamma, self.kbi.tolerance)
        self.assertAlmostEqual(val, 1.5523329679084736)


class KumagaiSetupFunctionsTest(unittest.TestCase):
    def setUp(self):
        self.bl = Locpot.from_file(bl_path)
        self.dl = Locpot.from_file(dl_path)
        self.bs = self.bl.structure
        self.ds = self.dl.structure

    def test_kumagai_init(self):
        angset, bohrset, vol, determ, invdiel = kumagai_init(
                self.bs, [[15, 0.1, -0.1], [0.1, 13, 0], [-0.1, 0, 20]])
        newangset = [list(row) for row in angset]
        newbohrset = [list(row) for row in bohrset]
        newinvdiel = [list(row) for row in invdiel]
        self.assertEqual(newangset, [[5.750183, 0, 0], [0, 5.750183, 0],
                                  [0, 0, 5.750183]])
        self.assertEqual(newbohrset, [[10.866120815099999, 0, 0],
                                      [0, 10.866120815099999, 0],
                                      [0, 0, 10.866120815099999]])
        self.assertAlmostEqual(vol, 1282.9909362724345)
        self.assertAlmostEqual(determ, 3899.6699999999969)
        tmpinvdiel = [[0.066672308169665628, -0.00051286390899742801, 0.00033336154084832821],
                      [-0.00051286390899742801, 0.076927022030069209, -0.0000025643195449871406],
                      [0.00033336154084832826, -0.0000025643195449871406, 0.050001666807704244]]
        self.assertAlmostEqual(newinvdiel, tmpinvdiel)

    def test_real_sum(self):
        a = self.bs.lattice.matrix[0]
        b = self.bs.lattice.matrix[1]
        c = self.bs.lattice.matrix[2]
        tmpdiel = [[15, 0.1, -0.1], [0.1, 13, 0], [-0.1, 0, 20]]
        val = real_sum(a, b, c, np.array([0.1, 0.1, 0.1]), -1, tmpdiel, 3, 1)
        self.assertAlmostEqual(val, -0.0049704211394050414)

    def test_disttrans(self):
        #not sure what best way to test dictionary is...
        pass

    def test_wigner_seitz_radius(self):
        self.assertAlmostEqual(wigner_seitz_radius(self.bs), 2.8750914999999999)

    def test_read_ES_avg(self):
        #not sure what best way to test dictionary is...
        pass

    def test_read_ES_avg_fromlocpot(self):
        #not sure what best way to test dictionary is...
        pass


class EnergyFunctionsTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_k_to_eV(self):
        g = [0.1, 0.2, 0.3]
        self.assertAlmostEqual(k_to_eV(g), 0.5333804)

    def test_eV_to_k(self):
        self.assertAlmostEqual(eV_to_k(1), 0.9681404248678961)


if __name__ == '__main__':
    unittest.main()
