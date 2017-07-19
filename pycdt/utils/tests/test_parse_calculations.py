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
from shutil import copyfile

from monty.serialization import dumpfn
from monty.json import MontyEncoder
from monty.tempfile import ScratchDir
from pymatgen import __file__ as initfilep
from pymatgen.io.vasp import Vasprun
from pycdt.utils.parse_calculations import PostProcess

pmgtestfiles_loc = os.path.join(os.path.split(os.path.split(initfilep)[0])[0], 'test_files')
# file_loc = os.path.abspath(os.path.join('..', '..', '..', 'test_files')) #Pycdt Testfiles

class PostProcessTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_parse_defect_calculations_AND_compile_all(self):
        #testing both parse defect_calculatiosn And the compile all methods because they both require a file structure...
        with ScratchDir('.'):
            #make a fake file structure to parse vaspruns and locpot paths
            os.mkdir('bulk')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml'), 'bulk/vasprun.xml')
            os.mkdir('bulk/LOCPOT') #locpot path just needs to exist..doesnt need to be real locpot file...
            bulktrans = {"supercell": [3, 3, 3], "defect_type": "bulk"}
            dumpfn(bulktrans, 'bulk/transformation.json', cls=MontyEncoder)

            os.mkdir('dielectric')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml.dfpt.ionic'), 'dielectric/vasprun.xml')

            vrobj = Vasprun(os.path.join(pmgtestfiles_loc, 'vasprun.xml'))
            lattobj = vrobj.final_structure.lattice
            energycompare = vrobj.final_energy

            os.mkdir('vac_1_As')
            os.mkdir('vac_1_As/charge_0')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml'), 'vac_1_As/charge_0/vasprun.xml')
            os.mkdir('vac_1_As/charge_0/LOCPOT') #locpot path just needs to exist
            transchg0 = {'charge': 0, 'supercell': [3, 3, 3], 'defect_type': 'vac_1_As',
                         'defect_supercell_site': vrobj.final_structure.sites[0]}
            dumpfn(transchg0, 'vac_1_As/charge_0/transformation.json', cls=MontyEncoder)

            os.mkdir('vac_1_As/charge_-1')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml.dfpt.unconverged'), #make this one unconverged...
                                    'vac_1_As/charge_-1/vasprun.xml')
            os.mkdir('vac_1_As/charge_-1/LOCPOT') #locpot path just needs to exist
            transchgm1 = {'charge': -1, 'supercell': [3, 3, 3], 'defect_type': 'vac_1_As',
                         'defect_supercell_site': vrobj.final_structure.sites[0]}
            dumpfn(transchgm1, 'vac_1_As/charge_-1/transformation.json', cls=MontyEncoder)

            os.mkdir('as_1_Cs_on_As')
            os.mkdir('as_1_Cs_on_As/charge_2')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml'), 'as_1_Cs_on_As/charge_2/vasprun.xml')
            os.mkdir('as_1_Cs_on_As/charge_2/LOCPOT') #locpot path just needs to exist
            transchg2 = {'charge': 0, 'supercell': [3, 3, 3], 'defect_type': 'as_1_Cs_on_As',
                         'defect_supercell_site': vrobj.final_structure.sites[1]}
            dumpfn(transchg2, 'as_1_Cs_on_As/charge_2/transformation.json', cls=MontyEncoder)

            #now test parse_defect_calculations
            pp = PostProcess('.')
            pdd = pp.parse_defect_calculations()
            self.assertEqual(pdd['bulk_entry'].energy, energycompare)
            self.assertEqual(len(pdd['bulk_entry'].structure), 25)
            self.assertEqual(pdd['bulk_entry'].data['locpot_path'], os.path.abspath('bulk/LOCPOT'))
            self.assertEqual(pdd['bulk_entry'].data['supercell_size'], [3, 3, 3])

            self.assertEqual(len(pdd['defects']), 2)
            self.assertEqual(pdd['defects'][0].entry.energy, energycompare)
            self.assertEqual(len(pdd['defects'][0].entry.structure), 25)
            self.assertEqual(pdd['defects'][0].entry.data['locpot_path'],
                             os.path.abspath('vac_1_As/charge_0/LOCPOT'))
            self.assertEqual(pdd['defects'][0].full_name, 'vac_1_As_0')
            self.assertFalse(pdd['defects'][0].multiplicity)
            self.assertEqual(list(pdd['defects'][0].site.coords),
                             list(vrobj.final_structure.sites[0].coords))
            self.assertEqual(list(pdd['defects'][1].site.coords),
                             list(vrobj.final_structure.sites[1].coords))
            self.assertEqual(pdd['defects'][0].supercell_size, [3, 3, 3])

            #now test compile_all quickly...
            ca = pp.compile_all()
            self.assertEqual(ca.keys(), ['epsilon', 'vbm', 'gap', 'defects', 'bulk_entry', 'mu_range'])
            answer = [[521.83587174, -0.00263523, 0.0026437],
                      [-0.00263523, 24.46276268, 5.381848290000001],
                      [0.0026437, 5.381848290000001, 24.42964103]]
            self.assertEqual(ca['epsilon'], answer)
            self.assertEqual(ca['vbm'], 1.5516000000000001)
            self.assertEqual(ca['gap'], 2.5390000000000001)
            self.assertEqual(len(ca['defects']), 2)
            self.assertEqual(ca['bulk_entry'].energy, energycompare)
            #INSERT a simpletest for mu_range...

    def test_get_vbm_bandgap(self):
        with ScratchDir('.'):
            os.mkdir('bulk')
            #first check the direct vasprun load method
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml'), 'bulk/vasprun.xml')
            pp = PostProcess('.')
            (testvbm, testgap) = pp.get_vbm_bandgap()
            self.assertEqual(testvbm, 1.5516000000000001)
            self.assertEqual(testgap, 2.5390000000000001)
            #secondly check a band gap pull
            pp = PostProcess('.', mpid='mp-2534') #GaAs mpid
            (testvbm_mpid, testgap_mpid) = pp.get_vbm_bandgap()
            self.assertEqual(testvbm_mpid, 2.6682)
            self.assertEqual(testgap_mpid, 0.18869999999999987)

    def test_get_chempot_limits(self):
        #Will wait to finish this until the ChemPotAnalyzer has been written...?
        #NOTE will want to small test to the compile_all test above when I DO implement something...
        pass

    def test_dielectric_calculation(self):
        with ScratchDir('.'):
            os.mkdir('dielectric')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml.dfpt.ionic'), 'dielectric/vasprun.xml')
            pp = PostProcess('.')
            eps = pp.parse_dielectric_calculation()
            answer = [[521.83587174, -0.00263523, 0.0026437],
                      [-0.00263523, 24.46276268, 5.381848290000001],
                      [0.0026437, 5.381848290000001, 24.42964103]]
            self.assertEqual(eps, answer)

    def test_compile_all(self):
        #this function calls all the other functions...to save time, using the first test to look at this test
        pass

if __name__ == '__main__':
    unittest.main()