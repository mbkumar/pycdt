# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "May 6, 2015"

import os
import unittest
from shutil import copyfile

from monty.serialization import dumpfn
from monty.json import MontyEncoder
from monty.tempfile import ScratchDir

from pymatgen import __file__ as initfilep
from pymatgen.core import Element
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest

from pycdt.utils.parse_calculations import PostProcess

pmgtestfiles_loc = os.path.join(
        os.path.split(os.path.split(initfilep)[0])[0], 'test_files')
file_loc = os.path.abspath(
        os.path.join(__file__, '..', '..', '..', '..', 'test_files')) #Pycdt Testfiles


class LegacyConversionTest(PymatgenTest):
    def test_convert_cd_to_de(self):
        #TODO: make this unit test
        pass

class SingleDefectParserTest(PymatgenTest):
    def test_from_paths(self):
        #TODO: make this unit test
        pass

    def test_freysoldt_loader(self):
        #TODO: make this unit test
        pass

    def test_kumagai_loader(self):
        #TODO: make this unit test
        pass

    def test_get_stdrd_metadata(self):
        #TODO: make this unit test
        pass

    def test_identify_defect_info(self):
        #TODO: make this unit test
        pass

    def test_get_bulk_gap_data(self):
        #TODO: make this unit test
        pass


class PostProcessTest(PymatgenTest):
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

            vrobj = Vasprun( os.path.join(pmgtestfiles_loc, 'vasprun.xml'))

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

            os.mkdir('sub_1_Cs_on_As')
            os.mkdir('sub_1_Cs_on_As/charge_2')
            copyfile(os.path.join(pmgtestfiles_loc, 'vasprun.xml'), 'sub_1_Cs_on_As/charge_2/vasprun.xml')
            os.mkdir('sub_1_Cs_on_As/charge_2/LOCPOT') #locpot path just needs to exist
            transchg2 = {'charge': 0, 'supercell': [3, 3, 3], 'defect_type': 'sub_1_Cs_on_As',
                         'defect_supercell_site': vrobj.final_structure.sites[1], 'substitution_specie': 'Cs'}
            dumpfn(transchg2, 'sub_1_Cs_on_As/charge_2/transformation.json', cls=MontyEncoder)

            #now test parse_defect_calculations
            pp = PostProcess('.')
            pdd = pp.parse_defect_calculations()
            self.assertEqual(pdd['bulk_entry'].energy, vrobj.final_energy)
            self.assertEqual(len(pdd['bulk_entry'].structure), 25)
            self.assertEqual(pdd['bulk_entry'].data['bulk_path'], './bulk')
            self.assertEqual(pdd['bulk_entry'].data['supercell_size'], [3, 3, 3])

            self.assertEqual(len(pdd['defects']), 2)
            self.assertEqual(pdd['defects'][0].energy, 0.)
            self.assertEqual(len(pdd['defects'][0].bulk_structure), 25)
            self.assertEqual(pdd['defects'][0].parameters['defect_path'],
                             './vac_1_As/charge_0')
            self.assertEqual(pdd['defects'][0].parameters['fldr_name']+'_'+str(pdd['defects'][0].charge),
                             'vac_1_As_0')
            self.assertEqual(pdd['defects'][0].multiplicity, 1)
            self.assertEqual(list(pdd['defects'][0].defect.site.coords),
                             list(vrobj.final_structure.sites[0].coords))
            self.assertEqual(pdd['defects'][0].parameters['supercell_size'], [3, 3, 3])

            self.assertEqual(list(pdd['defects'][1].defect.site.coords),
                             list(vrobj.final_structure.sites[1].coords))
            self.assertEqual( pdd['defects'][1].defect.site.specie.symbol,
                              'Cs')

            #now test compile_all quickly...
            ca = pp.compile_all()
            lk = sorted(list(ca.keys()))
            self.assertEqual(len(lk), 6)
            self.assertEqual(lk, sorted(['epsilon', 'vbm', 'gap', 'defects', 'bulk_entry', 'mu_range']))
            answer = [[521.83587174, -0.00263523, 0.0026437],
                      [-0.00263523, 24.46276268, 5.381848290000001],
                      [0.0026437, 5.381848290000001, 24.42964103]]
            self.assertEqual(ca['epsilon'], answer)
            self.assertEqual(ca['vbm'], 1.5516000000000001)
            self.assertEqual(ca['gap'], 2.5390000000000001)
            self.assertEqual(len(ca['defects']), 2)
            self.assertEqual(ca['bulk_entry'].energy, vrobj.final_energy)
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
        with ScratchDir('.'):
            #make a fake file structure to parse vaspruns
            os.mkdir('bulk')
            copyfile(os.path.join(file_loc, 'vasprun.xml_GaAs'), 'bulk/vasprun.xml')
            pp = PostProcess('.')
            gaas_cp = pp.get_chempot_limits()
            self.assertEqual(set([u'As-GaAs', u'Ga-GaAs']), set(gaas_cp.keys()))
            self.assertEqual([-4.6580705550000001,  -4.9884807425],
                             [gaas_cp['As-GaAs'][Element(elt)] for elt in ['As','Ga']])
            self.assertEqual([-6.609317857500001,  -3.03723344],
                             [gaas_cp['Ga-GaAs'][Element(elt)] for elt in ['As','Ga']])

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

if __name__ == '__main__':
    unittest.main()
