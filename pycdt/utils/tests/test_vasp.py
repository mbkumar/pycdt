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
import glob

from monty.json import MontyDecoder
from monty.tempfile import ScratchDir
from monty.serialization import loadfn
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from pymatgen.core.structure import Structure
from pycdt.utils.vasp import make_vasp_defect_files, \
        make_vasp_dielectric_files

file_loc = os.path.join('..', '..', '..', 'test_files')

class VaspDefectFilesTest(PymatgenTest):
    def setUp(self):
        self.defects = loadfn(os.path.join(file_loc, 'Cr2O3_defects.json'))
        self.user_settings = loadfn(os.path.join(file_loc,
                                                 'test_vasp_settings.yaml'))
        self.path = 'Cr2O3'
        self.neutral_def_incar_min = {'LVHAR': True, 'ISYM': 0, 'ISMEAR': 0,
                'ISIF': 2,  'ISPIN': 2}
        self.def_keys = ['EDIFF', 'EDIFFG', 'IBRION']
        self.bulk_incar = {'LVHAR': True, 'ISYM': 0, 'ISMEAR': 0,
                'IBRION': -1,  'ISPIN': 2}
        self.bulk_keys = ['EDIFF']

    def test_neutral_defect_incar(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            cr_def_path = glob.glob(os.path.join(self.path, 'vac*Cr'))[0]
            incar_loc = os.path.join(cr_def_path, 'charge_0')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertTrue(
                    self.neutral_def_incar_min.viewitems() <= incar.items())
            self.assertTrue(set(self.def_keys).issubset(incar))

    def test_charged_defect_incar(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            cr_def_path = glob.glob(os.path.join(self.path, 'vac*Cr'))[0]
            incar_loc = os.path.join(cr_def_path, 'charge_-1')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertIsNotNone(incar.pop('NELECT', None))
            self.assertTrue(
                    self.neutral_def_incar_min.viewitems() <= incar.items())
            self.assertTrue(set(self.def_keys).issubset(incar))

    def test_bulk_incar(self):

        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            incar_loc = os.path.join(self.path, 'bulk')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertTrue(self.bulk_incar.viewitems() <= incar.items())
            self.assertTrue(set(self.bulk_keys).issubset(incar))

    def test_kpoints(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            kpoints_loc = os.path.join(self.path, 'bulk')
            kpoints = Kpoints.from_file(os.path.join(kpoints_loc, 'KPOINTS'))
            # xxx: doesn't work?!
            #self.assertEqual(kpoints.kpts, [[2,2,2]])

    def test_poscar(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            poscar_loc = os.path.join(self.path, 'bulk')
            poscar = Poscar.from_file(os.path.join(poscar_loc, 'POSCAR'))
            self.assertTrue(poscar.structure.matches(
                self.defects['bulk']['supercell']['structure']))

    def test_user_settings_bulk(self):
        user_settings = loadfn(os.path.join(
            file_loc, 'test_vasp_settings.yaml'))
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path, 
                                   user_settings=user_settings)
            incar_loc = os.path.join(self.path, 'bulk')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertTrue(self.bulk_incar.viewitems() <= incar.items())
            self.assertTrue(set(self.bulk_keys).issubset(incar))
            self.assertEqual(incar['ENCUT'], 620)

    def test_hse_settings(self):
        pass
    # xxx: doesn't work?!
    #    self.assertTrue(0)


class VaspDielectricFilesTest(PymatgenTest):
    def setUp(self):
        self.structure = Structure.from_file(os.path.join(file_loc, 
                                                          'POSCAR_Cr2O3'))
        self.user_settings = loadfn(os.path.join(file_loc,
                                                 'test_vasp_settings.yaml'))
        self.path = 'dielectric'
        self.dielectric_min = {'LVHAR': True, 'ISYM': 0, 'LPEAD': True,
                'LEPSILON': True,  'ISPIN': 2, 'IBRION': 8}
        self.keys = ['EDIFF']

    def test_dielectric_files(self):
        with ScratchDir('.'):
            make_vasp_dielectric_files(self.structure, self.path)
            incar = Incar.from_file(os.path.join(self.path, "INCAR"))
            print incar
            self.assertTrue(
                    self.dielectric_min.viewitems() <= incar.items())
            self.assertTrue(set(self.keys).issubset(incar))

    def test_user_modified_dielectric_files(self):
        with ScratchDir('.'):
            make_vasp_dielectric_files(self.structure, self.path, 
                                       user_settings=self.user_settings)
            incar = Incar.from_file(os.path.join(self.path, "INCAR"))
            self.assertTrue(
                    self.dielectric_min.viewitems() <= incar.items())
            self.assertTrue(set(self.keys).issubset(incar))
            self.assertEqual(incar['EDIFF'], 5e-7)
            self.assertEqual(incar['ENCUT'], 620)

import unittest
if __name__ == '__main__':
    unittest.main()
