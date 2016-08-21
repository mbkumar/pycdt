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

from monty.json import MontyDecoder
from monty.tempfile import ScratchDir
from monty.serialization import loadfn
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from pycdt.utils.vasp import make_vasp_defect_files, \
        make_vasp_dielectric_files

file_loc = os.path.join('..', '..', '..', 'test_files')

class VaspDefectFilesTest(unittest.TestCase):
    def setUp(self):
        self.defects = loadfn(os.path.join(file_loc, 'Cr2O3_defects.json'))
        self.user_settings = loadfn(os.path.join(file_loc,
                                                 'test_vasp_settings.yaml'))
        self.path = 'Cr2O3'
        self.neutral_Cr_vac_incar = {}

    def test_neutral_defect_incar(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            incar_loc = os.path.join(self.path, 'vac_1_Cr', 'charge_0')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertEqual(incar, self.neutral_Cr_vac_incar)

    def test_charged_defect_incar(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            incar_loc = os.path.join(self.path, 'vac_1_Cr', 'charge_-3')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertIsNotNone(incar.pop('NELECT', None))
            self.assertEqual(incar, self.neutral_Cr_vac_incar)

    def test_bulk_incar(self):
        self.cr2o3_bulk_incar = {}
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            incar_loc = os.path.join(self.path, 'bulk')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertEqual(incar, self.cr2O3_bulk_incar)

    def test_kpoints(self):
        self.cr2o3_kpoints = {}
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            kpoints_loc = os.path.join(self.path, 'bulk')
            kpoints = Kpoints.from_file(os.path.join(kpoints_loc, 'KPOINTS'))
            self.assertEqual(kpoints, self.cr2O3_kpoints)

    def test_poscar(self):
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            poscar_loc = os.path.join(self.path, 'bulk')
            poscar = Poscar.from_file(os.path.join(poscar_loc, 'POSCAR'))
            self.assertEqual(
                poscar,
                Poscar(self.defects['bulk']['supercell']['structure']))

    def test_user_settings_bulk(self):
        self.cr2o3_bulk_incar = {}
        with ScratchDir('.'):
            make_vasp_defect_files(self.defects, self.path)
            incar_loc = os.path.join(self.path, 'bulk')
            incar = Incar.from_file(os.path.join(incar_loc, "INCAR"))
            self.assertEqual(incar, self.cr2O3_bulk_incar)
        pass

    def test_hse_settings(self):
        pass



class VaspDielectricFilesTest(unittest.TestCase):
    def setUp(self):
        pass

if __name__ == '__main__':
    unittest.main()
