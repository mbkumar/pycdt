# coding: utf-8

from __future__ import division

__author__ = "Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Development"
__date__ = "Jan 26, 2016"

import unittest
import os
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Locpot
from pycdt.corrections.finite_size_charge_correction import *


class EnergyConversionTest(unittest.TestCase):
    """
    Test for  k-vector to energy [eV] conversion and viceversa via hbar*k^2/2m
    """
    def setUp(self):
        pass
    def test_ev_to_k():
        pass
    def test_k_to_ev():
        pass


class ReciprocalVectorGenerationTest(unittest.TestCase):
    """
    Test for  k-vector to energy [eV] conversion and viceversa via hbar*k^2/2m
    """
    def setUp(self):
        self.struct = Structure.from_file('POSCAR')
        self.latt = self.struct.lattice
    def test_dg_vals(self):
        reci_latt = self.latt.reciprocal_lattice
        #print (reci_latt.abc)
        vol = self.latt.volume
        [a1, a2, a3] = self.latt.get_cartesian_coords(1)
        #print (a1, a2, a3)
        b1 = 2*np.pi/vol*np.cross(a2,a3)
        b2 = 2*np.pi/vol*np.cross(a3,a1)
        b3 = 2*np.pi/vol*np.cross(a1,a2)
        print (b1, b2, b3)
        norm = np.linalg.norm
        #print norm(b1), norm(b2), norm(b3)
        [bp1, bp2, bp3] = reci_latt.get_cartesian_coords(1)
        print bp1, bp2, bp3


class LocPotMeshTest(unittest.TestCase):
    """
    Test for  k-vector to energy [eV] conversion and viceversa via hbar*k^2/2m
    """
    def setUp(self):
        self.locpot = Locpot.from_file('LOCPOT')
        self.latt = self.locpot.structure.lattice
    def test_dim_vals(self):
        dim = self.locpot.dim
        print dim


class QModelTest(unittest.TestCase):
    """
    Test for the charge model used in correction 
    """
    def setUp(self):
        pass
    def test_rho_rec():
        pass
    def test_rho_rec_limit0():
        pass


class WignerSeitzRadiusTest(unittest.TestCase):
    """
    """
    def setUp(self):
        pass
    def test_radius():
        pass


class ChargeCorrectionTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_freysoldt_pc(self):
        pss

    def test_freysoldt_potalign(self):
        pass

    def test_kumagai_correction(self):
        pass

    def test_kumagai_pc(self):
        pass
    
    def test_kumagai_potalign(self):
        pass

    def test_markovpayne(self):
        pass

