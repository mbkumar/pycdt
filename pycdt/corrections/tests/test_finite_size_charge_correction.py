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

from pycdt.corrections.finite_size_charge_correction import *


def EnergyConversionTest(unittest.TestCase):
    """
    Test for  k-vector to energy [eV] conversion and viceversa via hbar*k^2/2m
    """
    def setUp(self):
        pass
    def test_ev_to_k():
        pass
    def test_k_to_ev():
        pass


def ReciprocalVectorGenerationTest(unittest.TestCase):
    """
    Test for  k-vector to energy [eV] conversion and viceversa via hbar*k^2/2m
    """
    def setUp(self):
        pass
    def test_gen():
        pass


def QModelTest(unittest.TestCase):
    """
    Test for the charge model used in correction 
    """
    def setUp(self):
        pass
    def test_rho_rec():
        pass
    def test_rho_rec_limit0():
        pass


def WignerSeitzRadiusTest(unittest.TestCase):
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

