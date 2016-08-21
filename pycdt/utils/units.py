"""
This module contains unit conversion constants and  functions
"""
__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import math
import numpy as np

# Define conversion_constants
hart_to_ev = 27.2114
ang_to_bohr = 1.8897
invang_to_ev = 3.80986
kb = 8.6173324e-5 #eV / K
hbar = 6.58211928e-16  #eV s

#no idea what is meaning for this constant. Used in defects_analyzer
conv = math.sqrt((9.1*1e-31)**3)*math.sqrt((1.6*1e-19)**3)/((1.05*1e-34)**3)

def k_to_eV(g):
    """
    Convert a k-vector to energy [eV] via hbar*k^2/2m
    Args:
        a: Reciprocal vector (units of 1/A).

    Returns:
        (double) Energy in eV
    """
    return invang_to_ev * np.dot(g,g)


def eV_to_k(energy):
    """
    Convert energy to reciprocal vector magnitude k via hbar*k^2/2m
    Args:
        a: Energy in eV.

    Returns:
        (double) Reciprocal vector magnitude (units of 1/Bohr).
    """
    return math.sqrt(energy/invang_to_ev) * ang_to_bohr