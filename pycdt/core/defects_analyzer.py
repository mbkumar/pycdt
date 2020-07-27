#!/usr/bin/env python

__author__ = "Geoffroy Hautier, Bharat Medasani,  Danny Broberg"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier, Bharat Medasani"
__email__ = "geoffroy@uclouvain.be, mbkumar@gmail.com"
__status__ = "Development"
__date__ = "November 4, 2012"

from math import sqrt, pi, exp
from collections import defaultdict 
from itertools import combinations

import os
import numpy as np

from pymatgen.core import Element
from pymatgen.core.structure import PeriodicSite, Structure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pycdt.corrections.finite_size_charge_correction import get_correction_freysoldt, get_correction_kumagai
from pycdt.utils.parse_calculations import SingleDefectParser
from pycdt.utils.units import kb, conv, hbar

import warnings
warnings.simplefilter('default')


def freysoldt_correction_from_paths( defect_file_path, bulk_file_path, dielectric,
                                     defect_charge, plot=False):
    """
    A function for performing the Freysoldt correction with a set of file paths.
    If this correction is used, please reference Freysoldt's original paper.
    doi: 10.1103/PhysRevLett.102.016402

    Does not require transformation.json file to exist in file path.

    :param defect_file_path (str): file path to defect folder of interest
    :param bulk_file_path (str): file path to bulk folder of interest
    :param dielectric (float or 3x3 matrix): Dielectric constant (or tensor) for the structure
    :param defect_charge (int): charge of defect structure of interest
    :param plot (bool): allow for plotting electrostatic potential
    :return:
        Dictionary of Freysoldt Correction for defect
    """
    sdp = SingleDefectParser.from_paths( defect_file_path, bulk_file_path, dielectric, defect_charge)
    _ = sdp.freysoldt_loader()
    plt_title = os.path.join( defect_file_path,
                             "{}_chg_{}".format(sdp.defect_entry.name, defect_charge)) if plot else None
    correction = get_correction_freysoldt(sdp.defect_entry,
                                          dielectric,
                                          title=plt_title)

    return correction

def kumagai_correction_from_paths( defect_file_path, bulk_file_path, dielectric,
                                   defect_charge, plot=False):
    """
    A function for performing the Kumagai correction with a set of file paths.
    If this correction is used, please reference Kumagai and Oba's original paper
    (doi: 10.1103/PhysRevB.89.195205) as well as Freysoldt's original
    paper (doi: 10.1103/PhysRevLett.102.016402

    Does not require transformation.json file to exist in file path.

    :param defect_file_path (str): file path to defect folder of interest
    :param bulk_file_path (str): file path to bulk folder of interest
    :param dielectric (float or 3x3 matrix): Dielectric constant (or tensor) for the structure
    :param defect_charge (int): charge of defect structure of interest
    :param plot (bool): allow for plotting electrostatic potential
    :return:
        Dictionary of Kumagai Correction for defect
    """
    sdp = SingleDefectParser.from_paths( defect_file_path, bulk_file_path, dielectric, defect_charge)
    _ = sdp.kumagai_loader()
    plt_title = os.path.join( defect_file_path,
                             "{}_chg_{}".format(sdp.defect_entry.name, defect_charge)) if plot else None
    correction = get_correction_kumagai(sdp.defect_entry,
                                          dielectric,
                                          title=plt_title)

    return correction



