from __future__ import division
"""
Computes the correction for surace energy error associated with vacancy 
formation energy by deducing the surface area of vacancy from the 
first principles calculation of vacancies with lda, pbe and pw91 
functionals.

The unit area correction computed is defined in Phys.Rev.B 73, 195123, 2006.
The surface area computation is defined in Phys.Rev.B 85, 144118, 2012.
"""

__author__ = "Bharat Medasani"

import math

from numpy import array, linalg, ones
from pymatgen.symmetry.analyzer import SymmetryAnalyzer
from pymatgen.io.vaspio_set import MPGGAVaspInputSet
from pymargen.io.vaspio.vasp_input import Potcar

ergpercmsq_to_evperangsq = 6.24150934e-5

se_corr = { 
        "LDA":{"A":448.454, "B":-55.845},
        "PW91":{"A":1577.2, "B":-231.29}, 
        "PBE":{"A":1193.7, "B":-174.37}}
Bohr_rad = 5.2917721092e-1

def unit_xc_correction(atoms, valence, volume, functional):
    """
    Computes the unit area surface energy correction for a functional
    Args:
        atoms: # of atoms in lattice
        valence: Valence of the atoms
        volume: Volume of the lattice.
        functional: Accepted values are PBE, PW91, and LDA
    Returns:
        Correction in ev/A^2
    """
    blk_electron_den = valence*atoms/volume
    rs = (3/(4*math.pi*blk_electron_den))**(1/3.0)
    rs = rs/Bohr_rad
    rspa = rs**(-5.0/2)
    rspb = rs**(-3.0/2)
    a = se_corr[functional]['A']
    b = se_corr[functional]['B']
    corr = (a*rspa + b*rspb) * ergpercmsq_to_evperangsq
    return corr

def correction(energy_dict, lattice_dict, valence_dict=None):
    """
    Computes the surface energy correction for a neutral vacancy in metal
    Args:
        energy_dict: Uncorrected vacancy formation energy for each functional 
        lattice_dict: Bulk structures for each functional
        valence_dict: # of valence electrons for each functional
            If not specified obtained from the potcar files. 
            (Needs vasp potcar files if not specified.)
    Returns:
        Corrected vacancy formation energy and surface area of vacancies
    """

    lda_uc_e0 = energy_dict['LDA']
    lda_struct = lattice_dict['LDA']
    lda_uc_struct = SymmetryAnalyzer(
            lda_struct).get_conventional_standard_structure()
    lda_volume = lda_uc_struct.volume/lda_uc_struct.num_sites
    lda_surfarea = lda_volume**(2.0/3)

    pbe_uc_e0 = energy_dict['PBE']
    pbe_struct = lattice_dict['PBE']
    pbe_uc_struct = SymmetryAnalyzer(
            pbe_struct).get_conventional_standard_structure()
    pbe_volume = pbe_uc_struct.volume/pbe_uc_struct.num_sites
    pbe_surfarea = pbe_volume**(2.0/3)

    pw91_uc_e0 = energy_dict['PW91']
    pw91_struct = lattice_dict['PW91']
    pw91_uc_struct = SymmetryAnalyzer(
            pw91_struct).get_conventional_standard_structure()
    pw91_volume = pw91_uc_struct.volume/pw91_uc_struct.num_sites
    pw91_surfarea = pw91_volume**(2.0/3)

    if valence_dict:
        lda_valence = valence_dict['LDA']
        pbe_valence = valence_dict['PBE']
        pw91_valence = valence_dict['PW91']
    else:
        try:
            mpvis = MPGGAVaspINputSet()
        except:
            raise ValueError('POTCAR not found. Supply valence of element')
        potcar = mpvis.get_potcar(pbe_uc_struct)
        pbe_valence = potcar[0].zval
        potcar_dict = potcar.to_dict
        potcar_dict.update({'functional':'LDA'})
        potcar = Potcar.from_dict(potcar_dict)
        lda_valence = potcar[0].zval
        potcar_dict.update({'functional':'PW91'})
        potcar = Potcar.from_dict(potcar_dict)
        pw91_valence = potcar[0].zval

    lda_xc_cor = correction(1, lda_valence, lda_volume, 'LDA')
    pbe_xc_cor = correction(1, pbe_valence, pbe_volume, 'PBE')
    pw91_xc_cor = correction(1, pw91_valence, pw91_volume, 'PW91')


    x = [lda_xc_cor*lda_surfarea, 
         pbe_xc_cor*pbe_surfarea, 
         pw91_xc_cor*pw91_surfarea]
    A = array([x, ones(3)])
    y = [lda_uc_e0, pbe_uc_e0, pw91_uc_e0]

    w = linalg.lstsq(A.T, y)[0]

    surf_err = {
            'vac_surface_area':{
                'LDA':w[0]*lda_surfarea,
                'PBE':w[0]*pbe_surfarea,
                'PW91':w[0]*pw91_surfarea
                },
            'vacancy_form_energy_corrected':{
                'LDA':lda_uc_e0-w[0]*lda_surfarea*lda_xc_cor,
                'PBE':pbe_uc_e0-w[0]*pbe_surfarea*pbe_xc_cor,
                'PW91':pw91_uc_e0-w[0]*pw91_surfarea*pw91_xc_cor
                }
            }
    return surf_err
                '

    

