"""
This module implements helper functions used in both Freysoldt and Kumagai
correction metodules 
"""
__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import sys
import math
import logging

import numpy as np
norm = np.linalg.norm

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


def cleanlat(dat):
    """
    # TODO: Identify better function and input variable name
    return lattice constants
    Args:
        dat: array of lattice vectors

    Returns:
        (double) Lattice constants (in same units as lattice vectors)
    """
    return list(map(norm, dat))


def genrecip(a1, a2, a3, encut):
    """
    Args:
        a1, a2, a3: lattice vectors in bohr
        encut: energy cut off in eV
    Returns:
        reciprocal lattice vectors with energy less than encut
    """
    vol = np.dot(a1, np.cross(a2, a3))  # 1/bohr^3
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  # units 1/bohr
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)

    # create list of recip space vectors that satisfy |i*b1+j*b2+k*b3|<=encut
    gcut = eV_to_k(encut)
    imax = int(math.ceil(gcut/min(map(norm, [b1, b2, b3]))))

    #recip = []
    for i in range(-imax, imax + 1):
        for j in range(-imax, imax + 1):
            for k in range(-imax, imax + 1):
                vec = i*b1 + j*b2 + k*b3
                en = invang_to_ev * (((1.0/ang_to_bohr) * norm(vec))**2)
                if (en <= encut and en != 0):
                    yield vec
                    #recip.append(vec)

    #return recip  #output is 1/bohr recip


def generate_reciprocal_vectors_squared(a1, a2, a3, encut):
    """
    Generate reciprocal vector magnitudes within the cutoff along the specied
    lattice vectors. 
    Args:
        a1: Lattice vector a (in Bohrs)
        a2: Lattice vector b (in Bohrs)
        a3: Lattice vector c (in Bohrs)
        encut: Reciprocal vector energy cutoff

    Returns:
        [[g1^2], [g2^2], ...] Square of reciprocal vectors (1/Bohr)^2 
        determined by a1, a2, a3 and whose magntidue is less than gcut^2.
    """
    vol = np.dot(a1, np.cross(a2, a3))  
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)

    # Max (i,j,k) that doesn't upset the condition |i*b1+j*b2+k*b3|<=gcut
    gcut = eV_to_k(encut)
    imax = int(math.ceil(gcut/min(norm(b1), norm(b2), norm(b3))))
    gcut2 = gcut * gcut
    #recip = []
    for i in range(-imax, imax+1):
        for j in range(-imax, imax+1):
            for k in range(-imax, imax+1):
                vec = i*b1 + j*b2 + k*b3
                vec2 = np.dot(vec,vec)
                if (vec2 <= gcut2 and vec2 != 0.0):
                    #recip.append(vec2)
                    yield vec2
    #return recip


def closestsites(struct_blk, struct_def, pos):
    """
    Returns closest site to the input position
    for both bulk and defect structures
    Args:
        struct_blk: Bulk structure
        struct_def: Defect structure
        pos: Position 
    Return: (site object, dist, index)
    """
    # input bulk and defect structures and get site that is nearest 
    # to the (cartesian) input position
    blk_close_sites = struct_blk.get_sites_in_sphere(pos, 5, include_index=True)
    blk_close_sites.sort(key=lambda x:x[1])
    def_close_sites = struct_def.get_sites_in_sphere(pos, 5, include_index=True)
    def_close_sites.sort(key=lambda x:x[1])

    return blk_close_sites[0], def_close_sites[0] 


def find_defect_pos(struct_blk, struct_def):
    """
    output cartesian coords of defect in bulk,defect cells.

    If vacancy defectpos=None, 
    if interstitial bulkpos=None, 
    if antisite/sub then both defined
    """
    if len(struct_blk.sites) > len(struct_def.sites):
        type_def = 'vacancy'
    elif len(struct_blk.sites) < len(struct_def.sites):
        type_def = 'interstitial'
    else:
        type_def = 'substitution' # also corresponds to antisite

    sitematching = []
    for site in struct_blk.sites:
        blksite, defsite = closestsites(struct_blk, struct_def, site.coords)
        if blksite[0].specie.symbol != defsite[0].specie.symbol:
            if type_def == 'vacancy':
                return blksite[0].coords, None
            elif type_def == 'interstitial':
                return None, defsite[0].coords
            else: #subs or antisite type
                return blksite[0].coords, defsite[0].coords
        sitematching.append([blksite[0], blksite[1], defsite[0], defsite[1]])

    if type_def == 'vacancy': 
        #in case site type is same for closest site to vacancy
        sitematching.sort(key=lambda x:x[3])
        vacant = sitematching[-1]
        return vacant[0].coords, None
    elif type_def == 'interstitial':
        #just in case site type is same for closest site to interstit
        sitematching.sort(key=lambda x:x[1])
        interstit = sitematching[-1]
        return  None, interstit[2].coords

    return None, None #if you get here there is an error