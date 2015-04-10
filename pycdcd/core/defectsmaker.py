#!/usr/bin/env python
from __future__ import division
"""
Code to generate charged defects structure.
Ideas from pydii's code and geoffroy's code are merged.
"""

__author__ = "Bharat Medasani, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com,geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

import copy

from pymatgen.core.structure import PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.point_defects import Vacancy

def get_sc_scale(inp_struct, final_site_no):
    lengths = inp_struct.lattice.abc
    no_sites = inp_struct.num_sites
    mult = (final_site_no/no_sites*lengths[0]*lengths[1]*lengths[2]) ** (1/3)
    num_mult = [int(round(mult/l)) for l in lengths]
    num_mult = [i if i > 0 else 1 for i in num_mult]
    sc = copy.deepcopy(inp_struct)
    sc.make_supercell(num_mult)
    if sc.num_sites > final_site_no:
        max_sc_dim = max(num_mult)
        i = num_mult.index(max_sc_dim)
        num_mult[i] -= 1
    return num_mult

def get_optimized_sc_scale(inp_struct, final_site_no):
    target_site = inp_struct.sites[0]
    dictio={}
    result=[]
    for k1 in range(1,4):
        for k2 in range(1,4):
            for k3 in range(1,4):
                struct = inp_struct.copy()
                struct.make_supercell([k1,k2,k3])
                if len(struct.sites) > final_site_no:
                    continue
                site_target=None
                index=None
                for i in range(struct.num_sites):
                    s=struct._sites[i]
                    if s.distance_from_point(target_site.coords)<0.001:
                        index=i
                min=1000.0
                for a in range(-1,1):
                    for b in range(-1,1):
                        for c in range(-1,1):
                            distance = struct.get_distance(index,index,(a,b,c))
                            if  distance < min and distance>0.00001:
                                min = distance
                min=round(min,3)
                if dictio.has_key(min):
                    if dictio[min]['num_sites'] > struct.num_sites:
                        dictio[min]['num_sites'] = struct.num_sites
                        dictio[min]['supercell'] = [k1,k2,k3]
                else:
                    dictio[min]={}
                    dictio[min]['num_sites'] = struct.num_sites
                    dictio[min]['supercell'] = [k1,k2,k3]
    min=-1.0
    biggest=None
    for c in dictio:
        if c>min:
            biggest=dictio[c]['supercell']
            min=c
    return biggest


class ChargedDefectsStructures(object):
    """
    A class to generate charged defective structures for use in first 
    principles supercell formalism. The standard defects such as antisites, 
    vacancies are generated.
    TODO: develop a better way to find interstitials
    """
    def __init__(self, structure, max_min_oxi=None, substitutions=None, 
                 oxi_states=None, cellmax=128, interstitial_sites=[], 
                 standardized=False):
        """
        Args:
            structure:
                the bulk structure
            max_min_oxi:
                The minimal and maximum oxidation state of each element as a 
                dict. For instance {"O":(-2,0)}
            substitutions:
                The allowed substitutions of elements as a dict. If not given, 
                intrinsic defects are computed. If given, intrinsic (e.g., 
                anti-sites) and extrinsic are considered explicitly specified. 
                Example: {"Co":["Zn","Mn"]} means Co sites can be substituted 
                by Mn or Zn.
            oxi_states:
                The oxidation state of the elements in the compound e.g. 
                {"Fe":2,"O":-2}. If not given, the oxidation state of each
                site is computed with bond valence sum. WARNING: Bond-valence 
                method can fail for mixed-valence compounds
            cellmax:
                Maximum number of atoms allowed in the supercell
            interstitials_sites:
                A list of PeriodicSites in the bulk structure on which we put 
                an interstitial
        """

        self.defects = []
        self.cellmax = cellmax
        self.struct = structure

        spa = SpacegroupAnalyzer(structure,symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        if standardized:
            struct = prim_struct
        else:
            struct = structure
        conv_prim_rat = int(struct.num_sites/prim_struct.num_sites)
        sc_scale = get_optimized_sc_scale(struct,cellmax)
        self.defects = {}
        sc = struct.copy()
        sc.make_supercell(sc_scale)
        self.defects['bulk'] = {'name':'bulk',
                'supercell':{'size':sc_scale,'structure':sc}}

        vacancies = []
        sub_defs = []

        vac = Vacancy(struct, {}, {})
        vac_scs = vac.make_supercells_with_defects(sc_scale)
        struct_species = struct.types_of_specie
        nb_per_elts = {e:0 for e in structure.composition.elements}

        for i in range(vac.defectsite_count()):
            vac_site = vac.get_defectsite(i)
            site_mult = vac.get_defectsite_multiplicity(i)
            site_mult = int(site_mult/conv_prim_rat)
            vac_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol
            vac_sc = vac_scs[i+1]

            list_charges=[]
            for c in range(max_min_oxi[vac_symbol][0], 
                    max_min_oxi[vac_symbol][1]+1):
                list_charges.append(-c)
            nb_per_elts[vac_specie] += 1

            vacancies.append({
                'name': vac_symbol+str(nb_per_elts[vac_specie])+"_vac",
                'unique_site': vac_site,
                'supercell':{'size':sc_scale,'structure':vac_sc},
                'charges':list_charges })

            # Substitutional defects generation
            if vac_symbol in substitutions:
                for subspecie_symbol in substitutions[vac_symbol]:
                    sub_sc = vac_sc.copy()
                    sub_sc.append(subspecie_symbol, vac_site.frac_coords)
                    sub_defs.append({
                        'name': vac_symbol+str(nb_per_elts[vac_specie])+ \
                                "_subst_"+subspecie_symbol,
                        'unique_site': vac_site,
                        'supercell':{'size':sc_scale,'structure':sub_sc},
                        'charges':[c-oxi_states[vac_symbol] for c in range(
                            max_min_oxi[subspecie_symbol][0],
                            max_min_oxi[subspecie_symbol][1]+1)]})

        self.defects['vacancies'] = vacancies 
        self.defects['substitutions'] = sub_defs

        #interstitials
        interstitials = []
        for elt in self.struct.composition.elements:
            count = 1
            for frac_coord in interstitial_sites:
                site = PeriodicSite(elt, frac_coord, structure.lattice)
                interstitials.append({
                    'name':elt.symbol+str(count)+"_inter",
                    'unique_site':site,
                    'supercell':{'size':s_size,
                        'structure':self.make_interstitial(site, sc_scale)},
                    'charges':[c for c in range(max_min_oxid[elt][0],
                        max_min_oxid[elt][1]+1)]})
                count = count+1
        self.defects['interstitials'] = interstitials

    
    def make_interstitial(self, target_site, sc_scale):
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        sc.append(target_site.specie, target_site.frac_coords)
        
        return sc
