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
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.point_defects import ValenceIonicRadiusEvaluator

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
    def __init__(self, structure, max_min_oxi={}, substitutions={}, 
                 oxi_states={}, cellmax=128, interstitial_sites=[], 
                 standardized=False):
        """
        Args:
            structure:
                the bulk structure
            max_min_oxi:
                The minimal and maximum oxidation state of each element as a 
                dict. For instance {"O":(-2,0)}. If not given, the oxi-states 
                of pymatgen are considered.
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
        as_defs = []
        sub_defs = []

        vac = Vacancy(struct, {}, {})
        vac_scs = vac.make_supercells_with_defects(sc_scale)
        struct_species = struct.types_of_specie
        nb_per_elts = {e:0 for e in structure.composition.elements}

        if not oxi_states:
            if len(struct_species) == 1:
                oxi_states = {struct.types_of_specie[0].symbol: 0}
            else:
                vir = ValenceIonicRadiusEvaluator(struct)
                oxi_states = vir.valences

        #if not oxi_states == 0 or len(max_min_oxi) == 0:
        #    if len(struct_species) == 1:
        #        struct_oxi = struct.copy()
        #        struct_oxi.add_oxidation_state_by_element(
        #            {struct.types_of_specie[0].symbol: 0})
        #    else:
        #        vba = BVAnalyzer()
        #        struct_oxi = vba.get_oxi_state_decorated_structure(struct)

        #if len(oxi_states) == 0:
        #    local_oxi_states = {}
        #    for s in struct_oxi:
        #        ele_sym = s.specie.element.symbol
        #        if ele_sym not in local_oxi_states.keys():
        #            local_oxi_states[ele_sym]=s.specie.oxi_state
        #else:
        #    local_oxi_states = oxi_states
        #if len(local_oxi_states) != len(struct_species):
        #    raise ValueError("Number of oxidation states does not"
        #                     " match number of specie types!")

        if not max_min_oxi:
            max_min_oxi = {}
            for s in struct_species:
                if isinstance(s, Specie):
                    el = s.element
                    max_min_oxi[el.symbol] = (el.min_oxidation_state, 
                            el.max_oxidation_state)
                elif isinstance(s, Element):
                    max_min_oxi[s.symbol] = (s.min_oxidation_state, 
                            s.max_oxidation_state)
                else:
                    continue
            #local_max_min_oxi = {}
            #for s in struct_oxi:
            #    spec = s.specie
            #    spec_oxi = spec.oxi_state
            #    ele = spec.element
            #    ele_sym = ele.symbol
            #    if ele_sym not in local_max_min_oxi.keys():
            #        oxi_min = ele.min_oxidation_state
            #        oxi_max = ele.max_oxidation_state
            #        n_oxi = len(ele.oxidation_states)
            #        if len(struct_species) > 1 and n_oxi > 1:
            #            oxis = list(sorted(ele.oxidation_states))
            #            if spec_oxi <= 0 and oxi_max > 0:
            #                for i in range(n_oxi):
            #                    oxi_max = oxis[n_oxi-1-i]
            #                    if oxi_max == spec_oxi:
            #                        break
            #                    elif oxi_max < spec_oxi:
            #                        raise ValueError("Unexpected oxidation"
            #                                         " state!")
            #            elif spec_oxi >= 0 and oxi_min < 0:
            #                for i in range(n_oxi):
            #                    oxi_min = oxis[i]
            #                    if oxi_min == spec_oxi:
            #                        break
            #                    elif oxi_min > spec_oxi:
            #                        raise ValueError("Unexpected oxidation"
            #                                         " state!")
            #        local_max_min_oxi[ele_sym]=(oxi_min, oxi_max)
        #else:
            #local_max_min_oxi = max_min_oxi
        #if len(local_max_min_oxi) != len(struct_species):
        #    raise ValueError("Number of ranges of oxidation states does"
        #                     " not match number of specie types!")

        for i in range(vac.defectsite_count()):
            vac_site = vac.get_defectsite(i)
            site_mult = vac.get_defectsite_multiplicity(i)
            site_mult = int(site_mult/conv_prim_rat)
            vac_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol
            vac_sc = vac_scs[i+1]
            vac_sc_site = list(set(vac_scs[0].sites) - set(vac_sc.sites))[0]

            list_charges=[]
            for c in range(max_min_oxi[vac_symbol][0]-1, 
                    max_min_oxi[vac_symbol][1]+2):
                list_charges.append(-c)
            nb_per_elts[vac_specie] += 1

            vacancies.append({
                'name': vac_symbol+str(nb_per_elts[vac_specie])+"_vac",
                'unique_site': vac_site,
                'supercell':{'size':sc_scale,'structure':vac_sc},
                'charges':list_charges })

            # Antisite defects generation
            for as_specie in set(struct_species)-set([vac_specie]):
                as_symbol = as_specie.symbol
                as_sc = vac_sc.copy()
                as_sc.append(as_symbol, vac_sc_site.frac_coords)
                oxi_min = min(max_min_oxi[as_symbol][0]-1,
                        max_min_oxi[vac_symbol][0]-1,0)
                oxi_max = max(max_min_oxi[as_symbol][1]+1,
                        max_min_oxi[vac_symbol][0]+1,1)
                as_defs.append({
                    'name': vac_symbol+str(nb_per_elts[vac_specie])+ \
                            "_subst_"+as_symbol,
                    'unique_site': vac_site,
                    'supercell':{'size':sc_scale,'structure':as_sc},
                    'charges':[c for c in range(oxi_min, oxi_max+1)]})

            # Substitutional defects generation
            if vac_symbol in substitutions:
                for subspecie_symbol in substitutions[vac_symbol]:
                    sub_sc = vac_sc.copy()
                    sub_sc.append(subspecie_symbol, vac_sc_site.frac_coords)
                    oxi_min = min(max_min_oxi[subspecie_symbol][0]-1,0)
                    oxi_max = max(max_min_oxi[subspecie_symbol][1]+1,1)
                    sub_defs.append({
                        'name': vac_symbol+str(nb_per_elts[vac_specie])+ \
                                "_subst_"+subspecie_symbol,
                        'unique_site': vac_site,
                        'supercell':{'size':sc_scale,'structure':sub_sc},
                        'charges':[c-oxi_states[vac_symbol] for c in range(
                            oxi_min, oxi_max+1)]})

        self.defects['vacancies'] = vacancies 
        self.defects['substitutions'] = sub_defs
        self.defects['substitutions'] += as_defs

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
                    'charges':[c for c in range(max_min_oxi[elt][0],
                        max_min_oxi[elt][1]+1)]})
                count = count+1
        self.defects['interstitials'] = interstitials

    
    def make_interstitial(self, target_site, sc_scale):
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        sc.append(target_site.specie, target_site.frac_coords)
        
        return sc
