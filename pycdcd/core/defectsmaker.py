#!/usr/bin/env python
from __future__ import divison

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
#from pymatgen.transformations.defect_transformation import \
#        VacancyTransformation,AntisiteDefectTransformation, \
#        SubstitutionDefectTransformation

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


class ChargedDefectsStructures(object):
    """
    A class to generate charged defective structures for use in first 
    principles supercell formalism. The standard defects such as antisites, 
    vacancies are generated.
    TODO: develop a better way to find interstitials
    """
    def __init__(self, structure, max_min_oxid=None, intrinsic=True,
                 allowed_subst=None, oxid_states=None, cellmax=128, 
                 interstitial_sites=[], standardized=False):
        """
        Args:
            structure:
                the bulk structure
            max_min_oxid:
                The minimal and maximum oxidation state of each element as a 
                dict. For instance {"O":(-2,0)}
            intrinsic_flag:
                If True, compute intrinsic defects, (vacancies and antisites)
            allowed_subst:
                The allowed substitutions of elements as a dict. If not given, 
                intrinsic defects are computed. If given, intrinsic (e.g., 
                anti-sites) and extrinsic are considered explicitly specified. 
                Example: {"Co":["Zn","Mn"]} means Co sites can be substituted 
                by Mn or Zn.
            oxid_states:
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

        spa = SpacegroupAnalyzer(structure,symprec=1e-2)
        self.struct = spa.get_symmetrized_structure()

        prim_struct = spa.get_primitive_standard_structure()
        if standardized:
            struct = prim_struct
        else:
            struct = structure
        conv_prim_rat = int(struct.num_sites/prim_struct.num_sites)
        sc_scale = get_sc_scale(struct,cellmax)
        self.defects = {}
        sc = struct.copy()
        sc.make_supercell(sc_scale)
        self.defects['bulk'] = {'name':'bulk',
                'supercell':{'size':sc_scale,'structure':sc}}

        if intrinsic:
            vacancies = []
            antisites = []

        vac = Vacancy(struct, {}, {})
        vac_scs = vac.make_supercells_with_defects(sc_scale)
        struct_species = struct.types_of_specie
        nb_per_elts = {e:0 for e in structure.composition.elements}
        for i in range(vac.defectsite_count):
            vac_site = vac.get_defectsite(i)
            site_mult = vac.get_defectsite_multiplicity(i)
            site_mult = int(site_mult/conv_prim_rat)
            vac_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol
            vac_sc = vac_scs[i+1]

            list_charges=[]
            for c in range(max_min_oxid[vac_symbol][0], 
                    max_min_oxid[vac_symbol][1]+1):
                list_charges.append(-c)
            nb_per_elts[vac_specie] += 1

            if intrinsci:
                vacancies.append({
                    'name': vac_symbol+str(nb_per_elts[vac_specie])+"_vac",
                    'unique_site': vac_site,
                    'supercell':{'size':sc_scale,'structure':vac_sc},
                    'charges':list_charges })

                # Antisite generation at all vacancy sites
                for specie in set(struct_species)-set([vac_specie]):
                    subspecie_symbol = specie.symbol
                    antisite_sc = vac_sc.copy()
                    antisite_sc.append(specie, vac_site.frac_coords)
                    antisites.append({
                        'name': vac_symbol+str(nb_per_elts[vac_specie])+ \
                                "_subst_"+subspecie_symbol,
                        'unique_site': vac_site,
                        'supercell':{'size':sc_scale,'structure':antiste_sc},
                        'charges':[c-oxid_states[vac_symbol] for c in range(
                            max_min_oxid[subspecie_symbol][0],
                            max_min_oxid[subspecie_symbol][1]+1)]})

        if intrinsic:
            self.defects['vacancies'] = vacancies 
            self.defects['antisites'] = antisites

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
