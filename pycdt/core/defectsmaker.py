# coding: utf-8
from __future__ import division
"""
Code to generate charged defects structure.
Ideas from pydii's code and geoffroy's code are merged.
"""

__author__ = "Bharat Medasani, Geoffroy Hautier, Danny Broberg"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com,geoffroy@uclouvain.be,dbroberg@berkeley.edu"
__status__ = "Development"
__date__ = "Janurary 6, 2016"

import copy

from monty.string import str2unicode
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
    print 'number of sites in bulk cell=',inp_struct.num_sites,\
	  '\nnumber of sites in final super cell=', final_site_no
    target_site = inp_struct.sites[0]
    dictio={}
    result=[]
    for k1 in range(1,6):
        for k2 in range(1,6):
            for k3 in range(1,6):
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
                 antisites_flag=True, standardized=False, 
                 charge_states='liberal'):
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
            antisites_flag: 
                If False, don't generate antisites
            charge_states:
                Options are 'liberal' and 'conservative'. If liberal is selected
                more charge states are computed
        """

        self.defects = []
        self.cellmax = cellmax
        self.substitutions = {}
        self.charge_states = charge_states
        for key,val in substitutions.items():
            self.substitutions[str2unicode(key)] = val

        spa = SpacegroupAnalyzer(structure,symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        if standardized:
            self.struct = prim_struct
        else:
            self.struct = structure

        struct_species = self.struct.types_of_specie
        if not oxi_states:
            if len(struct_species) == 1:
                oxi_states = {self.struct.types_of_specie[0].symbol: 0}
            else:
                vir = ValenceIonicRadiusEvaluator(self.struct)
                oxi_states = vir.valences
        self.oxi_states = {}
        for key,val in oxi_states.items():
            strip_key = ''.join([s for s in key if s.isalpha()])
            self.oxi_states[str2unicode(strip_key)] = val

        print 'oxidation states for bulk=',self.oxi_states


        conv_prim_rat = int(self.struct.num_sites/prim_struct.num_sites)
        sc_scale = get_optimized_sc_scale(self.struct,cellmax)
        self.defects = {}
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        self.defects['bulk'] = {
                'name': 'bulk',
                'supercell': {'size': sc_scale, 'structure': sc}}

        if not max_min_oxi: 
            max_min_oxi = {}
            for s in struct_species:
                if isinstance(s, Specie):
                    el = s.element
                elif isinstance(s, Element):
                    el = s
                else:
                    continue
                max_oxi = max(el.common_oxidation_states)
                min_oxi = min(el.common_oxidation_states)
                max_min_oxi[str2unicode(el.symbol)] = (min_oxi,max_oxi)
            for s, subspecies in self.substitutions.items():
                for subspecie in subspecies:
                    el = Element(subspecie)
                    max_oxi = max(el.common_oxidation_states)
                    min_oxi = min(el.common_oxidation_states)
                    max_min_oxi[str2unicode(el.symbol)] = (min_oxi,max_oxi)
        print 'max/min oxidation states=',max_min_oxi
        self.max_min_oxi = max_min_oxi
	
	if self.charge_states=='liberal': #check that all substitutions exist for all species 
		subelts=[]
		for s, subspecies in self.substitutions.items():
			for j in subspecies:
				subelts.append(j)
		subeltlis=list(set(subelts))
		warnlist=['\nWARNING - because of liberal setting, \
			will make sure all substitution elements are tried on each native element']
		for s,subspecies in self.substitutions.items():
			tmp=[]
			for j in subeltlis:
				if j not in subspecies:
					tmp.append(j)
					self.substitutions[s].append(j)
			if tmp:
				warnlist.append(str(tmp)+' added to substitution list of '+str(s))
		if len(warnlist)!=1:
			for j in warnlist: print j			
		print 'final sub dictionary=',self.substitutions,'\n'

        vacancies = []
        as_defs = []
        sub_defs = []

        vac = Vacancy(self.struct, {}, {})
        vac_scs = vac.make_supercells_with_defects(sc_scale)

        print 'oxidation states = ', self.oxi_states
        for i in range(vac.defectsite_count()):
            vac_site = vac.get_defectsite(i)
            site_mult = vac.get_defectsite_multiplicity(i)
            site_mult = int(site_mult/conv_prim_rat)
            vac_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol
            vac_sc = vac_scs[i+1]
            vac_sc_site = list(set(vac_scs[0].sites) - set(vac_sc.sites))[0]

            list_charges=[]
            print 'vac_symbol=', vac_symbol
            vac_oxi_state = self.oxi_states[str2unicode(vac_symbol)]
            if vac_oxi_state < 0:
                min_oxi = min(vac_oxi_state, self.max_min_oxi[vac_symbol][0])
		if self.charge_states=='liberal':	
			max_oxi = 2
		else:
                	max_oxi = 0
            elif vac_oxi_state > 0:
                max_oxi = max(vac_oxi_state, self.max_min_oxi[vac_symbol][1])
		if self.charge_states=='liberal':
			min_oxi = -2
		else:
                	min_oxi = 0
            for c in range(min_oxi, max_oxi+1):
                list_charges.append(-c)
            print 'charge states for ',vac_symbol,' vacancy =',list_charges

            vacancies.append({
                'name': "vac_{}_{}".format(i+1, vac_symbol),
                'unique_site': vac_site,
                'bulk_supercell_site': vac_sc_site,
                'defect_type': 'vacancy',
                'site_specie': vac_symbol,
                'site_multiplicity': site_mult,
                'supercell': {'size': sc_scale,'structure': vac_sc},
                'charges': list_charges })

            # Antisite defects generation

            if antisites_flag:
                for as_specie in set(struct_species)-set([vac_specie]):
                    as_symbol = as_specie.symbol
                    as_sc = vac_sc.copy()
                    as_sc.append(as_symbol, vac_sc_site.frac_coords)
                    if vac_oxi_state > 0:
                        oxi_max = max(self.max_min_oxi[as_symbol][1],0)
                        oxi_min = 0
                    else:
                        oxi_max = 0
                        oxi_min = min(self.max_min_oxi[as_symbol][0],0)
                    if self.charge_states=='liberal' and oxi_min==oxi_max:
                        if oxi_min - vac_oxi_state > 0:
                            charges = list(range(-1,oxi_min-vac_oxi_state+1))
                        else:
                            charges = list(range(oxi_min-vac_oxi_state-1,1))
                    else:
                        charges = [c - vac_oxi_state for c in range(
                            oxi_min, oxi_max+1)]
		    print 'charges for ',as_symbol,' on ', \
				vac_symbol,'=',charges

                    as_defs.append({
                        'name': "as_{}_{}_on_{}".format(
                            i+1, as_symbol, vac_symbol),
                        'unique_site': vac_site,
                        'bulk_supercell_site': vac_sc_site,
                        'defect_type': 'antisite',
                        'site_specie': vac_symbol,
                        'substitution_specie': as_symbol,
                        'site_multiplicity': site_mult,
                        'supercell': {'size': sc_scale,'structure': as_sc},
                        'charges': charges})

            # Substitutional defects generation
	    if vac_symbol in self.substitutions:
                for subspecie_symbol in self.substitutions[vac_symbol]:
                    sub_sc = vac_sc.copy()
                    sub_sc.append(subspecie_symbol, vac_sc_site.frac_coords)
                    if vac_oxi_state > 0:
                        oxi_max = max(self.max_min_oxi[subspecie_symbol][1],0)
                        oxi_min = 0
                    else:
                        oxi_max = 0
                        oxi_min = min(self.max_min_oxi[subspecie_symbol][0],0)
                    if self.charge_states=='liberal' and oxi_min==oxi_max:
                        if oxi_min - vac_oxi_state > 0:
                            charges = list(range(-1,oxi_min-vac_oxi_state+1))
                        else:
                            charges = list(range(oxi_min-vac_oxi_state-1,1))
                    else:
                        charges = [c - vac_oxi_state for c in range(
                            oxi_min, oxi_max+1)]
                    print 'charges for ',subspecie_symbol,' on',vac_symbol, \
				' substitution=',charges

                    sub_defs.append({
                        'name': "sub_{}_{}_on_{}".format(
                            i+1, subspecie_symbol, vac_symbol),
                        'unique_site': vac_site,
                        'bulk_supercell_site': vac_sc_site,
                        'defect_type':'antisite',
                        'site_specie':vac_symbol,
                        'substitution_specie':subspecie_symbol,
                        'site_multiplicity':site_mult,
                        'supercell':{'size':sc_scale,'structure':sub_sc},
                        'charges':charges})

        self.defects['vacancies'] = vacancies 
        self.defects['substitutions'] = sub_defs
        self.defects['substitutions'] += as_defs

        #interstitials
	#if (not self.interstitial_sites and self.charge_states=='liberal'):
		
	##THIS IS STRICTLY FOR DANNY TESTING INTERSTITIAL AUTOMATION
	#from pymatgen.analysis.defects.point_defects import Interstitial as Inter
	#from pymatgen.matproj.rest import MPRester
	#from pymatgen.core.periodic_table import Element
	#mp =MPRester()
	#struct=mp.get_structure_by_material_id('mp-2231') 
	#valdic={'Sn':2,'S':-2}
	#radi={}
	#for j in valdic.keys():
	#radi[j]=Element(j).atomic_radius
	#s=Inter(struct,valdic,radi)

        interstitials = []
        for elt in self.struct.composition.elements:
            count = 1
            for frac_coord in interstitial_sites:
                site = PeriodicSite(elt, frac_coord, structure.lattice)
                interstitials.append({
                    'name': elt.symbol+str(count)+"_inter",
                    'unique_site': site,
                    'supercell': {'size': s_size,
                        'structure': self.make_interstitial(site, sc_scale)},
                    'charges': [c for c in range(
                        max_min_oxi[elt][0], max_min_oxi[elt][1]+1)]})
                count = count+1
        self.defects['interstitials'] = interstitials

	print '\nNumber of jobs created:'
	tottmp=0
	for j in self.defects.keys():
		if j=='bulk':
			print 'bulk'
			tottmp+=1
		else:
			print j
			for lis in self.defects[j]:
				print '   ',lis['name'],'=',len(lis['charges'])
				tottmp+=len(lis['charges'])
	print 'Total (non dielectric) jobs created = ',tottmp,'\n'

    def make_defect_complexes(max_complex_size=0, include_vacancies=True):
        """
        Function to generate defect complexes
        Args:
            max_complex_size: max. number of defects in a complex.
                If zero, the max size possible is considered based 
                on no. of subsitutions
            include_vacancies: Include vacancies in the defect complex
        """
        if not max_complex_size:
            max_complex_size = len(self.defects['substitutions'])
            if include_vacancies:
                max_complex_size += 1

        complexes = []
        for size in range(2, max_complex_size+1):
            continue
            
    def make_interstitial(self, target_site, sc_scale):
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        sc.append(target_site.specie, target_site.frac_coords)
        
        return sc
