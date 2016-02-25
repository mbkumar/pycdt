# coding: utf-8
from __future__ import division

"""
Code to generate charged defects structure.
Ideas from pydii's code and geoffroy's code are merged.
"""

__author__ = "Bharat Medasani, Geoffroy Hautier, Danny Broberg," + \
        " Nils E. R. Zimmermann"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com, geoffroy@uclouvain.be," + \
        " dbroberg@berkeley.edu, n.zimmermann@tuhh.de"
__status__ = "Development"
__date__ = "Janurary 6, 2016"

import copy

from monty.string import str2unicode
from pymatgen.core.structure import PeriodicSite
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.point_defects import \
        ValenceIonicRadiusEvaluator
try:
    from pymatgen.analysis.defects.alt_interstitial_class import \
            StructureMotifInterstitial
    gen_inter = True
except:
    gen_inter = False

def get_sc_scale(inp_struct, final_site_no):

    """
    Get the scaling to generate supercells with atoms less than the
    final_site_no.
    """

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

    """
    Get the optimal scaling to generate supercells with atoms less than
    the final_site_no.
    """

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
                min_dist = 1000.0
                for a in range(-1,2):
                    for b in range(-1,2):
                        for c in range(-1,2):
                            distance = struct.get_distance(index,index,(a,b,c))
                            if  distance < min_dist and distance>0.00001:
                                min_dist = distance
                min_dist = round(min_dist, 3)
                if dictio.has_key(min_dist):
                    if dictio[min_dist]['num_sites'] > struct.num_sites:
                        dictio[min_dist]['num_sites'] = struct.num_sites
                        dictio[min_dist]['supercell'] = [k1,k2,k3]
                else:
                    dictio[min_dist]={}
                    dictio[min_dist]['num_sites'] = struct.num_sites
                    dictio[min_dist]['supercell'] = [k1,k2,k3]
    min_dist=-1.0
    biggest=None
    for c in dictio:
        if c>min_dist:
            biggest=dictio[c]['supercell']
            min_dist=c
    return biggest


class ChargedDefectsStructures(object):

    """
    A class to generate charged defective structures for use in first
    principles supercell formalism. The standard defects such as antisites
    and vacancies are generated.  Interstitial finding is also implemented
    (optional).
    """

    def __init__(self, structure, max_min_oxi={}, substitutions={},
                 oxi_states={}, cellmax=128, antisites_flag=True,
                 include_interstitials=False, intersites=[],
                 standardized=False, charge_states='liberal'):

        """
        Args:
            structure (Structure):
                the bulk structure.
            max_min_oxi (dict):
                The minimal and maximum oxidation state of each element as a
                dict. For instance {"O":(-2,0)}. If not given, the oxi-states
                of pymatgen are considered.
            substitutions (dict):
                The allowed substitutions of elements as a dict. If not given,
                intrinsic defects are computed. If given, intrinsic (e.g.,
                anti-sites) and extrinsic are considered explicitly specified.
                Example: {"Co":["Zn","Mn"]} means Co sites can be substituted
                by Mn or Zn.
            oxi_states (dict):
                The oxidation state of the elements in the compound e.g.
                {"Fe":2,"O":-2}. If not given, the oxidation state of each
                site is computed with bond valence sum. WARNING: Bond-valence
                method can fail for mixed-valence compounds.
            cellmax (int):
                Maximum number of atoms allowed in the supercell.
            antisites_flag (bool):
                If False, don't generate antisites.
            include_interstitials (bool):
                If true, do generate interstitial defect configurations
                (default: False).
            intersites ([PeriodicSite]):
                A list of PeriodicSites in the bulk structure on which we put
                an interstitial.  Note that you have to set flag
                include_interstitials to True in order to make use of this
                manual way of providing interstitial sites.
            standardized (bool):
                If True, use the primitive standard structure as unit cell
                for generating the defect configurations (default is False).
                The primitive standard structure is obtained from the
                SpacegroupAnalyzer class with a symprec of 0.01.
            charge_states (string):
                Options are 'liberal' and 'conservative'. If liberal is selected,
                more charge states are computed.
        """

        self.defects = []
        self.cellmax = cellmax
        self.substitutions = {}
        self.charge_states = charge_states
        for key,val in substitutions.items():
            self.substitutions[str2unicode(key)] = val

        spa = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        if standardized:
            self.struct = prim_struct
        else:
            self.struct = structure

        # If interstitials are provided as a list of PeriodicSites,
        # make sure that the lattice has not changed.
        if include_interstitials and intersites:
            smat = self.struct.lattice.matrix
            for intersite in intersites:
                imat = intersite.lattice.matrix
                for i1 in range(3):
                    for i2 in range(3):
                        if fabs(imat[i1][i2]-smat[i1][i2])/fabs(
                                imat[i1][i2]) > 1.0e-4:
                            raise RuntimeError("Discrepancy between lattices"
                                    " underlying the input interstitials and"
                                    " the bulk structure; possibly because of"
                                    " standardizing the input structure.")

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
	
	if self.charge_states=='liberal': 
        #check that all substitutions exist for all species 
		subelts=[]
		for s, subspecies in self.substitutions.items():
			for j in subspecies:
				subelts.append(j)
		subeltlis=list(set(subelts))
		warnlist=['\nWARNING - because of liberal setting, \
			all substitution elements are tried on each native element']
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
	#if (not self.interstitial_frac_coords and self.charge_states=='liberal'):
		
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
        # Find interstitial sites if intersites is empty,
        # thus, not providing any input fractional coordinates for
        # interstitial positions.
        # Note that we do want to use the smaller self.struct unit cell
        # and not the supercell sc to search for interstitial sites.
        # However, this approach tacitly assumes that both structures
        # directly map because the interstitial positions found will
        # be placed into the supercell sc.
        # We use the first element in the Composition object underlying
        # our input structure, but the result of the interstitial fin
        if include_interstitials:
            interstitials = []
            inter_types = []
            inter_cns = []
            inter_multi = []
            if not intersites and gen_inter:
                intersites = []
                smi = StructureMotifInterstitial(
                        self.struct,
                        self.struct.composition.elements[0].symbol,
                        dl=0.2)
                n_inters = len(smi.enumerate_defectsites())
                for i_inter in range(n_inters):
                    intersites.append(
                            smi.get_defectsite(i_inter))
                    inter_types.append(smi.get_motif_type(i_inter))
                    inter_cns.append(smi.get_coordinating_elements_cns(i_inter))
                    inter_multi.append(int(smi.get_defectsite_multiplicity(
                            i_inter)/conv_prim_rat))

            # For now, we focus on intrinsic interstitials;
            # extrinsic interstitials are, however,
            # not a dramatic extension.
            for elt in self.struct.composition.elements:
                for i_inter, intersite in enumerate(intersites):
                    if inter_types and inter_cns:
                        tmp_string = ""
                        for elem, cn in inter_cns[i_inter].items():
                            tmp_string = tmp_string + "_{}{}".format(elem, cn)
                        if tmp_string == "":
                            raise RuntimeError("no coordinating neighbors")
                        name = "inter_{}_{}_{}{}".format(i_inter+1, elt.symbol, inter_types[i_inter],
                                tmp_string)
                        site_mult = inter_multi[i_inter]

                    else:
                        name = "inter_{}_{}".format(i_inter+1, elt.symbol)
                        # This needs further attention at some point.
                        site_mult = int(1 / conv_prim_rat)

                    site = PeriodicSite(elt, intersite.frac_coords,
                            intersite.lattice)
                    site_sc = PeriodicSite(elt, site.coords, sc.lattice,
                            coords_are_cartesian=True)
                    sc_with_inter = sc.copy()
                    sc_with_inter.append(elt.symbol,
                        site_sc.frac_coords)

                    charges=[]
                    print 'inter_symbol=', elt
                    inter_oxi_state = self.oxi_states[str2unicode(elt.symbol)]
                    if inter_oxi_state < 0:
                        min_oxi = min(inter_oxi_state, self.max_min_oxi[elt.symbol][0])
                        if self.charge_states=='liberal':
                                max_oxi = 2
                        else:
                                max_oxi = 0
                    elif inter_oxi_state > 0:
                        max_oxi = max(inter_oxi_state, self.max_min_oxi[elt.symbol][1])
                        if self.charge_states=='liberal':
                                min_oxi = -2
                        else:
                                min_oxi = 0
                    for c in range(min_oxi, max_oxi+1):
                        charges.append(-c)
                    print 'charge states for ',elt.symbol,' interstitial =',charges

                    interstitials.append({
                            'name': name,
                            'unique_site': site,
                            'bulk_supercell_site': site_sc,
                            'defect_type': 'interstitial',
                            'site_specie': elt,
                            'site_multiplicity': site_mult,
                            'supercell': {'size': sc_scale, 'structure': sc_with_inter},
                            'charges': charges})


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

        """
        Function to generate a supercell that contains an
        interstitial site.
        Args:
            target_site (PeriodicSite): interstitial site
                to be inserted into a supercell of a
                copy of self.struct.
            sc_scale (3x3 matrix): supercell scaling matrix
                to be applied on the copy of self.struct.
        Returns:
            sc (Structure): supercell containing an
                interstitial site.
        """

        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        sc.append(target_site.specie, target_site.frac_coords)
        
        return sc
