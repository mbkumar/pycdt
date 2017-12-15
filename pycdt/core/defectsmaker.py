# coding: utf-8
from __future__ import division

"""
Code to generate charged defects structure files.
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
import abc
import math
import logging
from itertools import product


from monty.string import str2unicode
from monty.serialization import dumpfn
from pymatgen.core.structure import PeriodicSite
from pymatgen.core.periodic_table import Element, Specie, get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.analysis.local_env import ValenceIonicRadiusEvaluator as VIRE
from pymatgen.analysis.defects.point_defects import \
        StructureMotifInterstitial


def get_optimized_sc_scale(inp_struct, final_site_no):

    """
    Get the optimal scaling to generate supercells with atoms less than
    the final_site_no.
    """

    if final_site_no < len(inp_struct.sites):
        final_site_no = len(inp_struct.sites)

    dictio={}
    result=[]
    for k1 in range(1,6):
        for k2 in range(1,6):
            for k3 in range(1,6):
                struct = inp_struct.copy()
                struct.make_supercell([k1, k2, k3])
                if len(struct.sites) > final_site_no:
                    continue

                min_dist = 1000.0
                for a in range(-1,2):
                    for b in range(-1,2):
                        for c in range(-1,2):
                            try:
                                distance = struct.get_distance(0, 0, (a,b,c))
                            except:
                                print (a, b, c)
                                raise
                            if  distance < min_dist and distance>0.00001:
                                min_dist = distance
                min_dist = round(min_dist, 3)
                if min_dist in dictio:
                    if dictio[min_dist]['num_sites'] > struct.num_sites:
                        dictio[min_dist]['num_sites'] = struct.num_sites
                        dictio[min_dist]['supercell'] = [k1,k2,k3]
                else:
                    dictio[min_dist]={}
                    dictio[min_dist]['num_sites'] = struct.num_sites
                    dictio[min_dist]['supercell'] = [k1,k2,k3]
    min_dist = -1.0
    biggest = None
    for c in dictio:
        if c > min_dist:
            biggest = dictio[c]['supercell']
            min_dist = c
    if biggest is None or min_dist < 0.0:
        raise RuntimeError('could not find any supercell scaling vector')
    return biggest


class DefectCharger:
    __metaclass__ = abc.ABCMeta
    """
    Abstract base class to define the properties of a defect charge generator
    """
    def __init__(self, structure):
        pass

    @abc.abstractmethod
    def get_charges(self, defect_type, site_specie=None, sub_specie=None):
        """
        Based on the type of defect, site and substitution (if any) species
        the defect charge states are generated.
        Args:
            defect_type (str): Options are vacancy, antisite, substitution,
                               and interstitial
            site_specie: Specie on the host lattice site
                         For interstitials, use this
            sub_specie: Specie that is replacing the site specie.
                        For antisites and substitution defects

        """
        raise NotImplementedError


class DefectChargerSemiconductor(DefectCharger):
    """
    Determine oxidation states from bond valence method (unless oxidation states specified)
    Use BVM or specified oxidation states for site oxidation preference
    Use possible oxidation state list for sub_site elements
    """
    def __init__(self, structure, min_max_oxi=None, oxi_states=None):
        """
        Charge assignment based on the oxidation states referenced from 
        semiconductor database. Targetted materials are shallow and some 
        wideband semiconductors. For these systems, antisites are common and 
        their charge assignment for antisites follows vacancies
        Args: structure
            pymatgen structure object to determine the oxidation states
            min_max_oxi: any user specified min/max oxidation ranges for elements in structure
            oxi_states: any user specified oxidation states of elements in structure
        """
        min_max_oxi = min_max_oxi if min_max_oxi is not None else {}
        oxi_states = oxi_states if oxi_states is not None else {}

        struct_species = structure.types_of_specie
        if (len(struct_species) == 1) and struct_species[0].symbol not in oxi_states.keys():
            oxi_states[struct_species[0].symbol] = 0
        else:
            vir = VIRE(structure)
            for elt, oxi in vir.valences.items():
                strip_key = ''.join([s for s in elt if s.isalpha()])
                if strip_key not in oxi_states.keys():
                    oxi_states[strip_key] = oxi
        self.oxi_states = oxi_states

        self.min_max_oxi_bulk = [0,0]
        for s in struct_species:
            if isinstance(s, Specie):
                el = s.element
            elif isinstance(s, Element):
                el = s
            else:
                continue
            if el.symbol not in min_max_oxi.keys():
                max_oxi = max(el.common_oxidation_states)
                min_oxi = min(el.common_oxidation_states)
                min_max_oxi[el.symbol] = [min_oxi,max_oxi]
            if min_max_oxi[el.symbol][0] < self.min_max_oxi_bulk[0]:
                self.min_max_oxi_bulk[0] = min_max_oxi[el.symbol][0]
            if min_max_oxi[el.symbol][1] > self.min_max_oxi_bulk[1]:
                self.min_max_oxi_bulk[1] = min_max_oxi[el.symbol][1]
        self.min_max_oxi = min_max_oxi

    def get_charges(self, defect_type, site_specie, sub_specie=None):
        """
        Based on the type of defect, site and substitution (if any) species
        the defect charge states are generated.
        Args:
            defect_type (str): Options are vacancy, antisite, substitution,
                               and interstitial
            site_specie (str): Specie on the host lattice site
                         Use this for interstitials as well
            sub_specie (str): Specie that is replacing the site specie.
                        At present used for substitution and antisite defects
        """
        if site_specie not in self.min_max_oxi.keys():
            max_oxi = max(Element(site_specie).common_oxidation_states)
            min_oxi = min(Element(site_specie).common_oxidation_states)
            self.min_max_oxi[site_specie] = [min_oxi, max_oxi]
        if sub_specie:
            if sub_specie not in self.min_max_oxi.keys():
                max_oxi = max(Element(sub_specie).common_oxidation_states)
                min_oxi = min(Element(sub_specie).common_oxidation_states)
                self.min_max_oxi[sub_specie] = [min_oxi, max_oxi]
        if defect_type == 'vacancy':
            site_oxi = self.oxi_states[site_specie]
            if site_oxi:
                return list(range(-abs(site_oxi), abs(site_oxi) + 1))
            else:
                min_oxi = self.min_max_oxi[site_specie][0]
                max_oxi = self.min_max_oxi[site_specie][1]
                if abs(min_oxi) < abs(max_oxi):
                    return list(range(-abs(max_oxi), abs(max_oxi) + 1))
                else:
                    return list(range(-abs(min_oxi), abs(min_oxi) + 1))
        elif defect_type == 'antisite':
            return list(range(self.min_max_oxi_bulk[0],
                              self.min_max_oxi_bulk[1]-1))
        elif defect_type == 'substitution':
            oxi_sub_set = Element(sub_specie).common_oxidation_states
            oxi_site = self.oxi_states[site_specie]
            min_max_oxi_bulk_sub = [min(min(oxi_sub_set)-oxi_site,-1),
                                    max(max(oxi_sub_set)-oxi_site,1)]
            if (min_max_oxi_bulk_sub[1] - min_max_oxi_bulk_sub[0]) > 2:
                if min_max_oxi_bulk_sub[1] > 2:
                    return list(range(min_max_oxi_bulk_sub[0],
                                      min_max_oxi_bulk_sub[1]-2)) # if range exists, less likely to be a higher charge
                else:
                    return list(range(min_max_oxi_bulk_sub[0],
                                      2)) # extend to 2 if upper range does not already include this
            else:
                return list(range(min_max_oxi_bulk_sub[0]-1,
                                  min_max_oxi_bulk_sub[1]+2)) #likely upper bound is 0, so extend to 2
        elif defect_type == 'interstitial':
            if self.min_max_oxi[site_specie][0] > 0:
                min_oxi = 0
            else:
                min_oxi = self.min_max_oxi[site_specie][0]
            if self.min_max_oxi[site_specie][1] < 0:
                max_oxi = 0
            else:
                max_oxi = self.min_max_oxi[site_specie][1]
            return list(range(min_oxi, max_oxi+1))
        else:
            raise ValueError("Defect type not understood")


class DefectChargerInsulator(DefectCharger):
    """
    Conservative charge assignment based on the oxidation statess determined 
    by bond valence. Targetted materials are wideband semiconductors and 
    insulators. AxBy where A is cation and B is anion will have charge 
    assignments {A: [0:y], B:[-x:0]}. For these systems, antisites typically
    have very high formation energies and are ignored.
    """
    def __init__(self, structure):
        """
        Conservative defect charge generator based on the oxidation statess 
        determined by bond valence. Targetted materials are wideband 
        semiconductors and insulators. AxBy where A is cation and B is 
        anion will have charge assignments {A: [0:y], B:[-x:0]}. For these 
        systems, antisites typically have very high formation energies and 
        are ignored.
        Args:
            structure: pymatgen structure object 
        """
        struct_species = structure.types_of_specie
        if len(struct_species) == 1:
            oxi_states = {struct_species[0].symbol: 0}
        else:
            vir = VIRE(structure)
            oxi_states = vir.valences
        self.oxi_states = {}
        for key,val in oxi_states.items():
            strip_key = ''.join([s for s in key if s.isalpha()])
            self.oxi_states[str2unicode(strip_key)] = val

        self.min_max_oxi = {}
        for s in struct_species:
            if isinstance(s, Specie):
                el = s.element
            elif isinstance(s, Element):
                el = s
            else:
                continue
            max_oxi = max(el.common_oxidation_states)
            min_oxi = min(el.common_oxidation_states)
            self.min_max_oxi[str2unicode(el.symbol)] = (min_oxi,max_oxi)
        
    def get_charges(self, defect_type, site_specie=None, sub_specie=None):
        """
        Based on the type of defect, site and substitution (if any) species
        the defect charge states are generated.
        Args:
            defect_type (str): Options are vacancy, antisite, substitution,
                               and interstitial
            site_specie: Specie on the host lattice site 
                         For interstitials, use this
            sub_specie: Specie that is replacing the site specie.
                        For antisites and substitution defects
        """
        if defect_type == 'vacancy':
            vac_symbol = get_el_sp(site_specie).symbol
            vac_oxi_state = self.oxi_states[str2unicode(vac_symbol)]
            if vac_oxi_state < 0:
                min_oxi = max(vac_oxi_state, self.min_max_oxi[vac_symbol][0])
                max_oxi = 0
            elif vac_oxi_state > 0:
                max_oxi = min(vac_oxi_state, self.min_max_oxi[vac_symbol][1])
                min_oxi = 0
            else: # most probably single element
                oxi_states = get_el_sp(site_specie).common_oxidation_states
                min_oxi = min(oxi_states)
                max_oxi = max(oxi_states)
            return [-c for c in range(min_oxi, max_oxi+1)]

        elif defect_type == 'antisite':
            vac_symbol = get_el_sp(site_specie).symbol
            vac_oxi_state = self.oxi_states[str2unicode(vac_symbol)]
            as_symbol = get_el_sp(sub_specie).symbol
            if vac_oxi_state > 0:
                oxi_max = max(self.min_max_oxi[as_symbol][1],0)
                oxi_min = 0
            else:
                oxi_max = 0
                oxi_min = min(self.min_max_oxi[as_symbol][0],0)
            return [c - vac_oxi_state for c in range(
                        oxi_min, oxi_max+1)]
    
        elif defect_type == 'substitution':
            site_specie = get_el_sp(site_specie)
            sub_specie = get_el_sp(sub_specie)
            vac_symbol = site_specie.symbol
            vac_oxi_state = self.oxi_states[str2unicode(vac_symbol)]

            max_oxi_sub = max(sub_specie.common_oxidation_states)
            min_oxi_sub = min(sub_specie.common_oxidation_states)
            if vac_oxi_state > 0:
                if max_oxi_sub < 0:
                    raise ValueError("Substitution seems not possible")
                else:
                    if max_oxi_sub > vac_oxi_state:
                        return list(range(max_oxi_sub - vac_oxi_state + 1))
                    else:
                        return [max_oxi_sub - vac_oxi_state]
            else:
                if min_oxi_sub > 0:
                    raise ValueError("Substitution seems not possible")
                else:
                    if min_oxi_sub < vac_oxi_state:
                        return list(range(min_oxi_sub - vac_oxi_state, 1))
                    else:
                        return [min_oxi_sub - vac_oxi_state]
        
        elif defect_type == 'interstitial':
            site_specie = get_el_sp(site_specie)
            min_oxi = min(min(site_specie.common_oxidation_states), 0)
            max_oxi = max(max(site_specie.common_oxidation_states), 0)

            return list(range(min_oxi, max_oxi+1))


class DefectChargerUserCustom(DefectCharger):
    """
        NOTE from developers:
            This code is for full control over charging approach of defects
            but will be replaced with a future version.
            Code has no unit test and will not be maintained
            going forward (as of 12/15/2017).
            However, we are keeping function here to allow for
            current users to make use of it...

    Determine oxidation states from bond valence method 
    (unless oxidation states specified)
    Then ask user what charges they want
    """
    def __init__(self, structure, oxi_states=None):
        """
        Does initial Valence bond method to initialize oxidation states for user help
        Overwridden by oxi_state specified by
        Args: structure
            pymatgen structure object to determine the oxidation states
            min_max_oxi: any user specified min/max oxidation ranges for elements in structure
            oxi_states: any user specified oxidation states of elements in structure
        """
        oxi_states = oxi_states if oxi_states is not None else {}
        struct_species = structure.types_of_specie
        if (len(struct_species) == 1) and struct_species[0].symbol not in oxi_states.keys():
            oxi_states[struct_species[0].symbol] = 0
        else:
            vir = VIRE(structure)
            for elt, oxi in vir.valences.items():
                strip_key = ''.join([s for s in elt if s.isalpha()])
                if strip_key not in oxi_states.keys():
                    oxi_states[strip_key] = oxi
        self.oxi_states = oxi_states
        print ('\nThis is Full-User Charge Generation Mode.\n'
               'Options are: (1) Range mode (input min and max of each each type of defect) '
               'or (2) Individual mode (input each charge state you want for each defect)\n'
               '\nWhen finished with specific defect, press ENTER to continue.')
        rng_mod = raw_input('Please specify Range (R) or Individual (I) Mode:')
        if 'R' == rng_mod.upper()[0]:
            self.rangemode = True
        else:
            self.rangemode = False


    def get_charges(self, defect_type, site_specie=None, sub_specie=None):
        """
        Based on the type of defect, site and substitution (if any) species
        the defect charge states are generated.
        Args:
            defect_type (str): Options are vacancy, antisite, substitution,
                               and interstitial
            site_specie: Specie on the host lattice site
                         For interstitials, use this
            sub_specie: Specie that is replacing the site specie.
                        For antisites and substitution defects
        """
        if site_specie in self.oxi_states.keys():
            sitechg = self.oxi_states[site_specie]
        else:
            sitechg = False

        if sub_specie:
            if sub_specie in self.oxi_states.keys():
                subchg = self.oxi_states[sub_specie]
            else:
                subchg = False

        def get_users_charges():
            tmpchgs = raw_input('What charges would you like? : ')
            chgs = [int(c) for c in tmpchgs.split()]
            if self.rangemode:
                return list(range(chgs[0], chgs[-1]+1))
            else:
                return list(chgs)

        if defect_type == 'vacancy':
            if not sitechg:
                print (site_specie, defect_type, 
                       'charge suggestion unknown (specify oxidation states to get suggestion)')
            else:
                print (site_specie, defect_type, 'has charge =', -sitechg, 
                       'according to Simple Ionic Theory')
        elif defect_type in ['antisite','substitution']:
            nom = sub_specie+'_on_'+site_specie
            if not sitechg or not subchg:
                print (nom,defect_type,'charge suggestion unknown (specify oxidation states to get suggestion)')
            else:
                print (nom,defect_type,'has charge = ',subchg - sitechg,'according to Simple Ionic Theory')
        elif defect_type == 'interstitial':
            if not sitechg:
                print (site_specie,defect_type,'charge suggestion unknown (specify oxidation states to get suggestion)')
            else:
                print (site_specie,defect_type,'has charge = ',sitechg,'according to Simple Ionic Theory')
        outchgs = get_users_charges()
        print ('    Charges generated:',outchgs)
        return outchgs


class ChargedDefectsStructures(object):
    """
    A class to generate charged defective structures for use in first
    principles supercell formalism. The standard defects such as antisites
    and vacancies are generated.  Interstitial finding is also implemented
    (optional).
    """
    def __init__(self, structure,  max_min_oxi=None, substitutions=None, 
                 oxi_states=None, cellmax=128, antisites_flag=True, 
                 include_interstitials=False, interstitial_elements=None,
                 intersites=None, standardized=False, 
                 struct_type='semiconductor'):
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
            interstitial_elements ([str]):
                List of strings containing symbols of the elements that are
                to be considered for interstitial sites.  The default (None)
                triggers self-interstitial generation,
                given that include_interstitials is True.
            intersites ([PeriodicSite]):
                A list of PeriodicSites in the bulk structure on which we put
                interstitials.  Note that you still have to set flag
                include_interstitials to True in order to make use of this
                manual way of providing interstitial sites.
            standardized (bool):
                If True, use the primitive standard structure as unit cell
                for generating the defect configurations (default is False).
                The primitive standard structure is obtained from the
                SpacegroupAnalyzer class with a symprec of 0.01.
            struct_type (string):
                Options are 'semiconductor' and 'insulator'. If semiconductor 
                is selected, charge states based on database of semiconductors
                is used to assign defect charges. For insulators, defect 
                charges are conservatively assigned. 
        """
        max_min_oxi = max_min_oxi if max_min_oxi is not None else {}
        substitutions = substitutions if substitutions is not None else {}
        oxi_states = oxi_states if oxi_states is not None else {}
        interstitial_elements = interstitial_elements if interstitial_elements is not None else []
        intersites = intersites if intersites is not None else []

        self.defects = []
        self.cellmax = cellmax
        self.substitutions = {}
        self.struct_type = struct_type
        for key, val in substitutions.items():
            self.substitutions[str2unicode(key)] = val

        spa = SpacegroupAnalyzer(structure, symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        if standardized:
            self.struct = prim_struct
        else:
            self.struct = structure

        struct_species = self.struct.types_of_specie
        if self.struct_type == 'semiconductor':
            self.defect_charger = DefectChargerSemiconductor(self.struct,
                                                             min_max_oxi=max_min_oxi)
        elif self.struct_type == 'insulator':
            self.defect_charger = DefectChargerInsulator(self.struct)
        elif self.struct_type == 'manual':
            self.defect_charger = DefectChargerUserCustom(self.struct,
                                                          oxi_states=oxi_states)
        else:
            raise NotImplementedError
        
        if include_interstitials and interstitial_elements:
            for elem_str in interstitial_elements:
                if not Element.is_valid_symbol(elem_str):
                    raise ValueError("invalid interstitial element"
                            " \"{}\"".format(elem_str))

        conv_prim_rat = int(self.struct.num_sites/prim_struct.num_sites)
        sc_scale = get_optimized_sc_scale(self.struct,cellmax)
        self.defects = {}
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        self.defects['bulk'] = {
                'name': 'bulk',
                'supercell': {'size': sc_scale, 'structure': sc}}

        # If interstitials are provided as a list of PeriodicSites,
        # make sure that the lattice has not changed.
        if include_interstitials and intersites:
            smat = sc.lattice.matrix
            for intersite in intersites:
                imat = intersite.lattice.matrix
                print(str(imat))
                print(str(smat))
                for i1 in range(3):
                    for i2 in range(3):
                        if math.fabs(imat[i1][i2]-smat[i1][i2]) > 0.001:
                            raise RuntimeError("Discrepancy between lattices"
                                    " underlying the input interstitials and"
                                    " the bulk structure; possibly because of"
                                    " standardizing the input structure.")

        vacancies = []
        as_defs = []
        sub_defs = []

        vac = Vacancy(self.struct, {}, {})
        vac_scs = vac.make_supercells_with_defects(sc_scale)

        print("Setting up vacancies and antisites...")
        for i in range(vac.defectsite_count()):
            vac_site = vac.get_defectsite(i)
            site_mult = vac.get_defectsite_multiplicity(i)
            site_mult = int(site_mult/conv_prim_rat)
            vac_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol
            vac_sc = vac_scs[i+1]
            vac_sc_site = list(set(vac_scs[0].sites) - set(vac_sc.sites))[0]
            charges_vac = self.defect_charger.get_charges('vacancy',
                                                          vac_symbol)
            vacancies.append({
                'name': "vac_{}_{}".format(i+1, vac_symbol),
                'unique_site': vac_site,
                'bulk_supercell_site': vac_sc_site,
                'defect_type': 'vacancy',
                'site_specie': vac_symbol,
                'site_multiplicity': site_mult,
                'supercell': {'size': sc_scale,'structure': vac_sc},
                'charges': charges_vac})

            if antisites_flag:
                for as_specie in set(struct_species)-set([vac_specie]):
                    charges_as = self.defect_charger.get_charges(
                            'antisite', vac_symbol, as_specie.symbol)
                    as_symbol = as_specie.symbol
                    as_sc = vac_sc.copy()
                    as_sc.append(as_symbol, vac_sc_site.frac_coords)

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
                        'charges': charges_as})

            if vac_symbol in self.substitutions:
                for subspecie_symbol in self.substitutions[vac_symbol]:
                    sub_sc = vac_sc.copy()
                    sub_sc.append(subspecie_symbol, vac_sc_site.frac_coords)

                    charges_sub = self.defect_charger.get_charges(
                            'substitution', vac_symbol, subspecie_symbol)
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
                        'charges':charges_sub})

        self.defects['vacancies'] = vacancies 
        self.defects['substitutions'] = sub_defs
        self.defects['substitutions'] += as_defs

        if include_interstitials:
            print("Searching for interstitial sites (this can take awhile)...")
            interstitials = []
            inter_types = []
            inter_cns = []
            inter_multi = []
            if interstitial_elements:
                inter_elems = interstitial_elements
            else:
                inter_elems = [elem.symbol for elem in \
                        self.struct.composition.elements]
            if len(inter_elems) == 0:
                raise RuntimeError("empty element list for interstitials")
            if not intersites:
                intersites = []
                smi = StructureMotifInterstitial(self.struct, inter_elems[0],
                                                 dl=0.2)
                n_inters = len(smi.enumerate_defectsites())
                for i in range(n_inters):
                    intersites.append(smi.get_defectsite(i))
                    inter_types.append(smi.get_motif_type(i))
                    inter_cns.append(smi.get_coordinating_elements_cns(i))
                    inter_multi.append(int(
                        smi.get_defectsite_multiplicity(i)/conv_prim_rat))

            # Now set up the interstitials.
            for elt in inter_elems:
                for i, intersite in enumerate(intersites):
                    if inter_types and inter_cns:
                        tmp_string = ""
                        for elem, cn in inter_cns[i].items():
                            tmp_string = tmp_string + "_{}{}".format(elem, cn)
                        if tmp_string == "":
                            raise RuntimeError("no coordinating neighbors")
                        name = "inter_{}_{}_{}{}".format(
                                    i+1, elt, inter_types[i], tmp_string)
                        site_mult = inter_multi[i]

                    else:
                        name = "inter_{}_{}".format(i+1, elt)
                        # TODO: This needs better structuring
                        site_mult = int(1.0 / conv_prim_rat)

                    site = PeriodicSite(Element(elt), intersite.frac_coords,
                                        intersite.lattice)
                    site_sc = PeriodicSite(
                            Element(elt), site.coords, sc.lattice,
                            coords_are_cartesian=True)
                    sc_with_inter = sc.copy()
                    sc_with_inter.append(elt, site_sc.frac_coords)

                    charges_inter = self.defect_charger.get_charges(
                            'interstitial', elt)

                    interstitials.append({
                            'name': name,
                            'unique_site': site,
                            'bulk_supercell_site': site_sc,
                            'defect_type': 'interstitial',
                            'site_specie': Element(elt),
                            'site_multiplicity': site_mult,
                            'supercell': {'size': sc_scale,
                                          'structure': sc_with_inter},
                            'charges': charges_inter})

            self.defects['interstitials'] = interstitials

        print("\nNumber of jobs created:")
        tottmp=0
        for j in self.defects.keys():
            if j=='bulk':
                print("    bulk = 1")
                tottmp += 1
            else:
                print("    {}:".format(j))
                for lis in self.defects[j]:
                    print("        {} = {}".format(lis['name'],
                                                   len(lis['charges'])))
                    tottmp += len(lis['charges'])
        print("Total (non dielectric) jobs created = {}\n".format(tottmp))

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

    def to(self, outfile):
        dumpfn(self.defects, outfile)

    def get_n_defects_of_type(self, defect_type):
        """
        Get the number of defects of the given type.
        Args:
            defect_type (str): defect type (vacancies,
                interstitials, substitutions).
        Returns:
            n_defects (int): number of defects of given type.
        """
        try:
            return len(self.defects[defect_type])
        except:
            return 0

    def get_ith_supercell_of_defect_type(self, i, defect_type):
        """
        Get the (i-1)-th defect supercell with the given defect type.
        Args:
            i (int): index (starting from 0) of target supercell.
            defect_type (str): defect type (vacancies,
                interstitials, substitutions).
        Returns:
            sc (Structure): copy of the defect supercell.
        """
        return self.defects[defect_type][i]['supercell'][
            'structure'].copy()


