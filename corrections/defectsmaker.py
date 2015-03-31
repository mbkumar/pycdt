#!/usr/bin/env python

__author__ = "Bharat Medasani, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com,geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.structure_modifier import SupercellMaker, StructureEditor
from pymatgen.transformations.defect_transformation import \
        VacancyTransformation,AntisiteDefectTransformation, \
        SubstitutionDefectTransformation
from pymatgen.core.periodic_table import DummySpecie, Element
from pymatgen.core.structure import PeriodicSite


class ChargedDefectsMaker(object):
    """
    A class to generate charged defective structures for use in first 
    principles supercell formalism. The standard defects such as antisites, 
    vacancies are generated.
    TODO: develop a better way to find interstitials
    """
    def __init__(self, structure, max_min_oxid=None, allowed_subst=None, 
            oxid_states=None, cellmax=128, interstitial_sites=[]):
        """
        Args:
            structure:
                the bulk structure
            max_min_oxid:
                The minimal and maximum oxidation state of each element as a 
                dict. For instance {"O":(-2,0)}
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

        #if standardized:
        #    finder = SymmetryFinder(structure,symprec=1e-2)
        #    structure=finder.get_primitive_standard_structure()
        finder=SymmetryFinder(structure,symprec=1e-2)
        self.struct = finder.get_symmetrized_structure()
        self.defects = []
        self.size_limit = size_limit
        nb_per_elts = {e:0 for e in structure.composition.elements}
        for s in self.struct.equivalent_sites:
            list_charges=[]
            for c in range(max_min_oxid[s[0].specie][0], max_min_oxid[s[0].specie][1]+1):
                list_charges.append(-c)
            nb_per_elts[s[0].specie] = nb_per_elts[s[0].specie]+1
            self.defects.append({'short_name': s[0].specie.symbol+str(nb_per_elts[s[0].specie])+"_vac",
                'unique_sites': s[0],
                'supercells':[{'size': s_size,'structure':self.make_defect_cell_vacancy(s[0], s_size)}
                    for s_size in self.get_optimized_supercell(s[0])],
                'charges':list_charges})

            if s[0].specie in allowed_subst:
                for subst in allowed_subst[s[0].specie]:
                    self.defects.append({'short_name':s[0].specie.symbol+str(nb_per_elts[s[0].specie])+"_subst_"
                        +subst.symbol, 'unique_sites':s[0],
                        'supercells':[{'size':s_size,'structure':
                            self.make_defect_cell_intrinsic_subst_defects(s[0], subst, s_size)}
                            for s_size in self.get_optimized_supercell(s[0])],
                        'charges':[c-oxid_states[s[0].specie] for c in range(max_min_oxid[subst][0],
                            max_min_oxid[subst][1]+1)]})



                        #interstitials
        for elt in self.struct.composition.elements:
            count = 1
            for frac_coord in interstitial_sites:
                self.defects.append({'short_name':elt.symbol+str(count)+"_inter",
                    'unique_sites':PeriodicSite(elt, frac_coord, structure.lattice),
                    'supercells':[{'size':s_size,
                        'structure':self.make_defect_interstitial(
                            PeriodicSite(elt, frac_coord, structure.lattice), s_size)}
                        for s_size in self.get_optimized_supercell(s[0])],
                    'charges':[c for c in range(max_min_oxid[elt][0],max_min_oxid[elt][1]+1)]})
                count = count+1

    def make_defect_cell_vacancy(self, target_site, supercell):
        """
        Create a supercell for a vacancy
        TODO: This method needs to be changed to something smarter
        """
        k1=supercell[0]
        k2=supercell[1]
        k3=supercell[2]
        struct = SupercellMaker(self.struct, ((k1, 0, 0), (0, k2, 0), (0, 0, k3))).modified_structure
        editor = StructureEditor(struct)
        index = None
        for i in range(struct.num_sites):
            s = struct._sites[i]
            if s.distance_from_point(target_site.coords)<0.001:
                index=i
        editor.replace_site(index, DummySpecie())
        editor.remove_species([DummySpecie()])
        return editor.modified_structure
    
    def make_defect_interstitial(self, target_site, supercell):
        k1=supercell[0]
        k2=supercell[1]
        k3=supercell[2]
        struct=SupercellMaker(self.struct,((k1,0,0),(0,k2,0),(0,0,k3))).modified_structure
        editor=StructureEditor(struct)
        editor.append_site(target_site.specie, target_site.coords, coords_are_cartesian=True)
        
        return editor.modified_structure
        
    
    def make_defect_cell_intrinsic_subst_defects(self, target_site, subst, supercell):
        """
        subst defines elements that you allow on this site
        charge is the original oxidation state of the element on the target site
        """
        k1=supercell[0]
        k2=supercell[1]
        k3=supercell[2]
        struct=SupercellMaker(self.struct,((k1,0,0),(0,k2,0),(0,0,k3))).modified_structure
        editor=StructureEditor(struct)
        index=None
        for i in range(struct.num_sites):
            s=struct._sites[i]
            if s.distance_from_point(target_site.coords)<0.001:
                index=i
        editor.replace_site(index, subst)
        #editor.remove_species([DummySpecie()])
        return editor.modified_structure   
        
    def get_optimized_supercell(self, target_site, only_biggest=True):
        #print SupercellMaker(struct,((2,1,0),(0,1,0),(0,0,1))).modified_structure
        dictio={}
        result=[]
        for k1 in range(1,3):
            for k2 in range(1,3):
                for k3 in range(1,3):
                    struct=SupercellMaker(self.struct,((k1,0,0),(0,k2,0),(0,0,k3))).modified_structure
                    if len(struct.sites) > self.size_limit:
                        continue
                    #print struct
                    site_target=None
                    index=None
                    for i in range(struct.num_sites):
                        s=struct._sites[i]
                        #print s.coords
                        #print target_site.coords
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
                    print str(k1)+"x"+str(k2)+"x"+str(k3)+" "+str(struct.num_sites)+" "+str(min)
                    #print index
        if only_biggest==True:
            min=-1.0
            biggest=None
            for c in dictio:
                #print str(c)+" "+str(dictio[c]['num_sites'])+" "+str(dictio[c]['supercell'])
                if c>min:
                    biggest=dictio[c]['supercell']
                    min=c
            print min, biggest
            return [biggest]
        else:
            return [dictio[c]['supercell'] for c in dictio]

