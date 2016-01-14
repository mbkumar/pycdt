#made this for recreating SnS defects run in paper by Malone 2014
__author__ = 'Danny'

from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.matproj.rest import MPRester

import sys
sys.path.append("../")
from pycdcd.core.defectsmaker import ChargedDefectsStructures
from pycdcd.utils.vasp import make_vasp_defect_files, \
        make_vasp_dielectric_files
import sys
sys.path.append("../")

def example_maker(mpid='mp-2231',subs={}):
	print('connect to MP database')
	mp =MPRester()
	struct=mp.get_structure_by_material_id(mpid)  
	print('structure retrieved')
        make_vasp_dielectric_files(struct)
	make_vasp_defect_files(ChargedDefectsStructures(struct, substitutions=subs,charge_states='conservative').defects, struct.composition.reduced_formula)
	#make_vasp_defect_files(ChargedDefectsStructures(struct,charge_states='conservative').defects,struct.composition.reduced_formula)
	#make_vasp_defect_files(ChargedDefectsStructures(struct).defects,struct.composition.reduced_formula)
	#deflis=['N','P','As','Sb','O','Cl','Cu','Na','In','Cd','Zn']
        #make_vasp_defect_files(ChargedDefectsStructures(struct, substitutions={'Sn':deflis,'S':['N','P','As','Sb','O','Cl','Cu']}).defects, struct.composition.reduced_formula)

if __name__ == '__main__':
    mpvals=['mp-2133','mp-10695','mp-1190','mp-2176','mp-1132','mp-13031','mp-804']
    subvals=[{'Zn':['N','P'],'O':[]},{'Zn':[],'S':[]},{'Zn':[],'Se':[]},{'Zn':[],'Te':[]},
		{'Cd':[],'O':[]},{'Mg':[],'Se':[]},{'Ga':['C'],'N':[]}]
    for i in range(len(mpvals)):
	print '\n',i
    	example_maker(mpid=mpvals[i],subs=subvals[i])
