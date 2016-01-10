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

def example_maker():
	print('connect to MP database')
	mp =MPRester()
	struct=mp.get_structure_by_material_id('mp-2231')  #get SnS structure
	print('structure retrieved')
        make_vasp_dielectric_files(struct)
	#make_vasp_defect_files(ChargedDefectsStructures(struct,charge_states='conservative').defects,struct.composition.reduced_formula)
	#make_vasp_defect_files(ChargedDefectsStructures(struct).defects,struct.composition.reduced_formula)
	deflis=['N','P','As','Sb','O','Cl','Cu','Na','In','Cd','Zn']
        make_vasp_defect_files(ChargedDefectsStructures(struct, substitutions={'Sn':deflis,'S':['N','P','As','Sb','O','Cl','Cu']}).defects, struct.composition.reduced_formula)

if __name__ == '__main__':
    example_maker()
