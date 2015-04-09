#!/usr/bin/env python

"""
TODO create a VaspInputSet instead?
"""

__author__ = "Geoffroy Hautier, Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.io.vaspio.vasp_input import Kpoints
from pymatgen.io.vaspio_set import MPVaspInputSet
import json
import os

def make_vasp_defect_files(defects, path_base, user_settings=None, 
                           hse=False):
    """
    Generates VASP files for defect computations
    Args:
        defect_structs:
            the defects data as a dictionnary. Ideally this is generated
            from core.defectsmaker.ChargedDefectsStructures.
        path_base:
            where we write the files
        user_settings:
            Settings in dict format to override the defaults used in 
            generating vasp files. The format of the dictionary is
            {'incar':{...},
             'kpoints':...}
        hse:
            hse run or not
    """
    count=1
    bulk_sys = defects['bulk']
    comb_defs = reduce(lambda x,y: x+y, [
        defects[key] for key in defects if key != 'bulk'])

    for defect in comb_defs:
        for charge in defect['charges']:
            print charge
            #for s in defect['supercell']:
            s = defect['supercell']
            print type(s['structure'])
            dict_transf={
                    'defect_type': defect['name'], 
                    'defect_site': defect['unique_site'].as_dict, 
                    'charge': charge, 'supercell': s['size']}
            dict_params=MPVaspInputSet().get_all_vasp_input(s['structure'])
            incar=dict_params['INCAR']
            incar['IBRION']=2
            incar['ISIF']=0
            incar['ISPIN']=1
            incar['LWAVE']=False
            #incar['EDIFF']=0.0001
            incar['EDIFF']=0.001
            incar['ISMEAR']=0
            incar['SIGMA']=0.05
            incar['LVHAR']=True
            incar['LORBIT']=11
            incar['ALGO']="Fast"
            #incar['ALGO']="Normal"
            if hse == True:
                incar['LHFCALC']=True
                incar["ALGO"]="All"
                incar["HFSCREEN"]=0.2
                incar["PRECFOCK"]="Fast"
                incar["AEXX"]=0.45
            kpoint=Kpoints.monkhorst_automatic()
            path=os.path.join(path_base, 
                    str(site['short_name']),"charge"+str(charge))
            os.makedirs(path)
            #f=open(path+"/transformations.json",'w')
            #f.write(json.dumps(dict_transf))
            dumpfn(dict_transf,os.path.join(path,'transformations.json'))
            comp=s['structure'].composition
            sum_elec=0
            elts=set()
            for p in dict_params['POTCAR']:
                if p.element not in elts:
                    sum_elec+=comp.to_dict[p.element]*p.nelectrons
                    elts.add(p.element)
                #print p.element
                #print comp.to_dict[p.element]
                #print p.valence
            print sum_elec
            if charge != 0:
                incar['NELECT']=sum_elec-charge
            dict_params['POTCAR'].write_file(path+"POTCAR")
            incar.write_file(os.path.join(path,"INCAR"))
            kpoint.write_file(os.path.join(path,"KPOINTS"))
            dict_params['POSCAR'].write_file(os.path.join(path,"POSCAR"))
            #print Poscar(s['structure'])
            #Poscar(s['structure']).write_file(path+"POSCAR")
            count=count+1
            #f.close()
