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

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPVaspInputSet
from monty.serialization import loadfn, dumpfn
from monty.json import MontyDecoder, MontyEncoder

#import json
import os

def make_vasp_defect_files(defects, path_base, user_settings={}, hse=False):
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
            {'defects:{'INCAR':{...},'KPOINTS':{...},
             'bulk':{'INCAR':{...},'KPOINTS':{...}}
        hse:
            hse run or not
    """
    bulk_sys = defects['bulk']['supercell']
    comb_defs = reduce(lambda x,y: x+y, [
        defects[key] for key in defects if key != 'bulk'])

    for defect in comb_defs:
        #print type(defect)
        #print defect['charges']
        for charge in defect['charges']:
            s = defect['supercell']
            dict_transf = {
                    'defect_type': defect['name'], 
                    'defect_site': defect['unique_site'], 
                    'defect_supercell_site': defect['bulk_supercell_site'],
                    'charge': charge, 'supercell': s['size']}
            if 'substitution_specie' in  defect:
                dict_transf['substitution_specie'] = defect['substitution_specie']

            dict_params = MPVaspInputSet().get_all_vasp_input(s['structure'])
            incar = dict_params['INCAR']
            incar.update({
                'IBRION': 2, 'ISIF': 2, 'ISPIN': 2, 'LWAVE':False, 
                'EDIFF': 1e-5, 'EDIFFG': -1e-2, 'ISMEAR': 0, 'SIGMA': 0.05,
                'LVTOT': True, 'LVHAR': True, 'LORBIT': 11, 'ALGO': "Fast",
                'ISYM':0})
            #if hse == True:
            #    incar.update({
            #        'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2,
            #        "PRECFOCK": "Fast", 'NKRED': 2})#, "AEXX": 0.45})
            if user_settings:
                if 'INCAR' in user_settings.get('defects',None):
                    incar.update(user_settings['defects']['INCAR'])

            comp = s['structure'].composition
            sum_elec = 0
            elts = set()
            for p in dict_params['POTCAR']:
                if p.element not in elts:
                    sum_elec += comp.as_dict()[p.element]*p.nelectrons
                    elts.add(p.element)
            if charge != 0:
                incar['NELECT'] = sum_elec - charge

            kpoint=dict_params['KPOINTS'].monkhorst_automatic()

            path=os.path.join(path_base,defect['name'],"charge_"+str(charge))
            try:
                os.makedirs(path)
            except:
                pass
            incar.write_file(os.path.join(path,"INCAR"))
            kpoint.write_file(os.path.join(path,"KPOINTS"))
            dict_params['POSCAR'].write_file(os.path.join(path,"POSCAR"))
            dict_params['POTCAR'].write_file(os.path.join(path,"POTCAR"))
            dumpfn(dict_transf, os.path.join(path,'transformation.json'),
                   cls=MontyEncoder)
            if hse:
                incar.update({"LWAVE": True})
                incar.write_file(os.path.join(path,"INCAR.gga"))
                incar.update({
                    'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2, 
                    "PRECFOCK": "Fast", 'NKRED': 2}) #"AEXX": 0.45, 
                #incar.write_file(os.path.join(path,"INCAR.hse"))
                incar.write_file(os.path.join(path,"INCAR.hse1"))
                del incar['PRECFOCK']
                del incar['NKRED']
                incar.write_file(os.path.join(path,"INCAR.hse2"))

    # Generate bulk supercell inputs
    s = bulk_sys
    dict_transf = {'defect_type': 'bulk', 'supercell': s['size']}

    dict_params = MPVaspInputSet().get_all_vasp_input(s['structure'])
    incar = dict_params['INCAR']
    incar.update({
        'IBRION': -1, "NSW": 0, 'ISPIN': 2, 'LWAVE': False, 'EDIFF': 1e-5,
        'ISMEAR': 0, 'SIGMA': 0.05, 'LVTOT': True, 'LVHAR': True, 
        'ALGO': 'Fast', 'ISYM': 0})
    if user_settings:
        if 'INCAR' in user_settings.get('bulk',None):
            incar.update(user_settings['bulk']['INCAR'])
    #if hse == True:
    #    incar.update({
    #        'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2, "AEXX": 0.45, 
    #        "PRECFOCK": "Fast", 'NKRED': 2})
    kpoint=dict_params['KPOINTS'].monkhorst_automatic()
    path = os.path.join(path_base,'bulk')
    try:
        os.makedirs(path)
    except:
        pass
    kpoint.write_file(os.path.join(path,"KPOINTS"))
    dict_params['POSCAR'].write_file(os.path.join(path,"POSCAR"))
    dict_params['POTCAR'].write_file(os.path.join(path,"POTCAR"))
    dumpfn(dict_transf, os.path.join(path,'transformation.json'),
           cls=MontyEncoder)
    if hse:
        incar.update({"LWAVE": True})
        incar.write_file(os.path.join(path,"INCAR.gga"))
        incar.update({
            'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2, #"AEXX": 0.45, 
            "PRECFOCK": "Fast", 'NKRED': 2})
        incar.write_file(os.path.join(path,"INCAR.hse1"))
        del incar['PRECFOCK']
        del incar['NKRED']
        incar.write_file(os.path.join(path,"INCAR.hse2"))
    else:
        incar.write_file(os.path.join(path,"INCAR"))

def make_vasp_defect_files_dos(defects, path_base, user_settings={}, 
                           hse=False, dos_limits=(-1,7)):
    """
    Generates VASP files for defect computations which include dos 
    generation. Useful when the user don't want to use MPWorks for 
    dos calculations.
    Args:
        defects:
            the defects data as a dictionnary. Ideally this is generated
            from core.defectsmaker.ChargedDefectsStructures.
        path_base:
            where we write the files
        user_settings:
            Settings in dict format to override the defaults used in 
            generating vasp files. The format of the dictionary is
            {'defects:{'INCAR':{...},'KPOINTS':{...},
             'bulk':{'INCAR':{...},'KPOINTS':{...}}
        hse:
            hse run or not
        dos_limits:
            Lower and upper limits for dos plot as a tuple. The default
            (-1,7) should work for most of the cases.
    """
    bulk_sys = defects['bulk']['supercell']
    comb_defs = reduce(
            lambda x,y: x+y, 
            [defects[key] for key in defects if key != 'bulk'])

    for defect in comb_defs:
        print type(defect)
        print defect['charges']
        for charge in defect['charges']:
            s = defect['supercell']
            dict_transf = {
                    'defect_type': defect['name'], 
                    'defect_site': defect['unique_site'], 
                    'defect_supercell_site': defect['bulk_supercell_site'],
                    'charge': charge, 'supercell': s['size']}
            if 'substitution_specie' in  defect:
                dict_transf['substitution_specie'] = defect['substitution_specie']


            dict_params = MPVaspInputSet().get_all_vasp_input(s['structure'])
            incar = dict_params['INCAR']
            incar.update({
                'IBRION': 2, 'ISIF': 2, 'ISPIN': 2, 'LWAVE': False, 
                'EDIFF': 1e-5, 'EDIFFG': -1e-2, 'ISMEAR': 0, 'SIGMA': 0.05, 
                'LVTOT': True, 'LVHAR': True, 'LORBIT': 11, 'ALGO': "Fast",
                'ISYM': 0})
            #if hse == True:
            #    incar.update({
            #        'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2,
            #        "PRECFOCK": "Fast", 'NKRED': 2})#, "AEXX": 0.45})
            if user_settings:
                if 'INCAR' in user_settings.get('defects', None):
                    incar.update(user_settings['defects']['INCAR'])

            comp=s['structure'].composition
            sum_elec=0
            elts=set()
            for p in dict_params['POTCAR']:
                if p.element not in elts:
                    sum_elec += comp.as_dict()[p.element]*p.nelectrons
                    elts.add(p.element)
            if charge != 0:
                incar['NELECT'] = sum_elec-charge

            kpoint=dict_params['KPOINTS'].monkhorst_automatic()		#I think we need better than 222 grid. should use automatic_density with 1000 /atom?

            path = os.path.join(
                    path_base, defect['name'], "charge_"+str(charge))
            try:
                os.makedirs(path)
            except:
                pass
            if hse:
                incar.update({"LWAVE": True})
                incar.write_file(os.path.join(path,"INCAR.relax.gga"))
                incar.update({
                    'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2, #"AEXX": 0.45, 
                    "PRECFOCK": "Fast", 'NKRED': 2})
                incar.write_file(os.path.join(path,"INCAR.relax.hse1"))
                del incar['PRECFOCK']
                del incar['NKRED']
                incar.write_file(os.path.join(path,"INCAR.relax.hse2"))
            else:
                incar.write_file(os.path.join(path,"INCAR.relax"))
            kpoint.write_file(os.path.join(path, "KPOINTS"))
            dict_params['POSCAR'].write_file(os.path.join(path, "POSCAR.orig"))
            dict_params['POTCAR'].write_file(os.path.join(path, "POTCAR"))
            dumpfn(dict_transf, os.path.join(path, 'transformation.json'),
                   cls=MontyEncoder)

            # Write addition incar files for dos plots of defect levels
            del incar['NSW']
            del incar['ISIF']
            del incar['EDIFFG']
            del incar['LVHAR']
            del incar['LVTOT']
            incar['IBRION'] = -1
            incar['ICHARG'] = 1
            incar['EDIFF'] = 1e-6
            # High acc not required for chgcar
            incar.update({"PRECFOCK": "Fast", 'NKRED': 2}) 
            incar.write_file(os.path.join(path, "INCAR.static"))
            del incar['PRECFOCK']
            del incar['NKRED']
            incar['ICHARG'] = 11
            incar['NEDOS'] = int((dos_limits[1]-dos_limits[0])/0.006)
            incar['EMIN'] = dos_limits[0]
            incar['EMAX'] = dos_limits[1]
            incar.write_file(os.path.join(path, "INCAR.dos"))

    # Generate bulk supercell inputs
    s = bulk_sys
    dict_transf={'defect_type': 'bulk', 'supercell': s['size']}

    dict_params = MPVaspInputSet().get_all_vasp_input(s['structure'])
    incar = dict_params['INCAR']
    incar.update({
        'IBRION': -1, "NSW": 0, 'ISPIN': 2, 'LWAVE': False, 'EDIFF': 1e-5,
        'ISMEAR': 0, 'SIGMA': 0.05, 'LVTOT': True, 'LVHAR': True, 
        'ALGO': 'Fast', 'ISYM': 0})
    if user_settings:
        if 'INCAR' in user_settings.get('bulk', None):
            incar.update(user_settings['bulk']['INCAR'])
    #if hse == True:
    #    incar.update({
    #        'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2,
    #        "PRECFOCK": "Fast", "AEXX": 0.45})
    kpoint=dict_params['KPOINTS'].monkhorst_automatic()
    path=os.path.join(path_base, 'bulk')
    try:
        os.makedirs(path)
    except:
        pass
    if hse:
        incar.update({"LWAVE": True})
        incar.write_file(os.path.join(path,"INCAR.gga"))
        incar.update({
            'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2, #"AEXX": 0.45, 
            "PRECFOCK": "Fast", 'NKRED': 2})
        incar.write_file(os.path.join(path,"INCAR.hse1"))
        del incar['PRECFOCK']
        del incar['NKRED']
        incar.write_file(os.path.join(path,"INCAR.hse2"))
    else:
        incar.write_file(os.path.join(path,"INCAR"))
    kpoint.write_file(os.path.join(path, "KPOINTS"))
    dict_params['POSCAR'].write_file(os.path.join(path, "POSCAR"))
    dict_params['POTCAR'].write_file(os.path.join(path,"POTCAR"))
    dumpfn(dict_transf, os.path.join(path,'transformation.json'),
           cls=MontyEncoder)

def make_vasp_dielectric_files(struct, path=None, user_settings={}, 
        hse=False):
    """
    Generates VASP files for dielectric constant computations
    Args:
        struct:
            unitcell in pymatgen structure format 
        user_settings:
            Settings in dict format to override the defaults used in 
            generating vasp files. The format of the dictionary is
            {'INCAR':{...}, 'KPOINTS':{...}}
        hse:
            hse run or not
    """

    # Generate vasp inputs for dielectric constant

    dict_params = MPVaspInputSet().get_all_vasp_input(struct)
    incar = dict_params['INCAR']
    incar.update({
        'ISPIN': 1, 'LWAVE': False, 'EDIFF': 1e-5,
        'ISMEAR': -5, 'ALGO': 'Fast', 'ISIF': 2})
    incar.update({'IBRION': 8, 'LEPSILON': True, 'LPEAD': True})
    if 'INCAR' in user_settings:
        incar.update(user_settings['INCAR'])
    if 'NSW' in incar:
        del incar['NSW']
    if hse == True:
        incar.update({
            'LHFCALC': True, "ALGO": "All", "HFSCREEN": 0.2,
            "PRECFOCK": "Fast", 'NKRED': 2})#, "AEXX": 0.45})

    kpoints = dict_params['KPOINTS'].automatic_density(struct,1000,force_gamma=True)

    if not path:
        path_base = struct.composition.reduced_formula
        path = os.path.join(path_base, 'dielectric')
    os.makedirs(path)
    incar.write_file(os.path.join(path, "INCAR"))
    kpoints.write_file(os.path.join(path, "KPOINTS"))
    dict_params['POSCAR'].write_file(os.path.join(path, "POSCAR"))
    dict_params['POTCAR'].write_file(os.path.join(path, "POTCAR"))
