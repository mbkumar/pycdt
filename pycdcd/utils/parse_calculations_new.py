# coding: utf-8
"""
Parses the vasprun.xml files generated during VASP defect calculations
in conformation with DefectsAnalyzer written by Geoffroy.
"""
#from __future__ import unicode_literals
from __future__ import division

__author__ = "Bharat Medasani, Nils Zimmermann, Danny Broberg"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = 'mbkumar@gmail.com'
__date__  = "Sep 14, 2014"

import os
import sys
import glob
from argparse import ArgumentParser

from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.matproj.rest import MPRester
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pycdcd.corrections.defects_analyzer import ParsedDefect, DefectsAnalyzer


def parse_defect_calculations(root_fldr):
    """
    Parses the defect calculations as ComputedStructureEntries/ParsedDefects.
    Charge correction is missing in the first run.
    """
    parsed_defects = []
    subfolders = glob.glob(os.path.join(root_fldr,"bulk"))
    subfolders += glob.glob(os.path.join(root_fldr,"vac_*"))
    subfolders += glob.glob(os.path.join(root_fldr,"as_*"))
    subfolders += glob.glob(os.path.join(root_fldr,"sub_*"))

    def get_vr_and_check_locpot(fldr):
        vr_file = os.path.join(fldr,'vasprun.xml') 
        if not os.path.exists(vr_file):
            error_msg = ": Failure, vasprun.xml doesn't exist."
            return (None, error_msg) #Further processing is not useful

        try:
            vr = Vasprun(vr_file)
        except:
            error_msg = ": Failure, couldn't parse vaprun.xml file."
            return (None, error_msg)

        if not vr.converged:
            error_msg = ": Failure, Vasp calculation not converged."
            return (None, error_msg) # Further processing is not useful

        # Check if locpot exists
        locpot_file = os.path.join(fldr, 'LOCPOT')
        if not os.path.exists(locpot_file):
            error_msg = ": Failure, LOCPOT doesn't exist" 
            return (None, error_msg) #Further processing is not useful

        return (vr, None)

    for fldr in subfolders:
        fldr_name = os.path.split(fldr)[1]
        fldr_fields = fldr_name.split("_")
        if 'bulk' in fldr_fields:
            vr, error_msg = get_vr_and_check_locpot(fldr)
            if errro_msg:
                print (fldr_name, error_msg)
                print "Abandoing parsing of the calculations"
                break
            bulk_energy = vr.final_energy
            bulk_struct = vr.final_structure
            bulk_sites = vr.structures[-1].num_sites
            bulk_locpot_path = os.path.abspath(os.path.join(fldr,'LOCPOT'))
            bulk_entry = ComputedStructureEntry(bulk_struct, bulk_energy, 
                    data={'locpot_path':bulk_locpot_path})
        else:
            chrg_fldrs = glob.glob(os.path.join(fldr,'charge_*'))
            for chrg_fldr in chrg_fldrs:
                vr, error_msg = get_vr_and_check_locpot(chrg_fldr)
                if error_msg:
                    print (fldr_name, 'charge- ', chrg, error_msg)
                    print "But parsing of the rest of the calculations"
                    continue  # The successful calculations maybe useful

                #chrg = int(chrg_fldr_name.split("_")[1])
                trans_dict = loadfn('transformation.json')
                chrg = trans_dict['charge']
                site = trans_dict['defect_supercell_site']
                energy = vr.final_energy
                struct = vr.final_structure
                encut = vr.incar['ENCUT']
                locpot_path = os.path.abspath(os.path.join(chrg_fldr, 'LOCPOT'))
                parsed_defects.append(ParsedDefect(
                    ComputedStructureEntry(struct, energy, 
                        data={'locpot_path':locpot_path,'encut':encut}),
                    site_in_bulk=site, charge=chrg, name=fldr_name))

    else: 
        parsed_defects_data = {}
        parsed_defects_data['bulk_entry'] = bulk_entry 
        parsed_defects_data['defects'] = parsed_defects 
        return parsed_defects_data

    return {} # Return Null dict due to failure


def get_vbm(mpid, mapi_key=None):
    """
    Returns the valence band maxiumum (float) of the structure with
    MP-ID mpid.

    Args:
        mpid (str): MP-ID for which the valence band maximum is to
            be fetched from the Materials Project database
        mapi_key: Materials API key to access database
    """

    if not mapi_key:
        with MPRester() as mp:
            bs = mp.get_bandstructure_by_material_id(mpid)
    else:
        with MPRester(mapi_key) as mp:
            bs = mp.get_bandstructure_by_material_id(mpid)
    if  not bs:
        raise ValueError("Could not fetch band structure!")

    vbm_energy = bs.get_vbm()['energy']
    if not vbm_energy:
        vbm_energy = 0


def get_atomic_chempots(structure):
    """
    gets atomic chempots from MP database

    note: could also do this with mpid if that would be easier..
    """
    specs=list(set(structure.species))
    listspec=[i.symbol for i in specs]
    print 'look for atomic chempots relative to:',listspec
    
    mp=MPRester() 
    entries=mp.get_entries_in_chemsys(listspec)
    if len(listspec)==1:
        print 'this is elemental system! use bulk value.'
        vals=[i.energy_per_atom for i in entries]
        chempot={specs[0]:min(vals)}
        return chempot
    pd=PhaseDiagram(entries)
    print pd

    chemlims={}
    for i in specs:
        name=str(i)+'-limiting'
        tmpchemlims=PDAnalyzer(pd).get_chempot_range_stability_phase(self._structure.composition,i)
        chemlims[name]={str(i)+'-rich':{},str(i)+'-poor':{}}
        for j in tmpchemlims.keys():
            chemlims[name][str(i)+'-rich'][j]=tmpchemlims[j][1]
            chemlims[name][str(i)+'-poor'][j]=tmpchemlims[j][0]

    #make this less confusing for binary systems...
    if len(specs)==2:
	chemlims=chemlims[chemlims.keys()[0]]

    return chemlims 


def parse_dielectric_calculation(root_fldr):
    """
    Parses the "vasprun.xml" file in subdirectory "dielec" of root
    directory root_fldr and returns the dielectric tensor.

    Args:
        root_fldr (str):
            root directory where subdirectory "dielec" is expected
    Returns:
        eps (3x3 float matrix):
            dielectric tensor
    """

    vrun = Vasprun(os.path.join(root_fldr,"dielec","vasprun.xml"))
    eps_ion = vrun.epsilon_ionic
    eps_stat = vrun.epsilon_static

    eps = []
    for i in range(len(eps_ion)):
        eps.append([e[0]+e[1] for e in zip(eps_ion[i],eps_stat[i])])
    
    return eps

