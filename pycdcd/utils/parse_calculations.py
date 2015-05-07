# coding: utf-8
"""
Parses the vasprun.xml files generated during VASP defect calculations.
"""
#from __future__ import unicode_literals
from __future__ import division

__author__ = "Bharat Medasani"
__data__  = "Sep 14, 2014"

import os
import sys
import glob
from argparse import ArgumentParser

from pymatgen.matproj.rest import MPRester
from monty.serialization import dumpfn
from monty.json import MontyEncoder
from pymatgen.io.vaspio.vasp_output import Vasprun


def parse_defect_energy(structure, root_fldr):

    energy_dict = {}

    antisites = []
    vacancies = []
    substitutions = []
    subfolders += glob.glob(os.path.join(root_fldr,"bulk"))
    subfolders = glob.glob(os.path.join(root_fldr,"vac_*"))
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
            bulk_sites = vr.structures[-1].num_sites
            bulk_locpot_path = os.path.abspath(os.path.join(
                fldr_name,'LOCPOT'))
        else:
            chrg_fldrs = glob.glob(os.path.join(fldr,'charge_*'))
            energies = {}
            locpot_paths = {}
            for chrg_fldr in chrg_fldrs:
                chrg = int(chrg_fldr_name.split("_")[1])
                vr, error_msg = get_vr_and_check_locpot(chrg_fldr)
                if error_msg:
                    print (fldr_name, 'charge- ', chrg, error_msg)
                    print "But parsing of the rest of the calculations"
                    continue  # The successful calculations maybe useful
                energies[chrg] = vr.final_energy
                locpot_paths[chrg] = os.path.abspath(os.path.join(
                    chrg_fldr, 'LOCPOT')
            if 'vac' in fldr_fields:
                site_index = int(fldr_fields[1])
                site_multiplicity = int(fldr_fields[2].split("-")[1])
                site_specie = fldr_fields[3].split("-")[1]
                vacancies.append({'site_index':site_index,
                    'site_specie':site_specie,'energies':energies,
                    'site_multiplicity':site_multiplicity,
                    'locpot_paths':locpot_paths
                    })

            else:
                site_index = int(fldr_fields[1])
                site_multiplicity = int(fldr_fields[2].split("-")[1])
                site_specie = fldr_fields[3].split("-")[1]
                substitution_specie = fldr_fields[4].split("-")[1]
                dict_def = {'site_index':site_index,
                        'site_specie':site_specie,'energies':energies,
                        'substitution_specie':substitution_specie,
                        'site_multiplicity':site_multiplicity,
                        'locpot_paths':locpot_paths
                        }
                if 'as' in fldr_fields:
                    antisites.append(dict_def)
                else:
                    substitutions.append(dict_def)
    else:
        print "All calculations successful for ", mpid
        e0 = bulk_energy/bulk_sites*structure.num_sites
        for vac in vacancies:
            for charge, energy in vac['energies'].items():
                flip_energy = energy - bulk_energy
                vac['energies'][charge] = flip_energy
        vacancies.sort(key=lambda entry: entry['site_index'])
        for antisite in antisites:
            for charge, energy in antisite['energies'].items():
                flip_energy = energy - bulk_energy
                antisite['energies'][charge] = flip_energy
        antisites.sort(key=lambda entry: entry['site_index'])
        for sub in substitutions:
            for charge, energy in sub['energies'].items():
                flip_energy = energy - bulk_energy
                sub['energies'][charge] = flip_energy
        substitutions.sort(key=lambda entry: entry['site_index'])
        energy_dict[unicode(mpid)] = {u"structure":structure,
                'e0':e0,'vacancies':vacancies,'antisites':antisites,
                'substitutions':substitutions}
        return energy_dict

    return {} # Return Null dict due to failure

