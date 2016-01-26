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
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer

from pycdt.corrections.defects_analyzer import ParsedDefect, DefectsAnalyzer


class PostProcess(object):
    def __init__(self, root_fldr, mpid, mapi_key=None):
        """
        Post processing object for charged point-defect calculations.

        Args:
            root_fldr (str): path (relative) to directory
                in which data of charged point-defect calculations for
                a particular system are to be found;
            mpid (str): Materials Project ID of bulk structure; 
                format "mp-X", where X is an integer;
            mapi_key (str): Materials API key to access database.
        """
        self._root_fldr = root_fldr
        self._mpid = mpid
        self._mapi_key = mapi_key
        self._substitution_species = set()

    def parse_defect_calculations(self):
        """
        Parses the defect calculations as ComputedStructureEntries/ParsedDefects.
        Charge correction is missing in the first run.
        """
        parsed_defects = []
        subfolders = glob.glob(os.path.join(self._root_fldr,"bulk"))
        subfolders += glob.glob(os.path.join(self._root_fldr,"vac_*"))
        subfolders += glob.glob(os.path.join(self._root_fldr,"as_*"))
        subfolders += glob.glob(os.path.join(self._root_fldr,"sub_*"))

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

        def get_encut_from_potcar(fldr):
            potcar_file = os.path.join(fldr,'POTCAR')
            if not os.path.exists(potcar_file):
                error_msg = ": Failure, POTCAR file."
                return (None, error_msg) #Further processing is not useful

            try:
                potcar = Potcar.from_file(potcar_file)
            except:
                error_msg = ": Failure, couldn't read POTCAR file."
                return (None, error_msg)

            encut = max(ptcr_sngl.enmax for ptcr_sngl in potcar)
            return encut, None


        for fldr in subfolders:
            fldr_name = os.path.split(fldr)[1]
            fldr_fields = fldr_name.split("_")
            if 'bulk' in fldr_fields:
                vr, error_msg = get_vr_and_check_locpot(fldr)
                if error_msg:
                    print (fldr_name, error_msg)
                    print "Abandoing parsing of the calculations"
                    break
                bulk_energy = vr.final_energy
                bulk_struct = vr.final_structure
                bulk_sites = vr.structures[-1].num_sites
                bulk_locpot_path = os.path.abspath(os.path.join(fldr,'LOCPOT'))
                bulk_entry = ComputedStructureEntry(
                        bulk_struct, bulk_energy,
                        data={'locpot_path':bulk_locpot_path})
            else:
                chrg_fldrs = glob.glob(os.path.join(fldr,'charge*'))
                print ('charge folders', chrg_fldrs)
                for chrg_fldr in chrg_fldrs:
                    trans_dict = loadfn(
                            os.path.join(chrg_fldr, 'transformation.json'), 
                            cls=MontyDecoder)
                    chrg = trans_dict['charge']
                    vr, error_msg = get_vr_and_check_locpot(chrg_fldr)
                    if error_msg:
                        print (fldr_name, 'charge- ', chrg, error_msg)
                        print "But parsing of the rest of the calculations"
                        continue  # The successful calculations maybe useful
                    if 'substitution_specie' in trans_dict:
                        self._substitution_species.add(
                                trans_dict['substitution_specie'])
                        
                    site = trans_dict['defect_supercell_site']
                    energy = vr.final_energy
                    struct = vr.final_structure
                    try: 
                        encut = vr.incar['ENCUT'] 
                    except: # ENCUT not specified in INCAR. Read from POTCAR
                        try:
                            encut = get_encut_from_potcar(chrg_fldr)
                        except:
                            print (fldr_name, 'Not able to determine ENCUT') 
                            print "But parsing of the rest of the calculations"
                            continue

                    locpot_path = os.path.abspath(
                            os.path.join(chrg_fldr, 'LOCPOT'))
                    comp_data = {'locpot_path': locpot_path, 'encut': encut}
                    if 'substitution_specie' in trans_dict:
                        comp_data['substitution_specie'] = \
                                trans_dict['substitution_specie']
                    comp_def_entry = ComputedStructureEntry(
                            struct, energy, data=comp_data)
                    parsed_defects.append(
                            ParsedDefect( 
                                comp_def_entry, site_in_bulk=site, 
                                charge=chrg, name=fldr_name))

        else:
            parsed_defects_data = {}
            parsed_defects_data['bulk_entry'] = bulk_entry
            parsed_defects_data['defects'] = parsed_defects
            return parsed_defects_data

        return {} # Return Null dict due to failure

    def get_vbm_bandgap(self):
        """
        Returns the valence band maxiumum (float) of the structure with
        MP-ID mpid.

        Args:
            mpid (str): MP-ID for which the valence band maximum is to
                be fetched from the Materials Project database
        """

        if not self._mapi_key:
            with MPRester() as mp:
                bs = mp.get_bandstructure_by_material_id(self._mpid)
        else:
            with MPRester(self._mapi_key) as mp:
                bs = mp.get_bandstructure_by_material_id(self._mpid)
        if not bs:
            raise ValueError("Could not fetch band structure!")

        vbm = bs.get_vbm()['energy']
        if not vbm:
            vbm = 0
        bandgap = bs.get_band_gap()
        return (vbm, bandgap)

    def get_chempot_limits(self, structure=None):
        """
        Returns atomic chempots from mpid
        """
        if not structure:
            if not self._mapi_key:
                with MPRester() as mp:
                    structure = mp.get_structure_by_material_id(self._mpid)
            else:
                with MPRester(self._mapi_key) as mp:
                    structure = mp.get_structure_by_material_id(self._mpid)
            if  not structure:
                raise ValueError("Could not fetch structure for atomic chempots!")

        species = structure.types_of_specie
        species_symbol = [s.symbol for s in species]
        #print 'look for atomic chempots relative to:',species

        if not self._mapi_key:
            with MPRester() as mp:
                entries = mp.get_entries_in_chemsys(species_symbol)
        else:
            with MPRester(self._mapi_key) as mp:
                entries = mp.get_entries_in_chemsys(species_symbol)
        if  not entries:
            raise ValueError("Could not fetch entries for atomic chempots!")

        if len(species) == 1:
            print 'this is elemental system! use bulk value.'
            vals = [entry.energy_per_atom for entry in entries]
            chempot = {species[0]: min(vals)}
            return chempot

        pd = PhaseDiagram(entries)
        print pd

        chem_lims = {}
        print 'species=', species
        for specie in species:
            mu_lims = PDAnalyzer(pd).get_chempot_range_stability_phase(
                    structure.composition, specie)
            sp_symb = specie.symbol
            chem_lims[sp_symb] = {'rich': {},'poor': {}}
            for el in mu_lims.keys():
                chem_lims[sp_symb]['rich'][el.symbol] = mu_lims[el][1]
                chem_lims[sp_symb]['poor'][el.symbol] = mu_lims[el][0]

        # For substitution species the values are equal to single element vals
        for sub_el in self._substitution_species:
            print sub_el
            if not self._mapi_key:
                with MPRester() as mp:
                    entries = mp.get_entries_in_chemsys([sub_el])
            else:
                with MPRester(self._mapi_key) as mp:
                    entries = mp.get_entries_in_chemsys([sub_el])
            if  not entries:
                raise ValueError("Could not fetch entries for atomic chempots!")

            vals = [entry.energy_per_atom for entry in entries]
            for sp_symb in chem_lims:
                chem_lims[sp_symb]['rich'][sub_el] = min(vals)
                chem_lims[sp_symb]['poor'][sub_el] = min(vals)

        #make this less confusing for binary systems...
        if len(species) == 2:
            first_specie = sorted(chem_lims.keys())[0]
            for key in chem_lims.keys():
                if key is not first_specie:
                    del chem_lims[key]
            #chem_lims = chem_lims[chem_lims.keys()[0]]



        return chem_lims

    def parse_dielectric_calculation(self):
        """
        Parses the "vasprun.xml" file in subdirectory "dielectric" of 
        root directory root_fldr and returns the average of the trace
        of the dielectric tensor.

        Args:
            root_fldr (str):
                root directory where subdirectory "dielec" is expected
        Returns:
            eps (float):
                average of the trace of the dielectric tensor
        """

        try:
            vr = Vasprun(os.path.join(
                self._root_fldr,"dielectric","vasprun.xml"))
        except:
            print 'Parsing Dielectric calculation failed'
            return None

        eps_ion = vr.epsilon_ionic
        eps_stat = vr.epsilon_static

        eps = []
        for i in range(len(eps_ion)):
            eps.append([e[0]+e[1] for e in zip(eps_ion[i],eps_stat[i])])

        return eps

    def compile_all(self):
        """
        Run to get all post processing objects as dictionary

        note: still need to implement
            1) ability for substitutional atomic chempots
            2) incorporated charge corrections for defects
        """
        output = self.parse_defect_calculations()    
        output['epsilon'] = self.parse_dielectric_calculation()
        output['mu_range'] = self.get_chempot_limits()
        vbm,gap = self.get_vbm_bandgap()
        output['vbm'] = vbm
        output['gap'] = gap

        return output

