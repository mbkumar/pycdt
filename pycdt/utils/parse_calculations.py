# coding: utf-8
"""
Parses the computed data from VASP defect calculations.
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
import glob
import logging

from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.entries.computed_entries import ComputedStructureEntry

from pycdt.core.defects_analyzer import ComputedDefect 
from pycdt.core.chemical_potentials import ChemPotAnalyzer

class PostProcess(object):
    def __init__(self, root_fldr, mpid=None, mapi_key=None):
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
        #self._chem_pot_details = None

    def parse_defect_calculations(self):
        """
        Parses the defect calculations as ComputedDefects.
        Charge correction is missing in the first run.
        """
        logger = logging.getLogger(__name__)
        parsed_defects = []
        subfolders = glob.glob(os.path.join(self._root_fldr, "bulk"))
        subfolders += glob.glob(os.path.join(self._root_fldr, "vac_*"))
        subfolders += glob.glob(os.path.join(self._root_fldr, "as_*"))
        subfolders += glob.glob(os.path.join(self._root_fldr, "sub_*"))
        subfolders += glob.glob(os.path.join(self._root_fldr, "inter_*"))

        def get_vr_and_check_locpot(fldr):
            vr_file = os.path.join(fldr,'vasprun.xml')
            if not os.path.exists(vr_file):
                logger.warning("{} doesn't exit".format(vr_file))
                error_msg = ": Failure, vasprun.xml doesn't exist."
                return (None, error_msg) #Further processing is not useful

            try:
                vr = Vasprun(vr_file)
            except:
                logger.warning("Couldn't parse {}".format(vr_file))
                error_msg = ": Failure, couldn't parse vaprun.xml file."
                return (None, error_msg)

            if not vr.converged:
                logger.warning(
                    "Vasp calculation at {} not converged".format(fldr))
                error_msg = ": Failure, Vasp calculation not converged."
                return (None, error_msg) # Further processing is not useful

            # Check if locpot exists
            locpot_file = os.path.join(fldr, 'LOCPOT')
            if not os.path.exists(locpot_file):
                logger.warning("{} doesn't exit".format(locpot_file))
                error_msg = ": Failure, LOCPOT doesn't exist"
                return (None, error_msg) #Further processing is not useful

            return (vr, None)

        def get_encut_from_potcar(fldr):
            potcar_file = os.path.join(fldr,'POTCAR')
            if not os.path.exists(potcar_file):
                logger.warning("Not POTCAR in {} to parse ENCUT".format(fldr))
                error_msg = ": Failure, No POTCAR file."
                return (None, error_msg) #Further processing is not useful

            try:
                potcar = Potcar.from_file(potcar_file)
            except:
                logger.warning("Couldn't parse {}".format(potcar_file))
                error_msg = ": Failure, couldn't read POTCAR file."
                return (None, error_msg)

            encut = max(ptcr_sngl.enmax for ptcr_sngl in potcar)
            return (encut, None)


        for fldr in subfolders:
            fldr_name = os.path.split(fldr)[1]
            fldr_fields = fldr_name.split("_")
            if 'bulk' in fldr_fields:
                vr, error_msg = get_vr_and_check_locpot(fldr)
                if error_msg:
                    #print(fldr_name, error_msg)
                    logger.error("Abandoning parsing of the calculations")
                    break
                bulk_energy = vr.final_energy
                bulk_struct = vr.final_structure
                try: 
                    encut = vr.incar['ENCUT'] 
                except: # ENCUT not specified in INCAR. Read from POTCAR
                    encut, error_msg = get_encut_from_potcar(fldr)
                    if error_msg:
                        logger.error("Abandoning parsing of the calculations")
                        break
                        #raise AttributeError(error_msg)

                trans_dict = loadfn(
                            os.path.join(fldr, 'transformation.json'),
                            cls=MontyDecoder)
                supercell_size = trans_dict['supercell']

                bulk_locpot_path = os.path.abspath(os.path.join(fldr, 'LOCPOT'))
                bulk_entry = ComputedStructureEntry(
                        bulk_struct, bulk_energy,
                        data={'locpot_path': bulk_locpot_path,
                              'encut': encut,
                              'supercell_size': supercell_size})
            else:
                chrg_fldrs = glob.glob(os.path.join(fldr,'charge*'))
                for chrg_fldr in chrg_fldrs:
                    trans_dict = loadfn(
                            os.path.join(chrg_fldr, 'transformation.json'), 
                            cls=MontyDecoder)
                    chrg = trans_dict['charge']
                    supercell_size = trans_dict['supercell']
                    vr, error_msg = get_vr_and_check_locpot(chrg_fldr)
                    if error_msg:
                        logger.warning("Parsing the rest of the calculations")
                        continue
                    if 'substitution_specie' in trans_dict:
                        self._substitution_species.add(
                                trans_dict['substitution_specie'])
                    elif 'inter' in trans_dict['defect_type']:
                        #added because extrinsic interstitials don't have
                        # 'substitution_specie' character...
                        self._substitution_species.add(
                                trans_dict['defect_site'].specie.symbol)
                        
                    site = trans_dict['defect_supercell_site']
                    multiplicity = trans_dict.get('defect_multiplicity', None)
                    energy = vr.final_energy
                    struct = vr.final_structure
                    try: 
                        encut = vr.incar['ENCUT'] 
                    except: # ENCUT not specified in INCAR. Read from POTCAR
                        encut, error_msg = get_encut_from_potcar(chrg_fldr)
                        if error_msg:
                            logger.warning("Not able to determine ENCUT "
                                           "in {}".format(fldr_name))
                            logger.warning("Parsing the rest of the "
                                           "calculations")
                            continue

                    locpot_path = os.path.abspath(
                            os.path.join(chrg_fldr, 'LOCPOT'))
                    comp_data = {'locpot_path': locpot_path, 'encut': encut}
                    if 'substitution_specie' in trans_dict:
                        comp_data['substitution_specie'] = \
                                trans_dict['substitution_specie']
                    comp_def_entry = ComputedStructureEntry(
                            struct, energy, data=comp_data)
                    parsed_defects.append(ComputedDefect(
                            comp_def_entry, site_in_bulk=site,
                            multiplicity=multiplicity,
                            supercell_size=supercell_size,
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
        logger = logging.getLogger(__name__)
        if self._mpid is None:
                logger.warning(
                    'No mp-id provided, will fetch CBM/VBM details from the '
                    'bulk calculation.\nNote that it would be better to '
                    'perform real band structure calculation...')
                vr = Vasprun(os.path.join(self._root_fldr, 'bulk',
                                          'vasprun.xml'))
                bandgap = vr.eigenvalue_band_properties[0]
                vbm = vr.eigenvalue_band_properties[2]
        else:
            with MPRester(api_key=self._mapi_key) as mp:
                bs = mp.get_bandstructure_by_material_id(self._mpid)
            if not bs:
                logger.error("Could not fetch band structure!")
                raise ValueError("Could not fetch band structure!")

            vbm = bs.get_vbm()['energy']
            if not vbm:
                vbm = 0
            bandgap = bs.get_band_gap()['energy']

        return (vbm, bandgap)

    def get_chempot_limits(self, bulk_composition=None, chem_pot_details={}):
        """
        Returns atomic chempots from bulk_composition based on data in
        the materials project database. This is abstractly handled in the
        ChemPotAnalyzer so as to make beyond PBE-GGA easier to extend to

        chem_pot_details(dict): Is a dict from a previous query to MP
            database for chemical potentials It is not neccessary, but
            makes query faster because of less calls to MP database...
            typically will exist in defect_data.json
        """
        logger = logging.getLogger(__name__)
        if not bulk_composition:
            if not self._mpid:
                bulkvr = Vasprun(os.path.join(
                    self._root_fldr, "bulk", "vasprun.xml"))
                structure = bulkvr.final_structure
            else:
                with MPRester(api_key=self._mapi_key) as mp:
                    structure = mp.get_structure_by_material_id(self._mpid)
            if not structure:
                msg = "Could not fetch structure for atomic chempots!"
                logger.warning(msg)
                raise ValueError(msg)
            else:
                bulk_composition = structure.composition

        if chem_pot_details:
            cpa = ChemPotAnalyzer.from_dict(chem_pot_details)
        else:
            cpa = ChemPotAnalyzer(bulk_composition,
                                  sub_species=self._substitution_species)

        # Note that _substitution_species set is something that needs to be 
        # pre-loaded (as it is in the parse_defect_calculations attribute)
        if not self._mpid:
            bulkvr = Vasprun(os.path.join(
                self._root_fldr, "bulk", "vasprun.xml"))
            bulk_entry = bulkvr.get_computed_entry()
        else:
            bulk_entry = None # Use MP entry inside cpa

        chem_lims = cpa.analyze_GGA_chempots(bulk_entry=bulk_entry,
                                             mpid=self._mpid,
                                             mapi_key=self._mapi_key)

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
            logging.getLogger(__name__).warning(
                'Parsing Dielectric calculation failed')
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
        #entry data from MP query..could print for reducing future queries
        # output['chem_pot_details'] = self._chem_pot_details
        vbm,gap = self.get_vbm_bandgap()
        output['vbm'] = vbm
        output['gap'] = gap

        return output

