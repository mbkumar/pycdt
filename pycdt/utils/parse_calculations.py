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
import numpy as np

from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder

from pymatgen.core import PeriodicSite
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Potcar
# from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.analysis.defects.core import Vacancy, Substitution, Interstitial, DefectEntry

# from pycdt.core.defects_analyzer import ComputedDefect
from pycdt.core.chemical_potentials import MPChemPotAnalyzer


# below here is for single defect parser
from atomate.vasp.drones import VaspDrone

from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility


class SingleDefectParser(object):

    def __init__(self, path_to_defect, path_to_bulk, dielectric,
                 defect_charge = None, mpid = None,
                 compatibility=DefectCompatibility()):
        """
        Parse a defect object using many of the features of a standard Defect Builder.

        :param path_to_defect (str): path to defect file of interest
        :param path_to_bulk (str): path to bulk file of interest
        :param dielectric (float or 3x3 matrix): ionic + static contributions to dielectric constant
        :param defect_charge (int): Can specify defect charge. If not specified, with attempt to find
            charge based on standard Potcar settings, but this can sometimes lead to failure.
        """
        self.path_to_defect = path_to_defect
        self.path_to_bulk = path_to_bulk
        self.dielectric = dielectric
        self.parameters = {}
        self.defect_site = None
        self.defect_type = None
        self.defect_charge = defect_charge
        self.mpid = mpid
        self.compatibility = compatibility

    def get_defect_entry( self, run_compatibility=False):
        """
        :return: DefectEntry object with all possible defect meta data loaded
        """
        bulk_drone = VaspDrone( bandstructure_mode="auto", defect_wf_parsing=None)
        print("Parsing bulk file path")
        bulk_dict = bulk_drone.assimilate( self.path_to_bulk)
        print("Starting defect parsing...")
        defect_dict = bulk_drone.assimilate( self.path_to_defect) #just for grabbing initial structure
        print("\tInitial assimilation complete. Now run parser for "
              "defect-relevant metadata...")

        # add bulk simple properties
        bulk_energy = bulk_dict['output']['energy']
        bulk_sc_structure = Structure.from_dict( bulk_dict['input']['structure'])
        self.parameters.update({"dielectric": self.dielectric,
            'bulk_energy': bulk_energy, 'bulk_sc_structure': bulk_sc_structure})

        # add bulk data required for corrections and localization analysis
        bulklpt = bulk_dict['calcs_reversed'][0]['output']['locpot']
        axes = list(bulklpt.keys())
        axes.sort()
        self.parameters.update({'bulk_planar_averages': [bulklpt[ax] for ax in axes]})

        bulkoutcar = bulk_dict['calcs_reversed'][0]['output']['outcar']
        bulk_atomic_site_averages = bulkoutcar['electrostatic_potential']
        self.parameters.update({'bulk_atomic_site_averages': bulk_atomic_site_averages})

        # get functional / INCAR metadata
        self.get_stdrd_metadata( bulk_dict.copy(), defect_dict.copy())

        # identify defect site, structural information, and charge
        self.identify_defect_info( defect_dict.copy())

        # create initial defect object
        for_monty_defect = {"@module": "pymatgen.analysis.defects.core",
                            "@class": self.defect_type,
                            "charge": self.defect_charge,
                            "structure": self.parameters["bulk_sc_structure"].copy(),
                            "defect_site": self.defect_site}
        defect = MontyDecoder().process_decoded( for_monty_defect)
        test_defect_structure = defect.generate_defect_structure()
        if not StructureMatcher(primitive_cell=False, scale=False, attempt_supercell=False,
                                allow_subset=False).fit( test_defect_structure,
                                                         self.parameters['initial_defect_structure']):
            #NOTE: this does not insure that cartesian coordinates or indexing are identical
            raise ValueError("Error in defect matching!")

        # rerun parser with site information (initialized by identify_defect)
        defect_drone = VaspDrone( bandstructure_mode="auto", defect_wf_parsing=self.defect_site)
        defect_dict = defect_drone.assimilate( self.path_to_defect)

        # grab all additional metadata required for
        self.get_defect_metadata( defect_dict.copy())

        # get bulk bulk metadata
        self.get_bulk_gap_data( bulk_dict.copy())

        defect_entry = DefectEntry( defect, self.parameters['defect_energy'] - self.parameters['bulk_energy'],
                                    corrections = {}, parameters = self.parameters)

        if run_compatibility:
            print("Running compatibility...")
            defect_entry = self.compatibility.process_entry( defect_entry)

        return defect_entry

    def get_stdrd_metadata(self, bulk_dict, defect_dict):
        run_metadata = {}
        dincar = defect_dict["calcs_reversed"][0]["input"]["incar"]
        bincar = bulk_dict["calcs_reversed"][0]["input"]["incar"]
        d_kpoints = defect_dict['calcs_reversed'][0]['input']['kpoints']
        b_kpoints = bulk_dict['calcs_reversed'][0]['input']['kpoints']
        run_metadata.update( {"defect_incar": dincar, "bulk_incar": bincar, "defect_kpoints": d_kpoints,
                              "bulk_kpoints": b_kpoints})
        run_metadata.update( {"incar_calctype_summary":
                                  {k: dincar.get(k, None) if dincar.get(k) not in ['None', 'False', False] else None
                                   for k in ["LHFCALC", "HFSCREEN", "IVDW", "LUSE_VDW", "LDAU", "METAGGA"]
                                   }
                              }
                             )
        run_metadata.update({"potcar_summary":
                                 {'pot_spec': [potelt["titel"] for potelt in defect_dict['input']['potcar_spec']],
                                  'pot_labels': defect_dict['input']['pseudo_potential']['labels'][:],
                                  'pot_type': defect_dict['input']['pseudo_potential']['pot_type'],
                                  'functional': defect_dict['input']['pseudo_potential']['functional']}
                             })
        self.parameters.update({"run_metadata": run_metadata})

        return

    def identify_defect_info(self, defect_dict):
        bulk_sc_structure = self.parameters['bulk_sc_structure'].copy()

        # get defect structures (before and after)
        final_defect_structure = defect_dict['output']['structure']
        if type(final_defect_structure) != Structure:
            final_defect_structure = Structure.from_dict(final_defect_structure)

        # build initial_defect_structure from very first ionic relaxation step (ensures they are indexed the same]
        initial_defect_structure = Structure.from_dict(
            defect_dict['calcs_reversed'][-1]['output']['ionic_steps'][0]['structure'])
        if type(initial_defect_structure) != Structure:
            initial_defect_structure = Structure.from_dict(initial_defect_structure)

        self.parameters.update({'final_defect_structure': final_defect_structure,
                                'initial_defect_structure': initial_defect_structure})

        # get defect site and object information
        num_ids = len(initial_defect_structure)
        num_bulk = len(bulk_sc_structure)
        if num_ids == num_bulk - 1:
            self.defect_type = "Vacancy"
        elif num_ids == num_bulk + 1:
            self.defect_type = "Interstitial"
        elif num_ids == num_bulk:
            self.defect_type = "Substitution"
        else:
            raise ValueError("Could not identify defec type just from number of sites in structure: "
                             "{} in bulk vs. {} in defect?".format( num_ids, num_bulk ))

        if self.defect_charge is None:
            #TODO: run guessing routine for NELECT using POTCAR information in parameters metadata
            raise ValueError("Have not implemented defect charge routine yet")

        bulksites = [site.frac_coords for site in bulk_sc_structure]
        initsites = [site.frac_coords for site in initial_defect_structure]
        distmatrix = initial_defect_structure.lattice.get_all_distances(bulksites,
                                                                        initsites)
        min_dist_with_index = [[min(distmatrix[bulk_index]), int(bulk_index),
                                int(distmatrix[bulk_index].argmin())] for bulk_index in
                               range(len(distmatrix))]  # list of [min dist, bulk ind, defect ind]

        site_matching_indices = []
        poss_defect = []
        if self.defect_type in ["Vacancy", "Interstitial"]:
            for mindist, bulk_index, defect_index in min_dist_with_index:
                if mindist < 0.1:
                    site_matching_indices.append([bulk_index, defect_index])
                elif self.defect_type == "Vacancy":
                    poss_defect.append([bulk_index, bulksites[bulk_index][:]])

            if self.defect_type == "Interstitial":
                poss_defect = [[ind, fc[:]] for ind, fc in enumerate(initsites) \
                               if ind not in np.array(site_matching_indices)[:, 1]]

        elif self.defect_type == "Substitution":
            for mindist, bulk_index, defect_index in min_dist_with_index:
                species_match = bulk_sc_structure[bulk_index].specie == \
                                initial_defect_structure[defect_index].specie
                if mindist < 0.1 and species_match:
                    site_matching_indices.append([bulk_index, defect_index])

                elif not species_match:
                    poss_defect.append([defect_index, initsites[defect_index][:]])

        if len(poss_defect) == 1:
            defect_index_sc_coords = poss_defect[0][0]
            defect_frac_sc_coords = poss_defect[0][1]
        else:
            raise ValueError("Found {} possible defect sites when matching bulk and "
                             "defect structure".format(len(poss_defect)))

        if len(set(np.array(site_matching_indices)[:, 0])) != len(set(np.array(site_matching_indices)[:, 1])):
            raise ValueError("Error occured in site_matching routine. Double counting of site matching "
                              "occured:{}\nAbandoning structure parsing.".format(site_matching_indices))

        if self.defect_type == "Vacancy":
            self.defect_site = bulk_sc_structure[ defect_index_sc_coords]
        else:
            self.defect_site = initial_defect_structure[ defect_index_sc_coords]


        self.parameters.update( {'site_matching_indices': site_matching_indices,
                                 'defect_frac_sc_coords': defect_frac_sc_coords,
                                 'defect_index_sc_coords': defect_index_sc_coords})

        return

    def get_defect_metadata(self, defect_dict):

        defect_energy = defect_dict['output']['energy']
        self.parameters.update({'defect_energy': defect_energy})

        deflpt = defect_dict['calcs_reversed'][0]['output']['locpot']
        axes = list(deflpt.keys())
        axes.sort()
        defect_planar_averages = [deflpt[ax] for ax in axes]
        abc = self.parameters['initial_defect_structure'].lattice.abc
        axis_grid = []
        for ax in range(3):
            num_pts = len(defect_planar_averages[ax])
            axis_grid.append([i / num_pts * abc[ax] for i in range(num_pts)])
        self.parameters.update({'axis_grid': axis_grid,
                                'defect_planar_averages': defect_planar_averages})

        # grab Kumagai metadata (note that site matching was pulled during identify_defect_info)
        defoutcar = defect_dict['calcs_reversed'][0]['output']['outcar']
        defect_atomic_site_averages = defoutcar['electrostatic_potential']

        wz = self.parameters["initial_defect_structure"].lattice.get_wigner_seitz_cell()
        dist = []
        for facet in wz:
            midpt = np.mean(np.array(facet), axis=0)
            dist.append(np.linalg.norm(midpt))
        sampling_radius = min(dist)

        self.parameters.update( {'defect_atomic_site_averages': defect_atomic_site_averages,
                                 'sampling_radius': sampling_radius}
                                )

        # grab eigenvalue information for Band filling and localization analysis
        eigenvalues = defect_dict['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['eigenvalues']
        kpoint_weights = defect_dict['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['kpoint_weights']
        self.parameters.update({'eigenvalues': eigenvalues,
                                'kpoint_weights': kpoint_weights})

        # grab ks_delocalization information
        self.parameters.update({'defect_ks_delocal_data':
                                    defect_dict['calcs_reversed'][0]['output']['defect'].copy()})

        return

    def get_bulk_gap_data(self, bulk_dict):
        bulk_structure = self.parameters['bulk_sc_structure']

        if not self.mpid:
            try:
                with MPRester() as mp:
                    tmp_mplist = mp.get_entries_in_chemsys(list(bulk_structure.symbol_set))
                mplist = [ment.entry_id for ment in tmp_mplist if ment.composition.reduced_composition == \
                          bulk_structure.composition.reduced_composition]
            except:
                raise ValueError("Error with querying MPRester for {}".format( bulk_structure.composition.reduced_formula))

            mpid_fit_list = []
            for trial_mpid in mplist:
                with MPRester() as mp:
                    mpstruct = mp.get_structure_by_material_id(trial_mpid)
                if StructureMatcher(primitive_cell=True, scale=False, attempt_supercell=True,
                                    allow_subset=False).fit(bulk_structure, mpstruct):
                    mpid_fit_list.append( trial_mpid)

            if len(mpid_fit_list) == 1:
                self.mpid = mpid_fit_list[0]
                print("Single mp-id found for bulk structure:{}.".format( self.mpid))
            elif len(mpid_fit_list) > 1:
                num_mpid_list = [int(mp.split('-')[1]) for mp in mpid_fit_list]
                num_mpid_list.sort()
                self.mpid  = 'mp-'+str(num_mpid_list[0])
                print("Multiple mp-ids found for bulk structure:{}\nWill use lowest number mpid "
                      "for bulk band structure = {}.".format(str(mpid_fit_list), self.mpid))
            else:
                print("Could not find bulk structure in MP database after tying the "
                                  "following list:\n{}".format( mplist))
                self.mpid = None
        else:
            print("Manually fed mpid = {}".format( self.mpid))

        vbm, cbm, bandgap = None, None, None
        if 'task_level_metadata' not in self.parameters.keys():
            self.parameters['task_level_metadata'] = {}

        if self.mpid is not None:
            #TODO: NEED to be smarter about use of +U or HSE etc in MP gga band structure calculations...
            with MPRester(api_key=self._mapi_key) as mp:
                bs = mp.get_bandstructure_by_material_id(self._mpid)
            if bs:
                cbm = bs.get_cbm()['energy']
                vbm = bs.get_vbm()['energy']
                bandgap = bs.get_band_gap()['energy']
                self.parameters['task_level_metadata'].update({'MP_gga_BScalc_data':
                                                                   bs.get_band_gap().copy()})  # contains gap kpt transition

        if vbm is None or bandgap is None or cbm is None:
            if self.mpid:
                print('WARNING: Mpid {} was provided, but no bandstructure entry currently exists for it. '
                               'Reverting to use of bulk calculation.'.format( self._mpid))
            else:
                print( 'WARNING: No mp-id provided, will fetch CBM/VBM details from the '
                       'bulk calculation.')
            print('Note that it would be better to '
                  'perform real band structure calculation...')

            self.parameters['task_level_metadata'].update( {'MP_gga_BScalc_data': None}) #to signal no MP BS is used
            cbm = bulk_dict['output']['cbm']
            vbm = bulk_dict['output']['vbm']
            bandgap = bulk_dict['output']['bandgap']

        self.parameters.update( {'mpid': self.mpid, "cbm": cbm, "vbm": vbm, "gap": bandgap} )

        return


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

    def parse_defect_calculations(self):
        """
        Parses the defect calculations as ComputedDefects.
        Charge correction is missing in the first run.
        """
        logger = logging.getLogger(__name__)
        parsed_defects = []
        subfolders = glob.glob(os.path.join(self._root_fldr, "vac_*"))
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
                vr = Vasprun(vr_file, parse_potcar_file=False)
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

        # get bulk entry information first
        fldr = os.path.join(self._root_fldr, "bulk")
        vr, error_msg = get_vr_and_check_locpot(fldr)
        if error_msg:
            logger.error("Abandoning parsing of the calculations")
            return {}
        bulk_energy = vr.final_energy
        # bulk_struct = vr.final_structure
        try:
            encut = vr.incar['ENCUT']
        except:  # ENCUT not specified in INCAR. Read from POTCAR
            encut, error_msg = get_encut_from_potcar(fldr)
            if error_msg:
                logger.error("Abandoning parsing of the calculations")
                return {}

        trans_dict = loadfn(
            os.path.join(fldr, 'transformation.json'),
            cls=MontyDecoder)
        supercell_size = trans_dict['supercell']

        bulk_file_path = fldr
        # bulk_entry = ComputedStructureEntry(
        #     bulk_struct, bulk_energy,
        #     data={'bulk_path': fldr,
        #           'encut': encut,
        #           'supercell_size': supercell_size})

        # get defect entry information
        for fldr in subfolders:
            fldr_name = os.path.split(fldr)[1]
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
                defect_type = trans_dict.get('defect_type', None)
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

                comp_data = {'bulk_path': bulk_file_path,
                             'defect_path': chrg_fldr, 'encut': encut,
                             'fldr_name': fldr_name, 'supercell_size': supercell_size}
                if 'substitution_specie' in trans_dict:
                    comp_data['substitution_specie'] = \
                            trans_dict['substitution_specie']

                defect_dict = {'structure': struct, 'charge': chrg,
                               '@module': 'pymatgen.analysis.defects.core'
                               }
                defect_site = site
                if 'vac_' in defect_type:
                    defect_dict['@class'] = 'Vacancy'
                elif 'as_' in defect_type or 'sub_' in defect_type:
                    defect_dict['@class'] = 'Substitution'
                    substitution_specie = trans_dict['substitution_specie']
                    defect_site = PeriodicSite( substitution_specie, defect_site.frac_coords,
                                                defect_site.lattice, coords_are_cartesian=False)
                elif 'int_' in defect_type:
                    defect_dict['@class'] = 'Interstitial'
                else:
                    raise ValueError("defect type {} not recognized...".format(defect_type))

                defect_dict.update( {'defect_site': defect_site})
                defect = MontyDecoder().process_decoded( defect_dict)
                parsed_defects.append( DefectEntry( defect, energy - bulk_energy,
                                                    parameters=comp_data))

                #TODO: confirm not breaking anything by moving to DefectEntry approach above
                # comp_def_entry = ComputedStructureEntry(
                #         struct, energy, data=comp_data)
                # parsed_defects.append(ComputedDefect(
                #         comp_def_entry, site_in_bulk=site,
                #         multiplicity=multiplicity,
                #         supercell_size=supercell_size,
                #         charge=chrg, name=fldr_name))

        try:
            parsed_defects_data = {}
            # parsed_defects_data['bulk_entry'] = bulk_entry
            parsed_defects_data['defects'] = parsed_defects
            return parsed_defects_data
        except:
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
        vbm, bandgap = None, None

        if self._mpid is not None:
            with MPRester(api_key=self._mapi_key) as mp:
                bs = mp.get_bandstructure_by_material_id(self._mpid)
            if bs:
                vbm = bs.get_vbm()['energy']
                bandgap = bs.get_band_gap()['energy']

        if vbm is None or bandgap is None:
            if self._mpid:
                logger.warning('Mpid {} was provided, but no bandstructure entry currently exists for it. '
                               'Reverting to use of bulk calculation.'.format( self._mpid))
            else:
                logger.warning(
                    'No mp-id provided, will fetch CBM/VBM details from the '
                    'bulk calculation.')
            logger.warning('Note that it would be better to '
                           'perform real band structure calculation...')
            vr = Vasprun(os.path.join(self._root_fldr, 'bulk',
                                      'vasprun.xml'), parse_potcar_file=False)
            bandgap = vr.eigenvalue_band_properties[0]
            vbm = vr.eigenvalue_band_properties[2]

        return (vbm, bandgap)

    def get_chempot_limits(self):
        """
        Returns atomic chempots from bulk_composition based on data in
        the materials project database. This is abstractly handled in the
        ChemPotAnalyzer

        Note to user: If personal phase diagram desired,
            option exists in the pycdt.core.chemical_potentials to setup,
            run and parse personal phase diagrams for purposes of chemical potentials
        """
        logger = logging.getLogger(__name__)

        if self._mpid:
            cpa = MPChemPotAnalyzer( mpid = self._mpid, sub_species = self._substitution_species,
                                   mapi_key = self._mapi_key)
        else:
            bulkvr = Vasprun(os.path.join( self._root_fldr, "bulk", "vasprun.xml"), parse_potcar_file=False)
            if not bulkvr:
                msg = "Could not fetch computed entry for atomic chempots!"
                logger.warning(msg)
                raise ValueError(msg)
            cpa = MPChemPotAnalyzer( bulk_ce=bulkvr.get_computed_entry(),
                                   sub_species = self._substitution_species,
                                   mapi_key = self._mapi_key)

        chem_lims = cpa.analyze_GGA_chempots()

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
                self._root_fldr,"dielectric","vasprun.xml"), parse_potcar_file=False)
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
        vbm,gap = self.get_vbm_bandgap()
        output['vbm'] = vbm
        output['gap'] = gap

        return output

