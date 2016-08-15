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
import glob
import logging

from monty.serialization import loadfn, dumpfn
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer

from pycdt.core.defects_analyzer import ComputedDefect 


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
        Parses the defect calculations as ComputedStructureEntries/ComputedDefects.
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
                logger.warning("Vasp calculation at {} not converged".format(
                    fldr))
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

                bulk_locpot_path = os.path.abspath(os.path.join(fldr,'LOCPOT'))
                bulk_entry = ComputedStructureEntry(
                        bulk_struct, bulk_energy,
                        data={'locpot_path': bulk_locpot_path, 'encut': encut,
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
                        #print(fldr_name, 'charge- ', chrg, error_msg)
                        logger.warning("Parsing the rest of the calculations")
                        #print("But parsing of the rest of the calculations")
                        continue
                    if 'substitution_specie' in trans_dict:
                        self._substitution_species.add(
                                trans_dict['substitution_specie'])
                        
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
                            #print(fldr_name, 'Not able to determine ENCUT')
                            #print(error_msg)
                            #print("But parsing the rest of the calculations")
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
                            ComputedDefect( 
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
            if not self._mapi_key:
                with MPRester() as mp:
                    bs = mp.get_bandstructure_by_material_id(self._mpid)
            else:
                with MPRester(self._mapi_key) as mp:
                    bs = mp.get_bandstructure_by_material_id(self._mpid)
            if not bs:
                logger.error("Could not fetch band structure!")
                raise ValueError("Could not fetch band structure!")

            vbm = bs.get_vbm()['energy']
            if not vbm:
                vbm = 0
            bandgap = bs.get_band_gap()['energy']

        return (vbm, bandgap)

    def get_chempot_limits(self, structure=None):
        """
        TODO 1: get rid of MPRester pulling of structure. When you run this your structure should already be known
            It is an unneccssary complication for people who have structures that are not in the MP data base...
        TODO 2: allow for dependent elements to be used (i.e. PDA.get_chempot_range_stability_phase(target_comp, open_elt))
            (1 and 2 are related - since we want to allow for an option where we have the dependent chem pot limits)
        TODO 3: (when no mp-id present but composition is found) do structure check to see if input structure is
             identical to structure with identical composition in MP database

        Returns atomic chempots from mpid or structure input

        accounts for all different defect phases
        """
        logger = logging.getLogger(__name__)
        if not structure:
            if not self._mpid:
                    bulkvr = Vasprun(os.path.join(self._root_fldr,"bulk","vasprun.xml"))
                    structure = bulkvr.final_structure
            elif not self._mapi_key:
                with MPRester() as mp:
                    structure = mp.get_structure_by_material_id(self._mpid)
            else:
                with MPRester(self._mapi_key) as mp:
                    structure = mp.get_structure_by_material_id(self._mpid)
            if  not structure:
                msg = "Could not fetch structure for atomic chempots!"
                logger.warning(msg)
                raise ValueError(msg)

        bulk_species = structure.types_of_specie
        bulk_species_symbol = [s.symbol for s in bulk_species]
        bulk_composition = structure.composition

        def get_chempots_from_entries(structure, list_spec_symbol, comp, exceptions=[]):
            """
            outline for how to retrieve atomic chempots from Materials Project (MP) entries in a phase diagram (PD) object:
              1) check stability of computed entry w.r.t phase diagram generated from MP
              2) If stable then: a) if mp-id given and is in the stable entry list, proceed normally
                                 b) if mp-id not given, print 'congrats' and tell user about manual submission page on MP
                                    website, then manually insert the new object into the PD and generate chempots
                                        (note at this point there is no gurantee that the structure is actually unique,
                                          could be computational error making energy slightly lower than MP values)
              3) If not-stable then: a) check to see if composition exists among the structures in stable list of PD
                                     b) if a stable and identical composition exists in PD then
                                            print warning but continue as if it was stable (chem pots will depend on the stable phase)
                                     c) if no stable and identical composition exists in PD, then
                                            print warning and find facets that the composition is included in
            """
            if not self._mapi_key:
                with MPRester() as mp:
                    entries = mp.get_entries_in_chemsys(list_spec_symbol)
            else:
                with MPRester(self._mapi_key) as mp:
                    entries = mp.get_entries_in_chemsys(list_spec_symbol)
            if  not entries:
                msg = "Could not fetch entries for atomic chempots!"
                logger.warning(msg)
                raise ValueError(msg)

            chem_lims = {}
            pd = PhaseDiagram(entries)
            PDA = PDAnalyzer(pd)
            full_idlist = [i.entry_id for i in pd.qhull_entries]
            stable_idlist = [i.entry_id for i in pd.stable_entries]

            if self._mpid:
                if (self._mpid in full_idlist) and (self._mpid in stable_idlist):
                   logger.debug("Verified that mp-id is stable within "
                                "Materials Project {} phase diagram".format(
                                    '-'.join(list_spec_symbol)))
                   common_approach = True
                elif (self._mpid in full_idlist) and not (self._mpid in stable_idlist):
                    redcomp = comp.reduced_composition
                    common_approach = False
                    for i in pd.stable_entries:
                        if i.composition.reduced_composition==redcomp:
                            logger.warning(
                                "Input mp-id {} is unstable. Stable "
                                "composition mp-id found to be {}".format(
                                    self._mpid, i.entry_id))
                            logger.warning(
                                "Proceeding with atomic chemical potentials "
                                "with respect to stable phase.")
                            common_approach = True
                    if not common_approach:
                        logger.warning(
                            "Input mp-id {} is unstable. No stable structure "
                            "with same composition exists".format(self._mpid))
                        logger.warning(
                            "Proceeding with atomic chemical potentials "
                            "according to composition position within phase "
                            "diagram.")
                else:
                    logger.warning(
                        "Specified mp-id {} could not find it in MP phase "
                        "diagram. Reverting to assumption that mp-id is not "
                        "known.".format(self._mpid))
                    self._mpid = None
                    bulkvr = Vasprun(os.path.join(self._root_fldr, "bulk",
                                                  "vasprun.xml"))

            if not self._mpid:
                try:
                    ce = bulkvr.get_computed_entry()
                except:
                    bulkvr = Vasprun(os.path.join(self._root_fldr, "bulk",
                                                  "vasprun.xml"))
                    ce = bulkvr.get_computed_entry()
                decomp_en = round(PDA.get_decomp_and_e_above_hull(
                            ce, allow_negative=True)[1],4)
                redcomp = comp.reduced_composition
                stable_composition_exists = False
                for i in pd.stable_entries:
                    if i.composition.reduced_composition==redcomp:
                        stable_composition_exists = True

                if (decomp_en <= 0.) and stable_composition_exists:
                    #then stable and can proceed as normal
                    logger.debug(
                        "Bulk Computed Entry found to be stable with respect "
                        "to MP Phase Diagram. No mp-id specified, but found "
                        "stable MP composition to exist.")
                    common_approach = True
                elif (decomp_en <= 0.) and not stable_composition_exists:
                    logger.info(
                        "Bulk Computed Entry found to be stable with respect "
                        "to MP Phase Diagram.\nHowever, no stable entry with "
                        "this composition exists in the MP database!\nPlease "
                        "consider submitting the POSCAR to the MP xtaltoolkit,"
                        " so future users will know about this structure:"
                        " https://materialsproject.org/#apps/xtaltoolkit\n"
                        "Manually inserting structure into phase diagram and "
                        "proceeding as normal.")
                    entries.append(ce)
                    pd = PhaseDiagram(entries)
                    PDA = PDAnalyzer(pd)
                    common_approach = True
                elif stable_composition_exists:
                    logger.warning(
                        "Bulk Computed Entry not stable with respect to MP "
                        "Phase Diagram (e_above_hull = %f eV/atom), but found "
                        "stable MP composition to exist.\nProducing chemical "
                        "potentials with respect to stable phase.", decomp_en)
                    common_approach = True
                else:
                    logger.warning(
                        "Bulk Computed Entry not stable with respect to MP "
                        "Phase Diagram (e_above_hull = %f eV/atom) and no "
                        "stable structure with this composition exists in the "
                        "MP database.\nProceeding with atomic chemical "
                        "potentials according to composition position within "
                        "phase diagram.", decomp_en)
                    common_approach = False

            if common_approach:
                for facet in pd.facets:
                    fincomp = comp.reduced_composition
                    eltsinfac=[pd.qhull_entries[j].composition.reduced_composition for j in facet]
                    if fincomp in eltsinfac:
                        chempots = PDA.get_facet_chempots(facet)
                        if len(eltsinfac)!=1:
                            eltsinfac.remove(fincomp)
                        limnom=''
                        for sys in eltsinfac:
                            limnom+=str(sys.reduced_formula)+'-'
                        limnom=limnom[:-1]
                        if len(eltsinfac)==1:
                            limnom+='_rich'
                        print(limnom,chempots)
                        chemdict = {el.symbol:chempots[el] for el in pd.elements}
                        chem_lims[limnom]=chemdict
            else:
                #this uses basic form of creation of facets from
                # initialization of phase diagram object
                from scipy.spatial import ConvexHull
                tmpnew_qdata = list(pd.qhull_data)
                del tmpnew_qdata[-1]
                new_qdata = [[val[i] for i in range(len(val)-1)] for val in tmpnew_qdata]

                tmp_elts = [e for e in pd.elements]
                del tmp_elts[0]
                unstable_qdata_elt = [comp.get_atomic_fraction(el) for el in tmp_elts]
                new_qdata.append(unstable_qdata_elt)

                #take facets of composition space and see if the new
                # composition changes volume of facet
                facets = []
                for facet in pd.facets:
                    tmp_facet = [new_qdata[e] for e in facet]
                    prev_vol = ConvexHull(tmp_facet, qhull_options="QJ i").volume
                    tmp_facet.append(new_qdata[-1])
                    new_vol = ConvexHull(tmp_facet, qhull_options="QJ i").volume
                    if abs(prev_vol-new_vol) < 0.0001:
                        facets.append(facet)

                #now get chemical potentials
                for facet in facets:
                    chempots = PDA.get_facet_chempots(facet)
                    eltsinfac=[pd.qhull_entries[j].composition.reduced_composition for j in facet]
                    limnom=''
                    for sys in eltsinfac:
                        limnom+=str(sys.reduced_formula)+'-'
                    limnom=limnom[:-1]
                    if len(eltsinfac)==1:
                        limnom+='_rich'
                    print(limnom,chempots)
                    chemdict = {el.symbol:chempots[el] for el in pd.elements}
                    chem_lims[limnom]=chemdict

            return chem_lims

        #want to include all sub species within the chemical potential description
        for sub_el in self._substitution_species:
            if sub_el in bulk_species_symbol:
                continue
            else:
                from pymatgen.core import Element
                bulk_species.append(Element(sub_el))
                bulk_species_symbol.append(sub_el)

        chem_lims = get_chempots_from_entries(structure, bulk_species_symbol, bulk_composition)

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
        vbm,gap = self.get_vbm_bandgap()
        output['vbm'] = vbm
        output['gap'] = gap

        return output

