# coding: utf-8
"""
A class for performing analysis of chemical potentials with the grand
canonical linear programming approach
"""
from __future__ import division

__author__ = "Bharat Medasani, Nils Zimmermann, Danny Broberg"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = 'mbkumar@gmail.com'
__date__ = "Sep 14, 2014"

import os
import logging

from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer


class ChemPotAnalyzer(object):
    """
    Post processing for atomic chemical potentials used in defect
    calculations. (note this could be MP associated OR associated with 
    other inputs calculated by user?)

    Makes use of Materials Project pre-computed data to generate
    needed information for chem pots
        1) If using GGA-PBE vasp then can give numerical values for
        chem pots in different growth conditions
        2) If not using GGA-PBE Vasp then can give all needed structures 
        needed for computing chemical potentials
    """

    def __init__(self, bulk_composition, sub_species=set(), entries={}):
        """
        TODO: could have bulk entry object as input for faster parsing?

        Args:
            bulk_composition : Composition of bulk as a pymatgen Composition
                object. This and mapi_key are only actual required input for
                generating set of chemical potentials from Materials Project
                database
            subs_species : set of elemental species that are extrinsic to
                structure defaults to No substitutions needed.
            entries: pymatgen ComputedEntry objects to build phase diagram
        """
        self.bulk_composition = bulk_composition
        self.bulk_species_symbol = [s.symbol for s in bulk_composition.elements]
        self.sub_species = sub_species
        self.redcomp = bulk_composition.reduced_composition
        self.entries = entries
        self.bulk_ce = None  # could be improved by having direct loading...

    @staticmethod
    def from_dict(cls, d):
        """
        Create a ChemPotAnalyzer from a dictionary.
        Useful for setting up with previous chemical_data and reduces
        need for querying

        Args:
            d:  python dict
        """
        cpa = ChemPotAnalyzer(d['bulk_composition'],
                              sub_species=d['sub_species'],
                              entries=d['entries'])
        return cpa

    def as_dict(self):
        """""
        Json-serialization dict representation of the ChemPotAnalyzer.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             'bulk_composition': self.bulk_composition,
             'sub_species': self.sub_species,
             'bulk_species_symbol': self.bulk_species_symbol,
             'entries': self.entries}

        return d

    def get_mp_entries(self, mpid=None, mapi_key=None):
        """
        This queries MP database for computed entries according to
        input bulk and sub elements of interest

        Args:
            mpid (str): Structure id of the system in the MP databse.
            mapi_key (str): Materials API key to access database
                (if not in ~/.pmgrc.yaml already)
        """
        logger = logging.getLogger(__name__)

        with MPRester(api_key=mapi_key) as mp:
            self.entries['bulk_derived'] = mp.get_entries_in_chemsys(
                self.bulk_species_symbol)

            if mpid:
                self.bulk_ce = mp.get_entry_by_material_id(mpid)

        if not self.entries:
            msg = "Could not fetch bulk entries for atomic chempots!" \
                  "MPRester query error."
            logger.warning(msg)
            raise ValueError(msg)

        # now compile substitution entries
        self.entries['subs_set'] = dict()
        bulk_entry_set = [entry.entry_id for entry in
                          self.entries['bulk_derived']]
        for sub_el in self.sub_species:
            els = self.bulk_species_symbol + [sub_el]
            with MPRester(api_key=mapi_key) as mp:
                sub_entry_set = mp.get_entries_in_chemsys(els)
            if not sub_entry_set:
                msg = "Could not fetch sub entries for {} atomic chempots! " \
                      "Encountered MPRester query error".format(sub_el)
                logger.warning(msg)
                raise ValueError(msg)

            fin_sub_entry_set = []
            for entry in sub_entry_set:
                if entry.entry_id not in bulk_entry_set:
                    fin_sub_entry_set.append(entry)
            # All entries apart from the bulk entry set
            self.entries['subs_set'][sub_el] = fin_sub_entry_set

        return

    def get_mp_entries_from_symbols(self, list_spec_symbol, mapi_key=None):
        """
        Gets entries list from MP database based on entries of list_spec_symbol
        """
        logger = logging.getLogger(__name__)
        with MPRester(api_key=mapi_key) as mp:
            self.entries = mp.get_entries_in_chemsys(list_spec_symbol)
            self.recent_list_specs = list_spec_symbol
        if not self.entries:
            msg = "Could not fetch entries for atomic chempots! " \
                  "MPRester query error."
            logger.warning(msg)
            raise ValueError(msg)
        return

    def analyze_GGA_chempots(self, bulk_entry=None, mpid=None, mapi_key=None,
                             full_sub_approach=False):
        """
        For calculating GGA-PBE atomic chemical potentials by using
            Materials Project pre-computed data

        Args for input :
            (Among bulk_computed_entry, root_flr and mpid, only one of them is
            required)
            bulk_entry: Pymatgen ComputedStructureEntry object for
                bulk supercell
            mpid (str): Materials Project ID of bulk structure;
                format "mp-X", where X is an integer;
            mapi_key (str): Materials API key to access database
                (if not in ~/.pmgrc.yaml already)
            full_sub_approach: generate chemical potentials by looking at
                full phase diagram (setting to True is really NOT recommended
                if subs_species set has more than one element in it...)

        Outline for how this code retrieves atomic chempots from Materials
        Project (MP) entries in a phase diagram (PD) object:
          1) check stability of computed entry w.r.t phase diagram
                generated from MP
          2) If stable with respect to phase diagram, then:
             i) if mp-id given and it is in the stable entry list,
                proceed normally ("common approach")
             ii) if mp-id not given, prints message to user about possibility
                for manual submission page on MP website, then
                manually inserts the computed object into the local PD.
                Generate chempots from this PD (note at this point there 
                is no gurantee that the structure is actually unique, as 
                it could be a computational error in the DFT energy which 
                is slightly lower than MP calculated values)
          3) If not-stable with respect to phase diagram, then:
             i) check to see if composition exists among the structures
                in the stable list of the PD
             ii) if a stable and identical composition exists in PD
                then print warning but continue as if it were stable
                (chem pots will depend on the stable phase)
             iii) if no stable and identical composition exists in PD,
                then print warning and find facets that the bulk composition
                exists inside of

        Note on full_sub_approach:
            the default approach for subs is to only consider facets
            defined by N-2 phases with strictly elements from the BULK
            composition, and 1 sub_element(+possibly bulk-composition
            element) derived phases (along with the condition for all
            chemical potentials to be defined by the bulk entry, creating
            N equations to be solved for N atomic chemical potentials).
            This speeds up analysis SIGNFICANTLY when analyzing several
            substitutional species at once. It is essentially the
            assumption the majority of the elements in the total
            composition will be from the native species present rather
            than the sub species (a good approximation). If you would
            prefer to consider the full phase diagram (not recommended
            unless you have 1 or 2 substitutional defects), then set
            full_sub_approach to True.
        """
        if not self.entries:
            self.get_mp_entries(mapi_key=mapi_key)

        logger = logging.getLogger(__name__)
        # first get the computed entry
        if bulk_entry:
            self.bulk_ce = bulk_entry
        elif mpid:
            with MPRester(api_key=mapi_key) as mp:
                self.bulk_ce = mp.get_entry_by_material_id(mpid)
        else:
            msg = "No able to load computed entry. Cannot parse chemical " \
                  "potentials for job."
            logger.warning(msg)
            raise ValueError(msg)

        # figure out how system should be treated for chemical potentials
        # based on phase diagram
        entry_list = self.entries['bulk_derived']
        pd = PhaseDiagram(entry_list)
        pda = PDAnalyzer(pd)
        full_idlist = [i.entry_id for i in pd.qhull_entries]
        stable_idlist = [i.entry_id for i in pd.stable_entries]

        if mpid:
            if (mpid in full_idlist) and (mpid in stable_idlist):
                logger.debug("Verified that mp-id is stable within Materials "
                             "Project {} phase diagram".format('-'.join(
                                self.bulk_species_symbol)))
                common_approach = True
            elif (mpid in full_idlist) and not (mpid in stable_idlist):
                common_approach = False
                for i in pd.stable_entries:
                    if i.composition.reduced_composition == self.redcomp:
                        logger.warning(
                            "Input mp-id {} is unstable. Stable "
                            "composition mp-id found to be {}".format(
                                mpid, i.entry_id))
                        logger.warning(
                            "Proceeding with atomic chemical potentials "
                            "with respect to stable phase.")
                        common_approach = True
                if not common_approach:
                    logger.warning(
                        "Input mp-id {} is unstable. No stable structure "
                        "with same composition exists".format(mpid))
                    logger.warning(
                        "Proceeding with atomic chemical potentials "
                        "according to composition position within phase "
                        "diagram.")
            else:
                logger.warning(
                    "Specified mp-id {} not found in MP phase "
                    "diagram. Reverting to assumption that mp-id is not "
                    "known.".format(mpid))
                mpid = None

        if not mpid:
            decomp_en = round(pda.get_decomp_and_e_above_hull(
                                    self.bulk_ce, allow_negative=True)[1],
                              4)
            stable_composition_exists = False
            for i in pd.stable_entries:
                if i.composition.reduced_composition == self.redcomp:
                    stable_composition_exists = True

            if (decomp_en <= 0) and stable_composition_exists:
                # then stable and can proceed as normal
                logger.debug(
                    "Bulk Computed Entry found to be stable with respect "
                    "to MP Phase Diagram. No mp-id specified, but found "
                    "stable MP composition to exist.")
                common_approach = True
            elif (decomp_en <= 0) and not stable_composition_exists:
                logger.info(
                    "Bulk Computed Entry found to be stable with respect "
                    "to MP Phase Diagram.\nHowever, no stable entry with "
                    "this composition exists in the MP database!\nPlease "
                    "consider submitting the POSCAR to the MP xtaltoolkit,"
                    " so future users will know about this structure:"
                    " https://materialsproject.org/#apps/xtaltoolkit\n"
                    "Manually inserting structure into phase diagram and "
                    "proceeding as normal.")
                entry_list.append(self.bulk_ce)
                # pd = PhaseDiagram(self.entries)
                # PDA = PDAnalyzer(pd)
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

        if full_sub_approach:
            # Add all possible entries to entries list for phase diagram...
            # (not recommended)
            for sub, subentries in self.entries['subs_set'].items():
                for subentry in subentries:
                    entry_list.append(subentry)

        pd = PhaseDiagram(entry_list)
        pda = PDAnalyzer(pd)
        chem_lims = self.get_chempots_from_pda(
            pda, common_approach=common_approach)

        if not full_sub_approach:
            finchem_lims = {}  # this will be final chem_lims dictionary
            for key in chem_lims.keys():
                # TODO: Pretty sure this splitting up of strings was
                # TODO: neccessary earlier on, but probably can do this
                # TODO: in a prettier way now...
                face_list = key.split('-')
                blk, blknom, subnom = self.diff_bulk_sub_phases(face_list)
                finchem_lims[blknom] = {}
                finchem_lims[blknom] = chem_lims[key]

            # Now consider adding single elements to extend the phase diagram,
            # adding new additions to chemical potentials ONLY for the cases
            # where the phases in equilibria are those from the bulk phase
            # diagram. This is essentially the assumption that the majority of
            # the elements in the total composition will be from the native
            # species present rather than the sub species (a good approximation)
            for sub_el in self.sub_species:
                sub_specie_entries = entry_list[:]
                for entry in self.entries['subs_set'][sub_el]:
                    sub_specie_entries.append(entry)

                pd = PhaseDiagram(sub_specie_entries)
                pda = PDAnalyzer(pd)
                chem_lims = self.get_chempots_from_pda(pda)

                for key in chem_lims.keys():
                    face_list = key.split('-')
                    blk, blknom, subnom = self.diff_bulk_sub_phases(
                        face_list, sub_el=sub_el)
                    # if one less than number of bulk species then can be
                    # grouped with rest of structures
                    if len(blk)+1 == len(self.bulk_species_symbol):
                        if blknom not in finchem_lims.keys():
                            finchem_lims[blknom] = chem_lims[key]
                        else:
                            finchem_lims[blknom][sub_el] = \
                                chem_lims[key][sub_el]
                        if 'name-append' not in finchem_lims[blknom].keys():
                            finchem_lims[blknom]['name-append'] = subnom
                        else:
                            finchem_lims[blknom]['name-append'] += '-' + subnom
                    else:
                        # if chem pots determined by two (or more) sub-specie 
                        # containing phases, skip this facet!
                        continue
            chem_lims = finchem_lims.copy()

        return chem_lims

    def get_chempots_from_pda(self, pda, common_approach=True):
        # pass in a phase diagram Analyzer object and output chemical potentials
        # common_approach=True assumes that the bulk_entry exists while
        # common_approach = False determines chemical potentials based on the
        # facets that would contain the composition of interest
        chem_lims = {}
        pd = pda._pd
        if common_approach:
            for facet in pda._pd.facets:
                eltsinfac = [
                    pd.qhull_entries[j].composition.reduced_composition
                    for j in facet]
                if self.redcomp in eltsinfac:
                    chempots = pda.get_facet_chempots(facet)
                    if len(eltsinfac) != 1:
                        eltsinfac.remove(self.redcomp)
                    limnom = '-'.join(sys.reduced_formula for sys in eltsinfac)
                    if len(eltsinfac) == 1:
                        limnom += '_rich'
                    print(limnom, chempots)
                    chemdict = {
                        el.symbol: chempots[el] for el in pd.elements}
                    chem_lims[limnom] = chemdict
        else:
            # this uses basic form of creation of facets from initialization
            # of phase diagram object to find which facets of phase diagram
            # contain the composition of interest
            from scipy.spatial import ConvexHull
            tmpnew_qdata = list(pd.qhull_data)
            del tmpnew_qdata[-1]
            new_qdata = [[val[i] for i in range(len(val)-1)]
                         for val in tmpnew_qdata]

            tmp_elts = [e for e in pd.elements]
            del tmp_elts[0]
            unstable_qdata_elt = [
                self.bulk_composition.get_atomic_fraction(el)
                for el in tmp_elts]
            new_qdata.append(unstable_qdata_elt)

            # take facets of composition space and see if the new
            # composition changes volume of facet
            facets = []
            for facet in pd.facets:
                tmp_facet = [new_qdata[e] for e in facet]
                prev_vol = ConvexHull(
                    tmp_facet, qhull_options="QJ i").volume
                tmp_facet.append(new_qdata[-1])
                new_vol = ConvexHull(tmp_facet, qhull_options="QJ i").volume
                if abs(prev_vol-new_vol) < 0.0001:
                    facets.append(facet)

            # now get chemical potentials
            for facet in facets:
                mus = pda.get_facet_chempots(facet)
                eltsinfac = [
                    pd.qhull_entries[j].composition.reduced_composition
                    for j in facet]
                limnom = '-'.join(sys.reduced_formula for sys in eltsinfac)
                if len(eltsinfac) == 1:
                    limnom += '_rich'
                print(limnom, mus)
                chemdict = {el.symbol: mus[el] for el in pd.elements}
                chem_lims[limnom] = chemdict

        return chem_lims

    def diff_bulk_sub_phases(self, face_list, sub_el=None):
        # method for pulling out phases within a facet of a phase diagram
        # which may include a substitutional element...
        # face_list is an array of phases in a facet
        # sub_el is the element to look out for within the face_list array
        blk = []
        sub_spcs = []
        for face in face_list:
            if sub_el:
                if sub_el in face:
                    sub_spcs.append(face)
                else:
                    blk.append(face)
            else:
                blk.append(face)
        blk.sort()
        sub_spcs.sort()
        blknom = '-'.join(blk)
        subnom = '-'.join(sub_spcs)
        return blk, blknom, subnom

    def get_chempots_from_composition(self, mapi_key=None):
        """
        A simple method for getting GGA-PBE chemical potentials JUST
        from the composition information (Note: this only works if the
        composition already exists in the MP database)
        """
        if not self.entries:
            self.get_mp_entries(mapi_key=mapi_key)

        logger = logging.getLogger(__name__)
        # retrieve the most stable mp-id with the given composition
        lowest_energy_mpid = None
        lowest_energy = 1000.
        for i in self.entries['bulk_derived']:
            if (i.composition.reduced_composition == self.redcomp) \
                    and (i.energy_per_atom < lowest_energy):
                lowest_energy_mpid = i.entry_id
                lowest_energy = i.energy_per_atom

        if not lowest_energy_mpid:
            msg = "Not able to find an mpid for composition of interest. " \
                  "Cannot generate chempots without a computed entry."
            logger.warning(msg)
            raise ValueError(msg)
        else:
            mu = self.analyze_GGA_chempots(mpid=lowest_energy_mpid,
                                           mapi_key=mapi_key)

        return mu
