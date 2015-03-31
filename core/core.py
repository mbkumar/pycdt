#!/usr/bin/env python

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.structure import PeriodicSite
from pymatgen.util.io_utils import clean_json


class Defect(object):

    """
    an object with all the info concerning a defect computation: composition+structure, energy,
    correction on energy and name
    """

    def __init__(self, entry_defect, site_in_bulk, charge=0.0,
                 charge_correction=0.0, name=None):
        """
        Args:
            entry_defect: an Entry object corresponding to the defect
            charge: the charge of the defect
            charge_correction: some correction to the energy due to charge
            name: the name of the defect
        """

        self._entry = entry_defect
        self._site = site_in_bulk
        self._charge = charge
        self._charge_correction = charge_correction
        self._name = name
        self._full_name = self._name + "_" + str(charge)

    @property
    def to_dict(self):
        dictio = {}
        dictio['entry'] = self._entry.to_dict
        dictio['site'] = self._site.to_dict
        dictio['charge'] = self._charge
        dictio['charge_correction'] = self._charge_correction
        dictio['name'] = self._name
        dictio['full_name'] = self._full_name
        return clean_json(dictio)

    @staticmethod
    def from_dict(dictio):
        return Defect(ComputedEntry.from_dict(dictio['entry']), PeriodicSite.from_dict(dictio['site']),
                      charge=float(dictio['charge']), charge_correction=float(dictio['charge_correction']),
                      name=dictio['name'])