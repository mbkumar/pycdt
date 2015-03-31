# coding: utf-8

from __future__ import division

"""
This module defines the classes for representing defects.
"""

__author__ = "Geoffroy Hautier, Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani, Geoffroy Hautier"
__email__ = "mbkumar@gmail.com,geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.structure import PeriodicSite


class Defect(object):

    """
    Holds all the info concerning a defect computation: 
    composition+structure, energy, correction on energy and name
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

    def as_dict(self):
        return {'entry': self._entry.as_dict(),
                'site': self._site.as_dict(),
                'charge': self._charge,
                'charge_correction': self._charge_correction,
                'name': self._name,
                'full_name': self._full_name,
                '@module': self.__class__.__module__,
                '@class': self.__class__.__name__}

    @staticmethod
    def from_dict(cls, d):
        return Defect(ComputedEntry.from_dict(d['entry']), 
                      PeriodicSite.from_dict(d['site']),
                      charge=d.get('charge',0.0),
                      charge_correction=d.get('charge_correction',0.0),
                      name=d.get('name',None))
