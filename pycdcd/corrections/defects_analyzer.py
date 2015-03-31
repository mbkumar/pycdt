#!/usr/bin/env python


__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.periodic_table import Element
import math
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.util.io_utils import clean_json
from mpcollab.defects.core import Defect

#some constants
kb = 8.6173324e-5
hbar = 6.58211928e-16
conv = (math.sqrt((9.1*1e-31)**3)*math.sqrt((1.6*1e-19)**3))/((1.05*1e-34)**3)

class DefectsAnalyzer(object):
    """
    a class aimed at performing standard analysis of defects
    """
    def __init__(self, entry_bulk, e_vbm, mu_elts, band_gap):
        """
        Args:
            entry_bulk:
                the bulk data as an Entry
            e_vbm:
                the energy of the vbm (in eV)
            mu_elts:
                a dictionnary of {Element:value} giving the chemical
                potential of each element
            band_gap:
                the band gap (in eV)
        """
        self._entry_bulk = entry_bulk
        self._e_vbm = e_vbm
        self._mu_elts = mu_elts
        self._band_gap = band_gap
        self._defects = []
        self._formation_energies = []

    @property
    def to_dict(self):
        dictio = {}
        dictio['entry_bulk'] = self._entry_bulk.to_dict
        dictio['e_vbm'] = self._e_vbm
        dictio['mu_elts'] = self._mu_elts
        dictio['band_gap'] = self._band_gap
        dictio['defects'] = [d.to_dict for d in self._defects]
        dictio['formation_energies'] = self._formation_energies
        return clean_json(dictio)

    @staticmethod
    def from_dict(dictio):
        analyzer = DefectsAnalyzer(ComputedStructureEntry.from_dict(dictio['entry_bulk']), dictio['e_vbm'], {Element(e):dictio['mu_elts'][e] for e in dictio['mu_elts']}, dictio['band_gap'])
        for d in dictio['defects']:
            analyzer.add_defect(Defect.from_dict(d))
        return analyzer

    def add_defect(self, defect):
        """
        add a defect to the analyzer
        Args:
            defect:
                a Defect object
        """
        self._defects.append(defect)
        self._compute_form_en()

    def _get_all_defect_types(self):
        to_return = []
        for d in self._defects:
            if d._name not in to_return: to_return.append(d._name)
        return to_return

    def _compute_form_en(self):
        """
        compute the formation energies for all defects in the analyzer
        """
        self._formation_energies = []
        for d in self._defects:
            multiplier = None
            if math.floor((d._entry.composition.num_atoms + 1) / self._entry_bulk.composition.num_atoms) \
                    == (d._entry.composition.num_atoms+1) / self._entry_bulk.composition.num_atoms:
                multiplier = (d._entry.composition.num_atoms+1)/self._entry_bulk.composition.num_atoms
            elif math.floor((d._entry.composition.num_atoms-1)/self._entry_bulk.composition.num_atoms) \
                    == (d._entry.composition.num_atoms-1)/self._entry_bulk.composition.num_atoms:
                multiplier = (d._entry.composition.num_atoms-1)/self._entry_bulk.composition.num_atoms
            elif math.floor(d._entry.composition.num_atoms/self._entry_bulk.composition.num_atoms) \
                    == d._entry.composition.num_atoms/self._entry_bulk.composition.num_atoms:
                multiplier = d._entry.composition.num_atoms/self._entry_bulk.composition.num_atoms
            #go through each element in the defect and see how to "compensate" it with the chemical potentials
            mu_needed_coeffs = {}
            for elt in d._entry.composition.elements:
                if d._entry.composition[elt] > multiplier * self._entry_bulk.composition[elt]:
                    mu_needed_coeffs[elt] = -1.0
                if d._entry.composition[elt] < multiplier * self._entry_bulk.composition[elt]:
                    mu_needed_coeffs[elt] = 1.0

            sum_mus = 0.0
            for elt in mu_needed_coeffs:
                sum_mus = sum_mus + mu_needed_coeffs[elt] * self._mu_elts[elt]

            self._formation_energies.append(d._entry.energy - multiplier * self._entry_bulk.energy + sum_mus
                                            + d._charge * self._e_vbm + d._charge_correction)

    def correct_bg_simple(self, vbm_correct, cbm_correct):
        """
        correct the band gap in the analyzer.
        We assume the defects level remain the same when moving the band edges
        Args:
            vbm_correct:
                The correction on the vbm as a positive number. e.g., if the VBM
                goes 0.1 eV down vbm_correct=0.1
            cbm_correct:
                The correction on the cbm as a positive number. e.g., if the CBM
                goes 0.1 eV up cbm_correct=0.1

        """
        self._band_gap = self._band_gap + cbm_correct + vbm_correct
        self._e_vbm = self._e_vbm - vbm_correct
        self._compute_form_en()

    def correct_bg(self, dict_levels, vbm_correct, cbm_correct):
        """
        correct the band gap in the analyzer and make sure the levels move
        accordingly.
        There are two types of defects vbm_like and cbm_like and we need
        to provide a formal oxidation state
        The vbm-like will follow the vbm and the cbm_like the cbm. If nothing
        is specified the defect transition level does not move
        Args:
            dict_levels: a dictionnary of type {defect_name:
            {'type':type_of_defect,'q*':formal_ox_state}}
            Where type_of_defect is a string: 'vbm_like' or 'cbm_like'
        """


        self._band_gap = self._band_gap + cbm_correct + vbm_correct
        self._e_vbm = self._e_vbm - vbm_correct
        self._compute_form_en()
        for i in range(len(self._defects)):
            if not self._defects[i]._name in dict_levels:
                continue

            if dict_levels[self._defects[i]._name]['type'] == 'vbm_like':
                z = self._defects[i]._charge-dict_levels[self._defects[i]._name]['q*']
                self._formation_energies[i] = self._formation_energies[i] + z * vbm_correct
            if dict_levels[self._defects[i]._name]['type'] == 'cbm_like':
                z = dict_levels[self._defects[i]._name]['q*']-self._defects[i]._charge
                self._formation_energies[i] = self._formation_energies[i] + z * cbm_correct

    def _get_form_energy(self, ef, i):
        return self._formation_energies[i] + self._defects[i]._charge * ef

    def get_defects_concentration(self, temp=300, ef=0.0):
        """
        get the defect concentration for a temperature and Fermi level
        Args:
            temp:
                the temperature in K
            Ef:
                the fermi level in eV (with respect to the VBM)
        Returns:
            a list of dictionary of {'name':name of defect, 'charge':charge of defect
                                    'conc': concentration of defects in m-3}
        """
        conc=[]
        struct = SymmetryFinder(self._entry_bulk.structure,symprec=1e-1).get_symmetrized_structure()
        i = 0
        for d in self._defects:
            target_site=None
            #TODO make a better check this large tol. is weird
            for s in struct.sites:
                if abs(s.frac_coords[0]-d._site.frac_coords[0]) < 0.1 \
                    and abs(s.frac_coords[1]-d._site.frac_coords[1]) < 0.1 \
                    and abs(s.frac_coords[2]-d._site.frac_coords[2]) < 0.1:
                    target_site=s
            n = len(struct.find_equivalent_sites(target_site))*1e30/(struct.volume)
            conc.append({'name': d._name, 'charge': d._charge,
                         'conc': n*math.exp(-self._get_form_energy(ef, i)/(kb*temp))})
            i += 1
        return conc

    def _get_dos(self, e, m1, m2, m3, e_ext):
        return math.sqrt(2)/(math.pi**2*hbar**3)*math.sqrt(m1*m2*m3)*math.sqrt(e-e_ext)

    def _get_dos_fd_elec(self, e, ef, t, m1, m2, m3):
        return conv*(2.0/(math.exp((e-ef)/(kb*t))+1))*(math.sqrt(2)/(math.pi**2)) \
            * math.sqrt(m1 * m2 * m3)*math.sqrt(e - self._band_gap)

    def _get_dos_fd_hole(self, e, ef, t, m1, m2, m3):
        return conv*((math.exp((e-ef)/(kb*t))/(math.exp((e-ef)/(kb*t))+1))) * \
            (2.0 * math.sqrt(2)/(math.pi ** 2)) * math.sqrt(m1 * m2 * m3)*math.sqrt(-e)

    def _get_qd(self,ef,t):
        summation = 0.0
        for d in self.get_defects_concentration(t, ef):
            summation = summation + d['charge'] * d['conc']
        return summation

    def _get_qi(self, ef, t, m_elec, m_hole):
        from scipy import integrate
        return -integrate.quad(lambda e: self._get_dos_fd_elec(e, ef, t, m_elec[0], m_elec[1], m_elec[2]),
                               self._band_gap, 5+self._band_gap)[0] +\
                               integrate.quad(lambda e: self._get_dos_fd_hole(e, ef, t, m_hole[0], m_hole[1],
                                                                              m_hole[2]), -5, 0.0)[0]

    def _get_qtot(self, ef, t, m_elec, m_hole):
        return self._get_qd(ef, t) + self._get_qi(ef, t, m_elec, m_hole)

    def get_eq_ef(self, t, m_elec, m_hole):
        """
        access to equilibrium values of Fermi level and concentrations in defects and carriers
        obtained by self-consistent solution of charge balance + defect and carriers concentrations
        Args:
            t:
                temperature in K
            m_elec:
                electron effective mass as a list of 3 values (3 eigenvalues for the tensor)
            m_hole::
                hole effective mass as a list of 3 values (3 eigenvalues for the tensor)
        Returns:
            a dictionary with {'ef':eq fermi level,'Qi': the concentration of carriers
            (positive for holes, negative for e-) in m^-3,'conc': the concentration of defects
             as a list of dictionary
        """
        from scipy.optimize import bisect
        e_vbm = self._e_vbm
        e_cbm = self._e_vbm+self._band_gap
        ef = bisect(lambda e:self._get_qtot(e,t,m_elec,m_hole), 0, self._band_gap)
        return {'ef': ef, 'Qi': self._get_qi(ef, t, m_elec, m_hole),
                'QD': self._get_qd(ef,t), 'conc': self.get_defects_concentration(t, ef)}

    def get_non_eq_ef(self, tsyn, teq, m_elec, m_hole):
        """
        access to the non-equilibrium values of Fermi level and concentrations in defects and carriers
        obtained by self-consistent solution of charge balance + defect and carriers concentrations

        Implemented following Sun, R., Chan, M. K. Y., Kang, S., and Ceder, G. (2011). doi:10.1103/PhysRevB.84.035212

        Args:
            tsyn:
                the synthesis temperature in K
            teq:
                the temperature of use in K
            m_elec:
                electron effective mass as a list of 3 values (3 eigenvalues for the tensor)
            m_hole:
                hole effective mass as a list of 3 values (3 eigenvalues for the tensor)
        Returns:
            a dictionary with {'ef':eq fermi level,'Qi': the concentration of carriers
            (positive for holes, negative for e-) in m^-3,'conc': the concentration of defects
             as a list of dictionary
        """
        from scipy.optimize import bisect
        eqsyn = self.get_eq_ef(tsyn, m_elec, m_hole)
        cd = {}
        for c in eqsyn['conc']:
            if c['name'] in cd:
                cd[c['name']]=cd[c['name']]+c['conc']
            else:
                cd[c['name']]=c['conc']
        ef = bisect(lambda e:self._get_non_eq_qtot(cd, e, teq, m_elec, m_hole), -1.0, self._band_gap+1.0)
        return {'ef':ef, 'Qi':self._get_qi(ef, teq, m_elec, m_hole),'conc_syn': eqsyn['conc'],
                'conc':self._get_non_eq_conc(cd,ef,teq)}

    def _get_non_eq_qd(self, cd, ef, t):
        sum_tot = 0.0
        for n in cd:
            sum_d = 0.0
            sum_q = 0.0
            i = 0
            for d in self._defects:
                if d._name == n:
                    sum_d += math.exp(-self._get_form_energy(ef, i) / (kb * t))
                    sum_q += d._charge * math.exp(-self._get_form_energy(ef, i) / (kb * t))
                i += 1
            sum_tot += cd[n] * sum_q / sum_d
        return sum_tot

    def _get_non_eq_conc(self,cd,ef,t):
        sum_tot = 0.0
        res=[]
        for n in cd:
            sum_tot = 0
            i = 0
            for d in self._defects:
                if d._name == n:
                    sum_tot += math.exp(-self._get_form_energy(ef,i)/(kb*t))
                i += 1
            i=0
            for d in self._defects:
                if d._name == n:
                    res.append({'name':d._name,'charge':d._charge,
                                'conc':cd[n]*math.exp(-self._get_form_energy(ef,i)/(kb*t))/sum_tot})
                i += 1
        return res

    def _get_non_eq_qtot(self,cd,ef,t, m_elec, m_hole):
        return self._get_non_eq_qd(cd, ef, t)+self._get_qi(ef, t, m_elec, m_hole)