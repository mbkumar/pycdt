#!/usr/bin/env python

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.util.plotting_utils import get_publication_quality_plot


class DefectPlotter(object):
    """
    class performing all the typical plots from a defects study
    """

    def __init__(self, analyzer):
        self._analyzer = analyzer

    def get_plot_form_energy(self, xlim=None, ylim=None):
        """
        plots the typical formation energy vs Fermi energy plot
        Args:
            xlim:
                a tuple of (min, max) giving the edge of the x (fermi energy) axis
            ylim:
                a tuple of (min, max) giving the limits for the formation energy axis
        Returns:
            a matplotlib object

        """
        if xlim is None:
            xlim = (-0.5, self._analyzer._band_gap+1.0)
        max_lim = xlim[1]
        min_lim = xlim[0]
        nb_steps = 10000
        step = (max_lim - min_lim) / nb_steps
        x = [min_lim + step * i for i in range(nb_steps)]
        y = {}
        for t in self._analyzer._get_all_defect_types():
            y_tmp = []
            for x_step in x:
                min = 10000
                for i in range(len(self._analyzer._defects)):
                    if self._analyzer._defects[i]._name == t:
                        if self._analyzer._formation_energies[i] + self._analyzer._defects[i]._charge * x_step < min:
                            min = self._analyzer._formation_energies[i] + self._analyzer._defects[i]._charge * x_step
                y_tmp.append(min)
            y[t] = y_tmp
        plt = get_publication_quality_plot(12, 8)
        for c in y:
            plt.plot(x, y[c], linewidth=5)
        plt.plot([min_lim, max_lim], [0, 0], 'k-')
        plt.legend(y.keys())
        plt.axvline(x=0.0, linestyle='--', color='k', linewidth=3)
        if ylim is not None:
            plt.ylim(ylim)
        plt.xlabel("Fermi energy (eV)")
        plt.ylabel("Defect Formation Energy (eV)")
        return plt

    def plot_conc_temp(self, me=[1.0, 1.0, 1.0], mh=[1.0, 1.0, 1.0]):
        """
        plot the concentration of carriers vs temperature both in eq and non-eq after quenching at 300K
        Args:
            me:
                the effective mass for the electrons as a list of 3 eigenvalues
            mh:
                the effective mass for the holes as a list of 3 eigenvalues
        Returns;
            a matplotlib object

        """
        temps = [i*100 for i in range(3,20)]
        qi = []
        qi_non_eq = []
        for t in temps:
            qi.append(self._analyzer.get_eq_Ef(t,me,mh)['Qi']*1e-6)
            qi_non_eq.append(self._analyzer.get_non_eq_Ef(t,300,me,mh)['Qi']*1e-6)

        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        plt.xlabel("temperature (K)")
        plt.ylabel("carrier concentration (cm$^{-3}$)")
        plt.semilogy(temps,qi,linewidth=3.0)
        plt.semilogy(temps,qi_non_eq,linewidth=3)
        plt.legend(['eq','non-eq'])
        return plt

    def plot_carriers_ef(self, temp=300, me=[1.0, 1.0, 1.0], mh=[1.0, 1.0, 1.0]):
        """
        plot carrier concentration in function of the fermi energy
        Args:
            temp:
                temperature
            me:
                the effective mass for the electrons as a list of 3 eigenvalues
            mh:
                the effective mass for the holes as a list of 3 eigenvalues
        Returns:
            a matplotlib object
        """
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8)
        qi = []
        efs = []
        for ef in [x * 0.01 for x in range(0, 100)]:
            efs.append(ef)
            qi.append(self._analyzer.get_Qi(ef,temp,me,mh)*1e-6)
        plt.ylim([1e14, 1e22])
        return plt.semilogy(efs,qi)