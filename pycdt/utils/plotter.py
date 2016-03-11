#!/usr/bin/env python

__author__ = "Geoffroy Hautier, Bharat Medasani"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

from pymatgen.util.plotting_utils import get_publication_quality_plot


class DefectPlotter(object):
    """
    Class performing all the typical plots from a defects study
    """

    def __init__(self, analyzer):
        """
        Args:
            analyzer: DefectsAnalyzer object 
        """

        self._analyzer = analyzer

    def get_plot_form_energy(self, xlim=None, ylim=None):
        """
        Formation energy vs Fermi energy plot
        Args:
            xlim:
                Tuple (min,max) giving the range of the x (fermi energy) axis
            ylim:
                Tuple (min,max) giving the range for the formation energy axis
        Returns:
            a matplotlib object

        """
        if xlim is None:
            xlim = (-0.5, self._analyzer._band_gap+1.5)
        max_lim = xlim[1]
        min_lim = xlim[0]
        nb_steps = 10000
        step = (max_lim-min_lim) / nb_steps
        x = [min_lim+step*i for i in range(nb_steps)]
        y = {}
        for t in self._analyzer._get_all_defect_types():
            y_tmp = []
            for x_step in x:
                min = 10000
                for i, dfct in enumerate(self._analyzer._defects):
                    if dfct._name == t:
                        val = self._analyzer._formation_energies[i] + \
                                dfct._charge*x_step
                        if val < min:
                            min = val
                y_tmp.append(min)
            y[t] = y_tmp

        width,height = 12,8
        plt = get_publication_quality_plot(width, height)
        for c in y:
            plt.plot(x, y[c], linewidth=3)
        plt.plot([min_lim, max_lim], [0, 0], 'k-')

        def get_legends(types):
            legends = []
            for name in types:
                for dfct in self._analyzer._defects:
                    if name == dfct._name:
                        sub_str = '_{'+dfct.site.species_string+'}$'
                        if 'vac' in name:
                            base = '$Vac'
                        else: #'sub' in name or 'as' in name or 'antisite' in name:
                            base = '$'+'AS' # dfct.entry.data['substitution_specie'] --> KeyError
                        legend = base + sub_str
                        break
                legends.append(legend)
            return legends
            
        plt.legend(get_legends(y.keys()),fontsize=2*width)
        plt.axvline(x=0.0, linestyle='--', color='k', linewidth=3)
        plt.axvline(x=self._analyzer._band_gap, linestyle='--', color='k', linewidth=3)
        if ylim is not None:
            plt.ylim(ylim)
        plt.xlabel("Fermi energy (eV)",size=2*width)
        plt.ylabel("Defect Formation Energy (eV)",size=2*width)
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
            qi_non_eq.append(
                    self._analyzer.get_non_eq_Ef(t,300,me,mh)['Qi']*1e-6)

        plt = get_publication_quality_plot(12, 8)
        plt.xlabel("temperature (K)")
        plt.ylabel("carrier concentration (cm$^{-3}$)")
        plt.semilogy(temps, qi, linewidth=3.0)
        plt.semilogy(temps, qi_non_eq, linewidth=3)
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
        plt = get_publication_quality_plot(12, 8)
        qi = []
        efs = []
        for ef in [x * 0.01 for x in range(0, 100)]:
            efs.append(ef)
            qi.append(self._analyzer.get_Qi(ef,temp,me,mh)*1e-6)
        plt.ylim([1e14, 1e22])
        return plt.semilogy(efs, qi)
