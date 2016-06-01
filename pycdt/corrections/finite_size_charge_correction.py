"""
This module combines all the finite size supercell charge corrections
into one easy to use class with minimal inputs called ChargeCorrection

The methods implemented are
1) Freysoldt correction for isotropic systems. Includes:
       a) PC energy
       b) potential alignment by planar averaging.
2) Extended Freysoldt or Kumagai correction for anistropic systems. Includes:
       a) anisotropic PC energy
       b) potential alignment by atomic site averaging outside Wigner Seitz radius
3) Sxdefectalign wrapper - python interface to using robust C++ code written by
   Freysoldt et. al (in principle the results are identical to the output of our own
   Freysoldt python code)

If you use the corrections implemented in this module, cite
   Freysoldt, Neugebauer, and Van de Walle,
    Phys. Status Solidi B. 248, 1067-1076 (2011) for isotropic correction
   Kumagai and Oba, Phys. Rev. B. 89, 195205 (2014) for anisotropic correction
   in addition to the pycdt paper
"""

__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import os
import numpy as np
from pymatgen.io.vasp.outputs import Locpot


class ChargeCorrection(object):
    def __init__(self, dielectric_tensor, pure_locpot_path,
            defect_locpot_path, q, pure_locpot=None, defect_locpot=None, pos=None, energy_cutoff=520,
            madetol=0.0001, silence=False, optgamma=None, KumagaiBulk=None ):
        """
        Args:
            dielectric_tensor: Macroscopic dielectric tensor 
                 Include ionic also if defect is relaxed, othewise ion clamped.
                 Can be a matrix, array or scalar.
            pure_locpot_path: Bulk Locpot file path
            defect_locpot_path: Defect Locpot file path
            q: Charge associated with the defect (not homogen. background). Typically integer
            pos: fractional co-ordinates of defect position (just if you are doing sxdefectalign without kumagai or freysoldt first)
            energy_cutoff: Energy for plane wave cutoff (in eV).
                 If not given, Materials Project default 520 eV is used.
            madetol: Tolerance for convergence of energy terms in eV (double or float)
            silence: Flag for disabling/enabling  messages (Bool)
            q_model (QModel object): User defined charge for correction.
                 If not given, highly localized charge is assumed.
            optgamma: (For anisotropic code) If you have previously optimized gamma,
                put gamma value here for speed up calculation slightly.
            KumagaiBulk: This is a class object for Kumagai Anisotropic correction
                    that only needs to be calculated once for each bulk system looked at
        """
        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self._dielectricconst = float(dielectric_tensor)
            self._dieltens = np.diag(np.ones(3) * dielectric_tensor)
        else:
            self._dieltens = np.array(dielectric_tensor)
            self._dielectricconst = np.mean(np.diag(self._dieltens))
        if 'LOCPOT' not in pure_locpot_path:
            print 'pure LOCPOT not in path. appending it to path.'
            self._path_purelocpot = os.path.join(os.path.abspath(pure_locpot_path),'LOCPOT')
        else:
            self._path_purelocpot = os.path.abspath(pure_locpot_path)
        if 'LOCPOT' not in defect_locpot_path:
            print 'defect LOCPOT not in path. appending it to path.'
            self._path_deflocpot = os.path.join(os.path.abspath(defect_locpot_path),'LOCPOT')
        else:
            self._path_deflocpot = os.path.abspath(defect_locpot_path)

        if pure_locpot:
            self._purelocpot = pure_locpot   #actual pure locpot object
        else:
            self._purelocpot = self._path_purelocpot

        if defect_locpot:
            self._deflocpot = defect_locpot  #actual defect locpot object
        else:
            self._deflocpot = self._path_deflocpot
        self._madetol = madetol
        self._q = float(q)
        self._pos = pos
        self._encut = energy_cutoff
        self._silence = silence
        self._KumagaiBulk=KumagaiBulk
        if KumagaiBulk is None:
            self._optgamma=optgamma
        else:
            self._optgamma=KumagaiBulk.gamma

    def freysoldt(self, title=None, axis=0, partflag='All'):
        """
        Args:
            title: set if you want to plot the planar averaged potential
            axis: Specifies axis to average over (zero-defined)
            partflag: four options
                'pc' for just point charge correction, or
               'potalign' for just potalign correction, or
               'All' for both, or
               'AllSplit' for individual parts split up (form [PC,potterm,full])
        """

        from freysoldt_correction import FreysoldtCorrection

        s=FreysoldtCorrection(axis, self._dielectricconst, self._purelocpot,
            self._deflocpot, self._q, energy_cutoff=self._encut, madetol=self._madetol,
            silence=self._silence)

        if partflag in ['All','AllSplit']:
            nomtype='full correction'
            freyval=s.correction(title=title,partflag=partflag)
        elif partflag=='pc':
            nomtype='point charge correction'
            freyval=s.correction(title=title,partflag=partflag)
        elif partflag=='potalign':
            nomtype='potential alignment correction'
            freyval=s.correction(title=title,partflag=partflag)
        else:
            print partflag,' is incorrect potalign type. Must be "All","AllSplit", "pc", or "potalign".'
            return

        #if havent loaded locpots objects or positions of defect then save these for later
        # (good if you want to run freysoldt for multiple axes)
        if (type(s._purelocpot) is Locpot) and (type(self._purelocpot) is not Locpot):
            self._purelocpot=s._purelocpot
        if (type(s._deflocpot) is Locpot) and (type(self._deflocpot) is not Locpot):
            self._deflocpot=s._deflocpot

        print '\n Final Freysoldt',nomtype,'value is ',freyval

        return freyval

    def kumagai(self,title=None, partflag='All'):
        """
        Args:
            title: set if you want to plot the atomic site averaged potential
            partflag: four options
                'pc' for just point charge correction, or
               'potalign' for just pot. align correction, or
               'All' for both, or
               'AllSplit' for individual parts split up (form [PC,potterm,full])
        """

        from kumagai_correction import KumagaiBulkInit, KumagaiCorrection

        if self._KumagaiBulk is None:
             self._KumagaiBulk=KumagaiBulkInit(self._purelocpot, self._dieltens, encut=self._encut,
                tolerance=self._madetol, silence=self._silence, optgamma=self._optgamma)

        if (type(self._KumagaiBulk.bulk_locpot) is Locpot) and (type(self._purelocpot) is not Locpot):
             self._purelocpot=self._KumagaiBulk.bulk_locpot

        s=KumagaiCorrection(self._dieltens, self._purelocpot, self._deflocpot, self._q, gamma=self._KumagaiBulk.gamma,
                g_sum=self._KumagaiBulk.g_sum, energy_cutoff=self._encut, madetol=self._madetol, silence=self._silence)

        if partflag in ['All','AllSplit']:
            nomtype='full correction'
            kumval=s.correction(title=title,partflag=partflag)
        elif partflag=='pc':
            nomtype='point charge correction'
            kumval=s.correction(title=title,partflag=partflag)
        elif partflag=='potalign':
            nomtype='potential alignment correction'
            kumval=s.correction(title=title,partflag=partflag)

        if (type(s.locpot_blk) is Locpot) and (type(self._purelocpot) is not Locpot):
            self._purelocpot=s.locpot_blk
        if (type(s.locpot_def) is Locpot) and (type(self._deflocpot) is not Locpot):
            self._deflocpot=s.locpot_def

        print '\n Final Kumagai',nomtype,'value is ',kumval

        return kumval

    def sxdefect(self,lengths=None,pos=None, axiscalcs=[0],partflag='All',print_pot_flag='written'):
        """
        Args:
            lengths: for length conversion (makes calculation faster)
            pos: specify position for sxdefectalign code
            axiscalcs: Specifies axes to average over (zero-defined)
            partflag: four options
                'pc' for just point charge correction, or
               'potalign' for just potalign correction, or
               'All' for both, or
               'AllSplit' for individual parts split up (form [PC,potterm,full])
        """

        if not pos:
            pos=self._pos #already in fractional co-ordinates from earlier

        if lengths is None and type(self._purelocpot) is Locpot:
            lengths=self._purelocpot.structure.lattice.abc

        from sxdefect_correction import FreysoldtCorrection as SXD

        s=SXD(self._path_purelocpot, self._path_deflocpot, self._q, self._dielectricconst, pos,
                                        self._encut, axiscalcs = axiscalcs, lengths=lengths)

        if partflag in ['All','AllSplit']:
            nomtype='full correction'
        elif partflag=='pc':
            nomtype='point charge correction'
        elif partflag=='potalign':
            nomtype='potential alignment correction'

        sxvals=s.run_correction(print_pot_flag=print_pot_flag,partflag=partflag)

        print '\n Final Sxdefectalign ',nomtype,' correction value is ',sxvals

        return sxvals



if __name__ == '__main__':



