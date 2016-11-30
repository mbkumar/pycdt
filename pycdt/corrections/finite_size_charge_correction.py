"""
This module contains wrappers to the finite size supercell charge corrections
implemented 

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
from pycdt.corrections.kumagai_correction import KumagaiBulkInit, KumagaiCorrection

def get_correction_freysoldt(defect, bulk_entry, epsilon, title = None):
    """
    Function to compute the isotropic freysoldt correction for each defect.
    Args:
        defect: ComputedDefect object
        bulk_entry: ComputedStructureEntry corresponding to bulk OR bulk Locpot Object
        epsilon: dielectric constant
    """
    if type(bulk_entry) is Locpot:
        locpot_blk = bulk_entry
        locpot_path_blk = ''
    else:
        locpot_blk = None
        locpot_path_blk = bulk_entry.data['locpot_path']
    locpot_path_def = defect.entry.data['locpot_path']
    charge = defect.charge
    #frac_coords = defect.site.frac_coords  #maybe should be using this
    encut = defect.entry.data['encut']
    if not charge:
        print('charge is zero so charge correction is zero')
        return (0.,bulk_entry)

    #if either locpot is already loaded then load pure_locpot= or defect_locpot=
    # if you want to load position then can load it with pos=
    #if want to to change energy tolerance for correction convergence then change madetol= (default is 0.0001)
    # (for kumagai) if known optgamma, set optgamma=, if KumagaiBulk already initialized then set KumagaiBulk=
    corr_meth = ChargeCorrection(epsilon,
            locpot_path_blk, locpot_path_def, charge,
            pure_locpot = locpot_blk, #for quicker loading of bulk locpot objects...
            energy_cutoff = encut,
            silence=False)

    #could do an averaging over three axes but one axis works fine isotropic systems
    corr_val = corr_meth.freysoldt(title=title, axis=0, partflag='All')

    return (corr_val,corr_meth._purelocpot)


def get_correction_kumagai(defect, path_blk, bulk_init, bulk_locpot=None, 
                           title=None):
    """
    Function to compute the correction for each defect.
    Args:
        defect: ComputedDefect object
        path_blk: location to Bulk folder
        bulk_init: KumagainBulkInit class object
        bulk_locpot: BulkLocpot object 
                (if already loaded, otherwise will load from path_blk)
        type: 
            "freysoldt": Freysoldt correction for isotropic crystals
            "kumagai": modified Freysoldt or Kumagai for anisotropic crystals

    notes for ChargeCorrection class below:
    if either locpot already loaded then load pure_locpot= or defect_locpot=
    if you want to load position then can load it with pos=
    if want to to change energy tolerance for correction convergence then 
        change madetol= (default is 0.0001)
    (if known optgamma, set optgamma=, if KumagaiBulk already initialized 
    then set KumagaiBulk=
    """
    epsilon = bulk_init.epsilon
    charge = defect.charge
    encut = defect.entry.data['encut']

    outcar_path_blk = os.path.join(path_blk,'OUTCAR')
    locpot_path_def = defect.entry.data['locpot_path']
    dpat,dloc = os.path.split(locpot_path_def)
    outcar_path_def = os.path.join(dpat,'OUTCAR')

    if os.path.exists(outcar_path_blk) and os.path.exists(outcar_path_def):
        s = KumagaiCorrection(epsilon, charge, bulk_init.gamma, 
                bulk_init.g_sum, bulk_init.structure, defect.entry.structure, 
                energy_cutoff=encut, madetol=0.0001, 
                bulk_outcar=outcar_path_blk, defect_outcar=outcar_path_def)
    else:
        if not bulk_locpot:
            bulk_locpot = Locpot.from_file(path_blk)
        s = KumagaiCorrection(epsilon, charge, bulk_init.gamma,
                bulk_init.g_sum, bulk_init.structure, defect.entry.structure,
                energy_cutoff=encut, madetol=0.0001, 
                bulk_locpot=bulk_locpot, defect_locpot=locpot_path_def)

    kumval = s.correction(title=title, partflag='All')
    print('\n Kumagai Correction value is ', kumval)
    return kumval


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
            print('pure LOCPOT not in path. appending it to path.')
            self._path_purelocpot = os.path.join(os.path.abspath(pure_locpot_path),'LOCPOT')
        else:
            self._path_purelocpot = os.path.abspath(pure_locpot_path)
        if 'LOCPOT' not in defect_locpot_path:
            print('defect LOCPOT not in path. appending it to path.')
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

        from pycdt.corrections.freysoldt_correction import FreysoldtCorrection

        s=FreysoldtCorrection(axis, self._dielectricconst, self._purelocpot,
            self._deflocpot, self._q, energy_cutoff=self._encut, 
            madetol=self._madetol)

        freyval=s.correction(title=title,partflag=partflag)
        if partflag in ['All','AllSplit']:
            nomtype='full correction'
        elif partflag=='pc':
            nomtype='point charge correction'
        elif partflag=='potalign':
            nomtype='potential alignment correction'
        else:
            print(partflag,' is incorrect potalign type. Must be "All","AllSplit", "pc", or "potalign".')
            return

        #if havent loaded locpots objects or positions of defect then save these for later
        # (good if you want to run freysoldt for multiple axes)
        if (type(s._purelocpot) is Locpot) and (type(self._purelocpot) is not Locpot):
            self._purelocpot=s._purelocpot
        if (type(s._deflocpot) is Locpot) and (type(self._deflocpot) is not Locpot):
            self._deflocpot=s._deflocpot
        if self._pos is None: #want them in fractional coords
            self._pos = self._purelocpot.structure.lattice.get_fractional_coords(s._pos)

        print('\n Final Freysoldt',nomtype,'value is ',freyval)

        return freyval

    def kumagai(self,title=None, partflag='All', bulk_outcar_path=None, def_outcar_path=None):
        """
        Args:
            title: set if you want to plot the atomic site averaged potential
            partflag: four options
                'pc' for just point charge correction, or
               'potalign' for just pot. align correction, or
               'All' for both, or
               'AllSplit' for individual parts split up (form [PC,potterm,full])
            bulk_outcar_path: path to Bulk OUTCAR (quicker method for performing kumagai code)
            def_outcar_path: path to defect OUTCAR
        """
        
        if bulk_outcar_path is None:
            if type(self._purelocpot) is not Locpot:
                self._purelocpot = Locpot.from_file(self._purelocpot)
            s=KumagaiCorrection(self._dieltens, 
                    self._q, self._KumagaiBulk.gamma, self._KumagaiBulk.g_sum, 
                    self._purelocpot.structure, energy_cutoff=self._encut,
                    madetol=self._madetol, 
                    bulk_locpot=self._purelocpot, defect_locpot=self._deflocpot)
        else:
            if type(self._purelocpot) is not Locpot:
                self._purelocpot = Locpot.from_file(self._purelocpot)
            if type(self._deflocpot) is not Locpot:
                self._deflocpot = Locpot.from_file(self._deflocpot)

            s=KumagaiCorrection(self._dieltens, 
                    self._q, self._KumagaiBulk.gamma, self._KumagaiBulk.g_sum, 
                    self._purelocpot.structure, energy_cutoff=self._encut, 
                    madetol=self._madetol,
                    defstructure=self._deflocpot.structure,
                    bulk_outcar=bulk_outcar_path, defect_outcar=def_outcar_path)

        if partflag in ['All','AllSplit']:
            nomtype='full correction'
            kumval=s.correction(title=title,partflag=partflag)
        elif partflag=='pc':
            nomtype='point charge correction'
            kumval=s.correction(title=title,partflag=partflag)
        elif partflag=='potalign':
            nomtype='potential alignment correction'
            kumval=s.correction(title=title,partflag=partflag)

        print('\n Final Kumagai',nomtype,'value is ',kumval)

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
                                        self._encut, lengths=lengths)

        if partflag in ['All','AllSplit']:
            nomtype='full correction'
        elif partflag=='pc':
            nomtype='point charge correction'
        elif partflag=='potalign':
            nomtype='potential alignment correction'

        sxvals=s.run_correction(print_pot_flag=print_pot_flag, partflag=partflag)

        print('\n Final Sxdefectalign ',nomtype,' correction value is ',sxvals)

        return sxvals

