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

import os, shutil
import numpy as np
from pymatgen.io.vasp.outputs import Locpot
from pycdt.corrections.kumagai_correction import KumagaiBulkInit, KumagaiCorrection
from pycdt.corrections.freysoldt_correction import FreysoldtCorrection
from pycdt.corrections.sxdefect_correction import FreysoldtCorrection as SXD

def get_correction_freysoldt(defect, bulk_entry, epsilon, title = None,
                             partflag='All', averaging=False, defpos = None):
    """
    Function to compute the isotropic freysoldt correction for each defect.
    Args:
        defect: ComputedDefect object
        bulk_entry: ComputedStructureEntry corresponding to bulk OR bulk Locpot Object
        epsilon: dielectric constant
        title: decides whether to plot electrostatic potential plots or not...
            if None, no plot is printed, if a string,
            then the plot will include that string in it's label
        partflag: four options
                'pc' for just point charge correction, or
               'potalign' for just potalign correction, or
               'All' for both, or
               'AllSplit' for individual parts split up (form is [PC, potterm, full])
        averaging (bool): method for taking 3-axis average of freysoldt.
        defpos: (if known) defect position as a pymatgen Site object within bulk supercell
    """
    if partflag in ['All','AllSplit']:
        nomtype='full correction'
    elif partflag=='pc':
        nomtype='point charge correction'
    elif partflag=='potalign':
        nomtype='potential alignment correction'
    else:
        print(partflag,' is incorrect potalign type. Must be "All","AllSplit", "pc", or "potalign".')
        return

    if type(bulk_entry) is Locpot:
        locpot_blk = bulk_entry
    else:
        locpot_blk = bulk_entry.data['locpot_path']
    locpot_path_def = defect.entry.data['locpot_path']
    dpat, tmplocpotnom = os.path.split(locpot_path_def)
    charge = defect.charge
    encut = defect.entry.data['encut']

    if not charge:
        print('charge is zero so charge correction is zero')
        return (0.,bulk_entry)

    if not averaging:
        #single freysoldt correction along x-axis
        corr_meth = FreysoldtCorrection(0, epsilon, locpot_blk,
                                locpot_path_def, charge, energy_cutoff = encut, defect_position=defpos)
        corr_val = corr_meth.correction(title=title,partflag=partflag)
    else:
        #averagefreysoldtcorrection over threeaxes
        avgcorr = []
        fullcorrset = []
        locpot_def = locpot_path_def #start with path to defect
        for ax in range(3):
            corr_meth = FreysoldtCorrection(ax, epsilon, locpot_blk,
                                    locpot_def, charge, energy_cutoff = encut, defect_position=defpos)
            valset = corr_meth.correction(title=title+'ax'+str(ax+1), partflag=partflag)

            if partflag=='AllSplit':
                avgcorr.append(valset[1])
            else:
                avgcorr.append(valset)
            fullcorrset.append(valset)
            if title:
                homepat = os.path.abspath('.')
                src = os.path.join(homepat, title+'ax'+str(ax+1)+'FreyplnravgPlot.pdf')
                dst = os.path.join(dpat, title+'ax'+str(ax+1)+'FreyplnravgPlot.pdf')
                shutil.move(src, dst)
            locpot_blk = corr_meth._purelocpot #prevent reloading of locpots to save time
            locpot_def = corr_meth._deflocpot
        if partflag=='AllSplit':
            corr_val = [valset[0], np.mean(avgcorr), valset[0]+np.mean(avgcorr), fullcorrset]
        else:
            corr_val = np.mean(avgcorr)

    if partflag=='AllSplit':
        freyval = corr_val[2]
    else:
        freyval = corr_val

    print('\n Final Freysoldt',nomtype,'value is ',freyval)

    return (corr_val,corr_meth._purelocpot)


def get_correction_kumagai(defect, path_blk, bulk_init, bulk_locpot=None,
                           title=None, defpos = None):
    """
    Function to compute the Kumagai correction for each defect (modified freysoldt for anisotropic dielectric).
    NOTE that bulk_init class must be pre-instantiated to use this function
    Args:
        defect: ComputedDefect object
        path_blk: location to Bulk folder (only needed if bulk_locpot not supplied)
        bulk_init: KumagainBulkInit class object
            note this contains the dielectric tensor to be used...
        bulk_locpot: BulkLocpot object 
                (if already loaded, otherwise will load from path_blk)
        title: decides whether to plot electrostatic potential plots or not
            if None, no plot is printed, if a string, then the plot will include that string in it's label
        defpos: (if known) defect position as a pymatgen Site object within bulk supercell
    """
    epsilon = bulk_init.epsilon
    charge = defect.charge
    encut = defect.entry.data['encut']

    if not charge:
        print('charge is zero so charge correction is zero')
        return 0.

    outcar_path_blk = os.path.join(path_blk,'OUTCAR')
    locpot_path_def = defect.entry.data['locpot_path']
    dpat,dloc = os.path.split(locpot_path_def)
    outcar_path_def = os.path.join(dpat,'OUTCAR')

    if os.path.exists(outcar_path_blk) and os.path.exists(outcar_path_def):
        s = KumagaiCorrection(epsilon, charge, bulk_init.gamma, 
                bulk_init.g_sum, bulk_init.structure, defect.entry.structure, 
                energy_cutoff=encut, madetol=0.0001, 
                bulk_outcar=outcar_path_blk, defect_outcar=outcar_path_def, defect_position=defpos)
    else:
        if not bulk_locpot:
            bulk_locpot = Locpot.from_file(path_blk)
        s = KumagaiCorrection(epsilon, charge, bulk_init.gamma,
                bulk_init.g_sum, bulk_init.structure, defect.entry.structure,
                energy_cutoff=encut, madetol=0.0001, 
                bulk_locpot=bulk_locpot, defect_locpot=locpot_path_def, defect_position=defpos)

    kumval = s.correction(title=title, partflag='All')
    print('\n Kumagai Correction value is ', kumval)
    return kumval


def get_correction_sxdefect(path_def, path_blk, epsilon, pos, charge, title=None,
                            lengths=None, partflag='All', encut=520):
        """
            NOTE FROM DEVELOPERS:
            This is not unit tested and will not be maintained past 12/15/17.
            Code remaining here to allow for existing users to keep using it.

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
        #TODO: test this function

        if partflag in ['All','AllSplit']:
            nomtype='full correction'
        elif partflag=='pc':
            nomtype='point charge correction'
        elif partflag=='potalign':
            nomtype='potential alignment correction'
        else:
            print(partflag,' is incorrect potalign type. Must be "All","AllSplit", "pc", or "potalign".')
            return

        s = SXD(path_blk, path_def, charge, epsilon, pos, encut, lengths=lengths)

        if title:
            print_flag = 'plotfull'
        else:
            print_flag = 'none'
        sxvals=s.run_correction(print_pot_flag=print_flag, partflag=partflag)

        print('\n Final Sxdefectalign ',nomtype,' correction value is ',sxvals)

        return sxvals

