"""
This module computes finite size supercell charge corrections for
defects. The methods implemtend are
1) Freysoldt correction for isotropic systems.
   Freysoldt method includes
   a) PC energy
   b) potential alignment by planar averaging.
2) Extended Freysoldt or Kumagai correction for anistropic systems.
   Kumagai method includes
   a) anisotropic PC energy
   b) potential alignment by atomic site averaging at Wigner Seitz cell
      edge
If you use the corrections implemented in this module, cite
   Freysoldt, Neugebauer, and Van de Walle,
    Phys. Status Solidi B. 248, 1067-1076 (2011) for isotropic correction
   Kumagai and Oba, Phys. Rev. B. 89, 195205 (2014) for anisotropic correction
   in addition to the pycdt paper
"""
__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import sys
import math
import os

import numpy as np
#import matplotlib #cause nersc freaks out when you uncomment these...
#matplotlib.use('agg')
#import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Locpot
from pymatgen.core.lattice import Lattice

norm = np.linalg.norm  # define globally

# Define conversion_constants
hart_to_ev = 27.2114
ang_to_bohr = 1.8897

def k_to_eV(g):
    """
    Convert a k-vector to energy [eV] via hbar*k^2/2m
    Args:
        a: Reciprocal vector (units of 1/A).

    Returns:
        (double) Energy in eV
    """
    return 3.80986 * np.dot(g,g)

def eV_to_k(energy):
    """
    Convert energy to reciprocal vector magnitude k via hbar*k^2/2m
    Args:
        a: Energy in eV.

    Returns:
        (double) Reciprocal vector magnitude (units of 1/Bohr).
    """
    return math.sqrt(energy/3.80986)*1.8897

def cleanlat(dat):
    """
    return lattice constants
    Args:
        dat: array of lattice vectors

    Returns:
        (double) Lattice constants (in same units as lattice vectors)
    """
    return norm(dat[0]), norm(dat[1]), norm(dat[2])

def genrecip(a1, a2, a3, encut):
    # latt vectors in bohr, encut=eV
    # generate reciprocal lattice vectors with value less than encut
    # define recip vectors first, (units of 1/angstrom).
    vol = np.dot(a1, np.cross(a2, a3))  # 1/bohr^3
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  # units 1/bohr
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)
    recip = []
    flag = 0
    # create list of recip space vectors that satisfy |i*b1+j*b2+k*b3|<=encut
    #start by enumerating to find max i that doesn't upset the encut condition
    tol = 0
    while flag != 1:
        if 3.80986 * ((tol * (1 / ang_to_bohr) * min(norm(b1), norm(b2), norm(b3))) ** 2) < encut:
            #added the 1.8897 factor because the energy given converts 1/A to eV but b's in 1/bohr
            tol = tol + 1
        else:
            flag = 1
    #now look though all options for recip vectors to see what vectors are less than energy val
    for i in range(-tol, tol + 1):
        for j in range(-tol, tol + 1):
            for k in range(-tol, tol + 1):
                vec = (i * b1 + j * b2 + k * b3)
                en = 3.80986 * (((1 / ang_to_bohr) * norm(vec))** 2)
                #en = 3.80986 * (((ang_to_bohr) * norm(vec))** 2)  #isn't the line above  here what we want? Depends what units of recip vector are
                if (en <= encut and en != 0):
                    recip.append([i * b1[m] + j * b2[m] + k * b3[m] for m in range(3)])
    return recip  #output is 1/bohr recip

def generate_reciprocal_vectors_squared(a1, a2, a3, encut):
    """
    Generate reciprocal vector magnitudes within the cutoff along the specied
    lattice vectors. 
    Args:
        a1: Lattice vector a (in Bohrs)
        a2: Lattice vector b (in Bohrs)
        a3: Lattice vector c (in Bohrs)
        encut: Reciprocal vector energy cutoff

    Returns:
        [[g1^2], [g2^2], ...] Square of reciprocal vectors (1/Bohr)^2 
        determined by a1, a2, a3 and whose magntidue is less than gcut^2.
    """
    vol = np.dot(a1, np.cross(a2, a3))  
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)

    # Max (i,j,k) that doesn't upset the condition |i*b1+j*b2+k*b3|<=gcut
    gcut=eV_to_k(encut)
    max_index = int(math.ceil(gcut/min(norm(b1), norm(b2), norm(b3))))
    gcut2 = gcut*gcut
    recip = []
    for i in range(-max_index, max_index+1):
        for j in range(-max_index, max_index+1):
            for k in range(-max_index, max_index+1):
                vec = (i*b1 + j*b2 + k*b3)
                vec2 = np.dot(vec,vec)
                if (vec2 <= gcut2 and vec2 != 0.0):
                    recip.append(vec2)
    return recip


def closestsites(sb,sd,pos):
    #input bulk and defect structures and get site that is nearest to the (cartesian) input position
    bulkclosesites=sb.get_sites_in_sphere(pos,5)
    bulkclosesites.sort(key=lambda x:x[1])
    defclosesites=sd.get_sites_in_sphere(pos,5)
    defclosesites.sort(key=lambda x:x[1])
    return bulkclosesites[0],defclosesites[0] #returns closest (site object, dist) for both bulk and defect


def find_defect_pos(sb,sd):
    #Will output cartesian coords of defect in bulk,defect cells.
    #If vacancy defectpos=None, if interstitial bulkpos=None, if antisite/sub then both defined
    if len(sb.sites)>len(sd.sites):
        vactype=True
        interstittype=False
    elif len(sb.sites)<len(sd.sites):
        vactype=False
        interstittype=True
    else:
        vactype=False
        interstittype=False
    sitematching=[]
    for i in sb.sites:
        blksite,defsite=closestsites(sb,sd,i.coords)
        if vactype and blksite[0].specie.symbol != defsite[0].specie.symbol:
            return blksite[0].coords, None
        elif interstittype and blksite[0].specie.symbol != defsite[0].specie.symbol:
            return None, defsite[0].coords
        elif blksite[0].specie.symbol != defsite[0].specie.symbol: #subs or antisite type
            return blksite[0].coords, defsite[0].coords
        sitematching.append([blksite[0],blksite[1],defsite[0],defsite[1]])
    if vactype: #just in case site type is same for closest site to vacancy
        sitematching.sort(key=lambda x:x[3])
        vacant=sitematching[-1]
        return vacant[0].coords, None
    elif interstittype: #just in case site type is same for closest site to interstit
        sitematching.sort(key=lambda x:x[1])
        interstit=sitematching[-1]
        return  None, interstit[2].coords
    return None,None #if you get here there is an error


def kumagai_init(s1, dieltens, sil=True):
    angset = s1.lattice.get_cartesian_coords(1)
    if not sil:
        print 'defect lattice constants are (in angstroms)' + str(cleanlat(angset))
    [a1, a2, a3] = ang_to_bohr * angset  # convert to bohr
    bohrset = [a1, a2, a3]
    vol = np.dot(a1, np.cross(a2, a3))
    if not sil:
        print 'converted to bohr for atomic units, lat consts are:' + str(cleanlat([a1, a2, a3]))
    determ = np.linalg.det(dieltens)
    invdiel = np.linalg.inv(dieltens)
    if not sil:
        print 'inv dielectric tensor is ' + str(invdiel)
    return angset, bohrset, vol, determ, invdiel


def real_sum(a1, a2, a3, r, q, dieltens, gamma, tolerance, silence=True):
    invdiel = np.linalg.inv(dieltens)
    determ = np.linalg.det(dieltens)
    realpre = q / np.sqrt(determ)
    tolerance /= hart_to_ev

    #Real space sum by converging with respect to real space vectors
    #create list of real space vectors that satisfy |i*a1+j*a2+k*a3|<=N
    Nmaxlength = 40  #tolerance for stopping real space sum convergence
    N = 2
    r_sums = []
    while N < Nmaxlength:  
        r_sum = 0.0
        if norm(r):
            for i in range(-N, N + 1):
                for j in range(-N, N + 1):
                    for k in range(-N, N + 1):
                        r_vec = i*a1 + j*a2 + k*a3 - r
                        loc_res = np.dot(r_vec, np.dot(invdiel, r_vec))
                        nmr = math.erfc(gamma * np.sqrt(loc_res))
                        dmr = np.sqrt(determ * loc_res)
                        r_sum += nmr / dmr  
        else:
            for i in range(-N, N + 1):
                for j in range(-N, N + 1):
                    for k in range(-N, N + 1):
                        if i == j == k == 0:
                            continue
                        else:
                            r_vec = i*a1 + j*a2 + k*a3 
                            loc_res = np.dot(r_vec, np.dot(invdiel, r_vec))
                            nmr = math.erfc(gamma * np.sqrt(loc_res))
                            dmr = np.sqrt(determ * loc_res)
                            r_sum += nmr / dmr  
        r_sums.append([N, realpre * r_sum])

        if N == Nmaxlength-1:
            print('Direct part could not converge with real space ' +
                   'translation tolerance of {} for gamma {}'.format(
                       Nmaxlength-1, gamma))
            return
        elif len(r_sums) > 3:
            if abs(abs(r_sums[-1][1])-abs(r_sums[-2][1])) < tolerance:
                r_sum = r_sums[-1][1]
                if not silence:
                    print("gamma is {}".format(gamma))
                    print("convergence for real summatin term occurs at " + 
                           "step {}  where real sum is {}".format(
                               N,  r_sum * hart_to_ev))
                break

        N += 1
    return r_sum


def get_g_sum_at_r(g_sum, locpot_bulk, r):
    """
    Args:
        g_sum: Reciprocal summation calculated from reciprocal_sum method
        locpot_bulk: Bulk locpot 
        r: Position relative to defect (in cartesian coords)
    Returns:
        reciprocal summ value at g_sum[i_rx,j_ry,k_rz]
    """
    #abc=locpot_bulk.structure.lattice.abc
    fraccoord=locpot_bulk.structure.lattice.get_fractional_coords(r) #neccessary fix for non-cubic cells...
    # for i in range(3):
    #     r[i]=r[i]/abc[i] #translate to fractional coords for use in getgridind
    i, j, k = getgridind(locpot_bulk, fraccoord)
    return g_sum[i, j, k]


def anisotropic_madelung_potential(locpot_bulk, g_sum, r, dieltens, q,  gamma, tolerance,
                                  silence=True):
    """
    Compute the anisotropic Madelung potential at r not equal to 0.
    For r=(0,0,0) use anisotropic_pc_energy function
    Args:
        locpot_bulk: locpot object for bulk supercell
        g_sum: Precomputed reciprocal sum for all r_vectors
        r: r vector (in cartesian coordinates) relative to defect position. 
           Non zero r is expected
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        tolerance: Tolerance parameter for numerical convergence
        gamma (float): Convergence parameter 
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    structure = locpot_bulk.structure
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens, sil=silence)

    recippartreal = q * get_g_sum_at_r(g_sum, locpot_bulk, r)
    directpart = real_sum(a1, a2, a3, r, q, dieltens, gamma, tolerance,
                          silence)

    #now add up total madelung potential part with two extra parts:
    #self interaction term
    selfint = q * np.pi / (vol * (gamma ** 2))
    if not silence:
        print ('self interaction piece is {}'.format(selfint * hart_to_ev))

    #pot = hart_to_ev * (directpart + recippartreal - selfint)
    pot = (hart_to_ev/-q) * (directpart + recippartreal - selfint) #from danny: this is to have CORRECT conversion to V from atomic units for comparison with DFT data

    return pot


def anisotropic_pc_energy(structure, g_sum, dieltens, q, gamma, tolerance,
                          silence=True):
    """
    Compute the anistropic periodic point charge interaction energy.
    Args:
        structure: Bulk structure type
        g_sum : comes from KumagaiBulkInit class
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        gamma : convergence parameter optimized in KumagaiBulkInit class
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens, sil=silence)

    #g_sum = reciprocal_sum(locpot_bulk, dieltens, q, gamma, silence=silence)

    g_part = q*g_sum[0,0,0]
    r_part = real_sum(a1, a2, a3, [0,0,0], q, dieltens, gamma, tolerance, True)

    #now add up total madelung potential part with two extra parts:
    #self interaction term
    selfint = q*np.pi / (vol * (gamma**2))
    if not silence:
        print ('reciprocal piece is {}').format(g_part * hart_to_ev)
        print ('real piece is {}').format(r_part * hart_to_ev)
        print ('self interaction piece is {}'.format(selfint * hart_to_ev))

    #surface term (only for r not at origin)
    surfterm = 2*gamma*q / np.sqrt(np.pi*determ)
    if not silence:
        print ('surface term is {}'.format(surfterm * hart_to_ev))

    pc_energy = -q*0.5*hart_to_ev*(r_part + g_part - selfint - surfterm)
    if not silence:
        print ('Final PC Energy term is then ', pc_energy, ' (eV)')

    return pc_energy


def getgridind(locpot, r, gridavg=0.0):
    """
    Computes the index of a point, r, in the locpot grid
    Args:
        locpot: Can be any volumentric data object
        r: Relative co-ordinates with respect to abc lattice vectors
        gridavg: If you want to do atomic site averaging, set gridavg to the radius of the atom at r
    Returns:
        [i,j,k]: Indices as list
    TODO: Once final, remove the getgridind inside disttrans function
    """
    abc=locpot.structure.lattice.abc
    #abclats=locpot.structure.lattice.matrix #for getting correct distances in a non-orthogonal basis?
    grdind = []
    if gridavg:
        radvals=[] #radius in terms of indices
        dxvals=[]
    for i in range(3):
        if r[i] < 0:
            while r[i]<0:
                r[i] += 1
        elif r[i] >= 1:
            while r[i]>=1:
                r[i]-=1
        r[i] *= abc[i]
        x = locpot.get_axis_grid(i)
        dx = x[1] - x[0]
        x_rprojection_delta_abs = np.absolute(x - r[i])
        ind = np.argmin(x_rprojection_delta_abs)
        if x_rprojection_delta_abs[ind] > dx*1.1: #included this tolerance to avoid numerical errors that can occur...
            print i,ind,r
            print x_rprojection_delta_abs
            raise ValueError("Input position is not within the locpot grid")
        grdind.append(ind)
        if gridavg:
            radvals.append(int(np.ceil(gridavg/dx)))
            dxvals.append(dx)
    #grdind = [grdind] # To be consistent with non-zero gridavg #this breaks next part...
    grid_dim = locpot.dim
    if gridavg:
        grdindfull=[]
        for i in range(-radvals[0], radvals[0]+1):
            for j in range(-radvals[1], radvals[1]+1):
                for k in range(-radvals[2], radvals[2]+1):
                    dtoc = [i*dxvals[0], j*dxvals[1], k*dxvals[2]] #I wonder if this breaks if the trans vectors aren't orthogonal?
                    #dtoc = i*abclats[0]/grid_dim[0]+j*abclats[1]/grid_dim[1]+k*abclats[2]/grid_dim[2] #alternative for when trans vecs arent orthogonal?
                    if norm(dtoc) < gridavg:
                        ival = (i+grdind[0]) % grid_dim[0]
                        jval = (j+grdind[1]) % grid_dim[1]
                        kval = (k+grdind[2]) % grid_dim[2] #this was a big bug! last 2 was a 1...
                        grdindfull.append((ival,jval,kval))
        grdind=grdindfull

    return grdind


def disttrans(locpot_blk, locpot_def,  silence=False):
    """
    this is function for calculating distance from defect to each atom and finding NGX grid pts at each atom
    Args:
        locpot_blk: Bulk locpot object
        locpot_def: Defect locpot object
    """
    struct = locpot_blk.structure
    defstruct = locpot_def.structure

    #Find defect location in bulk and defect cells
    blksite,defsite=find_defect_pos(struct,defstruct)
    if blksite is None and defsite is None:
        print 'Error. Not able to determine defect site...'
        return
    if not silence:
        if blksite is None:
            print 'Found defect to be Interstitial type at ',defsite
        elif defsite is None:
            print 'Found defect to be Vacancy type at ',blksite
        else:
            print 'Found defect to be antisite/subsitution type at ',blksite,' in bulk, and ',defsite,' in defect cell'

    if blksite is None:
        blksite=defsite
    elif defsite is None:
        defsite=blksite

    def_fcoord = struct.lattice.get_fractional_coords(blksite)
    def_ccoord = blksite[:]
    defcell_def_fcoord = defstruct.lattice.get_fractional_coords(defsite)
    defcell_def_ccoord = defsite[:]

    if len(struct.sites)>=len(defstruct.sites):
        sitelist=struct.sites[:]
    else: #for interstitial list
        sitelist=defstruct.sites[:]


    #better image getter since pymatgen wasnt working
    def returnclosestr(vec):
        from operator import itemgetter
        listvals=[]
        abclats=defstruct.lattice.matrix
        trylist=[-1,0,1]
        for i in trylist:
            for j in trylist:
                for k in trylist:
                    transvec=i*abclats[0]+j*abclats[1]+k*abclats[2]
                    rnew=vec-(defcell_def_ccoord+transvec)
                    listvals.append([norm(rnew),rnew,transvec])
        listvals.sort(key=itemgetter(0))
        return listvals[0] #will return [dist,r to defect, and transvec for defect]

    grid_sites = {}  # dictionary with indices keys in order of structure list
    radlist = {'Li' : 0.97, 'O' : 0.72, 'Ti' : 1.28}  #this is a specific thing for Li2TiO3...
    for i in sitelist:
        if np.array_equal(i.coords,def_ccoord): #skip defect site
            print 'This is defect! Skipping ',i
            continue
        #radius=i.specie.atomic_radius
        radius=1
        blksite,defsite=closestsites(struct,defstruct,i.coords) #returns (site object, dist) for both bulk and defect cells

        blkindex=struct.index(blksite[0])
        #defindex=defstruct.index(defsite[0])

        cart_coord = blksite[0].coords
        frac_coord = blksite[0].frac_coords
        dcart_coord = defsite[0].coords
        dfrac_coord = defsite[0].frac_coords

        closeimage=returnclosestr(dcart_coord)
        cart_reldef=closeimage[1]
        defdist=closeimage[0]

        ## this seemed broken for whatever reason
        # dist, img = struct.lattice.get_distance_and_image(def_fcoord,
        #                                                   frac_coord)
        # defdist, defimg = defstruct.lattice.get_distance_and_image(defcell_def_fcoord,
        #                                                   dfrac_coord)
        # # cart_reldef = np.dot((frac_coord + img), defstruct.lattice._matrix) \
        # #               - def_ccoord
        # cart_reldef = np.dot((dfrac_coord + img), struct.lattice._matrix) \
        #           - defcell_def_ccoord
        if abs(norm(cart_reldef) - defdist) > 0.1:
            print('image locater issue encountered for site=', blkindex,
                  '(in def cell) distance should be ', defdist, ' but calculated to be ',
                  norm(cart_reldef))
            #return #dont want to break the code here, but want flag to exist...what to do?

        rad=radlist[i.species_string]
        if blkindex in grid_sites:
            print '(WARNING) index ',blkindex,' already exists in potinddict! overwriting information. '

        grid_sites[blkindex] = {'dist': defdist,'cart': dcart_coord,
                'cart_reldef': cart_reldef,
                #'defgrid':getgridind(locpot_def,  dfrac_coord, gridavg=float(radius)),
                'defgrid':getgridind(locpot_def,  dfrac_coord, gridavg=rad),
                #'bulkgrid': getgridind(locpot_blk,  frac_coord, gridavg=float(radius)),
                'bulkgrid': getgridind(locpot_blk,  frac_coord, gridavg=rad),
                'siteobj':[i.coords,i.frac_coords,i.species_string]}

    return grid_sites


def wigner_seitz_radius(structure):
    """
    Calculate the Wigner Seitz radius for the given structure.
    Args:
        s: Either structure or VolumetricData object
    """
    try:
        lat = Lattice(structure.lattice_vectors())
    except:
        lat = Lattice(structure.structure.lattice_vectors())

    wz = lat.get_wigner_seitz_cell()
    # wz is list of WS cell face vertices
    # make list of midpoints on faces of WS cell
    dist = []
    for facet in wz:
        midpt = np.mean(np.array(facet), axis=0)
        dist.append(norm(midpt))
    wsrad = min(dist)
    return wsrad


class KumagaiBulkInit(object):
    """
    Compute the anisotropic madelung potential array from the bulk 
    locpot. This helps in evaluating the bulk supercell related
    part once to speed up the calculations.
    """
    def __init__(self, bulk_locpot_path, epsilon, encut=520, tolerance=0.0001, 
                 silence=True, optgamma=False):
        """
        Args
            bulk_locpot_path: Path of bulk locpot file OR bulk Locpot Object
            epsilon: Dielectric tensor
            encut (float): Energy cutoff for optimal gamma
            tolerance (float): Accuracy parameter
            silence (bool): If True, intermediate steps are not printed
            optgamma: if you know optimized gamma already set this to it, otherwise it will be optimized
        """
        if type(bulk_locpot_path) is Locpot:
            self.bulk_locpot = bulk_locpot_path
        else:
            if not silence:
                print 'Loading bulk Locpot'
            self.bulk_locpot = Locpot.from_file(bulk_locpot_path)
        self.epsilon = epsilon
        self.encut = encut
        self.tolerance = tolerance
        self.silence = silence
        if not optgamma:
            self.gamma = self.find_optimal_gamma()
        else:
            self.gamma = optgamma
        self.g_sum = self.reciprocal_sum()

    def find_optimal_gamma(self):
        """
        Find optimal gamma by evaluating the brute force reciprocal
        summation and seeing when the values are on the order of 1
        this calculation is the anisotropic Madelung potential at r = (0, 0, 0)
        Args:
            structure: Bulk supercell structure
            dieltens: dielectric tensor
            madetol: Tolerance parameter for numerical convergence
            silence (bool): Verbosity flag. If False, messages are printed.
        """
        structure = self.bulk_locpot.structure
        dieltens = self.epsilon
        silence = self.silence
        madetol = self.tolerance

        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                structure, dieltens, sil=silence)
        optgam = None
        #do brute force recip summation
        def get_recippart(encut, gamma):
            recip = genrecip(a1, a2, a3, encut)
            recippart = 0.0
            for rec in recip:
                Gdotdiel = np.dot(rec, np.dot(dieltens, rec))
                summand = math.exp(-Gdotdiel / (4 * (gamma ** 2))) / Gdotdiel
                recippart += summand
            recippart *= 4*np.pi/vol # q is taken as 1
            return recippart, 0.0, len(recip)

        def do_summation(gamma):
            #do recip sum until it is bigger than 1eV
            #First do Recip space sum convergence with respect to encut for this gamma
            encut = 20  #start with small encut for expediency
            recippartreal1, recippartimag1, len_recip = get_recippart(encut, gamma)
            encut += 10
            recippartreal, recippartimag, len_recip = get_recippart(encut, gamma)
            converge = [recippartreal1, recippartreal]
            while abs(abs(converge[0]) - abs(converge[1])) * hart_to_ev > madetol:
                encut += 10
                recippartreal, recippartimag, len_recip = get_recippart(encut, gamma)
                converge.reverse()
                converge[1] = recippartreal
                if encut > self.encut:
                    raise ValueError(
                            'Optimal gamma not found at {} eV cutoff'.format(
                                self.encut))

            if abs(recippartimag) * hart_to_ev > madetol:
                if not silence:
                    print("Problem with convergence of imaginary part of recip sum."),
                    print("imag sum value is {} (eV)".format( recippartimag * hart_to_ev))
                return None, None
            if not silence:
                print('recip sum converged to {} (eV) at encut= {}'.format(
                            recippartreal * hart_to_ev, encut))
                print('Number of reciprocal vectors is {}'.format(len_recip))
            if (abs(converge[1]) * hart_to_ev < 1 and not optgam):
                if not silence:
                    print('Reciprocal summation value is less than 1 eV.')
                    print('This might lead to errors in the reciprocal summation.')
                    print('Changing gamma now.')
                return None, 'Try Again'

            return recippartreal, gamma

        #start with gamma s.t. gamma*L=5 (some paper said this is optimal)
        #optimizing gamma for the recip sum to improve convergence of calculation
        gamma = 5./(vol ** (1/3.))
        optimal_gamma_found = False
        while not optimal_gamma_found:
            recippartreal, optgamma = do_summation(gamma)
            if optgamma == gamma:
                print('optimized gamma found to be ', optgamma)
                optimal_gamma_found = True
            elif 'Try Again' in optgamma:
                gamma *= 1.5
            else:
                print('Had problem in gamma optimization process.')
                return None

            if gamma > 50:  #kind of an arbitrary number for cut off...
                print('WARNING. could not optimize gamma before gamma =', 50)
                return None

        return optgamma 

    def reciprocal_sum(self):
        """
        Compute the reciprocal summation in the anisotropic Madelung 
        potential.
        Args:
            locpot_bulk: Bulk calculation potential file
            dieltens: dielectric tensor
            silence (bool): Verbosity flag. If False, messages are printed.
            optgam (float): Convergence parameter (optional)
        TODO: Get the input to fft cut by half by using rfft instead of fft
        """
        locpot_bulk = self.bulk_locpot
        dieltens = self.epsilon 
        gamma = self.gamma
        silence = self.silence
        structure = locpot_bulk.structure

        if not silence:
            print 'calculating the reciprocal summation in Madeling potential'
        over_atob = 1.0/ang_to_bohr
        atob3=ang_to_bohr**3

        latt = structure.lattice
        vol = latt.volume*atob3 # in Bohr^3

        reci_latt = latt.reciprocal_lattice
        [b1, b2, b3] = reci_latt.get_cartesian_coords(1)
        b1 = np.array(b1)*over_atob # In 1/Bohr
        b2 = np.array(b2)*over_atob
        b3 = np.array(b3)*over_atob

        ndim = locpot_bulk.dim
        nx, ny, nz = ndim
        ind1 = np.arange(nx)
        for i in range(nx/2, nx):
            ind1[i] = i - nx
        ind2 = np.arange(ny)
        for i in range(ny/2, ny):
            ind2[i] = i - ny
        ind3 = np.arange(nz)
        for i in range(nz/2, nz):
            ind3[i] = i - nz

        g_array = np.zeros(ndim, np.dtype('c16'))
        gamm2 = 4*(gamma**2)
        for i in ind1:
            for j in ind2:
                for k in ind3:
                    g = i*b1 + j*b2 + k*b3
                    g_eps_g = np.dot(g, np.dot(dieltens, g))
                    if i == j == k == 0:
                        continue
                    else:
                        g_array[i,j,k] = math.exp(-g_eps_g/gamm2)/g_eps_g

        r_array = np.fft.fftn(g_array)
        over_vol = 4*np.pi/vol # Multiply with q later
        r_array *= over_vol
        r_arr_real = np.real(r_array)
        r_arr_imag = np.imag(r_array)

        if not silence:
            max_imag = r_arr_imag.max()
            print 'Max imaginary part found to be ', max_imag

        return r_arr_real

class KumagaiCorrection(object):
    """
    Class that implements the extended freysoldt correction developed 
    by Kumagai.
    """
    def __init__(self, dielectric_tensor, bulk_locpot,
            defect_locpot, q, gamma=None, g_sum=None,
            energy_cutoff=520, madetol=0.0001, silence=False, lengths=None):
        """
        Args:
            dielectric_tensor: Macroscopic dielectric tensor 
                 Include ionic also if defect is relaxed, othewise ion clamped.
                 Can be a matrix, array or scalar.
            bulk_locpot: Bulk Locpot file path OR Bulk Locpot Object
            defect_locpot: Defect Locpot file path OR defect Locpot Object
            q (int): Charge associated with the defect.
            gamma:  Convergence parameter. Obtained from KumagaiBulkPart
            g_sum: value that is dependent on the Bulk only. Obtained from KumagaiBulkPart
            energy_cutoff: Energy for plane wave cutoff (in eV).
                 If not given, Materials Project default 520 eV is used.
            madetol: Tolerance (double or float)
            silence: Flag for disabling/enabling  messages (Bool)
        """
        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self.dieltens = np.diag(
                np.ones(3) * dielectric_tensor)  #this would make kumagai correction pointless
        else:
            self.dieltens = np.array(dielectric_tensor)  #full dielectric tensor

        self.locpot_blk = bulk_locpot  #location of purelocpot OR locpot object
        self.locpot_def = defect_locpot  #location of defectlocpot OR locpotobject
        self.madetol = madetol #tolerance for convergence of energy terms in eV
        self.q = q  #charge of defect (not of the homogen. background)
        self.encut = energy_cutoff  #encut (eV) for calculation
        self.silence = silence  #for silencing printflags
        if not gamma:
            bi=KumagaiBulkInit(bulk_locpot, self.dieltens, encut=energy_cutoff, tolerance=madetol,
                 silence=silence)
            self.gamma=bi.gamma
            self.g_sum=bi.g_sum #even if you already put in g_sum, if gamma changed then you need to recompute g_sum
        else:
            self.gamma = gamma
            if g_sum is None:
                bi=KumagaiBulkInit(bulk_locpot, self.dieltens, encut=energy_cutoff, tolerance=madetol,
                        silence=silence,optgamma=gamma)
                self.g_sum = bi.g_sum
            else:
                self.g_sum = g_sum
        self.lengths=lengths #mainly used for quick plotting when structure isnt available...

    def correction(self, title=None, partflag='All'):
        """
        Computes the extended Freysoldt correction for anistropic systems 
        developed by Y. Kumagai and F. Oba (Ref: PRB 89, 195205 (2014)

        If you want a plot of potential averaging process set title to name of defect
        vb and vd are preloaded locpot objects for speeding this up.

        part flag can be 'pc' for just point charge correction,
               'potalign' for just potalign correction, 'All' for one combined correction,
              or 'AllSplit' for correction in form [PC,potterm,full]
        """
        if not self.silence:
            print 'This is Kumagai Correction.'
        if not self.q:
            if partflag=='AllSplit':
                return [0.,0.,0.]
            else:
                return 0.0

        if partflag!='potalign':
            energy_pc = self.pc()

        if partflag!='pc':
            if not self.silence:
                #print '\nPC calc done, correction =', round(energy_pc, 4),
                print 'Now run potential alignment script'
            #structure = self.locpot_blk.structure
            potalign = self.potalign(title=title)

        if not self.silence:
            print '\n\nKumagai Correction details:'
            if partflag!='potalign':
                print 'PCenergy = ', round(energy_pc, 5)
            if partflag!='pc':
                print 'potential alignment = ', round(potalign, 5)
            if partflag in ['All','AllSplit']:
                print 'TOTAL Kumagai correction = ', round(energy_pc - potalign, 5)

        if partflag=='pc':
            return round(energy_pc,5)
        elif partflag=='potalign':
            return round(potalign,5)
        elif partflag=='All':
            return round(energy_pc-potalign,5)
        else:
            return [round(energy_pc,5),round(potalign,5),round(energy_pc-potalign,5)]


    def pc(self):
        if not self.silence:
            print '\nrun Kumagai PC calculation'

        if not type(self.locpot_blk) is Locpot:
            if not self.silence:
                print 'Load bulk Locpot'
            self.locpot_blk=Locpot.from_file(self.locpot_blk)

        energy_pc = anisotropic_pc_energy(
                self.locpot_blk.structure, self.g_sum, self.dieltens, self.q,
                self.gamma, self.madetol, silence=self.silence)  #returns PCenergy in eV

        if not self.silence:
            print 'PC energy determined to be ', energy_pc, ' eV (', \
                    energy_pc/hart_to_ev, ' Hartree)'  #27.2114 eV/1 hartree

        return energy_pc  #PC energy in eV

    def potalign(self, title=None):
        """
        Potential alignment for Kumagai method
        Args:
            gamma: Convergence parameter
            locpot_blk: Bulk locpot object
            locpot_def: Defect locpot object
            title: Title for the plot. None will not generate the plot
        """
        if not self.silence:
            print ('\nrun Kumagai potential calculation (atomic site averaging)')

        if not type(self.locpot_blk) is Locpot:
            if not self.silence:
                print 'Load bulk Locpot'
            self.locpot_blk=Locpot.from_file(self.locpot_blk)

        if not type(self.locpot_def) is Locpot:
            if not self.silence:
                print 'Load defect Locpot'
            self.locpot_def=Locpot.from_file(self.locpot_def)

        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                self.locpot_blk.structure, self.dieltens, sil=self.silence)

        #this is to calculate distance matrix for plotting
        potinddict = disttrans(self.locpot_blk, self.locpot_def,silence=self.silence)

        #wsrad = wigner_seitz_radius(self.locpot_blk.structure) #fudge factor...not actually using wigner_seitz_radius
        wsrad = max(norm(a1),norm(a2),norm(a3))/2. #fudge factor...

        if not self.silence:
            print ('wsrad', wsrad)

        for i in potinddict.keys():
            if potinddict[i]['dist'] > wsrad:
                potinddict[i]['OutsideWS'] = True
            else:
                potinddict[i]['OutsideWS'] = False

        #calculate potential value at each atomic site
        puredat = self.locpot_blk.data["total"]
        defdat = self.locpot_def.data["total"]
        jup = 0
        for i in potinddict.keys():
            jup += 1
            if (not title and not potinddict[i]['OutsideWS']):
                #dont need to calculate inside WS if not printing plot
                continue
            if not self.silence:
                print '-------------------------------------'
                print "calculate alignment potential data for atom " + str(i) \
                      + " (dist=" + str(potinddict[i]['dist']) + ")"
            #dx, dy, dz = potinddict[i]['defgrid'] #really should change indices based on whether the atom moved or not...
            #dx, dy, dz = potinddict[i]['bulkgrid']
            #bx, by, bz = potinddict[i]['bulkgrid']
            #v_qb = defdat[dx][dy][dz] - puredat[bx][by][bz]
            #do averaging routine
            bulkvals=[]
            defvals=[]
            for u,v,w in potinddict[i]['bulkgrid']:
                bulkvals.append(puredat[u][v][w])
            for u,v,w in potinddict[i]['defgrid']:
                defvals.append(defdat[u][v][w])
            # print 'defdat val = ',np.mean(defvals)
            # print 'puredat val = ',np.mean(bulkvals)
            print 'defdat val = ',np.mean(defvals)
            print 'puredat val = ',np.mean(bulkvals)
            v_qb = np.mean(defvals) - np.mean(bulkvals)
            cart_reldef = potinddict[i]['cart_reldef'] #should change this to averaging pure
            # and def within a range then subtract

            v_pc = anisotropic_madelung_potential(self.locpot_blk, self.g_sum,
                    cart_reldef, self.dieltens, self.q, self.gamma,
                    self.madetol, silence=True)

            #this is fudge factor that I cant figure out... (Danny 3/21/16)
            #v_pc/=-self.q #fixed this is the function...
            v_qb*=-1 #change sign convention of electron charge

            potinddict[i]['Vpc'] = v_pc
            potinddict[i]['Vqb'] = v_qb
            if not self.silence:
                print 'Has anisotropic madelung potential =', v_pc
                print 'DFT bulk/defect difference = ', v_qb
                print 'atoms left to calculate = ' + str(len(potinddict.keys()) - jup) #this is broken if you arent plotting
        if not self.silence:
            print '--------------------------------------'

        #now parse and plot if neccessary
        if title:  #to make shading region prettier
            fullspecset = self.locpot_blk.structure.species
            specset = list(set(fullspecset))
            shade, forplot = {}, {}
            for i in specset:
                shade[i.symbol] = {'r': [], 'Vpc': [], 'Vqb': []}
                forplot[i.symbol] = {'r': [], 'Vpc': [], 'Vqb': [],'sites':[]}

        forcorrection = []
        for i in potinddict.keys():
            if (not title and not potinddict[i]['OutsideWS']):
                continue
            if potinddict[i]['OutsideWS']:
                forcorrection.append(potinddict[i]['Vqb'] - potinddict[i]['Vpc'])
                if title:
                    elt = fullspecset[i].symbol
                    shade[elt]['r'].append(potinddict[i]['dist'])
                    shade[elt]['Vpc'].append(potinddict[i]['Vpc'])
                    shade[elt]['Vqb'].append(potinddict[i]['Vqb'])
            if title:
                elt = fullspecset[i].symbol
                forplot[elt]['r'].append(potinddict[i]['dist'])
                forplot[elt]['Vpc'].append(potinddict[i]['Vpc'])
                forplot[elt]['Vqb'].append(potinddict[i]['Vqb'])
                forplot[elt]['sites'].append(potinddict[i]['siteobj'])

        potalign = np.mean(forcorrection)  #I think I might need to multiply this by the charge now since I am doing weird fudge factor?
        #note I am already returning a potalign correction that is multiplied by q...

        if title:
            if title!='written':
                import matplotlib.pyplot as plt
                plt.figure(2)
                plt.clf()
                collis = ['b', 'g', 'c', 'm', 'y', 'w', 'k']
                ylis = []
                rlis = []
                for i in range(len(forplot.keys())):
                    inkey = forplot.keys()[i]
                    for k in forplot[inkey]['r']:
                        rlis.append(k)
                    for k in ['Vqb', 'Vpc']:
                        for u in forplot[inkey][k]:
                            ylis.append(u)
                    plt.plot(forplot[inkey]['r'], forplot[inkey]['Vqb'], color=collis[i], marker='^', linestyle='None',
                             label=str(inkey) + ': V_{q/b}')
                    plt.plot(forplot[inkey]['r'], forplot[inkey]['Vpc'], color=collis[i], marker='o', linestyle='None',
                             label=str(inkey) + ': Vpc')
                full = []
                for i in forplot.keys():
                    for k in range(len(forplot[i]['Vpc'])):
                        full.append([forplot[i]['r'][k], forplot[i]['Vqb'][k] - forplot[i]['Vpc'][k]])
                realfull = sorted(full, key=lambda x: x[0])
                r, y = [], []
                for i in realfull:
                    r.append(i[0])
                    y.append(i[1])
                plt.plot(r, y, color=collis[-1], marker='x', linestyle='None', label='V_{q/b} - Vpc')
                plt.xlabel('Distance from defect (A)')
                plt.ylabel('Potential (V)')
                x = np.arange(wsrad, max(self.locpot_blk.structure.lattice.abc), 0.01)
                plt.fill_between(x, min(ylis) - 1, max(ylis) + 1, facecolor='red', alpha=0.15, label='sampling region')
                plt.axhline(y=potalign, linewidth=0.5, color='red', label='pot. align.')
                plt.legend(loc=8)
                plt.axhline(y=0, linewidth=0.2, color='black')
                plt.ylim([min(ylis) - .5, max(ylis) + .5])
                plt.xlim([0, max(rlis) + 3])

                plt.title(str(title) + ' atomic site potential plot')
                #plt.show()
                plt.savefig(str(title) + 'kumagaisiteavgPlot.png')
            else:
                from monty.serialization import dumpfn
                from monty.json import MontyEncoder
                forplot['EXTRA']={'wsrad':wsrad,'potalign':potalign}
                fname='KumagaiData.json'
                #fname='KumagaiData.dat'
                dumpfn(forplot, fname, cls=MontyEncoder)
                # with open(fname,'w') as f:
                #     f.write(str(forplot))

        if self.silence == False:
            print 'Atomic site method potential alignment term is ' + str(np.mean(forcorrection))
            print 'this yields total (q*align) Kumagai potential correction energy of ' \
                  + str(self.q * np.mean(forcorrection)) + ' (eV) '

        return self.q * np.mean(forcorrection)

    def plot_from_datfile(self,name='KumagaiData.json',title='default'):
        """
        Takes data file called 'name' and does plotting.
        Good for later plotting of locpot data after running run_correction()

        """
        if type(self.locpot_blk) is not Locpot and not self.lengths:
            self.locpot_blk=Locpot.from_file(self.locpot_blk) #do this for plotting axes sizes...

        from monty.serialization import loadfn
        from monty.json import MontyDecoder
        forplot=loadfn(name, cls=MontyDecoder)

        import matplotlib.pyplot as plt
        plt.figure()
        plt.clf()
        collis = ['b', 'g', 'c', 'm', 'y', 'w', 'k']
        ylis = []
        rlis = []
        for i in range(len(forplot.keys())):
            inkey = forplot.keys()[i]
            if inkey=='EXTRA':
                continue
            for k in forplot[inkey]['r']:
                rlis.append(k)
            for k in ['Vqb', 'Vpc']:
                for u in forplot[inkey][k]:
                    ylis.append(u)
            plt.plot(forplot[inkey]['r'], forplot[inkey]['Vqb'], color=collis[i], marker='^', linestyle='None',
                     label=str(inkey) + ': V_{q/b}')
            plt.plot(forplot[inkey]['r'], forplot[inkey]['Vpc'], color=collis[i], marker='o', linestyle='None',
                     label=str(inkey) + ': Vpc')
        full = []
        for i in forplot.keys():
            if i=='EXTRA':
                continue
            for k in range(len(forplot[i]['Vpc'])):
                full.append([forplot[i]['r'][k], forplot[i]['Vqb'][k] - forplot[i]['Vpc'][k]])
        realfull = sorted(full, key=lambda x: x[0])
        r, y = [], []
        for i in realfull:
            r.append(i[0])
            y.append(i[1])
        wsrad=forplot['EXTRA']['wsrad']
        potalign=forplot['EXTRA']['potalign']
        plt.plot(r, y, color=collis[-1], marker='x', linestyle='None', label='V_{q/b} - Vpc')
        plt.xlabel('Distance from defect (A)')
        plt.ylabel('Potential (V)')
        try:
            x = np.arange(wsrad, max(self.locpot_blk.structure.lattice.abc), 0.01)
        except:
            x = np.arange(wsrad, max(self.lengths), 0.01)
        plt.fill_between(x, min(ylis) - 1, max(ylis) + 1, facecolor='red', alpha=0.15, label='sampling region')
        plt.axhline(y=potalign, linewidth=0.5, color='red', label='pot. align.')
        plt.legend(loc=8)
        plt.axhline(y=0, linewidth=0.2, color='black')
        plt.ylim([min(ylis) - .5, max(ylis) + .5])
        plt.xlim([0, max(rlis) + 3])

        plt.title(str(title) + ' atomic site potential plot')
        #plt.show()
        plt.savefig(str(title) + 'kumagaisiteavgPlot.png')



if __name__ == '__main__':
    s = KumagaiCorrection(18.099,
        '../../../../Gavacm3testMachgcorr/LOCPOT_vref','../../../../Gavacm3testMachgcorr/LOCPOT_vdef',-3,
        gamma=4.160061)
    s.correction(title='Testing')
