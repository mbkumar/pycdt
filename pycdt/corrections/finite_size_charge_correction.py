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

from pymatgen.io.vaspio.vasp_output import Locpot
from pymatgen.core.lattice import Lattice
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import math

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


class QModel():
    """
    Model for the defect charge distribution.
    A combination of exponential tail and gaussian distribution is used
    (see Freysoldt (2011), DOI: 10.1002/pssb.201046289 )
    q_model(r) = q [x exp(-r/gamma) + (1-x) exp(-r^2/beta^2)]
            without normalization constants
    By default, gaussian distribution with 1 Bohr width is assumed.
    If defect charge is more delocalized, exponential tail is needed.
    """
    def __init__(self, beta=1.0, expnorm=0.0, gamma=1.0):
        """
        Args:
            beta: Gaussian decay constant. Default value is 1 Bohr.
                  When delocalized (eg. diamond), 2 Bohr is more appropriate.
            expnorm: Weight for the exponential tail in the range of [0-1].
                     Default is 0.0 indicating no tail .
                     For delocalized charges ideal value is around 0.54-0.6.
            gamma: Exponential decay constant
        """
        self.beta2 = beta * beta
        self.x = expnorm
        self.gamma2 = gamma * gamma
        if expnorm and not gamma:
            raise ValueError("Please supply exponential decay constant.")

    def rho_rec(self, g2):
        """
        Reciprocal space model charge value
        for input squared reciprocal vector.
        Args:
            g2: Square of reciprocal vector

        Returns:
            Charge density at the reciprocal vector magnitude
        """
        return self.x/np.sqrt(1+self.gamma2*g2) + \
               (1-self.x)*np.exp(-0.25*self.beta2*g2)

    def rho_rec_limit0(self):
        """
        Reciprocal space model charge value
        close to reciprocal vector 0 .
        rho_rec(g->0) -> 1 + rho_rec_limit0 * g^2
        """
        return -2*self.gamma2*self.x - 0.25*self.beta2*(1-self.x)


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

def reciprocal_sum(locpot_bulk, dieltens, q, gamma, silence=False):
    """
    Compute the reciprocal summation in the anisotropic Madelung 
    potential.
    Args:
        locpot_bulk: Bulk calculation potential file
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        silence (bool): Verbosity flag. If False, messages are printed.
        optgam (float): Convergence parameter (optional)
    TODO: Get the input to fft cut by half by using rfft instead of fft
    """
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
    over_vol = 4*np.pi*q/vol
    r_array *= over_vol
    r_arr_real = np.real(r_array)
    r_arr_imag = np.imag(r_array)

    #ind1 = np.arange(nx/2+1)
    #for i in range(nx/2, nx):
    #    ind1[i] = i - nx
    #ind2 = np.arange(ny/2+1)
    #for i in range(ny/2, ny):
    #    ind2[i] = i - ny
    #ind3 = np.arange(nz/2+1)
    #for i in range(nz/2, nz):
    #    ind3[i] = i - nz

    #g_array = np.zeros([nx/2+1,ny/2+1,nz/2+1], np.dtype('c16'))
    #for i in ind1:
    #    for j in ind2:
    #        for k in ind3:
    #            g = (i*b1 + j*b2 + k*b3)
    #            g_eps_g = np.dot(g, np.dot(dieltens, g))
    #            if i == j == k == 0:
    #                continue
    #            else:
    #                g_array[i,j,k] = math.exp(-g_eps_g/gamm2)/g_eps_g

    #r_array1 = np.fft.irfftn(g_array)
    #r_array1 *= over_vol
    #r_arr_real1 = np.real(r_array1)
    #r_arr_imag1 = np.imag(r_array1)
    
    #print ('r_arr_vals', r_arr_real[0,0,0], r_arr_real1[0,0,0])

    if not silence:
        max_imag = r_arr_imag.max()
        print 'Max imaginary part found to be ', max_imag

    #import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(np.arange(nx), np.arange(ny), r_arr_real[:,:,0])
    #plt.savefig('3dplot.png')
    return r_arr_real

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

def find_optimal_gamma(structure, dieltens, q, madetol, silence=True):
    """
    Find optimal gamma by evaluating the brute force reciprocal
    summation and seeing when the values are on the order of 1
    this calculation is the anisotropic Madelung potential at r = (0, 0, 0)
    Args:
        structure: Bulk supercell structure
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        madetol: Tolerance parameter for numerical convergence
        silence (bool): Verbosity flag. If False, messages are printed.
    """
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
        recippart *= 4*np.pi*q/vol
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
            if encut > 300:
                print('Problem, recip sum not converged at encut = 300eV')
                return

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

def get_g_sum_at_r(g_sum, locpot_bulk, r):
    """
    Args:
        g_sum: Reciprocal summation calculated from reciprocal_sum method
        locpot_bulk: Bulk locpot 
        r: Position relative to defect (in cartesian coords)
    Returns:
        reciprocal summ value at g_sum[i_rx,j_ry,k_rz]
    """
    abc=locpot_bulk.structure.lattice.abc
    for i in range(3):
        r[i]=r[i]/abc[i] #translate to fractional coords for use in getgridind1
    i, j, k = getgridind1(locpot_bulk, r)
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

    recippartreal = get_g_sum_at_r(g_sum, locpot_bulk, r)
    directpart = real_sum(a1, a2, a3, r, q, dieltens, gamma, tolerance,
                          silence)

    #now add up total madelung potential part with two extra parts:
    #self interaction term
    selfint = q * np.pi / (vol * (gamma ** 2))
    if not silence:
        print ('self interaction piece is {}'.format(selfint * hart_to_ev))

    pot = hart_to_ev * (directpart + recippartreal - selfint)

    return pot

def anisotropic_pc_energy(structure, g_sum, dieltens, q, gamma, tolerance,
                          silence=True):
    """
    Compute the anistropic periodic point charge interaction energy.
    Args:
        structure: Bulk structure type
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens, sil=silence)

    #g_sum = reciprocal_sum(locpot_bulk, dieltens, q, gamma, silence=silence)

    g_part = g_sum[0,0,0]
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

def getgridind1(locpot, r, gridavg=False):
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
    grdind = []
    if gridavg:
        radvals=[] #radius in terms of indices
        dxvals=[]
    for i in range(3):
        if r[i]<0:
            r[i]+=1
        r[i]*=abc[i]
        x = locpot.get_axis_grid(i)
        dx = x[1] - x[0]
        x_rprojection_delta_abs = np.absolute(x - r[i])
        ind = np.argmin(x_rprojection_delta_abs)
        if x_rprojection_delta_abs[ind] > dx:
            print i,ind,r
            print x_rprojection_delta_abs
            raise ValueError("Input position is not within the locpot grid")
        grdind.append(ind)
        if gridavg:
            radvals.append(int(gridavg/dx))
            dxvals.append(dx)
    if gridavg:
        grdindfull=[]
        for i in range(-radvals[0],radvals[0]+1):
            for j in range(-radvals[1],radvals[1]+1):
                for k in range(-radvals[2],radvals[2]+1):
                    dtoc=[float(i)*dxvals[0],float(j)*dxvals[1],float(k)*dxvals[2]]
                    if norm(dtoc)<gridavg:
                        ival=i+grdind[0]
                        jval=j+grdind[1]
                        kval=k+grdind[2]
                        if ival<0:
                            ival+=len(locpot.get_axis_grid(0))
                        if jval<0:
                            jval+=len(locpot.get_axis_grid(1))
                        if kval<0:
                            kval+=len(locpot.get_axis_grid(2))
                        grdindfull.append((ival,jval,kval))
        grdind=grdindfull

    return grdind

def disttrans1(locpot_blk, locpot_def, r_def, coords_are_cartesian=False):
    """
    this is function for calculating distance to each atom and finding NGX grid pts at each atom
    Args:
        locpot_blk: Bulk locpot object
        locpot_def: Defect locpot object
        r_def: Coordinates of defect (defaulted fractional)
        coords_are_cartesian (bool): Default is False
    """
    struct = locpot_blk.structure
    defstruct = locpot_def.structure

    #THIS IS ALL DANNY ADDED STUFF FOR FINDING DEFECT...NEED SMARTER WAY?
    #defect_index is where defect is (0 defined), if site >= to c then subtract one to get defect site...specific for vacancies
    from pymatgen.analysis.structure_matcher import StructureMatcher as SM
    smatch = SM(primitive_cell=False, allow_subset=True)
    matches = smatch.get_transformation(struct, defstruct)
    if not matches:
        print 'Could not match bulk and defect structures! (problem with structure id-ing)'
        return
    defect_index = -1
    for i in range(len(matches[2])):
        if matches[2][i] == None:  #don't use 'not matches[2][i]' because 0 might be in list
            if defect_index != -1:
                print 'More than one defect identified!(problem with structure id-ing)' \
                      ' Following is association matrix:', matches[2]
                return
            defect_index = i
            print 'Found defect index to be :',defect_index,'. This is site: ',struct.sites[defect_index-1]
    if defect_index == -1:
        print 'Could not find defect index! (problem with structure id-ing)'
        return

    grid_sites = {}  # dictionary with indices keys in order of structure list

    if coords_are_cartesian:
        def_fcoord = struct.lattice.get_fractional_coords(r_def)
        def_ccoord = r_def
    else:
        def_fcoord = r_def
        def_ccoord = struct.lattice.get_cartesian_coords(r_def)

    for i, site in enumerate(struct.sites): #go over defect sites?
        if i==defect_index:
            continue
        elif i<defect_index:
            id=i #id is defect index (zero defined)
        else:
            id=i-1 #this assumes vacancy for now...
        radius=site.specie.atomic_radius
        cart_coord = site.coords
        frac_coord = site.frac_coords
        dcart_coord = defstruct.sites[id].coords #0 defined
        dfrac_coord = defstruct.sites[id].frac_coords
        #dist, img = struct.lattice.get_distance_and_image(def_fcoord,
        #                                                  frac_coord)
        #I think we should be doing this for defect cell...
        dist, img = struct.lattice.get_distance_and_image(def_fcoord,
                                                          dfrac_coord)
        # cart_reldef = np.dot((frac_coord + img), struct.lattice._matrix) \
        #               - def_ccoord
        cart_reldef = np.dot((dfrac_coord + img), defstruct.lattice._matrix) \
                  - def_ccoord
        if abs(norm(cart_reldef) - dist) > 0.001:
            print('image locater issue encountered for site=', i, 
                  ' distance should be ', dist, ' but calculated to be ', 
                  norm(cart_reldef))
            return
        grid_sites[i] = {
                'dist': dist, 'cart': cart_coord, 'cart_reldef': cart_reldef,
                # cart_reldef=cartesian coord with defect as origin
                # should try and make this a range of grid indices in the future?
                'bulkgrid': getgridind1(locpot_blk,  frac_coord, gridavg=float(radius)), #note this is a list of indices now
                'defgrid' : getgridind1(locpot_def,  dfrac_coord, gridavg=float(radius))}

    #do we need this?
    # Danny: Yes this info might matter when we consider defects that relax and move away from where it was placed in the bulk,
    #       but I will leave it out for now. We will put it in the "grid_sites" dictionary if we use it.
    #grid_defect = getgridind1(def_ccoord)

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


# #reference sxdefectalign call for this example:
# '~/sxdefectalign --vasp -a1 --relative --pos 0.0,0.0,0.0 --charge 3 --ecut 38.2192757 --eps 18.099 --vref LOCPOT_vref --vdef LOCPOT_vdef'

class ChargeCorrection(object):
    def __init__(self, axis, dielectric_tensor, pure_locpot_path, 
            defect_locpot_path, q, pos,coords_are_cartesian=False, energy_cutoff=520,
            madetol=0.0001, silence=False, q_model=None, optgamma=None):
        """
        Args:
            axis: axis to do Freysoldt averaging over. Has no effect on 
                 Kumagai correction, so better move to Freysoldt potalign
            dielectric_tensor: Macroscopic dielectric tensor 
                 Include ionic also if defect is relaxed, othewise ion clamped.
                 Can be a matrix, array or scalar.
            pure_locpot_path: Bulk Locpot file path
            defect_locpot_path: Defect Locpot file path
            q: Charge associated with the defect. Typically integer
            pos: Position of the defect in the cell
            coords_are_cartesian: False if pos is in relative co-ordinates
            energy_cutoff: Energy for plane wave cutoff (in eV).
                 If not given, Materials Project default 520 eV is used.
            madetol: Tolerance (double or float)
            silence: Flag for disabling/enabling  messages (Bool)
            q_model (QModel object): User defined charge for correction.
                 If not given, highly localized charge is assumed.
            optgamma: (For anisotropic code) If you have previously optimized gamma,
                put gamma value here for speed up calculation slightly.
        """
        self._axis = axis  #needs to be zero defined (0,1,2); says which axis to do planar averaging on...
        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self._dielectricconst = float(dielectric_tensor)
            self._dieltens = np.diag(
                np.ones(3) * dielectric_tensor)  #this would make kumagai correction pointless
        else:
            self._dieltens = np.array(dielectric_tensor)  #full dielectric tensor
            self._dielectricconst = np.mean(np.diag(self._dieltens))  #take dielconstant to be average of trace
        self._purelocpot = pure_locpot_path  #location of purelocpot, could change so that this is locpot object?
        self._deflocpot = defect_locpot_path  #location of defectlocpot
        self._madetol = madetol #tolerance for convergence of energy terms in eV
        self._q = q  #charge of defect (not of the homogen. background)
        self._encut = energy_cutoff  #encut (eV) for calculation
        self._pos = pos
        self._coords_are_cartesian = coords_are_cartesian
        self._silence = silence  #for silencing printflags
        if not q_model:
            self._q_model = QModel()
        self._optgamma=optgamma
        self._g_sum_done=False


    def freysoldt_correction(self, vb=None, vd=None, title=None):
        #runs correction. input v1 and v2 locpot objects for speeding things up
        if not self._silence:
            print 'This is Freysoldt Correction.'
        if not self._q:
            return 0.0
        if not vb:
            if not self._silence:
                print 'Load bulk locpot'
            vb = Locpot.from_file(self._purelocpot)
        if not vd:
            if not self._silence:
                print 'Load defect locpot'
            vd = Locpot.from_file(self._deflocpot)
        if not self._silence:
            print '\nRun PC energy'
        energy_pc = self.freysoldt_pc(vb.structure)
        if not self._silence:
            print '\nPC calc done, correction =', round(energy_pc, 4)
            print 'Now run potenttial alignment script'
        potalign = self.freysoldt_potalign(v1=vb, v2=vd, title=title)
        if not self._silence:
            print '\n\nAlright so the corrections are:'
            print 'PCenergy = ', round(energy_pc, 5), '  potential alignment = ', round(potalign, 5)
            print 'TOTAL Freysoldt correction = ', round(energy_pc - potalign, 5)

        return round(energy_pc - potalign, 5)

    def freysoldt_pc(self, s1=None):
        #note this ony needs structural info
        # so s1=structure object speeds this calculation up alot
        if not s1:
            print 'load Pure locpot'
            s1 = Locpot.from_file(self._purelocpot).structure

        ap = s1.lattice.get_cartesian_coords(1)
        if self._silence == False:
            print 'run Freysoldt 2011 PC calculation (should be equivalent to sxdefectalign)'
            print 'defect lattice constants are (in angstroms)' + str(cleanlat(ap))
        [a1, a2, a3] = ang_to_bohr * ap
        if self._silence == False:
            print 'In atomic units, lat consts are (in bohr):' \
                  + str(cleanlat([a1, a2, a3]))
        vol = np.dot(a1, np.cross(a2, a3))  #vol in bohr^3

        #compute isolated energy
        step = 0.0001
        encut1 = 20.  #converge to some smaller encut [eV]
        flag = 0
        converge = []
        while (flag != 1):
            eiso = 1.
            gcut = eV_to_k(encut1)  #gcut is in units of 1/A
            g = step  #initalize
            while g < (gcut + step):
                #simpson integration
                eiso += 4*(self._q_model.rho_rec(g*g) ** 2)
                eiso += 2*(self._q_model.rho_rec((g+step) ** 2) ** 2)
                g += 2. * step
            eiso -= self._q_model.rho_rec(gcut ** 2) ** 2
            eiso *= (self._q ** 2) * step / (3 * round(np.pi, 6))
            converge.append(eiso)
            if len(converge) > 2:
                if abs(converge[-1] - converge[-2]) < self._madetol:
                    flag = 1
                elif encut1 > self._encut:
                    print 'Error encountered: Eiso did not converge before ' + str(self._encut) + ' eV'
                    sys.exit()
            encut1 += 20
        eiso = converge[-1]
        if self._silence == False:
            print 'Eisolated : ' + str(round(eiso, 5)) + ' converged at encut=' + str(encut1 - 20)

        #compute periodic energy;
        encut1 = 20.  #converge to some smaller encut
        flag = 0
        converge = []
        while flag != 1:
            eper = 0.0
            recip1 = generate_reciprocal_vectors_squared(a1, a2, a3, encut1)
            for g2 in recip1:
                eper += (self._q_model.rho_rec(g2) ** 2) / g2
            eper *= (self._q ** 2) * 2 * round(np.pi, 6) / vol
            eper += (self._q ** 2) * 4 * round(np.pi, 6) * self._q_model.rho_rec_limit0() / vol
            converge.append(eper)
            if len(converge) > 2:
                if abs(converge[-1] - converge[-2]) < self._madetol:
                    flag = 1
                elif encut1 > self._encut:
                    print 'Error encountered: Eper did not converge before ' + str(self._encut) + ' eV'
                    return
            encut1 += 20
        eper = converge[-1]

        if self._silence == False:
            print 'Eperiodic : ' + str(round(eper, 5)) + ' converged at encut=' + str(encut1 - 20)
            print 'difference (periodic-iso) is ' + str(round(eper - eiso, 6)) + ' hartree'
            print 'difference in (eV) is ' + str(round((eper - eiso) * hart_to_ev, 4))  #27.2114 eV/1 hartree
        PCfreycorr = round(((eiso - eper) / self._dielectricconst) * hart_to_ev, 6)  #converted to eV
        if self._silence == False:
            print 'Defect Correction without alignment (eV): ', PCfreycorr
        return PCfreycorr

    def freysoldt_potalign(self, v1=None, v2=None, title=None,  widthsample=1.):
        #NOTE this hasnt been coded for arbitrary defect position yet...
        #title is for name of plot, if you dont care about plot then leave it as None
        #widthsample is the width of the region in between defects where the potential alignment correction is averaged
        if not v1:
            print 'load pure locpot object'
            v1 = Locpot.from_file(self._purelocpot)
        if not v2:
            print 'load defect locpot object'
            v2 = Locpot.from_file(self._deflocpot)

        ind = []  #stores axes besides self._axis
        for i in range(3):
            if self._axis == i:
                continue
            else:
                ind.append(i)

        x = np.array(v1.get_axis_grid(self._axis))  #angstrom
        nx = len(x)
        print 'run Freysoldt potential alignment method'
        #perform potential alignment part
        pureavg = v1.get_average_along_axis(self._axis)  #eV
        defavg = v2.get_average_along_axis(self._axis)  #eV

        if not self._silence:
            print 'calculating lr part along planar avg axis'
        latt = v1.structure.lattice
        reci_latt = latt.reciprocal_lattice
        dg = reci_latt.abc[self._axis]
        v_G = np.empty(len(x), np.dtype('c16'))
        epsilon = self._dielectricconst
        # q needs to be that of the back ground? so added a (-) sign...
        v_G[0] = 4*np.pi * -self._q /epsilon * self._q_model.rho_rec_limit0()
        for i in range(1,nx):
            if (2*i < nx):
                g = i * dg
            else:
                g = (i-nx) * dg
            g2 = g*g
            v_G[i] = 4*np.pi / (epsilon * g2) * -self._q * self._q_model.rho_rec(g2)
        if not (nx % 2):
            v_G[nx/2] = 0
        v_R = np.fft.fft(v_G)
        v_R_imag = np.imag(v_R)
        v_R = np.real(v_R)* hart_to_ev
        v_R /= latt.volume

        max_imag_vr = v_R_imag.max()
        if abs(max_imag_vr) > self._madetol:
            print 'imaginary part found to be ', max_imag_vr, ' this is an issue'
            sys.exit()

        #now get correction and do plots
        short = (defavg - pureavg - v_R)
        checkdis = int((widthsample / 2) / x[1])  #index window for getting potential alignment correction
        mid = len(short) / 2
        tmppot = [short[i] for i in range(mid - checkdis, mid + checkdis)]

        C = np.mean(tmppot)
        Cquik = short[mid]  #this just uses mid point rather than average
        finalshift = [short[j] - C for j in range(len(v_R))]

        print 'short=',np.array(short)
        #get integrated final shift value
        vol=v1.structure.volume
        #alignmentterm_integrated=sum(finalshift) * x[1]/vol #not sure if should divide by vol
        alignmentterm_integrated=sum(finalshift)*x[1] #if locpot data is potential/unit vol this is right
        if not self._silence:
           print 'C value is averaged to be ' + str(C) + ' eV, '
           print 'compare with quicker C value (value halfway between defects) =',Cquik,' eV'
           print '\nCorresponds to an integrated alignment like term of ' + str(alignmentterm_integrated)+\
                    '\nIf this is not close to -C value then that should be sign of bad delocalization?'
           print 'will use eAlign term of ', -C, ' instead'
           print 'Pot. align correction is then (eV) : ' + str(-float(self._q) * float(C)) #= q*/Delta
        if title:
            plt.figure(1)
            plt.clf()
            plt.plot(x, v_R, c="green", zorder=1, label="long range from model")
            plt.plot(x, defavg - pureavg, c="red", label="DFT locpot diff")
            #plt.plot(x, short, c="purple", label="short range not shifted by C")
            plt.plot(x, finalshift, c="blue", label="short range (aligned)")
            plt.xlabel('planar average along axis ' + str(self._axis))
            plt.ylabel('Potential')
            plt.legend(loc=9)
            plt.axhline(y=0, linewidth=0.2, color='black')
            plt.axhline(y=C, linewidth=0.2, color='black')
            plt.title(str(title) + ' defect potential')
            plt.xlim(0,max(x))
            plt.savefig(str(title)+'FreyplnravgPlot.png')

        return -float(self._q)*C  #pot align energy correction (eV), add to the energy output of PCfrey

    def kumagai_correction1(self, locpot_blk=None, locpot_def=None, title=None):
        """
        Computes the extended Freysoldt correction for anistropic systems 
        developed by Y. Kumagai and F. Oba (Ref: PRB 89, 195205 (2014)

        If you want a plot of potential averaging process set title to name of defect
        vb and vd are preloaded locpot objects for speeding this up.
        """
        if not self._silence:
            print 'This is Kumagai Correction.'
        if not self._q:
            return 0.0
        if not locpot_blk:
            if not self._silence:
                print 'Load bulk locpot'
            locpot_blk = Locpot.from_file(self._purelocpot)

        structure = locpot_blk.structure
        if not self._optgamma:
            self._optgamma = find_optimal_gamma(structure, self._dieltens, self._q,
                                   self._madetol, silence=True)

        energy_pc = self.kumagai_pc1(locpot_blk)

        if not self._silence:
            #print '\nPC calc done, correction =', round(energy_pc, 4),
            print 'Now run potenttial alignment script'
        if not locpot_def:
            if not self._silence:
                print 'Load defect locpot'
            locpot_def = Locpot.from_file(self._deflocpot)

        potalign = self.kumagai_potalign1(locpot_blk=locpot_blk,
                                          locpot_def=locpot_def, title=title)
        if not self._silence:
            print '\n\nSummary of Kumagai corrections:'
            print 'PCenergy = ', round(energy_pc, 5), \
               '  potential alignment = ', round(potalign, 5)
            print 'TOTAL Kumagai correction = ', round(energy_pc - potalign, 5)
            return round(energy_pc - potalign, 5)
        return

    def kumagai_pc1(self, locpot_blk=None):
        #note that this ony needs structure info, not locpot info;
        # so s1=structure object speeds this calculation up alot
        if not locpot_blk:
            if not self._silence:
                print 'load structure from Pure locpot'
            locpot_blk = Locpot.from_file(self._purelocpot)
        if not self._silence:
            print '\nrun Kumagai PC calculation'

        if not self._optgamma:
            if not self._silence:
                print 'Optimize gamma first'
            self._optgamma = find_optimal_gamma(locpot_blk.structure,
                                self._dieltens, self._q, self._madetol, silence=True)

        if not self._g_sum_done:
            if not self._silence:
                print 'Calculating gsum'
            self._g_sum =  reciprocal_sum(locpot_blk, self._dieltens, self._q,
                                self._optgamma, silence=self._silence)
            self._g_sum_done=True

        energy_pc = anisotropic_pc_energy(
                locpot_blk.structure, self._g_sum, self._dieltens,  self._q,
                self._optgamma, self._madetol, silence=self._silence)  #returns PCenergy in eV

        if not self._silence:
            print 'PC energy determined to be ', energy_pc, ' eV (', \
                    energy_pc/hart_to_ev, ' Hartree)'  #27.2114 eV/1 hartree

        return energy_pc  #PC energy in eV

    def kumagai_potalign1(self, locpot_blk=None, locpot_def=None, title=None):
        """
        Potential alignment for Kumagai method
        Args:
            gamma: Convergence parameter
            locpot_blk: Bulk locpot object
            locpot_def: Defect locpot object
            title: Title for the plot. None will not generate the plot
        """
        if not self._silence:
            print ('\nrun Kumagai potential calculation (atomic site averaging)')
        if not locpot_blk:
            if not self._silence:
                print 'load pure locpot object'
            locpot_blk = Locpot.from_file(self._purelocpot)
        if not locpot_def:
            if not self._silence:
                print 'load defect locpot object'
            locpot_def = Locpot.from_file(self._deflocpot)

        if not self._optgamma:
            if not self._silence:
                print 'Optimizing gamma first'
            self._optgamma = find_optimal_gamma(locpot_blk.structure, self._dieltens, self._q,
                                   self._madetol, silence=True)

        if not self._g_sum_done:
            if not self._silence:
                print 'Calculating gsum'
            self._g_sum =  reciprocal_sum(locpot_blk, self._dieltens, self._q,
                                self._optgamma, silence=self._silence)
            self._g_sum_done=True

        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                locpot_blk.structure, self._dieltens, sil=self._silence)

        #this is to calculate distance matrix for plotting
        potinddict = disttrans1(locpot_blk, locpot_def, self._pos,
                coords_are_cartesian=self._coords_are_cartesian)

        wsrad = wigner_seitz_radius(locpot_blk.structure)
        if not self._silence:
            print ('wsrad', wsrad)
        
        for i in potinddict.keys():
            if potinddict[i]['dist'] > wsrad:
                potinddict[i]['OutsideWS'] = True
            else:
                potinddict[i]['OutsideWS'] = False

        #calculate potential value at each atomic site
        puredat = locpot_blk.data["total"]
        defdat = locpot_def.data["total"]
        jup = 0
        for i in potinddict.keys():
            jup += 1
            if (not title and not potinddict[i]['OutsideWS']):
                #dont need to calculate inside WS if not printing plot
                continue
            if not self._silence:
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
            print 'defdat val = ',np.mean(defvals)
            print 'puredat val = ',np.mean(bulkvals)
            v_qb = np.mean(defvals)-np.mean(bulkvals)
            cart_reldef = potinddict[i]['cart_reldef'] #should change this to averaging pure
            # and def within a range then subtract

            v_pc = anisotropic_madelung_potential(locpot_blk, self._g_sum,
                    cart_reldef, self._dieltens, self._q,self._optgamma,
                    self._madetol, silence=True)
            potinddict[i]['Vpc'] = v_pc
            potinddict[i]['Vqb'] = v_qb
            if not self._silence:
                print 'Has anisotropic madelung potential =', v_pc
                print 'DFT bulk/defect difference = ', v_qb
                print 'atoms left to calculate = ' + str(len(potinddict.keys()) - jup) #this is broken if you arent plotting
        if not self._silence:
            print '--------------------------------------'

        #now parse and plot if neccessary
        if title:  #to make shading region prettier
            fullspecset = locpot_blk.structure.species
            specset = list(set(fullspecset))
            shade, forplot = {}, {}
            for i in specset:
                shade[i.symbol] = {'r': [], 'Vpc': [], 'Vqb': []}
                forplot[i.symbol] = {'r': [], 'Vpc': [], 'Vqb': []}

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

        potalign = np.mean(forcorrection)

        if title:
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
                for k in range(len(forplot[inkey]['Vpc'])):
                    full.append([forplot[i]['r'][k], forplot[i]['Vqb'][k] - forplot[i]['Vpc'][k]])
            realfull = sorted(full, key=lambda x: x[0])
            r, y = [], []
            for i in realfull:
                r.append(i[0])
                y.append(i[1])
            plt.plot(r, y, color=collis[-1], marker='x', linestyle='None', label='V_{q/b} - Vpc')
            plt.xlabel('Distance from defect (A)')
            plt.ylabel('Potential (V)')
            x = np.arange(wsrad, max(locpot_blk.structure.lattice.abc), 0.01)
            plt.fill_between(x, min(ylis) - 1, max(ylis) + 1, facecolor='red', alpha=0.15, label='sampling region')
            plt.axhline(y=potalign, linewidth=0.5, color='red', label='pot. align.')
            plt.legend(loc=8)
            plt.axhline(y=0, linewidth=0.2, color='black')
            plt.ylim([min(ylis) - .5, max(ylis) + .5])
            plt.xlim([0, max(rlis) + 3])

            plt.title(str(title) + ' atomic site potential plot')
            #plt.show()
            plt.savefig(str(title) + 'kumagaisiteavgPlot.png')

        if self._silence == False:
            print 'Atomic site method potential alignment term is ' + str(np.mean(forcorrection))
            print 'this yields total (q*align) Kumagai potential correction energy of ' \
                  + str(self._q * np.mean(forcorrection)) + ' (eV) '

        return self._q * np.mean(forcorrection)


if __name__ == '__main__':
    s = ChargeCorrection(
            0, 18.099, 
            '../../../../Gavacm3testMachgcorr/LOCPOT_vref',
            '../../../../Gavacm3testMachgcorr/LOCPOT_vdef',
            -3, [0,0,0], silence=False,optgamma=2)

    #print 'load locpots'
    #vb=Locpot.from_file('../../Gavacm3testMachgcorr/LOCPOT_vref')
    #vb=Locpot.from_file('trans/LOCPOT_vref')
    #vd=Locpot.from_file('../../Gavacm3testMachgcorr/LOCPOT_vdef')


    #s.freysoldt_correction(title='Test')
    #s.freysoldt_pc()
    #s.freysoldt_pc(vb.structure)
    #s.freysoldt_potalign(v1=vb,v2=vd,title='Test')

    s.kumagai_correction1(title='Testing')
    #s.KumagaiPC(vb.structure)
    #s.Kumagaipotalign(vb,vd,optgam=1.848916,printflag='Gavac+3test')

