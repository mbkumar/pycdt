"""
This module computes finite size supercell charge corrections for
defects in anistropic systems using extended Freysoldt (or Kumagai) method 
developed by Kumagai and Oba.
Kumagai method includes
   a) anisotropic PC energy
   b) potential alignment by atomic site averaging at Wigner Seitz cell
      edge
If you use the corrections implemented in this module, cite
 a) Kumagai and Oba, Phys. Rev. B. 89, 195205 (2014) and
 b) Freysoldt, Neugebauer, and Van de Walle,
    Phys. Status Solidi B. 248, 1067-1076 (2011)  and
in addition to the pycdt paper
"""
__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import math
import logging

import numpy as np

from pymatgen.io.vasp.outputs import Locpot, Outcar
from pymatgen.core.lattice import Lattice

norm = np.linalg.norm

# Define conversion_constants
hart_to_ev = 27.2114
ang_to_bohr = 1.8897

# Define the logging stuff here
logging.basicConfig(filename='kumagai_debug.log',level=logging.DEBUG)


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
    """
    Args:
        a1, a2, a3: lattice vectors in bohr
        encut: energy cut off in eV
    Returns:
        reciprocal lattice vectors with energy less than encut
    """
    vol = np.dot(a1, np.cross(a2, a3))  # 1/bohr^3
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  # units 1/bohr
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)
    recip = []
    flag = 0

    # create list of recip space vectors that satisfy |i*b1+j*b2+k*b3|<=encut
    tol = 0
    while flag != 1:
        if 3.80986 * ((tol * (1 / ang_to_bohr) * min(norm(b1), norm(b2), norm(b3))) ** 2) < encut:
            #added the 1.8897 factor because the energy given converts 1/A to eV but b's in 1/bohr
            tol = tol + 1
        else:
            flag = 1

    for i in range(-tol, tol + 1):
        for j in range(-tol, tol + 1):
            for k in range(-tol, tol + 1):
                vec = (i * b1 + j * b2 + k * b3)
                en = 3.80986 * (((1 / ang_to_bohr) * norm(vec))** 2)
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


def kumagai_init(structure, dieltens, sil=True):
    angset = structure.lattice.get_cartesian_coords(1)

    dieltens = np.array(dieltens)
    if not len(dieltens.shape):
        dieltens = dieltens*np.identity(3)
    elif len(dieltens.shape) == 1:
        dieltens = np.diagflat(dieltens)

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


def get_g_sum_at_r(g_sum, structure, dim, r):
    """
    Args:
        g_sum: Reciprocal summation calculated from reciprocal_sum method
        structure: Bulk structure pymatgen object
        dim : ngxf dimension
        r: Position relative to defect (in cartesian coords)
    Returns:
        reciprocal summ value at g_sum[i_rx,j_ry,k_rz]
    """

    fraccoord=structure.lattice.get_fractional_coords(r)
    i, j, k = getgridind(structure, dim, fraccoord)

    return g_sum[i, j, k]


def anisotropic_madelung_potential(structure, dim, g_sum, r, dieltens, q,  gamma, tolerance,
                                  silence=True):
    """
    Compute the anisotropic Madelung potential at r not equal to 0.
    For r=(0,0,0) use anisotropic_pc_energy function
    Args:
        structure: Bulk pymatgen structure type
        dim : ngxf dimension
        g_sum: Precomputed reciprocal sum for all r_vectors
        r: r vector (in cartesian coordinates) relative to defect position. 
           Non zero r is expected
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        tolerance: Tolerance parameter for numerical convergence
        gamma (float): Convergence parameter 
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens, sil=silence)

    recippartreal = q * get_g_sum_at_r(g_sum, structure, dim, r)
    directpart = real_sum(a1, a2, a3, r, q, dieltens, gamma, tolerance,
                          silence)

    #now add up total madelung potential part with two extra parts:
    #self interaction term
    selfint = q * np.pi / (vol * (gamma ** 2))
    if not silence:
        print ('self interaction piece is {}'.format(selfint * hart_to_ev))

    #pot = (hart_to_ev/-q) * (directpart + recippartreal - selfint)
    pot = hart_to_ev * (directpart + recippartreal - selfint)  #reverted dividing by q to match kumagai data...

    return pot


def anisotropic_pc_energy(structure, g_sum, dieltens, q, gamma, tolerance,
                          silence=True):
    """
    Compute the anistropic periodic point charge interaction energy.
    Args:
        structure: Bulk pymatgen structure type
        g_sum : comes from KumagaiBulkInit class
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        gamma : convergence parameter optimized in KumagaiBulkInit class
        silence (bool): Verbosity flag. If False, messages are printed.
    """
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens, sil=silence)

    g_part = q*g_sum[0,0,0]
    r_part = real_sum(a1, a2, a3, [0,0,0], q, dieltens, gamma, tolerance, True)

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


def getgridind(structure, dim, r, gridavg=0.0):
    """
    Computes the index of a point, r, in the locpot grid
    Args:
        structure: Pymatgen structure object
        dim: dimension of FFT grid (NGXF dimension list in VASP)
        r: Relative co-ordinates with respect to abc lattice vectors
        gridavg: If you want to do atomic site averaging, set gridavg to the radius of the atom at r
    Returns:
        [i,j,k]: Indices as list
    TODO: Once final, remove the getgridind inside disttrans function
    """
    abc=structure.lattice.abc
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
        num_pts = dim[i]
        x = [now_num / float(num_pts) * abc[i] for now_num in range(num_pts)]
        dx = x[1] - x[0]
        x_rprojection_delta_abs = np.absolute(x - r[i])
        ind = np.argmin(x_rprojection_delta_abs)
        if x_rprojection_delta_abs[ind] > dx*1.1: #to avoid numerical errors
            print i,ind,r
            print x_rprojection_delta_abs
            raise ValueError("Input position is not within the locpot grid")
        grdind.append(ind)
        if gridavg:
            radvals.append(int(np.ceil(gridavg/dx)))
            dxvals.append(dx)

    if gridavg:
        grdindfull=[]
        for i in range(-radvals[0], radvals[0]+1):
            for j in range(-radvals[1], radvals[1]+1):
                for k in range(-radvals[2], radvals[2]+1):
                    dtoc = [i*dxvals[0], j*dxvals[1], k*dxvals[2]]
                    if norm(dtoc) < gridavg:
                        ival = (i+grdind[0]) % dim[0]
                        jval = (j+grdind[1]) % dim[1]
                        kval = (k+grdind[2]) % dim[2]
                        grdindfull.append((ival,jval,kval))
        grdind=grdindfull

    return grdind


def disttrans(struct, defstruct, dim, silence=False):
    """
    for calculating distance from defect to each atom and finding NGX grid pts at each atom
    Args:
        struct: Bulk structure object
        defstruct: Defect structure object
    """

    #Find defect location in bulk and defect cells
    blksite,defsite = find_defect_pos(struct,defstruct)
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
    for i in sitelist:
        if np.array_equal(i.coords,def_ccoord):
            print 'This is defect! Skipping ',i
            continue

        radius=1
        blksite,defsite=closestsites(struct,defstruct,i.coords)

        blkindex=struct.index(blksite[0])
        defindex=defstruct.index(defsite[0])

        cart_coord = blksite[0].coords
        frac_coord = blksite[0].frac_coords
        dcart_coord = defsite[0].coords
        dfrac_coord = defsite[0].frac_coords

        closeimage=returnclosestr(dcart_coord)
        cart_reldef=closeimage[1]
        defdist=closeimage[0]

        if abs(norm(cart_reldef) - defdist) > 0.1:
            print('image locater issue encountered for site=', blkindex,
                  '(in def cell) distance should be ', defdist, ' but calculated to be ',
                  norm(cart_reldef))
            #return #dont want to break the code here, but want flag to exist...what to do?

        rad = 1.0
        if blkindex in grid_sites:
            print '(WARNING) index ',blkindex,' already exists in potinddict! overwriting information. '

        grid_sites[blkindex] = {'dist': defdist,'cart': dcart_coord,
                'cart_reldef': cart_reldef,
                'defgrid':getgridind(struct, dim,  dfrac_coord, gridavg=rad),
                'bulkgrid': getgridind(defstruct, dim,  frac_coord, gridavg=rad),
                'siteobj':[i.coords,i.frac_coords,i.species_string],
                'bulk_site_index':blkindex, 'def_site_index':defindex}

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

    dist = []
    for facet in wz:
        midpt = np.mean(np.array(facet), axis=0)
        dist.append(norm(midpt))
    wsrad = min(dist)

    return wsrad


def read_ES_avg(location_outcar):
    #reads NGXF information and Electrostatic potential at each atomic site from VASP OUTCAR file
    with open(location_outcar,'r') as file:
        tmp_out_dat = file.read()

        out_dat = tmp_out_dat.split('\n')
        start_line = 0
        end_line = 0
        ngxf_line = 0
        for line_num, line in enumerate(out_dat):
            if "dimension x,y,z NGXF" in line:
                ngxf_line = line_num
            elif "average (electrostatic) potential at core" in line:
                start_line = line_num
                end_line = 0 #make sure to zero end_line if there are multiple electrostatic read outs
            elif start_line and not end_line and not len(line.split()):
                end_line = line_num

        ngxlineout = out_dat[ngxf_line].split()
        ngxf_dims = map(int, ngxlineout[3:8:2])

        rad_line = out_dat[start_line+1].split()
        radii = [float(rad) for rad in rad_line[5:]]

        ES_data = {'sampling_radii': radii, 'ngxf_dims': ngxf_dims}
        pot = []
        for line_num in range(start_line+3, end_line):
            line = out_dat[line_num].split()
            avg_es = map(float,line[1::2])
            pot += avg_es
        ES_data.update({'potential': pot})

        return ES_data

    return None


class KumagaiBulkInit(object):
    """
    Compute the anisotropic madelung potential array from the bulk 
    locpot. This helps in evaluating the bulk supercell related part 
    once to speed up the calculations.
    """
    def __init__(self, structure, dim, epsilon, encut=520, tolerance=0.0001,
                 silence=True, optgamma=False):
        """
        Args
            structure: 
                Pymatgen structure object of bulk cell
            dim: 
                Fine FFT grid dimensions as a list 
                For vasp this is NGXF grid dimensions
            epsilon: 
                Dielectric tensor
            encut (float): 
                Energy cutoff for optimal gamma
            tolerance (float): 
                Accuracy parameter
            silence (bool): 
                If True, intermediate steps are not printed
            optgamma: 
                if you know optimized gamma, give its value. 
                Otherwise it will be computed.
        """
        self.structure = structure
        self.dim = dim
        self.epsilon = epsilon
        self.encut = encut
        self.tolerance = tolerance
        self.silence = silence
        if not optgamma:
            self.gamma = self.find_optimal_gamma()
        else:
            self.gamma = optgamma
        self.g_sum = self.reciprocal_sum()
        logging.info('optimized gaama: %f', self.gamma)
        logging.info('g_sum: %f', self.g_sum)

    def find_optimal_gamma(self):
        """
        Find optimal gamma by evaluating the brute force reciprocal
        summation and seeing when the values are on the order of 1
        this calculation is the anisotropic Madelung potential at r = (0, 0, 0)
        Note this only requires the STRUCTURE not the LOCPOT object.
        """
        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                self.structure, self.epsilon, sil=self.silence)
        optgam = None

        #do brute force recip summation
        def get_recippart(encut, gamma):
            recip = genrecip(a1, a2, a3, encut)
            recippart = 0.0
            for rec in recip:
                Gdotdiel = np.dot(rec, np.dot(self.epsilon, rec))
                summand = math.exp(-Gdotdiel / (4 * (gamma ** 2))) / Gdotdiel
                recippart += summand
            recippart *= 4*np.pi/vol
            return recippart, 0.0, len(recip)

        def do_summation(gamma):
            #do recip sum until it is bigger than 1eV
            #First do Recip space sum convergence with respect to encut for this gamma
            encut = 20  #start with small encut for expediency
            recippartreal1, recippartimag1, len_recip = get_recippart(encut, gamma)
            encut += 10
            recippartreal, recippartimag, len_recip = get_recippart(encut, gamma)
            converge = [recippartreal1, recippartreal]
            while abs(abs(converge[0]) - abs(converge[1])) * hart_to_ev > self.tolerance:
                encut += 10
                recippartreal, recippartimag, len_recip = get_recippart(encut, gamma)
                converge.reverse()
                converge[1] = recippartreal
                if encut > self.encut:
                    raise ValueError(
                            'Optimal gamma not found at {} eV cutoff'.format(
                                self.encut))

            if abs(recippartimag) * hart_to_ev > self.tolerance:
                #if not self.silence:
                logging.error("Problem with convergence of imaginary part of recip sum."),
                logging.error("imag sum value is {} (eV)".format( recippartimag * hart_to_ev))
                return None, None
            #if not self.silence:
            logging.info('recip sum converged to {} (eV) at encut= {}'.format(
                        recippartreal * hart_to_ev, encut))
            logging.info('Number of reciprocal vectors is {}'.format(len_recip))
            if (abs(converge[1]) * hart_to_ev < 1 and not optgam):
                #if not self.silence:
                logging.info('Reciprocal summation value is less than 1 eV.')
                logging.info('This might lead to errors in the reciprocal summation.')
                logging.info('Changing gamma now.')
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
                logging.error('Had problem in gamma optimization process.')
                return None

            if gamma > 50:
                logging.error('Could not optimize gamma before gamma = %d', 50)
                return None

        return optgamma 

    def reciprocal_sum(self):
        """
        Compute the reciprocal summation in the anisotropic Madelung 
        potential.

        TODO: Get the input to fft cut by half by using rfft instead of fft
        """

        if not self.silence:
            print 'calculating the reciprocal summation in Madeling potential'
        over_atob = 1.0/ang_to_bohr
        atob3=ang_to_bohr**3

        latt = self.structure.lattice
        vol = latt.volume*atob3 # in Bohr^3

        reci_latt = latt.reciprocal_lattice
        [b1, b2, b3] = reci_latt.get_cartesian_coords(1)
        b1 = np.array(b1)*over_atob # In 1/Bohr
        b2 = np.array(b2)*over_atob
        b3 = np.array(b3)*over_atob

        nx, ny, nz = self.dim
        logging.debug('nx: %d, ny: %d, nz: %d', nx, ny, nz)
        ind1 = np.arange(nx)
        for i in range(nx/2, nx):
            ind1[i] = i - nx
        ind2 = np.arange(ny)
        for i in range(ny/2, ny):
            ind2[i] = i - ny
        ind3 = np.arange(nz)
        for i in range(nz/2, nz):
            ind3[i] = i - nz

        g_array = np.zeros(self.dim, np.dtype('c16'))
        gamm2 = 4*(self.gamma**2)
        for i in ind1:
            for j in ind2:
                for k in ind3:
                    g = i*b1 + j*b2 + k*b3
                    g_eps_g = np.dot(g, np.dot(self.epsilon, g))
                    if i == j == k == 0:
                        continue
                    else:
                        g_array[i,j,k] = math.exp(-g_eps_g/gamm2)/g_eps_g

        r_array = np.fft.fftn(g_array)
        over_vol = 4*np.pi/vol # Multiply with q later
        r_array *= over_vol
        r_arr_real = np.real(r_array)
        r_arr_imag = np.imag(r_array)

        #if not self.silence:
        max_imag = r_arr_imag.max()
        logging.debug('Max imaginary part found to be %f', max_imag)

        return r_arr_real


class KumagaiCorrection(object):
    """
    Class that implements the extended freysoldt correction developed 
    by Kumagai.
    """
    def __init__(self, dielectric_tensor, q, gamma, g_sum, bulk_structure,
            energy_cutoff=520, madetol=0.0001, silence=False,
            lengths=None,  defstructure=None, **kw):
        """
        Args:
            dielectric_tensor: Macroscopic dielectric tensor
                 Include ionic also if defect is relaxed, othewise ion clamped.
                 Can be a matrix array or scalar.
            q: Charge associated with the defect (not of the homogen. background). Typically integer
            gamma:  Convergence parameter. Obtained from KumagaiBulkPart
            g_sum: value that is dependent on the Bulk only. Obtained from KumagaiBulkPart
            energy_cutoff: Energy for plane wave cutoff (in eV).
                 If not given, Materials Project default 520 eV is used.
            madetol: Tolerance for convergence of energy terms in eV (double or float)
            silence: Flag for disabling/enabling  messages (Bool)
            lengths: Lengths of axes, for speeding up plotting slightly
            structure: bulk Pymatgen structure object. Need to specify this if using Outcar method for atomic site avg.
                (If you specify outcar files for bulk_file_path but dont specify structure then code will break)
                (TO DO: resolve this dumb dependency by being smarter about where structure comes from?)
            defstructure: defect Pymatgen structure object. only needed if using Outcar method...
            keywords:
                1) bulk_locpot: Bulk Locpot file path OR Bulk Locpot Object
                   defect_locpot: Defect Locpot file path or defect Locpot Object
                2) (Or) bulk_outcar:   Bulk Outcar file path OR Bulk Outcar Object
                   defect_outcar: Defect outcar file path or defect outcar Object
        """
        if not silence:
            print '\nThis is Anisotropic Freysoldt (Kumagai) Correction'
        if isinstance(dielectric_tensor, int) or \
                isinstance(dielectric_tensor, float):
            self.dieltens = np.identity(3) * dielectric_tensor
        else:
            self.dieltens = np.array(dielectric_tensor)

        if 'bulk_locpot' in kw:
            if isinstance(kw['bulk_locpot'], Locpot):
                self.locpot_blk = bulk_file_path
            else:
                self.locpot_blk = Locpot.from_file(kw['bulk_locpot'])
            if isinstance(kw['defect_locpot'], Locpot):
                self.locpot_def = kw['defect_locpot']
            else:
                self.locpot_def = Locpot.from_file(kw['defect_locpot'])
            self.dim = self.locpot_blk.dim

            self.outcar_blk = None
            self.outcar_def = None
            self.do_outcar_method = False

        if 'bulk_outcar' in kw:
            self.outcar_blk = kw['bulk_outcar']
            self.outcar_def = kw['defect_outcar']
            self.do_outcar_method = True
            self.locpot_blk = None
            self.locpot_def = None
            #this would be part where I read dims in from Outcar pymatgen attribute
            #for now use hack function
            tmpdict = read_ES_avg(self.outcar_blk)
            self.dim = tmpdict['ngxf_dims']

        self.madetol = madetol
        self.q = q
        self.encut = energy_cutoff
        self.silence = silence
        self.structure = bulk_structure
        self.defstructure = defstructure
        self.gamma = gamma
        self.g_sum = g_sum

        self.lengths=lengths

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
                print 'Now run potential alignment script'
            potalign = self.potalign(title=title)

        if not self.silence:
            print '\n\nKumagai Correction details:'
            if partflag!='potalign':
                print 'PCenergy (E_lat) = ', round(energy_pc, 5)
            if partflag!='pc':
                print 'potential alignment (-q*delta V) = ', round(potalign, 5)
            if partflag in ['All','AllSplit']:
                print 'TOTAL Kumagai correction = ', round(energy_pc + potalign, 5)

        if partflag=='pc':
            return round(energy_pc,5)
        elif partflag=='potalign':
            return round(potalign,5)
        elif partflag=='All':
            return round(energy_pc+potalign,5)
        else:
            return [round(energy_pc,5),round(potalign,5),round(energy_pc+potalign,5)]

    def pc(self):
        if not self.silence:
            print '\nrun Kumagai PC calculation'

        energy_pc = anisotropic_pc_energy(
                self.structure, self.g_sum, self.dieltens, self.q,
                self.gamma, self.madetol, silence=self.silence)

        if not self.silence:
            print 'PC energy determined to be ', energy_pc, ' eV (', \
                    energy_pc/hart_to_ev, ' Hartree)'

        return energy_pc

    def potalign(self, title=None):
        """
        Potential alignment for Kumagai method
        Args:
            title: Title for the plot. None will not generate the plot
        """
        if not self.silence:
            print ('\nrun Kumagai potential calculation (atomic site averaging)')

        if (not type(self.locpot_blk) is Locpot) and not self.structure: #if structure not specified the locpot_path is loaded
            if not self.silence:
                print 'Load bulk Locpot'
            self.locpot_blk=Locpot.from_file(self.locpot_blk)
            self.structure = self.locpot_blk.structure
        elif not self.structure:
            self.structure = self.locpot_blk.structure

        if (not type(self.locpot_def) is Locpot) and not self.defstructure:
            if not self.silence:
                print 'Load defect Locpot'
            self.locpot_def=Locpot.from_file(self.locpot_def)
            self.defstructure = self.locpot_def.structure
        elif not self.defstructure:
            self.defstructure = self.locpot_def.structure

        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                self.structure, self.dieltens, sil=self.silence)

        potinddict = disttrans(self.structure, self.defstructure, self.dim, silence=self.silence)

        minlat=min(norm(a1),norm(a2),norm(a3))
        lat_perc_diffs=[100*abs(norm(a1)-norm(lat))/minlat for lat in [a2,a3]]
        lat_perc_diffs.append(100*abs(norm(a2)-norm(a3))/minlat)
        if not all(i < 45 for i in lat_perc_diffs):
            print 'NOTICE! detected that cell was not very cubic. ' \
                  'Might want to be smarter in way you sample atoms outside wigner-seitz cell with Kumagai scheme'
        wsrad = wigner_seitz_radius(self.structure)
        if not self.silence:
            print ('wsrad', wsrad)

        for i in potinddict.keys():
            if potinddict[i]['dist'] > wsrad:
                potinddict[i]['OutsideWS'] = True
            else:
                potinddict[i]['OutsideWS'] = False

        if not self.do_outcar_method:
            puredat = self.locpot_blk.data["total"]
            defdat = self.locpot_def.data["total"]
        else:
            #note this is hack until we get OUTCAR object attribute working for ES potential
            puredat = read_ES_avg(self.outcar_blk)
            defdat = read_ES_avg(self.outcar_def)

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

            if not self.do_outcar_method:
                #NOTE this method needs to be improved a bit by specifying radius type based on ENAUG or something?

                ##single point routine
                #dx, dy, dz = potinddict[i]['defgrid']
                #dx, dy, dz = potinddict[i]['bulkgrid']
                #bx, by, bz = potinddict[i]['bulkgrid']
                #v_qb = defdat[dx][dy][dz] - puredat[bx][by][bz]

                #averaging point routine
                bulkvals=[]
                defvals=[]
                for u,v,w in potinddict[i]['bulkgrid']:
                    bulkvals.append(puredat[u][v][w])
                for u,v,w in potinddict[i]['defgrid']:
                    defvals.append(defdat[u][v][w])
                print 'defdat val = ',np.mean(defvals)
                print 'puredat val = ',np.mean(bulkvals)
                v_qb = np.mean(defvals) - np.mean(bulkvals)
            else:
                defindex = potinddict[i]['def_site_index'] #assuming this is zero defined...
                bulkindex = potinddict[i]['bulk_site_index']
                v_qb = defdat['potential'][defindex] - puredat['potential'][bulkindex]


            cart_reldef = potinddict[i]['cart_reldef']
            v_pc = anisotropic_madelung_potential(self.structure, self.dim, self.g_sum,
                    cart_reldef, self.dieltens, self.q, self.gamma,
                    self.madetol, silence=True)
            v_qb*=-1 #change charge sign convention

            potinddict[i]['Vpc'] = v_pc
            potinddict[i]['Vqb'] = v_qb
            if not self.silence:
                print 'Has anisotropic madelung potential =', v_pc
                print 'DFT bulk/defect difference = ', v_qb
                print 'atoms left to calculate = ' + str(len(potinddict.keys()) - jup)
        if not self.silence:
            print '--------------------------------------'

        if title:
            fullspecset = self.structure.species
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

        potalign = np.mean(forcorrection)

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
                             label=str(inkey) + ': $V_{q/b}$')
                    plt.plot(forplot[inkey]['r'], forplot[inkey]['Vpc'], color=collis[i], marker='o', linestyle='None',
                             label=str(inkey) + ': $V_{pc}$')
                full = []
                for i in forplot.keys():
                    for k in range(len(forplot[i]['Vpc'])):
                        full.append([forplot[i]['r'][k], forplot[i]['Vqb'][k] - forplot[i]['Vpc'][k]])
                realfull = sorted(full, key=lambda x: x[0])
                r, y = [], []
                for i in realfull:
                    r.append(i[0])
                    y.append(i[1])
                plt.plot(r, y, color=collis[-1], marker='x', linestyle='None', label='$V_{q/b} - V_{pc}$')
                plt.xlabel('Distance from defect (A)')
                plt.ylabel('Potential (V)')
                x = np.arange(wsrad, max(self.structure.lattice.abc), 0.01)
                plt.fill_between(x, min(ylis) - 1, max(ylis) + 1, facecolor='red', alpha=0.15, label='sampling region')
                plt.axhline(y=potalign, linewidth=0.5, color='red', label='pot. alignment')
                plt.legend()
                plt.axhline(y=0, linewidth=0.2, color='black')
                plt.ylim([min(ylis) - .5, max(ylis) + .5])
                plt.xlim([0, max(rlis) + 3])

                plt.title(str(title) + ' atomic site averaging potential plot')
                plt.savefig(str(title) + 'kumagaisiteavgPlot.pdf')
            else:
                from monty.serialization import dumpfn
                from monty.json import MontyEncoder
                forplot['EXTRA']={'wsrad':wsrad,'potalign':potalign}
                fname='KumagaiData.json'
                dumpfn(forplot, fname, cls=MontyEncoder)

        if self.silence == False:
            print 'Atomic site method potential alignment term is ' + str(np.mean(forcorrection))
            print 'this yields total (-q*align) Kumagai potential correction energy of ' \
                  + str(- self.q * np.mean(forcorrection)) + ' (eV) '

        return -self.q * np.mean(forcorrection)

    def plot_from_datfile(self,name='KumagaiData.json',title='default'):
        """
        Takes data file called 'name' and does plotting.
        Good for later plotting of locpot data after running run_correction()

        """
        if type(self.locpot_blk) is not Locpot and not self.lengths:
            self.locpot_blk=Locpot.from_file(self.locpot_blk)

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
