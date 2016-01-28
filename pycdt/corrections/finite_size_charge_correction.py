"""
chg corrections:
1) Freysoldt correction 
	a) PC energy
	b) potential alignment [axis averaging]
2) Kumagai correction
	a) anisotropic PC energy
	b) potential alignment [atomic site averaging outside WS cell]

3) simple Madelung correction? 
I think this might be possible with some simple pymatgen commands?

"""

__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com'

from pymatgen.io.vaspio.vasp_output import Locpot
from pymatgen.core.lattice import Lattice
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import math

norm = np.linalg.norm  # define globally

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
        if 3.80986 * ((tol * (1 / 1.8897) * min(norm(b1), norm(b2), norm(b3))) ** 2) < encut:
            #added the 1.8897 factor because the energy given converts 1/A to eV but b's in 1/bohr
            tol = tol + 1
        else:
            flag = 1
    #now look though all options for recip vectors to see what vectors are less than energy val
    for i in range(-tol, tol + 1):
        for j in range(-tol, tol + 1):
            for k in range(-tol, tol + 1):
                vec = (i * b1 + j * b2 + k * b3)
                en = 3.80986 * (((1 / 1.8897) * norm(vec))** 2)
                if (en <= encut and en != 0):
                    recip.append([i * b1[m] + j * b2[m] + k * b3[m] for m in range(3)])
    return recip  #output is 1/bohr recip

def generate_3d_reciprocal_vectors(a1, a2, a3, encut, gcutflag=False):
    # latt vectors in 1/bohr, encut=eV
    # generate reciprocal lattice vectors with value less than encut
    # define recip vectors first, (units of 1/angstrom).
    # gcut flag =True just quits and gives you gcut rather than computing all reciprocal lattice vectors
    vol = np.dot(a1, np.cross(a2, a3))  # 1/bohr^3
    b1 = (2 * np.pi / vol) * np.cross(a2, a3)  # units 1/bohr
    b2 = (2 * np.pi / vol) * np.cross(a3, a1)
    b3 = (2 * np.pi / vol) * np.cross(a1, a2)
    recip = []
    flag = 0
    # create list of recip space vectors that satisfy |i*b1+j*b2+k*b3|<=encut
    #start by enumerating to find max i that doesn't upset the encut condition
    tol = 0
    gcut = 0.  #this is max number of repetitions for minimum recip vector to upset encut condition
    #print 'value of smallest b vec is '+str(min(mod(b1),mod(b2),mod(b3)))
    while flag != 1:
        if 3.80986 * ((tol * (1 / 1.8897) * min(norm(b1), norm(b2), norm(b3))) ** 2) < encut:
            #added the 1.8897 factor because the energy given converts 1/A to eV but b'2 in 1/bohr
            tol = tol + 1
        else:
            #print ('tolerance is',tol)
            gcut = tol * min(norm(b1), norm(b2), norm(b3))  #1/bohr
            flag = 1
    if gcutflag:
        print ('gcut', gcut)
        return gcut
    print ('tol', tol)
    #now look though all options for recip vectors to see what vectors are less than energy val
    for i in range(-tol, tol + 1):
        for j in range(-tol, tol + 1):
            for k in range(-tol, tol + 1):
                vec = (i * b1 + j * b2 + k * b3)
                en = 3.80986 * (((1 / 1.8897) * norm(vec))** 2)
                if (en <= encut and en != 0):
                    recip.append([i * b1[m] + j * b2[m] + k * b3[m] for m in range(3)])
    return recip, gcut  #output is 1/bohr recip and 1/bohr gcut

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
        print 'defect lattice constants are (in angstroms)' + str(angset)
    [a1, a2, a3] = 1.8897 * angset  # convert to bohr
    bohrset = [a1, a2, a3]
    vol = np.dot(a1, np.cross(a2, a3))
    if not sil:
        print 'converted to bohr for atomic units, lat consts are:' + str([a1, a2, a3])
    determ = np.linalg.det(dieltens)
    invdiel = np.linalg.inv(dieltens)
    if not sil:
        print 'inv dielectric tensor is ' + str(invdiel)
    return angset, bohrset, vol, determ, invdiel


def get_anistropic_pc_energy(structure, dieltens, q, madetol, r, silence, 
        gamma=None):
    """
    Compute the anistropic point charge energy correction.
    Args:
        structure: pymatgen Structure 
        dieltens: dielectric tensor
        q: Point charge (in units of e+)
        madetol: Accuracy parameter 
        r: Position with respect to defect
        silence: Verbosity flag. If off, messages are printed.
        gamma: Convergence parameter (optional)
    """
    # if r=[0,0,0] return PCenergy, otherwise return the potential energy part
    # if gamma has already been optimized, then set optgam to the optimized gamma
    angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
            structure, dieltens, sil=silence)
    
    if not gamma: # What is the optimal gamma and how to infer it?
        gamma = 5./(vol ** (1/3.))
    if not silence:
        print ('Gamma not yet optimized for Kumagai calc. Initial gamma is ', 
                gamma)

    # for recip summation part
    if norm(r):  #produces potential term
        # Comment by Bharat to Danny: 
        # The following can be realized with FFT
        def get_recippart(encut, gamma):
            recip = genrecip(a1, a2, a3, encut)
            recippartreal, recippartimag = 0.0, 0.0
            for rec in recip:
                Gdotdiel = np.dot(rec, np.dot(dieltens, rec))
                Gdotr = np.dot(rec, r)
                summand = math.exp(-Gdotdiel / (4 * (gamma ** 2))) / Gdotdiel
                recippartreal += summand * math.cos(Gdotr)
                recippartimag += summand * math.sin(Gdotr)
            recippartreal *= 4*np.pi*q/vol
            recippartimag *= 4*np.pi*q/vol
            return recippartreal, recippartimag, len(recip)
    else:  #produces PC energy term
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
        directlist = []
        Nmaxlength = 40  #tolerance for stopping real space sum convergence
        N = 2
        realpre = q / np.sqrt(determ)
        #Real space sum by converging with respect to real space vectors
        #create list of real space vectors that satisfy |i*a1+j*a2+k*a3|<=N
        while N < Nmaxlength:  
            realvecsum = []
            if not norm(r):
                for i in range(-N, N + 1):
                    for j in range(-N, N + 1):
                        for k in range(-N, N + 1):
                            vec = i * a1 + j * a2 + k * a3
                            if norm(vec):
                                realvecsum.append(vec)
            else:  # if r!=0 then include zero in summation
                for i in range(-N, N + 1):
                    for j in range(-N, N + 1):
                        for k in range(-N, N + 1):
                            realvecsum.append(i * a1 + j * a2 + k * a3 - r)

            #calculation real summation up to N
            directpart = 0.0
            for i in realvecsum:
                local_response = np.dot(i, np.dot(invdiel, i))
                nmr = math.erfc(gamma * np.sqrt(local_response))
                dmr = np.sqrt(determ * local_response)
                directpart += nmr / dmr  
            directlist.append([N, realpre * directpart])

            if N == Nmaxlength-1:
                print('Direct part could not converge with real space ' +
                       'translation tolerance of {} for gamma {}'.format(
                           Nmaxlength-1, gamma))
                return
            elif len(directlist) > 3:
                if abs(abs(directlist[-1][1]) - abs(directlist[-2][1])) * 27.2114 < madetol:
                    directpart = directlist[-1][1]
                    if not silence:
                        print("gamma is {}".format(gamma))
                        print("convergence for direct term occurs at step " + 
                               "{}  where direct sum is {}".format(
                                   N,  directpart * 27.2114))
                        print('There are {} real vectors'.format(len(realvecsum)))
                    break

            N += 1

        #Recip space sum convergence with respect to encut
        encut = 20  #start with small encut for expediency
        recippartreal1, recippartimag1, len_recip = get_recippart(encut, gamma)
        encut += 10
        recippartreal, recippartimag, len_recip = get_recippart(encut, gamma)
        converge = [recippartreal1, recippartreal]
        while abs(abs(converge[0]) - abs(converge[1])) * 27.2114 > madetol:
            encut += 10
            recippartreal, recippartimag, len_recip = get_recippart(encut, gamma)
            converge.reverse()
            converge[1] = recippartreal
            if encut > 700: # Bharat: Why 700 eV?
                # Danny: Because if it isnt converged by then somethign is wrong...could do somethign smaller like 300eV?
                print('Problem, imaginary part not converged at encut = 700eV')
                return
        if abs(recippartimag) * 27.2114 > madetol:
            print("Problem with convergence of imaginary part of recip sum."),
            print("imag sum value is {} (eV)".format( recippartimag * 27.2114))
            return None, None, None
        if not silence:
            print('recip sum converged to {} (eV) at encut= {}'.format(
                        recippartreal * 27.2114, encut))
            print('Number of reciprocal vectors is {}'.format(len_recip))
            if (abs(converge[1]) * 27.2114 < 1 and not optgam):  
                #only optimize is flag set for optimizing routine
                print('Warning: reciprocal summation value is less than 1 eV.')
                print('Last recip sum value = {}.'.format(converge[1])) 
                print('This might lead to errors in the reciprocal summation.')
                print('Changing gamma now.')
                return None, None, 'Try Again'

        return directpart, recippartreal, gamma

    #start with gamma s.t. gamma*L=5 (some paper said this is optimal)
    if not optgam:
        #optimizing gamma for the recipt sum to improve convergence of calculation
        flag = 0
        while not flag:
            directpart, recippartreal, optgamma = do_summation(gamma)
            if optgamma == gamma:
                print('optimized gamma found to be ', optgamma)
                flag += 1
            elif 'Try Again' in optgamma:
                gamma *= 1.5
            else:
                print('Had problem in gamma optimization process.')
                return

            if gamma > 50:  #kind of an arbitrary number for cut off...
                print('WARNING. could not optimize gamma before gamma =', 50)
                return
    else:
        directpart, recippartreal, optgamma = do_summation(optgam)

    #now add up total madelung potential part with two extra parts:
    #self interaction term
    selfint = q * np.pi / (vol * (optgamma ** 2))
    if not silence:
        print ('self interaction piece is {}'.format(selfint * 27.2114))

    #surface term (only for r not at origin)
    if not norm(r):
        surfterm = 2 * optgamma * q / np.sqrt(np.pi * determ)
        if not silence:
            print ('surface term is {}'.format(surfterm * 27.2114))
        totalPC = -q * 0.5 * 27.2114 * (
                directpart + recippartreal - selfint - surfterm)
    else:
        #returns a potential, not an energy (no multiplying by -q/2)
        totalPC = (directpart + recippartreal - selfint) * 27.2114

    if not silence:
        print ('Final PC Energy term is then ', totalPC, ' (eV)')

    return totalPC, optgamma


def disttrans(s1, s2, c):
    """
    this is function for calculating distance to each atom and finding NGX grid pts at each atom
    input pure Locpot Object(s1), defect Locpot Object(s2)
        and zero-definedindex of where defect is in POSCAR (c)
    should amend to calculate potential out to radius of atom for averaging LocpotData?
    """
    struct = s1.structure
    defstruct = s2.structure
    griddict = {}  # dictionary with indices keys in order of structure list

    def getgridind(coord):
        # get grid index (if wanted to do a range then
        #      s1.structure.sites[i].specie.atomic_radius gives radius...
        grdind = []
        for j in range(3):
            axis = s1.get_axis_grid(j)
            ind = len(axis) + 1  # so function will fail if index not found
            diff = 10.0
            for k in range(len(axis)):
                if abs(coord[j] - axis[k]) < diff:  # minimize distance from pt to grid point
                    diff = abs(coord[j] - axis[k])
                    ind = k
            grdind.append(ind)
        return grdind

    def_fcoord = struct.sites[c].frac_coords
    def_ccoord = struct.sites[c].coords
    for i in range(len(struct.species)):
        if i != c:
            cart_coord = struct.sites[i].coords
            frac_coord = struct.sites[i].frac_coords
            dist, img = struct.lattice.get_distance_and_image(def_fcoord, frac_coord)
            cart_reldef = np.dot((frac_coord + img), struct.lattice._matrix) - def_ccoord
            if abs(norm(cart_reldef) - dist) > 0.001:
                print 'image locater issue encountered for site=', i, \
                    ' distance should be ', dist, ' but calculated to be ', norm(cart_reldef)
                return
            griddict[i] = {'dist': dist, 'cart': cart_coord, 'cart_reldef': cart_reldef,
                           # cart_reldef=cartesian coord with defect as origin
                           'bulkgrid': getgridind(
                               cart_coord)}  # should try and make this a range of grid indices in the future?

    # now get defect grid stuff
    if len(defstruct.species) == len(struct.species):  # antisites, subs
        for j in range(len(struct.species)):
            if j == c:
                continue
            griddict[j]['defind'] = j
            griddict[j]['defgrid'] = getgridind(defstruct.sites[j].coords)
    elif (len(defstruct.species) + 1) == len(struct.species):  # vacancies
        for j in range(len(struct.species)):
            if j == c:
                continue
            elif j > c:
                defind = j - 1
            else:
                defind = j
            griddict[j]['defind'] = defind
            griddict[j]['defgrid'] = getgridind(defstruct.sites[defind].coords)
    else:
        print 'defect type must be interstitial or something that this code isnt programmed to do yet.'
        return

    return griddict


def wigner_seitz_radius(s):
    """
    Calculate the Wigner Seitz radius for the given structure.
    Args:
        s: Either structure or VolumetricData object
    """
    try:
        lat = Lattice(s.structure.lattice_vectors())
    except:
        lat = Lattice(s.lattice_vectors())

    wz = lat.get_wigner_seitz_cell()
    # wz is list of WS cell face vertices
    # make list of midpoints on faces of WS cell
    dist = []
    for facet in wz:
        midpt = np.mean(np.array(facet), axis=0)
        dist.append(norm(midpt))
    wsrad = min(dist)
    return wsrad


# Add Specific Tool for Madelung Correction?


# #reference sxdefectalign call for this example:
# '~/sxdefectalign --vasp -a1 --relative --pos 0.0,0.0,0.0 --charge 3 --ecut 38.2192757 --eps 18.099 --vref LOCPOT_vref --vdef LOCPOT_vdef'

class ChargeCorrection(object):
    def __init__(self, axis, dielectric_tensor, pure_locpot_path, 
            defect_locpot_path, q, pos, energy_cutoff=520, madetol=0.0001, 
            silence=False, q_model=None):
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
            pos: Position of the defect in the cell (Frac or coord?)
            energy_cutoff: Energy for plane wave cutoff (in eV).
                 If not given, Materials Project default 520 eV is used.
            madetol: Tolerance (double or float)
            silence: Flag for disabling/enabling  messages (Bool)
            q_model (QModel object): User defined charge for correction.
                 If not given, highly localized charge is assumed.
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
        self._silence = silence  #for silencing printflags
        if not q_model:
            self._q_model = QModel()

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

        #lattice constants from pure cell
        def cleanlat(dat):
            return norm(dat[0]), norm(dat[1]), norm(dat[2])

        ap = s1.lattice.get_cartesian_coords(1)
        if self._silence == False:
            print 'run Freysoldt 2011 PC calculation (should be equivalent to sxdefectalign)'
            print 'defect lattice constants are (in angstroms)' + str(cleanlat(ap))
        [a1, a2, a3] = 1.8897 * ap
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
            print 'difference in (eV) is ' + str(round((eper - eiso) * 27.2114, 4))  #27.2114 eV/1 hartree
        PCfreycorr = round(((eiso - eper) / self._dielectricconst) * 27.2114, 6)  #converted to eV
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
        v_R = np.real(v_R)* 27.2114
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

    def kumagai_correction(self, title=None, vb=None, vd=None):
        #runs correction. If you want a plot of potential averaging process set title to name of defect
        #vb and vd are preloaded locpot objects for speeding this up.
        if not self._silence:
            print 'This is Kumagai Correction.'
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
        print '\nRun PC energy'
        energy_pc, optgam = self.kumagai_pc(vb.structure)
        print '\nPC calc done, correction =', round(energy_pc, 4), \
                ' optimized gamma found to be: ', round(optgam, 4)
        print 'Now run potenttial alignment script'
        potalign = self.kumagai_potalign(vb, vd, optgam=optgam, title=title)
        print '\n\nAlright so the corrections are:'
        print 'PCenergy = ', round(energy_pc, 5), '  potential alignment = ', round(potalign, 5)
        print 'TOTAL Kumagai correction = ', round(energy_pc - potalign, 5)
        return round(energy_pc - potalign, 5)

    def kumagai_pc(self, s1=None):
        #note that this ony needs structure info, not locpot info;
        # so s1=structure object speeds this calculation up alot
        if not s1:
            print 'load structure from Pure locpot'
            tmp = Locpot.from_file(self._purelocpot)
            s1 = tmp.structure
        print 'run Kumagai PC calculation'

        #get aniso PCenergy (equation 8 from kumagai paper)
        energy_pc, optgamma = get_anistropic_pc_energy(
                s1, self._dieltens,  self._q, 
                self._madetol, [0., 0., 0.], self._silence)  #returns PCenergy in eV

        if not self._silence:
            print 'PC energy determined to be ', energy_pc, ' eV (', \
                    energy_pc/27.2114, ' Hartree)'  #27.2114 eV/1 hartree
            print 'Optimized Gamma found to be ' + str(optgamma)

        return energy_pc, optgamma  #PC energy in eV

    def kumagai_potalign(self, v1=None, v2=None, optgam=None, title=None):
        """
        Potential alignment for Kumagai method
        Args:
            v1: Bulk locpot object
            v2: Defect locpot object
            optgam: None if gamma has not been optimized yet, otherwise feed it optimized value
            title: Title for the plot. None will not generate the plot
        """
        if not self._silence:
            print ('run Kumagai potential calculation (atomic site averaging)')
        if not v1:
            print 'load pure locpot object'
            v1 = Locpot.from_file(self._purelocpot)
        if not v2:
            print 'load defect locpot object'
            v2 = Locpot.from_file(self._deflocpot)

        angset, [a1, a2, a3], vol, determ, invdiel = kumagai_init(
                v1.structure, self._dieltens, sil=self._silence)

        #find location of defect
        from pymatgen.analysis.structure_matcher import StructureMatcher as SM

        smatch = SM(primitive_cell=False, allow_subset=True)
        matches = smatch.get_transformation(v1.structure, v2.structure)
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

        if defect_index == -1:
            print 'Could not find defect index! (problem with structure id-ing)'
            return

        #average potential part
        # create distance matrix for plotting
        # get a list of atoms outside of WS radius.
        potinddict = disttrans(v1, v2, defect_index)  #this is to calculate distance matrix for plotting
        wsrad = wigner_seitz_radius(v1)
        for i in potinddict.keys():
            if potinddict[i]['dist'] > wsrad:
                potinddict[i]['OutsideWS'] = True
            else:
                potinddict[i]['OutsideWS'] = False

        #calculate potential value at each atomic site
        puredat = v1.data["total"]
        defdat = v2.data["total"]
        jup = 0
        for i in potinddict.keys():
            jup += 1
            if (not title and not potinddict[i]['OutsideWS']):
                #dont need to calculate inside WS if not printing plot
                continue
            if not self._silence:
                print '-------------------------------------'
                print "calculate alignment potential data for atom " + str(i)
            dx, dy, dz = potinddict[i]['defgrid']
            bx, by, bz = potinddict[i]['bulkgrid']
            v_qb = defdat[dx][dy][dz] - puredat[bx][by][bz]  #should change this to averaging pure
            # and def within a range then subtract

            v_pc, gam1 = get_anistropic_pc_energy(v1.structure, self._dieltens, 
                    self._q, self._madetol, potinddict[i]['cart_reldef'], 
                    silence=True, optgam=optgam)
            potinddict[i]['Vpc'] = v_pc
            potinddict[i]['Vqb'] = v_qb
            if not self._silence:
                print 'Has anisotropic point charge energy = ', v_pc
                print 'DFT bulk/defect difference = ', v_qb
                print 'atoms left to calculate = ' + str(len(potinddict.keys()) - jup) #this is broken if you arent plotting
        if not self._silence:
            print '--------------------------------------'

        #now parse and plot if neccessary
        if title:  #to make shading region prettier
            fullspecset = v1.structure.species
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
            x = np.arange(wsrad, max(v1.structure.lattice.abc), 0.01)
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

    def madelung_corr(self, s1=None):
        #NOT DONE

        if not s1:
            print 'load pure file'
            s1 = Locpot.from_file(self._purelocpot)

            #could make this run similar to PCenergy function for Kumagai but change anisotropic terms to isotropic?

            ##or could do a pymatgen trick like:
            #from pymatgen.analysis.ewald import EwaldSummation as ES
            #ts=struct(s1.lattice,[s1.species[0]],[[0,0,0]])
            #ts.add_oxidation_state_by_site([self._q])
            #val=ES(ts)
            #MCen=(val.total_energy)/self._dielectricconst


            #s1prim = s1.structure.get_primitive_structure()
            #L = s1prim.lattice.volume ** 1 / 3.
            #gam = 5. / L
            #madelung = madelungconst(s1prim, gam)
            #Mcorr = self._q ** 2 * madelung / (2 * self._dielectricconst * L)
            #print 'Madelung leading order ES correction is then (units?): ' + Mcorr
            #return Mcorr

    def mp_corr(self, s1=None):
        """
        Markov Payne correction. 
        Compute alpha for each structure or use the given ne
        """
        #NOT DONE
        pass


if __name__ == '__main__':
    s = ChargeCorrection(
            0, 18.099, 
            '../../Gavacm3testMachgcorr/LOCPOT_vref',
            '../../Gavacm3testMachgcorr/LOCPOT_vdef',
            -3, silence=False)

    #print 'load locpots'
    #vb=Locpot.from_file('../../Gavacm3testMachgcorr/LOCPOT_vref')
    #vb=Locpot.from_file('trans/LOCPOT_vref')
    #vd=Locpot.from_file('../../Gavacm3testMachgcorr/LOCPOT_vdef')


    s.freysoldt_correction(title='Test')
    #s.freysoldt_pc()
    #s.freysoldt_pc(vb.structure)
    #s.freysoldt_potalign(v1=vb,v2=vd,title='Test')

    #s.kumagai_correction(title='Test')
    #s.KumagaiPC(vb.structure)
    #s.Kumagaipotalign(vb,vd,optgam=1.848916,printflag='Gavac+3test')

