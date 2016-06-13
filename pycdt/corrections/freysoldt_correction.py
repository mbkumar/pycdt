"""
This module is Freysoldt correction for isotropic systems
1) Freysoldt correction for isotropic systems.
includes
   a) PC energy
   b) potential alignment by planar averaging.
If you use the corrections implemented in this module, cite
   Freysoldt, Neugebauer, and Van de Walle, Phys. Rev. Lett. 102, 016402 (2009)
   [Optionally Phys. Status Solidi B. 248, 1067-1076 (2011) ]
   in addition to the pycdt paper
"""
__author__ = 'Danny Broberg, Bharat Medasani'
__email__ = 'dbroberg@gmail.com, mbkumar@gmail.com'

import sys
import math

import numpy as np

from pymatgen.io.vasp.outputs import Locpot
from pymatgen.core.structure import Structure

norm = np.linalg.norm

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
    """
    Args:
        a1, a2, a3: lattice vectors in bohr
        encut: energy cut off in eV
    Returns:
        reciprocal lattice vectors with energy less than encut
    """

    # define recip vectors first, (units of 1/angstrom).
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

    return recip

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

def closestsites(struct_blk, struct_def, pos):
    #input bulk and defect structures and get site that is nearest to the (cartesian) input position
    bulkclosesites = struct_blk.get_sites_in_sphere(pos, 5)
    bulkclosesites.sort(key=lambda x:x[1])
    defclosesites = struct_def.get_sites_in_sphere(pos, 5)
    defclosesites.sort(key=lambda x:x[1])

    return bulkclosesites[0],defclosesites[0] #returns closest (site object, dist) for both bulk and defect

def find_defect_pos(struct_blk, struct_def):
    """
    output cartesian coords of defect in bulk,defect cells.

    If vacancy defectpos=None, if interstitial bulkpos=None, if antisite/sub then both defined
    """
    if len(struct_blk.sites) > len(struct_def.sites):
        vactype = True
        interstittype = False
    elif len(struct_blk.sites) < len(struct_def.sites):
        vactype = False
        interstittype = True
    else:
        vactype = False
        interstittype = False

    sitematching = []
    for site in struct_blk.sites:
        blksite, defsite = closestsites(struct_blk, struct_def, site.coords)
        if vactype and blksite[0].specie.symbol != defsite[0].specie.symbol:
            return blksite[0].coords, None
        elif interstittype and blksite[0].specie.symbol != defsite[0].specie.symbol:
            return None, defsite[0].coords
        elif blksite[0].specie.symbol != defsite[0].specie.symbol: #subs or antisite type
            return blksite[0].coords, defsite[0].coords
        sitematching.append([blksite[0],blksite[1],defsite[0],defsite[1]])

    if vactype: #just in case site type is same for closest site to vacancy
        sitematching.sort(key=lambda x:x[3])
        vacant = sitematching[-1]
        return vacant[0].coords, None
    elif interstittype: #just in case site type is same for closest site to interstit
        sitematching.sort(key=lambda x:x[1])
        interstit = sitematching[-1]
        return  None, interstit[2].coords

    return None,None #if you get here there is an error


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


class FreysoldtCorrection(object):
    def __init__(self, axis, dielectricconst, pure_locpot_path,
            defect_locpot_path, q, energy_cutoff=520,
            madetol=0.0001, silence=False, q_model=None):
        """
        Args:
            axis: axis to do Freysoldt averaging over (zero-defined). Has no effect on
                 Kumagai correction, so better move to Freysoldt potalign
            dielectric_tensor: Macroscopic dielectric tensor 
                 Include ionic also if defect is relaxed, othewise ion clamped.
                 Can be a matrix array or scalar.
            pure_locpot_path: Bulk Locpot file path OR locpot object
            defect_locpot_path: Defect Locpot file path OR locpot object
            q: Charge associated with the defect (not of the homogen. background). Typically integer
            energy_cutoff: Energy for plane wave cutoff (in eV).
                 If not given, Materials Project default 520 eV is used.
            madetol: Tolerance for convergence of energy terms in eV (double or float)
            silence: Flag for disabling/enabling  messages (Bool)
            q_model (QModel object): User defined charge for correction.
                 If not given, highly localized charge is assumed.
        """
        self._axis = axis
        if isinstance(dielectricconst, int) or \
                isinstance(dielectricconst, float):
            self._dielectricconst = float(dielectricconst)
        else:
            self._dielectricconst = float(np.mean(np.diag(dielectricconst)))
        self._purelocpot = pure_locpot_path
        self._deflocpot = defect_locpot_path
        self._madetol = madetol
        self._q = q
        self._encut = energy_cutoff
        self._pos = None #code will determine positions of defect in bulk cell
        self._defpos = None #code will determine defect position in defect cell (after relaxation)
        self._silence = silence
        if not q_model:
            self._q_model = QModel()

    def correction(self, title=None, partflag='All'):
        """
        Args:
            title: set if you want to plot the planar averaged potential
            partflag: four options
                'pc' for just point charge correction, or
               'potalign' for just potalign correction, or
               'All' for both, or
               'AllSplit' for individual parts split up (form [PC,potterm,full])
        """
        if not self._silence:
            print 'This is Freysoldt Correction.'
        if not self._q:
            if partflag=='AllSplit':
                return [0.0,0.0,0.0]
            else:
                return 0.0

        if not type(self._purelocpot) is Locpot:
            if not self._silence:
                print 'Load bulk locpot'
            self._purelocpot = Locpot.from_file(self._purelocpot)

        if not self._silence:
            print '\nRun PC energy'
        if partflag!='potalign':
            energy_pc = self.pc()
            if not self._silence:
                print '\nPC calc done, correction =', round(energy_pc, 4)
                print 'Now run potenttial alignment script'

        if partflag!='pc':
            if not type(self._deflocpot) is Locpot:
                if not self._silence:
                    print 'Load defect locpot'
                self._deflocpot = Locpot.from_file(self._deflocpot)
            potalign = self.potalign(title=title)

        if not self._silence:
            print '\n\nFreysoldt Correction details:'
            if partflag!='potalign':
                print 'PCenergy (E_lat) = ', round(energy_pc, 5)
            if partflag!='pc':
                print 'potential alignment (-q*delta V) = ', round(potalign, 5)
            if partflag in ['All','AllSplit']:
                print 'TOTAL Freysoldt correction = ', round(energy_pc + potalign, 5)

        if partflag=='pc':
            return round(energy_pc,5)
        elif partflag=='potalign':
            return round(potalign,5)
        elif partflag=='All':
            return round(energy_pc+potalign,5)
        else:
            return [round(energy_pc,5),round(potalign,5),round(energy_pc+potalign,5)]

    def pc(self,struct=None):
        """
        Peform Electrostatic Correction
        note this ony needs structural info
        so struct input object speeds this calculation up
        equivalently fast if input Locpot is a locpot object
        """
        if type(struct) is Structure:
            s1=struct
        else:
            if not type(self._purelocpot) is Locpot:
                if not self._silence:
                    print 'load Pure locpot'
                self._purelocpot = Locpot.from_file(self._purelocpot)
            s1=self._purelocpot.structure

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
        encut1 = 20.  #converge to some smaller encut first [eV]
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
            print 'difference in (eV) is ' + str(round((eper - eiso) * hart_to_ev, 4))
        PCfreycorr = round(((eiso - eper) / self._dielectricconst) * hart_to_ev, 6)  #converted to eV
        if self._silence == False:
            print 'Defect Correction without alignment (eV): ', PCfreycorr

        return PCfreycorr

    def potalign(self, title=None,  widthsample=1., axis=None):
        """
        For performing planar averaging potential alignment

        Accounts for defects in arbitrary positions
        title is for name of plot, if you dont want a plot then leave it as None
        widthsample is the width of the region in between defects where the potential alignment correction is averaged
        axis allows you to override the axis setting of class
                (good for quickly plotting multiple axes without having to reload Locpot)
        """
        if not axis:
            axis=self._axis
        else:
            axis=axis

        if not type(self._purelocpot) is Locpot:
            if not self._silence:
                print 'load pure locpot object'
            self._purelocpot = Locpot.from_file(self._purelocpot)
        if not type(self._deflocpot) is Locpot:
            if not self._silence:
                print 'load defect locpot object'
            self._deflocpot = Locpot.from_file(self._deflocpot)

        #determine location of defects
        blksite,defsite=find_defect_pos(self._purelocpot.structure,self._deflocpot.structure)
        if blksite is None and defsite is None:
            print 'Error. Not able to determine defect site...'
            return
        if not self._silence:
            if blksite is None:
                print 'Found defect to be Interstitial type at ',defsite
            elif defsite is None:
                print 'Found defect to be Vacancy type at ',blksite
            else:
                print 'Found defect to be antisite/substitution type at ',blksite,' in bulk, and ',\
                        defsite,' in defect cell'

        #It is important to do planar averaging at same position, otherwise
        #you can get rigid shifts due to atomic changes at far away from defect
        #note these are cartesian co-ordinate sites...
        if defsite is None: #vacancies
            #self._defpos=blksite
            self._pos=blksite
        else: #all else, do w.r.t defect site
            #self._defpos=defsite
            self._pos=defsite

        ind = []
        for i in range(3):
            if axis == i:
                continue
            else:
                ind.append(i)

        x = np.array(self._purelocpot.get_axis_grid(axis))  #angstrom
        nx = len(x)
        print 'run Freysoldt potential alignment method'

        #perform potential alignment part
        pureavg = self._purelocpot.get_average_along_axis(axis)  #eV
        defavg = self._deflocpot.get_average_along_axis(axis)  #eV

        #now shift these planar averages to have defect at origin
        blklat=self._purelocpot.structure.lattice
        #deflat=self._deflocpot.structure.lattice
        axfracval=blklat.get_fractional_coords(self._pos)[axis]
        #axdefval=deflat.get_fractional_coords(self._defpos)[axis]
        axbulkval=axfracval*blklat.abc[axis]
        #axdefval*=deflat.abc[axis]
        if axbulkval<0:
            axbulkval += blklat.abc[axis]
        elif axbulkval > blklat.abc[axis]:
            axbulkval -= blklat.abc[axis]

        if axbulkval:
            for i in range(len(x)):
                if axbulkval<x[i]:
                    break
            rollind=len(x)-i
            pureavg=np.roll(pureavg,rollind)
            defavg = np.roll(defavg,rollind)
        # if axdefval:
        #     for i in range(len(x)):
        #         if axdefval<x[i]:
        #             break
        #     rollind=len(x)-i
        #     defavg=np.roll(defavg,rollind)


        if not self._silence:
            print 'calculating lr part along planar avg axis'
        latt = self._purelocpot.structure.lattice
        reci_latt = latt.reciprocal_lattice
        dg = reci_latt.abc[axis]
        dg/=ang_to_bohr #convert to bohr to do calculation in atomic units

        v_G = np.empty(len(x), np.dtype('c16'))
        epsilon = self._dielectricconst
        # q needs to be that of the back ground
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
        v_R /= (latt.volume*ang_to_bohr**3)
        v_R = np.real(v_R)* hart_to_ev

        max_imag_vr = v_R_imag.max()
        if abs(max_imag_vr) > self._madetol:
            print 'imaginary part found to be ', max_imag_vr, ' this is an issue'
            sys.exit()

        #now get correction and do plots
        short = (defavg - pureavg - v_R)
        checkdis = int((widthsample / 2) / (x[1]-x[0]))
        mid = len(short) / 2

        tmppot = [short[i] for i in range(mid - checkdis, mid + checkdis)]
        if not self._silence:
            print('shifted defect position on axis (',axbulkval,') to origin')
            print('means sampling region is (', x[mid-checkdis], ',', x[mid+checkdis], ')')

        C = -np.mean(tmppot)
        print 'C=',C
        finalshift = [short[j] + C for j in range(len(v_R))]
        v_R = [v_R[j]-C for j in range(len(v_R))]

        if not self._silence:
           print 'C value is averaged to be ' + str(C) + ' eV, '
           print 'Pot. align correction (-q*delta V) is then (eV) : ' + str(-float(self._q) * float(C))
        if title:
            if title!='written':
                import matplotlib.pyplot as plt
                plt.figure(1)
                plt.clf()
                plt.plot(x, v_R, c="green", zorder=1, label="long range from model")
                plt.plot(x, defavg - pureavg, c="red", label="DFT locpot diff")
                plt.plot(x, finalshift, c="blue", label="short range (aligned)")
                tmpx=[x[i] for i in range(mid - checkdis, mid + checkdis)]
                plt.fill_between(tmpx, -100, 100, facecolor='red', alpha=0.15, label='sampling region')
                plt.xlim(np.floor(x[0]),np.ceil(x[-1]))
                ymin=min(min(v_R),min(defavg - pureavg),min(finalshift))
                ymax=max(max(v_R),max(defavg - pureavg),max(finalshift))
                plt.ylim(np.floor(ymin),np.ceil(ymax))
                plt.xlabel('planar average along axis ' + str(axis+1)+' (Angstrom)')
                plt.ylabel('Potential (V)')
                plt.legend(loc=9)
                plt.axhline(y=0, linewidth=0.2, color='black')
                plt.title(str(title) + ' planar averaged electrostatic potential')
                plt.xlim(0,max(x))
                plt.savefig(str(title)+'FreyplnravgPlot.pdf')
            else:
                #this is because current uploading format for written is bad for numpy arrays...
                #Might want to update this in future so that dumping format is smarter than just dumping/loading a string
                xtmp=[]
                v_Rtmp=[]
                DFTdifftmp=[]
                finalshifttmp=[]
                for j in range(len(v_R)):
                    xtmp.append(x[j])
                    v_Rtmp.append(v_R[j])
                    DFTdifftmp.append(defavg[j]-pureavg[j])
                    finalshifttmp.append(finalshift[j])

                forplotting={'x':xtmp,'v_R':v_Rtmp,'DFTdiff':DFTdifftmp,'finalshift':finalshifttmp,'checkrange':[mid - checkdis,mid + checkdis]}
                fname='FreyAxisData.dat'
                with open(fname,'w') as f:
                    f.write(str(forplotting))

        return -float(self._q)*C  #pot align energy correction (eV), add to energy output of PCfrey

    def plot_from_datfile(self,name='FreyAxisData.dat',title='default'):
        """
        Takes data file called 'name' and does plotting.
        Good for later plotting of locpot data after running run_correction()
        """
        import ast
        with open(name,'r') as f:
            plotvals=f.read()
        plotvals=ast.literal_eval(plotvals) #converting string to dictionary

        x=plotvals['x']
        v_R=plotvals['v_R']
        DFTdiff=plotvals['DFTdiff']
        finalshift=plotvals['finalshift']
        check=plotvals['checkrange']

        import matplotlib.pyplot as plt
        plt.figure()
        plt.clf()
        plt.plot(plotvals['x'], plotvals['v_R'], c="green", zorder=1, label="long range from model")
        plt.plot(x, DFTdiff, c="red", label="DFT locpot diff")
        plt.plot(x, finalshift, c="blue", label="short range (aligned)")
        tmpx=[x[i] for i in range(check[0], check[1])]
        plt.fill_between(tmpx, -100, 100, facecolor='red', alpha=0.15, label='sampling region')
        plt.xlim(round(x[0]),round(x[-1]))
        ymin=min(min(v_R),min(DFTdiff),min(finalshift))
        ymax=max(max(v_R),max(DFTdiff),max(finalshift))
        plt.ylim(-0.2+ymin,0.2+ymax)
        plt.xlabel('planar average along axis ' + str(1))
        plt.ylabel('Potential')
        plt.legend(loc=9)
        plt.axhline(y=0, linewidth=0.2, color='black')
        plt.title(str(title) + ' defect potential')
        plt.xlim(0,max(x))
        plt.savefig(str(title)+'FreyplnravgPlot.png')



if __name__ == '__main__':
    s = FreysoldtCorrection(0,11.814,'../../bulk/LOCPOT', 'LOCPOT', 1)
    #s.correction(title='written',partflag='AllSplit')
    s.plot_from_datfile()
    # s = FreysoldtCorrection(0, 18.099, '../../../../Gavacm3testMachgcorr/LOCPOT_vref',
    #         '../../../../Gavacm3testMachgcorr/LOCPOT_vdef', -3)
    #s.correction(title='ThisIstest')


