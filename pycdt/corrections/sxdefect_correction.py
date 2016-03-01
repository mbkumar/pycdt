#!/usr/bin/env python

"""
This module is wraparound to sxdefectalign code by Freysoldt
Applies Freysoldt correction with planar averaged potential method
"""

__author__ = "Geoffroy Hautier, Danny Broberg"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Geoffroy Hautier/Danny Broberg"
__email__ = "geoffroy@uclouvain.be, dbroberg@berkeley.edu"
__status__ = "Development"
__date__ = "April 24, 2015"

import subprocess
import os
import numpy as np

from pymatgen.io.vasp.outputs import Locpot
from monty.tempfile import ScratchDir


class SxDefectCorrection(object):
    """
    This class applies the Freysoldt correction to remove electrostatic defect
    interaction contribution to energy and apply potential alignment.
    Ideally this sxdefectalign wrapper class should be replaced by a python
    code implementing a generalized anisotropic Freysoldt correction
    """
    
    def __init__(self, locpot_bulk, locpot_defect, charge, epsilon,
                 site_frac_coords, encut, axiscalcs=[0], lengths=None, name=''):
        """
        Args:
            locpot_bulk: 
                Location of LOCPOT of bulk (or actual Locpot object)
            locpot_defect: 
                Location of LOCPOT of defect  (or actual Locpot object)
            charge: 
                Charge relative to neutral defect (not relative to bulk)
            epsilon: 
                Dielectric constant obtained from relaxation run
            site_frac_coords: 
                Fractional coordinates of defect site as list
            encut: 
                Planewave basis energy cutoff used in VASP run (in eV)
            name:
                Name of the defect for writing/plotting files
            lengths:
                Length of lattice vectors.
        """
        self._locpot_bulk = os.path.abspath(locpot_bulk)
        self._locpot_defect = os.path.abspath(locpot_defect)
        self._charge = charge 
        self._epsilon = epsilon 
        self._frac_coords = site_frac_coords   
        self._encut = encut
        self._axiscalcs = axiscalcs
        if not lengths:
            if type(self._locpot_bulk) is Locpot:
                self._lengths=self._locpot_bulk.structure.lattice.abc
            else:
                struct=Locpot.from_file(locpot_bulk)
                self._lengths=struct.structure.lattice.abc
                print('had to import lengths, if want to speed up set lengths='+str(self._lengths))
        else:
            self._lengths = lengths
        
    def prepare_files(self):
        if  self._charge==0:
            print('defect has charge 0, so freysoldt correction is 0')
            return
        #self._locpot_bulk.write_file("LOCPOT_vref",True)   #from when the input is a locpot object
        path, name = os.path.split(self._locpot_bulk)
        mod_blk_locpot = os.path.join(path, name+'_vref')
        self.mod_bulk_locpot = mod_blk_locpot
        if not os.path.exists(mod_blk_locpot):
            print('prep pure Locpot')
            with open(self._locpot_bulk) as input:
                with open(mod_blk_locpot,'w') as output:
                    cnt = 0
                    for line in input:
                        cnt += 1
                        if cnt != 6:
                            output.write(line)
            #cmd="perl -n -e 'print if $. != 6' "+str(self._locpot_bulk)+" > "+mod_blk_locpot
            #print cmd
            #os.system(cmd)
        #self._locpot_defect.write_file("LOCPOT_vdef",True) #from when these are locpot objects
        path, name = os.path.split(self._locpot_defect)
        mod_defect_locpot = os.path.join(path, name+'_vdef')
        self.mod_defect_locpot = mod_defect_locpot
        if not os.path.exists(mod_defect_locpot):
            print('prep defect Locpot')
            with open(self._locpot_defect) as input:
                with open(mod_defect_locpot,'w') as output:
                    cnt = 0
                    for line in input:
                        cnt += 1
                        if cnt != 6:
                            output.write(line)
            #md="perl -n -e 'print if $. != 6' "+str(self._locpot_defect)+" > "+mod_def_locpot
            #rint cmd
            #s.system(cmd)
        print('locpots prepared for sxdefectalign')

    def plot_hartree_pot(self):
        #plot planar averages of bulk and defect (good for seeing global changes)
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials')
        bulkloc=Locpot.from_file(self._locpot_bulk)
        defloc=Locpot.from_file(self._locpot_defect)
        get_agrid = bulkloc.get_axis_grid
        get_baavg = bulkloc.get_average_along_axis
        get_daavg = defloc.get_average_along_axis
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            latt_len = self._lengths[axis]
            ax.plot(get_agrid(axis),get_baavg(axis),'r',
                    label="Bulk potential")
            ax.plot(get_agrid(axis), get_daavg(axis),'b',
                    label="Defect potential")
            ax.plot([self._frac_coords[axis]*latt_len], [0], 'or', 
                    markersize=4.0, label="Defect site")
            ax.set_ylabel("axis "+str(axis+1))
            if not axis:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        #plt.savefig('locpotavgplot.png')
        plt.show()

    def plot_hartree_pot_diff(self):
        #only plot the difference in planar averaged potentials
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potential difference')
        bulkloc=Locpot.from_file(self._locpot_bulk)
        defloc=Locpot.from_file(self._locpot_defect)
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            defect_axis = defloc.get_axis_grid(axis)
            defect_pot = defloc.get_average_along_axis(axis)
            pure_pot = bulkloc.get_average_along_axis(axis)
            latt_len = self._lengths[axis]
            ax.plot(defect_axis,defect_pot-pure_pot,'b',
                    label='Defect-Bulk difference')
            ax.plot([self._frac_coords[axis] * latt_len], [0], 'or',
                    markersize=4.0, label='Defect site')
            ax.set_ylabel("axis "+str(axis+1))
            if not axis:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        #plt.savefig('locpotdiffplot.png')
        plt.show()

    def plot_all_hartree_pot(self):
        #plot planar averaged locpots, along with the difference for each axis
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials and difference')
        bulkloc=Locpot.from_file(self._locpot_bulk)
        defloc=Locpot.from_file(self._locpot_defect)
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            defect_axis = defloc.get_axis_grid(axis)
            defect_pot = defloc.get_average_along_axis(axis)
            pure_pot = bulkloc.get_average_along_axis(axis)
            ax.plot(defect_axis, pure_pot, 'r', label="Bulk potential")
            ax.plot(defect_axis, defect_pot, 'b', label="Defect potential")
            ax.plot(defect_axis, defect_pot-pure_pot, 'k',
                    label='Defect-Bulk difference')
            ax.plot([self._frac_coords[axis]*self._lengths[axis]],\
                    [0], 'or', markersize=4.0, label='Defect site')
            ax.set_ylabel("axis "+str(axis+1))
            if not axis:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        #plt.savefig('locpotavgdiffplot.png')
        plt.show()

    def plot_pot_diff(self, align=None, print_pot_flag='written',name='default',widthsample=0.5):
        """
        Runs sxdefectalign and obtains three alignment constants for 
        each axis. To obtain Freysoldt correction call this function twice 
        (once to get alignment, and second time to get final correction 
        and plots). In the first call, use print_pot_flag=none.  
        Then in second call get corrections and leave vline-eV.dat 
        files for each axis and produce potential alignment correction
        or if  print flag= written then get different parts of vline-eV.dat 
        written to files for plotting (good for nersc or where you can't 
        plot easy) or if print flag=plotfull then show the full matplotlib 
        plot for all three axes at end of calculation. 
        encut is eV cutoff from VASP

        Would like final workflow to include the flag that I have put here 
        for when planar average varies by more than 0.2 eV around far region

        name is what to call the plot or written file...
        """
        if align is None:
            align=[0. for i in range(len(self._axiscalcs))]

        if not self._charge: #don't need charge correction if charge is zero
            zrochg=[]
            for i in range(2):
                zrochg.append([0 for j in  range(len(self._axiscalcs))])
            return zrochg # in form [[0,0,0],[0,0,0]] for totalPC and pot align terms on each axis
        #correction from output (should include alignment once alignment 
        # has been done)
        result = []   
        platy = []    #alignment terms for each axis
        # if want to plot right here, then build dictionary for storing 
        # planar average values of each axis
        if print_pot_flag == 'plotfull' or 'written':
            plotvals = {str(i):{} for i in self._axiscalcs}
        for axis in self._axiscalcs:
            print('do axis '+str(axis+1))
            #print self._frac_coords[1:]
            #relpos = (str(self._frac_coords)[1:])[:-1]
            #relpos = relpos.replace(" ",",")
            relpos = ",".join(str(i) for i in self._frac_coords)
            command = ['~/bin/sxdefectalign', '--vasp', '-a'+str(axis+1), 
                    '--relative', '--pos', relpos, 
                    '--charge', str(-self._charge), 
                    '--ecut', str(self._encut/13.6057), #eV to Ry for sxdefect 
                    '--eps', str(self._epsilon), 
                    '-C', str(-float(align[axis])), 
                    '--vref', self.mod_bulk_locpot,
                    '--vdef', self.mod_defect_locpot]
            print(command)

            #standard way of running NERSC commands.
            #in case NERSC (Hopper) has issues with python subprocess can use hack
            #p = subprocess.Popen(command, stdout=subprocess.PIPE, 
            #        stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            #output, err = p.communicate()
            #out = out.decode('utf-8')
            #err = err.decode('utf-8')
            #print 'output from sxdefectalign = ', str(out)
            #result.append(float(out[0].split("\n")[12].split()[3]))
            #print "chg correction is "+str(result[-1])

            #this is hack wrap-around for when subprocess doesn't work
            #(which is always an issue on hopper now...)
            cmd = ' '.join(command)
            cmd += ' > tmpoutput'
            os.system(cmd)
            
            with open('tmpoutput') as f:
                output = f.readlines()
            #print 'output from sxdefectalign = '+str(output)
            val =  output[-1].split()[3].strip()
            #result.append(float(output[-1].split()[3]))
            result.append(float(val))
            print("chg correction is "+str(result[-1]))
            os.remove('tmpoutput')

            x_lr, y_lr = [], []
            x, y = [], []
            x_diff, y_diff = [], []
            with open("vline-eV.dat",'r') as f_sr: #read in potential 
                for r in f_sr:
                    tmp = r.split("\t")
                    if(len(tmp)<3 and not r.startswith("&")):
                       x_lr.append(float(tmp[0])/1.889725989)   # to Angstrom
                       y_lr.append(float(tmp[1])) 
                    if len(tmp) > 2:
                        x.append(float(tmp[0])/1.889725989)     # to Angstrom
                        x_diff.append(float(tmp[0])/1.889725989)# to Angstrom
                        y.append(float(tmp[2].rstrip("\n")))
                        y_diff.append(float(tmp[1]))

            if print_pot_flag == 'none':
                if os.path.exists("vline-eV.dat"):
                    os.rename("vline-eV.dat","axis"+str(axis)+"vline-eV.dat")
            else:
                if os.path.exists("vline-eV.dat"):
                    os.remove("vline-eV.dat")

            # Extract potential alignment term averaging window of +/- 1 Ang 
            # around point halfway between neighboring defects
            latt_len = self._lengths[axis]
            if self._frac_coords[axis] >= 0.5:
                 platx = (self._frac_coords[axis]-0.5) * latt_len
            else:
                 platx = (self._frac_coords[axis]+0.5) * latt_len
            print("half way between defects is: ", platx)

            #sample short range potential in width of size widthsample
            radsample=widthsample/2.
            xmin = latt_len - (radsample-platx) if platx < radsample else platx-radsample
            xmax = radsample-(latt_len-platx) if platx > latt_len-radsample else radsample+platx
            print('means sampling region is (', xmin, ',', xmax, ')')

            tmpalign=[]
            if xmax < xmin:
                print('wrap around detected, special alignment needed')
                for i in range(len(x)):
                    if x[i]<xmax or x[i]>xmin:
                        tmpalign.append(y[i])
                    else:
                        continue
            else:
                for i in range(len(x)):
                    if (x[i]>xmin and x[i]<xmax):
                        tmpalign.append(y[i])
                    else:
                        continue

            print('alignment is ', -np.mean(tmpalign))
            platy.append(-np.mean(tmpalign))
            flag = 0
            for i in tmpalign:  #check to see if alignment region varies too much
                if np.abs(i-platy[-1])>0.2:
                    flag = 1
                else:
                    continue
            if flag != 0:
                print('Warning: potential aligned region varied by more ' + \
                      'than 0.2eV (in range of halfway between defects ' + \
                      '+/-1 \Angstrom). Might have issues with Freidel ' + \
                      'oscillations or atomic relaxation')

            if print_pot_flag == 'written' or print_pot_flag == 'plotfull':
                # #commenting this out in favor of just one output written file
                # def write_xy(x, y, fname):
                #     """
                #     Write the x, y vectors to file
                #     """
                #     with open(os.path.join('..',fname),'w') as f:
                #         for pair in zip(x, y):
                #             print >> f, ' '.join(str(i) for i in pair)
                #
                # name = self._name
                # charge = str(self._charge)
                # fname = '_'.join([name,charge,'xylong',str(axis)]) + '.dat'
                # write_xy(x_lr, y_lr, fname)
                # fname = '_'.join([name,charge,'xy',str(axis)]) + '.dat'
                # #fname = name + charge + 'xy' + str(axis) + '.dat'
                # write_xy(x, y, fname)
                # fname = '_'.join([name,charge,'xy',str(axis),'diff.dat'])
                # #fname = name + charge + 'xy' + str(axis) + 'diff.dat'
                # write_xy(x_diff, y_diff, fname)

                ##THIS is for outputting just ONE simpler file at end...
                if self._frac_coords[axis]:
                    axblkval = self._frac_coords[axis]*self._lengths[axis]
                    for i in range(len(x_lr)):
                        if axblkval< x_lr[i]:
                            print 'axblkval= ',axblkval,' i = ',i
                            break
                    #rollind=len(x_lr)-i
                    y_lr=np.roll(y_lr,-i)
                    y=np.roll(y,-i)
                    y_diff=np.roll(y_diff,-i)

                #this is because current uploading format for written is bad for numpy arrays...
                #Might want to update this in future so that dumping format is smarter than just dumping/loading a string
                if print_pot_flag == 'written':
                    y_lrtmp=[]
                    ytmp=[]
                    y_difftmp=[]
                    for j in range(len(x_lr)):
                        y_lrtmp.append(y_lr[j])
                        ytmp.append(y[j])
                        y_difftmp.append(y_diff[j])
                    y_lr = y_lrtmp
                    y = ytmp
                    y_diff = y_difftmp

                plotvals[str(axis)]['xylong'] = [x_lr,y_lr]
                plotvals[str(axis)]['xy'] = [x,y]
                plotvals[str(axis)]['xydiff'] = [x_diff,y_diff]


        if print_pot_flag == 'written': #then output a simple output written file
            fname='_'.join([str(name),'Frey',str(len(self._axiscalcs)),'AxisData'])+'.dat'
            with open(os.path.join('..',fname),'w') as f:
                f.write(str(plotvals))

        if print_pot_flag == 'plotfull':  #plot all three planar averaged potentials
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(15.0,12.0))
            print('plot full plot')
            for axis in self._axiscalcs:
                print('plot axis ',axis+1)
                ax = fig.add_subplot(len(self._axiscalcs), 1, axis+1)
                ax.set_ylabel('axis '+str(axis+1))
                #pylab.hold(True)
                vals_plot = plotvals[str(axis)]
                ax.plot(vals_plot['xy'][0], vals_plot['xy'][1])
                ax.plot(vals_plot['xydiff'][0], vals_plot['xydiff'][1], 'r')
                ax.plot(vals_plot['xylong'][0], vals_plot['xylong'][1], 'g')

                if not axis:
                    ax.set_title("Sxdefectalign electrostatic planar averaged potential")
                    ax.legend(['V_defect-V_ref-V_lr','V_defect-V_ref','V_lr'])

                latt_len = self._lengths[axis]
                fcoords = self._frac_coords[axis]
                # ax.plot([fcoords*latt_len], [0], 'or', markersize=4.0)
                # if fcoords >= 0.5:
                #     ax.plot([(fcoords-0.5)*latt_len],[0],'og',markersize=4.0)
                # else:
                #     ax.plot([(fcoords+0.5)*latt_len],[0],'og',markersize=4.0)
                #plt has been shited to be centered on defect
                ax.plot([fcoords], [0], 'or', markersize=4.0)
                ax.plot([(0.5)*latt_len],[0],'og',markersize=4.0)
                plt.axhline(y=0.0, linewidth=0.8, color='black')
            plt.savefig(os.path.join('..', "locpotgraph.png"))

        return [result,platy]

    def plot_from_datfile(self,name='default_Frey_3_AxisData.dat',widthsample=0.5):
        """
        Takes data file called 'name' and does plotting.
        Good for later plotting of locpot data after running run_correction()
        """
        import ast
        with open(name,'r') as f:
            plotvals=f.read()
        plotvals=ast.literal_eval(plotvals) #converting string to dictionary

        axiscalcs=plotvals.keys()
        axiscalcs.sort()
        axiscalcs=[int(i) for i in axiscalcs] #change to integer type

        #shifting values to origin for better looking plot; havent done this for plotfull option yet...
        import matplotlib.pyplot as plt
        #fig = plt.figure(figsize=(15.0,12.0))
        fig=plt.figure()
        plt.clf()
        print('plot full plot')
        for axis in axiscalcs:
            print('plot axis ',axis+1)
            ax = fig.add_subplot(len(axiscalcs), 1, axis+1)
            ax.set_ylabel('axis '+str(axis+1))
            #pylab.hold(True)
            vals_plot = plotvals[str(axis)]
            ax.plot(vals_plot['xylong'][0], vals_plot['xylong'][1], 'g')
            ax.plot(vals_plot['xydiff'][0], vals_plot['xydiff'][1], 'r')
            ax.plot(vals_plot['xy'][0], vals_plot['xy'][1],'b')

            if not axis:
                ax.set_title("Sxdefectalign electrostatic planar averaged potential")
                ax.legend(['long range from model','DFT locpot diff','short range (aligned)'],loc=9)

            latt_len = self._lengths[axis]
            fcoords = self._frac_coords[axis]
            #have moved defect to origin
            #ax.plot([fcoords*latt_len], [0], 'or', markersize=4.0)
            # if fcoords >= 0.5:
            #     ax.plot([(fcoords-0.5)*latt_len],[0],'og',markersize=4.0)
            # else:
            #     ax.plot([(fcoords+0.5)*latt_len],[0],'og',markersize=4.0)
            x=vals_plot['xy'][0]
            mid=len(x)/2
            checkdis=int((widthsample/2.)/(x[1]-x[0])) #for checking a window of 1 around defect
            tmpx=[x[i] for i in range(mid - checkdis, mid + checkdis)]
            ax.fill_between(tmpx, -100, 100, facecolor='red', alpha=0.15, label='sampling region')
            ax.set_xlim(round(x[0]),round(x[-1]))
            ymin=min(min(vals_plot['xy'][1]),min(vals_plot['xydiff'][1]),min(vals_plot['xylong'][1]))
            ymax=max(max(vals_plot['xy'][1]),max(vals_plot['xydiff'][1]),max(vals_plot['xylong'][1]))
            ax.set_ylim(-0.2+ymin,ymax+0.2)
            ax.axhline(y=0.0, linewidth=0.8, color='black')
        plt.savefig(os.path.join("Freylocpotgraph.png")) #do something smarter than this?

        return

    def run_correction(self, print_pot_flag='written', partflag='All'):
        """
        Runs all neccessary parts to get freysoldt corrections out with
        planar averaged potentials
        Change plot_pot_flag if you want plotted planar averages
        set transflag to True if you want to write flags

        part flag can be 'pc' for just returning point charge correction,
               'potalign' for just potalign correction, 'All' for one combined correction,
              or 'AllSplit' for correction in form [PC,potterm,full]
        """
        with ScratchDir('.'):
            self.prepare_files()
            s = self.plot_pot_diff(print_pot_flag='none')
            print('--')
            print('potential alignments determined to be: '+str(s[1]))
            print('get final correction terms')
            print('--')
            #NOTE this should not be neccessary...should be able to know answer after one sxdefectalign call...
            #To get locpot plots use print_pot_flag = 'written' or 'plotfull'
            vals = self.plot_pot_diff(align=s[1], print_pot_flag=print_pot_flag)
            print('vals is '+str(vals))
            for i in self._axiscalcs:
                if np.abs(vals[1][i]) > 0.0001:
                    print('PROBLEM! planar averaging didnt work. Issue ' +\
                            'with axis', str(i+1))
            print('--')
            print('Final Freysoldt sxdefectalign correction values are: ', \
                    str(vals[0]))

        if partflag=='pc':
            return s[0]
        elif partflag=='potalign':
            return [s[1][i]*self._charge for i in range(3)]
        elif partflag=='All':
            return vals[0]
        else:
            return [s[0],[round(s[1][i]*self._charge,5) for i in range(len(s[1]))],vals[0]]


if __name__ == '__main__':
    s=SxDefectCorrection('../../bulk/LOCPOT', 'LOCPOT', 1, 11.814,
                 [ 0.375,  0.125,  0.125], 520, lengths=[13.258098, 13.258098, 13.258098])
    #sxvals=s.run_correction(partflag='AllSplit')
    s.plot_from_datfile(name='default_Frey_1_AxisData.dat')
    #print sxvals
    #s.plot_from_datfile(name='default_Frey_1_AxisData.dat')
