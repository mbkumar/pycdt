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

from monty.tempfile import ScratchDir

class DefectCorrectionFreysoldt(object):

    """
    This class applies the Freysoldt correction to remove electrostatic defect
    interaction contribution to energy and apply potential alignment.
    Ideally this sxdefectalign wrapper class should be replaced by a python
    code implementing a generalized anisotropic Freysoldt correction
    """
    
    def __init__(self, locpot_bulk, locpot_defect, charge, epsilon, 
            frac_coords,encut):
        """
        locpot_bulk = LOCPOT of bulk as Locpot object
        locpot_defect=LOCPOT of defect as Locpot object
        """
        self._locpot_bulk=locpot_bulk
        self._locpot_defect=locpot_defect
        self._charge=charge #relative to neutral defect (not relative to bulk)
        self._epsilon=epsilon #this should include ionic relaxation effects
        self._frac_coords=frac_coords   #in form [0,0,0]
        self._encut=encut #encut value in eV from Vasp
        
    def prepare_files(self):
        if  self._charge==0:
            print 'defect has charge 0, so freysoldt correction is 0'
            return
        print 'Need to prep Locpots'
        self._locpot_bulk.write_file("LOCPOT_vref",True)
        self._locpot_defect.write_file("LOCPOT_vdef",True)
        print 'locpots prepared for sxdefectalign'

    def plot_hartree_pot(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials')
        get_agrid = self._locpot_bulk.get_axis_grid
        get_aavg = self._locpot_bulk.get_average_along_axis
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            latt_len = self._locpot_defect.structure.lattice._lengths[axis]
            ax.plot(get_agrid(axis),get_aavg(axis),'r',
                    label="Bulk potential")
            ax.plot(get_agrid(axis), get_aavg(axis),'b',
                    label="Defect potential")
            ax.plot([self._frac_coords[axis]*latt_len], [0], 'or', 
                    markersize=4.0, label="Defect site")
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        #plt.savefig('locpotavgplot.png')
        plt.show()

    def plot_hartree_pot_diff(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potential difference')
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            defect_axis=self._locpot_defect.get_axis_grid(axis)
            defect_pot=self._locpot_defect.get_average_along_axis(axis)
            pure_pot=self._locpot_bulk.get_average_along_axis(axis)
            latt_len = self._locpot_defect.structure.lattice._lengths[axis]
            ax.plot(defect_axis,defect_pot-pure_pot,'b',
                    label='Defect-Bulk difference')
            ax.plot([self._frac_coords[axis] * latt_len], [0], 'or',
                    markersize=4.0, label='Defect site')
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        #plt.savefig('locpotdiffplot.png')
        plt.show()

    def plot_all_hartree_pot(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials and difference')
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            defect_axis=self._locpot_defect.get_axis_grid(axis)
            defect_pot=self._locpot_defect.get_average_along_axis(axis)
            pure_pot=self._locpot_bulk.get_average_along_axis(axis)
            ax.plot(defect_axis,pure_pot,'r',label="Bulk potential")
            ax.plot(defect_axis,defect_pot,'b',label="Defect potential")
            ax.plot(defect_axis,defect_pot-pure_pot,'k',
                    label='Defect-Bulk difference')
            ax.plot([self._frac_coords[axis] * \
                    self._locpot_defect.structure.lattice._lengths[axis]],\
                    [0],'or',markersize=4.0,label='Defect site')
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        #plt.savefig('locpotavgdiffplot.png')
        plt.show()

    def plot_pot_diff(self, align=[0.0,0.0,0.0], print_pot_flag='written'):
        """
        runs sxdefectalign and obtains three alignment constants for 
        each axis. To do freysoldt correction you need to run this twice 
        (once to get alignment, and second time to get final correction 
        and plots). To do this just call the run_correction() method with
        print flag = none then just get corrections and leave vline-eV.dat 
        files for each axis and produce potential alignment correction
        or if  print flag= written then get different parts of vline-eV.dat 
        written to files for plotting (good for nersc or where you can't 
        plot easy) or if print flag=plotfull then show the full matplotlib 
        plot for all three axes at end of calculation encut is eV cutoff 
        from VASP

        Would like final workflow to include the flag that I have put here 
        for when planar average varies by more than 0.2 eV around far region
        """
        if self._charge==0: #don't need charge correction if charge is zero
            return [[0,0,0],[0,0,0]]
        #correction from output (should include alignment once alignment 
        # has been done)
        result = []   
        platy = []    #alignment terms for each axis
        # if want to plot right here, then build dictionary for storing 
        # planar average values of each axis
        if print_pot_flag == 'plotfull':  
            plotvals = {'0':{},'1':{},'2':{}}
        for axis in [0,1,2]:
            print 'do axis '+str(axis+1)
            relpos = (str(self._frac_coords)[1:])[:-1]
            relpos = relpos.replace(" ","")
            command = ['~/sxdefectalign', '--vasp', '-a'+str(axis+1), 
                    '--relative', '--pos', relpos, 
                    '--charge', str(-self._charge), 
                    '--ecut', str(self._encut/13.6057), #eV to Ry for sxdefect 
                    '--eps', str(self._epsilon), 
                    '-C', str(-float(align[axis])), 
                    '--vref', 'LOCPOT_vref', 
                    '--vdef', 'LOCPOT_vdef']
            print command

            #standard way of running NERSC commands.
            #in case NERSC (Hopper) has issues with python subprocess can use hack
            p = subprocess.Popen(command, stdout=subprocess.PIPE, 
                    stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            output, err = p.communicate()
            out = out.decode('utf-8')
            err = err.decode('utf-8')
            print 'output from sxdefectalign = ', str(out)
            result.append(float(out[0].split("\n")[12].split()[3]))
            print "chg correction is "+str(result[-1])

            ##this is hack wrap-around for when subprocess doesn't work
            ##(which is always an issue on hopper now...)
            #cmd=''
            #for i in range(len(command)):
            #    cmd+=command[i]+' '
            #cmd+=' > tmpoutput'
            #os.system(cmd)
            #output=[]
            #f=open('tmpoutput','r')
            #for r in f.readlines():
            #    output.append(r)
            #f.close()
            #print 'output from sxdefectalign = '+str(output)
            #result.append(float(output[-1].split()[3]))
            #print "chg correction is "+str(result[-1])
            #os.remove('tmpoutput')

            x_lr, y_lr = [], []
            x, y = [], []
            x_diff, y_diff = [], []
            with open("vline-eV.dat",'r') as f_sr: #read in potential 
                for r in f_sr.readlines(): 
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
            latt_len = self._locpot_defect.structure.lattice._lengths[axis]
            if self._frac_coords[axis] >= 0.5:
                 platx = (self._frac_coords[axis]-0.5) * latt_len
            else:
                 platx = (self._frac_coords[axis]+0.5) * latt_len
            print "half way between defects is: ", platx

            xmin = latt_len - (1-platx) if platx < 1 else platx-1
            xmax = 1-(latt_len-platx) if platx > latt_len-1 else 1+platx
            print 'means sampling region is (', xmin, ',', xmax, ')'

            tmpalign=[]
            if xmax < xmin:
                print 'wrap around detected, special alignment needed'
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

            print 'alignment is ', -np.mean(tmpalign)
            platy.append(-np.mean(tmpalign))
            flag = 0
            for i in tmpalign:  #check to see if alignment region varies too much
                if np.abs(i-platy[-1])>0.2:
                    flag = 1
                else:
                    continue
            if flag != 0:
                print 'Warning: potential aligned region varied by more ' + \
                      'than 0.2eV (in range of halfway between defects ' + \
                      '+/-1 \Angstrom). Might have issues with Freidel ' + \
                      'oscillations or atomic relaxation'

            if print_pot_flag == 'written':
                with open(os.path.join('..',"xylong"+str(axis)+".dat"),'w') as f:
                    for i in range(len(x_lr)):
                        f.write(str(x_lr[i])+" "+str(y_lr[i])+"\n")
                with open(os.path.join('..',"xy"+str(axis)+".dat"),'w') as f:
                    for i in range(len(x)):
                        f.write(str(x[i])+" "+str(y[i])+"\n")
                with open(os.path.join('..',"xy"+str(axis)+"diff.dat"),'w') as f:
                    for i in range(len(x_diff)):
                        f.write(str(x_diff[i])+" "+str(y_diff[i])+"\n")

            elif print_pot_flag == 'plotfull': #store data for plotting at end of all calcs
                plotvals[str(axis)]['xylong'] = [x_lr,y_lr]
                plotvals[str(axis)]['xy'] = [x,y]
                plotvals[str(axis)]['xydiff'] = [x_diff,y_diff]

        if print_pot_flag == 'plotfull':  #plot all three planar averaged potentials
            import matplotlib.pyplot as plt
            import pylab
            fig=plt.figure(figsize=(15.0,12.0))
            print 'plot full plot'
            for axis in [0,1,2]:
                print 'plot axis ',axis+1
                ax=fig.add_subplot(3,1,axis+1)
                ax.set_ylabel('axis '+str(axis+1))
                #pylab.hold(True)
                vals_plot = plotvals[str(axis)]
                ax.plot(vals_plot['xy'][0],vals_plot['xy'][1])
                ax.plot(vals_plot['xydiff'][0],vals_plot['xydiff'][1],'r')
                ax.plot(vals_plot['xylong'][0],vals_plot['xylong'][1],'g')

                if axis==0:
                    ax.set_title("Electrostatic planar averaged potential")
                    ax.legend(['V_defect-V_ref-V_lr','V_defect-V_ref','V_lr'])

                latt_len = self._locpot_defect.structure.lattice._lengths[axis]
                fcoords = self._frac_coords[axis]
                ax.plot([fcoords*latt_len], [0], 'or', markersize=4.0)
                if fcoords >= 0.5:
                    ax.plot([(fcoords-0.5)*latt_len],[0],'og',markersize=4.0)
                else:
                    ax.plot([(fcoords+0.5)*latt_len],[0],'og',markersize=4.0)
                plt.axhline(y=0.0,linewidth=0.8, color='black')
            plt.savefig(os.path.join('..',"locpotgraph.png"))

        return [result,platy]

    def run_correction(self):
        """
        Runs all neccessary parts to get freysoldt corrections out with
        planar averaged potentials
        Change plot_pot_flag if you want plotted planar averages
        """
        with ScratchDir('.'):
            self.prepare_files()
            s = self.plot_pot_diff(print_pot_flag='none')
            print '--'
            print 'potential alignments determined to be: '+str(s[1])
            print 'get final correction terms'
            print '--'
            #To get locpot plots use print_pot_flag = 'written' or 'plotfull'
            vals = self.plot_pot_diff(align=s[1], print_pot_flag='none')   
            print 'vals is '+str(vals)
            for i in range(3):
                if np.abs(vals[1][i]) > 0.0001:
                    print 'PROBLEM! planar averaging didnt work. Issue ' +\
                            'with axis', str(i+1)
            print '--'
            print 'Final Freysoldt sxdefectalign correction values are: ', \
                    str(vals[0])

        return vals[0]


if __name__ == '__main__':

    #for testing, example given here is Se_sn+1 in 72 atom cell
    #from pymatgen.io.vaspio.vasp_output import Locpot
    #print 'load locpots'
    #sdef=Locpot.from_file("LOCPOT")
    #spure=Locpot.from_file("../pure/LOCPOT")
    #s1=DefectCorrectionFreysoldt(spure,sdef,1,47.852,[0.0833330000000032,0.0307343554185392,0.3830636916206969],520)
    #s1.run_correction()

    #test=s1.plot_pot_diff(align=[0.0,0.0,0.0],print_pot_flag='none')
    #print 'output of plotting code is:'+str(test)
