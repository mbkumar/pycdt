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
__date__ = "April 10, 2015"

import subprocess
import os
import numpy as np


class DefectCorrectionFreysoldt(object):

    """
    This class apply the Freysoldt correction to remove electrostatic defect interaction and apply
    electostatic alignemnt
    Ideally this wrapper around sxdefectalign should be replaced by a python code implementing
    a generalized anisotropic Freysoldt correction

    locpot_bulk = LOCPOT of bulk as Locpot object
    locpot_defect=LOCPOT of defect as Locpot object
    """
    
    def __init__(self, locpot_bulk, locpot_defect, charge, epsilon, frac_coords,encut):
        self._locpot_bulk=locpot_bulk
        self._locpot_defect=locpot_defect
        self._charge=charge
        self._epsilon=epsilon
        self._frac_coords=frac_coords   #in form [0,0,0]
        self._encut=encut #encut value in eV from Vasp
        
        
    def prepare_files(self):
        if self._charge==0:
            print 'defect has charge 0! dont need freysoldt'
            return
        if not os.path.exists("tmp"):
            os.makedirs("tmp")
        if (os.path.exists("tmp/LOCPOT_vdef")==True and os.path.exists("tmp/LOCPOT_vref")==True):
            print 'locpot already written, good!'
        else:
            print 'Need to prep Locpots'
            self._locpot_bulk.write_file("tmp/LOCPOT_vref",True)
            self._locpot_defect.write_file("tmp/LOCPOT_vdef",True)
            print 'locpots prepared for sxdefectalign'

    def plot_hartree_pot(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials')
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            ax.plot(self._locpot_bulk.get_axis_grid(axis),self._locpot_bulk.get_average_along_axis(axis),'r',label="Bulk potential")
            ax.plot(self._locpot_defect.get_axis_grid(axis),self._locpot_defect.get_average_along_axis(axis),'b',label="Defect potential")
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0,label="Defect site")
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        pylab.show()

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
            ax.plot(defect_axis,defect_pot-pure_pot,'b',label='Defect-Bulk difference')
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0,label='Defect site')
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        pylab.show()

    def plot_all_hartree_pot(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials')
        for axis in [0,1,2]:
            ax = fig.add_subplot(3, 1, axis+1)
            defect_axis=self._locpot_defect.get_axis_grid(axis)
            defect_pot=self._locpot_defect.get_average_along_axis(axis)
            pure_pot=self._locpot_bulk.get_average_along_axis(axis)
            ax.plot(defect_axis,pure_pot,'r',label="Bulk potential")
            ax.plot(defect_axis,defect_pot,'b',label="Defect potential")
            ax.plot(defect_axis,defect_pot-pure_pot,'k',label='Defect-Bulk difference')
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0,label='Defect site')
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        pylab.show()

    def GEOFFplot_hartree_pot(self):
        """
        I really don't understand why we can't use the more simplified version I have done above...seems like the multiple stuff is all unnecessary
        """
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potentials')
        for axis in [0,1,2]:

            multiple = int(round(max(self._locpot_defect.get_axis_grid(axis))/max(self._locpot_bulk.get_axis_grid(axis))))
            supercell_axis = [self._locpot_bulk.get_axis_grid(axis)[j] 
                              + i * self._locpot_bulk.structure.lattice.abc[axis]
                            for i in range(0, multiple)
                            for j in range(0, len(self._locpot_bulk.get_axis_grid(axis)))]
            tmp_av=self._locpot_bulk.get_average_along_axis(axis)
            supercell_pot = [tmp_av[j]
                            for i in range(0, multiple)
                            for j in range(0, len(tmp_av))]
            print axis, min(self._locpot_bulk.get_axis_grid(axis)), max(self._locpot_bulk.get_axis_grid(axis))
            print axis, min(supercell_axis), max(supercell_axis)
            print axis, min(self._locpot_defect.get_axis_grid(axis)), max(self._locpot_defect.get_axis_grid(axis))

            ax = fig.add_subplot(3, 1, axis+1)
            ax.plot(supercell_axis,
                       supercell_pot,'r',label="Bulk potential")
            ax.plot(self._locpot_defect.get_axis_grid(axis),self._locpot_defect.get_average_along_axis(axis),'b',label="Defect potential")
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0,label="Defect site")
            ax.set_ylabel("axis "+str(axis+1))
            if axis==0:
                ax.legend()
        ax.set_xlabel("distance (Angstrom)")
        pylab.show()
    
    def GEOFFplot_hartree_pot_diff(self):
        """
        Again...multiple stuff seems unnecessarily complicated
        """
        import matplotlib.pyplot as plt
        fig=plt.figure()
        ax=fig.add_subplot(3,1,1)
        ax.set_title('Locpot planar averaged potential difference')
        for axis in [0,1,2]:
            multiple = int(round(max(self._locpot_defect.get_axis_grid(axis))/max(self._locpot_bulk.get_axis_grid(axis))))
            supercell_axis = [self._locpot_bulk.get_axis_grid(axis)[j] 
                              + i * self._locpot_bulk.structure.lattice.abc[axis]
                            for i in range(0, multiple)
                            for j in range(0, len(self._locpot_bulk.get_axis_grid(axis)))]
            tmp_av=self._locpot_bulk.get_average_along_axis(axis)
            supercell_pot = [tmp_av[j]
                            for i in range(0, multiple)
                            for j in range(0, len(tmp_av))]
            ax = fig.add_subplot(3, 1, axis+1)
            defect_axis=self._locpot_defect.get_axis_grid(axis)
            defect_pot=self._locpot_defect.get_average_along_axis(axis)
            ax.plot(supercell_axis,[defect_pot[i]-supercell_pot[i] for i in range(0,len(defect_pot))],'b')
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0)
            ax.set_ylabel("axis "+str(axis+1))
        ax.set_xlabel("distance (Angstrom)")
        pylab.show()
    
    def NEEDSFIXINGtest_rho(self,axis=2,alignment_sr=0.0):
        """
        This is supposed to be for testing if different values of rho produce different values for the freysoldt correction
        COULDBE CLEANED UP AT SOME POINT...Are we looking for differences to plot or numerical differences in energy corrections?
        """
        for expnorm in [0.0,0.25,0.5,0.75,1.0]:
            command=['~/sxdefectalign']
            #./sxdefectalign -C 1.6 --vasp --charge -2 --ecut 6 --vref LOCPOT_SnO --vdef LOCPOT_SnO_Vac_O_2
            command.append("--vasp")
            command.append("-a"+str(axis+1))
            command.append("--relative")
            command.append("--pos")
            relpos=(str(self._frac_coords)[1:])[:-1]
            relpos=relpos.replace(" ","")
            command.append(relpos)
            command.append("--expnorm")
            command.append(str(expnorm))
            command.append("--charge")
            command.append(str(self._charge))
            command.append("--ecut")
            command.append("6")
            command.append("--eps")
            command.append(str(self._epsilon))
            command.append("-C")
            command.append(str(alignment_sr))
            command.append("--vref")
            command.append("tmp/LOCPOT_vref")
            command.append("--vdef")
            command.append("tmp/LOCPOT_vdef")
            print command
            p = subprocess.call(command, stdout=subprocess.PIPE, close_fds=True)
            output=p.communicate()
            print output
            f_sr=open("vline-eV.dat",'r')
            x=[]
            y=[]
            for r in f_sr.readlines():
                tmp=r.split("\t")
                if len(tmp)>2:
                    x.append(float(tmp[0]))
                    y.append(float(tmp[2].rstrip("\n")))
            
            print x
            print y
            pylab.plot(x,y)
            f_sr.close()
            os.remove("vline-eV.dat")
        pylab.figure()
        pylab.show()
    
    def NEEDSFIXINGcompute_correction_rho(self,axis=2,alignment_sr=0.0,plot=False,expnorm=0.25):
        """
        Needs to be fixed like the 'plot_pot_diff' attribute, but not sure what is added from this? information about expnorm and rho...
        """
        command=['~/sxdefectalign']
        #./sxdefectalign -C 1.6 --vasp --charge -2 --ecut 6 --vref LOCPOT_SnO --vdef LOCPOT_SnO_Vac_O_2
        command.append("--vasp")
        command.append("-a"+str(axis+1))
        command.append("--relative")
        command.append("--pos")
        command.append(str(self._frac_coords))
        #command.append("[0,0,0]")
        command.append("--expnorm")
        command.append("0.5")
        command.append("--charge")
        command.append(str(self._charge))
        command.append("--ecut")
        command.append("20")
        command.append("--eps")
        command.append(str(self._epsilon))
        command.append("-C")
        command.append(str(alignment_sr))
        command.append("--vref")
        command.append("/home/geoffroy/tmp/LOCPOT_vref")
        command.append("--vdef")
        command.append("/home/geoffroy/tmp/LOCPOT_vdef")
        print command
        p = subprocess.Popen(command, stdout=subprocess.PIPE, close_fds=True)
        output=p.communicate()
        f_sr=open("vline-eV.dat",'r')
        x_lr=[]
        y_lr=[]
        x=[]
        y=[]
        if plot==True:
            for r in f_sr.readlines():
                tmp=r.split("\t")
                if(len(tmp)<3 and not r.startswith("&")):
                   x_lr.append(float(tmp[0]))
                   y_lr.append(float(tmp[1])) 
                if len(tmp)>2:
                    x.append(float(tmp[0]))
                    y.append(float(tmp[2].rstrip("\n")))
            pylab.figure()
            pylab.plot(x,y)
            #print [min(x),max(y)]
            #print [y[int(math.floor(len(y)/2))],y[int(math.floor(len(y)/2))]]
            #pylab.plot([min(x),max(y)],[y[int(math.floor(len(y)/2))],y[int(math.floor(len(y)/2))]],'k-')
            #print self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]*1.889725989
            pylab.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]*1.889725989],[0],'or',markersize=4.0)
            #pylab.plot(x_lr,y_lr)
            pylab.show()
            pylab.show()
        print output[0]
        print output[0].split("\n")[12].split()[3]
        return float(output[0].split("\n")[12].split()[3])

    def plot_pot_diff(self, align=[0.0,0.0,0.0],print_pot_flag='written'):
        """
        runs sxdefectalign and obtains three alignment constants for each axis. To do freysoldt correction you need to run this twice (once to get alignment, and second time to get final correction and plots)
        To do this just use the run_correction() attribute below
        print flag = none then just get corrections and leave vline-eV.dat files for each axis and produce potential alignment correction
            or if  print flag= written then get different parts of vline-eV.dat written to files for plotting (good for nersc or where you can't plot easy)
            or if print flag=plotfull then show the full matplotlib plot for all three axes at end of calculation
        encut is eV cutoff from VASP

        The planar alignment is obtained by looking at the average in  window of +/- 1 Angstrom around the furthest point from a defet (which is neccessary in ionic relaxation planar averages)
        Would like final workflow to include the flag that I have put here for when planar average varies by more than 0.2 eV around far region as possible flag for defect delocalization
        """
        if self._charge==0:
            return [[0,0,0],[0,0,0]]
        result=[]   #correction from output (will include alignment correction once alignment has been done)
        platy=[]    #alignment terms
        if print_pot_flag=='plotfull':  #if want to plot right here, then build dictionary for storing planar average values of each axis
            plotvals={}
            plotvals['0']={}
            plotvals['1']={}
            plotvals['2']={}
        for axis in [0,1,2]:
            print 'do axis '+str(axis+1)
            command=['~/sxdefectalign']
            command.append("--vasp")
            command.append("-a"+str(axis+1))
            command.append("--relative")
            command.append("--pos")
            relpos=(str(self._frac_coords)[1:])[:-1]
            relpos=relpos.replace(" ","")
            command.append(relpos)
            command.append("--charge")
            command.append(str(-self._charge))
            command.append("--ecut")
            command.append(str(self._encut/13.6057))  #converts eV to Ry for sxdefectalign
            command.append("--eps")
            command.append(str(self._epsilon))
            command.append("-C")
            command.append(str(-float(align[axis])))
            command.append("--vref")
            command.append("tmp/LOCPOT_vref")
            command.append("--vdef")
            command.append("tmp/LOCPOT_vdef")
            print command

            #for some reason subprocess command is having issues on NERSC/hopper due to python upgrades
            #p = subprocess.Popen(command, stdout=subprocess.PIPE, close_fds=True)
            #output=p.communicate()
            #print 'output from sxdefectalign = '+str(output)
            #result.append(float(output[0].split("\n")[12].split()[3]))
            #print "chg correction is "+str(result[-1])

            #this is rediculous wraparound to deal with subprocess not working...
            print 'running wrap around to subprocess command due to NERSC python issues'
            cmd=''
            for i in range(len(command)):
                cmd+=command[i]+' '
            cmd+=' > tmpoutput'
            os.system(cmd)
            output=[]
            f=open('tmpoutput','r')
            for r in f.readlines():
                output.append(r)
            f.close()
            #print 'output from sxdefectalign = '+str(output)
            result.append(float(output[-1].split()[3]))
            print "chg correction is "+str(result[-1])
            os.remove('tmpoutput')


            f_sr=open("vline-eV.dat",'r')
            x_lr=[]
            y_lr=[]
            x=[]
            y=[]
            x_diff=[]
            y_diff=[]
            for r in f_sr.readlines():  #read in different parts of potential from vline-eV.dat file
                #print r
                tmp=r.split("\t")
                if(len(tmp)<3 and not r.startswith("&")):
                   x_lr.append(float(tmp[0])/1.889725989)   #convert to Angstrom
                   y_lr.append(float(tmp[1])) 
                if len(tmp)>2:
                    x.append(float(tmp[0])/1.889725989)   #convert to Angstrom
                    x_diff.append(float(tmp[0])/1.889725989)   #convert to Angstrom
                    y.append(float(tmp[2].rstrip("\n")))
                    y_diff.append(float(tmp[1]))
            f_sr.close()

            if print_pot_flag=='none':
                if os.path.exists("vline-eV.dat"):
                    os.rename("vline-eV.dat","axis"+str(axis)+"vline-eV.dat")
            else:
                if os.path.exists("vline-eV.dat"):
                    os.remove("vline-eV.dat")

            #now extract potential alignment term averaging window of +/- 1 angstrom around point halfway between neighboring defects
            if self._frac_coords[axis]>=0.5:
                 platx=float(self._frac_coords[axis]-0.5)*float(self._locpot_defect.structure.lattice._lengths[axis])
            else:
                 platx=float(self._frac_coords[axis]+0.5)*float(self._locpot_defect.structure.lattice._lengths[axis])
            print "half way between defects is: "+str(platx)
            if platx<1:
                xmin=float(self._locpot_defect.structure.lattice._lengths[axis])-(1-platx)
            else:
                xmin=platx-1

            if platx>float(self._locpot_defect.structure.lattice._lengths[axis])-1:
                xmax=1-(float(self._locpot_defect.structure.lattice._lengths[axis])-platx)
            else:
                xmax=platx+1
            print 'means sampling region is ('+str(xmin)+','+str(xmax)+')'

            tmpalign=[]
            if xmax<xmin:
                print 'wrap around detected, special alignment needed'
                for i in range(len(x)):
                    if x[i]<xmax:
                        tmpalign.append(y[i])
                    elif x[i]>xmin:
                        tmpalign.append(y[i])
                    else:
                        continue
            else:
                for i in range(len(x)):
                    if (x[i]>xmin and x[i]<xmax):
                        tmpalign.append(y[i])
                    else:
                        continue

            print 'alignment is '+str(-np.mean(tmpalign))
            platy.append(-np.mean(tmpalign))
            flag=0
            for i in tmpalign:
                if np.abs(i-platy[-1])>0.2:
                    flag=1
                else:
                    continue
            if flag!=0:
                    print 'Warning: potential aligned region varied by more than 0.2eV (in range of halfway between defects +/-1 Angstrom). Might have issues with Freidel oscilattions or atomic relaxation'

            if print_pot_flag=='written':
                f=open("xylong"+str(axis)+".dat",'w')
                for i in range(len(x_lr)):
                        f.write(str(x_lr[i])+" "+str(y_lr[i])+"\n")
                f.close()
                f=open("xy"+str(axis)+".dat",'w')
                for i in range(len(x)):
                        f.write(str(x[i])+" "+str(y[i])+"\n")
                f.close()
                f=open("xy"+str(axis)+"diff.dat",'w')
                for i in range(len(x_diff)):
                        f.write(str(x_diff[i])+" "+str(y_diff[i])+"\n")
                f.close()
            elif print_pot_flag=='plotfull': #store data for end of all three calcs
                plotvals[str(axis)]['xylong']=[x_lr,y_lr]
                plotvals[str(axis)]['xy']=[x,y]
                plotvals[str(axis)]['xydiff']=[x_diff,y_diff]

        if print_pot_flag=='plotfull':  #plot all three
            import matplotlib.pyplot as plt
            import pylab
            fig=plt.figure(figsize=(15.0,12.0))
            print 'plot full plot'
            for axis in [0,1,2]:
                print axis+1
                ax=fig.add_subplot(3,1,axis+1)
                plt.ylabel('axis '+str(axis+1))
                pylab.hold(True)
                ax.plot(plotvals[str(axis)]['xy'][0],plotvals[str(axis)]['xy'][1])
                ax.plot(plotvals[str(axis)]['xydiff'][0],plotvals[str(axis)]['xydiff'][1],'r')
                ax.plot(plotvals[str(axis)]['xylong'][0],plotvals[str(axis)]['xylong'][1],'g')
                if axis==0:
                    ax.set_title("Electrostatic planar averaged potential")
                    ax.legend(['V_defect-V_ref-V_lr','V_defect-V_ref','V_lr'])
                ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0)
                if self._frac_coords[axis]>=0.5:
                        ax.plot([(self._frac_coords[axis]-0.5)*self._locpot_defect.structure.lattice._lengths[axis]],[0],'og',markersize=4.0)
                else:
                        ax.plot([(self._frac_coords[axis]+0.5)*self._locpot_defect.structure.lattice._lengths[axis]],[0],'og',markersize=4.0)
                plt.axhline(y=0.0,linewidth=0.8, color='black')
            plt.savefig("locpotgraph.png")

        return [result,platy]

    def run_correction(self):
        """
        This runs all neccessary parts to get freysoldt correction out with planar averaged potential
        """
        self.prepare_files()
        s=self.plot_pot_diff(align=[0.0,0.0,0.0],print_pot_flag='none') #first run sxdefect align on each axis to get alignment corrections for next step
        print '--'
        print 'alignments determined to be: '+str(s[1])
        print 'get final correction terms'
        print '--'
        vals=self.plot_pot_diff(align=s[1],print_pot_flag='none')   #if you want to get out locpot plots then change the print_pot_flag to 'written' or 'plotfull'
        print 'vals is '+str(vals)
        for i in range(3):
            if np.abs(vals[1][i])>0.0001:
                print 'problem! planar averaging didnt work? issue with axis '+str(i+1)
        print '---'
        print 'final Freysoldt sxdefectalign correction values are: '+str(vals[0])

        return vals[0]

    def GEOFFplot_pot_diff(self, align=0.0):
        """
        I don't think I need to keep this
        """
        import matplotlib.pyplot as plt

        fig=plt.figure()
        #./sxdefectalign -C 1.6 --vasp --charge -2 --ecut 6 --vref LOCPOT_SnO --vdef LOCPOT_SnO_Vac_O_2
        for axis in [0,1,2]:
            print axis
            command=['/home/geoffroy/research/TCO/defects/sxdefectalign']
            command.append("--vasp")
            command.append("-a"+str(axis+1))
            command.append("--relative")
            command.append("--pos")
            command.append(str(self._frac_coords))
            #command.append("[0.0,0.026,0.14]")
            command.append("--expnorm")
            command.append("0.0")
            command.append("--gamma")
            command.append("5.0")
            command.append("--average")
            command.append("4.0")
            command.append("--charge")
            command.append(str(-self._charge))
            command.append("--ecut")
            command.append("20")
            command.append("--eps")
            command.append(str(self._epsilon))
            command.append("-C")
            command.append(str(align))
            command.append("--vref")
            command.append("/home/geoffroy/tmp/LOCPOT_vref")
            command.append("--vdef")
            command.append("/home/geoffroy/tmp/LOCPOT_vdef")
            print command
            p = subprocess.Popen(command, stdout=subprocess.PIPE, close_fds=True)
            output=p.communicate()
            f_sr=open("vline-eV.dat",'r')
            x_lr=[]
            y_lr=[]
            x=[]
            y=[]
            x_diff=[]
            y_diff=[]

            for r in f_sr.readlines():
                #print r
                tmp=r.split("\t")
                if(len(tmp)<3 and not r.startswith("&")):
                   x_lr.append(float(tmp[0])/1.889725989)
                   y_lr.append(float(tmp[1]))
                if len(tmp)>2:
                    x.append(float(tmp[0])/1.889725989)
                    x_diff.append(float(tmp[0])/1.889725989)
                    y.append(float(tmp[2].rstrip("\n")))
                    y_diff.append(float(tmp[1]))
            f_sr.close()
            if os.path.exists("vline-eV.dat"):
                os.remove("vline-eV.dat")
            #pylab.figure()
            #pylab.plot(x,y)
            ax = fig.add_subplot(3, 1, axis + 1)
            #align in blue!
            pylab.hold(True)
            ax.plot(x,y)
            ax.plot(x_diff,y_diff,'r')
            ax.plot(x_lr,y_lr,'g')
            ax.legend(['V_defect-V_ref-V_lr','V_defect-V_ref','V_lr'])
            #print [min(x),max(y)]
            #print [y[int(math.floor(len(y)/2))],y[int(math.floor(len(y)/2))]]
            #pylab.plot([min(x),max(y)],[y[int(math.floor(len(y)/2))],y[int(math.floor(len(y)/2))]],'k-')
            #print self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]*1.889725989
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0)
            if self._frac_coords[axis]>=0.5:
                ax.plot([(self._frac_coords[axis]-0.5)*self._locpot_defect.structure.lattice._lengths[axis]],[0],'og',markersize=4.0)
            else:
                ax.plot([(self._frac_coords[axis]+0.5)*self._locpot_defect.structure.lattice._lengths[axis]],[0],'og',markersize=4.0)
            #pylab.plot(x_lr,y_lr)
        return plt

    def GEOFFcompute_correction(self,axis=2,alignment_sr=0.0,plot=False, beta=None):
        """
        I don't think I need to keep this
        """
        command=['/home/geoffroy/research/TCO/defects/sxdefectalign']
        #./sxdefectalign -C 1.6 --vasp --charge -2 --ecut 6 --vref LOCPOT_SnO --vdef LOCPOT_SnO_Vac_O_2
        command.append("--vasp")
        command.append("-a"+str(axis+1))
        command.append("--relative")
        command.append("--pos")
        command.append(str(self._frac_coords))
        #command.append("[0,0,0]")
        command.append("--expnorm")
        command.append("0.0")
        command.append("--gamma")
        command.append("10.0")
        command.append("--average")
        command.append("4.0")
        command.append("--charge")
        command.append(str(-self._charge))
        command.append("--ecut")
        command.append("20")
        command.append("--eps")
        command.append(str(self._epsilon))
        command.append("-C")
        command.append(str(alignment_sr))
        command.append("--vref")
        command.append("/home/geoffroy/tmp/LOCPOT_vref")
        command.append("--vdef")
        command.append("/home/geoffroy/tmp/LOCPOT_vdef")
        print command
        p = subprocess.Popen(command, stdout=subprocess.PIPE, close_fds=True)
        output=p.communicate()
        f_sr=open("vline-eV.dat",'r')
        x_lr=[]
        y_lr=[]
        x=[]
        y=[]
        x_diff=[]
        y_diff=[]
        #if plot==True:
        for r in f_sr.readlines():
            #print r
            tmp=r.split("\t")
            if(len(tmp)<3 and not r.startswith("&")):
               x_lr.append(float(tmp[0]))
               y_lr.append(float(tmp[1]))
            if len(tmp)>2:
                x.append(float(tmp[0]))
                x_diff.append(tmp[0])
                y.append(float(tmp[2].rstrip("\n")))
                y_diff.append(float(tmp[1]))
        return (float(output[0].split("\n")[13].split()[3]),[x_lr,y_lr],[x_diff,y_diff],[x,y])



##below is example code to run for getting correction term out of Se_sn+1 defect
#from pymatgen.io.vaspio.vasp_output import Locpot
#print 'load locpots'
#sdef=Locpot.from_file("LOCPOT")
#spure=Locpot.from_file("../pure/LOCPOT")
#print 'good, now run code'
#s1=DefectCorrectionFreysoldt(spure,sdef,1,47.852,[0.0833330000000032,0.0307343554185392,0.3830636916206969],520)   #this is for Se_sn+1
#s1.run_correction()
