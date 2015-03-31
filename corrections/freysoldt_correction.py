#!/usr/bin/env python

"""
This module provides classes to apply typical defect corrections
"""

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "November 4, 2012"

import pylab
import subprocess
import os
from pymatgen.electronic_structure.core import Spin


class DefectCorrectionFreysoldt(object):

    """
    This class apply the Freysoldt correction to remove electrostatic defect interaction and apply
    electostatic alignemnt
    Ideally this wrapper around sxdefectalign should be replaced by a python code implementing
    a generalized anisotropic Freysoldt correction
    """
    
    def __init__(self, locpot_bulk, locpot_defect, charge, epsilon, frac_coords):
        self._locpot_bulk=locpot_bulk
        self._locpot_defect=locpot_defect
        self._charge=charge
        self._epsilon=epsilon
        self._frac_coords=frac_coords
        
        
    def prepare_files(self):
        self._locpot_bulk.write_file("/home/geoffroy/tmp/LOCPOT_vref",True)
        self._locpot_defect.write_file("/home/geoffroy/tmp/LOCPOT_vdef",True)
    
    def plot_hartree_pot(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
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
                       supercell_pot,'r')
            ax.plot(self._locpot_defect.get_axis_grid(axis),self._locpot_defect.get_average_along_axis(axis),'b')
            ax.plot([self._frac_coords[axis]*self._locpot_defect.structure.lattice._lengths[axis]],[0],'or',markersize=4.0)
        pylab.show()
    
    def plot_hartree_pot_diff(self):
        import matplotlib.pyplot as plt
        fig=plt.figure()
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
        pylab.show()
    
    def test_rho(self,axis=2,alignment_sr=0.0):
        for expnorm in [0.0,0.25,0.5,0.75,1.0]:
            command=['/home/geoffroy/research/TCO/defects/sxdefectalign']
            #./sxdefectalign -C 1.6 --vasp --charge -2 --ecut 6 --vref LOCPOT_SnO --vdef LOCPOT_SnO_Vac_O_2
            command.append("--vasp")
            command.append("-a"+str(axis+1))
            command.append("--relative")
            command.append("--pos")
            command.append(str(self._frac_coords))
            #command.append("[0.0,0.25,0.0]")
            #command.append("--expnorm")
            #command.append(str(expnorm))
            command.append("--charge")
            command.append(str(self._charge))
            command.append("--ecut")
            command.append("6")
            command.append("--eps")
            command.append(str(self._epsilon))
            command.append("-C")
            command.append(str(alignment_sr))
            command.append("--vref")
            command.append("/home/geoffroy/tmp/LOCPOT_vref")
            command.append("--vdef")
            command.append("/home/geoffroy/tmp/LOCPOT_vdef")
            
            output = subprocess.call(command, stdout=subprocess.PIPE)
            #output=p.communicate()
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
    
    def compute_correction_rho(self,axis=2,alignment_sr=0.0,plot=False,expnorm=0.25):
        command=['/home/geoffroy/research/TCO/defects/sxdefectalign']
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
    
    
    def plot_pot_diff(self, align=0.0):
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
    
    def compute_correction(self,axis=2,alignment_sr=0.0,plot=False, beta=None):
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

        
    def plot_bulk_potential(self,axis=0):
        data=self._locpot_bulk.get_avg_potential_along_axis(axis)
        #print self._locpot_bulk.data[Spin.up]
        print len(self._locpot_bulk.data[Spin.up])
        print 
        pylab.plot(data[0],data[1])
        pylab.show()
        
    def plot_defect_potential(self,axis=0):
        data=self._locpot_defect.get_avg_potential_along_axis(axis)
        print data
        pylab.plot(data[0],data[1])
        pylab.show()
        
    def plot_compare_potentials(self,axis=0):
        locpot1=self._locpot_bulk
        data1=locpot1.get_avg_potential_along_axis(axis)
        #print len(data1)
        x_ticks1=[x * locpot1.poscar._struct.lattice.c/len(data1) for x in range(0, len(data1))]
        #for i in range(len(x_ticks)):
        #    print str(x_ticks[i])+" "+str(data[i][0])
        locpot_defect=self._locpot_defect
        data_defect=locpot_defect.get_avg_potential_along_axis(axis)
        locpot_bulk=self._locpot_bulk
        data_bulk=locpot_bulk.get_avg_potential_along_axis(axis)
        print locpot_defect.poscar.struct.lattice.abc[axis]
        print locpot_bulk.poscar.struct.lattice.abc[axis]
        print locpot_defect.poscar.struct.lattice.abc[axis]/locpot_bulk.poscar.struct.lattice.abc[axis]
        multiplier=int(round(locpot_defect.poscar.struct.lattice.abc[axis]/locpot_bulk.poscar.struct.lattice.abc[axis]))
        print multiplier
        f1=pylab.figure()
        pylab.plot(data_defect[0],data_defect[1])
        pylab.plot([x+i*locpot_bulk.poscar.struct.lattice.abc[axis] for i in range(multiplier) for x in data_bulk[0] ],[y for i in range(multiplier) for y in data_bulk[1] ],'r')
        from scipy import interpolate
        tck = interpolate.splrep([x+i*locpot_bulk.poscar.struct.lattice.abc[axis] for i in range(multiplier) for x in data_bulk[0] ],[y for i in range(multiplier) for y in data_bulk[1] ],s=0)
        f2=pylab.figure()
        pylab.plot(data_defect[0],data_defect[1]-interpolate.splev(data_defect[0],tck,der=0))
        pylab.show()
        