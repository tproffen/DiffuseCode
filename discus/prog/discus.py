#!/usr/bin/env python
from discuspy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
setup()

class structure(object):
    """
    The crystal structure.
    """
    def __init__(self,name="crystal",spacegroup="P1"):
        self.name=str(name)
        self.spcgr = spacegroup
        self.natoms=0
        self.pos=[]
        self.iscat=[]
        self.at_lis=dict()
    def __str__(self):
        return self.name
    def getName(self):
        return self.name
    def listAtoms(self):
        for i in range(self.natoms):
            print self.at_lis[self.iscat[i]], \
                self.pos[0][i],self.pos[1][i],self.pos[2][i]
    def getStructure(self):
        self.pos=get_cr_pos()
        self.natoms=get_natoms()
        self.name=get_crystal_name()
        self.spcgr=get_crystal_spcgr()
        self.at_lis=get_cr_at_lis()
        self.iscat=get_cr_iscat()
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('x label')
        ax.set_ylabel('y label')
        ax.set_zlabel('z label')
        c_table = ['k','b','g','r','c','m','y']
        c=[]
        for n in range(self.natoms):
            c.append(c_table[self.iscat[n]%8])
        ax.scatter(self.pos[0][:self.natoms],
                   self.pos[1][:self.natoms],
                   self.pos[2][:self.natoms],c=c)
        fig.show()
    def readCell(self,fname,dim1=1,dim2=1,dim3=1):
        read_cell(fname,dim1,dim2,dim3)
        self.getStructure()

class pdf(object):
    """
    PDF calculation.
    """
    def __init__(self,structure,name="PDF"):
        self.name      = name
        self.structure = structure
        self.pdf       = []
        self.r         = []
        self.rmax      = 50.0
        self.qmax      = 30.0
        self.deltar    = 0.01
        self.skal      = 1.0
        self.sigmaq    = 0.0
        self.xq        = 0.0
        self.rfmin     = 0.05
        self.rfmax     = 15.0
        self.delta     = 0.0
        self.rcut      = 0.0
        self.srat      = 1.0
        self.gamma     = 0.0
        self.qalp      = 0.0
        self.dnorm     = 1.0
        self.rho0      = 0.0
        self.sphere    = 0.0
        self.diam_poly = 0.0
        self.diam      = 0.0
        self.shape     = 0.0
        self.scale     = 1.0
        self.poly      = [0]*5
        self.bin       = 1
        self.finite    = 0
        self.radiation = 1
        self.gauss     = False
        self.d2d       = False
        self.lweights  = False
        self.lrho0     = True
        self.lexact    = False
        self.lrho0_rel = False        
        self.chem_period=[True]*3
    def __str__(self):
        return self.name+"\nStructure: "+self.structure.getName()
    def setup(self):
        pdf_setup()
    def calc(self):
        self.set_f()
        self.setup()
        pdf_determine()
        [self.r,self.pdf]=get_pdf()
    def show(self):
        self.setup()
        print 'Current PDF calculation settings :'
        print '  Maximum r [A]              : '+str(self.rmax)
        print '  Grid size DR [A]           : '+str(self.deltar) \
            + ' ('+str(self.bin)+')'
        print '  Radiation                  : '+str(getRadiation(self.radiation)) \
            + ' (at '+str(self.xq)+' A**-1)'
        print
        print '  Applied corrections        :'
        if self.qmax == 0:
            print '    Q termination (SINC)     : not applied'
        else:
            print '    Q termination (SINC)     : Qmax = '+str(self.qmax)+' A**-1'
        if self.sigmaq == 0:
            print '    Instrument resolution    : not applied'
        else:
            print '    Instrument resolution    : '+str(self.sigmaq)
        if self.chem_period[0] and self.chem_period[1] and self.chem_period[2]:
            if self.d2d:
                print '    Periodic boundaries      : applied - 2D'
            else:
                print '    Periodic boundaries      : applied - 3D'
            if self.finite == 0:
                print '    Particle size            : infinite'
            elif self.finite == 2:
                print '    Particle size is sphere  : '+str(self.pdf_sphere) \
                    +' A diameter'
            else:
                print '    Periodic boundaries      : illegal'
        else:
            print '    Periodic boundaries      : not applied'
            if pdf_finite == 0:
                print '    4 Pi Rho r correction    : none'
            elif pdf_finite == 2:
                print '    4 Pi Rho r correction    : '+str(self.sphere)
            elif pdf_finite == 1:
                print 'TO DO !!!!!'
            elif pdf_finite == 3:
                print '    4 Pi Rho r correction    : treated by tanh'
                print '       Particle diameter     : '+str(self.diam)
                print '       Particle shape param. : '+str(self.shape)
        print '    Weight correction        : '+str(self.scale)
        if self.gauss:
            print '    Convolution therm. Gauss.: applied'
            print '    Resolution broadening    : '+str(self.qalp)
            print '    Quad. correlation fac.   : '+str(self.delta)
            print '    Linear correlation fac.  : '+str(self.gamma)
            print '    PDF peak width ratio     : '+str(self.srat) \
                + ' below '+str(self.rcut)
        else:
            print '    Convolution therm. Gauss.: not applied'
            print '    Resolution broadening    : '+str(self.qalp)
        if self.lrho0:
            print '    Number density           : automatic'
        else:
            if self.lrho0_rel:
                print '    Number density           : '+str(self.rho0) \
                    + ' (Mode: relative)'
            else:
                print '    Number density           : '+str(self.rho0) \
                    + ' (Mode: absolute)'
        print '    Correction for RHO0      : '+str(self.dnorm)
        if self.lexact:
            print '    PDF calculation mode     : all atoms'
        else:
            print '    PDF calculation mode     : neighboring unit cells'
        print
        print '  Refinement settings        :'
        print '    Fit minimum r [a]        : '+str(self.rfmin) \
            + '(pt.'+str(int(self.rfmin/self.deltar))+')'
        print '    Fit maximum r [A]        : '+str(self.rfmax) \
            +'(pt.'+str(int(self.rfmax/self.deltar))+')'
        print '    Current scale factor     : '+str(self.skal) \
            +'(refined = '+'str(self.doskal)'+')'
        print
        print 'Selected atoms for PDF calculation :'
        print '  ...'
        print
        if self.lweights:
            print 'Weights for PDF calculation (user defined) :'
        else:
            print 'Weights for PDF calculation (internal) :'
        print '  ...'
    def show_f(self):
        pdf_show()
    def set_f(self):
        set_pdf(self.rmax,self.qmax,self.deltar,
                self.skal,self.sigmaq,self.xq,
                self.rfmin,self.rfmax,self.delta,
                self.rcut,self.srat,self.gamma,
                self.qalp,self.dnorm,self.rho0,
                self.sphere,self.diam_poly,self.diam,
                self.shape,self.scale,self.poly,
                self.bin,self.finite,self.radiation,
                self.gauss,self.d2d,self.lweights,
                self.lrho0,self.lexact,self.lrho0_rel)
    def plot(self):
        plt.plot(self.r,self.pdf)
        plt.ylabel('G(r) ($\AA^{-2}$)')
        plt.xlabel('r ($\AA$)')
        plt.show()
    def range(self,rmax,deltar):
        self.rmax=rmax
        self.deltar=deltar
        self.frmin=deltar
        self.rfmax=rmax
        self.bin=int(rmax/deltar)+1
    def set(self,choice,value):
        menu = {'rmax':self.rmax,
                'qmax':self.qmax,
                'deltar':self.deltar,
                'skal':self.skal,
                'qsig':self.sigmaq,
                'wavelgnth':self.xq,
                'rfmin':self.rfmin,
                'rfmax':self.rfmax,
                'delta':self.delta,
                'rcut':self.rcut,
                'srat':self.srat,
                'gamma':self.gamma,
                'qalp':self.qalp,
                'dnorm':self.dnorm,
                'rho0':self.rho0,
                'sphere':self.sphere,
                'diam_poly':self.diam_poly,
                'diam':self.diam,
                'shape':self.shape,
                'weight':self.scale,
                'poly':self.poly,
                'bin':self.bin,
                'finite':self.finite,
                'rad':self.radiation,
                'gauss':self.gauss,
                '2d':self.d2d,
                'lweights':self.lweights,
                'lrho0':self.lrho0,
                'lexact':self.lexact,
                'lrho0_rel':self.lrho0_rel}
        print menu[choice],choice,value
        menu[choice]=value
        

def getRadiation(r):
    rad = {1:"X-Rays", 2:"Neutron", 3:"Electrons"}
    if r in rad: 
        return rad[r]
    else:
        return "Invalid"


class powder(object):
    """
    Powder calculation.
    """
    def __init__(self,structure):
        self.structure=structure

class mmc(object):
    """
    MMC calculation.
    """
    def __init__(self,structure):
        self.structure=structure

class fourier(object):
    """
    Fourier calculation.
    """
    def __init__(self,structure):
        self.structure=structure
