#!/usr/bin/env python
from discuspy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from copy import copy
import numpy as np
setup()

class structure(object):
    """
    The crystal structure.
    """
    def __init__(self,name="crystal"):
        self.ncatoms = 1
        self.icc     = np.ones(3,dtype=np.int32)
        self.v       = 1.
        self.gten    = np.asfortranarray(np.identity(3,dtype=np.float32))
        self.scat    = None
        self.name    = str(name)
        self.spcgr   = "P1"
        self.at_lis  = None
        self.nscat   = 0
        self.natoms  = 0
        self.n_real_atoms = 0
        self.a0      = np.ones(3,dtype=np.float32)
        self.win     = np.ones(3,dtype=np.float32)*90
        self.dw      = None
        self.iscat   = None
        self.pos     = None

    def __str__(self):
        return self.name

    def getName(self):
        return self.name

    def listAtoms(self):
        print 'Name','x'.rjust(5),'y'.rjust(12),'z'.rjust(12),'B'.rjust(12)
        for i in range(self.natoms):
            print self.at_lis[self.iscat[i]].ljust(4), \
                '{:12.6f}'.format(self.pos[0][i]), \
                '{:12.6f}'.format(self.pos[1][i]), \
                '{:12.6f}'.format(self.pos[2][i]), \
                '{:12.6f}'.format(self.dw[self.iscat[i]])

    def getStructure(self):
        self.ncatoms,self.icc,self.v, \
            self.gten,self.scat,self.name, \
            self.spcgr,self.at_lis,self.nscat, \
            self.natoms,self.n_real_atoms,self.a0, \
            self.win,self.dw,self.iscat, \
            self.pos = get_cr()

    def setStructure(self):
        set_cr(self.ncatoms,self.icc,self.v, \
                   self.gten,self.scat,self.name, \
                   self.spcgr,self.at_lis,self.nscat, \
                   self.natoms,self.n_real_atoms,self.a0, \
                   self.win,self.dw,self.iscat, \
                   self.pos)

    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('x label')
        ax.set_ylabel('y label')
        ax.set_zlabel('z label')
        c_table = ['k','b','g','r','c','m','y']
        c=[]
        for n in range(self.natoms):
            c.append(c_table[self.iscat[n]%7])
        ax.scatter(self.pos[0][:self.natoms],
                   self.pos[1][:self.natoms],
                   self.pos[2][:self.natoms],c=c)
        fig.show()

    def readCell(self,fname,dim1=1,dim2=1,dim3=1):
        read_cell(fname,dim1,dim2,dim3)
        self.getStructure()

    def metric(self):
        print 'Lattice constants :'
        print '    '+'a'.ljust(11)+'b'.ljust(11)+'c'.ljust(11)\
            +'alpha'.ljust(11)+'beta'.ljust(11)+'gamma'.ljust(11)+'volume'
        print '  '+'{:9.5f}'.format(self.a0[0])+'  '+'{:9.5f}'.format(self.a0[1]) \
            +'  '+'{:9.5f}'.format(self.a0[2])+'  '+'{:9.5f}'.format(self.win[0]) \
            +'  '+'{:9.5f}'.format(self.win[1])+'  '+'{:9.5f}'.format(self.win[2]) \
            +'  '+'{:9.5f}'.format(self.v)
        print
        print 'Metric Tensor     :'
        for i in range(3):
            print '  '+'{:11.5f}'.format(self.gten[0][i]) \
                +'  '+'{:11.5f}'.format(self.gten[1][i]) \
                +'  '+'{:11.5f}'.format(self.gten[2][i])

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
        self.poly      = np.zeros(5,dtype=np.int32)
        self.bin       = 1
        self.finite    = 0
        self.radiation = 1
        self.poly_n    = 0
        self.lxray     = False
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

    def show(self,what='all'):
        if what == 'all':
            print 'Current PDF calculation settings'.ljust(30)+':'
            print '  Maximum r [A]'.ljust(30)+':',self.rmax
            print '  Grid size DR [A]'.ljust(30)+':',self.deltar,'('+str(self.bin)+')'
            print '  Radiation'.ljust(30)+':',getRadiation(self.radiation),\
                '(at '+str(self.xq)+' A**-1)'
            print
            print '  Applied corrections'.ljust(30)+':'
            if self.qmax == 0:
                print '    Q termination (SINC)'.ljust(30)+': not applied'
            else:
                print '    Q termination (SINC)'.ljust(30)+': Qmax =',self.qmax,' A**-1'
            if self.sigmaq == 0:
                print '    Instrument resolution'.ljust(30)+': not applied'
            else:
                print '    Instrument resolution'.ljust(30)+':',self.sigmaq
            if self.chem_period[0] and self.chem_period[1] and self.chem_period[2]:
                if self.d2d:
                    print '    Periodic boundaries'.ljust(30)+': applied - 2D'
                else:
                    print '    Periodic boundaries'.ljust(30)+': applied - 3D'
                if self.finite == 0:
                    print '    Particle size'.ljust(30)+': infinite'
                elif self.finite == 2:
                    print '    Particle size is sphere'.ljust(30)+':',self.sphere, \
                        'A diameter'
                else:
                    print '    Periodic boundaries'.ljust(30)+': illegal'
            else:
                print '    Periodic boundaries'.ljust(30)+': not applied'
                if self.finite == 0:
                    print '    4 Pi Rho r correction'.ljust(30)+': none'
                elif self.finite == 2:
                    print '    4 Pi Rho r correction'.ljust(30)+': ',self.sphere
                elif self.finite == 1:
                    print 'TO DO !!!!!'
                elif self.finite == 3:
                    print '    4 Pi Rho r correction'.ljust(30)+': treated by tanh'
                    print '       Particle diameter'.ljust(30)+':',self.diam
                    print '       Particle shape param.'.ljust(30)+':',self.shape
            print '    Weight correction'.ljust(30)+':',self.scale
            if self.gauss:
                print '    Convolution therm. Gauss.'.ljust(30)+': applied'
                print '    Resolution broadening'.ljust(30)+':',self.qalp
                print '    Quad. correlation fac.'.ljust(30)+':',self.delta
                print '    Linear correlation fac.'.ljust(30)+':',self.gamma
                print '    PDF peak width ratio'.ljust(30)+':',self.srat,'below',self.rcut
            else:
                print '    Convolution therm. Gauss.'.ljust(30)+': not applied'
                print '    Resolution broadening'.ljust(30)+':',self.qalp
            if self.lrho0:
                print '    Number density'.ljust(30)+': automatic'
            else:
                if self.lrho0_rel:
                    print '    Number density'.ljust(30)+':',self.rho0,'(Mode: relative)'
                else:
                    print '    Number density'.ljust(30)+':',self.rho0,'(Mode: absolute)'
            print '    Correction for RHO0'.ljust(30)+':',self.dnorm
            if self.lexact:
                print '    PDF calculation mode'.ljust(30)+': all atoms'
            else:
                print '    PDF calculation mode'.ljust(30)+': neighboring unit cells'
            print
            print '  Refinement settings'.ljust(30)+':'
            print '    Fit minimum r [A]'.ljust(30)+':',self.rfmin, \
                '(pt.'+str(int(self.rfmin/self.deltar))+')'
            print '    Fit maximum r [A]'.ljust(30)+':',self.rfmax, \
                '(pt.'+str(int(self.rfmax/self.deltar))+')'
            print '    Current scale factor'.ljust(30)+':',self.skal, \
                '(refined = '+'str(self.doskal)'+')'
            print
            print 'Selected atoms for PDF calculation :'
            print '  ...'
            print
            if self.lweights:
                print 'Weights for PDF calculation (user defined) :'
            else:
                print 'Weights for PDF calculation (internal) :'
            print '  ...'
        else:
            menu = {'rmax':       self.rmax,
                    'qmax':       self.qmax,
                    'deltar':     self.deltar,
                    'skal':       self.skal,
                    'qsig':       self.sigmaq,
                    'wavelength': self.xq,
                    'rfmin':      self.rfmin,
                    'rfmax':      self.rfmax,
                    'delta':      self.delta,
                    'rcut':       self.rcut,
                    'srat':       self.srat,
                    'gamma':      self.gamma,
                    'qalp':       self.qalp,
                    'dnorm':      self.dnorm,
                    'rho0':       self.rho0,
                    'sphere':     self.sphere,
                    'diam_poly':  self.diam_poly,
                    'diam':       self.diam,
                    'shape':      self.shape,
                    'weight':     self.scale,
                    'poly':       self.poly,
                    'bin':        self.bin,
                    'finite':     self.finite,
                    'rad':        self.radiation,
                    'gauss':      self.gauss,
                    '2d':         self.d2d,
                    'lweights':   self.lweights,
                    'lrho0':      self.lrho0,
                    'lexact':     self.lexact,
                    'lrho0_rel':  self.lrho0_rel}
            if menu.has_key(what):
                print what,'=',menu[what]

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
                self.poly_n)
        self.set_f_logical()

    def set_f_logical(self):
        set_pdf_logical(self.lxray,self.gauss,self.d2d,
                        self.lweights,self.lrho0,self.lexact,
                        self.lrho0_rel,self.chem_period)

    def plot(self):
        plt.plot(self.r,self.pdf)
        plt.ylabel('G(r) ($\AA^{-2}$)')
        plt.xlabel('r ($\AA$)')
        plt.show()

    def setBound(self,t='crystal',d='3D'):
        if t[0:3].lower()=='per':
            self.chem_period=[True]*3
            self.lexact=False
            if d.lower()=='2d':
                self.d2d=True
                print 'Setting PDF calculation to','periodic bound. 2D, unit cell','mode ..'
            else:
                self.d2d=False
                print 'Setting PDF calculation to','periodic bound. 3D, unit cell','mode ..'
        else:
            self.chem_period=[False]*3
            if t[0:2].lower()=='ex':
                self.lexact=True
                print 'Setting PDF calculation to','no periodic bound., exact','mode ..'
            else:
                self.lexact=False
                print 'Setting PDF calculation to','no periodic bound., unit cell','mode ..'

    def setFinite(self,t='periodic',p=0):
        self.finite=getFiniteNumber(t)
        if self.finite == getFiniteNumber('Sphere'):
            self.sphere = p
        elif self.finite == getFiniteNumber('Polygon'):
            self.diam_poly = p

    def setDensity(self,den):
        self.lrho0_rel=False
        if type(den) is str:
            if den[0:2].lower()=='au':
                self.lrho0=True
            else:
                print 'Did nothing!'
        else:
            self.rho0  = den
            self.lrho0 = False

    def setRDensity(self,den):
        self.lrho0_rel = True
        if type(den) is str:
            if den[0:2].lower()=='au':
                self.lrho0=True
            else:
                print 'Did nothing!'
        else:
            self.rho0  = den
            self.lrho0 = False

    def setRad(self,rad,wav=-1):
        if rad[0].lower()=='n':
            self.radiation=2
            #self.xq=wav
            self.lxray=False
        elif rad[0].lower()=='x':
            self.radiation=1
            if wav != -1:
                self.xq=wav
            self.lxray=True
        elif rad[0].lower()=='e':
            self.radiation=3
            if wav != -1:
                self.xq=wav
            self.lxray=True
        else:
            print "Invalid radiation"
            exit
        print "Setting radiation to",getRadiation(self.radiation),
        if self.lxray:
            print "with wavelength",self.xq,"A"

    def setRange(self,rmax,deltar):
        self.rmax   = rmax
        self.deltar = deltar
        self.rfmin  = deltar
        self.rfmax  = rmax
        self.bin    = int(rmax/deltar)+1

    def setFRange(self,min,max):
        self.rfmin = min
        self.rfmax = max

    def setPolygon(self,p):
        self.poly_n = len(p)
        self.poly[:len(p)] = p

    def setSRatio(self,srat,rcut):
        self.srat = srat
        self.rcut = rcut

    def setThermal(self,gauss):
        if gauss[0:3].lower()=='gau':
            self.gauss=True
        else:
            self.gauss=False

    def set(self,choice,value):
        #print choice,menu[choice]
        #menu[choice]=value
        pass
        

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

def getRadiation(r):
    rad = {1:"X-Rays", 2:"Neutrons", 3:"Electrons"}
    if rad.has_key(r): 
        return rad[r]
    else:
        return "Invalid"

def getFiniteNumber(f):
    if f[0:3].lower()=='pol':
        return 1
    elif f[0:3].lower()=='sph':
        return 2
    elif f[0:3].lower()=='tan':
        return 3
    else:
        return 0

def getFinite(f):
    finite = {0:"Periodic", 1:"Polygon", 2:"Sphere", 3:"tanh"}
    if finite.has_key(f):
        return finite[f]
    else:
        return "Invalid"

