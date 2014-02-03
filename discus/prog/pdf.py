#!/usr/bin/env python
from discuspy import *
import numpy as np
from func import getRadiation, getRadiationPar
import matplotlib.pyplot as plt

class PDF(object):
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
        self.radiation,self.lxray = getRadiationPar(rad)
        if self.radiation != 2:
            if wav != -1:
                self.xq=wav
        print "Setting radiation to",getRadiation(self.radiation),
        if self.lxray:
            print "with wavelength",self.xq,"A"

    def setRange(self,rmax=None,deltar=None):
        if rmax==None or deltar==None:
            print "To change enter setRange( 'Rmax' , 'delta R' )"
            print 'Current:'
        else:
            self.rmax   = rmax
            self.deltar = deltar
            self.rfmin  = deltar
            self.rfmax  = rmax
            self.bin    = int(rmax/deltar)+1
        print '  Maximum r [A]:',self.rmax
        print '  Grid size DR [A]:',self.deltar,'('+str(self.bin)+')'

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

    def show(self,what='all'):
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

