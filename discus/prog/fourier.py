#!/usr/bin/env python
from discuspy import *
import numpy as np
from func import getRadiation, getRadiationPar
import matplotlib.pyplot as plt

class Fourier(object):
    """
    Fourier calculation.
    """
    def __init__(self,structure,name="Fourier"):
        self.name       = name
        self.structure  = structure
        self.dsi        = np.empty((121,121),dtype=np.float32,order='F')
        self.fave       = 0.0
        self.nlots      = 1
        self.ilots      = 1 # OFF = 1, BOX = 2, ELI = 3
        self.ls_xyz     = [5,5,5]
        self.four_mode  = 1
        self.lperiod    = True
        self.inc        = [121,121,1]
        self.ano        = False
        self.ldbw       = False
        self.lxray      = True
        self.radiation  = 1
        self.eck        = np.zeros((4,3),dtype=np.float32,order='F')
        self.eck[1,0]   = self.eck[2,1] = 5.0
        self.vi         = np.zeros((3,3),dtype=np.float32,order='F')
        self.vi[0,0]    = self.vi[1,1] = 0.05
        self.rlambda    = 0.709260
        self.ltop       = False

    def __str__(self):
        return self.name+"\nStructure: "+self.structure.getName()

    def ll(self,e11,e12,e13):
        self.eck[0,:] = [e11, e12, e13]
        for i in range(3):
            self.vi[0,i]  = (self.eck[1,i]-self.eck[0,i])/max(1,self.inc[0]-1)
            self.vi[1,i]  = (self.eck[2,i]-self.eck[0,i])/max(1,self.inc[1]-1)
            self.vi[2,i]  = (self.eck[3,i]-self.eck[0,i])/max(1,self.inc[2]-1)

    def lr(self,e21,e22,e23):
        self.eck[1,:] = [e21, e22, e23]
        for i in range(3):
            self.vi[0,i]  = (self.eck[1,i]-self.eck[0,i])/max(1,self.inc[0]-1)

    def ul(self,e31,e32,e33):
        self.eck[2,:] = [e31, e32, e33]
        for i in range(3):
            self.vi[1,i]  = (self.eck[2,i]-self.eck[0,i])/max(1,self.inc[1]-1)

    def tl(self,e41,e42,e43):
        self.eck[3,:] = [e41, e42, e43]
        for i in range(3):
            self.vi[2,i]  = (self.eck[3,i]-self.eck[0,i])/max(1,self.inc[2]-1)
        self.ltop = True

    def na(self,inc1):
        self.inc[0]=inc1
        for i in range(3):
            self.vi[0,i]  = (self.eck[1,i]-self.eck[0,i])/max(1,self.inc[0]-1)

    def no(self,inc2):
        self.inc[1]=inc2
        for i in range(3):
            self.vi[1,i]  = (self.eck[2,i]-self.eck[0,i])/max(1,self.inc[1]-1)

    def layer(self,e11,e12,e13,e21,e22,e23,e31,e32,e33,inc1,inc2):
        self.inc[0]=inc1
        self.inc[1]=inc2
        self.eck[0,:] = [e11, e12, e13]
        self.eck[1,:] = [e21, e22, e23]
        self.eck[2,:] = [e31, e32, e33]
        for i in range(3):
            self.vi[0,i]  = (self.eck[1,i]-self.eck[0,i])/max(1,self.inc[0]-1)
            self.vi[1,i]  = (self.eck[2,i]-self.eck[0,i])/max(1,self.inc[1]-1)
            self.vi[2,i]  = (self.eck[3,i]-self.eck[0,i])/max(1,self.inc[2]-1)

    def get_dsi(self):
        self.dsi=get_diffuse_dsi(self.inc[0],self.inc[1],self.inc[2])

    def lots(self,t,lx=5,ly=5,lz=5,n=1,p='yes'):
        self.nlots  = n
        self.ls_xyz = [lx,ly,lz]
        if t[0].lower() == 'b':
            self.ilots = 2 # Box
        elif t[0].lower() == 'e':
            self.ilots = 3 # Ellipsoid
        if p[0].lower() == 'y':
            self.lperiod = True
        else:
            self.lperiod = False

    def lots_off():
        self.nlots = 1 # Off
        self.ilots = 1

    def set_aver(self,ave):
        self.fave = ave*0.01

    def setRad(self,rad,wav=-1):
        self.radiation,self.lxray = getRadiationPar(rad)
        if wav != -1:
            self.rlambda=wav
        print "Setting radiation to",getRadiation(self.radiation),
        print "with wavelength",self.rlambda,"A"

    def temp(self,use='None'):
        """
        Usage: temp(<'use'|'ignore'>)
        
        """
        if use==None:
            print '   Temp. factors:',
            if self.ldbw:
                print 'used'
            else:
                print 'ignored'
        elif use[0:2].lower()=='us': # Use
            self.ldbw = True
        elif use[0:2].lower()=='ig': # Ignore
            self.ldbw = False

    def disp(self,anom='p'):
        """
        print "Usage: disp(<'anom'|'off'>)"
        
        """
        if anom[0].lower()=='a':
            self.ano=True
        elif anom[0].lower()=='o':
            self.ano=False
        else:
            print '   Anomalous scat.    :',
            if self.ano:
                print 'used'
            else:
                print 'ignored'

    def delf(self):
        pass

    def set_f(self):
        set_fourier(self.fave,self.nlots,self.ilots,self.ls_xyz,\
                        self.inc,self.radiation,self.eck,self.vi,self.rlambda)

    def set_f_logical(self):
        set_diffuse_logical(self.lperiod,self.lxray)

    def run(self):
        self.set_f()
        self.set_f_logical()
        set_ltop(self.ltop)
        #dlink()
        four_run()
        self.get_dsi()

    def plot(self,lim=-1):
        fig = plt.figure()
        imgplot = plt.imshow(self.dsi)
        if lim > 0:
            imgplot.set_clim(0,lim)
        fig.show()

    def show_f(self):
        four_show_f()

    def show(self):
        if self.fave==0:
            print ' Fourier technique    : turbo Fourier'
        else:
            print ' Fourier technique    : turbo Fourier, minus <F> (based on',\
                self.fave*100,'% of cryst.)'
        if self.four_mode==1:
            print ' Fourier calculated by: atom form factors'
        elif self.four_mode==0:
            print ' Fourier calculated by: object form factors'
        if self.ilots==0:
            print '   Fourier volume     : complete crystal'
        elif self.ilots==1:
            print '   Fourier volume     :',self.nlots,' box shaped lots'
            print '   Lot size           :',self.ls_xyz[0],' x ',self.ls_xyz[1],\
                ' x ',self.ls_xyz[2],' unit cells (periodic boundaries = ',self.lperiod,')'
        elif self.ilots==2:
            print '   Fourier volume     :',self.nlots,' ellipsoid shaped lots'
            print '   Lot size           :',self.ls_xyz[0],' x ',self.ls_xyz[1],\
                ' x ',self.ls_xyz[2],' unit cells (periodic boundaries = ',self.lperiod,')'
        print '  Radiation           :',getRadiation(self.radiation),\
            '(at '+str(self.rlambda)+' A**-1)'
        print '   Temp. factors      :',
        if self.ldbw:
            print 'used'
        else:
            print 'ignored'
        print '   Anomalous scat.    : ',
        if self.ano:
            print 'used'
        else:
            print 'ignored'
        print
        print ' Reciprocal layer     :'
        print '   lower left  corner :',self.eck[0]
        print '   lower right corner :',self.eck[1]
        print '   upper left  corner :',self.eck[2]
        print '   top   left  corner :',self.eck[3]
        print
        print '   hor. increment     :',self.vi[0],' -> ','A**-1'
        print '   vert. increment    :',self.vi[1],' -> ','A**-1'
        print '   top   increment    :',self.vi[2],' -> ','A**-1'
        print '   # of points        :',self.inc[0],'x',self.inc[1],'x',self.inc[2],'(h,k,l)'
        print '   Angle Ratio Aver '
