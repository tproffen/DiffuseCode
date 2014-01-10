#!/usr/bin/env python
import numpy as np
from func import getRadiation, getRadiationPar
import matplotlib.pyplot as plt

class Powder(object):
    """
    Powder calculation.
    """
    def __init__(self,structure,name="Powder"):
        self.name          = name
        self.structure     = structure
        self.powder        = []
        self.r             = []
        self.MAXPKT        = 1
        self.axis          = 1
        self.four_mode     = 0
        self.four_type     = 0
        self.pow_four_vers = 4
        self.pow_lp        = 1
        self.pow_l_all     = True
        self.tthmin        = 0.1
        self.tthmax        = 40.0
        self.deltatth      = 0.05
        self.qmin          = 0.2
        self.qmax          = 7.0
        self.deltaq        = 0.0001
        self.ds_max        = 0.0001
        self.ds_min        = 0.0001
        self.delta         = 0.0
        self.lp_fac        = 0.88
        self.lp_ang        = 20.0
        self.lp_cos        = 0.936
        self.back          = [0.0]*6
        self.scale         = 1.0
        self.hkl_max       = [4.0]*3
        self.hkl_del       = [0.05]*3
        self.hkl_shift     = [0.00]*3
        self.pref          = False
        self.pref_type     = 1
        self.pref_g1       = 0.0
        self.pref_g2       = 0.0
        self.pref_hkl      = [0., 0., 1.]
        self.profile       = 2
        self.pow_pr_par    = 0
        self.pow_fwhm      = 0.01
        self.pow_eta       = 0.5
        self.pow_etax      = 0.0
        self.pow_u         = 0.0
        self.pow_v         = 0.0
        self.pow_w         = 0.05
        self.pow_p1        = 0.0
        self.pow_p2        = 0.0
        self.pow_p3        = 0.0
        self.pow_p4        = 0.0
        self.pow_width     = 20.0
        self.size_of       = 0
        self.radiation     = 1
        self.lxray         = True
        self.ano           = False
        
    def __str__(self):
        return self.name+"\nStructure: "+self.structure.getName()
    
    def calc(self):
        pass
    
    def show(self):
        pass
    
    def show_f(self):
        pow_show()
    
    def set_f(self):
        pass
    
    def set_f_logical(self):
        pass
    
    def plot(self):
        plt.plot(self.r,self.powder)
        plt.ylabel('Intensity')
        plt.xlabel('???2theta (1/$\AA$)???')
        plt.show()
    
    def setRad(self,rad):
        self.radiation,self.lxray = getRadiationPar(rad)
        print "Setting radiation to",getRadiation(self.radiation)

