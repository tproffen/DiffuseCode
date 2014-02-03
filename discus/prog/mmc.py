#!/usr/bin/env python
from discuspy import *
import numpy as np

class MMC(object):
    """
    MMC calculation.
    """
    def __init__(self,structure):
        self.structure=structure



class vector(object):
    def __init__(self,iv,is1,is2,dx,dy,dz):
        self.iv  = iv
        self.is1 = is1
        self.is2 = is2
        self.dx  = dx
        self.dy  = dy
        self.dz  = dz

class neighbours(object):
    def __init__(self):
        self.neig = []
        self.target = 'corr'
    def add_neig(self):
        pass
    def del_neig(self):
        pass
    def set_target(self,ic,what,etc):
        pass
