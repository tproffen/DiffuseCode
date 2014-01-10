#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Structure(object):
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

    def readStru(self,fname):
        read_stru(fname)
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

