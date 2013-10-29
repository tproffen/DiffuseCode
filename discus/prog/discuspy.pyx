import numpy as np

cdef extern int cr_natoms
cdef extern int nmax

cdef extern:
    void setup_c()

cdef extern:
    void discus_loop()

cdef extern:
    void get_cr_pos_c(float*,int)


def setup():
    setup_c()

def interactive():
    discus_loop()

def get_natoms():
    return cr_natoms

def print_natoms():
    print cr_natoms

def get_nmax():
    return nmax

def print_nmax():
    print nmax

def set_nmax(n):
    global nmax
    nmax = n

def get_cr_pos():
    out = np.empty((3,nmax),dtype=np.float32)
    out = np.asfortranarray(out)
    cdef float [::1,:] out2 = out
    get_cr_pos_c(&out2[0][0],nmax)
    return out

