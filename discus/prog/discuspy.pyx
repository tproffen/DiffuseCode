import numpy as np

cdef extern int cr_natoms
cdef extern int nmax
cdef extern float pdf_deltar
cdef extern float pdf_skal
cdef extern float pdf_rfmin
cdef extern float pdf_rfmax
cdef extern int pdf_ndat

cdef extern void setup_c()
cdef extern void discus_loop()
cdef extern void get_cr_pos_c(float*,int)
cdef extern void set_cr_pos_c(float*,int)
cdef extern void get_pdf_c(double*,int,int)

def setup():
    setup_c()

def interactive():
    discus_loop()

def get_natoms():
    return cr_natoms

def set_natoms(n):
    global cr_natoms
    cr_natoms=n

def get_nmax():
    return nmax

def set_nmax(n):
    global nmax
    nmax = n

def get_cr_pos():
    foo = np.empty((3,nmax),dtype=np.float32)
    foo = np.asfortranarray(foo)
    cdef float [::1,:] bar = foo
    get_cr_pos_c(&bar[0][0],nmax)
    return foo

def set_cr_pos(pos_in):
    cdef float [::1,:] cr_pos_in = np.asfortranarray(pos_in)
    set_cr_pos_c(&cr_pos_in[0][0],nmax)

def get_pdf_ndat():
    return pdf_ndat

def get_pdf():
    nmi = int(pdf_rfmin/pdf_deltar)
    nma = int(pdf_rfmax/pdf_deltar)
    r=[i*pdf_deltar for i in range(nmi,nma)]
    pdf = np.empty(len(r),dtype=np.float)
    cdef double [:] pdf_c = pdf
    get_pdf_c(&pdf_c[0],nmi,len(r))
    return r,pdf
