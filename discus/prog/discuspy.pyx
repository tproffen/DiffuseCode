import numpy as np
#from libcpp cimport bool as c_bool

### config_mod ###
cdef extern int maxscat
cdef extern int nmax
### end config_mod ###

### crystal_mod ###
cdef extern int cr_natoms
cdef extern int cr_nscat
### endcrystal_mod ###

### pdf_mod ###
cdef extern int pdf_maxscat
cdef extern int pdf_maxdat
cdef extern int pdf_maxbnd
cdef extern int pdf_nscat
cdef extern int pdf_ndat
cdef extern int pdf_nbnd
cdef extern float pdf_rmax
cdef extern float pdf_qmax
cdef extern float pdf_deltar
cdef extern float pdf_skal
cdef extern float pdf_sigmaq
cdef extern float pdf_xq
cdef extern float pdf_rfmin
cdef extern float pdf_rfmax
cdef extern float pdf_delta
cdef extern float pdf_rcut
cdef extern float pdf_srat
cdef extern float pdf_gamma
cdef extern float pdf_qalp
cdef extern float pdf_dnorm
cdef extern float pdf_rho0
cdef extern float pdf_sphere
cdef extern float pdf_diam_poly
cdef extern float pdf_diam
cdef extern float pdf_shape
cdef extern float pdf_scale
cdef extern float pdf_poly[5]
cdef extern int pdf_bin
cdef extern int pdf_finite
cdef extern int pdf_radiation
#cdef extern bool pdf_lxray
#cdef extern bool pdf_gauss
#cdef extern bool pdf_2d
#cdef extern bool pdf_lweights
#cdef extern bool pdf_lrho0
#cdef extern bool pdf_lexact
#cdef extern bool pdf_lrho0_rel
### end pdf_mod ###


cdef extern void setup_c "setup" ()
cdef extern void discus_loop()
cdef extern void get_cr_pos_c "get_cr_pos" (float*,int)
cdef extern void set_cr_pos_c "set_cr_pos" (float*,int)
cdef extern void get_cr_iscat_c "get_cr_iscat" (int*,int)
cdef extern void get_crystal_name_f "get_crystal_name" (char*)
cdef extern void get_crystal_spcgr_f "get_crystal_spcgr" (char*)
cdef extern void get_pdf_c "get_pdf" ( double*,int,int)
cdef extern void pdf_determine_c (bool)
cdef extern void pdf_show_c(char*,int)
cdef extern void pdf_setup_c "pdf_setup" ()
cdef extern void rese_cr_f "rese_cr" ()
cdef extern void read_cell_f "read_cell" (char*,int,int,int,int)
cdef extern void get_atom_type(char*,int)
cdef extern void alloc_pdf_f ()
cdef extern void set_pdf_logical_f "set_pdf_logical" (bool,bool,bool,bool,bool,bool,bool)

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
    foo = np.empty((3,nmax),dtype=np.float32,order='F')
    #foo = np.asfortranarray(foo)
    cdef float [::1,:] bar = foo
    get_cr_pos_c(&bar[0][0],nmax)
    return foo

def set_cr_pos(pos_in):
    cdef float [::1,:] cr_pos_in = np.asfortranarray(pos_in)
    set_cr_pos_c(&cr_pos_in[0][0],nmax)

def get_cr_iscat():
    foo = np.empty((nmax),dtype=np.int32,order='F')
    cdef int [:] bar = foo
    get_cr_iscat_c(&bar[0],nmax)
    return foo

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

def pdf_determine():
    #cdef c_bool a = True
    pdf_determine_c(True)

def pdf_show():
    cdef char* what = "PDF"
    pdf_show_c(&what[0],len(what))

def pdf_setup():
    pdf_setup_c()

def get_crystal_name():
    cdef char[81] c_str
    get_crystal_name_f(&c_str[0])
    cdef bytes py_str = c_str
    return py_str

def get_crystal_spcgr():
    cdef char[81] c_str
    get_crystal_spcgr_f(&c_str[0])
    cdef bytes py_str = c_str
    return py_str

def read_cell(fname,dim1,dim2,dim3):
    cdef char* c_str = fname
    cdef int fname_length = len(fname)
    read_cell_f(&c_str[0],fname_length,dim1,dim2,dim3)

def get_cr_at_lis():
    out = dict()
    cdef char[5] c_str
    for i in range(maxscat+1):
        get_atom_type(&c_str[0],i)
        out[i]=c_str
    return out

def set_pdf(rmax,qmax,deltar,
    skal,sigmaq,xq,
    rfmin,rfmax,delta,
    rcut,srat,gamma,
    qalp,dnorm,rho0,
    sphere,diam_poly,diam,
    shape,scale,poly,
    bin,finite,radiation):
    global pdf_nscat,pdf_ndat,pdf_nbnd,pdf_rmax,pdf_qmax,pdf_deltar,pdf_skal,pdf_sigmaq,pdf_xq,pdf_rfmin,pdf_rfmax,pdf_delta,pdf_rcut,pdf_srat,pdf_gamma,pdf_qalp,pdf_dnorm,pdf_rho0,pdf_sphere,pdf_diam_poly,pdf_diam,pdf_shape,pdf_scale,pdf_poly,pdf_bin,pdf_finite,pdf_radiation
    pdf_rmax=rmax
    pdf_qmax=qmax
    pdf_deltar=deltar
    pdf_skal=skal
    pdf_sigmaq=sigmaq
    pdf_xq=xq
    pdf_rfmin=rfmin
    pdf_rfmax=rfmax
    pdf_delta=delta
    pdf_rcut=rcut
    pdf_srat=srat
    pdf_gamma=gamma
    pdf_qalp=qalp
    pdf_dnorm=dnorm
    pdf_rho0=rho0
    pdf_sphere=sphere
    pdf_diam_poly=diam_poly
    pdf_diam=diam
    pdf_shape=shape
    pdf_scale=scale
    ###pdf_poly=poly
    pdf_bin=bin
    pdf_finite=finite
    pdf_radiation=radiation
    if bin>pdf_maxdat:
        pdf_nscat = max(pdf_nscat, cr_nscat, pdf_maxscat, maxscat)
        pdf_ndat  = max(pdf_ndat , pdf_bin , pdf_maxdat)
        pdf_nbnd  = max(pdf_nbnd ,           pdf_maxbnd)
        alloc_pdf_f()

def set_pdf_logical(lxray,gauss,d2d,
                    lweights,lrho0,lexact,
                    lrho0_rel):
    set_pdf_logical_f(lxray,gauss,d2d,
                      lweights,lrho0,lexact,
                      lrho0_rel)
