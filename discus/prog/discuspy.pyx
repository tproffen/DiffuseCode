import numpy as np
#from libcpp cimport bool as c_bool

### config_mod ###
cdef extern int maxscat
cdef extern int nmax
# end config_mod #

### crystal_mod ###
cdef extern int cr_ncatoms
cdef extern int cr_icc[3]
cdef extern float cr_v
cdef extern float cr_gten[3][3]
cdef extern int cr_nscat
cdef extern float cr_a0[3]
cdef extern float cr_win[3]
cdef extern int cr_natoms
cdef extern int cr_n_real_atoms
# end crystal_mod #

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
cdef extern int pdf_poly_n
#cdef extern bool pdf_lxray
#cdef extern bool pdf_gauss
#cdef extern bool pdf_2d
#cdef extern bool pdf_lweights
#cdef extern bool pdf_lrho0
#cdef extern bool pdf_lexact
#cdef extern bool pdf_lrho0_rel
# end pdf_mod #

### powder_mod ###
cdef extern int POW_MAXPKT
cdef extern int pow_axis
cdef extern int pow_four_mode
cdef extern int pow_four_type
cdef extern int pow_four_vers
cdef extern int pow_lp
#cdef extern bool pow_l_all
cdef extern float pow_tthmin
cdef extern float pow_tthmax
cdef extern float pow_deltatth
cdef extern float pow_qmin
cdef extern float pow_qmax
cdef extern float pow_deltaq
cdef extern float pow_ds_max
cdef extern float pow_ds_min
cdef extern float pow_delta
cdef extern float pow_lp_fac
cdef extern float pow_lp_ang
cdef extern float pow_lp_cos
cdef extern float pow_back[6]
cdef extern float pow_scale
cdef extern float pow_hkl_max[3]
cdef extern float pow_hkl_del[3]
cdef extern float pow_hkl_shift[3]
#cdef extern bool pow_pref
cdef extern int pow_pref_type
cdef extern float pow_pref_g1
cdef extern float pow_pref_g2
cdef extern float pow_pref_hkl[3]
cdef extern int pow_profile
cdef extern int pow_pr_par
cdef extern float pow_fwhm
cdef extern float pow_eta
cdef extern float pow_etax
cdef extern float pow_u
cdef extern float pow_v
cdef extern float pow_w
cdef extern float pow_p1
cdef extern float pow_p2
cdef extern float pow_p3
cdef extern float pow_p4
cdef extern float pow_width
cdef extern int pow_size_of
# end powder_mod #

### diffuse_mod ###
cdef extern int dif_maxat
cdef extern int dif_maxscat
cdef extern float fave
cdef extern int nlots
cdef extern int ilots
cdef extern int ls_xyz[3]
cdef extern int nxat
cdef extern int inc[3]
cdef extern int diff_radiation
cdef extern float eck[4][3]
cdef extern float vi[3][3]
cdef extern float rlambda
# end diffuse_mod #


### routines ###
cdef extern void setup_c "setup" ()
cdef extern void discus_loop()
# end routines #

### crystal routines ###
cdef extern void get_cr_pos_c "get_cr_pos" (float*,int)
cdef extern void set_cr_pos_c "set_cr_pos" (float*,int)
cdef extern void get_cr_dw_c "get_cr_dw" (float*,int)
cdef extern void get_cr_scat_c "get_cr_scat" (float*,int)
cdef extern void set_cr_scat_c "set_cr_scat" (float*,int)
cdef extern void get_cr_iscat_c "get_cr_iscat" (int*,int)
cdef extern void get_crystal_name_f "get_crystal_name" (char*)
cdef extern void get_crystal_spcgr_f "get_crystal_spcgr" (char*)
cdef extern void rese_cr_f "rese_cr" ()
cdef extern void read_cell_f "read_cell" (char*,int,int,int,int)
cdef extern void read_stru_f "read_stru" (char*,int)
cdef extern void get_atom_type (char*,int)
cdef extern void set_atom_type (char*,int)
# end crystal routines #

### pdf routines ###
cdef extern void get_pdf_c "get_pdf" (double*,int,int)
cdef extern void pdf_determine_c (bool)
cdef extern void pdf_show_c(char*,int)
cdef extern void pdf_setup_c "pdf_setup" ()
cdef extern void alloc_pdf_f ()
cdef extern void set_pdf_logical_f "set_pdf_logical" (bool,bool,bool,
                                                      bool,bool,bool,
                                                      bool,bool,bool,
                                                      bool)
# end pdf routines #

### powder routines ###
cdef extern void get_powder_c "get_powder" (double*,int,int)
cdef extern void powder_run_c ()
cdef extern void pow_show_c "pow_show" ()
# end powder routines #

### fourier routines ###
cdef extern void get_diffuse_dsi_c "get_diffuse_dsi" (float*, int)
cdef extern void four_show_c (bool)
cdef extern void four_run_f ()
cdef extern void set_diffuse_logical_f "set_diffuse_logical" (bool,bool)
cdef extern void set_ltop_f "set_ltop" (bool)
cdef extern void dlink_f ()
# end fourier routines #

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

def get_cr():
    return cr_ncatoms,get_cr_icc(),cr_v, \
        get_cr_gten(),get_cr_scat(),get_crystal_name().strip(), \
        get_crystal_spcgr(),get_cr_at_lis(),cr_nscat, \
        cr_natoms,cr_n_real_atoms,get_cr_a0(), \
        get_cr_win(),get_cr_dw(),get_cr_iscat(), \
        get_cr_pos()

def set_cr(ncatoms,icc,v, \
               gten,scat,name, \
               spcgr,at_lis,nscat, \
               natoms,n_real_atoms,a0, \
               win,dw,iscat, \
               pos):
    global cr_ncatoms,cr_icc,cr_v,cr_gten,cr_nscat,cr_a0,cr_win,cr_natoms,cr_n_real_atoms
    cr_ncatoms  = ncatoms
    set_cr_icc(icc)
    cr_v        = v
    
    set_cr_gten(gten)
    set_cr_scat(scat,nscat)
    #name
    
    #spcgr
    #at_lis
    #set_cr_at_lis(at_lis)
    cr_nscat    = nscat
    
    cr_natoms   = natoms
    cr_n_real_atoms = n_real_atoms
    set_cr_a0(a0)
    
    set_cr_win(win)
    #dw
    #iscat
    
    set_cr_pos(pos,natoms)
    

def get_cr_icc():
    icc = np.empty(3,dtype=np.int32)
    cdef int [:] icc2 = icc
    icc2[:] = cr_icc
    return icc

def set_cr_icc(icc):
    global cr_icc
    for i in range(3):
        cr_icc[i] = icc[i]

def get_cr_a0():
    a0 = np.empty(3,dtype=np.float32)
    cdef float [:] a02 = a0
    a02[:] = cr_a0
    return a0

def set_cr_a0(a0):
    global cr_a0
    for i in range(3):
        cr_a0[i] = a0[i]

def get_cr_win():
    win = np.empty(3,dtype=np.float32)
    cdef float [:] win2 = win
    win2[:] = cr_win
    return win

def set_cr_win(win):
    global cr_win
    for i in range(3):
        cr_win[i] = win[i]

def get_cr_gten():
    gten = np.empty((3,3),dtype=np.float32,order='F')
    cdef float [:,:] gten2 = gten
    gten2[:][:] = cr_gten
    return gten

def set_cr_gten(gten):
    global cr_gten
    gten = np.asfortranarray(gten)
    for i in range(3):
        for j in range(3):
            cr_gten[i][j] = gten[i][j]

def get_cr_pos():
    foo = np.empty((3,nmax),dtype=np.float32,order='F')
    cdef float [::1,:] bar = foo
    get_cr_pos_c(&bar[0][0],nmax)
    return foo

def set_cr_pos(pos_in,n):
    cdef float [::1,:] cr_pos_in = np.asfortranarray(pos_in)
    set_cr_pos_c(&cr_pos_in[0][0],n)

def get_cr_dw():
    foo = np.empty((maxscat+1),dtype=np.float32,order='F')
    cdef float [:] bar = foo
    get_cr_dw_c(&bar[0],maxscat+1)
    return foo

def get_cr_scat():
    foo = np.empty((11,cr_nscat+1),dtype=np.float32,order='F')
    cdef float [::1,:] bar = foo
    get_cr_scat_c(&bar[0][0],cr_nscat+1)
    return foo

def set_cr_scat(scat,nscat):
    #foo = np.empty((11,cr_nscat+1),dtype=np.float32,order='F')
    cdef float [::1,:] cr_scat = scat
    set_cr_scat_c(&cr_scat[0][0],nscat+1)

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

def read_stru(fname):
    cdef char* c_str = fname
    cdef int fname_length = len(fname)
    read_stru_f(&c_str[0],fname_length)

def get_cr_at_lis():
    out = dict()
    cdef char[5] c_str
    for i in range(maxscat+1):
        get_atom_type(&c_str[0],i)
        out[i]=c_str.strip()
    return out

#def set_cr_at_lis(at_lis):
#    #cdef char[5] c_str
#    for i in at_lis:
#        set_atom_type(at_lis[i],i)

def set_pdf(rmax,qmax,deltar,
            skal,sigmaq,xq,
            rfmin,rfmax,delta,
            rcut,srat,gamma,
            qalp,dnorm,rho0,
            sphere,diam_poly,diam,
            shape,scale,poly,
            bin,finite,radiation,
            poly_n):
    global pdf_nscat,pdf_ndat,pdf_nbnd,pdf_rmax,pdf_qmax,pdf_deltar,pdf_skal,pdf_sigmaq,pdf_xq,pdf_rfmin,pdf_rfmax,pdf_delta,pdf_rcut,pdf_srat,pdf_gamma,pdf_qalp,pdf_dnorm,pdf_rho0,pdf_sphere,pdf_diam_poly,pdf_diam,pdf_shape,pdf_scale,pdf_poly,pdf_bin,pdf_finite,pdf_radiation,pdf_poly_n
    pdf_rmax      = rmax
    pdf_qmax      = qmax
    pdf_deltar    = deltar
    pdf_skal      = skal
    pdf_sigmaq    = sigmaq
    pdf_xq        = xq
    pdf_rfmin     = rfmin
    pdf_rfmax     = rfmax
    pdf_delta     = delta
    pdf_rcut      = rcut
    pdf_srat      = srat
    pdf_gamma     = gamma
    pdf_qalp      = qalp
    pdf_dnorm     = dnorm
    pdf_rho0      = rho0
    pdf_sphere    = sphere
    pdf_diam_poly = diam_poly
    pdf_diam      = diam
    pdf_shape     = shape
    pdf_scale     = scale
    #poly      = np.zeros(5,dtype=np.int32)
    #cdef int [:] poly2 = poly
    #pdf_poly = poly2
    pdf_bin       = bin
    pdf_finite    = finite
    pdf_radiation = radiation
    pdf_poly_n    = poly_n
    if bin>pdf_maxdat:
        pdf_nscat = max(pdf_nscat, cr_nscat, pdf_maxscat, maxscat)
        pdf_ndat  = max(pdf_ndat , pdf_bin , pdf_maxdat)
        pdf_nbnd  = max(pdf_nbnd ,           pdf_maxbnd)
        alloc_pdf_f()

def set_pdf_logical(lxray,gauss,d2d,
                    lweights,lrho0,lexact,
                    lrho0_rel,chem_period):
    chem_period1 = chem_period[0]
    chem_period2 = chem_period[1]
    chem_period3 = chem_period[2]
    set_pdf_logical_f(lxray,gauss,d2d,
                      lweights,lrho0,lexact,
                      lrho0_rel,
                      chem_period1,chem_period2,chem_period3)

def get_diffuse_dsi(i1,i2,i3):
    dsi = np.empty((i3,i2,i1),dtype=np.float32,order='F')
    cdef float [::1,:,:] dsi_c = dsi
    get_diffuse_dsi_c(&dsi_c[0][0][0],i1*i2*i3)
    return dsi

def four_show_f():
    four_show_c(True)

def set_fourier(f_fave,f_nlots,f_ilots,f_ls_xyz,f_inc,radiation,f_eck,f_vi,f_rlambda):
    global fave,nlots,ilots,ls_xyz,inc,diff_radiation,eck,vi,rlambda
    fave              = f_fave
    nlots             = f_nlots
    ilots             = f_ilots
    diff_radiation    = radiation
    rlambda           = f_rlambda
    for i in range(3):
        ls_xyz[i]     = f_ls_xyz[i]
    for i in range(3):
        inc[i]        = f_inc[i]
    for i in range(4):
        for j in range(3):
            eck[i][j] = f_eck[i][j]
    for i in range(3):
        for j in range(3):
            vi[i][j]  = f_vi[i][j]

def set_diffuse_logical(lperiod,lxray):
    set_diffuse_logical_f(lperiod,lxray)

def set_ltop(ltop):
    set_ltop_f(ltop)

def four_run():
    four_run_f()

def dlink():
    dlink_f()
    
def set_structure():
    pass

def pow_show():
    pow_show_c()

def pow_calc():
    powder_run_c()
