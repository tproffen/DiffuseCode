MODULE pdf_mod
USE iso_c_binding
!+
!     This file contains variables for pdf routines
!-
USE config_mod
!
SAVE
!
INTEGER, PARAMETER  ::  PDF_BACK_PERIOD  = 0
INTEGER, PARAMETER  ::  PDF_BACK_POLY    = 1
INTEGER, PARAMETER  ::  PDF_BACK_SPHERE  = 2
INTEGER, PARAMETER  ::  PDF_BACK_TANH    = 3
!
INTEGER , PARAMETER  :: PDF_RAD_XRAY = 1
INTEGER , PARAMETER  :: PDF_RAD_NEUT = 2
INTEGER , PARAMETER  :: PDF_RAD_ELEC = 3
!
INTEGER(C_INT), BIND(C) ::  PDF_MAXSCAT      = 1
INTEGER(C_INT), BIND(C) ::  PDF_MAXDAT       = 1
INTEGER(C_INT), BIND(C) ::  PDF_MAXBND       = 1
!
INTEGER(C_INT), BIND(C) ::  pdf_nscat = 1
INTEGER(C_INT), BIND(C) ::  pdf_ndat  = 1
INTEGER(C_INT), BIND(C) ::  pdf_nbnd  = 1
!
REAL*8 , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_calc   ! (MAXDAT)
REAL*8 , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_corr   ! (MAXDAT)
INTEGER, DIMENSION(:,:,:),ALLOCATABLE  ::  pdf_temp   ! (MAXDAT,0:MAXSCAT,0:MAXSCAT)
REAL   , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_obs    ! (MAXDAT)
REAL   , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_wic    ! (MAXDAT)
REAL   , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_sinc   ! (2*MAXDAT)
REAL   , DIMENSION(:,:  ),ALLOCATABLE  ::  pdf_weight ! (0:PDF_MAXSCAT,0:PDF_MAXSCAT)
!
!REAL*8 , DIMENSION(MAXDAT)             ::  pdf_calc   ! (MAXDAT)
!REAL*8 , DIMENSION(MAXDAT)             ::  pdf_corr   ! (MAXDAT)
!INTEGER, DIMENSION(MAXDAT,0:MAXSCAT,0:MAXSCAT) ::  pdf_temp   ! (MAXDAT,0:MAXSCAT,0:MAXSCAT)
!REAL   , DIMENSION(MAXDAT)             ::  pdf_obs    ! (MAXDAT)
!REAL   , DIMENSION(MAXDAT)             ::  pdf_wic    ! (MAXDAT)
!REAL   , DIMENSION(2*MAXDAT)           ::  pdf_sinc   ! (2*MAXDAT)
!REAL   , DIMENSION(0:MAXSCAT,0:MAXSCAT) ::  pdf_weight ! (0:PDF_MAXSCAT,0:PDF_MAXSCAT)
REAL(C_FLOAT), BIND(C) ::  pdf_rmax   = 50.00
REAL(C_FLOAT), BIND(C) ::  pdf_qmax   = 30.00
REAL(C_FLOAT), BIND(C) ::  pdf_deltar =  0.01
REAL(C_FLOAT), BIND(C) ::  pdf_skal   =  1.00
REAL(C_FLOAT), BIND(C) ::  pdf_sigmaq =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_xq     =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_rfmin  =  0.05
REAL(C_FLOAT), BIND(C) ::  pdf_rfmax  = 15.00
REAL(C_FLOAT), BIND(C) ::  pdf_delta  =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_rcut   =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_srat   =  1.00
REAL(C_FLOAT), BIND(C) ::  pdf_gamma  =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_qalp   =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_dnorm  =  1.00
REAL(C_FLOAT), BIND(C) ::  pdf_rho0   =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_sphere =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_diam_poly =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_diam   =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_shape  =  0.00
REAL(C_FLOAT), BIND(C) ::  pdf_scale  =  1.00
REAL(C_FLOAT), BIND(C) ::  pdf_poly(5)=  0.00
!
INTEGER, DIMENSION(:,:), ALLOCATABLE  :: pdf_bnd     ! (3,-MAXBND:2*MAXBND)
!INTEGER, DIMENSION(3,-MAXBND:2*MAXBND):: pdf_bnd     ! (3,-MAXBND:2*MAXBND)
INTEGER(C_INT), BIND(C) ::  pdf_bin    = 1
INTEGER(C_INT), BIND(C) ::  pdf_finite = PDF_BACK_PERIOD
INTEGER(C_INT), BIND(C) ::  pdf_poly_n = 0
INTEGER             ::  pdf_sel_prop(0:1) = 0
!                                                
INTEGER(C_INT), BIND(C) ::  pdf_radiation = PDF_RAD_XRAY
INTEGER             ::  pdf_power     = 4
LOGICAL             ::  pdf_lxray  = .false.
LOGICAL             ::  pdf_gauss  = .false.
LOGICAL             ::  pdf_2d     = .false.
LOGICAL, DIMENSION(:),ALLOCATABLE  ::  pdf_allowed_i ! (0:PDF_MAXSCAT)
LOGICAL, DIMENSION(:),ALLOCATABLE  ::  pdf_allowed_j ! (0:PDF_MAXSCAT)
!LOGICAL, DIMENSION(0:    MAXSCAT)  ::  pdf_allowed_i ! (0:PDF_MAXSCAT)
!LOGICAL, DIMENSION(0:    MAXSCAT)  ::  pdf_allowed_j ! (0:PDF_MAXSCAT)
LOGICAL             ::  pdf_ldata     = .false.
LOGICAL             ::  pdf_lweights  = .false.
LOGICAL             ::  pdf_lrho0     = .true.
LOGICAL             ::  pdf_lexact    = .false.
LOGICAL             ::  pdf_lrho0_rel = .false.
INTEGER             ::  pdf_size_of   = 0
!
END MODULE pdf_mod
