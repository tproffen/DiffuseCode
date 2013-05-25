MODULE pdf_mod
!+
!     This file contains COMMON variables for pdf routines
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
INTEGER             ::  PDF_MAXSCAT      = 1
INTEGER             ::  PDF_MAXDAT       = 1
INTEGER             ::  PDF_MAXBND       = 1
!
INTEGER             ::  pdf_nscat = 1
INTEGER             ::  pdf_ndat  = 1
INTEGER             ::  pdf_nbnd  = 1
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
REAL                ::  pdf_rmax   = 50.00
REAL                ::  pdf_qmax   = 30.00
REAL                ::  pdf_deltar =  0.01
REAL                ::  pdf_skal   =  1.00
REAL                ::  pdf_sigmaq =  0.00
REAL                ::  pdf_xq     =  0.00
REAL                ::  pdf_rfmin  =  0.05
REAL                ::  pdf_rfmax  = 15.00
REAL                ::  pdf_delta  =  0.00
REAL                ::  pdf_rcut   =  0.00
REAL                ::  pdf_srat   =  1.00
REAL                ::  pdf_gamma  =  0.00
REAL                ::  pdf_qalp   =  0.00
REAL                ::  pdf_dnorm  =  1.00
REAL                ::  pdf_rho0   =  0.00
REAL                ::  pdf_sphere =  0.00
REAL                ::  pdf_diam_poly =  0.00
REAL                ::  pdf_diam   =  0.00
REAL                ::  pdf_shape  =  0.00
REAL                ::  pdf_scale  =  1.00
REAL                ::  pdf_poly(5)=  0.00
!
INTEGER, DIMENSION(:,:), ALLOCATABLE  :: pdf_bnd     ! (3,-MAXBND:2*MAXBND)
!INTEGER, DIMENSION(3,-MAXBND:2*MAXBND):: pdf_bnd     ! (3,-MAXBND:2*MAXBND)
INTEGER             ::  pdf_bin    = 1
INTEGER             ::  pdf_finite = PDF_BACK_PERIOD
INTEGER             ::  pdf_poly_n = 0
INTEGER             ::  pdf_sel_prop(0:1) = 0
!                                                
INTEGER             ::  pdf_radiation = PDF_RAD_XRAY
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
!     COMMON /pdfbl/ pdf_calc,pdf_rmax,pdf_lxray,pdf_qmax,              &
!    &               pdf_deltar,pdf_bin,pdf_sinc,pdf_obs,               &
!    &               pdf_wic,pdf_skal,pdf_sigmaq,pdf_rfmin,             &
!    &               pdf_rfmax,pdf_weight,pdf_xq,pdf_bnd,               &
!    &               pdf_allowed_i,pdf_allowed_j,pdf_gauss,             &
!    &               pdf_delta,pdf_rcut,pdf_srat,pdf_dnorm,             &
!    &               pdf_sphere,                                        &
!    &               pdf_corr,pdf_2d,pdf_rho0,pdf_gamma,                &
!    &               pdf_diam,pdf_shape,pdf_ldata,pdf_qalp,             &
!    &               pdf_scale,pdf_finite,pdf_poly,pdf_poly_n,          &
!    &               pdf_lweights,pdf_lrho0_rel,pdf_temp,               &
!    &               pdf_lrho0,pdf_lexact,                              &
!    &               pdf_sel_prop
!
END MODULE pdf_mod
