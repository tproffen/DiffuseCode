MODULE pdf_mod
!+
!     This file contains variables for pdf routines
!-
USE discus_config_mod
USE precision_mod
!
SAVE
!
INTEGER, PARAMETER  ::  PDF_BACK_PERIOD  = 0
INTEGER, PARAMETER  ::  PDF_BACK_SPHERE  = 1
INTEGER, PARAMETER  ::  PDF_BACK_POLY    = 2
INTEGER, PARAMETER  ::  PDF_BACK_TANH    = 3
!
INTEGER , PARAMETER  :: PDF_RAD_XRAY = 1
INTEGER , PARAMETER  :: PDF_RAD_NEUT = 2
INTEGER , PARAMETER  :: PDF_RAD_ELEC = 3
!
INTEGER, PARAMETER :: PDF_DO_CALC = 0 ! Calculate PDF from structure
INTEGER, PARAMETER :: PDF_DO_FIT  = 1 ! Refine    PDF
INTEGER, PARAMETER :: PDF_DO_SHOW = 2 ! Show PDF, mode not yet known
!
INTEGER             ::  PDF_MAXSCAT      = 1
INTEGER             ::  PDF_MAXSITE      = 1
INTEGER             ::  PDF_MAXDAT       = 1
INTEGER             ::  PDF_MAXBND       = 1
INTEGER             ::  PDF_MAXTEMP      = 1
INTEGER             ::  PDF_MAXSINCC     = 2**12+1
!
INTEGER             ::  pdf_nscat = 1
INTEGER             ::  pdf_nsite = 1
INTEGER             ::  pdf_ndat  = 1
INTEGER             ::  pdf_nbnd  = 1
INTEGER             ::  pdf_ntemp = 1
!
REAL(PREC_DP) , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_calc   ! (MAXDAT)
REAL(PREC_DP) , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_ppp    ! (MAXDAT)
REAL(PREC_DP) , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_corr   ! (MAXDAT)
INTEGER, DIMENSION(:,:,:,:),ALLOCATABLE  ::  pdf_temp   ! (MAXTEMP,0:MAXSCAT,0:MAXSCAT)
REAL(kind=PREC_DP)   , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_obs    ! (MAXDAT)
REAL(kind=PREC_DP)   , DIMENSION(  :  ),ALLOCATABLE  ::  pdf_wic    ! (MAXDAT)
REAL(PREC_DP), DIMENSION(  :  ),ALLOCATABLE  ::  pdf_sincc  ! (2*MAXDAT)
REAL(kind=PREC_DP)   , DIMENSION(:,:  ),ALLOCATABLE  ::  pdf_weight ! (0:PDF_MAXSCAT,0:PDF_MAXSCAT)
REAL(PREC_DP), DIMENSION(  :  ),ALLOCATABLE  ::  pdf_exp    ! (4000)
!
REAL(kind=PREC_DP)                ::  pdf_rmin   =  0.001   ! Minimum distance to calculate, internal value
REAL(kind=PREC_DP)                ::  pdf_rminu  =  0.01    ! Minimum distance to calculate, user value
REAL(kind=PREC_DP)                ::  pdf_rmax   = 50.00    ! Maximum distance to calculate, internal value
REAL(kind=PREC_DP)                ::  pdf_rmaxu  = 50.00    ! Maximum distance to calculate, user value
REAL(kind=PREC_DP)                ::  pdf_qmax   = 30.00
REAL(kind=PREC_DP)                ::  pdf_deltari=  0.001   !0.001   ! internal delta r
REAL(kind=PREC_DP)                ::  pdf_deltar =  0.001   !0.001   ! internal delta r
REAL(kind=PREC_DP)                ::  pdf_deltars=  0.0005  !0.0005
REAL(kind=PREC_DP)                ::  pdf_deltaru=  0.01    ! User supplied delta r
INTEGER             ::  pdf_us_int =  1       ! Ratio user steps to internal steps
INTEGER             ::  pdf_mode   =  PDF_DO_CALC ! PDF mode, default to 'calc'
REAL(kind=PREC_DP)                ::  pdf_skal   =  1.00
REAL(kind=PREC_DP)                ::  pdf_sigmaq =  0.00
REAL(kind=PREC_DP)                ::  pdf_xq     =  0.00
REAL(kind=PREC_DP)                ::  pdf_rfmin  =  0.01    ! distance range for refinement, internal value
REAL(kind=PREC_DP)                ::  pdf_rfmax  =  0.01    ! distance range for refinement, internal value
REAL(kind=PREC_DP)                ::  pdf_rfminu =  0.01    ! distance range for refinement, user value
REAL(kind=PREC_DP)  ::  pdf_rfmaxu =  0.01    ! distance range for refinement, user value
REAL(kind=PREC_DP)                ::  pdf_rfminf =  0.01    ! distance range for refinement, file value
REAL(kind=PREC_DP)                ::  pdf_rfmaxf =  0.01    ! distance range for refinement, file value
REAL(kind=PREC_DP)  ::  pdf_cquad_a=  0.00    ! r^-2 dependet sigm correction for atoms
REAL(kind=PREC_DP)                ::  pdf_rcut   =  0.00
REAL(kind=PREC_DP)                ::  pdf_srat   =  1.00
REAL(kind=PREC_DP)  ::  pdf_clin_a =  0.00    ! r^-1 dependent sigma correction for atoms
REAL(kind=PREC_DP)                ::  pdf_qalp   =  0.00
REAL(kind=PREC_DP)                ::  pdf_dnorm  =  1.00
REAL(kind=PREC_DP)                ::  pdf_rho0   =  0.00
REAL(kind=PREC_DP)                ::  pdf_sphere =  0.00
REAL(kind=PREC_DP)                ::  pdf_diam_poly =  0.00
REAL(kind=PREC_DP)                ::  pdf_diam   =  0.00
REAL(kind=PREC_DP)                ::  pdf_shape  =  0.00
REAL(kind=PREC_DP)                ::  pdf_scale  =  1.00
REAL(kind=PREC_DP)                ::  pdf_poly(5)=  0.00
!
INTEGER, DIMENSION(:,:), ALLOCATABLE  :: pdf_bnd     ! (3,-MAXBND:2*MAXBND)
INTEGER             ::  pdf_bin    = 1
INTEGER             ::  pdf_finite = PDF_BACK_PERIOD
INTEGER             ::  pdf_poly_n = 0
INTEGER             ::  pdf_sel_prop(0:1) = 0
!                                                
INTEGER             ::  pdf_radiation = PDF_RAD_XRAY
INTEGER             ::  pdf_power     = 4
INTEGER             ::  pdf_nmol      = 0   ! pdf_temp dimension if molecules are relevant
LOGICAL             ::  pdf_success= .FALSE.
LOGICAL             ::  pdf_lxray  = .false.
LOGICAL             ::  pdf_gauss  = .false.
LOGICAL             ::  pdf_gauss_init  = .true.
REAL(PREC_DP), PARAMETER ::  pdf_gauss_step = 0.0005d0
LOGICAL             ::  pdf_2d     = .false.
LOGICAL, DIMENSION(:),ALLOCATABLE  ::  pdf_allowed_i ! (0:PDF_MAXSCAT)
LOGICAL, DIMENSION(:),ALLOCATABLE  ::  pdf_allowed_j ! (0:PDF_MAXSCAT)
LOGICAL, DIMENSION(:),ALLOCATABLE  ::  pdf_lsite_i   ! (0:PDF_MAXSCAT)
LOGICAL, DIMENSION(:),ALLOCATABLE  ::  pdf_lsite_j   ! (0:PDF_MAXSCAT)
INTEGER, DIMENSION(:),ALLOCATABLE  ::  pdf_has_atom  ! (0:PDF_MAXSCAT)
INTEGER, DIMENSION(:,:),ALLOCATABLE  ::  pdf_look_mol ! (0:PDF_MAXSCAT)
REAL(kind=PREC_DP),    DIMENSION(:)  , ALLOCATABLE ::  pdf_bvalue_mole ! effective mol bvalues
REAL(kind=PREC_DP),    DIMENSION(:)  , ALLOCATABLE ::  pdf_clin_mole ! linear correction mol
REAL(kind=PREC_DP),    DIMENSION(:)  , ALLOCATABLE ::  pdf_cqua_mole ! quadratic correction mol
!LOGICAL, DIMENSION(0:    MAXSCAT)  ::  pdf_allowed_i ! (0:PDF_MAXSCAT)
!LOGICAL, DIMENSION(0:    MAXSCAT)  ::  pdf_allowed_j ! (0:PDF_MAXSCAT)
LOGICAL             ::  pdf_ldata     = .false.
LOGICAL             ::  pdf_lweights  = .false.
LOGICAL             ::  pdf_lrho0     = .true.
LOGICAL             ::  pdf_lexact    = .false.
LOGICAL             ::  pdf_lrho0_rel = .false.
INTEGER             ::  pdf_size_of   = 0
!
LOGICAL               :: pdf_refine_scale   = .false.
LOGICAL               :: pdf_refine_density = .false.
LOGICAL, DIMENSION(6) :: pdf_refine_lattice = .false.
!
END MODULE pdf_mod
