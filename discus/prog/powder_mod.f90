MODULE powder_mod
!+
!     variables needed for the powder diffraction
!-
USE discus_config_mod
USE precision_mod
!
SAVE
!
INTEGER             :: POW_MAXPKT          = 1
!
INTEGER, PARAMETER  :: POW_FOURIER         = 0
INTEGER, PARAMETER  :: POW_STACK           = 1
!
INTEGER, PARAMETER  :: POW_COMPL           = 0
INTEGER, PARAMETER  :: POW_DEBYE           = 1
!
INTEGER, PARAMETER  :: POW_PREF_RIET       = 1
INTEGER, PARAMETER  :: POW_PREF_MARCH      = 2
!
INTEGER, PARAMETER  :: POW_PROFILE_GAUSS   = 1
INTEGER, PARAMETER  :: POW_PROFILE_PSVGT   = 2
!
INTEGER, PARAMETER  :: POW_AXIS_DSTAR      = 0
INTEGER, PARAMETER  :: POW_AXIS_Q          = 1
INTEGER, PARAMETER  :: POW_AXIS_TTH        = 2
!
INTEGER, PARAMETER  :: POW_PROFILE_PAR_TTH = 0
INTEGER, PARAMETER  :: POW_PROFILE_PAR_Q   = 1
!
INTEGER, PARAMETER  :: POW_LP_NONE         = 0
INTEGER, PARAMETER  :: POW_LP_BRAGG        = 1
INTEGER, PARAMETER  :: POW_LP_NEUT         = 2
INTEGER, PARAMETER  :: POW_LP_SYNC         = 3
INTEGER, PARAMETER  :: POW_LP_CORRE        = 4
!
INTEGER                  :: pow_axis       = POW_AXIS_Q
INTEGER                  :: pow_npkt       = 1           ! Actual number of powder data points
!
INTEGER                  :: pow_four_mode  = 0
INTEGER                  :: pow_four_type  = POW_COMPL
!INTEGER                  :: pow_four_vers  = POW_DEBYE
!
INTEGER                  :: pow_lp         = POW_LP_BRAGG
!
LOGICAL                  :: pow_l_all      = .true.
LOGICAL                  :: pow_qtthmin    = .TRUE.    ! User  provided Qmin(==true) TTHmin(=false)
LOGICAL                  :: pow_qtthmax    = .TRUE.    ! User  provided Qmax(==true) TTHmax(=false)
LOGICAL                  :: pow_deltaqtth  = .TRUE.    ! User  provided Qstp(==true) TTHstp(=false)
LOGICAL                  :: pow_qtthzero   = .TRUE.    ! User  provided Qzero(==true) TTHzero(=false)
!
REAL(KIND=PREC_DP)       :: pow_tthzero    =  0.0D0
REAL(KIND=PREC_DP)       :: pow_tthmin     =  0.1D0
REAL(KIND=PREC_DP)       :: pow_tthmax     = 40.0D0
REAL(KIND=PREC_DP)       :: pow_deltatth   =  0.05D0
REAL(KIND=PREC_DP)       :: pow_qzero      =  0.0D0
REAL(KIND=PREC_DP)       :: pow_qmin       =  0.2D0
REAL(KIND=PREC_DP)       :: pow_qmax       =  7.0D0
REAL(KIND=PREC_DP)       :: pow_deltaq     =  0.001D0
REAL(KIND=PREC_DP)       :: pow_qmin_c     =  0.2D0    ! Limits used in actual calculation
REAL(KIND=PREC_DP)       :: pow_qmax_c     =  7.0D0    ! May vary due to FWHM, zero point
REAL(KIND=PREC_DP)       :: pow_deltaq_c   =  0.001D0
REAL(KIND=PREC_DP)       :: pow_qmin_u     =  0.2D0    ! Temporary stored user limits in case of corrlin, corrquad
REAL(KIND=PREC_DP)       :: pow_qmax_u     =  7.0D0    ! "         Actual calculation will proceed to 
REAL(KIND=PREC_DP)       :: pow_deltaq_u   =  0.0001D0 ! "         qmax * 1.5
REAL(KIND=PREC_DP)       :: pow_tthmax_buf =  1.0D0    ! additional buffer 
REAL(KIND=PREC_DP)       :: pow_qmax_buf   =  0.5D0    ! for convolutions
REAL(kind=PREC_DP)       :: pow_ds_max     =  0.0001
REAL(kind=PREC_DP)       :: pow_ds_min     =  0.0001
REAL(kind=PREC_DP)       :: pow_delta      =  0.0
REAL(kind=PREC_DP)       :: pow_lp_fac     =  0.88
REAL(kind=PREC_DP)       :: pow_lp_ang     = 20.0
REAL(kind=PREC_DP)       :: pow_lp_cos     =  0.936
!
integer                  :: pow_npkt_u     =  1
LOGICAL                  :: pow_lperiod    = .FALSE.
REAL(kind=PREC_DP)       :: pow_period     =  0.000
REAL(kind=PREC_DP)       :: pow_period_cut =  0.800
!
INTEGER                  :: pow_nback      = 0
REAL(kind=PREC_DP)   , DIMENSION(0:5)  :: pow_back       = 0.0
REAL(kind=PREC_DP)       :: pow_scale      = 1.0
!
REAL(kind=PREC_DP)   , DIMENSION(3)    :: pow_hkl_max    = 4.0
REAL(kind=PREC_DP)   , DIMENSION(3)    :: pow_hkl_del    = 0.05
REAL(kind=PREC_DP)   , DIMENSION(3)    :: pow_hkl_shift  = 0.00
!
LOGICAL                  :: pow_pref       = .false.
INTEGER                  :: pow_pref_type  = POW_PREF_RIET
REAL(kind=PREC_DP)       :: pow_pref_g1    = 0.0
REAL(kind=PREC_DP)       :: pow_pref_g2    = 0.0
REAL(kind=PREC_DP), DIMENSION(3)    :: pow_pref_hkl   = (/0.0D0, 0.0D0, 1.0D0/)
!
INTEGER                  :: pow_profile    = POW_PROFILE_PSVGT
INTEGER                  :: pow_pr_par     =  0
REAL(kind=PREC_DP)       :: pow_fwhm       =  0.01
REAL(kind=PREC_DP)       :: pow_eta        =  0.5
REAL(kind=PREC_DP)       :: pow_eta_l      =  0.0
REAL(kind=PREC_DP)       :: pow_eta_q      =  0.0
REAL(kind=PREC_DP)       :: pow_u          =  0.0
REAL(kind=PREC_DP)       :: pow_v          =  0.0
REAL(kind=PREC_DP)       :: pow_w          =  0.05
REAL(kind=PREC_DP)   , DIMENSION(4,3)  :: pow_asym       =  0.0D0
REAL(kind=PREC_DP)       :: pow_width      = 20.0D0
!
REAL(kind=PREC_DP)       :: pow_ka21       =  0.0d0
LOGICAL                  :: pow_ka21_u     =  .FALSE.
!
!
logical                               :: pow_l_partial = .FALSE.   ! Default to full powder/PDF
logical , dimension(:,:), allocatable :: pow_do_partial    ! 0:maxscat, 0:maxscat)
integer , dimension(:  ), allocatable :: pow_nn_partial    ! 0:maxscat, 0:maxscat)
!
REAL(KIND=PREC_DP), DIMENSION(:)    , ALLOCATABLE :: pow_qsp     !  (0:POW_MAXPKT)
REAL(KIND=PREC_DP), DIMENSION(:)    , ALLOCATABLE :: pow_f2aver  !  (0:POW_MAXPKT)
REAL(KIND=PREC_DP), DIMENSION(:)    , ALLOCATABLE :: pow_faver2  !  (0:POW_MAXPKT)
REAL(KIND=PREC_DP), DIMENSION(:,:)  , ALLOCATABLE :: pow_f2      !  (0:POW_MAXPKT, 0:MAXSCAT)
REAL(KIND=PREC_DP), DIMENSION(:)    , ALLOCATABLE :: pow_fu      !  (0:POW_MAXPKT, 0:MAXSCAT)
REAL(KIND=PREC_DP), DIMENSION(:)    , ALLOCATABLE :: pow_conv    !  (0:POW_MAXPKT)
REAL(KIND=PREC_DP), DIMENSION(:)    , ALLOCATABLE :: pow_sq      !  (0:POW_MAXPKT)
INTEGER                            :: pow_nreal  = 0
INTEGER                            :: pow_ncreal = 0
REAL(kind=PREC_DP)                 :: pow_u2aver = 0.0
!
INTEGER                  :: pow_size_of  = 0 ! Bytes allocated for powder
!
END MODULE powder_mod
