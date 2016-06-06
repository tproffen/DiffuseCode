MODULE powder_mod
!+
!     variables needed for the powder diffraction
!-
USE discus_config_mod
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
INTEGER, PARAMETER  :: POW_LONG            = 2
INTEGER, PARAMETER  :: POW_FAST            = 3
INTEGER, PARAMETER  :: POW_HIST            = 4
INTEGER, PARAMETER  :: POW_NEW             = 5
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
INTEGER                  :: pow_four_vers  = POW_HIST
!
INTEGER                  :: pow_lp         = POW_LP_BRAGG
!
LOGICAL                  :: pow_l_all      = .true.
!
REAL                     :: pow_tthmin     =  0.1
REAL                     :: pow_tthmax     = 40.0
REAL                     :: pow_deltatth   =  0.05
REAL                     :: pow_qmin       =  0.2
REAL                     :: pow_qmax       =  7.0
REAL                     :: pow_deltaq     =  0.0001
REAL                     :: pow_ds_max     =  0.0001
REAL                     :: pow_ds_min     =  0.0001
REAL                     :: pow_delta      =  0.0
REAL                     :: pow_lp_fac     =  0.88
REAL                     :: pow_lp_ang     = 20.0
REAL                     :: pow_lp_cos     =  0.936
!
INTEGER                  :: pow_nback      = 0
REAL   , DIMENSION(0:5)  :: pow_back       = 0.0
REAL                     :: pow_scale      = 1.0
!
REAL   , DIMENSION(3)    :: pow_hkl_max    = 4.0
REAL   , DIMENSION(3)    :: pow_hkl_del    = 0.05
REAL   , DIMENSION(3)    :: pow_hkl_shift  = 0.00
!
LOGICAL                  :: pow_pref       = .false.
INTEGER                  :: pow_pref_type  = POW_PREF_RIET
REAL                     :: pow_pref_g1    = 0.0
REAL                     :: pow_pref_g2    = 0.0
REAL   , DIMENSION(3)    :: pow_pref_hkl   = (/0., 0., 1./)
!
INTEGER                  :: pow_profile    = POW_PROFILE_PSVGT
INTEGER                  :: pow_pr_par     =  0
REAL                     :: pow_fwhm       =  0.01
REAL                     :: pow_eta        =  0.5
REAL                     :: pow_etax       =  0.0
REAL                     :: pow_u          =  0.0
REAL                     :: pow_v          =  0.0
REAL                     :: pow_w          =  0.05
REAL                     :: pow_p1         =  0.0
REAL                     :: pow_p2         =  0.0
REAL                     :: pow_p3         =  0.0
REAL                     :: pow_p4         =  0.0
REAL                     :: pow_width      = 20.0
!
!
REAL   (KIND=KIND(0.0D0)), DIMENSION(:), ALLOCATABLE :: pow_qsp     !  (0:POW_MAXPKT)
REAL   (KIND=KIND(0.0D0)), DIMENSION(:), ALLOCATABLE :: pow_f2aver  !  (0:POW_MAXPKT)
REAL   (KIND=KIND(0.0D0)), DIMENSION(:), ALLOCATABLE :: pow_faver2  !  (0:POW_MAXPKT)
INTEGER                            :: pow_nreal
REAL                               :: pow_u2aver
!
INTEGER                  :: pow_size_of  = 0 ! Bytes allocated for powder
!
END MODULE powder_mod
