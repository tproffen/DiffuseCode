MODULE rmc_mod
!+
!     Variables for REVERSE-MONTE-CARLO level
!-
USE config_mod
!
SAVE
!
INTEGER                 :: RMC_MAX_PLANES    = 1
INTEGER                 :: RMC_MAXSCAT       = 1
INTEGER                 :: RMC_MAX_Q         = 1
INTEGER                 :: RMC_MAX_SQ        = 1
INTEGER                 :: RMC_MAX_LOTS      = 1
INTEGER                 :: RMC_MAX_SYM       = 48  ! Maximum no of sym op in rec space
!
INTEGER, PARAMETER      :: rmc_mode_shift    = 1
INTEGER, PARAMETER      :: rmc_mode_swchem   = 2
INTEGER, PARAMETER      :: rmc_mode_swdisp   = 3
INTEGER, PARAMETER      :: rmc_mode_external = 4
INTEGER, PARAMETER      :: RMC_N_MOVE        = 4
!
INTEGER, PARAMETER      :: rmc_data_nipl = 1
INTEGER, PARAMETER      :: rmc_data_pgm  = 2
!
INTEGER, PARAMETER      :: rmc_wic_eins = 1
INTEGER, PARAMETER      :: rmc_wic_sqrt = 2
INTEGER, PARAMETER      :: rmc_wic_log  = 3
INTEGER, PARAMETER      :: rmc_wic_lin  = 4
INTEGER, PARAMETER      :: rmc_wic_qua  = 5
INTEGER, PARAMETER      :: rmc_wic_inv  = 6
INTEGER, PARAMETER      :: rmc_wic_isq  = 7
INTEGER, PARAMETER      :: rmc_wic_dat  = 8
!
INTEGER, PARAMETER      :: rmc_local_all     = 1
INTEGER, PARAMETER      :: rmc_local_loc     = 2
INTEGER, PARAMETER      :: rmc_local_locsite = 3
INTEGER, PARAMETER      :: rmc_local_site    = 4
INTEGER, PARAMETER      :: rmc_local_conn    = 5
!
INTEGER , PARAMETER     :: RMC_RAD_XRAY = 1
INTEGER , PARAMETER     :: RMC_RAD_NEUT = 2
INTEGER , PARAMETER     :: RMC_RAD_ELEC = 3
!
INTEGER, DIMENSION(RMC_N_MOVE)          ::  rmc_move_local = 1   !=rmc_local_all  (RMC_N_MOVE)
REAL   , DIMENSION(RMC_N_MOVE)          ::  rmc_move_prob  = 0.0 ! (RMC_N_MOVE)
REAL   , DIMENSION(RMC_N_MOVE)          ::  rmc_move_cprob = 0.0 ! (RMC_N_MOVE)
!
LOGICAL, DIMENSION(:), ALLOCATABLE                 :: rmc_allowed   ! (0:RMC_MAXSCAT)
!
CHARACTER (LEN=80), DIMENSION(:), ALLOCATABLE     :: rmc_fname      ! (RMC_MAX_PLANES)
CHARACTER (LEN= 4), DIMENSION(:), ALLOCATABLE     :: rmc_lambda     ! (RMC_MAX_PLANES)
REAL              , DIMENSION(:,:), ALLOCATABLE   :: rmc_xy         ! (4,RMC_MAX_PLANES)
REAL   , DIMENSION(:), ALLOCATABLE                :: rmc_rlambda    ! (RMC_MAX_PLANES)
REAL   , DIMENSION(:), ALLOCATABLE                :: rmc_skal       ! (RMC_MAX_PLANES)
REAL   , DIMENSION(:), ALLOCATABLE                :: rmc_back       ! (RMC_MAX_PLANES)
REAL   , DIMENSION(:), ALLOCATABLE                :: rmc_chi2       ! (RMC_MAX_PLANES)
REAL   , DIMENSION(:), ALLOCATABLE                :: rmc_wtot       ! (RMC_MAX_PLANES)
INTEGER, DIMENSION(:), ALLOCATABLE                :: offq           ! (RMC_MAX_PLANES+1)
INTEGER, DIMENSION(:), ALLOCATABLE                :: rmc_wic_typ    ! (RMC_MAX_PLANES)
INTEGER, DIMENSION(:), ALLOCATABLE                :: rmc_nsym       ! (RMC_MAX_PLANES)
INTEGER, DIMENSION(:), ALLOCATABLE                :: rmc_constrain  ! (RMC_MAX_PLANES)
INTEGER, DIMENSION(:,:), ALLOCATABLE              :: rmc_num        ! (2,RMC_MAX_PLANES)
INTEGER, DIMENSION(:), ALLOCATABLE                :: rmc_radiation  ! (RMC_MAX_PLANES)
INTEGER, DIMENSION(:), ALLOCATABLE                :: rmc_power      ! (RMC_MAX_PLANES)
LOGICAL, DIMENSION(:), ALLOCATABLE                :: rmc_lxray      ! (RMC_MAX_PLANES)
LOGICAL, DIMENSION(:), ALLOCATABLE                :: rmc_ano        ! (RMC_MAX_PLANES)
LOGICAL, DIMENSION(:), ALLOCATABLE                :: rmc_ldbw       ! (RMC_MAX_PLANES)
!
REAL   , DIMENSION(:,:,:,:), ALLOCATABLE          :: rmc_eck        ! (3,3,RMC_MAX_SYM,RMC_MAX_PLANES)
REAL   , DIMENSION(:,:,:,:), ALLOCATABLE          :: rmc_vi         ! (3,2,RMC_MAX_SYM,RMC_MAX_PLANES)
INTEGER, DIMENSION(:,:)    , ALLOCATABLE          :: offsq          ! (RMC_MAX_PLANES+1,RMC_MAX_SYM) 

COMPLEX, DIMENSION(:,:,:)  , ALLOCATABLE          :: rcfact         ! (0:CFPKT, DEF_MAXSCAT, RMC_MAX_PLANES)
REAL   , DIMENSION(:,:)    , ALLOCATABLE          :: rmc_maxmove    ! (3,0:DEF_MAXSCAT)
REAL   , DIMENSION(:,:)    , ALLOCATABLE          :: rmc_mindist    ! (DEF_MAXSCAT,DEF_MAXSCAT)
!
COMPLEX, DIMENSION(:,:)    , ALLOCATABLE           :: rmc_csf       ! (RMC_MAX_SQ, RMC_MAX_LOTS)
COMPLEX, DIMENSION(:,:)    , ALLOCATABLE           :: rmc_csf_new   ! (RMC_MAX_SQ, RMC_MAX_LOTS)
INTEGER, DIMENSION(:)      , ALLOCATABLE           :: ristl         ! (RMC_MAX_SQ)
INTEGER, DIMENSION(:,:)    , ALLOCATABLE           :: rmc_lots_orig ! (3,RMC_MAX_LOTS)
!
REAL   , DIMENSION(:)      , ALLOCATABLE           :: rmc_int       ! (RMC_MAX_Q)
REAL   , DIMENSION(:)      , ALLOCATABLE           :: rmc_wic       ! (RMC_MAX_Q)
!
!
CHARACTER (LEN=80)                                 :: rmc_lname
!
REAL                                               :: rmc_mindist_max
REAL                                               :: rmc_ave
REAL                                               :: rmc_sigma
REAL                                               :: rmc_qmin,rmc_qmax
REAL                                               :: rmc_llim,rmc_ulim
!
INTEGER                                            :: rmc_nplane
INTEGER                                            :: rmc_data
INTEGER, DIMENSION(3)                              :: rmc_csize
INTEGER                                            :: rmc_nlots,rmc_ilots
INTEGER                                            :: rmc_maxcyc,rmc_display
INTEGER                                            :: rmc_mode,rmc_local
INTEGER, DIMENSION(0:1)                            :: rmc_sel_prop
!
LOGICAL                                            :: rmc_doskal,rmc_doback
LOGICAL                                            :: rmc_calc_f,rmc_log,rmc_dosym,rmc_nosym
LOGICAL                                            :: rmc_ranloc,rmc_sel_atom
!
INTEGER                                            :: rmc_size_of  ! Bytes allocates for rmc
INTEGER                                            :: rmc_n_sym    ! Actual number of symmetry operations
INTEGER                                            :: rmc_n_qxy    ! Size of RMC_MAX_Q arrays
INTEGER                                            :: rmc_n_sq     ! Size of RMC_MAX_Q arrays* Number of Symmetry
!     
END MODULE rmc_mod
