MODULE mmc_mod
!+
!     Variables for MONTE-CARLO multi level
!-
USE discus_config_mod
USE precision_mod
!
SAVE
!
INTEGER, PARAMETER  ::  MC_N_ENERGY        = 11
INTEGER, PARAMETER  ::  MC_MULTI_ENERGY    =  1
INTEGER, PARAMETER  ::  MC_SINGLE_ENERGY   =  0
!
INTEGER, PARAMETER  ::  MC_N_MOVE          =  4
!
INTEGER, PARAMETER  ::  MMC_SELECT_RANDOM  =  0
INTEGER, PARAMETER  ::  MMC_SELECT_ALL     =  0
!
INTEGER, PARAMETER  ::  MMC_MAX_ATOM       =  2    ! Number of selected atoms
!
INTEGER             ::  MMC_MAX_ANGLES     =  1
!
INTEGER             ::  MMC_MAX_CENT       =  2
!
INTEGER, PARAMETER  ::  MC_MOVE_NONE       =  0
INTEGER, PARAMETER  ::  MC_MOVE_SWCHEM     =  1
INTEGER, PARAMETER  ::  MC_MOVE_SWDISP     =  2
INTEGER, PARAMETER  ::  MC_MOVE_DISP       =  3
INTEGER, PARAMETER  ::  MC_MOVE_INVDISP    =  4
!
INTEGER, PARAMETER  ::  MMC_C_XYZ    =  0
INTEGER, PARAMETER  ::  MMC_C_RADIUS =  1
!
INTEGER, PARAMETER  ::  MMC_L_CELLS  =  0
INTEGER, PARAMETER  ::  MMC_L_ATOMS  =  1
!
INTEGER             ::  MMC_MAX_CORR       =  1
INTEGER             ::  MMC_MAX_SCAT       =  1
INTEGER             ::  MMC_MAX_SITE       =  1
INTEGER             ::  MMC_LENN_CORR      =  0
INTEGER             ::  MMC_LENN_SCAT      =  0
INTEGER             ::  MMC_BUCK_CORR      =  0
INTEGER             ::  MMC_BUCK_SCAT      =  0
INTEGER             ::  MMC_REP_CORR       =  0
INTEGER             ::  MMC_REP_SCAT       =  0
!
INTEGER             ::  MMC_H_NNNN         =  3
!
INTEGER             ::  mmc_move           =  0
INTEGER ::  mmc_select_mode
!
INTEGER, DIMENSION(MC_N_MOVE)          ::  mmc_local      = 1   !=rmc_local_all  (MC_N_MOVE)
REAL   , DIMENSION(MC_N_MOVE)          ::  mmc_move_prob  = 0.0 ! (MC_N_MOVE)
REAL   , DIMENSION(MC_N_MOVE)          ::  mmc_move_cprob = 0.0 ! (MC_N_MOVE)
!
INTEGER, DIMENSION(0:MC_N_ENERGY)      ::  n_e_av_p       = 0   ! (0:MC_N_ENERGY)
INTEGER, DIMEnSION(0:MC_N_ENERGY)      ::  n_e_av_m       = 0   ! (0:MC_N_ENERGY)
INTEGER, DIMEnSION(0:MC_N_ENERGY)      ::  n_e_av_z       = 0   ! (0:MC_N_ENERGY)
REAL   , DIMEnSION(0:MC_N_ENERGY)      ::  e_aver_p       = 0.0 ! (0:MC_N_ENERGY)
REAL   , DIMEnSION(0:MC_N_ENERGY)      ::  e_aver_m       = 0.0 ! (0:MC_N_ENERGY)
!
INTEGER ::  mmc_n_angles        =  0
LOGICAL ::  mmc_l_constrains    = .false.
LOGICAL ::  mmc_sel_atom        = .true.
INTEGER ::  mmc_sel_prop(0:1)   =  0
INTEGER ::  mmc_constrain_type
REAL    ::  mmc_c_min(3)        = -1.E10
REAL    ::  mmc_c_max(3)        =  1.E10
REAL    ::  mmc_c_rad
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE ::  mmc_nvec        ! (CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
REAL , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_target_corr ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_depth       ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_ach_corr    ! (CHEM_MAX_COR,0:MC_N_ENERGY,-1:DEF_MAXSCAT,-1:DEF_MAXSCAT)
REAL , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_ach_sigm    ! (CHEM_MAX_COR,0:MC_N_ENERGY,-1:DEF_MAXSCAT,-1:DEF_MAXSCAT)
!REAL,DIMENSION(:,:,:,:,:), ALLOCATABLE ::  mmc_vec         ! (4,12,CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL , DIMENSION(:,:)    , ALLOCATABLE ::  mmc_const       ! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
REAL , DIMENSION(:,:)    , ALLOCATABLE ::  mmc_cfac        ! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
!
REAL, DIMENSION(:,:)     , ALLOCATABLE :: mmc_h_diff             ! history of achieved correlation differences
INTEGER                                :: mmc_h_number= 0        ! Number of achieved targets
INTEGER                                :: mmc_h_ctarg = 0        ! Number of current  target 
INTEGER                                :: mmc_h_index = 0        ! Current cycle entry in achieved history
INTEGER                                :: mmc_h_ncycl = 0        ! Number of feedback cycles achieved
REAL, DIMENSION(:)       , ALLOCATABLE :: mmc_h_targ             ! Target values
REAL, DIMENSION(:)       , ALLOCATABLE :: mmc_h_aver             ! average changes from cycle to cycle
REAL, DIMENSION(:,:)     , ALLOCATABLE :: mmc_h_maxd             ! Maximum change from cycle to cycle
REAL                                   :: mmc_h_conv_m = 0.090   ! convergence Maximum difference to target over last cycles
REAL                                   :: mmc_h_conv_c = 0.050   ! convergence Maximum change in difference over last cycles
REAL                                   :: mmc_h_conv_a = 0.001   ! convergence average change in difference over last cycles
LOGICAL                                :: mmc_h_stop   = .TRUE.  ! stop upon cycles==F or convergence==T
!
LOGICAL, DIMENSION(:,:)    , ALLOCATABLE ::  mmc_cor_energy! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_pair      ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:),       ALLOCATABLE ::  mmc_latom     ! (0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:),       ALLOCATABLE ::  mmc_lsite     ! (0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:),       ALLOCATABLE ::  mmc_allowed   ! (0:DEF_MAXSCAT)
!
!
INTEGER, DIMENSION(:), ALLOCATABLE     ::  mmc_angles     ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL   , DIMENSION(:), ALLOCATABLE     ::  mmc_target_angl! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL   , DIMENSION(:), ALLOCATABLE     ::  mmc_depth_angl ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL   , DIMENSION(:), ALLOCATABLE     ::  mmc_ach_angl   ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL   , DIMENSION(:), ALLOCATABLE     ::  mmc_ang_sigm   ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
!
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_a      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_b      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_m      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_n      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_a      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_b      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_c      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_m      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL                                   ::  mmc_rep_low = 1e9
!
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_a     ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_rho   ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_b     ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_rmin  ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_atmin ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
LOGICAL ::  mmc_l_limited
INTEGER ::  mmc_l_type
INTEGER ::  mmc_l_center(3)
INTEGER ::  mmc_l_extend(3)
INTEGER ::  mmc_l_lower
INTEGER ::  mmc_l_upper
INTEGER(KIND=PREC_INT_LARGE) ::  mmc_no_valid = 1000
INTEGER ::  mmc_size_of = 0
!
END MODULE mmc_mod
