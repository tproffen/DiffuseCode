MODULE mmc_mod
!+
!     Variables for MONTE-CARLO multi level
!-
USE discus_config_mod
USE mc_mod
!
USE precision_mod
!
SAVE
!
LOGICAL, PARAMETER  ::  MMC_CLASSIC        = .TRUE.
LOGICAL, PARAMETER  ::  MMC_GROWTH         = .FALSE.
!
!
!INTEGER, PARAMETER  ::  MC_N_ENERGY        = 12
INTEGER, PARAMETER  ::  MC_MULTI_ENERGY    =  1
INTEGER, PARAMETER  ::  MC_SINGLE_ENERGY   =  0
!
INTEGER, PARAMETER  ::  MC_N_MOVE          =  8
!
INTEGER, PARAMETER  ::  MMC_SELECT_RANDOM  =  0
INTEGER, PARAMETER  ::  MMC_SELECT_ALL     =  0
!
INTEGER, PARAMETER  ::  MMC_MAX_ATOM       =  4    ! Number of selected atoms
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
INTEGER, PARAMETER  ::  MC_MOVE_ROTATE     =  5
INTEGER, PARAMETER  ::  MC_MOVE_SWNEIG     =  6
INTEGER, PARAMETER  ::  MC_MOVE_EXNEIG     =  7
INTEGER, PARAMETER  ::  MC_MOVE_SWVALUE    =  8
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
INTEGER             ::  MMC_H_NNNN         =  5
!
LOGICAL             ::  mmc_style          = MMC_IS_ATOM
!
!LOGICAL             ::  mmc_algo           = MMC_GROWTH         ! MMC Algorithm == CLASSIC or GROWTH
LOGICAL             ::  mmc_algo           = MMC_CLASSIC        ! MMC Algorithm == CLASSIC or GROWTH
REAL(kind=PREC_DP)  ::  mmc_g_rate         = 0.05     ! Use unlikely atoms in growth select
REAL(kind=PREC_DP)  ::  mmc_g_bad          = 0.05     ! Accept bad moved   in growth
REAL(kind=PREC_DP)  ::  mmc_g_neut         = 0.25     ! Accept neutral moved   in growth
INTEGER             ::  mmc_move           =  0
INTEGER ::  mmc_select_mode
!
logical             :: mmc_out_feed  = .true.         ! Screen output during feedbacks
logical             :: mmc_out_final = .true.         ! Screen after loop
!
INTEGER, DIMENSION(MC_N_MOVE)          ::  mmc_local      = 1   !=rmc_local_all  (MC_N_MOVE)
REAL(kind=PREC_DP)   , DIMENSION(MC_N_MOVE)          ::  mmc_move_prob  = 0.0 ! (MC_N_MOVE)
REAL(kind=PREC_DP)   , DIMENSION(MC_N_MOVE)          ::  mmc_move_cprob = 0.0 ! (MC_N_MOVE)
!
INTEGER, DIMENSION(0:MC_N_ENERGY)      ::  n_e_av_p       = 0   ! (0:MC_N_ENERGY)
INTEGER, DIMEnSION(0:MC_N_ENERGY)      ::  n_e_av_m       = 0   ! (0:MC_N_ENERGY)
INTEGER, DIMEnSION(0:MC_N_ENERGY)      ::  n_e_av_z       = 0   ! (0:MC_N_ENERGY)
REAL(kind=PREC_DP)   , DIMEnSION(0:MC_N_ENERGY)      ::  e_aver_p       = 0.0 ! (0:MC_N_ENERGY)
REAL(kind=PREC_DP)   , DIMEnSION(0:MC_N_ENERGY)      ::  e_aver_m       = 0.0 ! (0:MC_N_ENERGY)
!
INTEGER ::  mmc_n_angles        =  0
LOGICAL ::  mmc_l_constrains    = .false.
LOGICAL ::  mmc_sel_atom        = .true.
INTEGER ::  mmc_sel_prop(0:1)   =  0
INTEGER ::  mmc_constrain_type
REAL(kind=PREC_DP)    ::  mmc_c_min(3)        = -1.E10
REAL(kind=PREC_DP)    ::  mmc_c_max(3)        =  1.E10
REAL(kind=PREC_DP)    ::  mmc_c_rad
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE ::  mmc_nvec        ! (CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
integer            , dimension(:,:,:,:), allocatable ::  mmc_ach_pairs   ! Initial number of pairs
integer            , dimension(:,:,:,:), allocatable ::  mmc_ini_pairs   ! Initial number of pairs
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_target_corr ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_depth       ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE, target ::  mmc_ach_corr    ! (CHEM_MAX_COR,0:MC_N_ENERGY,-1:DEF_MAXSCAT,-1:DEF_MAXSCAT)
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE, target ::  mmc_ini_corr    ! (CHEM_MAX_COR,0:MC_N_ENERGY,-1:DEF_MAXSCAT,-1:DEF_MAXSCAT)
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_ach_sigm    ! (CHEM_MAX_COR,0:MC_N_ENERGY,-1:DEF_MAXSCAT,-1:DEF_MAXSCAT)
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_ini_sigm    ! (CHEM_MAX_COR,0:MC_N_ENERGY,-1:DEF_MAXSCAT,-1:DEF_MAXSCAT)
REAL(kind=PREC_DP) , DIMENSION(:,:)    , ALLOCATABLE ::  mmc_const       ! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
REAL(kind=PREC_DP) , DIMENSION(:,:)    , ALLOCATABLE ::  mmc_cfac        ! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
logical, DIMENSION(:,:)    , ALLOCATABLE ::  mmc_lfeed       ! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
REAL(kind=PREC_DP) , DIMENSION(:)      , ALLOCATABLE ::  mmc_depth_def   ! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:)  , ALLOCATABLE ::  mmc_pid_diff    ! PID Difference term
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:,:), ALLOCATABLE ::  mmc_pid_inte    ! PID Integral term
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:,:), ALLOCATABLE ::  mmc_pid_deri    ! PID Differential term
REAL(kind=PREC_DP) , DIMENSION(3,2)                    ::  mmc_pid_pid     ! PID PID parameters, Average of (diff, inte,dervi)
real(kind=PREC_DP) , dimension(:)        , allocatable ::  mmc_comp        ! Crystal composition
integer                                  :: mmc_pid_pid_n    ! Number of correlations contributing
real(kind=PREC_DP) :: mmc_pid_change = 0.0
REAL(kind=PREC_DP) , DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_pre_corr    ! Previous correlation 
INTEGER, DIMENSION(:,:,:), ALLOCATABLE ::  mmc_pneig       ! (0:DEF_MAXSCAT, 0:DEF_MAXSCAT, 1:CHEM_MAX_COR)
!
logical                                :: mmc_feed_auto=.FALSE.  ! Set Feedback cycles automatically
REAL(kind=PREC_DP), DIMENSION(:,:)     , ALLOCATABLE :: mmc_h_diff             ! history of achieved correlation differences
INTEGER                                :: mmc_h_number= 0        ! Number of achieved targets
INTEGER                                :: mmc_h_ctarg = 0        ! Number of current  target 
INTEGER                                :: mmc_h_index = 0        ! Current cycle entry in achieved history
INTEGER                                :: mmc_m_index = 0        ! Maximum number of feedback cycles so far
INTEGER                                :: mmc_h_ncycl = 0        ! Number of feedback cycles achieved
REAL(kind=PREC_DP), DIMENSION(:)       , ALLOCATABLE :: mmc_h_targ             ! Target values
REAL(kind=PREC_DP), DIMENSION(:)       , ALLOCATABLE :: mmc_h_aver             ! average changes from cycle to cycle
REAL(kind=PREC_DP), DIMENSION(:)       , ALLOCATABLE :: mmc_h_aver_r           ! average changes from cycle to cycle
REAL(kind=PREC_DP), DIMENSION(:,:)     , ALLOCATABLE :: mmc_h_maxd             ! Maximum change from cycle to cycle
REAL(kind=PREC_DP)                                   :: mmc_h_conv_m = 1.0E10  ! convergence Maximum difference to target over last cycles
REAL(kind=PREC_DP)                                   :: mmc_h_conv_r = 1.0E10  ! convergence Maximum difference to target over last cycles
REAL(kind=PREC_DP)                                   :: mmc_h_conv_c = 1.0E10  ! convergence Maximum change in difference over last cycles
REAL(kind=PREC_DP)                                   :: mmc_h_conv_a = 1.0E10  ! convergence average change in difference over last cycles
LOGICAL                                :: mmc_h_stop   = .TRUE.  ! stop upon cycles==F or convergence==T
LOGICAL                                :: mmc_h_log    = .FALSE. ! Screen log on / off 
INTEGER                                :: mmc_h_nfeed  = 0       ! Number of feed back this run
!
LOGICAL, DIMENSION(:,:)    , ALLOCATABLE ::  mmc_cor_energy! (0:CHEM_MAX_COR,0:MC_N_ENERGY)
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE ::  mmc_pair      ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
INTEGER, DIMENSION(:,:,:  ), ALLOCATABLE ::  mmc_left      ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT)
INTEGER, DIMENSION(:,:,:  ), ALLOCATABLE ::  mmc_right     ! (CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:),       ALLOCATABLE ::  mmc_latom     ! (0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:),       ALLOCATABLE ::  mmc_lsite     ! (0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:),       ALLOCATABLE ::  mmc_allowed   ! (0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:,:),     ALLOCATABLE ::  mmc_allowed_site   ! (0:DEF_MAXSCAT)
!
!
INTEGER, DIMENSION(:), ALLOCATABLE     ::  mmc_angles     ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE     ::  mmc_target_angl! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE     ::  mmc_depth_angl ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE     ::  mmc_ach_angl   ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE     ::  mmc_ang_sigm   ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE     ::  mmc_ini_angl   ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE     ::  mmc_ini_sang   ! (CHEM_MAX_COR*MMC_MAX_ANGLES)
!
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_a      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_b      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_m      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_len_n      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_a      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_b      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_c      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_rep_m      ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP)                                   ::  mmc_rep_low = 1e9
!
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_a     ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_rho   ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_b     ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_rmin  ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE    ::  mmc_buck_atmin ! (0:CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
LOGICAL ::  mmc_l_limited
INTEGER ::  mmc_l_type
INTEGER ::  mmc_l_center(3)
INTEGER ::  mmc_l_extend(3)
INTEGER ::  mmc_l_lower
INTEGER ::  mmc_l_upper
INTEGER(KIND=PREC_INT_LARGE) ::  mmc_no_valid = 1000
!
END MODULE mmc_mod
