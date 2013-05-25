MODULE mmc_mod
!+
!     Variables for MONTE-CARLO multi level
!-
USE config_mod
!
SAVE
!
INTEGER, PARAMETER  ::  MC_N_ENERGY        =  9
INTEGER, PARAMETER  ::  MC_MULTI_ENERGY    =  1
INTEGER, PARAMETER  ::  MC_SINGLE_ENERGY   =  0
!
INTEGER, PARAMETER  ::  MC_N_MOVE          =  4
!
INTEGER, PARAMETER  ::  MMC_SELECT_RANDOM  =  0
INTEGER, PARAMETER  ::  MMC_SELECT_ALL     =  0
!
INTEGER             ::  MMC_MAX_ANGLES     =  1
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
INTEGER             ::  MMC_LENN_CORR      =  0
INTEGER             ::  MMC_LENN_SCAT      =  0
INTEGER             ::  MMC_BUCK_CORR      =  0
INTEGER             ::  MMC_BUCK_SCAT      =  0
INTEGER             ::  MMC_REP_CORR       =  0
INTEGER             ::  MMC_REP_SCAT       =  0
!
INTEGER ::  mmc_move
INTEGER ::  mmc_select_mode
!
INTEGER, DIMENSION(MC_N_MOVE)          ::  mmc_local     ! (MC_N_MOVE)
REAL   , DIMENSION(MC_N_MOVE)          ::  mmc_move_prob ! (MC_N_MOVE)
REAL   , DIMENSION(MC_N_MOVE)          ::  mmc_move_cprob! (MC_N_MOVE)
!
INTEGER, DIMENSION(0:MC_N_ENERGY)      ::  n_e_av_p      ! (0:MC_N_ENERGY)
INTEGER, DIMEnSION(0:MC_N_ENERGY)      ::  n_e_av_m      ! (0:MC_N_ENERGY)
INTEGER, DIMEnSION(0:MC_N_ENERGY)      ::  n_e_av_z      ! (0:MC_N_ENERGY)
REAL   , DIMEnSION(0:MC_N_ENERGY)      ::  e_aver_p      ! (0:MC_N_ENERGY)
REAL   , DIMEnSION(0:MC_N_ENERGY)      ::  e_aver_m      ! (0:MC_N_ENERGY)
!
INTEGER ::  mmc_n_angles
LOGICAL ::  mmc_l_constrains
LOGICAL ::  mmc_sel_atom
INTEGER ::  mmc_sel_prop(0:1)
INTEGER ::  mmc_constrain_type
REAL    ::  mmc_c_min(3)
REAL    ::  mmc_c_max(3)
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
LOGICAL ::  mmc_cor_energy(0:CHEM_MAX_COR,0:MC_N_ENERGY)
LOGICAL ::  mmc_pair(CHEM_MAX_COR,0:MC_N_ENERGY,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE     ::  mmc_latom       ! (0:DEF_MAXSCAT)
LOGICAL ::  mmc_allowed(0:DEF_MAXSCAT)
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
INTEGER ::  mmc_size_of = 0
!
!     COMMON /mmcbl/ mmc_angles,mmc_move,mmc_local,mmc_cor_energy,      &
!    &               mmc_const,mmc_allowed,                             &
!    &               mmc_latom,mmc_sel_prop,                            &
!    &               mmc_md_s,mmc_md_nr,                                &
!    &               mmc_cfac,mmc_ach_corr,mmc_ach_sigm,                &
!    &               mmc_target_corr,                                   &
!    &               mmc_pair,mmc_depth,mmc_depth_angl,                 &
!    &               mmc_n_angles,mmc_target_angl,mmc_ach_angl,         &
!    &               mmc_ang_sigm,                                      &
!    &               mmc_nvec,mmc_vec,mmc_move_prob,                    &
!    &               mmc_move_cprob,mmc_l_constrains,                   &
!    &               mmc_constrain_type,                                &
!    &               mmc_c_min,mmc_c_max,mmc_c_rad,                     &
!    &               mmc_l_limited,                                     &
!    &               mmc_l_type,mmc_l_center,mmc_l_extend,              &
!    &               mmc_l_lower,mmc_l_upper,                           &
!    &               mmc_len_a,mmc_len_b,mmc_len_m,mmc_len_n,           &
!    &               mmc_buck_a,mmc_buck_rho,mmc_buck_b,                &
!    &               mmc_buck_rmin,mmc_buck_atmin
!
!     COMMON  /mmcav/ n_e_av_p,n_e_av_m,n_e_av_z,e_aver_p,e_aver_m
!
END MODULE mmc_mod
