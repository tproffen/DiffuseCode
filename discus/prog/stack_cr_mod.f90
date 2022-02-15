MODULE stack_cr_mod
!-
!
!     Definition of all variables for the microdomain - crystal
!+
USE discus_config_mod
USE stack_mod
!
use precision_mod
!
SAVE
!
INTEGER                             :: ST_MMAX = 1
INTEGER                             :: ST_MAX_SCAT  = 1
!
CHARACTER (LEN=80)                             :: st_name
CHARACTER (LEN=16)                             :: st_spcgr
CHARACTER (LEN=16)                             :: st_spcgr_set
CHARACTER (LEN= 3)                             :: st_set
CHARACTER (LEN=4 ), DIMENSION(:), ALLOCATABLE  :: st_at_lis ! (0:ST_MAX_SCAT)
CHARACTER (LEN=4 ), DIMENSION(:), ALLOCATABLE  :: sa_at_lis ! (0:ST_MAX_SCAT)
integer :: st_iset 
!
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: st_scat  ! (9,0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE :: st_dw    ! (  0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE :: st_occ   ! (  0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE :: sa_dw    ! (  0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE :: sa_occ   ! (  0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE :: st_delfr ! (  0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE :: st_delfi ! (  0:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: st_pos   ! (3,1:ST_MMAX)
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: sa_pos   ! (3,1:ST_MAX_SCAT)
REAL(kind=PREC_DP), DIMENSION(3)                :: st_a0
REAL(kind=PREC_DP), DIMENSION(3)                :: st_win
REAL(kind=PREC_DP)                              :: st_v
REAL(kind=PREC_DP), DIMENSION(3,3)              :: st_gten
REAL(kind=PREC_DP), DIMENSION(3,3,3)            :: st_eps
REAL(kind=PREC_DP), DIMENSION(3)                :: st_ar
REAL(kind=PREC_DP), DIMENSION(3)                :: st_wrez
REAL(kind=PREC_DP)                              :: st_vr
REAL(kind=PREC_DP), DIMENSION(3,3)              :: st_rten
REAL(kind=PREC_DP), DIMENSION(3,3,3)            :: st_reps
REAL(kind=PREC_DP), DIMENSION(3,2)              :: st_dim
!
INTEGER, parameter  ::  st_GEN_ADD_MAX = 192
INTEGER             ::  st_gen_add_n
INTEGER             ::  st_gen_add_power(st_GEN_ADD_MAX)
REAL(kind=PREC_DP)  ::  st_gen_add(4,4,0:st_GEN_ADD_MAX)
!
INTEGER, parameter  ::  st_SYM_ADD_MAX = 192
INTEGER             ::  st_sym_add_n
INTEGER             ::  st_sym_add_power(st_sym_ADD_MAX)
REAL(kind=PREC_DP)  ::  st_sym_add(4,4,0:st_sym_ADD_MAX)
!
INTEGER                              :: sa_natoms
INTEGER                              :: st_nscat
INTEGER                              :: sa_nscat
INTEGER , DIMENSION(:), ALLOCATABLE  :: st_iscat  ! (1:ST_MMAX)
INTEGER , DIMENSION(:), ALLOCATABLE  :: st_mole   ! (1:ST_MMAX)
INTEGER , DIMENSION(:,:), ALLOCATABLE:: st_surf   ! (1:ST_MMAX)
REAL(kind=PREC_DP)    , DIMENSION(:,:), ALLOCATABLE:: st_magn   ! (1:ST_MMAX)
INTEGER , DIMENSION(:), ALLOCATABLE  :: st_prop   ! (1:ST_MMAX)
INTEGER , DIMENSION(:), ALLOCATABLE  :: sa_iscat  ! (1:ST_MAX_SCAT)
INTEGER , DIMENSION(:), ALLOCATABLE  :: sa_prop   ! (1:ST_MAX_SCAT)
INTEGER                              :: st_ndel
INTEGER                              :: st_natoms
INTEGER                              :: st_spcgr_ianz
INTEGER                              :: st_spcgr_para
LOGICAL                              :: st_magnetic
INTEGER                              :: st_cr_size_of
!
END MODULE stack_cr_mod
