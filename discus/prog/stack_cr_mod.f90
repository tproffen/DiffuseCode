MODULE stack_cr_mod
!-
!
!     Definition of all variables for the microdomain - crystal
!+
USE discus_config_mod
USE stack_mod
!
SAVE
!
INTEGER                             :: ST_MMAX = 1
INTEGER                             :: ST_MAX_SCAT  = 1
!
CHARACTER (LEN=80)                             :: st_name
CHARACTER (LEN=16)                             :: st_spcgr
CHARACTER (LEN= 3)                             :: st_set
CHARACTER (LEN=4 ), DIMENSION(:), ALLOCATABLE  :: st_at_lis ! (0:ST_MAX_SCAT)
CHARACTER (LEN=4 ), DIMENSION(:), ALLOCATABLE  :: sa_at_lis ! (0:ST_MAX_SCAT)
!
REAL   , DIMENSION(:,:), ALLOCATABLE :: st_scat  ! (9,0:ST_MAX_SCAT)
REAL   , DIMENSION(  :), ALLOCATABLE :: st_dw    ! (  0:ST_MAX_SCAT)
REAL   , DIMENSION(  :), ALLOCATABLE :: st_occ   ! (  0:ST_MAX_SCAT)
REAL   , DIMENSION(  :), ALLOCATABLE :: sa_dw    ! (  0:ST_MAX_SCAT)
REAL   , DIMENSION(  :), ALLOCATABLE :: sa_occ   ! (  0:ST_MAX_SCAT)
REAL   , DIMENSION(  :), ALLOCATABLE :: st_delfr ! (  0:ST_MAX_SCAT)
REAL   , DIMENSION(  :), ALLOCATABLE :: st_delfi ! (  0:ST_MAX_SCAT)
REAL   , DIMENSION(:,:), ALLOCATABLE :: st_pos   ! (3,1:ST_MMAX)
REAL   , DIMENSION(:,:), ALLOCATABLE :: sa_pos   ! (3,1:ST_MAX_SCAT)
REAL   , DIMENSION(3)                :: st_a0
REAL   , DIMENSION(3)                :: st_win
REAL                                 :: st_v
REAL   , DIMENSION(3,3)              :: st_gten
REAL   , DIMENSION(3,3,3)            :: st_eps
REAL   , DIMENSION(3)                :: st_ar
REAL   , DIMENSION(3)                :: st_wrez
REAL                                 :: st_vr
REAL   , DIMENSION(3,3)              :: st_rten
REAL   , DIMENSION(3,3,3)            :: st_reps
REAL   , DIMENSION(3,2)              :: st_dim
!
INTEGER                              :: sa_natoms
INTEGER                              :: st_nscat
INTEGER                              :: sa_nscat
INTEGER , DIMENSION(:), ALLOCATABLE  :: st_iscat  ! (1:ST_MMAX)
INTEGER , DIMENSION(:), ALLOCATABLE  :: st_mole   ! (1:ST_MMAX)
INTEGER , DIMENSION(:,:), ALLOCATABLE:: st_surf   ! (1:ST_MMAX)
INTEGER , DIMENSION(:), ALLOCATABLE  :: st_prop   ! (1:ST_MMAX)
INTEGER , DIMENSION(:), ALLOCATABLE  :: sa_iscat  ! (1:ST_MAX_SCAT)
INTEGER , DIMENSION(:), ALLOCATABLE  :: sa_prop   ! (1:ST_MAX_SCAT)
INTEGER                              :: st_ndel
INTEGER                              :: st_natoms
INTEGER                              :: st_spcgr_ianz
INTEGER                              :: st_spcgr_para
INTEGER                              :: st_cr_size_of
!
END MODULE stack_cr_mod
