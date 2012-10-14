MODULE chem_mod
!+
!
!
!     This file contains COMMON variables for chem routines
!-
USE config_mod
!
SAVE
!
INTEGER, PARAMETER      :: CHEM_NONE  = 0
INTEGER, PARAMETER      :: CHEM_DIST  = 1
INTEGER, PARAMETER      :: CHEM_VEC   = 2
INTEGER, PARAMETER      :: CHEM_ANG   = 3
INTEGER, PARAMETER      :: CHEM_ENVIR = 4
INTEGER, PARAMETER      :: CHEM_RANGE = 5
INTEGER, PARAMETER      :: CHEM_CON   = 6
!
!
INTEGER                 :: CHEM_MAXAT_CELL   = 1
INTEGER                 :: CHEM_MAX_AVE_ATOM = 1
INTEGER                 :: CHEM_MAX_VEC      = 1
INTEGER                 :: CHEM_MAX_ANG      = 1
INTEGER                 :: CHEM_MAX_RAN      = 0
INTEGER                 :: CHEM_MAX_CON      = 1
INTEGER                 :: CHEM_MAX_ENV      = 1
!
CHARACTER(LEN=80)                        :: chem_fname     = 'blen.xy'
INTEGER                                  :: chem_bin       = 601
INTEGER                                  :: chem_ncor      = 1
REAL, DIMENSION(2)                       :: chem_blen_cut  = (/1.5,  7.5/)
REAL, DIMENSION(2)                       :: chem_bang_cut  = (/0.0,180.0/)
LOGICAL                                  :: chem_quick     = .true.
LOGICAL                                  :: chem_cluster   = .false.
LOGICAL, DIMENSION(3)                    :: chem_period    = .true.
LOGICAL                                  :: chem_sel_atom  = .true.
!
! Vector definitions
INTEGER, DIMENSION(:,:), ALLOCATABLE     :: chem_cvec        !  (5,CHEM_MAX_VEC)
INTEGER, DIMENSION(:,:), ALLOCATABLE     :: chem_use_vec     !  (CHEM_MAX_VEC,CHEM_MAX_COR)
!
! Range definitions
INTEGER, DIMENSION(:,:)  , ALLOCATABLE   :: chem_use_ran     !  (CHEM_MAX_RAN,CHEM_MAX_COR)
INTEGER, DIMENSION(:  ,:), ALLOCATABLE   :: chem_cran_cent   !  (0:CHEM_MAX_ATOM,CHEM_MAX_RAN)
INTEGER, DIMENSION(:  ,:), ALLOCATABLE   :: chem_cran_neig   !  (0:CHEM_MAX_ATOM,CHEM_MAX_RAN)
INTEGER, DIMENSION(:)    , ALLOCATABLE   :: chem_cran_nuvw   !  (CHEM_MAX_RAN)
INTEGER, DIMENSION(:)    , ALLOCATABLE   :: chem_cran_nshort !  (CHEM_MAX_RAN)
REAL   , DIMENSION(:,:,:), ALLOCATABLE   :: chem_cran_uvw    !  (3,48,CHEM_MAX_RAN)
REAL   , DIMENSION(:)    , ALLOCATABLE   :: chem_cran_sig    !  (CHEM_MAX_RAN)
REAL   , DIMENSION(:)    , ALLOCATABLE   :: chem_cran_wsig   !  (CHEM_MAX_RAN)
REAL   , DIMENSION(:)    , ALLOCATABLE   :: chem_cran_rmax   !  (CHEM_MAX_RAN)
REAL   , DIMENSION(:)    , ALLOCATABLE   :: chem_cran_rmin   !  (CHEM_MAX_RAN)
LOGICAL, DIMENSION(:)    , ALLOCATABLE   :: chem_cran_cang   !  (CHEM_MAX_RAN)
LOGICAL, DIMENSION(:)    , ALLOCATABLE   :: chem_cran_lsym   !  (CHEM_MAX_RAN)
LOGICAL, DIMENSION(:)    , ALLOCATABLE   :: chem_cran_short  !  (CHEM_MAX_RAN)
!
! Angle definitions
INTEGER, DIMENSION(:,:), ALLOCATABLE     :: chem_cwin        !  (9,CHEM_MAX_ANG)
INTEGER, DIMENSION(:,:), ALLOCATABLE     :: chem_use_win     !  (CHEM_MAX_ANG,CHEM_MAX_COR)
!
! Connectivity definitions
INTEGER, DIMENSION(:,:), ALLOCATABLE     :: chem_ccon        !  (2,CHEM_MAX_ANG)
INTEGER, DIMENSION(:,:), ALLOCATABLE     :: chem_use_con     !  (CHEM_MAX_CON,CHEM_MAX_COR)
!
! Environment definitions
!
INTEGER, DIMENSION(CHEM_MAX_BIN)                       :: chem_hist        !  (CHEM_MAX_BIN)
!
! Average crystal structure
INTEGER, DIMENSION(:)  , ALLOCATABLE                   :: chem_ave_n       !  (MAXAT_CELL)
INTEGER, DIMENSION(:,:), ALLOCATABLE                   :: chem_ave_iscat   !  (MAXAT_CELL,CHEM_MAX_ATOM)
REAL   , DIMENSION(:,:), ALLOCATABLE                   :: chem_ave_pos     !  (3,MAXAT_CELL)
REAL   , DIMENSION(:,:), ALLOCATABLE                   :: chem_ave_sig     !  (3,MAXAT_CELL)
REAL   , DIMENSION(:,:), ALLOCATABLE                   :: chem_ave_bese    !  (MAXAT_CELL,CHEM_MAX_ATOM)
!
! Displacement correlations
REAL, DIMENSION(:,:,:), ALLOCATABLE                    :: chem_disp_ave    !  (CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
REAL, DIMENSION(:,:,:), ALLOCATABLE                    :: chem_disp_sig    !  (CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
! Environment interactions
INTEGER, DIMENSION(:,:), ALLOCATABLE                   :: chem_use_env     !  (CHEM_MAX_ENV,CHEM_MAX_COR)
INTEGER, DIMENSION(:,:), ALLOCATABLE                   :: chem_cenv        !  (0:MAX_ATOM_ENV,CHEM_MAX_ENV)
INTEGER, DIMENSION(:)  , ALLOCATABLE                   :: chem_env_neig    !  (CHEM_MAX_ENV)
REAL,    DIMENSION(:)  , ALLOCATABLE                   :: chem_rmax_env    !  (CHEM_MAX_ENV)
REAL,    DIMENSION(:)  , ALLOCATABLE                   :: chem_rmin_env    !  (CHEM_MAX_ENV)
!
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_nran        !  (CHEM_MAX_COR)
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_nvec        !  (CHEM_MAX_COR)
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_ncon        !  (CHEM_MAX_COR)
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_nwin        !  (CHEM_MAX_COR)
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_ctyp        !  (CHEM_MAX_COR)
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_nnei        !  (CHEM_MAX_COR)
INTEGER, DIMENSION(CHEM_MAX_COR)                       :: chem_nenv        !  (CHEM_MAX_COR)
!REAL, DIMENSION(3,12,CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT) :: chem_vect_ave    !  (3,12,CHEM_MAX_COR,0:MAXSCAT,0:MAXSCAT)
!REAL, DIMENSION(3,12,CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT) :: chem_vect_sig    !  (3,12,CHEM_MAX_COR,0:MAXSCAT,0:MAXSCAT)
REAL, DIMENSION(3,48,CHEM_MAX_COR)                     :: chem_neig        !  (3,48,CHEM_MAX_COR)
REAL, DIMENSION(3,2,CHEM_MAX_COR)                      :: chem_dir         !  (3,2,CHEM_MAX_COR)
REAL, DIMENSION(CHEM_MAX_COR)                          :: chem_rmax        !  (CHEM_MAX_COR)
REAL, DIMENSION(CHEM_MAX_COR)                          :: chem_rmin        !  (CHEM_MAX_COR)
REAL, DIMENSION(CHEM_MAX_COR)                          :: chem_freq_sigma  !  (CHEM_MAX_COR)
REAL, DIMENSION(CHEM_MAX_COR)                          :: chem_wink_sigma  !  (CHEM_MAX_COR)
LOGICAL, DIMENSION(CHEM_MAX_COR)                       :: chem_cang        !  (CHEM_MAX_COR)
LOGICAL, DIMENSION(CHEM_MAX_COR)                       :: chem_ldall       !  (CHEM_MAX_COR)
!
!
INTEGER :: chem_aver_size_of = 0
INTEGER :: chem_ang_size_of  = 0
INTEGER :: chem_disp_size_of = 0
INTEGER :: chem_env_size_of  = 0
INTEGER :: chem_vec_size_of  = 0
INTEGER :: chem_con_size_of  = 0
INTEGER :: chem_ran_size_of  = 0

!
!     COMMON      /chbl/ chem_fname,chem_hist,                          &
!    &                   chem_blen_cut,chem_bang_cut,chem_bin,          &
!    &                   chem_quick,chem_freq_sigma,chem_wink_sigma,    &
!    &                   chem_ave_pos,chem_ave_sig,chem_ave_n,          &
!    &                   chem_ave_iscat,chem_ave_bese,chem_cang,        &
!    &                   chem_cvec,chem_nvec,chem_neig,                 &
!    &                   chem_cran_uvw,chem_cran_sig,chem_cran_wsig,    &
!    &                   chem_cran_cent,chem_cran_neig,chem_cran_cang,  &
!    &                   chem_cran_lsym,chem_cran_short,                &
!    &                   chem_cran_nshort,                              &
!    &                   chem_nran,chem_use_ran,chem_cran_nuvw,         &
!    &                   chem_ctyp,chem_use_vec,chem_ncor,chem_dir,     &
!    &                   chem_nnei,chem_rmax,chem_rmin,                 &
!    &                   chem_cran_rmax,chem_cran_rmin,                 &
!    &                   chem_disp_ave,chem_disp_sig,                   &
!    &                   chem_vect_ave,chem_vect_sig,                   &
!    &                   chem_period,                                   &
!    &                   chem_sel_atom,chem_ldall,                      &
!    &                   chem_cwin,chem_nwin,chem_use_win,              &
!    &                   chem_nenv,chem_use_env,chem_cenv,              &
!    &                   chem_env_neig,chem_rmin_env,chem_rmax_env,     &
!    &                   chem_cluster
!
!     COMMON      /bv/      bv_index,bv_r0,bv_b
!
END MODULE chem_mod
