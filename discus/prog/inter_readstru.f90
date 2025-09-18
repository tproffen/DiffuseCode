MODULE inter_readstru
!
USE precision_mod
IMPLICIT NONE
!
INTEGER rd_NMAX 
INTEGER rd_MAXSCAT 
!
CHARACTER (LEN=PREC_STRING)           :: rd_strucfile 
CHARACTER (LEN=  80)                  :: rd_cr_name 
CHARACTER (LEN=  16)                  :: rd_cr_spcgr 
REAL(kind=PREC_DP)  , DIMENSION(3)    :: rd_cr_a0
REAL(kind=PREC_DP)  , DIMENSION(3)    :: rd_cr_ar
REAL(kind=PREC_DP)  , DIMENSION(3)    :: rd_cr_win
REAL(kind=PREC_DP)  , DIMENSION(3)    :: rd_cr_wrez
REAL(kind=PREC_DP)  , DIMENSION(3,3)  :: rd_cr_gten
REAL(kind=PREC_DP)  , DIMENSION(3,3)  :: rd_cr_rten
INTEGER             , DIMENSION(3)    :: rd_cr_icc
INTEGER             , DIMENSION(0:1)  :: rd_cr_sel_prop
REAL(kind=PREC_DP)  , DIMENSION(3,3)  :: rd_cr_fmat
REAL(kind=PREC_DP)  , DIMENSION(3,3)  :: rd_cr_gmat
REAL(kind=PREC_DP)  , DIMENSION(3,3,3):: rd_cr_eps
REAL(kind=PREC_DP)  , DIMENSION(3,3,3):: rd_cr_reps
REAL(kind=PREC_DP)  , DIMENSION(4,4)  :: rd_tran_g
REAL(kind=PREC_DP)  , DIMENSION(4,4)  :: rd_tran_gi
REAL(kind=PREC_DP)  , DIMENSION(4,4)  :: rd_tran_f
REAL(kind=PREC_DP)  , DIMENSION(4,4)  :: rd_tran_fi
REAL(kind=PREC_DP)                    :: rd_cr_vol
REAL(kind=PREC_DP)                    :: rd_cr_vr
INTEGER                               :: rd_cr_spcgrno
INTEGER                               :: rd_cr_syst
INTEGER                               :: rd_cr_natoms
INTEGER                               :: rd_cr_ncatoms
INTEGER                               :: rd_cr_n_real_atoms
INTEGER                               :: rd_cr_nscat 
LOGICAL :: rd_cr_acentric
LOGICAL :: rd_cr_cartesian
LOGICAL :: rd_cr_newtype
REAL(kind=PREC_DP)  , DIMENSION(:),   ALLOCATABLE :: rd_cr_dw     ! (0:MAXSCAT) 
REAL(kind=PREC_DP)  , DIMENSION(:),   ALLOCATABLE :: rd_cr_occ    ! (0:MAXSCAT) 
CHARACTER (LEN=   4), DIMENSION(:),   ALLOCATABLE :: rd_cr_at_lis ! (0:MAXSCAT) 
CHARACTER (LEN=   4), DIMENSION(:),   ALLOCATABLE :: rd_cr_as_lis ! (0:MAXSCAT) 
CHARACTER (LEN=   4), DIMENSION(:),   ALLOCATABLE :: rd_cr_at_equ ! (0:MAXSCAT) 
LOGICAL             , DIMENSION(:),   ALLOCATABLE :: rd_cr_scat_equ ! (0:MAXSCAT) 
LOGICAL             , DIMENSION(:),   ALLOCATABLE :: rd_cr_scat_int ! (0:MAXSCAT) 
LOGICAL             , DIMENSION(:),   ALLOCATABLE :: rd_cr_delf_int ! (0:MAXSCAT) 
REAL(kind=PREC_DP)  , DIMENSION(:),   ALLOCATABLE :: rd_cr_delfi    ! (0:MAXSCAT) 
REAL(kind=PREC_DP)  , DIMENSION(:),   ALLOCATABLE :: rd_cr_delfr    ! (0:MAXSCAT) 
REAL(kind=PREC_DP)  , DIMENSION(:,:), ALLOCATABLE :: rd_cr_scat   ! (11,0_MAXSCAT)
REAL(kind=PREC_DP)  , DIMENSION(:,:), ALLOCATABLE :: rd_cr_pos    ! (1:3,1:NMAX)
INTEGER             , DIMENSION(:,:), ALLOCATABLE :: rd_cr_iscat  ! (1:NMAX)
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_cr_prop   ! (1:NMAX)
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_cr_mole   ! (1:NMAX)
INTEGER             , DIMENSION(:,:), ALLOCATABLE :: rd_cr_surf   ! (1:NMAX)
REAL(kind=PREC_DP)  , DIMENSION(:,:), ALLOCATABLE :: rd_cr_magn   ! (1:NMAX)
real(kind=PREC_DP)  , dimension(  :), allocatable :: rd_cr_valu   ! (1:NMAX)
REAL(kind=PREC_DP)  , DIMENSION(3,2)  :: rd_cr_dim
REAL(kind=PREC_DP)  , DIMENSION(3,2)  :: rd_cr_dim0
INTEGER                               :: rd_as_natoms 
CHARACTER (LEN=   4), DIMENSION(:),   ALLOCATABLE :: rd_as_at_lis ! (0:MAXSCAT) 
REAL(kind=PREC_DP)  , DIMENSION(:),   ALLOCATABLE :: rd_as_dw     ! (0:MAXSCAT) 
REAL(kind=PREC_DP)  , DIMENSION(:,:), ALLOCATABLE :: rd_as_pos    ! (3, MAXSCAT) 
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_as_iscat  ! (1:MAXSCAT)
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_as_prop   ! (1:MAXSCAT)
INTEGER             , DIMENSION(3)    :: rd_sav_ncell ! (3) 
LOGICAL                               :: rd_sav_r_ncell 
INTEGER                               :: rd_sav_ncatoms 
INTEGER                               :: rd_spcgr_ianz 
INTEGER                               :: rd_spcgr_para 
!                                                                       
! SAVE FLAGS
!
LOGICAL                           ::  rd_sav_scat      = .true.
LOGICAL                           ::  rd_sav_adp       = .true.
LOGICAL                           ::  rd_sav_gene      = .true.
LOGICAL                           ::  rd_sav_symm      = .true.
LOGICAL                           ::  rd_sav_obje      = .true.
LOGICAL                           ::  rd_sav_doma      = .true.
LOGICAL                           ::  rd_sav_mole      = .true.
LOGICAL                           ::  rd_sav_w_ncell   = .true.
LOGICAL                           ::  rd_sav_prop      = .true.
INTEGER, DIMENSION(0:1)           ::  rd_sav_sel_prop  = 0
INTEGER                           ::  rd_n_latom       = 0
LOGICAL, DIMENSION(:),ALLOCATABLE ::  rd_sav_latom ! (0:MAXSCAT)

!
END MODULE inter_readstru
