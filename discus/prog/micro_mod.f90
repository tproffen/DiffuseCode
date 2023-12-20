MODULE micro_mod
!-
!     Definition of all variables for the microdomain - crystal
!+
USE discus_config_mod
use precision_mod
!
SAVE
!
INTEGER      :: MK_MAX_SCAT
INTEGER      :: MK_MAX_ATOM
!
CHARACTER(LEN=80)  :: mk_name  = ' '
CHARACTER(LEN=16)  :: mk_spcgr = 'P1'
CHARACTER(LEN=16)  :: mk_spcgr_set = 'P1'
CHARACTER(LEN= 3)  :: mk_set   = 'abc'
INTEGER            :: mk_iset  =  1
CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE ::  mk_at_lis  ! (0:MK_MAX_SCAT)
!
REAL(kind=PREC_DP)            , DIMENSION(:), ALLOCATABLE ::  mk_dw      ! (0:MK_MAX_SCAT)
integer                       , DIMENSION(:), ALLOCATABLE ::  mk_is_sym  ! (0:ncatoms)
REAL(kind=PREC_DP)            , DIMENSION(:), ALLOCATABLE ::  mk_occ     ! (0:MK_MAX_SCAT)
REAL(kind=PREC_DP)            , DIMENSION(:,:), ALLOCATABLE ::  mk_anis    ! (0:MK_MAX_SCAT)
REAL(kind=PREC_DP)            , DIMENSION(:,:), ALLOCATABLE   :: mk_anis_full    ! (0:MK_MAX_SCAT)
REAL(kind=PREC_DP)            , DIMENSION(:,:,:), ALLOCATABLE :: mk_prin    ! (0:MK_MAX_SCAT)
REAL(kind=PREC_DP)  :: mk_a0(3)    = (/1.0D0, 1.0D0, 1.0D0/)
REAL(kind=PREC_DP)  :: mk_win(3)   = (/90.0D0, 90.0D0, 90.0D0/)
REAL(kind=PREC_DP)  :: mk_dim(3,2) = 0.0D0
!
INTEGER             :: ik,il
INTEGER             :: mk_nscat    = 0
INTEGER             :: mk_natoms   = 0
INTEGER             :: mk_nanis    = 0   ! Number of ADPs
!
INTEGER             :: mk_spcgr_ianz = 1,mk_spcgr_para = 0
INTEGER             :: mk_spcgr_no = 1
INTEGER, PARAMETER  ::  mk_GEN_ADD_MAX  =  192
INTEGER, PARAMETER  ::  mk_SYM_ADD_MAX  =  192
!
INTEGER             ::  mk_gen_add_n                  = 0
INTEGER             ::  mk_gen_add_power(mk_GEN_ADD_MAX) = 1
REAL(kind=PREC_DP)  ::  mk_gen_add(4,4,0:mk_GEN_ADD_MAX) = &
     RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,mk_GEN_ADD_MAX)/),(/4,4,mk_GEN_ADD_MAX+1/))
!
INTEGER             ::  mk_sym_add_n                  = 0
INTEGER             ::  mk_sym_add_power(mk_SYM_ADD_MAX) = 1
REAL(kind=PREC_DP)  ::  mk_sym_add(4,4,0:mk_SYM_ADD_MAX) = &
     RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,mk_SYM_ADD_MAX)/),(/4,4,mk_SYM_ADD_MAX+1/))
!
LOGICAL             :: mk_infile_internal = .FALSE.  ! File with domain atoms in stored internally
INTEGER             :: mk_iatom           = 1        ! Current atom number in internal file
!
END MODULE micro_mod
