MODULE micro_mod
!-
!     Definition of all variables for the microdomain - crystal
!+
USE config_mod
!
SAVE
!
INTEGER      :: MK_MAX_SCAT
INTEGER      :: MK_MAX_ATOM
!
CHARACTER(LEN=80) mk_name
CHARACTER(LEN=16) mk_spcgr
CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE ::  mk_at_lis  ! (0:MK_MAX_SCAT)
!
REAL            , DIMENSION(:), ALLOCATABLE ::  mk_dw      ! (0:MK_MAX_SCAT)
REAL                :: mk_a0(3)
REAL                :: mk_win(3)
REAL                :: mk_dim(3,2)
!
INTEGER             :: ik,il
INTEGER           :: mk_nscat
INTEGER             :: mk_natoms
!
INTEGER             :: mk_spcgr_ianz,mk_spcgr_para
INTEGER             :: mk_spcgr_no
INTEGER, PARAMETER  ::  mk_GEN_ADD_MAX  =  14
INTEGER, PARAMETER  ::  mk_SYM_ADD_MAX  =  14
!
INTEGER             ::  mk_gen_add_n                  = 0
INTEGER             ::  mk_gen_add_power(mk_GEN_ADD_MAX) = 1
REAL                ::  mk_gen_add(4,4,0:mk_GEN_ADD_MAX) = &
     RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,mk_GEN_ADD_MAX)/),(/4,4,mk_GEN_ADD_MAX+1/))
!
INTEGER             ::  mk_sym_add_n                  = 0
INTEGER             ::  mk_sym_add_power(mk_SYM_ADD_MAX) = 1
REAL                ::  mk_sym_add(4,4,0:mk_SYM_ADD_MAX) = &
     RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,mk_SYM_ADD_MAX)/),(/4,4,mk_SYM_ADD_MAX+1/))
!
LOGICAL             :: mk_infile_internal  ! File with domain atoms in stored internally
INTEGER             :: mk_iatom            ! Current atom number in internal file
INTEGER             :: mic_size_of
!
END MODULE micro_mod
