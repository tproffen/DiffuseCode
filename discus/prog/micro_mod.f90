MODULE micro_mod
!-
!     Definition of all variables for the microdomain - crystal
!+
USE config_mod
!
SAVE
!
CHARACTER*80 mk_name
CHARACTER*16 mk_spcgr
CHARACTER*4  mk_at_lis(0:MK_MAX_SCAT)
!
REAL         mk_dw(0:MK_MAX_SCAT)
REAL         mk_a0(3)
REAL         mk_win(3)
REAL         mk_dim(3,2)
!
INTEGER*4    mk_nscat
INTEGER      mk_natoms
!
INTEGER      mk_spcgr_ianz,mk_spcgr_para
!
END MODULE micro_mod
