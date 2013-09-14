MODULE unitcell_mod
!
!
!     This file contains subroutines for:
!     input conditions for the unit cell files
!
!     gen_sta             Sequence status of generators
!
!*****7*****************************************************************
!
SAVE
!
INTEGER, PARAMETER  ::     GEN_CENTER  =  0
INTEGER, PARAMETER  ::     GEN_SYMM    =  1
!
INTEGER             ::  gen_sta = GEN_SYMM
!
END MODULE unitcell_mod
