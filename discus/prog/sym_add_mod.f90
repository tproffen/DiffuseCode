MODULE sym_add_mod
!
!
!     This file contains subroutines for:
!     include file for additional symmetry matrizes 
!
!     sym_add_n           Number of different additional symmetries 
!     sym_add (4,4,ng)    List of all additional symmetries
!
!*****7*****************************************************************
!
SAVE
!
INTEGER, PRIVATE    :: ik
INTEGER, PRIVATE    :: il
INTEGER, PARAMETER  ::  SYM_ADD_MAX  =  14
!
INTEGER                               ::  sym_add_n     = 0
INTEGER, DIMENSION(SYM_ADD_MAX)       ::  sym_add_power = 1
REAL   , DIMENSION(4,4,0:SYM_ADD_MAX) ::  sym_add       = &
         RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,SYM_ADD_MAX)/),(/4,4,SYM_ADD_MAX+1/))
!
END MODULE sym_add_mod
