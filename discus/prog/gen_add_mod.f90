MODULE gen_add_mod
!
!
!     This file contains subroutines for:
!     include file for additional generator matrizes 
!
!     gen_add_n           Number of different additional generators 
!     gen_add (4,4,ng)    List of all additional generators
!
!*****7*****************************************************************
!
INTEGER, PRIVATE    ::  ik
INTEGER, PRIVATE    ::  il
INTEGER, PARAMETER  ::  GEN_ADD_MAX  =  14
!
INTEGER             ::  gen_add_n                  = 0
INTEGER             ::  gen_add_power(GEN_ADD_MAX) = 1
REAL                ::  gen_add(4,4,0:GEN_ADD_MAX) = &
         RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,GEN_ADD_MAX)/),SHAPE(gen_add ))
!
!     COMMON /gen_add/ gen_add_n,gen_add,gen_add_power
!
END MODULE gen_add_mod
