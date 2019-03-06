MODULE surface_mod
!+
!
!     Surface related definitions
!-
!
SAVE
!
INTEGER, PARAMETER  ::  SURF_MAXTYPE  =  5
INTEGER, PARAMETER  ::  SURF_NONE     =  0
INTEGER, PARAMETER  ::  SURF_PLANE    =  1
INTEGER, PARAMETER  ::  SURF_SPHERE   =  2
INTEGER, PARAMETER  ::  SURF_CYLINDER =  3
INTEGER, PARAMETER  ::  SURF_EDGE     =  4
INTEGER, PARAMETER  ::  SURF_CORNER   =  5
!
INTEGER, PARAMETER  ::  SURF_SURFACE  =  0
INTEGER, PARAMETER  ::  SURF_DOMAIN   =  1
!
INTEGER             ::  SURF_MAXSCAT  =  0
!
REAL                ::  SURF_DIST_DEF =  2.55
!
REAL, DIMENSION(:), ALLOCATABLE  ::  surf_ex_dist   !(0:MAXSCAT)
REAL, DIMENSION(:), ALLOCATABLE  ::  surf_in_dist   !(0:MAXSCAT)
!
INTEGER             ::  surf_size_of
!
END MODULE surface_mod
