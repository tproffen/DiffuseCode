MODULE surface_mod
!+
!
!     Surface related definitions
!-
!
SAVE
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
