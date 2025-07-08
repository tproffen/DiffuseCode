MODULE surface_mod
!+
!
!     Surface related definitions
!-
!
use precision_mod
!
SAVE
!
INTEGER, PARAMETER  ::  SURF_MAXTYPE  =  7
INTEGER, PARAMETER  ::  SURF_NONE     =  0
INTEGER, PARAMETER  ::  SURF_PLANE    =  1
INTEGER, PARAMETER  ::  SURF_SPHERE   =  2
INTEGER, PARAMETER  ::  SURF_CYLINDER =  3
INTEGER, PARAMETER  ::  SURF_EDGE     =  4
INTEGER, PARAMETER  ::  SURF_CORNER   =  5
INTEGER, PARAMETER  ::  SURF_LOCAL    =  6
!
INTEGER, PARAMETER  ::  SURF_SURFACE  =  0
INTEGER, PARAMETER  ::  SURF_DOMAIN   =  1
!
INTEGER             ::  SURF_MAXSCAT  =  0
!
REAL(kind=PREC_DP)  ::  SURF_DIST_DEF =  2.55
!
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE  ::  surf_ex_dist   !(0:MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE  ::  surf_in_dist   !(0:MAXSCAT)
!
logical           , DIMENSION(:), ALLOCATABLE  ::  surf_original  !(0:MAXSCAT)  Old atoms
logical           , DIMENSION(:), ALLOCATABLE  ::  surf_replace   !(0:MAXSCAT)  to be replaced by these
integer                                        ::  surf_n_repl    ! number of replacement options
!
LOGICAL             ::  surf_local_new = .TRUE.
logical             ::  surf_boundary 
!
!
END MODULE surface_mod
