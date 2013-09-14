MODULE config_mod
!
IMPLICIT NONE
PUBLIC
SAVE
!
!#######################################################################
!       MAXSCAT Configuration file
!#######################################################################
!
!     MAXDAT      : Maximum number of PDF data points
!     MAXDSET     : Maximum number of PDF data sets
!     MAXELEM     : Maximum number of element
!     MAXPARA     : Maximum number of fit PARAMETERs
!
!#######################################################################
!
   INTEGER, PARAMETER :: MAXDAT  = 10010
   INTEGER, PARAMETER :: MAXDSET =     2
   INTEGER, PARAMETER :: MAXELEM =    20
   INTEGER, PARAMETER :: MAXPARA =     2
!
!
END MODULE config_mod
