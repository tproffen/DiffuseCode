MODULE wink_mod
!+
!     This file the definitions for radian, 2 pi 
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   REAL, PARAMETER :: pi  = 3.1415926535897932384626433832795028841971693993751
   REAL, PARAMETER :: zpi = 2.0 * pi
   REAL, PARAMETER :: fpi = 4.0 * pi
   REAL, PARAMETER :: rad = pi/180.
!
!
END MODULE wink_mod
