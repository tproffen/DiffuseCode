MODULE wink_mod
!+
!     This file the definitions for radian, 2 pi 
!-
   USE precision_mod
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   REAL(PREC_DP), PARAMETER :: pi  = 3.1415926535897932384626433832795028841971693993751D0
   REAL(PREC_DP), PARAMETER :: zpi = 2.0D0 * pi
   REAL(PREC_DP), PARAMETER :: fpi = 4.0D0 * pi
   REAL(PREC_DP), PARAMETER :: rad = pi/180.D0
!
!
END MODULE wink_mod
