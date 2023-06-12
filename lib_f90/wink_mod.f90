module wink_mod
!+
!     This file the definitions for radian, 2 pi 
!-
use precision_mod
implicit none
public
save
!
real(PREC_DP), parameter :: pi  = 3.1415926535897932384626433832795028841971693993751D0
real(PREC_DP), parameter :: zpi = 2.0D0 * pi
real(PREC_DP), parameter :: fpi = 4.0D0 * pi
real(PREC_DP), parameter :: rad = pi/180.D0
!
!
end module wink_mod
