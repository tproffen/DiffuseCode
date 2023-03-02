module lib_functions_mod
!
!  Low level function of generic use
!+
!
contains
!
!###############################################################################
!
function sinc(x)  result(val)
!-
!  sinc function:   sin(x)/x
!+
use precision_mod
!
implicit none
!
real(kind=PREC_DP) :: val
!
real(kind=PREC_DP), intent(in) :: x
!
real(kind=PREC_DP), parameter :: TOL=1.0D-7
!
if(abs(x)<TOL) then
   val = 1.0D0
else
   val = sin(x)/x
endif
!
end function sinc
!
!###############################################################################
!
end module lib_functions_mod
