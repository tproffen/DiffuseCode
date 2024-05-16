module super_waves_mod
!
contains
!
!*******************************************************************************
!
real(kind=PREC_DP    ) function sup_fun_sine(amp, average, arg1, arg2)
!-
!  Sine type wave function : amp * sin(2PI*arg) + average
!+
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: amp
real(kind=PREC_DP), intent(in) :: average
real(kind=PREC_DP), intent(in) :: arg1
real(kind=PREC_DP), intent(in) :: arg2
!
sup_fun_sine = amp*sin(zpi*arg1) + average
!
end function sup_fun_sine
!
!*******************************************************************************
!
real(kind=PREC_DP    ) function sup_fun_cren(amp, average, arg1, arg2)
!-
!  Crennel type wave function : amp * box(frac(arg)) + average
!+
!
use lib_functions_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: amp
real(kind=PREC_DP), intent(in) :: average
real(kind=PREC_DP), intent(in) :: arg1
real(kind=PREC_DP), intent(in) :: arg2
!
real(kind=PREC_DP), parameter  :: TOL=1.0D-8
real(kind=PREC_DP)             :: arg
!
arg = frac(arg1)
if(arg<-TOL       ) then
   arg = frac(arg + 1.0_PREC_DP)
endif
!
if(arg  < arg2 .or. arg > 1.0_PREC_DP-TOL-arg2) then
   sup_fun_cren = average + amp
else
   sup_fun_cren = average - amp
endif
!
end function sup_fun_cren
!
!*******************************************************************************
!
end module super_waves_mod
