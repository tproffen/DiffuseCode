module recipro_mod
!+
!
!     This file contains variables for inverse Fourier and
!     Patterson input
!-
use precision_mod
!
save
!
integer, parameter  ::  REC_MAX_SYM  =  48
!
integer                             :: rec_n_sym
real(kind=PREC_DP), dimension(4,4,REC_MAX_SYM) :: rec_sym
!
end module recipro_mod
