module superspace_mod
!-
!  Variables for superspace applications
!+
!
use precision_mod
!
integer, parameter :: SUP_DISP = 0
integer, parameter :: SUP_DENS = 1
!
integer, parameter :: SUP_SINE = 1
integer, parameter :: SUP_USER = 2
!
character(len=PREC_STRING) :: sup_file
character(len=4)  , dimension(:)  , allocatable :: sup_atom   ! "Atom" name
character(len=4)  , dimension(:)  , allocatable :: sup_repl   ! Replace site by this atom type
integer           , dimension(:)  , allocatable :: sup_irepl  ! Replace site by this atom type
integer           , dimension(:)  , allocatable :: sup_char   ! Character DISPLACEMENT / DENSITY
integer           , dimension(:)  , allocatable :: sup_func   ! Function  : sine, user
real(kind=PREC_DP), dimension(:,:), allocatable :: sup_ampl   ! Displacement amplitude vector
real(kind=PREC_DP), dimension(:)  , allocatable :: sup_phase  ! Phase in range [0,1]
real(kind=PREC_DP), dimension(:,:), allocatable :: sup_prob   ! Density propabilities [min,max]
!
real(kind=PREC_DP), dimension(3) :: sup_qvec   ! Satellite vector
!
end module superspace_mod
