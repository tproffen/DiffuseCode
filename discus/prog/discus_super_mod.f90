module superspace_mod
!-
!  Variables for superspace applications
!+
!
use precision_mod
!
integer, parameter :: SUP_DISP = 1
integer, parameter :: SUP_DENS = 2
!
integer, parameter :: SUP_SINE = 1
integer, parameter :: SUP_CREN = 2
integer, parameter :: SUP_TRI  = 3
integer, parameter :: SUP_USER = 4
!
integer                                             :: sup_nwaves  = 1  ! Number of simultaneous waves
integer                                             :: sup_ngroups = 1  ! Number of pseudo atom groups
integer                                             :: sup_current = 1  ! Current wave number
character(len=PREC_STRING), dimension(:), allocatable:: sup_file   ! Input file
integer           , dimension(:)      , allocatable :: sup_group  ! Current group of pseudoatoms @ current wave
character(len=4)  , dimension(:,:,:,:), allocatable :: sup_atom   ! "Atom" name
!haracter(len=4)  , dimension(:,:)    , allocatable :: sup_repl   ! Replace site by this atom type
integer           , dimension(:,:)    , allocatable :: sup_irepl  ! Replace site by this many different atom types
integer           , dimension(:,:)    , allocatable :: sup_char   ! Character DISPLACEMENT / DENSITY
integer           , dimension(:,:)    , allocatable :: sup_func   ! Function  : sine, user
real(kind=PREC_DP), dimension(:,:)    , allocatable :: sup_func_p ! Function  : Additional parameter
logical           , dimension(:,:,:)  , allocatable :: sup_old    ! Density: old atom types
integer           , dimension(:,:,:)  , allocatable :: sup_new    ! Density: new atom types
real(kind=PREC_DP), dimension(:,:,:,:), allocatable :: sup_ampl   ! Displacement amplitude vector
real(kind=PREC_DP), dimension(:,  :,:), allocatable :: sup_phase  ! Phase in range [0,1]
real(kind=PREC_DP), dimension(:,:,:)  , allocatable :: sup_prob   ! Density propabilities [min,max]
!
real(kind=PREC_DP), dimension(:,:)    , allocatable :: sup_qvec   ! Satellite vector
!
end module superspace_mod
