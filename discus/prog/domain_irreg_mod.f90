module domain_irreg_mod
!-
!  Variables for the irreguarly shaped domeins
!+
!*******************************************************************************
!
use precision_mod
!
integer                                         :: ngroup  ! Group size
integer                                         :: mgroup  ! Group size short
integer,            dimension(:),   allocatable :: short   ! List of atoms in the current domain
integer                                         :: icent   ! Central atom of group
integer                                         :: icent_cr! Central atom in crystal
real(kind=PREC_DP), dimension(3)                :: origin  ! Origin of this domain
integer           , dimension(3)                :: orig_off  ! Origin of this domain
real(kind=PREC_DP), dimension(:,:), allocatable :: offs    ! Offsets from group origin
real(kind=PREC_DP), dimension(:,:), allocatable :: coor    ! Coordinates relative to group origin
real(kind=PREC_DP), dimension(3,2)              :: maxdim  ! Dimensions of this domain
!
end module domain_irreg_mod
