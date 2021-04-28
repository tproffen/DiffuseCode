module privat_mod
!-
!  A module in which private calculations can be performed
!  The actual file "private.f90" is maintained as this template only. 
!  Your private modifications may get lost if you keep "private.f90"
!  within DiffuseCode/discus/prog only.
!
!  Compare this module to any of the build in ones to create your own 
!  stuff.
!+
contains
!
!*******************************************************************************
!
subroutine do_private(line)
!-
!  Do a private structure modification.
!  The list mod used modules serves as example only.
!+
use crystal_mod
use atom_env_mod
use chem_mod
use do_find_mod
!
use precision_mod
use param_mod
!
implicit none
!
character(len=*), intent(inout) :: line
!
! With further use of the modules, you can perform private calculations / 
! modifications of the crystal.
!
end subroutine do_private
!
!*******************************************************************************
!
end module privat_mod
