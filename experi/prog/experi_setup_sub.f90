module experi_setup_sub_mod
!
contains
!
!*******************************************************************************
!
subroutine experi_set_sub
!-
! Sets the specific EXPERI interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!+
!
use set_sub_generic_mod
!
implicit none
!
interface
   subroutine experi_mache_kdo (line, lend, length)
!                                                                       
   character(len=*), intent(inout) :: line
   logical         , intent(  out) :: lend
   integer         , intent(inout) :: length
!
   end subroutine experi_mache_kdo
end interface
!
p_mache_kdo => experi_mache_kdo

end subroutine experi_set_sub
!
!*******************************************************************************
!
end module experi_setup_sub_mod
