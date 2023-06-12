module variable_test
!
use variable_mod
!
implicit none
!
contains
!
subroutine variable_exist(c_temp,l_temp, c_type, l_exist, l_type, var_no)
!
character (LEN=*), intent(in)  :: c_temp    ! Name to be tested
integer          , intent(in)  :: l_temp    ! length of input name
integer          , INTENT(in)  :: c_type    ! Type of input variable
logical          , intent(out) :: l_exist   ! True if exists
logical          , intent(out) :: l_type    ! True if correct type
integer          , intent(out) :: var_no    ! variable number if exists
!
integer :: i
!
l_exist = .false.
l_type  = .false.
var_no  = 0
!
search: do i=1,var_num
   if( c_temp(1:l_temp) == var_name(i)) then
      l_exist = .true.
      var_no  = i
      if( c_type == var_type(i)) then
         l_type = .true.
      endif
      exit search
   endif
enddo search
!
end subroutine variable_exist
!
end module variable_test
