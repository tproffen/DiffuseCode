module lib_forpython_mod
!-
!  Handle generic interaction with forpy_mod
!+
use errlist_mod
use forpy_mod
!
implicit none
!
logical :: forpy_isactive = .false.
!
private
public  forpy_active  !  logical function; TRUE if forpy_initialize was done
public  forpy_start   !  Subroutine      ; Initialize forpy_mod
!
contains
!
!*******************************************************************************
!
logical function forpy_active()
   forpy_active = forpy_isactive
end function forpy_active
!
!*******************************************************************************
!
subroutine forpy_start()
!-
!  Initialize forpy 
!-
if(.not. forpy_active()) then
   ier_num = forpy_initialize()
   forpy_isactive = .true.
endif
!
end subroutine forpy_start
!
!*******************************************************************************
!
subroutine forpy_finish()
!-
! Deactivate forpy_mod
!+
if(forpy_active()) then
   call forpy_finalize
   forpy_isactive = .false.
endif
!
end subroutine forpy_finish
!
!*******************************************************************************
!
end module lib_forpython_mod
