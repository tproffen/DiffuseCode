module lib_forpython_mod
!-
!  Handle generic interaction with forpy_mod
!+
use forpy_mod
!
implicit none
!
logical :: forpy_isactive = .false.
!
private
public  forpy_active  !  logical function; TRUE if forpy_initialize was done
public  forpy_start   !  Subroutine      ; Initialize forpy_mod
public  forpy_finish  !  Subroutine      ; Finalize  forpy_mod
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
subroutine forpy_start(ier_num)
!-
!  Initialize forpy 
!-
!
implicit none
!
integer, intent(out) :: ier_num
!
if(.not. forpy_active()) then
   ier_num = forpy_initialize()
   forpy_isactive = .true.
endif
!
end subroutine forpy_start
!
!*******************************************************************************
!
subroutine set_module_path(module_path, ier_num)!  , python_script)
!-
!  Add path to python modules
!+
!
implicit none
!
character(len=*), intent(in) :: module_path
integer, intent(out) :: ier_num
!character(len=*), intent(in) :: module_name
!type(module_py) , intent(in) :: python_script
!
type(list)      :: paths_to_module    ! python script path
!type(object)    :: return_value       ! forpy return value
!
ier_num = get_sys_path(paths_to_module)
ier_num = paths_to_module%append(module_path)
!ier_num = import_py(python_script, module_name)
!
end subroutine set_module_path
!
!*******************************************************************************
!
subroutine forpy_finish()
!-
! Deactivate forpy_mod
!+
!
implicit none
!
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
