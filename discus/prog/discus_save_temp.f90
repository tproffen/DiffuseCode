module save_temp_mod
!-
!  Routines to save /restore the current structure
!+
contains
!
!*******************************************************************************
!
subroutine save_temp(save_file)
!-
!  Save the current structure into the internal file "save_file"
!+
!
use crystal_mod
use discus_allocate_appl_mod
use prop_para_func
use save_menu
use discus_save_mod
!
use precision_mod
!
implicit none
!
character(len=*), intent(inout) :: save_file
!
character(len=PREC_STRING) :: int_save_file
character(len=PREC_STRING) :: line
integer                    :: length
!
call save_store_setting             ! Backup user "save" setting
call save_default_setting           ! Default to full saving
line       = 'ignore, all'          ! Ignore all properties
length     = 11
call property_select(line, length, sav_sel_prop)
line       = 'ignore, all'          ! Ignore all properties for global as well
length     = 11
call property_select(line, length,  cr_sel_prop)
!
if(save_file(1:8)/='internal') then
  int_save_file = 'internal.' // save_file(1:len_trim(save_file))
else
  int_save_file = save_file(1:len_trim(save_file))
endif
!
call alloc_unitcell(cr_ncatoms)
call save_internal(int_save_file)        !     thus this file name is unique
!
end subroutine save_temp
!
!*******************************************************************************
!
subroutine restore_temp(save_file)
!-
! Read the saved backup structure
!+
!
use read_internal_mod
use save_menu
!
!use errlist_mod
use lib_errlist_func
use precision_mod
!
implicit none
!
character(len=*), intent(inout) :: save_file
!
character(len=PREC_STRING) :: int_save_file
integer, parameter :: MAXMASK=4
logical, dimension(0:MAXMASK) :: uni_mask
!
uni_mask(0)   = .false.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
!
CALL errlist_save                   ! Keep error status 

CALL save_restore_setting
CALL no_error
!
if(save_file(1:8)/='internal') then
  int_save_file = 'internal.' // save_file(1:len_trim(save_file))
else
  int_save_file = save_file(1:len_trim(save_file))
endif
!
CALL readstru_internal(MAXMASK, int_save_file, uni_mask)   ! Read  core file
CALL errlist_restore                ! Restore error status
!
end subroutine restore_temp
!
!*******************************************************************************
!
end module save_temp_mod
