module lib_config
!-
!  Variables and routines for general config
!+
!
use precision_mod
!
private
public generic_read_config
public generic_write_config
public generic_get_interval
public generic_set_interval
!
character(len=PREC_STRING) :: config_name = 'DISCUS.config'
integer :: os_update_interval = -1
!
!*******************************************************************************
!
contains
!
!*******************************************************************************
!
subroutine generic_make_dir(discus_dir, discus_dir_l)
!-
!  Test and create $HOME/.DISCUS
!+
!
use errlist_mod
use lib_errlist_func
use precision_mod
use sys_compiler
!
implicit none
!
character(len=*), intent(in) :: discus_dir
integer         , intent(in) :: discus_dir_l
!
integer, parameter :: MAXW = 1
!
character(len=PREC_STRING)               :: line
character(len=PREC_STRING), dimension(1) :: cpara
character(len=PREC_STRING)               :: message
integer                   , dimension(1) :: lpara
integer :: ios
integer :: ier_cmd
logical :: lexist
!
call no_error
cpara(1) = discus_dir
lpara(1) = discus_dir_l
lexist = .false.
call sys_inquire_directory(MAXW, cpara, lpara, lexist)
if(.not.lexist) then                    ! Directory did not exist, create
   line = 'mkdir -p ' // discus_dir(1:discus_dir_l)
   call execute_command_line(line, wait=.true.,                                  &
        cmdstat=ier_cmd, cmdmsg=message, exitstat=ios)
   if(ios/=0) then
      ier_num    = -2
      ier_typ    = ER_IO
      ier_msg(1) = 'Check access rights to your HOME directory'
   endif
endif
!
end subroutine generic_make_dir
!
!*******************************************************************************
!
subroutine generic_read_config(discus_dir, discus_dir_l)
!-
!  Read the config file DISCUS.config  in $HOME/.DISCUS
!+
!
use blanks_mod
use errlist_mod
use precision_mod
use sys_compiler
!
implicit none
!
character(len=*), intent(in) :: discus_dir
integer         , intent(in) :: discus_dir_l
!
integer, parameter :: IRD=37
!
character(len=PREC_STRING)               :: string
character(len=PREC_STRING)               :: config_file
integer :: ios
integer :: length
logical :: lexist
!
call generic_make_dir(discus_dir, discus_dir_l)
if(ier_num/=0)  return
!
config_file = discus_dir(1:discus_dir_l) // config_name(1:len_trim(config_name))
inquire(file=config_file, exist=lexist)
!
if(lexist) then                      ! OLd config file exists, read
   open(unit=ird, file=config_file, status='old', iostat=ios) 
   loop_read: do
      read(ird,'(a)',iostat=ios) string
      if(ios/=0) exit loop_read
      length = len_trim(string)
      call rem_bl(string,length)
      cond_os_upd:if(string == '[OsUpdates]') then   ! Operating system info
         loop_os: do
            read(ird,'(a)',iostat=ios) string
            if(ios/=0) exit loop_read
            if(string(1:1)=='[') exit cond_os_upd
            if(string(1:9) == 'Interval=') then
               read(string(10:13),'(i4)',iostat=ios) os_update_interval

            endif
            if(ios/=0) exit loop_read
         enddo loop_os
      endif cond_os_upd
   enddo loop_read
   close(unit=IRD)
else                                 ! no config file create a new
   call generic_write_config(discus_dir, discus_dir_l)
   if(ier_num/=0)  return
endif
!
end subroutine generic_read_config
!
!*******************************************************************************
!
subroutine generic_write_config(discus_dir, discus_dir_l)
!-
!  Write the config file DISCUS.config  in $HOME/.DISCUS
!+
!
use errlist_mod
use precision_mod
use sys_compiler
!
implicit none
!
character(len=*), intent(in) :: discus_dir
integer         , intent(in) :: discus_dir_l
!
integer, parameter :: IWR=38
!integer, parameter :: MAXW = 1
!
character(len=PREC_STRING)               :: config_file
!character(len=PREC_STRING), dimension(1) :: cpara
!integer                   , dimension(1) :: lpara
integer :: ios
!
call generic_make_dir(discus_dir, discus_dir_l)
if(ier_num/=0)  return
!
config_file = discus_dir(1:discus_dir_l) // config_name(1:len_trim(config_name))
!
open(unit=IWR, file=config_file, status='unknown', iostat=ios) 
!
write(IWR,'(a)') '[OsUpdates]'
if(os_update_interval>0) then
   write(IWR,'(a,i4.4)') 'Interval=', os_update_interval
else
   write(IWR,'(a)') 'Interval=-001'
endif
close(unit=IWR)
!
end subroutine generic_write_config
!
!*******************************************************************************
!
integer function generic_get_interval()
!
! returns the os_update_interval
!+
generic_get_interval = os_update_interval
!
end function generic_get_interval
!
!*******************************************************************************
!
subroutine generic_set_interval(i)
!
! returns the os_update_interval
!+
implicit none
integer, intent(in) :: i
!
os_update_interval = i
!
end subroutine generic_set_interval
!
!*******************************************************************************
!
end module lib_config
