module kuplot_exit_mod
!
contains
!
!*****7**************************************************************** 
!
SUBROUTINE kuplot_do_exit 
!                                                                       
!       Things to do when KUPLOT exits finally 
!                                                                       
USE kuplot_config 
USE kuplot_mod 
!
USE errlist_mod 
USE exit_mod
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
!------ call system wide exit routine                                   
!                                                                       
CALL exit_all 
!                                                                       
!------ close PGPLOT devices                                            
!                                                                       
CALL PGEND 
!                                                                       
END SUBROUTINE kuplot_do_exit                        
!
!*****7**************************************************************** 
!
SUBROUTINE kuplot_exit 
!                                                                       
!       Things to do when KUPLOT exits finally 
!                                                                       
USE kuplot_mod
use kuplot_save_mod
!
use errlist_mod
use precision_mod
use lib_global_flags_mod
use berechne_mod
!                                                                       
IMPLICIT none 
character(len=4)           :: form
character(len=PREC_STRING) :: string
integer                    :: kid, indiv
integer                    :: ik
integer                    :: ianz
integer, parameter         :: MAXW = 4
integer                    :: length       ! String length
logical                    :: l_m999
real(kind=PREC_DP), dimension(MAXW) :: werte
integer :: ier_cmd
integer :: exit_msg
character(len=PREC_STRING) :: message
!
if(iz==1) return   ! No data immediate return
!
if(lib_global_flags(1) == 1) then   ! Save last data set
   string = 'mkdir -p DISCUS_SUITE_DERIVATIVES'
   call execute_command_line(string, CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg, wait=.true.)
   string='REF_KID'
   length = 7
   kid = nint(berechne(string, length))
   string='REF_INDIV'
   length = 9
   indiv = nint(berechne(string, length))
   string = ' '
   write(string,'(a,i4.4,a)') 'DISCUS_SUITE_DERIVATIVES/data.', kid! , '.', indiv
   length = 35
   ik     = iz - 1
   form   = 'H5'
   ianz   = 0
   werte  = 0.0_PREC_DP
   l_m999 = .FALSE.
   call check_form (ik, form, ianz, werte, MAXW)
   call do_save(ik, string, form, ianz, werte, MAXW, l_m999)
endif
!                                                                       
END SUBROUTINE kuplot_exit
!
SUBROUTINE kuplot_sigint
!
!     Handle KUPLOT specific part of a CTRL-C interrupt.
!     CALLED within standalone KUPLOT only
!     This subroutine calls all KUPLOT specific emergency handlers
!     which can also be called from the SUITE
!
USE exit_mod
IMPLICIT NONE
CHARACTER(LEN=1) :: dummy
WRITE(*,*)
WRITE(*,*) ' EMERGENCY Shutdown with USER CTRL-C Interrupt'
!
CALL exit_all
!
WRITE(*,*) 
WRITE(*,*) ' KUPLOT closed by User Request CTRL-C '
WRITE(*,*) ' For final close down hit ENTER key'
READ(*,'(a)') dummy
STOP        ! Terminate program
!
END SUBROUTINE kuplot_sigint
!
end module kuplot_exit_mod
