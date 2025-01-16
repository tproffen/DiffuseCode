module experi_loop_mod
!-
! Main EXPERI Loop module, runs the loop to execute commands
!-
contains
!
!*******************************************************************************
!
subroutine experi_loop
!
USE class_macro_internal    ! Macro tree routines
use errlist_mod             ! Variables ier_*
use doact_mod               ! spropmt, loop active
use do_if_mod               ! perform do loop, if condition
use lib_errlist_func        ! subroutine errlist
use lib_macro_func          ! Macro routines
use prompt_mod              ! program prompt, io variables
use precision_mod           ! Variable precision
use sup_mod                 ! Generic support, get_cmd
!
implicit none
!
character(len=PREC_STRING) :: line    ! User command line
character(len=PREC_STRING) :: string  ! User argument line
character(len=4)           :: command ! User types this command
integer :: length                     ! Command line length
integer :: lcmd                       ! Command length
integer :: lp                         ! argument line length
logical :: lend                       ! end flag return to SUITE
!
lend = .FALSE.
!                                                                       
!------ This is the main loop: reading commands ..                      
!
loop_main: do
   if(lend) exit loop_main
!
   call get_cmd(line, length, command, lcmd, string, lp, prompt)
   if(ier_num==0 .and. length==0) cycle loop_main             ! No/empty command
   if((line(1:1)=='#' .or. line(1:1)=='!')) cycle loop_main   ! Comment, ignore
!
!  Handle error messages
!
   cond_err: if(ier_num/=0) then
      if(ier_num ==-9.and. ier_typ==ER_IO)  then
         write(output_io, 8000)
         write(output_io, 9000)
         stop
      endif
!
!     Regular error!
!
      call errlist
      cond_sta: if(ier_sta/=ER_S_LIVE) then
         cond_macro: if(lmakro .OR. lmakro_error) then
            if(sprompt /= 'experi') then
               ier_num = -9
               ier_typ = ER_COMM
               exit loop_main
            endif
         else cond_macro
            if(lmacro_close) then
               call macro_close(-1)
               lmakro_error = .FALSE.
               PROMPT_STATUS = PROMPT_ON
               sprompt = ' '
            endif
         endif cond_macro
      endif cond_sta
   endif cond_err
!
   if(line(1:3)=='do ' .or. line(1:2) == 'if') then
      call do_loop(line, lend, length)
   else
      call experi_mache_kdo(line, lend, length)
   endif
!
enddo loop_main
!
!
8000 FORMAT(' ****EXIT**** Input error on normal read        ',        &
     &       '        ****',a1/)
9000 FORMAT(' ****EXIT**** EXPERI  terminated by error status',        &
     &       '        ****',a1/)
!
end subroutine experi_loop
!
!*******************************************************************************
!
end module experi_loop_mod
