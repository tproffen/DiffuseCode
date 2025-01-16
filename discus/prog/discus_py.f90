MODULE discus
!
PRIVATE
PUBLIC  get, get_r, get_i, send,send_i, send_r
PUBLIC  interactive
PUBLIC  command
!
INTERFACE get
   MODULE PROCEDURE get_i, get_r
END INTERFACE get
INTERFACE send
   MODULE PROCEDURE send_i, send_r
END INTERFACE send
!
CONTAINS
!
SUBROUTINE interactive ()
!
!  Generic interface routine to start an interactive discus session
!  from the host system
!
USE prompt_mod
USE discus_setup_mod
USE discus_setup_sub_mod
USE discus_loop_mod
USE set_sub_generic_mod
!
IMPLICIT none 
!
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL discus_setup
ENDIF
CALL discus_set_sub
lstandalone = .false.
CALL discus_loop
!                                                                       
END SUBROUTINE interactive                            
!
!
!
SUBROUTINE command (incomming, ier_status)
!
! Interface routine to execute the command on the incomming line
! Control is returned to the host system
! Currently commands 'do' 'enddo' and 'if', 'elseif', 'else'
! 'endif' cannot be used. 
! A command '@macro.mac' is fine, and the macro may even contain
! loops, and if blocks. 
! If an error occurs inside a sub-menu, discus will start an 
! interactive session at this point
! Commands that branch into sub-menus cause an interactive section.
! 
! 
USE discus_setup_mod
USE discus_setup_sub_mod
!use discus_mache_kdo_mod
USE errlist_mod
USE class_macro_internal
USE lib_errlist_func
USE lib_length
USE lib_macro_func
USE precision_mod
USE prompt_mod
USE set_sub_generic_mod
USE sup_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN   ) :: incomming
INTEGER         , INTENT(  OUT) :: ier_status
!
CHARACTER(LEN=MAX(PREC_STRING.LEN(incomming)))  :: line
CHARACTER(LEN=MAX(PREC_STRING.LEN(incomming)))  :: zeile
CHARACTER(LEN=   4)  :: befehl
INTEGER              :: laenge
INTEGER              :: lbef
INTEGER              :: lp
LOGICAL              :: lend
!
!
!EXTERNAL             :: discus_mache_kdo
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL discus_setup
ENDIF
CALL diffev_set_sub
lend = .false.
!
laenge = len_str(incomming)     ! 
IF ( laenge > LEN(line)) THEN    ! Excessively long string, refuse
   ier_status = -1
   RETURN
ENDIF
line   = incomming              ! Make a local working copy
!
CALL discus_mache_kdo (line, lend, laenge)  ! Execute initial command
!
IF( ier_num /= 0 ) THEN         ! Handle error messages
   CALL errlist
   ier_status = -1
   CALL macro_close(-1)
   CALL no_error
ELSE
   main: DO WHILE( lmakro )     ! Initial command was a macro, run this macro
      CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt)
      ok: IF (ier_num.eq.0.and.laenge.gt.0) then 
!                                                                       
!     - If not a comment continue                                       
!                                                                       
         IF (.not. (line (1:1) .eq.'#'.or.line (1:1) .eq.'!') ) then 
!                                                                       
!     - execute command                                                 
!                                                                       
            IF (line (1:3) .eq.'do '.OR.line (1:2) .eq.'if') then 
               CALL do_loop (line, lend, laenge) 
            ELSE 
               CALL discus_mache_kdo (line, lend, laenge) 
            ENDIF 
         ENDIF 
      ENDIF ok
!                                                                       
!     - Handle error message                                            
!                                                                       
      IF( ier_num /= 0 ) THEN
         CALL errlist
         ier_status = -1
         CALL macro_close
         CALL no_error
         EXIT main
      ELSE
         ier_status = 0
      ENDIF
   ENDDO main
ENDIF
!
END SUBROUTINE command
!
!  INCLUDE the generic send and get routines from lib_f90
!  These allow to send/get sections of i[] and r[].
!  As these are identical to all programs, the source 
!  code is in lib_f90. As I want to have these routines 
!  to be part of this module, its easiest to include
!  the source code instead of adding another file to 
!  the f2py command
!
INCLUDE 'send_get.f90'
!
END MODULE discus
