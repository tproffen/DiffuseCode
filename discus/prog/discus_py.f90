MODULE discus_interface
!
PUBLIC discus
PUBLIC discus_run
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
USE discus_loop_mod
!
IMPLICIT none 
!
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL discus_setup
ENDIF
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
USE errlist_mod
USE macro_mod
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN   ) :: incomming
INTEGER         , INTENT(  OUT) :: ier_status
!
CHARACTER(LEN=1024)  :: line
CHARACTER(LEN=1024)  :: zeile
CHARACTER(LEN=   4)  :: befehl
INTEGER              :: laenge
INTEGER              :: lbef
INTEGER              :: lp
LOGICAL              :: lend
!
INTEGER              :: len_str
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL discus_setup
ENDIF
lend = .false.
!
laenge = len_str(incomming)     ! 
IF ( laenge > LEN(line)) THEN    ! Excessively long string, refuse
   ier_status = -1
   RETURN
ENDIF
line   = incomming              ! Make a local working copy
!
CALL mache_kdo (line, lend, laenge)  ! Execute initial command
!
IF( ier_num /= 0 ) THEN         ! Handle error messages
   CALL errlist
   ier_status = -1
   CALL macro_close
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
               CALL mache_kdo (line, lend, laenge) 
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
END MODULE discus_interface
