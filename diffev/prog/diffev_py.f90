MODULE diffev_interface
!
PUBLIC diffev
PUBLIC diffev_run
!
CONTAINS
!
SUBROUTINE interactive ()
!
!  Generic interface routine to start an interactive diffev session
!  from the host system
!
USE prompt_mod
USE diffev_setup_mod
USE diffev_loop_mod
!
IMPLICIT none 
!
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL diffev_setup
ENDIF
lstandalone = .false.
CALL diffev_loop
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
! If an error occurs inside a sub-menu, diffev will start an 
! interactive session at this point
! Commands that branch into sub-menus cause an interactive section.
! 
! 
USE diffev_setup_mod
USE run_mpi_mod
USE diffev_mpi_mod
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
INTEGER, PARAMETER   :: master = 0 ! MPI ID of MASTER process
!
INTEGER              :: len_str
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL diffev_setup
ENDIF
master_slave: IF ( run_mpi_myid == master ) THEN ! MPI master or standalone
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
ELSEIF ( run_mpi_active ) THEN master_slave
   CALL run_mpi_slave     ! MPI slave, standalone never
   CALL run_mpi_finalize  ! Always end MPI, if slave is finished
ELSE master_slave
   ier_num = -23          ! MPI returned a slave, BUT MPI is not active ?
ENDIF master_slave
!
IF ( lend ) THEN
   IF( run_mpi_myid == master .and. run_mpi_active ) THEN
CONTINUE
      CALL run_mpi_finalize
   ENDIF
ENDIF
!
END SUBROUTINE command
!
END MODULE diffev_interface
