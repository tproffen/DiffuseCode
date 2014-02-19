MODULE diffev_loop_mod
!
CONTAINS
!
SUBROUTINE diffev_loop
!                                                                       
USE diffev_mpi_mod
USE do_exit_mod
USE run_mpi_mod
USE doact_mod
USE errlist_mod 
USE learn_mod 
USE macro_mod 
USE prompt_mod 
!                                                                       
IMPLICIT none 
!
!*****7*****************************************************************
!                                                                       
!     Main loop for DIFFEV                                           
!                                                                       
!     This is the main program for DIFFEV. It sets up most              
!     variables and calls the loop interpreting the commands.           
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@mail.uni-erlangen.de)      
!                                                                       
!*****7*****************************************************************
!
!                                                                       
CHARACTER (LEN=1024)           :: line, zeile 
CHARACTER (LEN=4)              :: befehl 
LOGICAL                        :: lend    = .false.
INTEGER                        :: laenge, lp, lbef 
!
INTEGER, PARAMETER             :: master = 0 ! MPI ID of MASTER process
!                                                                       
with_mpi_error: IF ( ier_num == 0 ) THEN             ! No MPI error
   master_slave: IF ( run_mpi_myid == master ) THEN  ! MPI master or stand alone
!                                                                       
      CALL no_error 
!                                                                       
!------ This is the main loop: reading commands ..                      
!                                                                       
      main: DO while (.not.lend) 
         CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt) 
      no_cmd: IF (ier_num.eq.0.and.laenge.gt.0) then 
!                                                                       
!     - If not a comment continue                                       
!                                                                       
            IF (.not. (line (1:1) .eq.'#'.or.line (1:1) .eq.'!') ) then 
!                                                                       
!  - execute command                                                 
!                                                                       
               IF (line (1:3) .eq.'do '.OR.line (1:2) .eq.'if') then 
                  CALL do_loop (line, lend, laenge) 
               ELSE 
                  CALL mache_kdo (line, lend, laenge) 
               ENDIF 
            ENDIF 
         ENDIF no_cmd
!                                                                       
!     - Handle error message                                            
!                                                                       
         IF (ier_num.ne.0) then 
            CALL errlist 
            IF (ier_sta.ne.ER_S_LIVE) then 
               IF (lmakro.and.ier_sta.ne.ER_S_LIVE) then 
                  CALL macro_close 
                  prompt_status = PROMPT_ON 
               ENDIF 
               lblock = .false. 
               CALL no_error 
            ENDIF 
         ENDIF 
      ENDDO main
!                                                                       
!
   ELSEIF(run_mpi_active) THEN  master_slave
      CALL RUN_MPI_SLAVE  ! MPI slave, standalone never
   ELSE master_slave
      ier_num = -23       ! Mpi returned a slave ID, but MPI is not active !?!
      ier_typ = ER_APPL
   ENDIF master_slave
ENDIF with_mpi_error
!
CALL RUN_MPI_FINALIZE
!
IF ( lstandalone ) THEN
   CALL do_exit
ENDIF
!                                                                       
END SUBROUTINE diffev_loop
!*****7*****************************************************************
END MODULE diffev_loop_mod
