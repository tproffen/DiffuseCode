MODULE diffev_loop_mod
!
CONTAINS
!
SUBROUTINE diffev_loop
!                                                                       
USE diffev_mpi_mod
!use diffev_mache_kdo_mod
USE create_trial_mod
USE gen_mpi_mod
USE doact_mod
USE errlist_mod 
USE learn_mod 
USE class_macro_internal 
USE mpi_slave_mod
USE diffev_random
USE do_if_mod
USE lib_errlist_func
USE lib_macro_func
USE precision_mod
USE prompt_mod 
USE set_sub_generic_mod
USE sup_mod
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
!     Authors : R.B. Neder  (reinhard.neder@fau.de)      
!                                                                       
!*****7*****************************************************************
!
!                                                                       
CHARACTER (LEN=PREC_STRING)    :: line, zeile 
CHARACTER (LEN=4)              :: befehl 
LOGICAL                        :: lend   = .false.
INTEGER                        :: laenge, lp, lbef 
!
INTEGER, PARAMETER             :: master = 0 ! MPI ID of MASTER process
!
lend = .false.                                       ! Always initialize the loop
with_mpi_error: IF ( ier_num == 0 ) THEN             ! No MPI error
   master_slave: IF ( gen_mpi_myid == master ) THEN  ! MPI master or stand alone
!                                                                       
      CALL no_error 
!     INQUIRE(FILE='GENERATION', EXIST=lexist)
!     IF(lexist) THEN
!        CALL read_genfile
!     ENDIF
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
                  CALL diffev_mache_kdo (line, lend, laenge) 
               ENDIF 
            ENDIF 
         ENDIF no_cmd
!                                                                       
!     - Handle error message                                            
!                                                                       
         IF (ier_num.ne.0) then 
            IF( ier_num ==-9.and. ier_typ==ER_IO) THEN
               write(output_io, 8000)
               write(output_io, 9000)
               stop
            ENDIF
               IF(mpi_active .AND. ier_sta == ER_S_EXIT) THEN  ! Error while MPI is on
                  ier_sta = ER_S_LIVE              ! Fake Error status to prevent stop
                  CALL errlist                     ! but get error message
                  ier_sta = ER_S_EXIT              ! Signal EXIT back to SUITE
                  ier_num = -9                     ! Signal error condition to SUITE
                  ier_typ = ER_COMM
                  EXIT main                        ! Now terminate program gracefully
               ENDIF
               CALL errlist
               IF (ier_sta.ne.ER_S_LIVE) then 
                  IF (lmakro .OR. lmakro_error) then 
                     IF(sprompt /= 'diffev') THEN
                        ier_num = -9
                        ier_typ = ER_COMM
                        EXIT main
                     ELSE
                        IF(lmacro_close) THEN
                           CALL macro_close(-1)
                           lmakro_error = .FALSE.
                           PROMPT_STATUS = PROMPT_ON
                           sprompt = ' '
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
         ENDIF 
!
!        If loop was run from a non interactive remote and we
!        are no longer inside a makro, return after this command
         IF(.NOT. linteractive .AND. .NOT. lmakro) RETURN
      ENDDO main
!                                                                       
      IF(lend) THEN
         CALL diffev_best_macro
      ENDIF
!
   ELSEIF(gen_mpi_active) THEN  master_slave
      CALL RUN_MPI_SLAVE  ! MPI slave
   ELSE master_slave
      ier_num = -23       ! Mpi returned a slave ID, but MPI is not active !?!
      ier_typ = ER_APPL
   ENDIF master_slave
ENDIF with_mpi_error
!
8000 format(' ****EXIT**** Input error on normal read        ',        &
     &       '        ****',a1/)
9000 format(' ****EXIT**** DIFFEV  terminated by error status',        &
     &       '        ****',a1/)
!
END SUBROUTINE diffev_loop
!*****7*****************************************************************
END MODULE diffev_loop_mod
