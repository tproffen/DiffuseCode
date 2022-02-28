MODULE refine_loop_mod
!
CONTAINS
!
SUBROUTINE refine_loop
!
!use refine_mache_kdo_mod
USE diffev_mpi_mod
USE gen_mpi_mod
USE mpi_slave_mod
!
USE class_macro_internal 
USE doact_mod
USE do_if_mod
USE errlist_mod 
USE learn_mod 
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
!     Main loop for REFINE                                           
!                                                                       
!     This is the main program for REFINE. It sets up most              
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
                  CALL refine_mache_kdo (line, lend, laenge) 
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
            IF(lstandalone) THEN
               CALL errlist 
               IF (ier_sta.ne.ER_S_LIVE) then 
                  IF (lmakro.and.ier_sta.ne.ER_S_LIVE.AND.lmacro_close) then 
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
                  ENDIF 
                  lblock = .false. 
                  CALL no_error 
               ENDIF 
            ELSE
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
                     IF(sprompt /= 'refine') THEN
                        ier_num = -9
                        ier_typ = ER_COMM
                        EXIT main
                     ELSE
                        IF(lmacro_close) THEN
                           CALL macro_close
                           lmakro_error = .FALSE.
                           PROMPT_STATUS = PROMPT_ON
                           sprompt = ' '
                        ENDIF 
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
   ELSEIF(gen_mpi_active) THEN  master_slave
      CALL RUN_MPI_SLAVE  ! MPI slave, standalone never
   ELSE master_slave
      ier_num = -23       ! Mpi returned a slave ID, but MPI is not active !?!
      ier_typ = ER_APPL
   ENDIF master_slave
ENDIF with_mpi_error
!
8000 FORMAT(' ****EXIT**** Input error on normal read        ',        &
     &       '        ****',a1/)
9000 FORMAT(' ****EXIT**** REFINE  terminated by error status',        &
     &       '        ****',a1/)
!
END SUBROUTINE refine_loop
!
!*****7*****************************************************************
!
END MODULE refine_loop_mod
