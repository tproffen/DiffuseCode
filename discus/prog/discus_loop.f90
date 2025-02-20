MODULE discus_loop_mod
!
CONTAINS
!
!*****7*****************************************************************
!                                                                       
recursive SUBROUTINE discus_loop
!                                                                       
USE discus_exit_mod
!use discus_mache_kdo_mod
USE doact_mod 
USE errlist_mod 
USE learn_mod 
USE class_macro_internal 
USE mpi_slave_mod
USE do_if_mod
USE lib_errlist_func
USE lib_macro_func
USE precision_mod
USE prompt_mod 
USE sup_mod
!
IMPLICIT none 
!
!*****7*****************************************************************
!                                                                       
!     Main program for DISCUS                                           
!                                                                       
!     This is the main program for DISCUS. It sets up most              
!     variables and calls the loop interpreting the commands.           
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@fau.de)      
!               Th. Proffen (tproffen@ornl.gov)                         
!                                                                       
!*****7*****************************************************************
!                                                                       
CHARACTER (LEN=PREC_STRING) :: line   = ' '
CHARACTER (LEN=PREC_STRING) :: zeile  = ' '
CHARACTER (LEN=4)    :: befehl = ' '
LOGICAL              :: lend 
INTEGER              :: laenge = 1
INTEGER              :: lp     = 1
INTEGER              :: lbef   = 1
!                                                                       
!EXTERNAL             :: discus_mache_kdo    ! Declare DISCUS copy of mache_kdo
!                                                                       
lend    = .false. 
!                                                                       
!------ This is the main loop: reading commands ..                      
!                                                                       
main: DO while (.not.lend) 
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
            CALL do_loop (line, lend, laenge) !, discus_mache_kdo) 
         ELSE 
            CALL discus_mache_kdo (line, lend, laenge) 
         ENDIF 
      ENDIF 
   ENDIF ok
!                                                                       
!     - Handle error message                                            
!                                                                       
fehler: IF (ier_num.ne.0) then 
      IF( ier_num ==-9.and. ier_typ==ER_IO) THEN  
         WRITE(output_io, 8000)
         WRITE(output_io, 9000)
         stop
      ENDIF
      IF(lstandalone) THEN
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro.AND. lmacro_close) then 
               CALL macro_close(-1)
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
         IF (ier_sta /= ER_S_LIVE) THEN 
            IF (lmakro .OR.  lmakro_error) THEN 
               IF(sprompt /= 'discus') THEN
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
   ENDIF fehler
!
!  If loop was run from a non interactive remote and we
!  are no longer inside a makro, return after this command
   IF(.NOT. linteractive .AND. .NOT. lmakro) RETURN
ENDDO main
!                                                                       
IF ( lstandalone ) THEN
   CALL discus_do_exit 
ENDIF
!
8000 format(' ****EXIT**** Input error on normal read        ',        &
     &       '        ****',a1/)
9000 format(' ****EXIT**** Program terminated by error status',        &
     &       '        ****',a1/)
!                                                                       
END SUBROUTINE discus_loop
!                                                                       
!*****7*****************************************************************
!
END MODULE discus_loop_mod
