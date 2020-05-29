MODULE kuplot_loop_mod
!
CONTAINS
!
      SUBROUTINE kuplot_loop
!                                                                       
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal
      USE mpi_slave_mod
      USE do_if_mod
USE lib_errlist_func
USE lib_macro_func
      USE prompt_mod 
      USE sup_mod
!
      IMPLICIT none 
!*****7*****************************************************************
!       This is the universal plot program KUPLOT. It sets up most      
!     variables and calls the loop interpreting the commands.           
!*****7*****************************************************************
!                                                                       
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line, zeile 
      CHARACTER(4) befehl 
      LOGICAL lend 
      INTEGER lbef, lp, ll 
!
      EXTERNAL  :: kuplot_mache_kdo
!                                                                       
      lend = .false.
!                                                                       
!                                                                       
      main: DO WHILE (.NOT.lend) 
      CALL get_cmd (line, ll, befehl, lbef, zeile, lp, prompt) 
      ok: IF (ier_num.eq.0.and.ll.gt.0) THEN 
         IF (.NOT.(line==' ' .OR. line(1:1)=='#' .OR. line(1:1)=='!')) THEN
!                                                                       
!------ --- Execute command                                             
!                                                                       
         IF (befehl (1:3) .eq.'do '.or.befehl (1:2) .eq.'if') THEN 
            CALL do_loop (line, lend, ll) !, kuplot_mache_kdo) 
         ELSE 
            CALL kuplot_mache_kdo (line, lend, ll) !, previous) 
         ENDIF 
         ENDIF 
      ENDIF ok
!                                                                       
!     - Handle error message                                            
!                                                                       
      fehler: IF (ier_num.ne.0) THEN 
         IF( ier_num ==-9.and. ier_typ==ER_IO) THEN
            write(output_io, 8000)
            write(output_io, 9000)
            stop
         ENDIF
         IF(lstandalone) THEN
            CALL errlist 
            IF (ier_sta.ne.ER_S_LIVE) THEN 
               IF (lmakro) THEN 
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
            IF (ier_sta.ne.ER_S_LIVE) THEN 
               IF (lmakro .OR. lmakro_error) THEN 
                  IF(sprompt /= 'kuplot') THEN
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
      ENDIF  fehler
!
!        If loop was run from a non interactive remote and we
!        are no longer inside a makro, return after this command
         IF(.NOT. linteractive .AND. .NOT. lmakro) RETURN
      ENDDO main
!                                                                       
!------ END of PROGRAM                                                  
!                                                                       
      IF(lstandalone) THEN
         CALL kuplot_do_exit
      ENDIF
!                                                                       
8000 format(' ****EXIT**** Input error on normal read        ',        &
     &       '        ****',a1/)
9000 format(' ****EXIT**** KUPLOT  terminated by error status',        &
     &       '        ****',a1/)
!                                                                       
      END SUBROUTINE kuplot_loop
END MODULE kuplot_loop_mod
