MODULE kdo_all_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE kdo_all (bef, lbef, zei, lc) 
!                                                                       
!     This part contains all the program independent commands.          
!                                                                       
!     USE calc_expr_mod
!
USE appl_env_mod
USE arrays_mod
USE ber_params_mod
USE define_variable_mod
USE do_eval_mod
USE do_show_mod
USE do_set_mod
USE do_wait_mod
USE errlist_mod 
USE fput_mod
USE get_params_mod
USE learn_mod 
USE lib_learn
USE lib_do_operating_mod
USE lib_echo
USE lib_help
USE lib_macro_func
USE precision_mod
USE prompt_mod 
USE sockets_mod
USE random_state_mod
USE str_comp_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL, PARAMETER :: NO_DIFFEV = .FALSE.
      INTEGER, PARAMETER :: MAXPAR = 2
!                                                                       
      CHARACTER (LEN=*) , INTENT(IN)    :: bef 
      INTEGER           , INTENT(IN)    :: lbef 
      CHARACTER (LEN=*) , INTENT(INOUT) :: zei 
      INTEGER           , INTENT(INOUT) :: lc
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zei))) :: command 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zei))) :: cpara (maxpar) 
      INTEGER lpara (maxpar) 
      INTEGER ianz 
      REAL(KIND=PREC_DP):: werte (maxpar) 
!                                                                       
!                                                                       
!     change working directory                                          
!                                                                       
      IF (str_comp (bef, 'cd', 2, lbef, 2) ) THEN 
         CALL do_chdir (zei, lc, .true.) 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
      ELSEIF (str_comp (bef, 'continue', 2, lbef, 8) ) THEN 
         CALL macro_continue (zei, lc) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
      ELSEIF (str_comp (bef, 'echo', 2, lbef, 4) ) THEN 
         CALL echo (zei, lc) 
!                                                                       
!------ Echo a string to *
!                                                                       
      ELSEIF (str_comp (bef, 'flush', 3, lbef, 5) ) THEN 
         CALL do_flush(zei, lc) 
!                                                                       
!     Evaluate an expression, just for interactive check 'eval'         
!                                                                       
      ELSEIF (str_comp (bef, 'eval', 2, lbef, 4) ) THEN 
         CALL do_eval (zei, lc, .TRUE.) 
!                                                                       
!------ IO commands                                                     
!                                                                       
      ELSEIF (str_comp (bef, 'bye', 2, lbef, 3) .and.lconn) THEN 
         CALL socket_close (s_conid) 
         lconn = .false. 
!                                                                       
      ELSEIF (str_comp (bef, 'fopen', 2, lbef, 5) ) THEN 
         CALL do_fopen (zei, lc) 
      ELSEIF (str_comp (bef, 'fclose', 2, lbef, 6) ) THEN 
         CALL do_fclose (zei, lc) 
      ELSEIF (str_comp (bef, 'fend', 2, lbef, 4) ) THEN 
         CALL do_fend (zei, lc) 
      ELSEIF (str_comp (bef, 'fexist', 2, lbef, 6) ) THEN 
         CALL do_fexist (zei, lc,.TRUE.) 
      ELSEIF (str_comp (bef, 'fget', 2, lbef, 4) ) THEN 
         CALL do_fget (zei, lc) 
      ELSEIF (str_comp (bef, 'fput', 2, lbef, 4) ) THEN 
         CALL do_fput (zei, lc) 
      ELSEIF (str_comp (bef, 'fsub', 2, lbef, 4) ) THEN 
         CALL do_fgetsub (zei, lc) 
      ELSEIF (str_comp (bef, 'fformat', 2, lbef, 7) ) THEN 
         CALL do_fformat (zei, lc) 
!                                                                       
!------ Online help                                                     
!                                                                       
      ELSEIF (linteractive .AND. str_comp (bef, 'help', 2, lbef, 4) .OR.   &
                                 str_comp (bef, '?   ', 1, lbef, 4) ) THEN                                            
         IF (lc.gt.0) THEN 
            lc = lc + 7 
            WRITE (command, '(a6,1x,a)') pname, zei (1:lc) 
         ELSE 
            lc = 7 
            command = pname 
         ENDIF 
         CALL do_hel (command, lc) 
!                                                                       
!------ Online manual                                                     
!
      ELSEIF (linteractive .AND. str_comp (bef, 'manual', 2, lbef, 6) ) THEN
        IF (lc.gt.0) THEN
           WRITE (command, '(a)') zei (1:lc)
        ELSE
           WRITE (command, '(a,a)') 'section:',pname
           lc = LEN_TRIM(pname)
        ENDIF
        CALL do_manual(command, lc)
!                                                                       
!------ start learning a macro 'learn'                                  
!                                                                       
      ELSEIF (str_comp (bef, 'lear', 3, lbef, 4) ) THEN 
         CALL start_learn (zei, lc) 
!                                                                       
!------ end learning a macro 'lend'                                     
!                                                                       
      ELSEIF (str_comp (bef, 'lend', 3, lbef, 4) ) THEN 
         CALL ende_learn 
!                                                                       
!------ Multiply two user variable arrays 'matmul'                      
!                                                                       
      ELSEIF (str_comp (bef, 'matmul', 4, lbef, 6) ) THEN 
         CALL arr_matmul(zei,lc) 
!                                                                       
!------ Add      two user variable arrays 'matadd'                      
!                                                                       
      ELSEIF (str_comp (bef, 'matadd', 4, lbef, 6) ) THEN 
         CALL arr_matadd(zei,lc) 
!                                                                       
!------ Calculate determinant             'detmat'                      
!                                                                       
      ELSEIF (str_comp (bef, 'detmat', 4, lbef, 6) ) THEN 
         CALL arr_detmat(zei,lc) 
!                                                                       
!------ Calculate inverse matrix          'invmat'                      
!                                                                       
      ELSEIF (str_comp (bef, 'invmat', 4, lbef, 6) ) THEN 
         CALL arr_invmat(zei,lc) 
!                                                                       
!------ Calculate transpose matrix          'mattrans'                      
!                                                                       
      ELSEIF (str_comp (bef, 'mattrans', 4, lbef, 8) ) THEN 
         CALL arr_transpose(zei,lc) 
      ELSEIF (str_comp (bef, 'show' , 2, lbef, 4) ) THEN 
         CALL get_params (zei, ianz, cpara, lpara, maxpar, lc) 
         CALL do_show_generic (cpara, lpara, MAXPAR)
!                                                                       
!     Reset the seed for the random number generator 'seed'             
!                                                                       
      ELSEIF (str_comp (bef, 'seed', 3, lbef, 4) ) THEN 
         CALL do_seed (zei, lc) 
!                                                                       
!                                                                       
!------ Various settings 'set'                                          
!                                                                       
      ELSEIF (str_comp (bef, 'set', 3, lbef, 3) ) THEN 
         CALL do_set (zei, lc) 
!     elseif(str_comp(bef,'memory',3,lbef,6)) THEN
!        CALL get_params (zei, ianz, cpara, lpara, maxpar, lc) 
!        CALL ber_params (ianz, cpara, lpara, werte, maxpar) 
!        if(ier_num == 0) THEN
!        write(cpara(1),'(I10)') NINT(werte(1))
!        ENDIF
!        call memory_message(cpara(1)(1:20))
!                                                                       
!------ Sleep fo a while 'sleep'                                        
!                                                                       
      ELSEIF (str_comp (bef, 'sleep', 2, lbef, 5) ) THEN 
         CALL get_params (zei, ianz, cpara, lpara, maxpar, lc) 
         IF (ier_num.eq.0) THEN 
            CALL ber_params (ianz, cpara, lpara, werte, maxpar) 
            IF (ier_num.eq.0) THEN 
               IF (nint (werte (1) ) .gt.0) THEN 
                  CALL do_sleep (nint (werte (1) ) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
!------ Handling a scocket connection                                   
!                                                                       
      ELSEIF (str_comp (bef, 'socket', 3, lbef, 6) ) THEN 
         CALL do_socket (zei, lc) 
!                                                                       
!------ Operating System Kommandos 'syst'                               
!                                                                       
      ELSEIF (str_comp (bef, 'system', 2, lbef, 6) ) THEN 
         command = ' ' 
         IF (zei.ne.' ') THEN 
            command (1:lc) = zei (1:lc) 
            CALL do_operating (command, lc) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!
!------ mount a WSL drive
!
      ELSEIF (str_comp (bef, 'mount', 3, lbef, 5) ) THEN 
         CALL mount_drive(zei)
!
!------ umount a WSL drive
!
      ELSEIF (str_comp (bef, 'umount', 3, lbef, 6) ) THEN 
         CALL umount_drive(zei)
!
!------ Update Discus
!
      ELSEIF (str_comp (bef, 'update', 3, lbef, 6) ) THEN 
         CALL lib_f90_update_discus
!                                                                       
!------ definition of variables                                         
!                                                                       
      ELSEIF (str_comp (bef, 'variable', 3, lbef, 8) ) THEN 
         CALL define_variable (zei, lc, NO_DIFFEV) 
!                                                                       
!------ Wait for user input 'wait'                                      
!                                                                       
      ELSEIF (str_comp (bef, 'wait', 3, lbef, 4) ) THEN 
         CALL do_input (zei, lc) 
!                                                                       
!------ Unknown command                                                 
!                                                                       
      ELSE 
         ier_num = - 8 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
END SUBROUTINE kdo_all                        
!
!*****7*****************************************************************
!
END MODULE kdo_all_mod
