MODULE transform_menu
!
CONTAINS
!
!+                                                                      
!     Generalized unit cell transformation operations:                  
!     The user provides the transformation matrix between new and       
!     old coordinate system and the program transforms the old          
!     atom coordinates into the new ones.                               
!     Individual direct and reciprocal space transformation can be      
!     calculated as well.                                               
!                                                                       
!*****7*****************************************************************
SUBROUTINE transform 
!-                                                                      
!     Main menu for generalized transformation operations               
!+                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE metric_mod
USE modify_mod
USE discus_show_menu
USE transfrm_mod 
USE trans_sup_mod
!
USE ber_params_mod
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE class_macro_internal
USE prompt_mod 
USE sup_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MIN_PARA = 20  ! A command requires at leaset these no of parameters
INTEGER maxw 
LOGICAL, PARAMETER :: lold = .FALSE. 
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara ! (MAXSCAT) 
INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara ! (MAXSCAT)
REAL(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte ! (MAXSCAT) 
!
CHARACTER(LEN=8) :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=PREC_STRING) :: line, zeile
CHARACTER(LEN=PREC_STRING) :: tran_hfile 
INTEGER :: lp, length, lbef 
INTEGER :: tran_hlen 
INTEGER :: indxg, ianz, i, j 
LOGICAL :: lend, lchange, lscreen 
LOGICAL :: lspace
REAL(KIND=PREC_DP), DIMENSION(3), PARAMETER :: NULLV = (/ 0.00, 0.00, 0.00 /)
REAL(KIND=PREC_DP), DIMENSION(3) :: u
REAL(KIND=PREC_DP)               :: dd
REAL(KIND=PREC_DP), DIMENSION(4) :: hkl !(4) 
!                                                                       
!                                                                       
!                                                                       
DATA lchange / .TRUE. / 
DATA lscreen / .TRUE. / 
!                                                                       
maxw = MAX(MIN_PARA,MAXSCAT+1)
lend = .FALSE. 
CALL no_error 
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/tran' 
!                                                                       
main: DO WHILE (.NOT.lend) 
   CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
   IF (ier_num == 0) THEN 
      IF (line /= ' '      .AND. line(1:1) /= '#' .AND. &
          line /= CHAR(13) .AND. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
         indxg = index (line, '=') 
         IF (indxg /= 0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) )    &
                       .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                       .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                                   str_comp (befehl, '?   ',   2, lbef, 4) )    &
                       .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
            CALL do_math (line, indxg, length) 
         ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1)  == '@') THEN 
                  IF (length.ge.2) THEN 
                     line(1:length-1) = line(2:length)
                     line(length:length) = ' '
                     length = length - 1
                     CALL file_kdo(line, length)
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----list asymmetric unit 'asym'                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) THEN 
                  CALL show_asym 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) THEN 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----list atoms present in the crystal 'chem'                      
!                                                                       
         ELSEIF (str_comp (befehl, 'chemstry', 2, lbef, 8) ) THEN 
                  CALL show_chem 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
         ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) THEN 
                  CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .TRUE. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
         ELSEIF(str_comp(befehl, 'help', 2, lbef, 4) .OR.         &
                str_comp(befehl, '?   ', 1, lbef, 4)      ) THEN                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 12 
               CALL do_hel ('discus tran '//zeile, lp) 
            ENDIF 
!                                                                       
!     ----reset transformation 'reset'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) THEN 
                  CALL tran_reset
         ELSE
!--------------------------------------------------------------------------------
!---------TRANSFORM specific commands
!--------------------------------------------------------------------------------
!
            IF( cr_nscat > TRAN_MAXSCAT .OR. cr_ncatoms>TRAN_MAXSITE) THEN
               CALL alloc_transfrm ( cr_nscat, cr_ncatoms )
               IF ( ier_num < 0 ) THEN
                  RETURN
               ENDIF
            ENDIF
!                                                                       
!     ----a(new) in terms of old axes                                   
!                                                                       
            IF(str_comp (befehl, 'anew', 2, lbef, 4) ) THEN 
               ier_num = - 6 
               ier_typ = ER_COMM 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF(ier_num == 0.AND.ianz == 3) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF(ier_num == 0) THEN 
                     u = werte(1:3)
                     lspace = .TRUE.
                     dd = do_blen(lspace, u, NULLV) 
                     IF(ier_num==0) THEN
                        DO j = 1, 3 
                           tran_g (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_G 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
               ENDIF 
!                                                                       
!     ----a(old) in terms of new axes                                   
!                                                                       
            ELSEIF (str_comp(befehl, 'aold', 2, lbef, 4) ) THEN 
               ier_num = - 6 
               ier_typ = ER_COMM 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num == 0.AND.ianz == 3) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF(ier_num == 0) THEN 
                     u = werte(1:3)
                     lspace = .TRUE.
                     dd = do_blen(lspace, u, NULLV) 
                     IF(ier_num==0) THEN
                        DO j = 1, 3 
                           tran_gi (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_GI 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
               ENDIF 
!                                                                       
!     ----b(new) in terms of old axes                                   
!                                                                       
            ELSEIF (str_comp (befehl, 'bnew', 2, lbef, 4) ) THEN 
               ier_num = - 6 
               ier_typ = ER_COMM 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num == 0.AND.ianz == 3) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num == 0) THEN 
                     u = werte(1:3)
                     lspace = .TRUE.
                     dd = do_blen(lspace, u, NULLV) 
                     IF(ier_num==0) THEN
                        DO j = 1, 3 
                           tran_g (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_G 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
               ENDIF 
!                                                                       
!     ----b(old) in terms of new axes                                   
!                                                                       
            ELSEIF (str_comp (befehl, 'bold', 2, lbef, 4) ) THEN 
               ier_num = - 6 
               ier_typ = ER_COMM 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num == 0.AND.ianz == 3) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num == 0) THEN 
                     u = werte(1:3)
                     lspace = .TRUE.
                     dd = do_blen(lspace, u, NULLV) 
                     IF(ier_num==0) THEN
                        DO j = 1, 3 
                           tran_gi (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_GI 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
               ENDIF 
!                                                                       
!     ----calculate a single transformation operation 'c2new'           
!                                                                       
            ELSEIF (str_comp (befehl, 'c2new', 2, lbef, 5) ) THEN 
               IF (lchange) THEN 
                  CALL tran_setup 
               ENDIF 
               IF(ier_num==0) THEN
                  lchange = .FALSE. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz == 3) THEN 
                        cpara (4) = 'd' 
                        lpara (4) = 1 
                     ENDIF 
                     ianz = 3 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO i = 1, 3 
                           hkl(i) = werte (i) 
                        ENDDO 
                        IF(str_comp(cpara(4), 'd', 1, lpara(4), 1)) THEN                                          
                           hkl (4) = 1.0 
                           CALL tran_ca (hkl, tran_f, lscreen) 
                        ELSEIF(str_comp(cpara(4), 'r', 1, lpara(4), 1) ) THEN                                     
                           hkl (4) = 0.0 
                           CALL tran_ca (hkl, tran_g, lscreen) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     ----calculate a single transformation operation 'c2old'           
!                                                                       
               ELSEIF (str_comp (befehl, 'c2old', 2, lbef, 5) ) THEN 
                  IF (lchange) THEN 
                     CALL tran_setup 
                  ENDIF 
                  IF(ier_num==0) THEN
                  lchange = .FALSE. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz == 3) THEN 
                        cpara (4) = 'd' 
                        lpara (4) = 1 
                     ENDIF 
                     ianz = 3 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO i = 1, 3 
                        hkl (i) = werte (i) 
                        ENDDO 
                        IF (str_comp (cpara (4) , 'd', 1, lpara (4) , 1)&
                        ) THEN                                          
                           hkl (4) = 1.0 
                           CALL tran_ca (hkl, tran_fi, lscreen) 
                        ELSEIF (str_comp (cpara (4) , 'r', 1, lpara (4) &
                        , 1) ) THEN                                     
                           hkl (4) = 0.0 
                           CALL tran_ca (hkl, tran_gi, lscreen) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
                  ENDIF 
!                                                                       
!     ----c(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'cnew', 2, lbef, 4) ) THEN 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_g (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_G 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----c(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'cold', 2, lbef, 4) ) THEN 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_gi (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_GI 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Deselect which atoms are included in the wave 'dese'          
!                                                                       
               ELSEIF (str_comp (befehl, 'deselect', 1, lbef, 8) ) THEN 
                  CALL atom_select (zeile, lp, 0, TRAN_MAXSCAT, tran_latom, &
                  tran_lsite, 0, TRAN_MAXSITE,                              &
                  tran_sel_atom, lold,.FALSE.)
!
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num == 0) THEN 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num == 0) THEN 
!                       IF (werte (1)  ==  - 1) THEN 
!                          DO i = 0, cr_nscat 
!                          tran_latom (i) = .FALSE. 
!                          ENDDO 
!                       ELSE 
!                          DO i = 1, ianz 
!                          tran_latom (NINT (werte (i) ) ) = .FALSE. 
!                          ENDDO 
!                       ENDIF 
!                    ENDIF 
!                 ENDIF 
!                                                                       
!     ----transforme a list of reflections to new base 'h2new'          
!                                                                       
               ELSEIF (str_comp (befehl, 'h2new', 2, lbef, 5) ) THEN 
                  IF (lchange) THEN 
                     CALL tran_setup 
                  ENDIF 
                  IF(ier_num==0) THEN
                  lchange = .FALSE. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz.ge.1) THEN 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 1)                                        
                        IF (ier_num == 0) THEN 
                           tran_hfile = cpara (1) 
                           tran_hlen = lpara (1) 
                           CALL tran_hkl (tran_hfile, tran_hlen, tran_g) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
                  ENDIF 
!                                                                       
!     ----transforme a list of reflections to old base 'h2old'          
!                                                                       
               ELSEIF (str_comp (befehl, 'h2old', 2, lbef, 5) ) THEN 
                  IF (lchange) THEN 
                     CALL tran_setup 
                  ENDIF 
                  IF(ier_num==0) THEN
                  lchange = .FALSE. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz == 1) THEN 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 1)                                        
                        IF (ier_num == 0) THEN 
                           tran_hfile = cpara (1) 
                           tran_hlen = lpara (1) 
                           CALL tran_hkl (tran_hfile, tran_hlen,        &
                           tran_gi)                                     
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
                  ENDIF 
!                                                                       
!     ----Select range of atoms within crystal to be included 'incl'    
!                                                                       
               ELSEIF (str_comp (befehl, 'include', 1, lbef, 7) ) THEN 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz == 2) THEN 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num == 0) THEN 
                           tran_start = NINT (werte (1) ) 
                           tran_end = NINT (werte (2) ) 
                        ENDIF 
                     ELSEIF (ianz == 1) THEN 
                        IF(str_comp(cpara(1), 'all', 1, lpara(1), 3) ) THEN                                       
                           tran_start = 1 
                           tran_end = - 1 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!------ ----Select the location of the new origin in old                
!                  coordinates 'onew'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'onew', 2, lbef, 4) ) THEN 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz == 3) THEN 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num == 0) THEN 
                           tran_oold = .TRUE. 
                           DO i = 1, 3 
                           tran_orig (i) = werte (i) 
                           ENDDO 
                           lchange = .TRUE. 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Select the location of the old origin in new                  
!                  coordinates 'oold'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'oold', 2, lbef, 4) ) THEN 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz == 3) THEN 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num == 0) THEN 
                           tran_oold = .FALSE. 
                           DO i = 1, 3 
                           tran_orig (i) = werte (i) 
                           ENDDO 
                           lchange = .TRUE. 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----run transformation 'run'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'run', 2, lbef, 4) ) THEN 
                  IF (lchange) THEN 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .FALSE. 
                  IF(ier_num==0) CALL tran_op (.TRUE.)
!                                                                       
!     ----Select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'select', 3, lbef, 6) ) THEN 
                  CALL atom_select (zeile, lp, 0, TRAN_MAXSCAT, tran_latom, &
                  tran_lsite, 0, TRAN_MAXSITE,                              &
                  tran_sel_atom, lold,.TRUE.)
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num == 0) THEN 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num == 0) THEN 
!                       IF (werte (1)  ==  - 1) THEN 
!                          DO i = 0, cr_nscat 
!                          tran_latom (i) = .TRUE. 
!                          ENDDO 
!                       ELSE 
!                          DO i = 1, ianz 
!                          tran_latom (NINT (werte (i) ) ) = .TRUE. 
!                          ENDDO 
!                       ENDIF 
!                    ENDIF 
!                 ENDIF 
!                                                                       
!     ----Set parameters for various commands 'set'                     
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0) THEN 
                     IF (ianz.gt.0) THEN 
                        IF (str_comp (cpara (1) , 'deltahkl', 1, lbef,  &
                        8) ) THEN                                       
                           cpara (1) = '0.0' 
                           lpara (1) = 3 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num == 0) THEN 
                              tran_deltahkl = werte (2) 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
                  IF (lchange) THEN 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .FALSE. 
                  call transfrm_show
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) THEN 
                  IF (zeile /= ' ') THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
                  CALL do_input (zeile, lp) 
!                                                                       
!     ----x(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'xnew', 2, lbef, 4)            &
               .OR.str_comp (befehl, 'asnew', 3, lbef, 5) ) THEN        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_f (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_F 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----x(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'xold', 2, lbef, 4)            &
               .OR.str_comp (befehl, 'asold', 3, lbef, 5) ) THEN        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_fi (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_FI 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----y(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'ynew', 2, lbef, 4)            &
               .OR.str_comp (befehl, 'bsnew', 3, lbef, 5) ) THEN        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_f (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_F 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----y(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'yold', 2, lbef, 4)            &
               .OR.str_comp (befehl, 'bsold', 3, lbef, 5) ) THEN        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_fi (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_FI 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----z(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'znew', 2, lbef, 4)            &
               .OR.str_comp (befehl, 'csnew', 3, lbef, 5) ) THEN        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_f (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_F 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----z(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'zold', 2, lbef, 4)            &
               .OR.str_comp (befehl, 'csold', 3, lbef, 5) ) THEN        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num == 0.AND.ianz == 3) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num == 0) THEN 
                        DO j = 1, 3 
                        tran_fi (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_FI 
                        lchange = .TRUE. 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
   IF (ier_num /= 0) THEN 
      lchange = .TRUE. 
      CALL errlist 
      IF (ier_sta /= ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in transform menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close(-1)
                  prompt_status = PROMPT_ON 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (lblock) THEN 
            ier_num = - 11 
            ier_typ = ER_COMM 
            prompt_status = PROMPT_ON 
            prompt = orig_prompt
            RETURN 
         ENDIF 
         CALL no_error 
         lmakro_error = .FALSE.
         sprompt = ' '
      ENDIF 
   ENDIF 
ENDDO  main
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE transform                      
!
!*****7*************************************************************************
!
subroutine transfrm_show
!
use crystal_mod
use transfrm_mod
!
use prompt_mod
!
implicit none
!
integer :: i, j
!
WRITE (output_io, 3000) 
WRITE (output_io, 3010) ( (tran_g (i, j), j = 1, 3), i = 1, 3)                                             
WRITE (output_io, 3020) ( (tran_gi(i, j), j = 1, 3), i = 1, 3)                                             
WRITE (output_io, 3030) ( (tran_f (i, j), j = 1, 4), i = 1, 3)                                             
WRITE (output_io, 3040) ( (tran_fi(i, j), j = 1, 4), i = 1, 3)                                             
WRITE (output_io, 3050) ( (tran_f (i, j), j = 1, 3), i = 1, 3)                                             
WRITE (output_io, 3060) ( (tran_fi(i, j), j = 1, 3), i = 1, 3)                                             
WRITE (output_io, 3061) ( (tran_g (i, j), j = 1, 3), i = 1, 3)                                             
WRITE (output_io, 3062) ( (tran_gi(i, j), j = 1, 3), i = 1, 3)                                             
!                                                                       
IF (tran_det /= 0) THEN 
   WRITE (output_io, 3100) tran_det, 1. / tran_det 
ELSE 
   WRITE (output_io, 3200) 
ENDIF 
!                                                                       
WRITE (output_io, 3090) 
WRITE (output_io, 3091) 
DO i = 0, cr_nscat 
   IF (tran_latom (i) ) THEN 
      WRITE (output_io, 3092) i, cr_at_lis (i) 
   ENDIF 
ENDDO 
IF (tran_end ==  - 1) THEN 
   WRITE (output_io, 3080) 
ELSE 
   WRITE (output_io, 3081) tran_start, tran_end 
ENDIF 
!                                                                       
 3000 FORMAT    (20x,'    Unit cell transformations '/                  &
     &                  20x,' ==============================='//)       
 3010 FORMAT    (                                                       &
     &           ' ( a(new) ) = ( ',2(F9.5,','),f9.5,' )   ( a(old) )'/ &
     &           ' ( b(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(old) )'/ &
     &           ' ( c(new) ) = ( ',2(F9.5,','),f9.5,' )   ( c(old) )'/)
 3020 FORMAT    (                                                       &
     &           ' ( a(old) ) = ( ',2(F9.5,','),f9.5,' )   ( a(new) )'/ &
     &           ' ( b(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(new) )'/ &
     &           ' ( c(old) ) = ( ',2(F9.5,','),f9.5,' )   ( c(new) )'/)
 3030 FORMAT    (                                                       &
     &           ' ( x(new) ) = ( ',2(F9.5,','),f9.5,' )   ( x(old) )', &
     &           '   (',f9.5,')'/                                       &
     &           ' ( y(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(old) )', &
     &           ' + (',f9.5,')'/                                       &
     &           ' ( z(new) ) = ( ',2(F9.5,','),f9.5,' )   ( z(old) )', &
     &           '   (',f9.5,')'/)                                      
 3040 FORMAT    (                                                       &
     &           ' ( x(old) ) = ( ',2(F9.5,','),f9.5,' )   ( x(new) )', &
     &           '   (',f9.5,')'/                                       &
     &           ' ( y(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(new) )', &
     &           ' + (',f9.5,')'/                                       &
     &           ' ( z(old) ) = ( ',2(F9.5,','),f9.5,' )   ( z(new) )', &
     &           '   (',f9.5,')'/)                                      
 3050 FORMAT    (                                                       &
     &           ' ( a*(new)) = ( ',2(F9.5,','),f9.5,' )   ( a*(old))'/ &
     &           ' ( b*(new)) = ( ',2(F9.5,','),f9.5,' ) * ( b*(old))'/ &
     &           ' ( c*(new)) = ( ',2(F9.5,','),f9.5,' )   ( c*(old))'/)
 3060 FORMAT    (                                                       &
     &           ' ( a*(old)) = ( ',2(F9.5,','),f9.5,' )   ( a*(new))'/ &
     &           ' ( b*(old)) = ( ',2(F9.5,','),f9.5,' ) * ( b*(new))'/ &
     &           ' ( c*(old)) = ( ',2(F9.5,','),f9.5,' )   ( c*(new))'/)
 3061 FORMAT    (                                                       &
     &           ' ( h(new) ) = ( ',2(F9.5,','),f9.5,' )   ( h(old) )'/ &
     &           ' ( k(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( k(old) )'/ &
     &           ' ( l(new) ) = ( ',2(F9.5,','),f9.5,' )   ( l(old) )'/)
 3062 FORMAT    (                                                       &
     &           ' ( h(old) ) = ( ',2(F9.5,','),f9.5,' )   ( h(new) )'/ &
     &           ' ( k(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( k(new) )'/ &
     &           ' ( l(old) ) = ( ',2(F9.5,','),f9.5,' )   ( l(new) )'/)
 3100 FORMAT    (                                                       &
     &           '   V(new)   = ',G13.5E3,' * V(old)'/                     &
     &           '   V(old)   = ',G13.5E3,' * V(new)'/)                    
 3200 FORMAT    (                                                       &
     &           '   Determinant of transformation is zero'/)           
 3080 FORMAT    (' Range of atoms from to    : All atoms included') 
 3081 FORMAT    (' Range of atoms from to    : ',2(2x,i9)) 
 3090 FORMAT    (' selected atoms    :') 
 3091 FORMAT    ('                      type name') 
 3092 FORMAT    (16x,2x,i8,1x,a4) 
!
end subroutine transfrm_show
!
!*****7*************************************************************************
!
SUBROUTINE tran_setup 
!-                                                                      
!     Performs the generalized symmetry operation                       
!     See Sands, D.E. Vectors and Tensors in Crystallography Chapt. 4.7 
!+                                                                      
USE discus_config_mod 
USE transfrm_mod 
!
use matrix_mod
USE errlist_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER :: i, j 
REAL(kind=PREC_DP), DIMENSION(3,3) ::  a !(3, 3) 
!                                                                       
!     initialize fourth columns and rows                                
!                                                                       
DO i = 1, 4 
      tran_g (i, 4) = 0.0 
      tran_g (4, i) = 0.0 
      tran_gi (i, 4) = 0.0 
      tran_gi (4, i) = 0.0 
      tran_f (i, 4) = 0.0 
      tran_f (4, i) = 0.0 
      tran_fi (i, 4) = 0.0 
      tran_fi (4, i) = 0.0 
ENDDO 
tran_g (4, 4) = 1.0 
tran_gi (4, 4) = 1.0 
tran_f (4, 4) = 1.0 
tran_fi (4, 4) = 1.0 
!                                                                       
!     Matrix g was defined                                              
!                                                                       
IF (tran_inp == TRAN_INP_G) THEN 
   DO i = 1, 4 
      DO j = 1, 4 
         tran_gi (i, j) = tran_g (i, j) 
      ENDDO 
   ENDDO 
!  CALL invmat4 (tran_gi) 
   CALL matinv  (tran_g , tran_gi) 
   IF(ier_num/=0) THEN
      ier_msg(2) = 'Check transformation matrices'
      RETURN
   ENDIF
   DO i = 1, 3 
      DO j = 1, 3 
         tran_f (j, i) = tran_gi (i, j) 
         tran_fi (j, i) = tran_g (i, j) 
      ENDDO 
   ENDDO 
!                                                                       
!     Matrix gi was defined                                             
!                                                                       
ELSEIF (tran_inp == TRAN_INP_GI) THEN 
   DO i = 1, 4 
      DO j = 1, 4 
         tran_g (i, j) = tran_gi (i, j) 
      ENDDO 
   ENDDO 
!  CALL invmat4 (tran_g) 
   CALL matinv  (tran_gi, tran_g) 
   DO i = 1, 3 
      DO j = 1, 3 
         tran_f (j, i) = tran_gi (i, j) 
         tran_fi (j, i) = tran_g (i, j) 
      ENDDO 
   ENDDO 
!                                                                       
!     Matrix f was defined                                              
!                                                                       
ELSEIF (tran_inp == TRAN_INP_F) THEN 
   DO i = 1, 4 
      DO j = 1, 4 
         tran_fi (i, j) = tran_f (i, j) 
      ENDDO 
   ENDDO 
!  CALL invmat4 (tran_fi) 
   CALL matinv  (tran_f , tran_fi) 
   DO i = 1, 3 
      DO j = 1, 3 
         tran_g (j, i) = tran_fi (i, j) 
         tran_gi (j, i) = tran_f (i, j) 
      ENDDO 
   ENDDO 
!                                                                       
!     Matrix fi was defined                                             
!                                                                       
ELSEIF (tran_inp == TRAN_INP_FI) THEN 
   DO i = 1, 4 
      DO j = 1, 4 
         tran_f (i, j) = tran_fi (i, j) 
      ENDDO 
   ENDDO 
!  CALL invmat4 (tran_f) 
   CALL matinv  (tran_fi, tran_f) 
   DO i = 1, 3 
      DO j = 1, 3 
         tran_g (j, i) = tran_fi (i, j) 
         tran_gi (j, i) = tran_f (i, j) 
      ENDDO 
   ENDDO 
ENDIF 
!                                                                       
!     Origin was defined in terms of the old coordinates                
!                                                                       
IF (tran_oold) THEN 
   j = 4 
   DO i = 1, 3 
         tran_fi (i, 4) = tran_orig (i) 
   ENDDO 
   DO i = 1, 3 
      DO j = 1, 3 
         tran_f (i, 4) = tran_f (i, 4) - tran_f (i, j) * tran_fi (j, 4) 
      ENDDO 
   ENDDO 
ELSE 
   j = 4 
   DO i = 1, 3 
      tran_f (i, 4) = tran_orig (i) 
   ENDDO 
   DO i = 1, 3 
      DO j = 1, 3 
         tran_fi (i, 4) = tran_fi (i, 4) - tran_fi (i, j) * tran_f (j,  &
         4)                                                             
      ENDDO 
   ENDDO 
ENDIF 
!                                                                       
DO i = 1, 3 
   DO j = 1, 3 
      a (i, j) = tran_g (i, j) 
   ENDDO 
ENDDO 
tran_det = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2) ) &
         + a(2, 1) * (a(3, 2) * a(1, 3) - a(1, 2) * a(3, 3) ) &
         + a(3, 1) * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2) )                 
!                                                                       
END SUBROUTINE tran_setup                     
!*****7*****************************************************************
!
SUBROUTINE tran_op (lout)
!-                                                                      
!     Performs the actual transformation operation.                     
!     All atoms of the structure, all symmetry elements are transformed 
!     into the new cell.                                                
!+                                                                      
USE discus_config_mod 
use chem_mod
USE crystal_mod 
USE metric_mod
USE spcgr_apply, ONLY: setup_lattice
USE update_cr_dim_mod
USE transfrm_mod 
!
USE errlist_mod 
use precision_mod
!                                                                       
IMPLICIT none 
!
logical, intent(in) :: lout  ! Output of new lattice
!                                                                       
INTEGER :: i, j 
INTEGER :: i1, i2, i3 
INTEGER :: i_start, i_end 
LOGICAL :: lspace 
REAL(kind=PREC_DP) ::  usym (4), ures (4) 
REAL(kind=PREC_DP) :: werte (5) 
REAL(kind=PREC_DP) :: u (3), v (3), w (3) 
!                                                                       
DATA usym / 0.0, 0.0, 0.0, 1.0 / 
DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = tran_start 
      i_end = tran_end 
      IF (tran_end ==  - 1) i_end = cr_natoms 
!                                                                       
!     Apply transformation operation to all atoms within selected range 
!                                                                       
      DO i = i_start, i_end 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (tran_latom (cr_iscat(1,i) ) ) THEN 
!                                                                       
!     ----Copy atom to temporary place                                  
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) 
         ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
         usym (4) = 1.0 
!        CALL trans (usym, tran_f, ures, 4) 
         ures = matmul(tran_f, usym)
!                                                                       
!     ----Replace original atom by its image                            
!                                                                       
         DO j = 1, 3 
         cr_pos (j, i) = ures (j) 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
!     If all atoms of the crystal have been selected, calculate new     
!     lattice constants and symmetry operations                         
!                                                                       
      IF (i_start == 1.AND.i_end == cr_natoms) THEN 
         v (1) = 0.0 
         v (2) = 0.0 
         v (3) = 0.0 
         lspace = .TRUE. 
         DO i = 1, 3 
         DO j = 1, 3 
         usym (j) = tran_g (i, j) 
         ENDDO 
         cr_a0 (i) = sqrt (skalpro (usym, usym, cr_gten) ) 
!                                                                       
         i1 = mod (i - 1, 3) + 1 
         i2 = mod (i, 3) + 1 
         i3 = mod (i + 1, 3) + 1 
         DO j = 1, 3 
         u (j) = tran_g (i1, j) 
         w (j) = tran_g (i2, j) 
         ENDDO 
         cr_win (i3) = do_bang (lspace, u, v, w) 
         ENDDO 
         CALL tran_sym 
         cr_spcgr = 'P1' 
         cr_spcgr_set = 'P1' 
         cr_set    = 'abc'
         cr_iset   =  1
         cr_spcgrno = 1 
         spcgr_ianz = 0 
         spcgr_para = 1 
!        lout = .TRUE. 
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian, cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
         CALL update_cr_dim 
!
   chem_period = .false.     ! Expect periodic boundary conditions
   chem_quick = .false.      ! and fast lookup mode to be wrong
   cr_icc = 1
!
ENDIF 
!                                                                       
END SUBROUTINE tran_op                        
!
!*****7*****************************************************************
!
      SUBROUTINE tran_sym 
!-                                                                      
!     Transforms the symmetry operations                                
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE gen_add_mod 
      USE sym_add_mod 
      USE transfrm_mod 
      USE unitcell_mod 
!                                                                       
      USE errlist_mod 
use precision_mod
      USE prompt_mod 
!
      IMPLICIT none 
!                                                                       
      INTEGER :: i, j, k, n, nn=0, jp, l 
      LOGICAL lequal 
!                                                                       
      REAL(kind=PREC_DP) ::  mat (4, 4) !, arr (4, 4) 
      REAL(kind=PREC_DP) ::  eps 
!                                                                       
!                                                                       
      DATA eps / 0.001D0 / 
!                                                                       
!------ The original symmetry elements are going to be stored as        
!       additional                                                      
!     generators. The original additional generators are shifted up.    
!                                                                       
      IF (gen_add_n + generspcgr (0, cr_spcgrno) .gt.gen_add_MAX) THEN 
         ier_typ = ER_APPL 
         ier_num = - 73 
         RETURN 
      ENDIF 
!                                                                       
!------ Shift original additional generators by generspcgr(0,cr_spcgrno)
!       entries                                                         
!                                                                       
      k = generspcgr (0, cr_spcgrno) 
      DO n = gen_add_n, 1, - 1 
      nn = n + k 
      DO i = 1, 4 
      DO j = 1, 4 
      gen_add (i, j, nn) = gen_add (i, j, n) 
      ENDDO 
      ENDDO 
      gen_add_power (nn) = gen_add_power (n) 
      ENDDO 
!                                                                       
!     copy space group generators into additional generators            
!                                                                       
      DO n = 1, generspcgr (0, cr_spcgrno) 
      IF (gen_sta == GEN_SYMM) THEN 
         nn = generspcgr (n, cr_spcgrno) 
      ELSEIF (gen_sta == GEN_CENTER) THEN 
         nn = generspcgr_center (n, cr_spcgrno) 
      ENDIF 
      DO i = 1, 4 
      DO j = 1, 4 
      gen_add (i, j, n) = generators (i, j, nn) 
      ENDDO 
      ENDDO 
      gen_add_power (n) = generpower (nn) 
      ENDDO 
      gen_add_n = gen_add_n + generspcgr (0, cr_spcgrno) 
!                                                                       
!     transform the additional generators                               
!                                                                       
      DO n = 1, gen_add_n 
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j) = gen_add (i, j, n) 
      ENDDO 
      ENDDO 
!     CALL matmul4 (arr, tran_f, mat) 
!     CALL matmul4 (mat, arr, tran_fi) 
      mat = matmul( matmul(tran_f, mat), tran_fi)
      DO i = 1, 4 
      DO j = 1, 4 
      gen_add (i, j, n) = mat (i, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     transform the additional symmetry elements                        
!                                                                       
      DO n = 1, sym_add_n 
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j) = sym_add (i, j, n) 
      ENDDO 
      ENDDO 
!     CALL matmul4 (arr, tran_f, mat) 
!     CALL matmul4 (mat, arr, tran_fi) 
      mat = matmul( matmul(tran_f, mat), tran_fi)
      DO i = 1, 4 
      DO j = 1, 4 
      sym_add (i, j, n) = mat (i, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     Check transformation of primitive translations, if it results in  
!     non-integer translation, add as additional generator.             
!                                                                       
      k = gen_add_n 
      DO n = 1, 3 
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j) = 0.0 
      ENDDO 
      mat (i, i) = 1.0 
      ENDDO 
      mat (n, 4) = 1.0 
!     CALL matmul4 (arr, tran_f, mat) 
!     CALL matmul4 (mat, arr, tran_fi) 
      mat = matmul( matmul(tran_f, mat), tran_fi)
      IF(ABS(mod(mat(1, 4), 1.D0) ) .gt. eps .OR.  &
         ABS(mod(mat(2, 4), 1.D0) ) .gt. eps .OR.  &
         ABS(mod(mat(3, 4), 1.D0) ) .gt. eps      ) THEN       
         IF (gen_add_n + 1.gt.gen_add_MAX) THEN 
            ier_typ = ER_APPL 
            ier_num = - 73 
            RETURN 
         ENDIF 
!                                                                       
!     ----Determine power of generator                                  
!                                                                       
!                                                                       
!     ----After this loop j holds the power,which transforms the        
!         generator into a primitive translation                        
!                                                                       
!                                                                       
!     ----This loop transforms the translation into 0<=t<1              
!                                                                       
         DO i = 1, 3 
         mat(i, 4) = mod(mat(i, 4), 1.0D0) 
         IF(mat(i, 4) .lt.0.0) mat(i, 4) = mat(i, 4) + 1.0 
         ENDDO 
         jp = 1 
         DO WHILE (ABS (jp * mat (1, 4) - NINT (jp * mat (1, 4) ) )     &
         .gt.eps.OR.ABS (jp * mat (2, 4) - NINT (jp * mat (2, 4) ) )    &
         .gt.eps.OR.ABS (jp * mat (3, 4) - NINT (jp * mat (3, 4) ) )    &
         .gt.eps)                                                       
         jp = jp + 1 
         ENDDO 
!                                                                       
!     ----Compare to previous generators originating from primitive     
!         translations. If identical, skip this one.                    
!                                                                       
         lequal = .TRUE. 
         DO i = k + 1, gen_add_n 
         DO j = 1, 3 
         lequal = lequal.AND. (ABS (gen_add (j, 4, i) - mat (j, 4) )    &
         .lt.eps)                                                       
         ENDDO 
         ENDDO 
         IF (.NOT.lequal.OR.k == gen_add_n) THEN 
            gen_add_n = gen_add_n + 1 
            gen_add_power (gen_add_n) = 1 
            IF (jp.gt.1) jp = jp - 1 
            gen_add_power (gen_add_n) = jp 
            DO i = 1, 4 
            DO j = 1, 4 
            gen_add (i, j, gen_add_n) = mat (i, j) 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     Shift the generators resulting from the primitive translations    
!     up front                                                          
!                                                                       
      n = gen_add_n - k 
      DO i = 1, n 
      DO j = 1, 4 
      DO k = 1, 4 
      mat (j, k) = gen_add (j, k, gen_add_n) 
      ENDDO 
      ENDDO 
      jp = gen_add_power (gen_add_n) 
      DO l = gen_add_n, 2, - 1 
      DO j = 1, 4 
      DO k = 1, 4 
      gen_add (j, k, l) = gen_add (j, k, l - 1) 
      ENDDO 
      ENDDO 
      gen_add_power (l) = gen_add_power (l - 1) 
      ENDDO 
      DO j = 1, 4 
      DO k = 1, 4 
      gen_add (j, k, 1) = mat (j, k) 
      ENDDO 
      ENDDO 
      gen_add_power (1) = jp 
      ENDDO 
!                                                                       
      END SUBROUTINE tran_sym                       
!*****7*****************************************************************
      SUBROUTINE tran_hkl (infile, infile_l, matrix) 
!-                                                                      
!     Transforms a list of reflections, read from file.                 
!+                                                                      
      USE discus_config_mod 
      USE transfrm_mod 
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE precision_mod
      USE prompt_mod 
USE support_mod
!
IMPLICIT none 
!
CHARACTER(len=*)  , intent(in) :: infile 
INTEGER           ,  intent(in) :: infile_l
REAL(kind=PREC_DP), intent(in) ::  matrix (4, 4) 
!                                                                       
      INTEGER :: ird, iwr, irs 
      INTEGER :: i, j 
      LOGICAL :: lread, ltest 
!                                                                       
      CHARACTER(len=20) :: inten 
      CHARACTER(LEN=PREC_STRING) :: outfile 
      CHARACTER(LEN=PREC_STRING) :: restfile 
      INTEGER :: hkl (3) 
REAL(kind=PREC_DP) ::  usym (4), ures (4), utest (3) 
!                                                                       
      DATA ird, iwr, irs / 7, 8, 9 / 
!                                                                       
!     Open input file, and the two output files                         
!                                                                       
      lread = .TRUE. 
      CALL oeffne (ird, infile, 'unknown') 
      IF (ier_num /= 0) THEN 
         CLOSE (ird) 
         RETURN 
      ENDIF 
      lread = .FALSE. 
      outfile = infile (1:infile_l) //'.trans' 
      CALL oeffne (iwr, outfile, 'unknown') 
      IF (ier_num /= 0) THEN 
         CLOSE (ird) 
         CLOSE (iwr) 
         RETURN 
      ENDIF 
      lread = .FALSE. 
      restfile = infile (1:infile_l) //'.rest' 
      CALL oeffne (irs, restfile, 'unknown') 
      IF (ier_num /= 0) THEN 
         CLOSE (ird) 
         CLOSE (iwr) 
         CLOSE (irs) 
         RETURN 
      ENDIF 
!                                                                       
!-----      Loop over all reflections in input file                     
!                                                                       
cond_read: do
      READ (ird, 2000, end = 999) hkl, inten 
      DO i = 1, 3 
      usym (i) = REAL(hkl (i) ) 
      ENDDO 
      usym (4) = 0.0 
!                                                                       
!-----      -- Apply transformation operation                           
!                                                                       
!     CALL trans (usym, matrix, ures, 4) 
      ures = matmul(matrix, usym)
!                                                                       
!     -- Write new integer reflections to output file, non-integer      
!        reflections to restfile                                        
!                                                                       
      ltest = .TRUE. 
      DO i = 1, 3 
      utest (i) = ABS (REAL(NINT (ures (i) ) ) - ures (i) ) 
      ltest = ltest.AND. (utest (i) .lt.tran_deltahkl) 
      ENDDO 
      IF (ltest) THEN 
         WRITE (iwr, 2000) (NINT (ures (j) ), j = 1, 3), inten 
      ELSE 
         WRITE (irs, 3000) (ures (j), j = 1, 3), inten, (NINT (usym (j) &
         ), j = 1, 3)                                                   
      ENDIF 
enddo cond_read
!
  999 CONTINUE 
      CLOSE (ird) 
      CLOSE (iwr) 
      CLOSE (irs) 
!                                                                       
!                                                                       
 2000 FORMAT    (3i4,a20) 
 3000 FORMAT    (3f8.3,1x,a20,1x,3i4) 
      END SUBROUTINE tran_hkl                       
!
!*******************************************************************************
!
SUBROUTINE tran_reset
!
USE discus_allocate_appl_mod
USE transfrm_mod
!
IMPLICIT NONE
!
INTEGER :: ik
!
CALL alloc_transfrm(1, 1)
!
TRAN_MAXSCAT = 1
TRAN_MAXSITE = 1
!
tran_start   =  1
tran_end     = -1
tran_inp     = TRAN_INP_G
IF(ALLOCATED(tran_latom)) tran_latom(:) = .FALSE.  ! (0:TRAN_MAXSCAT)
IF(ALLOCATED(tran_lsite)) tran_lsite(:) = .TRUE.   ! (0:TRAN_MAXSCAT)
!
tran_oold      = .TRUE.
tran_sel_atom  = .TRUE.
tran_orig(3)   = 0.0
tran_det       = 1.0
tran_g   (:,:) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_g ))
tran_gi  (:,:) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_gi))
tran_f   (:,:) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_f ))
tran_fi  (:,:) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_fi))
tran_deltahkl  = 0.001
!
END SUBROUTINE tran_reset
!
!*******************************************************************************
!
END MODULE transform_menu
