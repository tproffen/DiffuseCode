MODULE chem_menu
!
CONTAINS
!*****7*****************************************************************
!                                                                       
SUBROUTINE chem 
!-                                                                      
!     This sublevel contains all routines to analyse a structure        
!     line calculation of bondlength distributions, the average         
!     structure and correlations.                                       
!                                                                       
!     Note: Some variables are used in the MC section as well           
!           and settings might be overwritten.                          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_aver_mod
      USE chem_symm_mod
      USE celltoindex_mod
      USE get_iscat_mod
      USE modify_mod
!
      USE build_name_mod
      USE calc_expr_mod
      USE doact_mod 
      USe do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE get_params_mod
      USE learn_mod 
      USE class_macro_internal 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
      USE param_mod 
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      USE string_convert_mod
      USE sup_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
CHARACTER(len=13) :: befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(LEN=PREC_STRING) line, zeile, cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw), wwerte (maxw), wwwerte (maxw) 
      REAL(KIND=PREC_DP) :: uwerte (maxw) 
      INTEGER lbeg (3) 
      INTEGER lpara (maxw), lp, length 
      INTEGER indxg, ianz, lbef, iianz, jjanz 
      INTEGER kkanz 
      LOGICAL lout , lsite
!                                                                       
!                                                                       
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/chem' 
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
             line == char(13) .or. line(1:1) == '!'  ) GOTO 10
!                                                                       
!------ search for "="                                                  
!                                                                       
indxg = index (line, '=') 
IF (indxg /= 0 .AND. .NOT. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
               .AND. .NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
               .AND. .NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                            str_comp (befehl, '?   ', 2, lbef, 4) )    &
               .AND. INDEX(line, '==') == 0                        ) THEN
            CALL do_math (line, indxg, length) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            line(1:length-1) = line(2:length)
            line(length:length) = ' '
            length = length - 1
            CALL file_kdo(line, length)
!
!------ command 'symmetry'
!
         ELSEIF (str_comp (befehl, 'apply_symmetry', 2, lbef, 13) ) then
            CALL chem_symm(zeile, lp)
!                                                                       
!------ Calculate average structure and sigmas 'aver'                   
!                                                                       
         ELSEIF (str_comp (befehl, 'aver', 2, lbef, 4) ) then 
            IF (.not.chem_sel_atom) then 
               ier_num = - 22 
               ier_typ = ER_CHEM 
            ELSE 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF(ier_num == 0) THEN
                  IF(IANZ==0 .or. str_comp (cpara(1), 'one', 2, lpara(1), 3)) THEN
                     lsite = .true.
                     CALL chem_aver (.true., lsite ) 
                  ELSEIF(IANZ== 1 .and. str_comp (cpara(1), 'ind', 2, lpara(1), 3)) THEN
                     lsite = .false.
                     CALL chem_aver (.true., lsite ) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF
            ENDIF 
!                                                                       
!------ Calculate bond-angle distribution 'bang'                        
!                                                                       
         ELSEIF (str_comp (befehl, 'bang', 2, lbef, 4) ) then 
            IF (.not.chem_sel_atom) then 
               ier_num = - 22 
               ier_typ = ER_CHEM 
               RETURN 
            ENDIF 
!                                                                       
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.ge.3) then 
                  iianz = 1 
                  jjanz = 1 
                  kkanz = 1 
                  CALL get_iscat (iianz, cpara, lpara, werte, maxw,     &
                  .false.)                                              
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  IF (ier_num.eq.0) then 
                     CALL get_iscat (jjanz, cpara, lpara, wwerte, maxw, &
                     .false.)                                           
                     IF (ier_num.eq.0) then 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (kkanz, cpara, lpara, uwerte,    &
                        maxw, .false.)                                  
                        IF (ier_num.eq.0) then 
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL do_build_name (ianz, cpara, lpara,      &
                           wwwerte, maxw, 1)                            
                           chem_fname = cpara (1) (1:lpara(1))
                           IF (ier_num.eq.0) then 
                              CALL chem_bang (iianz, jjanz, kkanz,      &
                              werte, wwerte, uwerte, maxw)              
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!------ Calculate bond-length distribution 'blen'                       
!                                                                       
         ELSEIF (str_comp (befehl, 'blen', 2, lbef, 4) ) then 
            IF (.not.chem_sel_atom) then 
               ier_num = - 22 
               ier_typ = ER_CHEM 
               RETURN 
            ENDIF 
!                                                                       
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.ge.2) then 
                  iianz = 1 
                  jjanz = 1 
                  CALL get_iscat (iianz, cpara, lpara, werte, maxw,     &
                  .false.)                                              
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  IF (ier_num.eq.0) then 
                     CALL get_iscat (jjanz, cpara, lpara, wwerte, maxw, &
                     .false.)                                           
                     IF (ier_num.eq.0) then 
                        IF (ianz.gt.1) then 
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL do_build_name (ianz, cpara, lpara,      &
                           wwwerte, maxw, 1)                            
                           chem_fname = cpara (1) (1:lpara(1))
                        ELSE 
                           chem_fname = 'discus.blen' 
                        ENDIF 
                        IF (ier_num.eq.0) then 
                           IF (chem_cluster) then 
                              CALL chem_blen_cluster (iianz, jjanz,     &
                              werte, wwerte, maxw)                      
                           ELSE 
                              CALL chem_blen (iianz, jjanz, werte,      &
                              wwerte, maxw)                             
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!------ 'bval' calculates the bond valence sum for a given atom         
!                                                                       
         ELSEIF (str_comp (befehl, 'bval', 2, lbef, 4) ) then 
            CALL chem_bval (zeile, lp) 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
            CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
            CALL echo (zeile, lp) 
!                                                                       
!------ Show relative amounts of elements 'elem'                        
!                                                                       
         ELSEIF (str_comp (befehl, 'element', 2, lbef, 7) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) then 
               CALL do_cap (cpara (1) ) 
               lout = (cpara (1) (1:2) .eq.'ON') 
            ELSE 
               lout = .true. 
            ENDIF 
!                                                                       
            IF (chem_sel_atom) then 
               CALL chem_elem (lout) 
            ELSE 
               CALL chem_mole (lout) 
            ENDIF 
!                                                                       
!------ 'env' finds neighbours around given atom, site or position      
!                                                                       
         ELSEIF (str_comp (befehl, 'env', 2, lbef, 3) ) then 
            CALL chem_env (zeile, lp) 
!                                                                       
!------ 'mode' toggles between working with molecules and atoms         
!                                                                       
         ELSEIF (str_comp (befehl, 'mode', 2, lbef, 4) ) then 
            CALL chem_mode (zeile, lp) 
!                                                                       
!------ 'neig' finds neighbours according to given neighbour definition 
!                                                                       
         ELSEIF (str_comp (befehl, 'neig', 2, lbef, 4) ) then 
            CALL chem_nei (zeile, lp) 
!                                                                       
!------ Evaluate an expression 'eval'                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
            CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!------ calculate correlation field                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'field', 2, lbef, 5) ) then 
            CALL chem_corr_field (zeile, lp) 
!                                                                       
!------ calculate neighbour probabilities 'corr'                        
!                                                                       
         ELSEIF (str_comp (befehl, 'corr', 3, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.3) then 
               IF (chem_ctyp (1) .eq.CHEM_NONE) then 
                  ier_num = - 11 
                  ier_typ = ER_CHEM 
               ELSE 
                  IF (str_comp (cpara (1) , 'occ', 1, lpara (1) , 3) )  &
                  then                                                  
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     lbeg (1) = - 1 
                     lbeg (2) = 0 
                     lbeg (3) = 0 
!                                                                       
                     IF (chem_sel_atom) then 
                        CALL chem_corr_occ (cpara, lpara, maxw,   &
                        .true., lbeg)                                   
                     ELSE 
                        CALL chem_corr_occ_mol (ianz, cpara, lpara,     &
                        maxw, .true., lbeg)                             
                     ENDIF 
                  ELSEIF (str_comp (cpara (1) , 'disp', 1, lpara (1) ,  &
                  4) ) then                                             
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     IF (chem_sel_atom) then 
                        CALL chem_corr_dis (ianz, cpara, lpara, maxw,   &
                        .true., lbeg)                                   
                     ELSE 
                        CALL chem_corr_dis_mol (ianz, cpara, lpara,     &
                        maxw, .true., lbeg)                             
                     ENDIF 
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ calculate angles for given neighbours                           
!                                                                       
         ELSEIF (str_comp (befehl, 'angle', 3, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.ne.0) return 
            IF (ianz.ge.2) then 
               IF (chem_ctyp (1) .eq.CHEM_NONE) then 
                  ier_num = - 11 
                  ier_typ = ER_CHEM 
               ELSE 
                  IF (ianz.gt.2) then 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 3)                                           
                  ENDIF 
                  IF (chem_sel_atom) then 
                     CALL chem_angle (ianz, cpara, lpara, werte, wwerte,&
                     wwwerte, maxw, .true.)                             
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ calculate distortions for given neighbours                      
!                                                                       
         ELSEIF (str_comp (befehl, 'disp', 3, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.ne.0) return 
            IF (ianz.ge.2) then 
               IF (chem_ctyp (1) .eq.CHEM_NONE) then 
                  ier_num = - 11 
                  ier_typ = ER_CHEM 
               ELSE 
                  IF (ianz.gt.2) then 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 3)                                           
                  ENDIF 
                  IF (chem_sel_atom) then 
                     CALL chem_disp (ianz, cpara, lpara, werte, maxw,   &
                     .true.)                                            
                  ELSE 
                     CALL chem_disp_mol (ianz, cpara, lpara, werte,     &
                     maxw, .true.)                                      
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ 'homo' checks how homogeneous the crystal is                    
!                                                                       
         ELSEIF (str_comp (befehl, 'homo', 2, lbef, 4) ) then 
            CALL chem_homo (zeile, lp) 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.&
              str_comp (befehl, '?   ', 1, lbef, 4) ) then
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 12 
               CALL do_hel ('discus chem '//zeile, lp) 
            ENDIF 
!                                                                       
!------ reset all parameters for 'chem' section: 'reset'                   
!                                                                       
         ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) then 
            CALL chem_reset
!                                                                       
!------ set most parameters for 'chem' section: 'set'                   
!                                                                       
         ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
            CALL chem_set (zeile, lp) 
!                                                                       
!------ show parameters 'show'                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            zeile = 'all' 
            CALL chem_show (zeile) 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) then 
            IF (zeile.ne.' ') then 
               CALL do_operating (zeile (1:lp), lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ calc cell-atomindex and atomindex-cell                          
!                                                                       
         ELSEIF (str_comp (befehl, 'trans', 2, lbef, 5) ) then 
            CALL chem_trans(zeile,lp)
!                                                                       
!------ Waiting for user input                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
            CALL do_input (zeile, lp) 
!                                                                       
!------ no command found                                                
!                                                                       
         ELSE 
            ier_num = - 8 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
!------ any errors ?                                                    
!                                                                       
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in chemistry menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN 
               ELSE
                  IF(lmacro_close) THEN
                    CALL macro_close 
                    prompt_status = PROMPT_ON 
                 ENDIF 
               ENDIF 
            ENDIF 
            IF (lblock) then 
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
      GOTO 10 
!                                                                       
 9999 CONTINUE 
      prompt = orig_prompt
      END SUBROUTINE chem                           
!*****7*****************************************************************
      SUBROUTINE chem_show (cmd) 
!+                                                                      
!     show current parameters                                           
!-                                                                      
      USE discus_config_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) cmd 
      INTEGER i, j, ic 
      LOGICAL lnone 
!                                                                       
      CALL do_cap (cmd) 
!                                                                       
!------ General settings                                                
!                                                                       
      IF (cmd (1:2) .eq.'AL'.or.cmd (1:2) .eq.'GE') then 
         WRITE (output_io, 1000) chem_quick, chem_period 
         IF (chem_sel_atom) then 
            WRITE (output_io, 1050) 'atoms' 
         ELSE 
            WRITE (output_io, 1050) 'molecules' 
         ENDIF 
         IF (ilots.eq.LOT_OFF) then 
            WRITE (output_io, 1100) 
         ELSEIF (ilots.eq.LOT_BOX) then 
            WRITE (output_io, 1110) nlots 
            WRITE (output_io, 1130) (ls_xyz (i), i = 1, 3), lperiod 
         ELSEIF (ilots.eq.LOT_ELI) then 
            WRITE (output_io, 1120) nlots 
            WRITE (output_io, 1130) (ls_xyz (i), i = 1, 3), lperiod 
         ENDIF 
      ENDIF 
!                                                                       
!------ Bond length distribution settings                               
!                                                                       
      IF (cmd (1:2) .eq.'AL'.or.cmd (1:2) .eq.'BO') then 
         WRITE (output_io, 2000) chem_blen_cut 
         WRITE (output_io, 2100) chem_bin 
      ENDIF 
!                                                                       
!------ Correlation determination settings                              
!                                                                       
      IF (cmd (1:2) .eq.'AL'.or.cmd (1:2) .eq.'CO') then 
         WRITE (output_io, 3000) 
         WRITE (output_io, 3010) 
!                                                                       
         DO ic = 1, chem_ncor 
!                                                                       
         IF (chem_ctyp (ic) .eq.CHEM_NONE) then 
            WRITE (output_io, 3050) 
!                                                                       
         ELSEIF (chem_ctyp (ic) .eq.CHEM_DIST) then 
            WRITE (output_io, 3100) ic, (chem_neig (i, 1, ic), i = 1, 3)&
            , chem_freq_sigma (ic), chem_wink_sigma (ic), chem_rmin (ic)&
            , chem_rmax (ic), chem_cang (ic), chem_nnei (ic)            
         ELSEIF (chem_ctyp (ic) .eq.CHEM_VEC) then 
            WRITE (output_io, 3200) ic, (chem_use_vec (i, ic), i = 1,   &
            chem_nvec (ic) )                                            
         ELSEIF (chem_ctyp (ic) .eq.CHEM_ANG) then 
            WRITE (output_io, 3300) ic, (chem_use_win (i, ic), i = 1,   &
            chem_nwin (ic) )                                            
         ELSEIF (chem_ctyp (ic) .eq.CHEM_ENVIR) then 
            WRITE (output_io, 3400) ic, (chem_use_env (i, ic), i = 1,   &
            chem_nenv (ic) )                                            
         ELSEIF (chem_ctyp (ic) .eq.CHEM_RANGE) then 
            WRITE (output_io, 3500) ic, (chem_use_ran (i, ic), i = 1,   &
            chem_nran (ic) )                                            
         ELSEIF (chem_ctyp (ic) .eq.CHEM_CON) then 
            WRITE (output_io, 3600) ic, (chem_use_con (i, ic), i = 1,   &
            chem_ncon (ic) )                                            
         ENDIF 
         ENDDO 
!                                                                       
!------ - Show defined interaction vectors                              
!                                                                       
         lnone = .true. 
         WRITE (output_io, 4000) 
         DO i = 1, CHEM_MAX_VEC 
         IF (chem_cvec (1, i) .ne. - 9999) then 
            WRITE (output_io, 4010) i, (chem_cvec (j, i), j = 1, 5) 
            lnone = .false. 
         ENDIF 
         ENDDO 
         IF (lnone) WRITE (output_io, 3050) 
!                                                                       
!------ - Show defined interaction angles                               
!                                                                       
         lnone = .true. 
         WRITE (output_io, 4500) 
         DO i = 1, CHEM_MAX_ANG 
         IF (chem_cwin (1, i) .ne. - 9999) then 
            WRITE (output_io, 4510) i, (chem_cwin (j, i), j = 1, 9) 
            lnone = .false. 
         ENDIF 
         ENDDO 
         IF (lnone) WRITE (output_io, 3050) 
!                                                                       
!------ - Show defined interaction environments                         
!                                                                       
         lnone = .true. 
         WRITE (output_io, 4600) 
         DO i = 1, CHEM_MAX_ENV 
         IF (chem_cenv (0, i) .ne. - 9999) then 
            WRITE (output_io, 4610) i, chem_cenv (0, i), chem_rmin_env (&
            i), chem_rmax_env (i)                                       
            WRITE (output_io, 4620) (chem_cenv (j, i), j = 1,           &
            chem_env_neig (i) )                                         
            lnone = .false. 
         ENDIF 
         ENDDO 
         IF (lnone) WRITE (output_io, 3050) 
!                                                                       
!------ - Show defined interaction ranges                               
!                                                                       
         lnone = .true. 
         WRITE (output_io, 4700) 
         DO i = 1, CHEM_MAX_RAN 
         IF (chem_cran_uvw (1, 1, i) .ne. - 9999) then 
            WRITE (output_io, 4710) i, (chem_cran_uvw (j, 1, i),        &
            j = 1, 3), chem_cran_sig (i), chem_cran_wsig (i)            
            IF (chem_cran_cent (0, i) .eq. - 1) then 
               WRITE (output_io, 4720) 
            ELSE 
               WRITE (output_io, 4730) (chem_cran_cent (j, i), j = 1,   &
               chem_cran_cent (0, i) )                                  
            ENDIF 
            IF (chem_cran_neig (0, i) .eq. - 1) then 
               WRITE (output_io, 4740) 
            ELSE 
               WRITE (output_io, 4750) (chem_cran_neig (j, i), j = 1,   &
               chem_cran_neig (0, i) )                                  
            ENDIF 
            IF (chem_cran_short (i) ) then 
               WRITE (output_io, 4760) chem_cran_nshort (i) 
            ENDIF 
            lnone = .false. 
         ENDIF 
         ENDDO 
         IF (lnone) WRITE (output_io, 3050) 
!                                                                       
!------ - Show defined displacement directions                          
!                                                                       
         lnone = .true. 
         WRITE (output_io, 5000) 
         DO i = 1, chem_ncor 
         IF (chem_ldall (i) ) then 
            WRITE (output_io, 5005) i 
            lnone = .false. 
         ELSE 
            IF (chem_dir (1, 1, i) .ne. - 9999) then 
               WRITE (output_io, 5010) i, (chem_dir (j, 1, i), j = 1, 3)&
               , (chem_dir (j, 2, i), j = 1, 3)                         
               lnone = .false. 
            ENDIF 
         ENDIF 
         ENDDO 
         IF (lnone) WRITE (output_io, 3050) 
!                                                                       
!------ - Show defined interaction connectivities                              
!                                                                       
         lnone = .true. 
         WRITE (output_io, 6000) 
         DO i = 1, CHEM_MAX_CON 
         IF (chem_ccon (1, i) .ne. - 9999) then 
            WRITE (output_io, 6010) i, (chem_ccon (j, i), j = 1, 2), &
            chem_cname(i)(1:chem_cname_l(i))
            lnone = .false. 
         ENDIF 
         ENDDO 
         IF (lnone) WRITE (output_io, 3050) 
      ENDIF 
!                                                                       
 1000 FORMAT ('    Neighbour determination mode   : quick = ',L1,/      &
     &        '    Periodic boundaries (x,y,z)    : ',3(L1,1x))         
 1050 FORMAT ('    Current operation mode         : ',A) 
 1100 FORMAT ('    Sample volume for homo-check   : complete crystal') 
 1110 FORMAT ('    Sample volume for homo-check   : ',I4,               &
     &        ' box shaped lots')                                       
 1120 FORMAT ('    Sample volume for homo-check   : ',I4,               &
     &        ' ellipsoid shaped lots')                                 
 1130 FORMAT ('    Lot size                       : ',I3,' x ',I3,      &
     &        ' x ',I3,' unit cells ',/,                                &
     &        '    Periodic boundaries            : ',L1)               
 2000 FORMAT ('    Allowed bond length range      : ',                  &
     &        F6.3,' A to ',F6.3,' A')                                  
 2100 FORMAT ('    # points for histogramms       : ',I6,/) 
 3000 FORMAT ('    Neighbour definitions          : ',/) 
 3010 FORMAT ('        #  Mode  Neighbour (or vec.)  fsig  wsig',       &
     &        '   rmin[A] rmax[A] ang sym',/,7x,67('-'))                
 3050 FORMAT (7x,'** none defined **') 
 3100 FORMAT (7x,i2,'  neig',2x,3(f5.2,1x),2x,2(f5.2,1x),               &
     &                   1x,2(f7.3,1x),2x,l1,2x,i2)                     
 3200 FORMAT (7x,i2,'  vec ',2x,50i4) 
 3300 FORMAT (7x,i2,'  ang ',2x,50i4) 
 3400 FORMAT (7x,i2,'  env ',2x,50i4) 
 3500 FORMAT (7x,i2,'  ran ',2x,50i4) 
 3600 FORMAT (7x,i2,'  con ',2x,50i4) 
 4000 FORMAT (/,'    Defined correlation vectors    : ',/) 
 4010 FORMAT ('       Correlation vector ',I3,'      : ',               &
     &        'Site',I3,' -> site',I3,', neig =',3I3)                   
 4500 FORMAT (/,'    Defined correlation angles     : ',/) 
 4510 FORMAT ('       Correlation angle  ',I3,'      : ',               &
     &        'Site',I3,' -> site',I3,', neig =',3I3,/,                 &
     &        44x,' -> site',I3,', neig =',3I3)                         
 4600 FORMAT (/,'    Defined correlation environments: ',/) 
 4610 FORMAT ('       Correlation envir  ',I3,'      : ',               &
     &        'Zentral Atom ',I4,/                                      &
     &        37x,'r min, r max : ',f8.3,2x,f8.3)                       
 4620 FORMAT (   37x,'neighbours   : ',15I4,/) 
 4700 FORMAT (/,'    Defined correlation ranges: ',/) 
 4710 FORMAT ('       Correlation range ',i3,': ',                      &
     &        3(F7.3,2x),'+-',F7.3,' A  +-',F7.3,'ø')                   
 4720 FORMAT ('            central atom    ',' ALL') 
 4730 FORMAT ('            central atom    ',50I4) 
 4740 FORMAT ('            neigh.  atom    ',' ALL') 
 4750 FORMAT ('            neigh.  atom    ',50I4) 
 4760 FORMAT ('            Only the nearest',  I4,                      &
     &        ' atoms are considered neighbours')                       
 5000 FORMAT (/,'    Defined displacement direc.    : ',/) 
 5005 FORMAT ('       Displacement corr. #',I3,'     : ',               &
     &        'all directions')                                         
 5010 FORMAT ('       Displacement corr. #',I3,'     : ',               &
     &        'A:',3(f5.2,1x),' B:',3(f5.2,1x))                         
 6000 FORMAT (/,'    Defined correlation connectiv. : ',/) 
 6010 FORMAT ('       Correlation connec ',I3,'      : ',               &
     &        'Type',I3,' -> def.',I3,1X, A)
      END SUBROUTINE chem_show                      
!*****7*****************************************************************
      SUBROUTINE chem_set (zeile, lp) 
!+                                                                      
!     sets most parameters for 'chem' section                           
!-                                                                      
      USE discus_config_mod 
USE discus_allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))) :: cpara (maxw) 
      REAL(KIND=PREC_DP):: werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER indxx, indxy, indxz 
INTEGER :: n_bin ! Dummy for histogramm allocation
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.2) then 
            CALL do_cap (cpara (1) ) 
!                                                                       
!------ --- 'set neig': setting correlation determination method        
!                                                                       
            IF (cpara (1) (1:2) .eq.'NE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_neig (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set cryst': setting crystal dimensions                     
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'CR') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.4) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     cr_icc (1) = nint (werte (1) ) 
                     cr_icc (2) = nint (werte (2) ) 
                     cr_icc (3) = nint (werte (3) ) 
                     cr_ncatoms = nint (werte (4) ) 
                     chem_purge = .FALSE.     ! Crystal dimension should allow periodic boundary
                  ELSE
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     --- 'set lots': sample volume (lots)                              
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'LOT') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               IF (ianz.eq.1) then 
                  CALL do_cap (cpara (1) ) 
                  IF (cpara (1) (1:1) .eq.'O') then 
                     ilots = LOT_OFF 
                     nlots = 1 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSEIF (ianz.eq.6) then 
                  CALL do_cap (cpara (1) ) 
                  IF (cpara (1) (1:1) .eq.'B') then 
                     ilots = LOT_BOX 
                  ELSEIF (cpara (1) (1:1) .eq.'E') then 
                     ilots = LOT_ELI 
                  ELSE 
                     ier_num = - 2 
                     ier_typ = ER_FOUR 
                  ENDIF 
                  IF (ier_num.eq.0) then 
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     CALL ber_params (ianz - 1, cpara, lpara, werte,    &
                     maxw)                                              
                     IF (ier_num.eq.0) then 
                        ls_xyz (1) = nint (werte (1) ) 
                        ls_xyz (2) = nint (werte (2) ) 
                        ls_xyz (3) = nint (werte (3) ) 
                        nlots = nint (werte (4) ) 
                        CALL do_cap (cpara (5) ) 
                        lperiod = (cpara (5) (1:1) .eq.'Y') 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- set mode: set calculation mode to quick/exact               
!                                                                       
     ELSEIF(cpara(1)(1:2) == 'MO') THEN 
        CALL do_cap (cpara (2) ) 
        IF(chem_purge .AND. (cpara(2)(1:3) ==  'QUI' .OR. &
                             cpara(3)(1:3) ==  'PER'     )) THEN
           ier_num = -31
           ier_typ = ER_CHEM
           ier_msg(1) = "Use >set crystal< in chem to define "
           ier_msg(2) = "Number of unit cells and atoms per unit cell"
           ier_msg(3) = "Or read a new cell/structure"
        ELSE
           chem_quick   = (cpara(2)(1:3) .eq.'QUI') 
           chem_cluster = (cpara(2)(1:3) .eq.'CLU') 
           IF(ianz >= 3) THEN 
              CALL do_cap(cpara(3) ) 
              IF(cpara(3)(1:3)  == 'PER') THEN 
                 IF(ianz == 4) THEN 
                    CALL do_cap(cpara(4) ) 
                    indxx = INDEX (cpara(4) , 'X') 
                    indxy = INDEX (cpara(4) , 'Y') 
                    indxz = INDEX (cpara(4) , 'Z') 
                    chem_period(1) = indxx.gt.0 
                    chem_period(2) = indxy.gt.0 
                    chem_period(3) = indxz.gt.0 
                 ELSE 
                    chem_period(1) = .TRUE. 
                    chem_period(2) = .TRUE. 
                    chem_period(3) = .TRUE. 
                 ENDIF 
              ELSE 
                 chem_period(1) = .FALSE. 
                 chem_period(2) = .FALSE. 
                 chem_period(3) = .FALSE. 
              ENDIF 
           ELSE 
              chem_period(1) = .FALSE. 
              chem_period(2) = .FALSE. 
              chem_period(3) = .FALSE. 
           ENDIF 
        ENDIF 
!                                                                       
!------ --- set bin: set number of points for histogramm binning        
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'BI') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               IF (ianz.eq.1) then 
                  n_bin = NINT(werte(1))
                  IF(n_bin >=  CHEM_MAX_BIN) THEN
                     CALL alloc_chem_hist(n_bin)
                  ENDIF
                  chem_bin = n_bin
!                 ELSE 
!                    ier_num = - 1 
!                    ier_typ = ER_CHEM 
!                 ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  chem_bin = 1
               ENDIF 
!                                                                       
!------ --- set blen: set allowed range for bond-length histogramm      
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'BL') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               IF (ianz.eq.2) then 
                  IF (werte (1) .gt.0.01D0.and.werte (2) .gt.werte (1) )  &
                  then                                                  
                     chem_blen_cut (1) = werte (1) 
                     chem_blen_cut (2) = werte (2) 
                  ELSE 
                     ier_num = - 3 
                     ier_typ = ER_CHEM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- set bang: set allowed range for bond-angle histogramm       
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'BA') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               IF (ianz.eq.2) then 
                  IF (werte (1) .gt.0.00D0.and.werte (2) .gt.werte (1) )  &
                  then                                                  
                     chem_bang_cut (1) = werte (1) 
                     chem_bang_cut (2) = werte (2) 
                  ELSE 
                     ier_num = - 25 
                     ier_typ = ER_CHEM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- set ang: sets correlation angles                            
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'AN') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_angle (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- set vec: sets correlation vectors                           
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'VE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_vec (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- set env: sets correlation environment                       
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'EN') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_envir (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- unknown subcommand entered                                  
!                                                                       
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
      END SUBROUTINE chem_set                       
!*****7*****************************************************************
SUBROUTINE chem_set_vec (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     'set vec' is processed here                                       
!-                                                                      
USE discus_allocate_appl_mod 
USE discus_config_mod 
USE crystal_mod 
USE chem_mod 
USE ber_params_mod
USE errlist_mod 
USE precision_mod
USE str_comp_mod
IMPLICIT none 
!                                                                       
       
!                                                                       
INTEGER                           , INTENT(IN)    :: ianz 
INTEGER                           , INTENT(IN)    :: maxw 
CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
REAL(KIND=PREC_DP),DIMENSION(MAXW), INTENT(INOUT) :: werte
INTEGER          , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!                                                                       
INTEGER     :: is1, is2, iv
INTEGER     :: n_vec  ! Dummy for allocations
INTEGER     :: n_cor  ! Dummy for allocations
!                                                                       
iv = 0
!                                                                       
IF (str_comp (cpara (1) , 'reset', 2, lpara (1) , 5) ) then 
   chem_cvec           =     0 ! all elements (i,j)
   chem_cvec    (1, :) = -9999 ! column (1,*)
   chem_use_vec        =     1 ! all elements (i,j)
!        DO i = 1, CHEM_MAX_VEC 
!        chem_cvec (1, i) = - 9999 
!        DO j = 2, 5 
!        chem_cvec (j, i) = 0 
!        ENDDO 
!        DO j = 1, CHEM_MAX_COR 
!        chem_use_vec (i, j) = 1 
!        ENDDO 
!        ENDDO 
   RETURN 
ENDIF 
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
IF (ianz.eq.6) then 
   iv  = nint (werte (1) ) 
   is1 = nint (werte (2) ) 
   is2 = nint (werte (3) ) 
!
!  allocate vectors
!
   IF (iv > 0 ) THEN
      IF (iv > CHEM_MAX_VEC .OR. CHEM_MAX_COR>UBOUND(chem_nvec,1)) THEN
         n_vec = CHEM_MAX_VEC + 10   ! Increment size by 10
         n_cor = CHEM_MAX_COR
         CALL alloc_chem_vec ( n_vec , n_cor )
         IF (ier_num /= 0) RETURN
      ENDIF
   ENDIF
   IF (iv.le.0.or.iv.gt.CHEM_MAX_VEC) then 
      ier_num = - 9 
      ier_typ = ER_CHEM 
   ELSE 
      IF (is1.le.0.or.is1.gt.cr_ncatoms.or.       &
         is2.le.0.or.is2.gt.cr_ncatoms    ) then                                                            
         ier_num = - 10 
         ier_typ = ER_CHEM 
      ELSE 
         chem_cvec (1, iv) = int( werte (2) )
         chem_cvec (2, iv) = int( werte (3) )
         chem_cvec (3, iv) = int( werte (4) )
         chem_cvec (4, iv) = int( werte (5) )
         chem_cvec (5, iv) = int( werte (6) )
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!
chem_ndef(CHEM_VEC) = max(chem_ndef(CHEM_VEC), iv)
!                                                                       
END SUBROUTINE chem_set_vec                   
!*****7*****************************************************************
SUBROUTINE chem_set_con (ianz, cpara, lpara, werte, maxw) 
!+
!     'set con' is processed here                                       
!-
   USE discus_allocate_appl_mod 
   USE discus_config_mod 
   USE crystal_mod 
   USE chem_mod 
   USE conn_sup_mod
   USE get_iscat_mod
   USE modify_mod
   USE ber_params_mod
   USE errlist_mod 
   USE get_params_mod
USE lib_errlist_func
USE precision_mod
USE str_comp_mod
   IMPLICIT none 
!
!
   INTEGER                           , INTENT(INOUT) :: ianz 
   INTEGER                           , INTENT(IN)    :: maxw 
   CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
   REAL(KIND=PREC_DP),DIMENSION(MAXW), INTENT(INOUT) :: werte
   INTEGER          , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!                                                                       
   CHARACTER (LEN=256)  :: c_name   ! Connectivity name
   INTEGER              :: c_name_l ! connectivity name length
   INTEGER     :: is1, ino, iianz, j
   INTEGER     :: iv = 0 ! con number in chemistry
   INTEGER     :: is     ! Loop over atom types
   INTEGER     :: j_con  ! Number of connectivities for an atom type
   INTEGER     :: n_con  ! Dummy for allocations
   INTEGER     :: n_cor  ! Dummy for allocations
   LOGICAL     :: lold   ! Atom types have to be present
!                                                                       
!                                                                       
   lold = .true.
   IF (str_comp (cpara (1) , 'reset', 2, lpara (1) , 5) ) then 
      chem_ccon           =     0 ! all elements (i,j)
      chem_ccon    (1, :) = -9999 ! column (1,*)
      chem_use_con        =     1 ! all elements (i,j)
      iv = 0                      ! no more con's
      RETURN 
   ENDIF 
!                                                                       
   IF(str_comp (cpara (1) , 'all', 2, lpara (1) , 3) ) THEN  ! Use all connectivities
      atom_types: DO is = 1, cr_nscat                    ! Loop over all atom types
         j_con = get_connectivity_numbers(is)            ! How many conns does this type have?
         connecs: DO j=1, j_con                          ! Loop over all connectivities 
            ino = j                                      ! identity may change ino
            c_name   = ' '                               ! Force search by number
            c_name_l = 1
            CALL get_connectivity_identity( is, ino, c_name, c_name_l)
            iv = iv + 1                                  ! autoincrement con.nr.
            IF (iv > CHEM_MAX_CON .OR. CHEM_MAX_COR>UBOUND(chem_ncon,1)) THEN
               n_con = CHEM_MAX_CON + 10                 ! Increment size by 10
               n_cor = CHEM_MAX_COR
               CALL alloc_chem_con ( n_con , n_cor )
               IF (ier_num /= 0) RETURN
            ENDIF
            chem_ccon (1, iv) = is
            chem_ccon (2, iv) = ino
            chem_cname  (iv)  = c_name
            chem_cname_l(iv)  = c_name_l
         ENDDO connecs
      ENDDO atom_types
      continue
   ELSE                                                  ! use all connectivities ???
   iianz = 1
   CALL ber_params (iianz, cpara, lpara, werte, maxw) 
   IF (ier_num.ne.0) return 
   nparams: IF (ianz.eq.3) then 
      iv  = nint (werte (1) )                            ! Set mmc connectivity number
      CALL del_params (1, ianz, cpara, lpara, maxw)      ! Remove Para 1
      CALL ber_params (iianz, cpara, lpara, werte, maxw) ! try to calc atom type
      IF (ier_num.ne.0) THEN                             ! Error must be atom name
         CALL get_iscat (iianz, cpara, lpara, werte, maxw, lold) 
         is1 = nint (werte (1) ) 
         CALL no_error
      ELSE                                               ! Success set to value
         is1 = nint (werte (1) )
      ENDIF 
      CALL del_params (1, ianz, cpara, lpara, maxw)      ! Remove Atom type
      CALL ber_params (iianz, cpara, lpara, werte, maxw) ! try to calc connect number
      IF (ier_num.ne.0) THEN                             ! Error must be a name
         c_name   = cpara(1)
         c_name_l = lpara(1)
         ino      = 0
         CALL no_error
      ELSE                                               ! Success set to value
         ino = nint (werte (1) ) 
         c_name   = ' '
         c_name_l = 1
      ENDIF
      CALL get_connectivity_identity( is1, ino, c_name, c_name_l)
!
!     allocate list of mmc connectivities
!
      IF (iv > 0 ) THEN
         IF (iv > CHEM_MAX_CON .OR. CHEM_MAX_COR>UBOUND(chem_ncon,1)) THEN
            n_con = CHEM_MAX_CON + 10   ! Increment size by 10
            n_cor = CHEM_MAX_COR
            CALL alloc_chem_con ( n_con , n_cor )
            IF (ier_num /= 0) RETURN
         ENDIF
      ENDIF
      IF (iv.le.0.or.iv.gt.CHEM_MAX_CON) then 
         ier_num = -28 
         ier_typ = ER_CHEM 
      ELSE 
         IF (is1.lt.0.or.is1.gt.cr_nscat) THEN
            ier_num = - 10 
            ier_typ = ER_CHEM 
         ELSE 
            chem_ccon (1, iv) = is1
            chem_ccon (2, iv) = ino
            chem_cname  (iv)  = c_name
            chem_cname_l(iv)  = c_name_l
         ENDIF 
      ENDIF 
   ELSE nparams
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF nparams
   ENDIF                              ! use all connectivities ???
!
!
chem_ndef(CHEM_CON) = max(chem_ndef(CHEM_CON), iv)
!                                                                       
END SUBROUTINE chem_set_con                   
!*****7*****************************************************************
SUBROUTINE chem_set_ranges (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     'set range' is processed here                                     
!-                                                                      
USE discus_allocate_appl_mod 
USE discus_config_mod 
USE crystal_mod 
USE chem_mod 
USE get_iscat_mod
USE metric_mod
USE modify_mod
USE rmc_symm_mod
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER                           , INTENT(INOUT) :: ianz 
INTEGER                           , INTENT(IN)    :: maxw 
CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
REAL(KIND=PREC_DP),DIMENSION(MAXW), INTENT(INOUT) :: werte
INTEGER          , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!                                                                       
INTEGER, PARAMETER  :: max_uvw = 48 
!                                                                       
INTEGER             :: iv, i, j
INTEGER             :: janz 
INTEGER             :: n_ran  ! Dummy for allocations
INTEGER             :: n_cor  ! Dummy for allocations
LOGICAL,PARAMETER   :: lold = .true.
LOGICAL             :: lacentric 
!                                                                       
REAL(KIND=PREC_DP) :: u (3), v (3) 
REAL(KIND=PREC_DP) :: uvw (4, max_uvw) 
REAL(KIND=PREC_DP) :: uvw_mat (4, 4, max_uvw) 
!                                                                       
iv = 0
!                                                                       
main: IF (str_comp (cpara (1) , 'reset', 2, lpara (1) , 5) ) then 
         chem_cran_uvw           = 0      ! (i,j,k)
         chem_cran_uvw (1, 1, :) = - 9999 
         chem_cran_sig           = 0.00   ! (i)
         chem_cran_wsig          = 0.00   ! (i)
         chem_use_ran            = 1      ! (i,j)
         chem_cran_cent (0, :)   = 0      ! (0,i)
         chem_cran_neig (0, :)   = 0      ! (0,i)
         RETURN 
ELSE main
!
   second: IF (str_comp (cpara (1) , 'direc', 2, lpara (1) , 5) ) then 
!                                                                       
!---- ---Define the direction, sigma of direction and sigma of angle    
!                                                                       
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      IF (ianz.gt.3) THEN 
         iv = nint (werte (1) ) 
!
!  allocate ranges
!
            IF (iv > 0 ) THEN
               IF (iv > CHEM_MAX_RAN .OR. CHEM_MAX_COR>UBOUND(chem_nran,1)) THEN
                  n_ran = CHEM_MAX_RAN + 10   ! Increment size by 10
                  n_cor = CHEM_MAX_COR
                  CALL alloc_chem_vec ( n_ran , n_cor )
                  IF (ier_num /= 0) RETURN
               ENDIF
            ENDIF
!                                                                       
            chem_cran_cang (iv) = .false. 
            chem_cran_lsym (iv) = .false. 
!                                                                       
!     ----Is last parameter "sym" or "nosym" ??                         
!                                                                       
            IF (str_comp(cpara(ianz),'sym',3,lpara(ianz),3) ) THEN
               chem_cran_lsym (iv) = .true. 
               ianz = ianz - 1 
            ELSEIF (str_comp (cpara (ianz) , 'nosym', 3, lpara (ianz) , &
            5) ) THEN
               chem_cran_lsym (iv) = .false. 
               ianz = ianz - 1 
            ENDIF 
!                                                                       
            chem_cran_uvw (1, 1, iv) = werte (2) 
            chem_cran_uvw (2, 1, iv) = werte (3) 
            chem_cran_uvw (3, 1, iv) = werte (4) 
            chem_cran_sig (iv) = werte (5) 
            IF (ianz.eq.6) THEN 
               chem_cran_wsig (iv) = werte (6) 
               chem_cran_cang (iv) = .true. 
            ENDIF 
!                                                                       
            DO i = 1, 3 
            u (i) = 0.0 
            v (i) = chem_cran_uvw (i, 1, iv) 
            ENDDO 
            chem_cran_rmax (iv) = do_blen (.true., u, v) 
            chem_cran_rmin (iv) = chem_cran_rmax (iv) - chem_cran_sig ( &
            iv) / 2.0                                                   
            chem_cran_rmax (iv) = chem_cran_rmax (iv) + chem_cran_sig ( &
            iv) / 2.0                                                   
!                                                                       
!------ - get symmetrically equivalent directions if needed             
!                                                                       
            IF (chem_cran_lsym (iv) ) then 
               DO i = 1, 3 
               uvw (i, 1) = chem_cran_uvw (i, 1, iv) 
               ENDDO 
               uvw (4, 1) = 0.0 
               DO i = 1, 4 
               DO j = 1, 4 
               uvw_mat (i, j, 1) = 0.0 
               ENDDO 
               uvw_mat (i, i, 1) = 1.0 
               ENDDO 
!                                                                       
               CALL rmc_symmetry (chem_cran_nuvw (iv), uvw, uvw_mat,    &
               max_uvw, .true., lacentric)                              
!                                                                       
               DO i = 1, chem_cran_nuvw (iv) 
               DO j = 1, 3 
               chem_cran_uvw (j, i, iv) = uvw (j, i) 
               ENDDO 
               ENDDO 
!                                                                       
            ELSE 
               chem_cran_nuvw (iv) = 1 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
   ELSEIF (str_comp (cpara (1) , 'central', 2, lpara (1) , 7) ) then 
!                                                                       
!     --Define the central atom(s)                                      
!                                                                       
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         janz = 1 
         CALL ber_params (janz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         IF (ianz.ge.2) then 
            iv = nint (werte (1) ) 
!
!  allocate ranges
!
            IF (iv > 0 ) THEN
               IF (iv > CHEM_MAX_RAN .OR. CHEM_MAX_COR>UBOUND(chem_nran,1)) THEN
                  n_ran = CHEM_MAX_RAN + 10   ! Increment size by 10
                  n_cor = CHEM_MAX_COR
                  CALL alloc_chem_vec ( n_ran , n_cor )
                  IF (ier_num /= 0) RETURN
               ENDIF
            ENDIF
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
            IF (ier_num.ne.0) return 
            IF (NINT(werte(1)) ==  -1) then 
               chem_cran_cent (0, iv) = - 1 
            ELSE 
               chem_cran_cent (0, iv) = ianz 
               DO i = 1, ianz 
               chem_cran_cent (i, iv) = nint (werte (i) ) 
               ENDDO 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
   ELSEIF (str_comp (cpara (1) , 'neig', 2, lpara (1) , 4) ) then 
!                                                                       
!     --Define the neighboring atom(s)                                  
!                                                                       
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         janz = 1 
         CALL ber_params (janz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         IF (ianz.ge.2) then 
            iv = nint (werte (1) ) 
!
!  allocate ranges
!
            IF (iv > 0 ) THEN
               IF (iv > CHEM_MAX_RAN .OR. CHEM_MAX_COR>UBOUND(chem_nenv,1)) THEN
                  n_ran = CHEM_MAX_RAN + 10   ! Increment size by 10
                  n_cor = CHEM_MAX_COR
                  CALL alloc_chem_vec ( n_ran , n_cor )
                  IF (ier_num /= 0) RETURN
               ENDIF
            ENDIF
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
            IF(ier_num.ne.0) return 
            IF(NINT(werte(1))  ==  -1) then 
               chem_cran_neig (0, iv) = - 1 
            ELSE 
               chem_cran_neig (0, iv) = ianz 
               DO i = 1, ianz 
               chem_cran_neig (i, iv) = nint (werte (i) ) 
               ENDDO 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
   ELSEIF (str_comp (cpara (1) , 'short', 2, lpara (1) , 5) ) then 
!                                                                       
!     --Define the number of the closest atoms to be considered         
!                                                                       
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         janz = 1 
         CALL ber_params (janz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         IF (ianz.eq.2) then 
            iv = nint (werte (1) ) 
!
!  allocate ranges
!
            IF (iv > 0 ) THEN
               IF (iv > CHEM_MAX_RAN .OR. CHEM_MAX_COR>UBOUND(chem_nvec,1)) THEN
                  n_ran = CHEM_MAX_RAN + 10   ! Increment size by 10
                  n_cor = CHEM_MAX_COR
                  CALL alloc_chem_vec ( n_ran , n_cor )
                  IF (ier_num /= 0) RETURN
               ENDIF
            ENDIF
            IF (str_comp (cpara (2) , 'none', 2, lpara (1) , 4) ) then 
               chem_cran_short (iv) = .false. 
               chem_cran_nshort (iv) = - 1 
            ELSE 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               chem_cran_short (iv) = .true. 
               chem_cran_nshort (iv) = nint (werte (ianz) ) 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
   ELSE second
         ier_num = - 6 
         ier_typ = ER_COMM 
         ier_msg (1) = cpara (1) 
   ENDIF second
ENDIF  main
!
chem_ndef(CHEM_RANGE) = max(chem_ndef(CHEM_RANGE), iv)
!                                                                       
      END SUBROUTINE chem_set_ranges                
!*****7*****************************************************************
SUBROUTINE chem_set_angle (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     'set angle' is processed here                                     
!-                                                                      
USE discus_allocate_appl_mod 
USE discus_config_mod 
USE crystal_mod 
USE chem_mod 
USE ber_params_mod
USE errlist_mod 
USE precision_mod
USE str_comp_mod
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER                           , INTENT(IN)    :: ianz 
INTEGER                           , INTENT(IN)    :: maxw 
CHARACTER (LEN=*), DIMENSION(maxw), INTENT(INOUT) :: cpara
REAL(KIND=PREC_DP),DIMENSION(MAXW), INTENT(INOUT) :: werte
INTEGER          , DIMENSION(MAXW), INTENT(INOUT) :: lpara
!                                                                       
INTEGER            :: is1, is2, is3, iv
!
INTEGER     :: n_ang  ! Dummy for allocations
INTEGER     :: n_cor  ! Dummy for allocations
!                                                                       
iv = 0
!                                                                       
IF (str_comp (cpara (1) , 'reset', 2, lpara (1) , 5) ) then 
   chem_cwin        =     0 ! All elements (i,j)
   chem_cwin (1, :) = -9999 ! Row (1,*)
   chem_use_win     =     1 ! All elements  (i,j)
!        DO i = 1, CHEM_MAX_ANG 
!        chem_cwin (1, i) = - 9999 
!        DO j = 2, 5 
!        chem_cwin (j, i) = 0 
!        ENDDO 
!        DO j = 1, CHEM_MAX_COR 
!        chem_use_win (i, j) = 1 
!        ENDDO 
!        ENDDO 
   RETURN 
ENDIF 
!                                                                       
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
IF (ianz.eq.10) then 
   iv = nint (werte (1) ) 
   is1 = nint (werte (2) ) 
   is2 = nint (werte (3) ) 
   is3 = nint (werte (7) ) 
!
!  allocate angles
!
   IF (iv > 0 ) THEN
      IF (iv > CHEM_MAX_ANG .OR. CHEM_MAX_COR>UBOUND(chem_nwin,1)) THEN
         n_ang = CHEM_MAX_ANG + 10   ! Increment size by 10
         n_cor = CHEM_MAX_COR
         call alloc_chem_ang ( n_ang , n_cor )
         IF (ier_num /= 0) RETURN
      ENDIF
   ENDIF
   IF (iv.le.0.or.iv.gt.CHEM_MAX_ANG) then 
      ier_num = - 24 
      ier_typ = ER_CHEM 
   ELSE 
      IF (is1.le.0.or.is1.gt.cr_ncatoms.or.  &
          is2.le.0.or.is2.gt.cr_ncatoms.or.  &
          is3.le.0.or.is3.gt.cr_ncatoms) then                           
         ier_num = - 10 
         ier_typ = ER_CHEM 
      ELSE 
         chem_cwin (1, iv) = int( werte (2) )
         chem_cwin (2, iv) = int( werte (3) )
         chem_cwin (3, iv) = int( werte (4) )
         chem_cwin (4, iv) = int( werte (5) )
         chem_cwin (5, iv) = int( werte (6) )
         chem_cwin (6, iv) = int( werte (7) )
         chem_cwin (7, iv) = int( werte (8) )
         chem_cwin (8, iv) = int( werte (9) )
         chem_cwin (9, iv) = int( werte (10) )
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!
chem_ndef(CHEM_ANG) = max(chem_ndef(CHEM_ANG), iv)
!                                                                       
      END SUBROUTINE chem_set_angle                 
!*****7*****************************************************************
      SUBROUTINE chem_set_envir (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     'set environment' is processed here                               
!-                                                                      
      USE discus_allocate_appl_mod 
      USE discus_config_mod 
USE atom_env_mod
      USE crystal_mod 
      USE chem_mod 
      USE get_iscat_mod
      USE modify_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, janz, iv, i, j
      INTEGER            :: n_atom
      INTEGER            :: n_env
      INTEGER            :: n_cor
!                                                                       
!                                                                       
      IF (str_comp (cpara (1) , 'reset', 2, lpara (1) , 5) ) then 
        chem_cenv           =     0 ! all elements (i,j)
        chem_cvec    (1, :) = -9999 ! column (1,*)
        chem_use_env        =     1 ! all elements (i,j)
        chem_env_neig       =     0 ! all elements (i)
!        DO i = 1, CHEM_MAX_ENV 
!        chem_cenv (0, i) = - 9999 
!        DO j = 1, MAX_ATOM_ENV 
!        chem_cenv (j, i) = 0 
!        ENDDO 
!        DO j = 1, CHEM_MAX_COR 
!        chem_use_env (i, j) = 1 
!        chem_env_neig (j) = 0 
!        ENDDO 
!        ENDDO 
         RETURN 
      ENDIF 
!                                                                       
!                                                                       
      janz = 3 
      CALL ber_params (janz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
      iv = nint (werte (1) ) 
!
!  allocate vectors
!
      IF (iv > 0 ) THEN
         IF (iv > CHEM_MAX_ENV .OR. CHEM_MAX_COR>UBOUND(chem_nenv,1)) THEN
            n_atom= MAX(n_atom,MAX_ATOM_ENV)
            n_env = CHEM_MAX_ENV + 10   ! Increment size by 10
            n_cor = CHEM_MAX_COR
            CALL alloc_chem_env ( n_atom , n_env, n_cor )
            IF (ier_num /= 0) RETURN
         ENDIF
      ENDIF
!
      IF (iv.le.0.or.iv.gt.CHEM_MAX_ENV) then 
         ier_num = - 26 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
!                                                                       
      chem_rmin_env (iv) = werte (2) 
      chem_rmax_env (iv) = werte (3) 
      CALL del_params (3, ianz, cpara, lpara, maxw) 
      chem_env_neig (iv) = ianz - 1 
      janz = ianz 
      DO i = 1, janz 
      j = 1 
      CALL get_iscat (j, cpara, lpara, werte, maxw, .false.) 
      chem_cenv (i - 1, iv) = nint (werte (1) ) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      ENDDO 
!
chem_ndef(CHEM_ENVIR) = max(chem_ndef(CHEM_ENVIR), iv)
!                                                                       
      END SUBROUTINE chem_set_envir                 
!*****7*****************************************************************
!
SUBROUTINE chem_set_neig (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     Command 'set neig' processed here                                 
!-                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE atom_env_mod
USE chem_mod 
USE metric_mod
USE rmc_sup_mod
USE rmc_symm_mod
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE string_convert_mod
use take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(INOUT) :: ianz 
INTEGER, INTENT(IN)    :: maxw 
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: cpara !(maxw) 
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara !(maxw) 
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(INOUT) :: werte !(maxw) 
!                                                                       
INTEGER, PARAMETER :: max_uvw = 48
      INTEGER i, j, iv , ii
INTEGER :: n_cor    ! Dummy number of correlations
INTEGER :: n_ang    ! Dummy number of vectors
INTEGER :: n_con    ! Dummy number of vectors
INTEGER :: n_env    ! Dummy number of vectors
INTEGER :: n_ran    ! Dummy number of vectors
INTEGER :: n_vec    ! Dummy number of vectors
INTEGER :: n_atom   ! Dummy number of vectors
      REAL(KIND=PREC_DP) :: u (3), v (3) 
      REAL(KIND=PREC_DP) :: uvw (4, max_uvw) 
      REAL(KIND=PREC_DP) :: uvw_mat (4, 4, max_uvw) 
      LOGICAL lacentric, csym 
      LOGICAL :: success = .TRUE.
      LOGICAL :: grand   = .FALSE.
integer :: iianz
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_NUM     = 1
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'number' /
DATA loname /  6       /
opara  =  (/ '0.0000'  /)   ! Always provide fresh default values
lopara =  (/  6        /)
owerte =  (/  0.0      /)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
if(lpresent(O_NUM)) then
  if(opara(O_NUM) == 'next') then
     chem_ncor = chem_ncor + 1
  else
     iianz = 1
     call ber_params(iianz, opara, lopara, owerte, NOPTIONAL) 
     chem_ncor = nint(owerte(1))
  endif
else
  chem_ncor = max(1,chem_ncor)             ! Old style without "number:" 
endif
   IF(chem_ncor >= CHEM_MAX_COR) THEN      ! Need to allocate more correlations
      n_cor = CHEM_MAX_COR + 10
      CALL alloc_chem_correlation(n_cor)
!     IF(MAXVAL(chem_nvec) > 0) THEN       ! Need to allocate vectors
      if(chem_ndef(CHEM_VEC)>0) then       ! Need to allocate vectors
         n_vec = CHEM_MAX_VEC
         CALL alloc_chem_vec(n_vec, n_cor)
      ENDIF
!     IF(MAXVAL(chem_ncon) > 0) THEN       ! Need to allocate vectors
      if(chem_ndef(CHEM_CON)>0) then       ! Need to allocate Connectivities
         n_con = CHEM_MAX_CON
         CALL alloc_chem_con(n_con, n_cor)
      ENDIF
      IF(MAXVAL(chem_nwin) > 0) THEN       ! Need to allocate vectors
         n_ang = CHEM_MAX_ANG
         CALL alloc_chem_ang(n_ang, n_cor)
      ENDIF
      IF(MAXVAL(chem_nran) > 0) THEN       ! Need to allocate vectors
         n_ran = CHEM_MAX_RAN
         n_atom= MAX(n_atom,MAX_ATOM_ENV)
         CALL alloc_chem_ran(n_ran, n_cor, n_atom)
      ENDIF
      IF(MAXVAL(chem_nenv) > 0) THEN       ! Need to allocate vectors
         n_env = CHEM_MAX_ENV
         n_atom= MAX(n_atom,MAX_ATOM_ENV)
         CALL alloc_chem_env ( n_atom , n_env, n_cor )
      ENDIF
   ENDIF
!
!                                                                       
!                                                                       
      CALL do_cap (cpara (1) ) 
!                                                                       
!------ set neig,add : Add neighbour definition to list                 
!                                                                       
IF(cpara(1)(1:2) == 'AD' .AND. ianz == 1) THEN 
   IF(chem_ncor==1) THEN                   ! This might be erroneous first 'add' command
      IF(MAXVAL(chem_nvec)==0 .AND. MAXVAL(chem_ncon)==0 .AND.         &
         MAXVAL(chem_nwin)==0 .AND.                                    &
         MAXVAL(chem_nran)==0 .AND. MAXVAL(chem_nenv)==0       ) THEN
         RETURN                            ! Ignore erroneous first 'add' command
      ENDIF 
   ENDIF 
   IF(chem_ncor >= CHEM_MAX_COR) THEN      ! Need to allocate more correlations
      n_cor = CHEM_MAX_COR + 10
      CALL alloc_chem_correlation(n_cor)
      IF(MAXVAL(chem_nvec) > 0) THEN       ! Need to allocate vectors
         n_vec = CHEM_MAX_VEC
         CALL alloc_chem_vec(n_vec, n_cor)
      ENDIF
      IF(MAXVAL(chem_ncon) > 0) THEN       ! Need to allocate vectors
         n_con = CHEM_MAX_CON
         CALL alloc_chem_con(n_con, n_cor)
      ENDIF
      IF(MAXVAL(chem_nwin) > 0) THEN       ! Need to allocate vectors
         n_ang = CHEM_MAX_ANG
         CALL alloc_chem_ang(n_ang, n_cor)
      ENDIF
      IF(MAXVAL(chem_nran) > 0) THEN       ! Need to allocate vectors
         n_ran = CHEM_MAX_RAN
         n_atom= MAX(n_atom,MAX_ATOM_ENV)
         CALL alloc_chem_ran(n_ran, n_cor, n_atom)
      ENDIF
      IF(MAXVAL(chem_nenv) > 0) THEN       ! Need to allocate vectors
         n_env = CHEM_MAX_ENV
         n_atom= MAX(n_atom,MAX_ATOM_ENV)
         CALL alloc_chem_env ( n_atom , n_env, n_cor )
      ENDIF
   ENDIF
!
   IF (chem_ncor.lt.CHEM_MAX_COR) THEN 
      chem_ncor = chem_ncor + 1 
   ELSE 
      ier_num = - 12 
      ier_typ = ER_CHEM 
   ENDIF 
!                                                                       
!------ set neig,rese : Reset neighbour list                            
!                                                                       
ELSEIF (cpara (1) (1:2) .eq.'RE') then 
         chem_ncor = 0 
         chem_nvec (1) = 0 
         chem_ncon (1) = 0 
         chem_nwin (1) = 0 
         chem_nran (1) = 0 
         chem_nenv (1) = 0 
!                                                                       
!------ set neig,vec : Define neighbour via vectors                     
!                                                                       
      ELSEIF (cpara (1) (1:2) .eq.'VE'.and.ianz.gt.1) then 
         CALL chem_set_nei_range (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         DO i = 1, ianz 
         iv = nint (werte (i) ) 
         IF (iv.gt.0.and.iv.le.CHEM_MAX_VEC) then 
            IF (chem_cvec (1, iv) .ne. - 9999) then 
               chem_use_vec (i, chem_ncor) = iv 
            ELSE 
               ier_num = - 9 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 9 
            ier_typ = ER_CHEM 
         ENDIF 
         ENDDO 
         IF (ier_num.ne.0) return 
         chem_nvec (chem_ncor) = ianz 
         chem_ctyp (chem_ncor) = CHEM_VEC 
!                                                                       
!------ set neig,con : Define neighbour via connectivity                     
!                                                                       
      ELSEIF (cpara (1) (1:2) .eq.'CO'.and.ianz.gt.1) then 
!                                           ! Try if params match connectivity names
         grand = .TRUE.
         search: DO ii=2, ianz              ! Loop over all params
            i = ii-1                        ! need to subtract 1 as names start in cpara(2)
            success = .FALSE.
            cons: DO j = 1, CHEM_MAX_CON    ! Have do do the full loop exit if -9999 is found
               IF(chem_ccon(1,j)==-9999) EXIT cons   ! no more cons
               IF(cpara(ii)==chem_cname(j)) THEN   ! need exact match
                  iv = j                    ! lets use same variable name as in old style
                  chem_use_con (i, chem_ncor) = iv 
                  success = .TRUE.          ! Found a match for current parameter
                  EXIT cons                 ! we are done for this parameter
               ENDIF
            ENDDO cons
            grand = grand .AND. success
         ENDDO search
         IF(grand) THEN 
            chem_ncon (chem_ncor) = ianz -1
            chem_ctyp (chem_ncor) = CHEM_CON 
            RETURN
         ELSE
!                                           ! try old style
         CALL chem_set_nei_range (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         DO i = 1, ianz 
         iv = nint (werte (i) ) 
         IF (iv.ge.1.and.iv.le.CHEM_MAX_CON) then 
            IF (chem_ccon (1, iv) .ne. - 9999) then 
               chem_use_con (i, chem_ncor) = iv 
            ELSE 
               ier_num = - 9 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 9 
            ier_typ = ER_CHEM 
         ENDIF 
         ENDDO 
         IF (ier_num.ne.0) return 
         chem_ncon (chem_ncor) = ianz 
         chem_ctyp (chem_ncor) = CHEM_CON 
         ENDIF 
!                                                                       
!------ set neig,ran : Define neighbour via ranges                      
!                                                                       
      ELSEIF (cpara (1) (1:2) .eq.'RA'.and.ianz.gt.1) then 
         CALL chem_set_nei_range (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         DO i = 1, ianz 
         iv = nint (werte (i) ) 
         IF (iv.gt.0.and.iv.le.CHEM_MAX_RAN) then 
            IF (chem_cran_uvw (1, 1, iv) .ne. - 9999) then 
               chem_use_ran (i, chem_ncor) = iv 
            ELSE 
               ier_num = - 9 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 9 
            ier_typ = ER_CHEM 
         ENDIF 
         ENDDO 
         IF (ier_num.ne.0) return 
         chem_nran (chem_ncor) = ianz 
         chem_ctyp (chem_ncor) = CHEM_RANGE 
!                                                                       
!------ set neig,ang : Define neighbour via angles                      
!                                                                       
      ELSEIF (cpara (1) (1:2) .eq.'AN'.and.ianz.gt.1) then 
         CALL chem_set_nei_range (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         DO i = 1, ianz 
         iv = nint (werte (i) ) 
         IF (iv.gt.0.and.iv.le.CHEM_MAX_ANG) then 
            IF (chem_cwin (1, iv) .ne. - 9999) then 
               chem_use_win (i, chem_ncor) = iv 
            ELSE 
               ier_num = - 9 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 9 
            ier_typ = ER_CHEM 
         ENDIF 
         ENDDO 
         IF (ier_num.ne.0) return 
         chem_nwin (chem_ncor) = ianz 
         chem_ctyp (chem_ncor) = CHEM_ANG 
!                                                                       
!------ set neig,env : Define neighbour via environment                 
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         
      ELSEIF (cpara (1) (1:3) .eq.'ENV'.and.ianz.gt.1) then 
         CALL chem_set_nei_range (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         DO i = 1, ianz 
         iv = nint (werte (i) ) 
         IF (iv.gt.0.and.iv.le.CHEM_MAX_ENV) then 
            IF (chem_cenv (0, iv) .ne. - 9999) then 
               chem_use_env (i, chem_ncor) = iv 
            ELSE 
               ier_num = - 9 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 9 
            ier_typ = ER_CHEM 
         ENDIF 
         ENDDO 
         IF (ier_num.ne.0) return 
         chem_nenv (chem_ncor) = ianz 
         chem_ctyp (chem_ncor) = CHEM_ENVIR 
!                                                                       
!------ set neig,dis : Define neighbour via distance                    
!                                                                       
      ELSEIF (cpara (1) (1:3) .eq.'DIS'.and.ianz.ge.5) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (ier_num.ne.0) return 
!                                                                       
         chem_cang (chem_ncor) = .false. 
         chem_ctyp (chem_ncor) = CHEM_DIST 
!                                                                       
         csym = .false. 
         IF (ianz.eq.5.or.ianz.eq.6) then 
            CALL do_cap (cpara (ianz) ) 
            csym = (cpara (ianz) (1:1) .eq.'S') 
            ianz = ianz - 1 
         ENDIF 
!                                                                       
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
!                                                                       
         chem_neig (1, 1, chem_ncor) = werte (1) 
         chem_neig (2, 1, chem_ncor) = werte (2) 
         chem_neig (3, 1, chem_ncor) = werte (3) 
         chem_freq_sigma (chem_ncor) = werte (4) 
!                                                                       
         IF (ianz.eq.5) then 
            chem_cang (chem_ncor) = .true. 
            chem_wink_sigma (chem_ncor) = werte (5) 
         ENDIF 
!                                                                       
         DO i = 1, 3 
         u (i) = 0.0 
         v (i) = chem_neig (i, 1, chem_ncor) 
         ENDDO 
         chem_rmax (chem_ncor) = do_blen (.true., u, v) 
         chem_rmin (chem_ncor) = chem_rmax (chem_ncor) -                &
         chem_freq_sigma (chem_ncor) / 2.0                              
         chem_rmax (chem_ncor) = chem_rmax (chem_ncor) +                &
         chem_freq_sigma (chem_ncor) / 2.0                              
!                                                                       
!------ - get symmetrically equivalent directions if needed             
!                                                                       
         IF (csym) then 
            DO i = 1, 3 
            uvw (i, 1) = chem_neig (i, 1, chem_ncor) 
            ENDDO 
            uvw (4, 1) = 0.0 
            DO i = 1, 4 
            DO j = 1, 4 
            uvw_mat (i, j, 1) = 0.0 
            ENDDO 
            uvw_mat (i, i, 1) = 1.0 
            ENDDO 
!                                                                       
            CALL rmc_symmetry (chem_nnei (chem_ncor), uvw, uvw_mat,     &
            max_uvw, .true., lacentric)                                 
!                                                                       
            DO i = 1, chem_nnei (chem_ncor) 
            DO j = 1, 3 
            chem_neig (j, i, chem_ncor) = uvw (j, i) 
            ENDDO 
            ENDDO 
!                                                                       
         ELSE 
            chem_nnei (chem_ncor) = 1 
         ENDIF 
!                                                                       
!------ set neig,dir: sets directions for disp. correlations            
!                                                                       
      ELSEIF (cpara (1) (1:3) .eq.'DIR') then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (cpara (1) (1:3) .eq.'all') then 
            chem_ldall (chem_ncor) = .true. 
         ELSE 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
!                                                                       
            IF (ianz.eq.3) then 
               chem_ldall (chem_ncor) = .false. 
               chem_dir (1, 1, chem_ncor) = werte (1) 
               chem_dir (2, 1, chem_ncor) = werte (2) 
               chem_dir (3, 1, chem_ncor) = werte (3) 
               chem_dir (1, 2, chem_ncor) = werte (1) 
               chem_dir (2, 2, chem_ncor) = werte (2) 
               chem_dir (3, 2, chem_ncor) = werte (3) 
            ELSEIF (ianz.eq.6) then 
               chem_ldall (chem_ncor) = .false. 
               chem_dir (1, 1, chem_ncor) = werte (1) 
               chem_dir (2, 1, chem_ncor) = werte (2) 
               chem_dir (3, 1, chem_ncor) = werte (3) 
               chem_dir (1, 2, chem_ncor) = werte (4) 
               chem_dir (2, 2, chem_ncor) = werte (5) 
               chem_dir (3, 2, chem_ncor) = werte (6) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
!------ unknown subcommand                                              
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_set_neig                  
!*****7*****************************************************************
      SUBROUTINE chem_set_nei_range (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     Parameters for command 'set neig' processed here                  
!-                                                                      
      USE discus_config_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
!                                                                       
      INTEGER i 
      INTEGER istart, iend 
      LOGICAL lrange 
!                                                                       
!                                                                       
      IF (ianz.lt.2) return 
!                                                                       
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      lrange = str_comp (cpara (ianz) , 'range', 2, lpara (ianz) , 5) 
      IF (lrange) then 
         ianz = ianz - 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         istart = nint (werte (1) ) 
         iend = nint (werte (2) ) 
         DO i = 1, iend-istart + 1 
         werte (i) = istart + i - 1 
         ENDDO 
         ianz = iend-istart + 1 
      ELSE 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_set_nei_range             
!*****7*****************************************************************
      SUBROUTINE chem_env (line, laenge) 
!-                                                                      
!     Finds the environment around an atom, site or position            
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE atom_name 
      USE chem_mod 
      USE celltoindex_mod
      USE do_find_mod
      USE modify_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      CHARACTER(9) at_name_i 
      INTEGER lpara (maxw) 
      INTEGER i, ianz, laenge 
      INTEGER iatom, isite, icell (3) 
      REAL(KIND=PREC_DP) ::  werte (maxw), dummy(1)
      REAL(KIND=PREC_DP) ::  pos (3) 
      REAL(KIND=PREC_DP) ::  radius
      LOGICAL latom 
!                                                                       
!                                                                       
      latom = .false. 
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Given atom index                                                
!                                                                       
      IF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) then 
         IF (ianz.eq.3) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            iatom = nint (werte (1) ) 
            IF (iatom.gt.0.and.iatom.le.cr_natoms) then 
               radius = werte (2) 
               CALL indextocell (iatom, icell, isite) 
               DO i = 1, 3 
               pos (i) = cr_pos (i, iatom) 
               ENDDO 
               latom = .true. 
            ELSE 
               ier_num = - 10 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ Given unit cell and site                                        
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'site', 1, lpara (1) , 4) ) then 
         IF (ianz.eq.6) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            DO i = 1, 3 
            icell (i) = nint (werte (i) ) 
            ENDDO 
            isite = nint (werte (4) ) 
            CALL celltoindex (icell, isite, iatom) 
            IF (iatom.gt.0.and.iatom.le.cr_natoms) then 
               radius = werte (5) 
               DO i = 1, 3 
               pos (i) = cr_pos (i, iatom) 
               ENDDO 
               latom = .true. 
            ELSE 
               ier_num = - 10 
               ier_typ = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ Given position                                                  
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'pos', 1, lpara (1) , 3) ) then 
         IF (ianz.eq.5) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            DO i = 1, 3 
            pos (i) = werte (i) 
            ENDDO 
            radius = werte (4) 
            latom = .false. 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ No valid subcommand                                             
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!---------------------------------------------------------------------  
!------ Find neighbours                                                 
!---------------------------------------------------------------------  
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
      dummy = - 1 
      CALL do_find_env (1, dummy, 1, pos, 0.01D0, radius, chem_quick,     &
      chem_period)                                                      
      IF (ier_num.ne.0) return 
!                                                                       
      IF (atom_env (0) .gt.0) then 
         WRITE (output_io, 1000) chem_quick, chem_period 
         WRITE (output_io, 1100) (pos (i), i = 1, 3), radius 
         IF (latom) then 
            at_name_i = at_name (cr_iscat (iatom) ) 
            WRITE (output_io, 1200) (icell (i), i = 1, 3), isite, iatom,&
            at_name_i                                                   
         ENDIF 
         WRITE (output_io, 1300) 
         DO ianz = 1, atom_env (0) 
         at_name_i = at_name (cr_iscat (atom_env (ianz) ) ) 
         CALL indextocell (atom_env (ianz), icell, isite) 
         WRITE (output_io, 1400) ianz, atom_env (ianz), at_name_i,      &
         atom_dis (ianz), (cr_pos (i, atom_env (ianz) ), i = 1, 3),     &
         icell, isite                                                   
         res_para (ianz) = REAL(atom_env (ianz) ) 
         ENDDO 
         res_para (0) = atom_env (0) 
      ENDIF 
!                                                                       
 1000 FORMAT (  ' Found neighbours (quick more = ',L1,                  &
     &          ' periodic boundaries (x,y,z) = ',3(L1,1x),') ')        
 1100 FORMAT (  '   Position of search center   : ',3(F9.3,2X),/        &
     &          '   Search radius [A]           : ',F9.3)               
 1200 FORMAT (  '   Unit cell and site          : ',3(I4,1X),'/',I4,/   &
     &          '   Atom index                  : ',I9,' (',A9,')')     
 1300 FORMAT (/,4X,'#',3X,'atom',4X,'name',5X,'dist [A]',12x,'pos',     &
     &        16x,'cell',5X,'site',/,3X,74('-'))                        
 1400 FORMAT (3X,I3,1X,I6,3X,A9,1X,F7.2,1X,3(F7.2,1X),                  &
     &        1X,3(I3,1X),1X,I3)                                        
!                                                                       
      END SUBROUTINE chem_env                       
!*****7*****************************************************************
      SUBROUTINE chem_bval (line, laenge) 
!-                                                                      
!     Calculates the bond valence sum for specified site/atom           
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE do_find_mod
      USE atom_env_mod 
      USE atom_name 
      USE chem_mod  
      USE celltoindex_mod
      USE modify_mod
      USE element_data_mod
      USE bv_data_mod
!
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      CHARACTER(9) at_name_i 
      CHARACTER (LEN=4)  :: el_name
      INTEGER lpara (maxw) 
      INTEGER i, iii, ianz, laenge 
      INTEGER iatom, isite, icell (3) 
      INTEGER      :: ie_1, ie_2 
      INTEGER      :: ilook  ! Lookup entry in bv_index_table
      REAL(KIND=PREC_DP) ::  werte (maxw), dummy (1)
      REAL(KIND=PREC_DP) ::  pos (3) 
      REAL(KIND=PREC_DP) ::  radius, bval
      LOGICAL latom 
!                                                                       
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Given atom index                                                
!                                                                       
      IF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) then 
         IF (ianz.eq.3) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            iatom = nint (werte (1) ) 
            IF (iatom.gt.0.and.iatom.le.cr_natoms) then 
               radius = werte (2) 
               CALL indextocell (iatom, icell, isite) 
               DO i = 1, 3 
               pos (i) = cr_pos (i, iatom) 
               ENDDO 
               latom = .true. 
            ELSE 
               ier_typ = - 10 
               ier_num = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ Given unit cell and site                                        
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'site', 1, lpara (1) , 4) ) then 
         IF (ianz.eq.6) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            DO i = 1, 3 
            icell (i) = nint (werte (i) ) 
            ENDDO 
            isite = nint (werte (4) ) 
            CALL celltoindex (icell, isite, iatom) 
            IF (iatom.gt.0.and.iatom.le.cr_natoms) then 
               radius = werte (5) 
               DO i = 1, 3 
               pos (i) = cr_pos (i, iatom) 
               ENDDO 
               latom = .true. 
            ELSE 
               ier_typ = - 10 
               ier_num = ER_CHEM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ No valid subcommand                                             
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!---------------------------------------------------------------------  
!------ Find neighbours and calculate BV sum                            
!---------------------------------------------------------------------  
!                                                                       
      ier_num = - 75 
      ier_typ = ER_APPL 
      ie_1    = 0
      el_name = cr_at_lis (cr_iscat (iatom))
      CALL symbf ( el_name, ie_1)
      IF (ier_num.ne.0) return 
!                                                                       
      dummy = - 1 
      CALL do_find_env (1, dummy, 1, pos, 0.01D0, radius, chem_quick,     &
      chem_period)                                                      
      IF (ier_num.ne.0) return 
!                                                                       
      IF (atom_env (0) .gt.0) then 
         bval = 0.0 
         DO ianz = 1, atom_env (0) 
            ier_num = - 75 
            ier_typ = ER_APPL 
            ie_2    = 0
            el_name = cr_at_lis (cr_iscat (iatom))
            CALL symbf ( el_name, ie_2)
            CALL bv_lookup ( ie_1, ie_2, ilook )
            IF (ie_2 == 0 .or. ilook == 0 ) THEN
               ier_num = - 75 
               ier_typ = ER_APPL 
               RETURN 
            ELSE 
               iii  = bv_target(ilook)
               bval = bval + exp ( (bv_r0 (iii) - res_para (ianz) )        &
               / bv_b (iii) )
            ENDIF 
         ENDDO 
         res_para (0) = 1 
         res_para (1) = bval 
         at_name_i = at_name (cr_iscat (iatom) ) 
         WRITE (output_io, 1000) at_name_i, bval 
      ELSE 
         ier_num = - 8 
         ier_typ = ER_CHEM 
      ENDIF 
!                                                                       
 1000 FORMAT    (' Bond valence sum for ',a9,' : ',f7.4) 
!                                                                       
      END SUBROUTINE chem_bval                      
!*****7*****************************************************************
      SUBROUTINE chem_mode (line, laenge) 
!-                                                                      
!     Sets molecule/atom mode                                           
!+                                                                      
      USE discus_config_mod 
      USE chem_mod 
      USE molecule_mod 
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
USE str_comp_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, laenge 
!                                                                       
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (str_comp (cpara (1) , 'mol', 1, lpara (1) , 3) ) then 
         IF (mole_num_mole.gt.0) then 
            chem_sel_atom = .false. 
            WRITE (output_io, 1000) 'molecules' 
         ELSE 
            ier_num = - 20 
            ier_typ = ER_CHEM 
         ENDIF 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) then 
         chem_sel_atom = .true. 
         WRITE (output_io, 1000) 'atoms' 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     ( ' New operation mode : ',a,' ...') 
      END SUBROUTINE chem_mode                      
!*****7*****************************************************************
      SUBROUTINE chem_nei (line, laenge) 
!-                                                                      
!     Finds the environment for given neighbour definition              
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE atom_name 
      USE chem_mod 
      USE molecule_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE lib_f90_allocate_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw! , maxatom 
      PARAMETER (maxw = 3) 
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      CHARACTER(9) at_name_i 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_DP) :: pos (3, 0:MAX_ATOM_ENV) 
      INTEGER lpara (maxw) 
      INTEGER ind (0:MAX_ATOM_ENV), iatom, imol 
      INTEGER i, j, n, ic, ianz, laenge 
      INTEGER :: n_res
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ - Atoms                                                         
!                                                                       
         IF (chem_sel_atom) then 
            iatom = nint (werte (1) ) 
            IF (iatom.le.0.or.iatom.gt.cr_natoms) then 
               ier_num = - 10 
               ier_typ = ER_CHEM 
            ENDIF 
!                                                                       
            ic = nint (werte (2) ) 
            IF (ic.le.0.or.ic.gt.chem_ncor) then 
               ier_num = - 14 
               ier_typ = ER_CHEM 
            ENDIF 
!                                                                       
            IF (chem_ctyp (1) .eq.CHEM_NONE) then 
               ier_num = - 11 
               ier_typ = ER_CHEM 
            ENDIF 
            IF (ier_num.ne.0) return 
!                                                                       
            CALL chem_neighbour (iatom, ic, ind, pos, n, MAX_ATOM_ENV) 
!                                                                       
            IF (n.gt.0) then 
               at_name_i = at_name (cr_iscat (iatom) ) 
               WRITE (output_io, 1000) at_name_i, (cr_pos (i, iatom),   &
               i = 1, 3), iatom, ic                                     
               WRITE (output_io, 1100) 
               DO i = 1, n 
               at_name_i = at_name (cr_iscat (ind (i) ) ) 
               WRITE (output_io, 1200) i, at_name_i, ind (i), (cr_pos ( &
               j, ind (i) ), j = 1, 3)                                  
               ENDDO 
!                                                                       
            ELSE 
               ier_num = - 8 
               ier_typ = ER_CHEM 
            ENDIF 
!                                                                       
!------ - Molecules                                                     
!                                                                       
         ELSE 
            imol = nint (werte (1) ) 
            IF (imol.le.0.or.imol.gt.mole_num_mole) then 
               ier_num = - 63 
               ier_typ = ER_APPL 
            ENDIF 
!                                                                       
            ic = nint (werte (2) ) 
            IF (ic.le.0.or.ic.gt.chem_ncor) then 
               ier_num = - 14 
               ier_typ = ER_CHEM 
            ENDIF 
!                                                                       
            IF (chem_ctyp (1) .eq.CHEM_NONE) then 
               ier_num = - 11 
               ier_typ = ER_CHEM 
            ENDIF 
            IF (ier_num.ne.0) return 
!                                                                       
            CALL chem_neighbour_mol (imol, ic, ind, n, MAX_ATOM_ENV) 
!                                                                       
            IF (n.gt.0) then 
               WRITE (output_io, 2000) (cr_pos (i, mole_cont (mole_off (&
               imol) + 1) ), i = 1, 3), imol, mole_cont (mole_off (imol)&
               + 1), ic                                                 
               WRITE (output_io, 2100) 
               DO i = 1, n 
               WRITE (output_io, 2200) i, mole_type (ind (i) ), ind (i),&
               cr_pos (1, mole_cont (mole_off (ind (i) ) + 1) ),        &
               cr_pos (2, mole_cont (mole_off (ind (i) ) + 1) ),        &
               cr_pos (3, mole_cont (mole_off (ind (i) ) + 1) )         
               ENDDO 
            ELSE 
               ier_num = - 19 
               ier_typ = ER_CHEM 
            ENDIF 
         ENDIF 
!                                                                       
!------ Wrong number of parameters                                      
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!------ store results in res_para                                       
!                                                                       
      IF (n.gt.MAXPAR_RES) then 
         n_res = MAX(n, MAXPAR_RES, MAX_ATOM_ENV)
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
!        ier_typ = ER_CHEM 
!        ier_num = - 2 
!     ELSE 
         res_para (0) = n 
         DO i = 1, n 
         res_para (i) = REAL(ind (i) ) 
         ENDDO 
!     ENDIF 
!                                                                       
 1000 FORMAT (  ' Neighbours found in defined direction ',/             &
     &          '   Center atom name            : ',a9,/                &
     &          '   Position of center atom     : ',3(f8.3,2x),/        &
     &          '   Center atom index           : ',i7,/                &
     &          '   Neighbour definition no.    : ',i7,/)               
 1100 FORMAT (  4x,'#',6x,'name',7x,'index',16x,'position',/            &
     &          3x,(58('-')))                                           
 1200 FORMAT (  2x,i3,3x,a9,3x,i9,4x,3(f8.3,2x)) 
!                                                                       
 2000 FORMAT (  ' Neighbours found in defined direction ',/             &
     &          '   Origin of center molecule   : ',3(f8.3,2x),/        &
     &          '   Center molecule index       : ',i7,/                &
     &          '   Atom index on origin of mol.: ',i7,/                &
     &          '   Neighbour definition no.    : ',i7,/)               
 2100 FORMAT (  4x,'#',6x,'type',7x,'index',17x,'origin',/              &
     &          3x,(58('-')))                                           
 2200 FORMAT (  2x,i3,3x,i9,3x,i9,4x,3(f8.3,2x)) 
      END SUBROUTINE chem_nei                       
!*****7*****************************************************************
      SUBROUTINE chem_homo (line, lp) 
!-                                                                      
!     Check homogeniety of crstal                                       
!+                                                                      
      USE discus_config_mod 
      USE chem_mod 
      USE get_iscat_mod
      USE modify_mod
      USE build_name_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 15) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw), catom (2) 
      REAL(KIND=PREC_DP) :: werte (maxw), wwerte (maxw) 
      INTEGER lpara (maxw), latom (2) 
      INTEGER lp, ianz, iianz 
      LOGICAL locc 
!                                                                       
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.lt.3) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!------ Check concentration                                             
!                                                                       
      IF (str_comp (cpara (1) , 'occ', 1, lpara (1) , 3) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         iianz = 1 
         IF (chem_sel_atom) then 
            CALL get_iscat (iianz, cpara, lpara, werte, maxw, .false.) 
         ELSE 
            iianz = 1 
            CALL ber_params (1, cpara, lpara, werte, maxw) 
         ENDIF 
         IF(ier_num.ne.0) return 
         IF(werte(1) .ge. 0.D0) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, wwerte, maxw, 1) 
            CALL chem_homo_occ (cpara (1), iianz, werte, maxw) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_CHEM 
         ENDIF 
!                                                                       
!------ Check correlations                                              
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'cor', 1, lpara (1) , 3) ) then 
         IF (ianz.ge.5) then 
            CALL do_cap (cpara (2) ) 
            locc = (cpara (2) (1:3) .eq.'OCC') 
            catom (1) = cpara (3) 
            latom (1) = lpara (3) 
            catom (2) = cpara (4) 
            latom (2) = lpara (4) 
            CALL del_params (4, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, wwerte, maxw, 1) 
            CALL chem_homo_corr (cpara (1), catom, latom, locc) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_homo                      
!*****7*****************************************************************
      SUBROUTINE chem_homo_occ (fname, ianz, werte, maxw) 
!-                                                                      
!     Calculates concentration distribution for given atom typ          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE fourier_sup, ONLY: four_csize, four_ranloc
      USE molecule_mod 
      USE modify_func_mod
      USE errlist_mod 
USE lib_length
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER ( * ) fname 
      INTEGER ianz, maxw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      INTEGER itot, ico 
      INTEGER csize (3), lbeg (3) 
      INTEGER i, imol, ibin, il, ia, iup 
      REAL(kind=PREC_DP) :: c, sc, scc, ave_c, sig_c 
!                                                                       
!     LOGICAL atom_allowed, chem_inlot 
!     LOGICAL chem_inlot 
!                                                                       
!------ Some setup                                                      
!                                                                       
      WRITE (output_io, 500) 
      IF (ilots.eq.LOT_OFF) then 
         WRITE (output_io, 600) 
      ELSEIF (ilots.eq.LOT_BOX) then 
         WRITE (output_io, 610) nlots 
         WRITE (output_io, 630) (ls_xyz (i), i = 1, 3), lperiod 
      ELSEIF (ilots.eq.LOT_ELI) then 
         WRITE (output_io, 620) nlots 
         WRITE (output_io, 630) (ls_xyz (i), i = 1, 3), lperiod 
      ENDIF 
!                                                                       
      CALL four_csize (cr_icc, csize, lperiod, ls_xyz) 
!                                                                       
      DO i = 1, CHEM_MAX_BIN 
      chem_hist (i) = 0 
      ENDDO 
!                                                                       
      c = 0.0 
      sc = 0.0 
      scc = 0.0 
!                                                                       
      iup = nlots / 10 
!                                                                       
!------ Loop over all lots                                              
!                                                                       
      DO il = 1, nlots 
      CALL four_ranloc (csize, lbeg) 
!                                                                       
!------ - Get concentration in selected lot                             
!                                                                       
      itot = 0 
      ico = 0 
!                                                                       
!------ - Atom mode                                                     
!                                                                       
      IF (chem_sel_atom) then 
         DO ia = 1, cr_natoms 
         IF (chem_inlot (ia, lbeg) ) then 
            itot = itot + 1 
            IF (atom_allowed (ia, werte, ianz, maxw) ) ico = ico + 1 
         ENDIF 
         ENDDO 
!                                                                       
!------ - Molecule mode                                                 
!                                                                       
      ELSE 
         DO ia = 1, mole_num_mole 
         imol = mole_cont (mole_off (ia) + 1) 
         IF (chem_inlot (imol, lbeg) ) then 
            itot = itot + 1 
            IF (mole_type (ia) .eq.nint (werte (1) ) ) ico = ico + 1 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      IF (itot.ne.0) then 
         c = REAL(ico) / REAL(itot) 
         sc = sc + c 
         scc = scc + c * c 
      ENDIF 
!                                                                       
      ibin = int ( (chem_bin - 1) * c) + 1 
      chem_hist (ibin) = chem_hist (ibin) + 1 
!                                                                       
      IF (nlots.gt.20) then 
         IF (mod (il, iup) .eq.0) WRITE (output_io, 900) il 
      ENDIF 
      ENDDO 
!                                                                       
!------ Write histogram file                                            
!                                                                       
      WRITE (output_io, 1000) fname (1:len_str (fname) ) 
!                                                                       
      CALL oeffne(37, fname, 'unknown') 
      IF (ier_num.ne.0) return 
      DO i = 1, chem_bin 
      WRITE (37, 2000) REAL(i - 1) / REAL(chem_bin - 1), chem_hist (&
      i)                                                                
      ENDDO 
      CLOSE (37) 
!                                                                       
!------ Output average concentration and sigma                          
!                                                                       
      ave_c = sc / REAL(nlots) 
      sig_c = scc / REAL(nlots) - ave_c * ave_c 
!                                                                       
      WRITE (output_io, 1200) ave_c, sig_c 
!                                                                       
  500 FORMAT    (  ' Computing concentration distribution ...') 
  600 FORMAT     ( '   Sample volume            : complete crystal') 
  610 FORMAT     ( '   Sample volume            : ',I4,                 &
     &                      ' box shaped lots')                         
  620 FORMAT     ( '   Sample volume            : ',I4,                 &
     &                      ' ellipsoid shaped lots')                   
  630 FORMAT     ( '   Lot size                 : ',I3,' x ',I3,        &
     &                    ' x ',I3,' unit cells ',/,                    &
     &                    '   Periodic boundaries      : ',L1,/)        
  900 FORMAT    (  '   Finished lot #',i5,' ...') 
 1000 FORMAT    (/,'   Saving histogram to file : ',a,/) 
 1200 FORMAT    (  '   Average concentration    : ',f6.4,' +- ',f6.4) 
 2000 FORMAT    (f8.3,1x,i12) 
!                                                                       
      END SUBROUTINE chem_homo_occ                  
!*****7*****************************************************************
      SUBROUTINE chem_homo_corr (fname, catom, latom, locc) 
!-                                                                      
!     Calculates correlation distribution                               
!+                                                                      
      USE discus_config_mod 
use discus_allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE fourier_sup, ONLY: four_csize, four_ranloc
      USE mc_mod 
      USE errlist_mod 
USE lib_length
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER ( * ) fname, catom (2) 
      INTEGER latom (2) 
      LOGICAL locc 
!                                                                       
      INTEGER csize (3), lbeg (3) 
      INTEGER i, ibin, il, iup, ic, l 
      REAL(KIND=PREC_DP) :: sc (CHEM_MAX_COR), scc (CHEM_MAX_COR) 
      REAL(KIND=PREC_DP) :: ave_c, sig_c 
!                                                                       
!                                                                       
!------ Some setup                                                      
!                                                                       
call alloc_mo_ach(chem_ncor)
      WRITE (output_io, 500) 
      IF (ilots.eq.LOT_OFF) then 
         WRITE (output_io, 600) 
      ELSEIF (ilots.eq.LOT_BOX) then 
         WRITE (output_io, 610) nlots 
         WRITE (output_io, 630) (ls_xyz (i), i = 1, 3), lperiod 
      ELSEIF (ilots.eq.LOT_ELI) then 
         WRITE (output_io, 620) nlots 
         WRITE (output_io, 630) (ls_xyz (i), i = 1, 3), lperiod 
      ENDIF 
!                                                                       
      CALL four_csize (cr_icc, csize, lperiod, ls_xyz) 
!                                                                       
      DO i = 1, CHEM_MAX_BIN 
      chem_hist (i) = 0 
      ENDDO 
!                                                                       
      DO i = 1, CHEM_MAX_COR 
      sc (i) = 0.0 
      scc (i) = 0.0 
      ENDDO 
!                                                                       
      iup = nlots / 10 
!                                                                       
!------ Loop over all lots                                              
!                                                                       
      DO il = 1, nlots 
      CALL four_ranloc (csize, lbeg) 
!                                                                       
!------ - Get correlation in selected lot                               
!                                                                       
      IF (locc) then 
         l = 2 
         IF (chem_sel_atom) then 
            CALL chem_corr_occ (catom, latom, 2, .false., lbeg) 
         ELSE 
            CALL chem_corr_occ_mol (l, catom, latom, 2, .false., lbeg) 
         ENDIF 
      ELSE 
         l = 2 
         IF (chem_sel_atom) then 
            CALL chem_corr_dis (l, catom, latom, 2, .false., lbeg) 
         ELSE 
            CALL chem_corr_dis_mol (l, catom, latom, 2, .false., lbeg) 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) return 
!                                                                       
      DO ic = 1, chem_ncor 
      sc (ic) = sc (ic) + mo_ach_corr (ic) 
      scc (ic) = scc (ic) + mo_ach_corr (ic) **2 
!                                                                       
      ibin = int ( (mo_ach_corr (ic) + 1) * chem_bin / 2.0) + 1 
      chem_hist (ibin) = chem_hist (ibin) + 1 
      ENDDO 
!                                                                       
      IF (nlots.gt.20) then 
         IF (mod (il, iup) .eq.0) WRITE (output_io, 900) il 
      ENDIF 
      ENDDO 
!                                                                       
!------ Write histogram file                                            
!                                                                       
      WRITE (output_io, 1000) fname (1:len_str (fname) ) 
!                                                                       
      CALL oeffne (37, fname, 'unknown') 
      IF (ier_num.ne.0) return 
      DO i = 1, chem_bin 
      WRITE (37, 2000) (2. * REAL(i - 1) / REAL(chem_bin - 1) )     &
      - 1., chem_hist (i)                                               
      ENDDO 
      CLOSE (37) 
!                                                                       
!------ Output average correlations and sigma                           
!                                                                       
      DO ic = 1, chem_ncor 
      ave_c = sc (ic) / REAL(nlots) 
      sig_c = scc (ic) / REAL(nlots) - ave_c * ave_c 
      WRITE (output_io, 1200) ic, ave_c, sig_c 
      ENDDO 
!                                                                       
  500 FORMAT    (  ' Computing correlation distribution ...') 
  600 FORMAT     ( '   Sample volume            : complete crystal') 
  610 FORMAT     ( '   Sample volume            : ',I4,                 &
     &                      ' box shaped lots')                         
  620 FORMAT     ( '   Sample volume            : ',I4,                 &
     &                      ' ellipsoid shaped lots')                   
  630 FORMAT     ( '   Lot size                 : ',I3,' x ',I3,        &
     &                    ' x ',I3,' unit cells ',/,                    &
     &                    '   Periodic boundaries      : ',L1,/)        
  900 FORMAT    (  '   Finished lot #',i5,' ...') 
 1000 FORMAT    (/,'   Saving histogram to file : ',a,/) 
 1200 FORMAT    (  '   Average correlation #',i3,' : ',f6.4,' +- ',f6.4) 
 2000 FORMAT    (f8.3,1x,i12) 
!                                                                       
      END SUBROUTINE chem_homo_corr                 
!*****7*****************************************************************
      LOGICAL function chem_inlot (ia, lbeg) 
!+                                                                      
!     Checks if given atom 'ia' is in given lot 'lbeg'                  
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE diffuse_mod 
use precision_mod
!     USE modify_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP) :: xtest (3), x0 
      INTEGER cr_end 
      INTEGER lbeg (3), iz (3), izmin, izmax 
      INTEGER ia, is, i 
      chem_inlot = .false.
!                                                                       
      cr_end = cr_ncatoms * cr_icc (1) * cr_icc (2) * cr_icc (3)        &
      + 1                                                               
!                                                                       
!------ We are not using lots                                           
!                                                                       
      IF (ilots.eq.LOT_OFF) then 
         chem_inlot = .true. 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (ilots.eq.LOT_BOX) then 
         IF (ia.lt.cr_end) then 
            CALL indextocell (ia, iz, is) 
         ELSE 
            DO i = 1, 3 
            iz (i) = int (cr_pos (i, ia) - cr_dim0 (i, 1) ) + 1 
            ENDDO 
         ENDIF 
!                                                                       
         chem_inlot = .true. 
!                                                                       
         DO i = 1, 3 
         izmin = lbeg (i) 
         izmax = lbeg (i) + ls_xyz (i) - 1 
         chem_inlot = chem_inlot.and. (iz (i) .ge.izmin) .and. (iz (i)  &
         .le.izmax)                                                     
!                                                                       
         IF (.not.chem_inlot.and. (izmax.gt.cr_icc (i) ) ) then 
            izmin = 1 
            izmax = izmax - cr_icc (i) 
            chem_inlot = chem_inlot.and. (iz (i) .ge.izmin) .and. (iz ( &
            i) .le.izmax)                                               
         ENDIF 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (ilots.eq.LOT_ELI) then 
         IF (ia.lt.cr_end) then 
            CALL indextocell (ia, iz, is) 
         ELSE 
            DO i = 1, 3 
            iz (i) = int (cr_pos (i, ia) + cr_dim0 (i, 1) ) - 1 
            ENDDO 
         ENDIF 
!                                                                       
         chem_inlot = .true. 
!                                                                       
         DO i = 1, 3 
         x0 = REAL(ls_xyz (i) ) / 2.0 
         izmin = lbeg (i) 
         izmax = lbeg (i) + ls_xyz (i) - 1 
         xtest (i) = (iz (i) - izmin - x0 + 0.5) **2 / x0**2 
         chem_inlot = chem_inlot.and. (iz (i) .ge.izmin) .and. (iz (i)  &
         .le.izmax)                                                     
!                                                                       
         IF (.not.chem_inlot.and. (izmax.gt.cr_icc (i) ) ) then 
            izmin = 1 
            izmax = izmax - cr_icc (i) 
            chem_inlot = chem_inlot.and. (iz (i) .ge.izmin) .and. (iz ( &
            i) .le.izmax)                                               
         ENDIF 
         ENDDO 
!                                                                       
         chem_inlot = chem_inlot.and. ( (xtest (1) + xtest (2) + xtest (&
         3) ) .le.1.0)                                                  
!                                                                       
      ENDIF 
      END FUNCTION chem_inlot                       
!*****7*****************************************************************
      SUBROUTINE chem_corr_field (line, lp) 
!-                                                                      
!     Calculates correlation field                                      
!+                                                                      
      USE discus_config_mod 
use discus_allocate_appl_mod
      USE chem_mod 
      USE mc_mod 
      USE build_name_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
USE support_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 15) 
      INTEGER maxval 
      PARAMETER (maxval = 1000) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: fname, catom (2) 
      INTEGER back_cvec (5, chem_max_vec) 
      INTEGER lpara (maxw) 
      INTEGER latom (2), lbeg (3), lname 
      INTEGER xmin, xmax, ymin, ymax 
      INTEGER l, i, ix, iy, ianz, lp 
      REAL(KIND=PREC_DP):: werte (maxw)
      REAL(KIND=PREC_DP):: nv (3) 
      REAL(KIND=PREC_DP):: cval (maxval) 
      REAL(KIND=PREC_DP):: back_neig (3, 48, CHEM_MAX_COR) 
      REAL(KIND=PREC_DP):: back_rmin (CHEM_MAX_COR) 
      REAL(KIND=PREC_DP):: back_rmax (CHEM_MAX_COR) 
      LOGICAL locc 
!                                                                       
call alloc_mo_ach(chem_ncor)
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.lt.6) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!------ Get correlation mode                                            
!                                                                       
      IF (str_comp (cpara (1) , 'occ', 1, lpara (1) , 3) ) then 
         locc = .true. 
      ELSEIF (str_comp (cpara (1) , 'dis', 1, lpara (1) , 3) ) then 
         locc = .false. 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!------ Save atom names                                                 
!                                                                       
      catom (1) = cpara (2) (1:lpara (2) ) 
      catom (2) = cpara (3) (1:lpara (3) ) 
      latom (1) = lpara (2) 
      latom (2) = lpara (3) 
!                                                                       
!------ Build filename                                                  
!                                                                       
      CALL del_params (3, ianz, cpara, lpara, maxw) 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      fname = cpara (1) 
      lname = lpara (1) 
!                                                                       
!------ Get calculation range                                           
!                                                                       
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ianz.eq.2) then 
         xmin = nint (werte (1) ) 
         xmax = nint (werte (2) ) 
         ymin = 0 
         ymax = 0 
      ELSEIF (ianz.eq.4) then 
         xmin = nint (werte (1) ) 
         xmax = nint (werte (2) ) 
         ymin = nint (werte (3) ) 
         ymax = nint (werte (4) ) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
!-------------------------------------------------------------------    
!------ Here starts the calculation of the correlation field            
!-------------------------------------------------------------------    
!                                                                       
      CALL oeffne (37, fname, 'unknown') 
      IF (ier_num.ne.0) return 
      WRITE (output_io, 1000) fname (1:lname), catom (1) (1:latom (1) ),&
      catom (2) (1:latom (2) )                                          
!                                                                       
!------ Save current 'neig' settings                                    
!                                                                       
      CALL chem_save_neig (back_cvec, back_neig, back_rmin, back_rmax,  &
      chem_cvec, chem_neig, chem_rmin, chem_rmax, CHEM_MAX_COR, CHEM_MAX_VEC)
!                                                                       
!------ We have a 1D file                                               
!                                                                       
      IF (ymin.eq.0.and.ymax.eq.0) then 
         CALL chem_inc_check (1) 
         IF (ier_num.ne.0) goto 9999 
         lbeg (1) = - 1 
         lbeg (2) = 0 
         lbeg (3) = 0 
!                                                                       
         DO ix = xmin, xmax 
         CALL chem_inc_neig (ix, 0, back_cvec, back_neig, nv) 
         IF (locc) then 
            l = 2 
            IF (chem_sel_atom) then 
               CALL chem_corr_occ (catom, latom, 2, .false., lbeg) 
            ELSE 
               CALL chem_corr_occ_mol (l, catom, latom, 2, .false.,     &
               lbeg)                                                    
            ENDIF 
         ELSE 
            l = 2 
            IF (chem_sel_atom) then 
               CALL chem_corr_dis (l, catom, latom, 2, .false., lbeg) 
            ELSE 
               CALL chem_corr_dis_mol (l, catom, latom, 2, .false.,     &
               lbeg)                                                    
            ENDIF 
         ENDIF 
         IF (ier_num.ne.0) goto 9999 
         WRITE (output_io, 1100) ix, 0, nv, mo_ach_corr (1) 
         WRITE (37, * ) ix, mo_ach_corr (1) 
         ENDDO 
!                                                                       
!------ 2D (NIPL) file                                                  
!                                                                       
      ELSE 
         CALL chem_inc_check (2) 
         IF (ier_num.ne.0) goto 9999 
         lbeg (1) = - 1 
         lbeg (2) = 0 
         lbeg (3) = 0 
!                                                                       
         WRITE (37, * ) (xmax - xmin + 1), (ymax - ymin + 1) 
         WRITE (37, * ) REAL(xmin), REAL(xmax), REAL(ymin),       &
         REAL(ymax)                                                   
         DO iy = ymin, ymax 
         DO ix = xmin, xmax 
         CALL chem_inc_neig (ix, iy, back_cvec, back_neig, nv) 
         IF (locc) then 
            l = 2 
            IF (chem_sel_atom) then 
               CALL chem_corr_occ (catom, latom, 2, .false., lbeg) 
            ELSE 
               CALL chem_corr_occ_mol (l, catom, latom, 2, .false.,     &
               lbeg)                                                    
            ENDIF 
         ELSE 
            l = 2 
            IF (chem_sel_atom) then 
               CALL chem_corr_dis (l, catom, latom, 2, .false., lbeg) 
            ELSE 
               CALL chem_corr_dis_mol (l, catom, latom, 2, .false.,     &
               lbeg)                                                    
            ENDIF 
         ENDIF 
         IF (ier_num.ne.0) goto 9999 
         WRITE (output_io, 1100) ix, iy, nv, mo_ach_corr (1) 
         cval (ix - xmin + 1) = mo_ach_corr (1) 
         ENDDO 
         WRITE (37, * ) (cval (i), i = 1, xmax - xmin + 1) 
         ENDDO 
      ENDIF 
!                                                                       
 9999 CONTINUE 
!                                                                       
!------ Restore 'neig' settings                                         
!                                                                       
      CALL chem_save_neig (chem_cvec, chem_neig, chem_rmin, chem_rmax,  &
      back_cvec, back_neig, back_rmin, back_rmax, CHEM_MAX_COR, CHEM_MAX_VEC)
!                                                                       
      CLOSE (37) 
!                                                                       
 1000 FORMAT (' Calculating correlation field ... ',/                   &
     &        '   Outputfile           : ',a,/,                         &
     &        '   Atoms/Molecules used : ',a4,' - ',a4,//,              &
     &        6x,'x',7x,'y',16x,'Neighbours',12x,'Correlation',/,       &
     &        3x,64('-'))                                               
 1100 FORMAT (3x,i5,3x,i5,4x,3(f8.2,2x),3x,f8.4) 
      END SUBROUTINE chem_corr_field                
!*****7*****************************************************************
      SUBROUTINE chem_save_neig (a_cvec, a_neig, a_rmin, a_rmax, b_cvec,&
      b_neig, b_rmin, b_rmax, CHEM_MAX_COR, CHEM_MAX_VEC)
!-                                                                      
!     Saves/restores current CHEM neighbour settings                    
!+                                                                      
      USE discus_config_mod 
use precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: CHEM_MAX_COR
INTEGER, INTENT(IN) :: CHEM_MAX_VEC
!                                                                       
      INTEGER a_cvec (5, chem_max_vec) 
      INTEGER b_cvec (5, chem_max_vec) 
      REAL(kind=PREC_DP) ::  a_neig (3, 48, CHEM_MAX_COR) 
      REAL(kind=PREC_DP) :: b_neig (3, 48, CHEM_MAX_COR) 
      REAL(kind=PREC_DP) :: a_rmin (CHEM_MAX_COR), a_rmax (CHEM_MAX_COR) 
      REAL(kind=PREC_DP) :: b_rmin (CHEM_MAX_COR), b_rmax (CHEM_MAX_COR) 
!                                                                       
      INTEGER i, j, k 
!                                                                       
      DO i = 1, chem_max_vec 
      DO j = 1, 5 
      a_cvec (j, i) = b_cvec (j, i) 
      ENDDO 
      ENDDO 
!                                                                       
      DO i = 1, CHEM_MAX_COR 
      a_rmin (i) = b_rmin (i) 
      a_rmax (i) = b_rmax (i) 
!                                                                       
      DO j = 1, 3 
      DO k = 1, 48 
      a_neig (j, k, i) = b_neig (j, k, i) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE chem_save_neig                 
!*****7*****************************************************************
      SUBROUTINE chem_inc_neig (ix, iy, cvec, neig, nv) 
!-                                                                      
!     Updates neigbouring distances for correlation field               
!     calculations.                                                     
!+                                                                      
      USE discus_config_mod 
      USE chem_mod 
      USE metric_mod
      USE errlist_mod 
use precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j, ix, iy, iv 
      INTEGER cvec (5, chem_max_vec) 
      REAL(kind=PREC_DP) ::  neig (3, 48, CHEM_MAX_COR) 
      REAL(KIND=PREC_DP) :: nv (3), u (3), v (3) 
!                                                                       
!                                                                       
!------ Mode VECTOR                                                     
!                                                                       
      IF (chem_ctyp (1) .eq.CHEM_VEC) then 
         DO i = 1, 3 
         DO iv = 1, chem_nvec (1) 
         chem_cvec (i + 2, chem_use_vec (iv, 1) ) = ix * cvec (i + 2,   &
         chem_use_vec (iv, 1) ) + iy * cvec (i + 2, chem_use_vec (iv, 2)&
         )                                                              
         ENDDO 
         nv (i) = chem_cvec (i + 2, chem_use_vec (1, 1) ) 
         ENDDO 
!                                                                       
!------ Mode DISTANCE                                                   
!                                                                       
      ELSEIF (chem_ctyp (1) .eq.CHEM_DIST) then 
         DO i = 1, 3 
         DO j = 1, chem_nnei (1) 
         chem_neig (i, j, 1) = REAL(ix) * neig (i, j, 1) + REAL(iy) &
         * neig (i, j, 2)                                               
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         u (i) = 0.0 
         v (i) = chem_neig (i, 1, 1) 
         nv (i) = chem_neig (i, 1, 1) 
         ENDDO 
         chem_rmax (1) = do_blen (.true., u, v) 
         chem_rmin (1) = chem_rmax (1) - chem_freq_sigma (1) / 2.0 
         chem_rmax (1) = chem_rmax (1) + chem_freq_sigma (1) / 2.0 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_inc_neig                  
!*****7*****************************************************************
      SUBROUTINE chem_inc_check (ic) 
!-                                                                      
!     Checks input for correlation field                                
!+                                                                      
      USE discus_config_mod 
      USE chem_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, ic 
!                                                                       
      DO i = 1, ic 
      IF (chem_ctyp (i) .eq.CHEM_NONE) then 
         ier_num = - 11 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
      IF (ic.eq.2.and.chem_ctyp (1) .ne.chem_ctyp (2) ) then 
         ier_num = - 17 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
!                                                                       
      IF (ic.eq.2.and.chem_nvec (1) .ne.chem_nvec (2) ) then 
         ier_num = - 18 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_inc_check                 
!*****7*****************************************************************
      SUBROUTINE chem_disp (ianz, cpara, lpara, werte, maxw, lout) 
!+                                                                      
!     Calculates distortions within the crystal                         
!-                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE atom_name 
      USE chem_mod 
      USE get_iscat_mod
      USE metric_mod
      USE mc_mod 
      USE modify_mod
      USE modify_func_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
!     INTEGER maxatom 
!                                                                       
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL(KIND=PREC_DP):: werte (maxw) 
      LOGICAL lout 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: fname 
      CHARACTER(9) at_name_i, at_name_j 
      INTEGER atom (0:MAX_ATOM_ENV), natom 
      INTEGER iianz, jjanz, i, j, k, is, js, ic 
      INTEGER bl_anz (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) ::  bl_sum (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) ::  bl_s2 (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) ::  patom (3, 0:MAX_ATOM_ENV) 
      REAL(kind=PREC_DP) ::  u (3), v (3), d (3), di 
      REAL(KIND=PREC_DP) :: wwerte (maxw) 
      LOGICAL lfile 
!                                                                       
!     allocate displacement arrays
!
      CALL alloc_chem_disp(CHEM_MAX_COR, MAXSCAT)
!
      fname = cpara (3) (1:lpara (3) ) 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) cpara (1) (1:lpara (1) ), cpara (2)    &
         (1:lpara (2) )                                                 
         IF (ianz.gt.2) then 
            CALL oeffne (37, fname, 'unknown') 
         ENDIF 
      ENDIF 
!                                                                       
      lfile = lout.and. (ianz.gt.2) 
      if(ianz>2) ianz = ianz - 1  ! Set ianz to actual atom type list
!                                                                       
      iianz = 1 
      jjanz = 1 
      CALL get_iscat (iianz, cpara, lpara, werte, maxw, .false.) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL get_iscat (jjanz, cpara, lpara, wwerte, maxw, .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ - reset counters, vectors, ..                                   
!                                                                       
      DO i = 0, cr_nscat 
      DO j = 0, cr_nscat 
      bl_sum (i, j) = 0.0 
      bl_s2 (i, j) = 0.0 
      bl_anz (i, j) = 0 
      ENDDO 
      ENDDO 
!                                                                       
!------ - calculate distortions                                         
!                                                                       
      DO i = 1, cr_natoms 
      IF (atom_allowed (i, werte, ianz, maxw) ) then 
         CALL chem_neighbour (i, ic, atom, patom, natom, MAX_ATOM_ENV) 
         IF (natom.gt.0) then 
            DO j = 1, natom 
            IF (atom_allowed (atom (j), wwerte, ianz, maxw) ) then 
               DO k = 1, 3 
               u (k) = cr_pos (k, i) 
               v (k) = patom (k, j) 
               d (k) = v (k) - u (k) 
               ENDDO 
               di = do_blen (.true., u, v) 
               is = cr_iscat (i) 
               js = cr_iscat (atom (j) ) 
               IF (lfile) WRITE (37, 3000) d, is, js 
               bl_sum (is, js) = bl_sum (is, js) + di 
               bl_s2 (is, js) = bl_s2 (is, js) + di**2 
               bl_anz (is, js) = bl_anz (is, js) + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ - WRITE results and save to res_para block                      
!                                                                       
      DO i = 0, cr_nscat 
      DO j = i, cr_nscat 
      IF (bl_anz (i, j) .ne.0.or.bl_anz (j, i) .ne.0) then 
         chem_disp_ave (ic, i, j) = (bl_sum (i, j) + bl_sum (j, i) )    &
         / (bl_anz (i, j) + bl_anz (j, i) )                             
         chem_disp_sig (ic, i, j) = (bl_s2 (i, j) + bl_s2 (j, i) )      &
         / (bl_anz (i, j) + bl_anz (j, i) )                             
         chem_disp_sig (ic, i, j) = (chem_disp_sig (ic, i, j) - (       &
         chem_disp_ave (ic, i, j) **2) )                                
         IF (chem_disp_sig (ic, i, j) .gt.0) then 
            chem_disp_sig (ic, i, j) = sqrt (chem_disp_sig (ic, i, j) ) 
         ELSE 
            chem_disp_sig (ic, i, j) = 0.0 
         ENDIF 
         IF (lout.and.bl_anz (i, j) .ne.0) then 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            WRITE (output_io, 2000) ic, at_name_i, at_name_j,           &
            chem_disp_ave (ic, i, j), chem_disp_sig (ic, i, j), bl_anz (&
            i, j) + bl_anz (j, i)                                       
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (ier_num.ne.0) return 
      ENDDO 
!                                                                       
      IF (lfile) close (37) 
!                                                                       
 1000 FORMAT (  ' Calculating distortions ',/,                          &
     &          '    Atom types : A = ',A4,' and B = ',A4,' ',//,       &
     &          '    Neig.  Atom A      Atom B       distance',         &
     &          '   sigma     # pairs',/,4x,60('-'))                    
 2000 FORMAT (4x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3000 FORMAT (3(f12.5,1x),3x,2(i3,1x)) 
!                                                                       
      END SUBROUTINE chem_disp                      
!*****7*****************************************************************
      SUBROUTINE chem_disp_mol (ianz, cpara, lpara, werte, maxw, lout) 
!+                                                                      
!     Calculates distortions within the crystal                         
!     Molecules version ..                                              
!-                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE molecule_mod 
      USE ber_params_mod
      USE errlist_mod 
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxmol 
!                                                                       
      PARAMETER (maxmol = 48) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL(KIND=PREC_DP) ::  werte (maxw) 
      LOGICAL lout 
!                                                                       
      INTEGER imol, jmol, i, j, k, it1, it2, is, js, ic 
      INTEGER mol (0:maxmol), nmol 
      INTEGER bl_anz (mole_max_type, mole_max_type) 
      REAL(kind=PREC_DP) :: bl_sum (mole_max_type, mole_max_type) 
      REAL(kind=PREC_DP) :: bl_s2 (mole_max_type, mole_max_type) 
      REAL(KIND=PREC_DP) :: u (3), v (3), d (3), di 
      LOGICAL lfile 
!                                                                       
!     allocate displacement arrays
!
      CALL alloc_chem_disp(CHEM_MAX_COR, MAXSCAT)
!                                                                       
!                                                                       
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ check input                                                     
!                                                                       
      it1 = nint (werte (1) ) 
      it2 = nint (werte (2) ) 
!                                                                       
      IF (it1.le.0.or.it1.gt.mole_num_type.or.it2.le.0.or.it2.gt.mole_nu&
     &m_type) then                                                      
         ier_num = - 64 
         ier_typ = ER_APPL 
      ENDIF 
      IF (ier_num.ne.0) return 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) it1, it2 
         IF (ianz.gt.2) then 
            CALL oeffne (37, cpara (3) , 'unknown') 
         ENDIF 
      ENDIF 
!                                                                       
      lfile = lout.and. (ianz.gt.2) 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ - reset counters, vectors, ..                                   
!                                                                       
      DO i = 1, mole_num_type 
      DO j = 1, mole_num_type 
      bl_sum (i, j) = 0.0 
      bl_s2 (i, j) = 0.0 
      bl_anz (i, j) = 0 
      ENDDO 
      ENDDO 
!                                                                       
!------ - calculate distortions                                         
!                                                                       
      DO i = 1, mole_num_mole 
      imol = mole_cont (mole_off (i) + 1) 
      IF (it1.eq.mole_type (i) .or.it2.eq.mole_type (i) ) then 
         CALL chem_neighbour_mol (i, ic, mol, nmol, maxmol) 
         IF (nmol.gt.0) then 
            DO j = 1, nmol 
            jmol = mole_cont (mole_off (mol (j) ) + 1) 
            IF (it1.eq.mole_type (mol (j) ) .or.it2.eq.mole_type (mol ( &
            j) ) ) then                                                 
               DO k = 1, 3 
               u (k) = cr_pos (k, imol) 
               v (k) = cr_pos (k, jmol) 
               d (k) = u (k) - v (k) 
               ENDDO 
               IF (lfile) WRITE (37, 3000) d 
               di = do_blen (.true., u, v) 
               is = mole_type (i) 
               js = mole_type (mol (j) ) 
               bl_sum (is, js) = bl_sum (is, js) + di 
               bl_s2 (is, js) = bl_s2 (is, js) + di**2 
               bl_anz (is, js) = bl_anz (is, js) + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
      DO i = 1, mole_num_type 
      DO j = i, mole_num_type 
      IF (bl_anz (i, j) .ne.0.or.bl_anz (j, i) .ne.0) then 
         chem_disp_ave (ic, i, j) = (bl_sum (i, j) + bl_sum (j, i) )    &
         / (bl_anz (i, j) + bl_anz (j, i) )                             
         chem_disp_sig (ic, i, j) = (bl_s2 (i, j) + bl_s2 (j, i) )      &
         / (bl_anz (i, j) + bl_anz (j, i) )                             
         chem_disp_sig (ic, i, j) = (chem_disp_sig (ic, i, j) - (       &
         chem_disp_ave (ic, i, j) **2) )                                
         IF (chem_disp_sig (ic, i, j) .gt.0) then 
            chem_disp_sig (ic, i, j) = sqrt (chem_disp_sig (ic, i, j) ) 
         ELSE 
            chem_disp_sig (ic, i, j) = 0.0 
         ENDIF 
         IF (lout) then 
            WRITE (output_io, 2000) ic, i, j, chem_disp_ave (ic, i, j), &
            chem_disp_sig (ic, i, j), bl_anz (i, j) + bl_anz (j, i)     
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (ier_num.ne.0) return 
      ENDDO 
!                                                                       
      IF (lfile) close (37) 
!                                                                       
 1000 FORMAT ( ' Calculating distortions ',/,                           &
     &         '    Molecule types : A = ',I4,' and B = ',I4,' ',//,    &
     &         '    Neig.  mole A      mole B       distance',          &
     &         '   sigma     # pairs',/,4x,60('-'))                     
 2000 FORMAT (4x,i3,3x,i9,3x,i9,5x,f7.3,3x,f7.3,3x,i8) 
 3000 FORMAT (3(f12.5,1x)) 
!                                                                       
      END SUBROUTINE chem_disp_mol                  
!*****7*****************************************************************
      SUBROUTINE chem_angle (ianz, cpara, lpara, werte, uerte, verte,   &
      maxw, lout)                                                       
!+                                                                      
!     Calculates angular distortions within the crystal                 
!-                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE atom_name 
      USE chem_mod 
      USE metric_mod
      USE mc_mod  
      USE modify_mod
      USE modify_func_mod
      USE errlist_mod 
      USE get_iscat_mod
      USE get_params_mod
      USE lib_f90_allocate_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
!     INTEGER maxatom 
!                                                                       
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_DP) :: uerte (maxw) 
      REAL(KIND=PREC_DP) :: verte (maxw) 
      LOGICAL lout 
!                                                                       
      CHARACTER(9) at_name_i, at_name_j 
      CHARACTER(9) name_1, name_2, name_3 
      INTEGER atom (0:MAX_ATOM_ENV), natom 
      INTEGER iianz, i, j, k, is, js, ic 
      INTEGER jjanz, kkanz 
      INTEGER lname_1, lname_2, lname_3 
      INTEGER ba_anz (0:maxscat, 0:maxscat) 
      INTEGER :: n_res
      REAL(kind=PREC_DP) :: ba_sum (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_s2 (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: patom (3, 0:MAX_ATOM_ENV) 
      REAL(kind=PREC_DP) :: u (3), v (3), w (3), wi 
      LOGICAL lfile 
!                                                                       
!     allocate displacement arrays
!
      CALL alloc_chem_disp(CHEM_MAX_COR, MAXSCAT)
!                                                                       
      res_para (0) = 0.0 
!                                                                       
      iianz = 1 
      CALL get_iscat (iianz, cpara, lpara, werte, maxw, .false.) 
      name_1 = cpara (1) 
      lname_1 = lpara (1) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      jjanz = 1 
      CALL get_iscat (jjanz, cpara, lpara, uerte, maxw, .false.) 
      name_2 = cpara (1) 
      lname_2 = lpara (1) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      kkanz = 1 
      CALL get_iscat (kkanz, cpara, lpara, verte, maxw, .false.) 
      name_3 = cpara (1) 
      lname_3 = lpara (1) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) name_1 (1:lname_1), name_2 (1:lname_2),&
         name_3 (1:lname_3)                                             
         IF (ianz.gt.1) then 
            CALL oeffne (37, cpara (2) , 'unknown') 
         ENDIF 
      ENDIF 
!                                                                       
      lfile = lout.and. (ianz.gt.1) 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ - reset counters, vectors, ..                                   
!                                                                       
      DO i = 0, cr_nscat 
      DO j = 0, cr_nscat 
      ba_sum (i, j) = 0.0 
      ba_s2 (i, j) = 0.0 
      ba_anz (i, j) = 0 
      ENDDO 
      ENDDO 
!                                                                       
!------ - calculate angles                                              
!                                                                       
      DO i = 1, cr_natoms 
      IF (atom_allowed (i, werte, iianz, maxw) ) then 
         CALL chem_neighbour (i, ic, atom, patom, natom, MAX_ATOM_ENV) 
         IF (natom.gt.1) then 
            DO k = 1, 3 
            u (k) = patom (k, 1) 
            ENDDO 
            DO j = 2, natom - 1, 2 
            IF (atom (j) .ne.i) then 
               IF (atom_allowed (atom (j), uerte, jjanz, maxw) ) then 
                  DO k = 1, 3 
                  v (k) = patom (k, j) 
                  ENDDO 
                  IF (atom (j + 1) .ne.atom (j) ) then 
                     IF (atom_allowed (atom (j + 1), verte, kkanz, maxw)&
                     ) then                                             
                        DO k = 1, 3 
                        w (k) = patom (k, j + 1) 
                        ENDDO 
                        wi = do_bang (.true., v, u, w) 
                        is = cr_iscat (atom (j) ) 
                        js = cr_iscat (atom (j + 1) ) 
                        IF (lfile) WRITE (37, 3000) u, is, js 
                        ba_sum (is, js) = ba_sum (is, js) + wi 
                        ba_s2 (is, js) = ba_s2 (is, js) + wi**2 
                        ba_anz (is, js) = ba_anz (is, js) + 1 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
      DO i = 0, cr_nscat 
      DO j = i, cr_nscat 
      IF (ba_anz (i, j) .ne.0.or.ba_anz (j, i) .ne.0) then 
         chem_disp_ave (ic, i, j) = (ba_sum (i, j) + ba_sum (j, i) )    &
         / (ba_anz (i, j) + ba_anz (j, i) )                             
         chem_disp_sig (ic, i, j) = (ba_s2 (i, j) + ba_s2 (j, i) )      &
         / (ba_anz (i, j) + ba_anz (j, i) )                             
         chem_disp_sig (ic, i, j) = (chem_disp_sig (ic, i, j) - (       &
         chem_disp_ave (ic, i, j) **2) )                                
         IF (chem_disp_sig (ic, i, j) .gt.0) then 
            chem_disp_sig (ic, i, j) = sqrt (chem_disp_sig (ic, i, j) ) 
         ELSE 
            chem_disp_sig (ic, i, j) = 0.0 
         ENDIF 
         IF (lout) then 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            WRITE (output_io, 2000) ic, at_name_i, at_name_j,           &
            chem_disp_ave (ic, i, j), chem_disp_sig (ic, i, j), ba_anz (&
            i, j) + ba_anz (j, i)                                       
         ENDIF 
!                                                                       
!------ ------- Save results to res[i]                                  
!                                                                       
         IF(res_para(0)+2 > MAXPAR_RES) THEN
            n_res = MAX(NINT(res_para(0)+2), INT(MAXPAR_RES*1.1+10), MAX_ATOM_ENV)
            CALL alloc_param(n_res)
            MAXPAR_RES = n_res
         ENDIF
         res_para (0) = res_para (0) + 2 
         IF (res_para (0) .lt.MAXPAR_RES) then 
            res_para (NINT(res_para (0)) - 1) = chem_disp_ave (ic, i, j)
            res_para (NINT(res_para (0)) ) = chem_disp_sig (ic, i, j) 
         ELSE 
            ier_typ = ER_CHEM 
            ier_num = - 2 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (ier_num.ne.0) return 
      ENDDO 
!                                                                       
      IF (lfile) close (37) 
!                                                                       
 1000 FORMAT ( ' Calculating angles ',/,                                &
     &         '    Zentral atom   = ',A4,/,                            &
     &         '    Atom types : A = ',A4,' and B = ',A4,' ',//,        &
     &         '    Neig.  Atom A      Atom B       angle   ',          &
     &         '   sigma     # pairs',/,4x,60('-'))                     
 2000 FORMAT (4x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3000 FORMAT (3(f12.5,1x),3x,2(i3,1x)) 
!                                                                       
      END SUBROUTINE chem_angle                     
!*****7*****************************************************************
      SUBROUTINE chem_corr_dis (ianz, cpara, lpara, maxw, lout, lbeg) 
!+                                                                      
!     Calculates displacement correlations within the crystal           
!       according to: cij = <x(i)x(j)>/sqrt(<x(i)**2><x(j)**2>)         
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE chem_mod 
      USE chem_aver_mod
      USE celltoindex_mod
      USE get_iscat_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod   
      USE modify_mod   
      USE modify_func_mod
      USE errlist_mod 
      USE get_params_mod
      USE lib_f90_allocate_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER                             ,INTENT(INOUT) :: ianz 
      INTEGER                             ,INTENT(IN   ) :: maxw 
      CHARACTER (LEN=*), DIMENSION(1:maxw), INTENT(IN) ::  cpara
      INTEGER          , DIMENSION(1:maxw), INTENT(IN) ::  lpara
      INTEGER          , DIMENSION(1:3   ), INTENT(IN) ::  lbeg 
      LOGICAL                             , INTENT(IN) ::  lout 
!                                                                       
      INTEGER maxww!, maxatom 
!                                                                       
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
!     INTEGER ianz, maxw 
!     CHARACTER ( * ) cpara (maxw) 
!     INTEGER lpara (maxw), lbeg (3) 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))), DIMENSION(:), ALLOCATABLE :: ccpara
      INTEGER         , DIMENSION(:), ALLOCATABLE :: llpara
!     LOGICAL lout
!                                                                       
      INTEGER atom (0:MAX_ATOM_ENV), natom 
      INTEGER icc (3), jcc (3) 
      INTEGER i, j, ii, is, js, ic, iianz, jjanz, nn 
      INTEGER :: n_res
      REAL(kind=PREC_DP) :: patom (3, 0:MAX_ATOM_ENV) 
      REAL(KIND=PREC_DP):: werte (MAXSCAT), wwerte (MAXSCAT) 
      REAL(KIND=PREC_DP):: idir (3), jdir (3), di (3), dj (3) 
      REAL(kind=PREC_DP) :: rdi, rdj, dpi, dpj 
      REAL(kind=PREC_DP) :: xij, xi2, xj2 
      LOGICAL lvalid 
      LOGICAL lauto     ! catch autocorrelation of atom to itself
!                                                                       
      lauto = .true.
!                                                                       
!------ writing output line                                             
!                                                                       
      maxww = MAXSCAT
      IF (lout) then 
         WRITE (output_io, 1000) cpara(1)(1:lpara(1)), &
                                 cpara(2)(1:lpara (2))                                                 
      ENDIF 
!
      ALLOCATE(ccpara(1:maxw))
      ALLOCATE(llpara(1:maxw))
      llpara(1) = MAX(1,MIN(lpara(1),4))
      llpara(2) = MAX(1,MIN(lpara(2),4))
      ccpara(1) = cpara(1)(1:llpara(1))
      ccpara(2) = cpara(2)(1:llpara(2))
!
!                                                                       
      iianz = 1 
      jjanz = 1 
      CALL get_iscat (iianz, ccpara, llpara, werte, maxww, .false.) 
      CALL del_params (1, ianz, ccpara, llpara, maxw) 
      CALL get_iscat (jjanz, ccpara, llpara, wwerte, maxww, .false.) 
      IF (ier_num.ne.0) THEN
         DEALLOCATE(ccpara)
         DEALLOCATE(llpara)
         return 
      ENDIF
!                                                                       
      CALL chem_aver (.false., .true.) 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      corr_loop: DO ic = 1, chem_ncor 
!                                                                       
!DBGDISPwrite(*,*) ' DIR ', chem_dir(1, 1, ic), .not.chem_ldall(ic), chem_ldall(ic)
         IF (chem_dir(1, 1, ic) ==  -9999..and..not.chem_ldall(ic)) THEN
            ier_num = - 15 
            ier_typ = ER_CHEM 
            DEALLOCATE(ccpara)
            DEALLOCATE(llpara)
            RETURN 
         ENDIF 
!                                                                       
         nn = 0 
         xij = 0.0 
         xi2 = 0.0 
         xj2 = 0.0 
         DO i = 1, 3 
            idir (i) = chem_dir (i, 1, ic) 
            jdir (i) = chem_dir (i, 2, ic) 
         ENDDO 
!                                                                       
!------ - calculate correlations                                        
!                                                                       
         rdi = skalpro (idir, idir, cr_gten) 
         rdj = skalpro (jdir, jdir, cr_gten) 
         IF (rdi.gt.0.0) rdi = sqrt (rdi) 
         IF (rdj.gt.0.0) rdj = sqrt (rdj) 
!DBGDISPwrite(*,*) ' CHEM CORR idir ', idir(:), rdi, chem_ldall (ic)
!DBGDISPwrite(*,*) ' CHEM CORR jdir ', jdir(:), rdj
!                                                                       
!DBGDISPwrite(*,*) ' ACCUM  , dj(1),           dpj,             dpj**2,                 nn, xij,              rdj'
         atom_loop: DO i = 1, cr_natoms 
!                                                                       
!------ --- Check if selected atom is valid                             
!                                                                       
      lvalid = atom_allowed (i, werte, iianz, maxww) 
      IF (lbeg (1) .gt.0) then 
         lvalid = lvalid.and.chem_inlot (i, lbeg) 
      ENDIF 
!                                                                       
      IF (lvalid) then 
         CALL chem_neighbour (i, ic, atom, patom, natom, MAX_ATOM_ENV) 
         IF (natom.gt.0) then 
            CALL indextocell (i, icc, is) 
            DO j = 1, 3 
               di(j) = cr_pos(j,i) - chem_ave_pos(j,is) -&
                       REAL(icc(j) - 1) - cr_dim0 (j,1)                                    
            ENDDO 
!                                                                       
            IF (chem_ldall (ic) ) then 
               DO j = 1, 3 
                  jdir (j) = di (j) 
               ENDDO 
               rdj = skalpro (jdir, jdir, cr_gten) 
               IF (rdj.gt.0.0) then 
                  rdj = sqrt (rdj) 
               ELSE 
                  rdj = 1.0 
               ENDIF 
               dpi = 1.0 
            ELSE 
               dpi = skalpro (di, idir, cr_gten) / rdi 
            ENDIF 
!                                                                       
            DO j = 1, natom 
               lauto = lauto .AND. cr_iscat(i)==cr_iscat(atom(j))
            IF (atom_allowed (atom (j), wwerte, jjanz, maxww) ) then 
               CALL indextocell (atom (j), jcc, js) 
               DO ii = 1, 3 
                  dj(ii) = cr_pos(ii,atom (j) ) - chem_ave_pos (ii,js) &
                          - REAL(jcc(ii) - 1) - cr_dim0(ii,1)                 
               ENDDO 
               dpj = skalpro (dj, jdir, cr_gten) / rdj 
               xij = xij + dpi * dpj 
               xi2 = xi2 + dpi**2 
               xj2 = xj2 + dpj**2 
               nn = nn + 1 
!DBGDISPwrite(*,*) ' ACCUM ', dj(1), dpj, dpj**2, nn, xij, rdj
            ENDIF 
            ENDDO 
         ENDIF 
         ENDIF 
         ENDDO atom_loop
!DBGDISPwrite(*,*) ' RDI RDJ ', rdi, rdj
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
!DBGDISPwrite(*,*) ' FINAL ', xij, xi2, xj2, xij / sqrt (xi2 * xj2)
         IF (nn.ne.0) then 
            xij = xij / REAL(nn) 
            xi2 = xi2 / REAL(nn) 
            xj2 = xj2 / REAL(nn) 
!                                                                       
!DBGDISPwrite(*,*) ' FINAL ', xij, xi2, xj2, xij / sqrt (xi2 * xj2)
            IF (xi2.ne.0.and.xj2.ne.0.0) then 
               mo_ach_corr (ic) = xij / sqrt (xi2 * xj2) 
            ELSE 
               mo_ach_corr (ic) = 0.0 
            ENDIF 
         ELSEIF(lauto) THEN   ! Got autocorrelation of atom to itself
           mo_ach_corr (ic) = 1.0
         ENDIF
         IF (nn>0 .or. (nn==0.AND.lauto)) THEN
!                                                                       
            IF (lout) then 
               IF (chem_ldall (ic) ) then 
                  WRITE(output_io,1100) ic,nn,mo_ach_corr(ic) 
               ELSE 
                  WRITE(output_io,1110) ic,idir,jdir,nn,mo_ach_corr(ic)                                                      
               ENDIF 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 8 
            ier_typ = ER_CHEM 
         ENDIF 
         IF (ier_num.ne.0) THEN
            DEALLOCATE(ccpara)
            DEALLOCATE(llpara)
            RETURN
         ENDIF
      ENDDO   corr_loop
!                                                                       
!------ Save results to res[i]                                          
!                                                                       
      IF(chem_ncor     > MAXPAR_RES) THEN
         n_res = MAX(chem_ncor ,MAXPAR_RES, MAX_ATOM_ENV)
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
      res_para (0) = chem_ncor 
      DO ic = 1, chem_ncor 
         IF (res_para (0) <=  MAXPAR_RES) then 
            res_para (ic) = mo_ach_corr (ic) 
         ELSE 
            ier_typ = ER_CHEM 
            ier_num = - 2 
         ENDIF 
      ENDDO 
!
      DEALLOCATE(ccpara)
      DEALLOCATE(llpara)
!                                                                       
 1000 FORMAT     (  ' Calculating correlations ',/,                     &
     &          '    Atom type: A = ',A4,' B = ',A4,                    &
     &        //,4x,'Neig.',5x,'Displacement A',6x,'Displacement B',    &
     &           6x,'# pairs',5x,'correlation',/,4x,73('-'))            
 1100 FORMAT     (5x,i3,16x,'all directions',15x,i8,5x,f7.4) 
 1110 FORMAT     (5x,i3,4x,3(f5.2,1x),2x,3(f5.2,1x),3x,i8,5x,f7.4) 
!                                                                       
      END SUBROUTINE chem_corr_dis                  
!*****7*****************************************************************
      SUBROUTINE chem_corr_dis_mol (ianz, cpara, lpara, maxw, lout,     &
      lbeg)                                                             
!+                                                                      
!     Calculates displacement correlations within the crystal           
!       according to: cij = <x(i)x(j)>/sqrt(<x(i)**2><x(j)**2>)         
!     Molecule version ...                                              
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE chem_mod 
      USE chem_aver_mod
      USE celltoindex_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod   
      USE modify_mod   
      USE molecule_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE lib_f90_allocate_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxww, maxmol 
!                                                                       
      PARAMETER (maxww = 3) 
      PARAMETER (maxmol = 48) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw), lbeg (3) 
      LOGICAL lout 
!                                                                       
      INTEGER mol (0:maxmol), nmol 
      INTEGER icc (3), jcc (3) 
      INTEGER i, j, ii, is, js, ic, nn 
      INTEGER it1, it2, imol, jmol 
      INTEGER :: n_res
      REAL(KIND=PREC_DP) :: werte (maxww) 
      REAL(KIND=PREC_DP) ::  idir (3), jdir (3), di (3), dj (3) 
      REAL(kind=PREC_DP) :: rdi, rdj, dpi, dpj 
      REAL(kind=PREC_DP) :: xij, xi2, xj2 
      LOGICAL lvalid 
!                                                                       
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ check input                                                     
!                                                                       
      it1 = nint (werte (1) ) 
      it2 = nint (werte (2) ) 
!                                                                       
      IF (it1.le.0.or.it1.gt.mole_num_type.or.it2.le.0.or.it2.gt.mole_nu&
     &m_type) then                                                      
         ier_num = - 64 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      CALL chem_aver (.false., .true.) 
!                                                                       
!------ writing output line                                             
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) it1, it2 
      ENDIF 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
      IF (chem_dir (1, 1, ic) .eq. - 9999..and..not.chem_ldall (ic) )   &
      then                                                              
         ier_num = - 15 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
!                                                                       
      nn = 0 
      xij = 0.0 
      xi2 = 0.0 
      xj2 = 0.0 
      DO i = 1, 3 
      idir (i) = chem_dir (i, 1, ic) 
      jdir (i) = chem_dir (i, 2, ic) 
      ENDDO 
!                                                                       
!------ - calculate correlations                                        
!                                                                       
      rdi = skalpro (idir, idir, cr_gten) 
      rdj = skalpro (jdir, jdir, cr_gten) 
      IF (rdi.gt.0.0) rdi = sqrt (rdi) 
      IF (rdj.gt.0.0) rdj = sqrt (rdj) 
!                                                                       
      DO i = 1, mole_num_mole 
      imol = mole_cont (mole_off (i) + 1) 
!                                                                       
!------ --- Check if selected atom is valid                             
!                                                                       
      lvalid = mole_type (i) .eq.it1.or.mole_type (i) .eq.it2 
      IF (lbeg (1) .gt.0) then 
         lvalid = lvalid.and.chem_inlot (imol, lbeg) 
      ENDIF 
!                                                                       
      IF (lvalid) then 
         CALL chem_neighbour_mol (i, ic, mol, nmol, maxmol) 
         IF (nmol.gt.0) then 
            CALL indextocell (imol, icc, is) 
            DO j = 1, 3 
            di (j) = cr_pos (j, imol) - chem_ave_pos (j, is) - REAL(  &
            icc (j) - 1) - cr_dim0 (j, 1)                               
            ENDDO 
!                                                                       
            IF (chem_ldall (ic) ) then 
               DO j = 1, 3 
               jdir (j) = di (j) 
               ENDDO 
               rdj = skalpro (jdir, jdir, cr_gten) 
               IF (rdj.gt.0.0) then 
                  rdj = sqrt (rdj) 
               ELSE 
                  rdj = 1.0 
               ENDIF 
               dpi = 1.0 
            ELSE 
               dpi = skalpro (di, idir, cr_gten) / rdi 
            ENDIF 
!                                                                       
            DO j = 1, nmol 
            IF (mole_type (mol (j) ) .eq.it1.or.mole_type (mol (j) )    &
            .eq.it2) then                                               
               jmol = mole_cont (mole_off (mol (j) ) + 1) 
               CALL indextocell (jmol, jcc, js) 
               DO ii = 1, 3 
               dj (ii) = cr_pos (ii, jmol) - chem_ave_pos (ii, js)      &
               - REAL(jcc (ii) - 1) - cr_dim0 (ii, 1)                 
               ENDDO 
               dpj = skalpro (dj, jdir, cr_gten) / rdj 
               xij = xij + dpi * dpj 
               xi2 = xi2 + dpi**2 
               xj2 = xj2 + dpj**2 
               nn = nn + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
      IF (nn.ne.0) then 
         xij = xij / REAL(nn) 
         xi2 = xi2 / REAL(nn) 
         xj2 = xj2 / REAL(nn) 
!                                                                       
         IF (xi2.ne.0.and.xj2.ne.0.0) then 
            mo_ach_corr (ic) = xij / sqrt (xi2 * xj2) 
         ELSE 
            mo_ach_corr (ic) = 0.0 
         ENDIF 
!                                                                       
         IF (lout) then 
            IF (chem_ldall (ic) ) then 
               WRITE (output_io, 1100) ic, nn, mo_ach_corr (ic) 
            ELSE 
               WRITE (output_io, 1110) ic, idir, jdir, nn, mo_ach_corr (&
               ic)                                                      
            ENDIF 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 8 
         ier_typ = ER_CHEM 
      ENDIF 
      IF (ier_num.ne.0) return 
      ENDDO 
!                                                                       
!------ Save results to res[i]                                          
!                                                                       
      IF(chem_ncor     > MAXPAR_RES) THEN
         n_res = MAX(chem_ncor ,MAXPAR_RES, MAX_ATOM_ENV)
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
      res_para (0) = chem_ncor 
      DO ic = 1, chem_ncor 
      IF (res_para (0) <=  MAXPAR_RES) then 
         res_para (ic) = mo_ach_corr (ic) 
      ELSE 
         ier_typ = ER_CHEM 
         ier_num = - 2 
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT     (  ' Calculating correlations ',/,                     &
     &          '    Molecule type: A = ',I4,' B = ',I4,                &
     &        //,4x,'Neig.',5x,'Displacement A',6x,'Displacement B',    &
     &           6x,'# pairs',5x,'correlation',/,4x,73('-'))            
 1100 FORMAT     (5x,i3,16x,'all directions',15x,i8,5x,f7.4) 
 1110 FORMAT     (5x,i3,4x,3(f5.2,1x),2x,3(f5.2,1x),3x,i8,5x,f7.4) 
!                                                                       
      END SUBROUTINE chem_corr_dis_mol              
!*****7*****************************************************************
      SUBROUTINE chem_corr_occ (cpara, lpara, maxw, lout, lbeg) 
!+                                                                      
!     Calculates occupational correlations within the crystal           
!       according to: cij = (Pij-T**2)/T(1-T).                          
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
use discus_allocate_appl_mod
USE atom_env_mod
      USE chem_mod 
      USE get_iscat_mod
      USE mc_mod 
      USE mmc_mod   
      USE modify_mod
      USE modify_func_mod
!
      USE errlist_mod 
      USE lib_f90_allocate_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxww!, maxatom 
!                                                                       
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw), lbeg (3) 
      LOGICAL lout 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: ccpara (MAXSCAT) 
      INTEGER llpara (maxw) 
!                                                                       
      INTEGER atom (0:MAX_ATOM_ENV), natom 
      INTEGER i, j, is, js, ic, iianz, jjanz 
      INTEGER nneig 
      INTEGER :: n_res
      REAL(kind=PREC_DP) :: patom (3, 0:MAX_ATOM_ENV) 
      REAL(KIND=PREC_DP) :: werte (MAXSCAT), wwerte (MAXSCAT) 
      REAL(kind=PREC_DP) :: pneig (2, 2) 
      REAL(kind=PREC_DP) :: pro00, pro01, pro11 
      REAL(kind=PREC_DP) :: thet 
      LOGICAL lvalid 
!                                                                       
!     LOGICAL atom_allowed, chem_inlot 
!     LOGICAL chem_inlot 
!                                                                       
      maxww = MAXSCAT
      iianz = 1 
      ccpara (1) = cpara (1) 
      llpara (1) = lpara (1) 
      CALL get_iscat (iianz, ccpara, llpara, werte, maxww, .false.) 
!                                                                       
      jjanz = 1 
      ccpara (1) = cpara (2) 
      llpara (1) = lpara (2) 
      CALL get_iscat (jjanz, ccpara, llpara, wwerte, maxww, .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ check input                                                     
!                                                                       
      IF (int (werte (1) ) .eq. - 1.or.int (wwerte (1) ) .eq. - 1) then 
         ier_num = - 6 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
!                                                                       
      IF (int (werte (1) ) .eq.int (wwerte (1) ) ) then 
         ier_num = - 7 
         ier_typ = ER_CHEM 
         RETURN 
      ENDIF 
!                                                                       
!------ write output                                                    
!                                                                       
      is = int (werte (1) ) 
      js = int (wwerte (1) ) 
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) cr_at_lis (is), cr_at_lis (js) 
      ENDIF 
call alloc_mo_ach(chem_ncor)
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ - reset counters, vectors, ..                                   
!                                                                       
      DO i = 1, 2 
      DO j = 1, 2 
      pneig (i, j) = 0.0 
      ENDDO 
      ENDDO 
!                                                                       
!------ - calculate correlations                                        
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!------ --- Check if selected atom is valid                             
!                                                                       
      lvalid = .true. 
      IF (lbeg (1) .gt.0) then 
         lvalid = lvalid.and.chem_inlot (i, lbeg) 
      ENDIF 
!                                                                       
      IF (lvalid) then 
         is = - 1 
         IF (atom_allowed (i, werte, iianz, maxww) ) is = 1 
         IF (atom_allowed (i, wwerte, jjanz, maxww) ) is = 2 
         IF (is.gt.0) then 
            CALL chem_neighbour (i, ic, atom, patom, natom, MAX_ATOM_ENV) 
            IF (natom.gt.0) then 
               DO j = 1, natom 
               js = - 1 
               IF (atom_allowed (atom (j), werte, iianz, maxww) ) js =  &
               1                                                        
               IF (atom_allowed (atom (j), wwerte, jjanz, maxww) ) js = &
               2                                                        
               IF (js.gt.0) then 
                  pneig (is, js) = pneig (is, js) + 1.0 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
      nneig = nint (pneig (1, 1) + pneig (2, 2) + pneig (1, 2) + pneig (&
      2, 1) )                                                           
      IF (nneig.ne.0) then 
         pro00 = pneig (1, 1) / REAL(nneig) 
         pro01 = (pneig (1, 2) + pneig (2, 1) ) / REAL(nneig) 
         pro11 = pneig (2, 2) / REAL(nneig) 
         thet = 0.5 * (2.0 * pneig (1, 1) + pneig (1, 2) + pneig (2, 1) &
         ) / REAL(nneig)                                              
         IF (thet.ne.0.and.thet.ne.1.0) then 
            mo_ach_corr (ic) = (pro00 - thet**2) / (thet * (1 - thet) ) 
         ELSE 
            mo_ach_corr (ic) = 0.0 
         ENDIF 
!                                                                       
         IF (lout) WRITE (output_io, 1100) ic, 100. * pro00, 100. *     &
         pro01, 100. * pro11, nneig, mo_ach_corr (ic)                   
!                                                                       
      ELSE 
         IF (lout) WRITE (output_io, 1200) ic
!         ier_num = - 8 
!         ier_typ = ER_CHEM 
      ENDIF 
      IF (ier_num.ne.0) return 
      ENDDO 
!                                                                       
!------ Save results to res[i]                                          
!                                                                       
      IF(chem_ncor     > MAXPAR_RES) THEN
         n_res = MAX(chem_ncor ,MAXPAR_RES, MAX_ATOM_ENV)
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
      res_para (0) = chem_ncor 
      DO ic = 1, chem_ncor 
      IF (res_para (0) <=  MAXPAR_RES) then 
         res_para (ic) = mo_ach_corr (ic) 
      ELSE 
         ier_typ = ER_CHEM 
         ier_num = - 2 
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT ( ' Calculating correlations ',/,                          &
     &         '    Atom types : A = ',A4,' and B = ',A4,' ',//,        &
     &         4x,'Neig.',5x,'AA',9x,'AB',9x,'BB',9x,                   &
     &         '# pairs    correlation',/,4x,65('-'))                   
 1100 FORMAT (5x,i3,3x,3(f6.2,' % ',2x),1x,i8,6x,f7.4) 
 1200 FORMAT (5x,i3,4x,'No pairs in this neighborhood')
!                                                                       
      END SUBROUTINE chem_corr_occ                  
!*****7*****************************************************************
      SUBROUTINE chem_corr_occ_mol (ianz, cpara, lpara, maxw, lout,     &
      lbeg)                                                             
!+                                                                      
!     Calculates occupational correlations within the crystal           
!       according to: cij = (Pij-T**2)/T(1-T).                          
!     Molecules version ...                                             
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE chem_mod 
      USE mc_mod 
      USE mmc_mod   
      USE molecule_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE lib_f90_allocate_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxmol, maxww 
      PARAMETER (maxmol = 48) 
      PARAMETER (maxww = 3) 
!                                                                       
      INTEGER ianz, maxw 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw), lbeg (3) 
      LOGICAL lout 
!                                                                       
      INTEGER mol (0:maxmol), nmol 
      INTEGER i, j, is, js, ic 
      INTEGER imol, it1, it2 
      INTEGER nneig 
      INTEGER :: n_res
      REAL(KIND=PREC_DP) :: werte (maxww) 
      REAL(kind=PREC_DP) :: pneig (2, 2) 
      REAL(kind=PREC_DP) :: pro00, pro01, pro11 
      REAL(kind=PREC_DP) :: thet 
      LOGICAL lvalid 
!                                                                       
!     LOGICAL chem_inlot 
!                                                                       
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ check input                                                     
!                                                                       
      IF (ianz.eq.2) then 
         it1 = nint (werte (1) ) 
         it2 = nint (werte (2) ) 
!                                                                       
      IF (it1.le.0.or.it1.gt.mole_num_type.or.it2.le.0.or.it2.gt.mole_nu&
     &m_type) then                                                      
            ier_num = - 64 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
         IF (it1.eq.it2) then 
            ier_num = - 21 
            ier_typ = ER_CHEM 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      IF (ier_num.ne.0) return 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) it1, it2 
      ENDIF 
!                                                                       
!------ loop over all defined correlations                              
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ - reset counters, vectors, ..                                   
!                                                                       
      DO i = 1, 2 
      DO j = 1, 2 
      pneig (i, j) = 0.0 
      ENDDO 
      ENDDO 
!                                                                       
!------ - calculate correlations                                        
!                                                                       
      DO i = 1, mole_num_mole 
      imol = mole_cont (mole_off (i) + 1) 
!                                                                       
!------ --- Check if selected atom is valid                             
!                                                                       
      lvalid = .true. 
      IF (lbeg (1) .gt.0) then 
         lvalid = lvalid.and.chem_inlot (imol, lbeg) 
      ENDIF 
!                                                                       
      IF (lvalid) then 
         is = - 1 
         IF (mole_type (i) .eq.it1) is = 1 
         IF (mole_type (i) .eq.it2) is = 2 
         IF (is.gt.0) then 
            CALL chem_neighbour_mol (i, ic, mol, nmol, maxmol) 
            IF (nmol.gt.0) then 
               DO j = 1, nmol 
               js = - 1 
               IF (mole_type (mol (j) ) .eq.it1) js = 1 
               IF (mole_type (mol (j) ) .eq.it2) js = 2 
               IF (js.gt.0) then 
                  pneig (is, js) = pneig (is, js) + 1.0 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ - write results and save to res_para block                      
!                                                                       
      nneig = nint (pneig (1, 1) + pneig (2, 2) + pneig (1, 2) + pneig (&
      2, 1) )                                                           
      IF (nneig.ne.0) then 
         pro00 = pneig (1, 1) / REAL(nneig) 
         pro01 = (pneig (1, 2) + pneig (2, 1) ) / REAL(nneig) 
         pro11 = pneig (2, 2) / REAL(nneig) 
         thet = 0.5 * (2.0 * pneig (1, 1) + pneig (1, 2) + pneig (2, 1) &
         ) / REAL(nneig)                                              
         IF (thet.ne.0.and.thet.ne.1.0) then 
            mo_ach_corr (ic) = (pro00 - thet**2) / (thet * (1 - thet) ) 
         ELSE 
            mo_ach_corr (ic) = 0.0 
         ENDIF 
!                                                                       
         IF (lout) WRITE (output_io, 1100) ic, 100. * pro00, 100. *     &
         pro01, 100. * pro11, nneig, mo_ach_corr (ic)                   
!                                                                       
      ELSE 
         IF (lout) WRITE (output_io, 1200) ic
!         ier_num = - 19 
!         ier_typ = ER_CHEM 
      ENDIF 
      IF (ier_num.ne.0) return 
      ENDDO 
!                                                                       
!------ Save results to res[i]                                          
!                                                                       
      IF(chem_ncor     > MAXPAR_RES) THEN
         n_res = MAX(chem_ncor ,MAXPAR_RES, MAX_ATOM_ENV)
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
      res_para (0) = chem_ncor 
      DO ic = 1, chem_ncor 
      IF (res_para (0) <=  MAXPAR_RES) then 
         res_para (ic) = mo_ach_corr (ic) 
      ELSE 
         ier_typ = ER_CHEM 
         ier_num = - 2 
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT (  ' Calculating correlations ',/,                         &
     &      '    Molecule types : A = ',I4,' and B = ',I4,' ',//,       &
     &      4x,'Neig.',5x,'AA',9x,'AB',9x,'BB',9x,                      &
     &      '# pairs    correlation',/,4x,65('-'))                      
 1100 FORMAT (5x,i3,3x,3(f6.2,' % ',2x),1x,i8,6x,f7.4) 
 1200 FORMAT (5x,i3,4x,'No pairs in this neighborhood')
!                                                                       
      END SUBROUTINE chem_corr_occ_mol              
!*****7*****************************************************************
      SUBROUTINE chem_neighbour (jatom, ic, iatom, patom, natom, maxw) 
!+                                                                      
!     Determine neighbours from given atom index 'jatom'.               
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE check_bound_mod
      USE chem_mod 
      USE celltoindex_mod
      USE do_find_mod
      USE metric_mod
      USE modify_mod
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      INTEGER iatom (0:maxw), jatom, natom, ic 
      REAL(kind=PREC_DP) :: patom (3, 0:maxw) 
!                                                                       
REAL(KIND=PREC_DP):: u (3), v (3), w (3), uu (3) 
REAL(kind=PREC_DP), dimension(3) :: offset !(3)
REAL(KIND=PREC_DP) :: dummy(1) 
      INTEGER jcell (3), icell (3), isite, jsite 
      INTEGER i, j, k, ii, iv, katom 
      LOGICAL lok 
!                                                                       
      natom = 0 
      dummy = - 1 
                                                                        
!                                                                       
!------ ------------------------------------------------------          
!------ Mode distance                                                   
!------ ------------------------------------------------------          
!                                                                       
      IF (chem_ctyp (ic) .eq.CHEM_DIST) then 
         DO j = 1, 3 
         u (j) = cr_pos (j, jatom) 
         w (j) = 0.0 
         ENDDO 
         CALL do_find_env (1, dummy, 1, u, chem_rmin (ic), chem_rmax (  &
         ic), chem_quick, chem_period)                                  
         IF (atom_env (0) .ne.0) then 
            DO j = 1, atom_env (0) 
            IF (chem_cang (ic) ) then 
!                                                                       
!------ ------- Check neighbouring angle and symmetry                   
!                                                                       
               DO k = 1, chem_nnei (ic) 
               DO ii = 1, 3 
               v (ii) = atom_pos (ii, j) - cr_pos (ii, jatom) 
               uu (ii) = chem_neig (ii, k, ic) 
               ENDDO 
               IF (abs (do_bang (.true., v, w, uu) )                    &
               .lt.chem_wink_sigma (ic) ) then                          
                  natom = natom + 1 
                  IF (natom.le.maxw) then 
                     iatom (natom) = atom_env (j) 
                     patom (1, natom) = atom_pos (1, j) 
                     patom (2, natom) = atom_pos (2, j) 
                     patom (3, natom) = atom_pos (3, j) 
                  ELSE 
                     ier_num = - 23 
                     ier_typ = ER_CHEM 
                     RETURN 
                  ENDIF 
               ENDIF 
               ENDDO 
!                                                                       
!------ ------- Check just distance                                     
!                                                                       
            ELSE 
               natom = natom + 1 
               IF (natom.le.maxw) then 
                  iatom (natom) = atom_env (j) 
                  patom (1, natom) = atom_pos (1, j) 
                  patom (2, natom) = atom_pos (2, j) 
                  patom (3, natom) = atom_pos (3, j) 
               ELSE 
                  ier_num = - 23 
                  ier_typ = ER_CHEM 
                  RETURN 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode vector                                                     
!------ ------------------------------------------------------          
!                                                                       
      ELSEIF (chem_ctyp (ic) .eq.CHEM_VEC) then 
         CALL indextocell (jatom, jcell, jsite) 
         DO i = 1, chem_nvec (ic) 
         iv = chem_use_vec (i, ic) 
         IF (jsite.eq.chem_cvec (1, iv) ) then 
            icell (1) = jcell (1) + chem_cvec (3, iv) 
            icell (2) = jcell (2) + chem_cvec (4, iv) 
            icell (3) = jcell (3) + chem_cvec (5, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cvec (2, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
         ELSEIF (jsite.eq.chem_cvec (2, iv) ) then 
            icell (1) = jcell (1) - chem_cvec (3, iv) 
            icell (2) = jcell (2) - chem_cvec (4, iv) 
            icell (3) = jcell (3) - chem_cvec (5, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cvec (1, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode angular neighbours                                         
!------ ------------------------------------------------------          
!                                                                       
      ELSEIF (chem_ctyp (ic) .eq.CHEM_ANG) then 
         CALL indextocell (jatom, jcell, jsite) 
         DO i = 1, chem_nwin (ic) 
         iv = chem_use_win (i, ic) 
         IF (jsite.eq.chem_cwin (1, iv) ) then 
            natom = natom + 1 
            iatom (natom) = jatom 
            patom (1, natom) = cr_pos (1, jatom) 
            patom (2, natom) = cr_pos (2, jatom) 
            patom (3, natom) = cr_pos (3, jatom) 
            icell (1) = jcell (1) + chem_cwin (3, iv) 
            icell (2) = jcell (2) + chem_cwin (4, iv) 
            icell (3) = jcell (3) + chem_cwin (5, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (2, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
            icell (1) = jcell (1) + chem_cwin (7, iv) 
            icell (2) = jcell (2) + chem_cwin (8, iv) 
            icell (3) = jcell (3) + chem_cwin (9, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (6, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
         ELSEIF (jsite.eq.chem_cwin (2, iv) ) then 
            icell (1) = jcell (1) - chem_cwin (3, iv) 
            icell (2) = jcell (2) - chem_cwin (4, iv) 
            icell (3) = jcell (3) - chem_cwin (5, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (1, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
            natom = natom + 1 
            iatom (natom) = jatom 
            patom (1, natom) = cr_pos (1, jatom) 
            patom (2, natom) = cr_pos (2, jatom) 
            patom (3, natom) = cr_pos (3, jatom) 
            icell (1) = jcell (1) - chem_cwin (3, iv) + chem_cwin (7,   &
            iv)                                                         
            icell (2) = jcell (2) - chem_cwin (4, iv) + chem_cwin (8,   &
            iv)                                                         
            icell (3) = jcell (3) - chem_cwin (5, iv) + chem_cwin (9,   &
            iv)                                                         
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (6, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
         ELSEIF (jsite.eq.chem_cwin (6, iv) ) then 
            icell (1) = jcell (1) - chem_cwin (7, iv) 
            icell (2) = jcell (2) - chem_cwin (8, iv) 
            icell (3) = jcell (3) - chem_cwin (9, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (1, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
            icell (1) = jcell (1) - chem_cwin (7, iv) + chem_cwin (3,   &
            iv)                                                         
            icell (2) = jcell (2) - chem_cwin (8, iv) + chem_cwin (4,   &
            iv)                                                         
            icell (3) = jcell (3) - chem_cwin (9, iv) + chem_cwin (5,   &
            iv)                                                         
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (2, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom = natom + 1 
               iatom (natom) = katom 
               patom (1, natom) = cr_pos (1, katom) + offset (1) 
               patom (2, natom) = cr_pos (2, katom) + offset (2) 
               patom (3, natom) = cr_pos (3, katom) + offset (3) 
            ENDIF 
            natom = natom + 1 
            iatom (natom) = jatom 
            patom (1, natom) = cr_pos (1, jatom) 
            patom (2, natom) = cr_pos (2, jatom) 
            patom (3, natom) = cr_pos (3, jatom) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_neighbour                 
!*****7*****************************************************************
      SUBROUTINE chem_neighbour_mol (jmol, ic, imol, nmol, maxww) 
!+                                                                      
!     Determine neighbours from given molecule index 'jmol'.            
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE chem_mod 
      USE molecule_mod 
      USE rmc_mod 
      USE errlist_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxww!, maxw 
!     PARAMETER (maxw = chem_max_neig) 
!                                                                       
      REAL(kind=PREC_DP) :: patom (3, 0:MAX_ATOM_ENV) 
      INTEGER imol (0:maxww), jmol, nmol, ic 
      INTEGER iatom (0:MAX_ATOM_ENV), jatom, natom 
      INTEGER i, im 
      LOGICAL no 
!                                                                       
      nmol = 0 
!                                                                       
!------ Get origin of selected molecule and find neighbours             
!                                                                       
      jatom = mole_cont (mole_off (jmol) + 1) 
      CALL chem_neighbour (jatom, ic, iatom, patom, natom, MAX_ATOM_ENV) 
!                                                                       
!------ Convert found atoms to molecules                                
!                                                                       
      IF (natom.gt.0) then 
         DO i = 1, natom 
         no = .true. 
         im = 0 
         DO while (no.and.im.lt.mole_num_mole) 
         im = im + 1 
         no = (iatom (i) .ne.mole_cont (mole_off (im) + 1) ) 
         ENDDO 
!                                                                       
         IF (.not.no) then 
            nmol = nmol + 1 
            IF (nmol.le.MAX_ATOM_ENV) then 
               imol (nmol) = im 
            ELSE 
               ier_num = - 23 
               ier_typ = ER_CHEM 
               RETURN 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE chem_neighbour_mol             
!*****7*****************************************************************
!
SUBROUTINE chem_blen (iianz, jjanz, werte, wwerte, maxw) 
!+                                                                      
!     Calculate bond length distribution within crystal                 
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE atom_name 
USE atom_env_mod 
USE chem_mod 
use discus_output_save_mod
USE do_find_mod
USE modify_mod
USE modify_func_mod
USE errlist_mod 
USE lib_f90_allocate_mod
USE param_mod 
USE precision_mod
USE prompt_mod 
!                                                                       
implicit none 
!                                                                       
integer, intent(in) :: iianz
integer, intent(in) :: jjanz
integer, intent(in) :: MAXW
real(kind=PREC_DP), dimension(maxw), intent(in) :: werte
real(kind=PREC_DP), dimension(maxw), intent(in) :: wwerte
!                                                                       
INTEGER i, j, k, l, is, js, ibin 
INTEGER bl_anz (0:maxscat, 0:maxscat), btot, btottot 
INTEGER :: n_res
REAL(kind=PREC_DP) :: u (3) 
REAL(KIND=PREC_DP) :: bl_min (0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP) :: bl_max (0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP) :: bl_sum (0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP) :: bl_s2 (0:maxscat, 0:maxscat) 
REAL(KIND=PREC_DP) :: bl_ave, bl_sig 
logical, dimension(0:maxscat, 0:maxscat) :: lprint   ! Print this pair
!                                                                       
CHARACTER(9) at_name_i 
CHARACTER(9) at_name_j 
!
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt
INTEGER                         :: all_status
!     LOGICAL atom_allowed 
!                                                                       
!------ write output                                                    
!                                                                       
WRITE (output_io, 500) chem_blen_cut, chem_fname, chem_bin 
!                                                                       
!------ Reset arrays                                                    
!                                                                       
bl_sum = 0.0
bl_s2  = 0.0
bl_anz = 0
bl_min = 0.0
bl_max = 0.0
chem_hist = 0
lprint = .TRUE.     ! Print all pairs
!      DO i = 0, cr_nscat 
!      DO j = 0, cr_nscat 
!      bl_sum (i, j) = 0.0 
!      bl_s2 (i, j) = 0.0 
!      bl_anz (i, j) = 0 
!      bl_min (i, j) = 999.0 
!      bl_max (i, j) = 0.0 
!      ENDDO 
!      ENDDO 
!      DO i = 1, CHEM_MAX_BIN 
!      chem_hist (i) = 0 
!      ENDDO 
!                                                                       
!------ Start the calculation                                           
!                                                                       
loop_atoms: DO i = 1, cr_natoms 
   IF (atom_allowed (i, werte, iianz, maxw) ) then 
!     DO j = 1, 3 
!        u (j) = cr_pos (j, i) 
!     ENDDO 
      u = cr_pos(1:3,i)
      CALL do_find_env(jjanz, wwerte, maxw, u, chem_blen_cut (1),               &
                       chem_blen_cut (2), chem_quick, chem_period)
      IF(ier_num.ne.0) return 
      DO k = 1, atom_env(0) 
         is = cr_iscat(i) 
         js = cr_iscat(atom_env (k) ) 
         bl_sum(is, js) = bl_sum(is, js) + res_para(k) 
         bl_s2 (is, js) = bl_s2 (is, js) + res_para(k) **2 
         bl_anz(is, js) = bl_anz(is, js) + 1 
         bl_min(is, js) = min(bl_min(is, js), res_para(k)) 
         bl_max(is, js) = max(bl_max(is, js), res_para(k)) 
                                                                        
         ibin = int((res_para(k) - chem_blen_cut(1)) * chem_bin /               &
                    (chem_blen_cut(2) - chem_blen_cut(1)))        + 1                 
         chem_hist (ibin) = chem_hist(ibin) + 1 
      ENDDO 
   ENDIF 
ENDDO  loop_atoms
!                                                                       
!------ write output                                                    
!                                                                       
l = 0 
btottot = 0 
n_res = (cr_nscat+2)*(cr_nscat+1)/2
IF(n_res                           > MAXPAR_RES) THEN
   n_res = MAX(n_res     ,MAXPAR_RES, MAX_ATOM_ENV)
   CALL alloc_param(n_res)
   MAXPAR_RES = n_res
ENDIF
!                                                                       
loop_fst: DO i = 0, cr_nscat  
   IF(-1==NINT(werte(1))  .OR. i==NINT(werte(1))  ) THEN
!     DO j = i, cr_nscat 
      loop_scd: DO j = 0, cr_nscat 
         if(lprint(i,j)) then               ! this pair still needs to be printed
         IF(-1==NINT(wwerte(1)) .OR. j==NINT(wwerte(1)) ) THEN
            IF (bl_anz (i, j) .ne.0.or.bl_anz (j, i) .ne.0) then 
               bl_ave = (bl_sum(i, j) + bl_sum(j, i)) / (bl_anz(i, j) + bl_anz(j, i))
               bl_sig = (bl_s2 (i, j) + bl_s2 (j, i)) / (bl_anz(i, j) + bl_anz(j, i))
               bl_sig = bl_sig - (bl_ave**2) 
               IF (bl_sig.gt.0.0) then 
                  bl_sig = sqrt (bl_sig) 
               ELSE 
                  bl_sig = 0.0 
               ENDIF 
            ELSE 
               bl_ave = 0.0 
               bl_sig = 0.0 
            ENDIF 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
!                                                                       
            IF (i.ne.j) then 
               btot = bl_anz (i, j) + bl_anz (j, i) 
            ELSE 
               btot = bl_anz (i, i) 
            ENDIF 
            btottot = btottot + btot 
!                                                                       
            WRITE (output_io, 1000) at_name_i, at_name_j, bl_ave, bl_sig,          &
            min(bl_min(i, j), bl_min(j, i)), max(bl_max(i, j), bl_max(j, i)), btot
!                                                                       
            IF ( (l + 3) <=  MAXPAR_RES) then 
               res_para (l + 1) = bl_ave 
               res_para (l + 2) = bl_sig 
               res_para (l + 3) = REAL(btot) 
               l = l + 3 
            ELSE 
               ier_num = - 2 
               ier_typ = ER_CHEM 
            ENDIF 
         ENDIF 
            lprint(j,i) = .FALSE.        ! Flag transposed pair as done
         ENDIF 
      ENDDO  loop_scd
   ENDIF 
ENDDO loop_fst
!                                                                       
      res_para (0) = l 
      WRITE (output_io, 1100) btottot 
!                                                                       
!------ write histogramm Allow optional write into kuplot
!                                                                       
      ALLOCATE(xwrt(1:chem_bin), stat=all_status)
      ALLOCATE(ywrt(1:chem_bin), stat=all_status)
      DO i = 1, chem_bin
         xwrt(i) = chem_blen_cut(1) + (chem_blen_cut(2) -chem_blen_cut(1)) * (i-1)/chem_bin
      ENDDO
      ywrt(1:chem_bin) = chem_hist(1:chem_bin)
      CALL output_save_file_1d( chem_fname, chem_bin, xwrt, ywrt, 0 )
      DEALLOCATE(xwrt, stat=all_status)
      DEALLOCATE(ywrt, stat=all_status)
!     OPEN (unit = 43, file = chem_fname, status = 'unknown',iostat=ios, &
!           IOMSG=message) 
!     IF(ios/=0) THEN
!        ier_num = -2
!        ier_typ = ER_IO
!        ier_msg(1)(1:60) = chem_fname(1:60)
!        ier_msg(3) = message(1:80)
!        RETURN
!     ENDIF
!     DO i = 1, chem_bin 
!     WRITE (43, 5000) chem_blen_cut (1) + (chem_blen_cut (2) -         &
!     chem_blen_cut (1) ) * (i - 1) / chem_bin, chem_hist (i)           
!     ENDDO 
!     CLOSE (43) 
!                                                                       
  500 FORMAT (' Calculating bond-length distibution',/,                 &
     &        '    Allowed range : ',F6.2,' A to ',F6.2,                &
     &        '  A / File : ',A12,' (',I4,' pts)',/)                    
 1000 FORMAT ('    ',A9,'- ',A9,': d =',F7.3,' +- ',F6.3,' A ',         &
     &        '(Min =',F7.3,', Max =',F7.3,')',/,                       &
     &        49x,'(Pairs = ',i18,')')                                  
 1100 FORMAT (49x,'(Total = ',i18,')') 
!5000 FORMAT (F8.3,I12) 
      END SUBROUTINE chem_blen                      
!*****7*****************************************************************
      SUBROUTINE chem_blen_cluster (iianz, jjanz, werte, wwerte, maxw) 
!+                                                                      
!     Calculate bond length distribution within crystal                 
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE chem_mod 
use discus_output_save_mod
      USE metric_mod
      USE modify_func_mod
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw, iianz, jjanz 
      REAL(KIND=PREC_DP) :: werte (maxw), wwerte (maxw) 
!                                                                       
      INTEGER i, j, k, ibin
      REAL(kind=PREC_DP) :: u (3), v (3), dist 
!
      REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt
      REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt
      INTEGER                         :: all_status
!                                                                       
!                                                                       
!------ write output                                                    
!                                                                       
      WRITE (output_io, 400) 
      WRITE (output_io, 500) chem_blen_cut, chem_fname, chem_bin 
!                                                                       
!------ Reset arrays                                                    
!                                                                       
      DO i = 1, CHEM_MAX_BIN 
      chem_hist (i) = 0 
      ENDDO 
!                                                                       
!------ Start the calculation                                           
!                                                                       
      DO i = 1, cr_natoms - 1 
      IF (atom_allowed (i, werte, iianz, maxw) ) then 
         DO k = 1, 3 
         u (k) = cr_pos (k, i) 
         ENDDO 
         DO j = i + 1, cr_natoms 
         IF (atom_allowed (j, wwerte, jjanz, maxw) ) then 
            DO k = 1, 3 
            v (k) = cr_pos (k, j) 
            ENDDO 
            dist = do_blen (.true., u, v) 
!                                                                       
            ibin = int ( (dist - chem_blen_cut (1) ) * chem_bin /       &
            (chem_blen_cut (2) - chem_blen_cut (1) ) ) + 1              
            IF (0.lt.ibin.and.ibin.le.CHEM_MAX_BIN) then 
               chem_hist (ibin) = chem_hist (ibin) + 1 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
!------ write histogramm                                                
!                                                                       
      ALLOCATE(xwrt(1:chem_bin), stat=all_status)
      ALLOCATE(ywrt(1:chem_bin), stat=all_status)
      DO i = 1, chem_bin
         xwrt(i) = chem_blen_cut(1) + (chem_blen_cut(2) -chem_blen_cut(1)) * (i-1)/chem_bin
      ENDDO
      ywrt(1:chem_bin) = chem_hist(1:chem_bin)
      CALL output_save_file_1d( chem_fname, chem_bin, xwrt, ywrt, 0 )
      DEALLOCATE(xwrt, stat=all_status)
      DEALLOCATE(ywrt, stat=all_status)
!     OPEN (unit = 43, file = chem_fname, status = 'unknown', &
!           IOSTAT=ios, IOMSG=message) 
!     IF(ios/=0) THEN
!        ier_num = -2
!        ier_typ = ER_IO
!        ier_msg(1)(1:60) = chem_fname(1:60)
!        ier_msg(3) = message(1:80)
!        RETURN
!     ENDIF
!     DO i = 1, chem_bin 
!     WRITE (43, 5000) chem_blen_cut (1) + (chem_blen_cut (2) -         &
!     chem_blen_cut (1) ) * (i - 1) / chem_bin, chem_hist (i)           
!     ENDDO 
!     CLOSE (43) 
!                                                                       
  400 FORMAT (/,' Calculating bond-length distibution (Mode: CLUSTER)') 
  500 FORMAT (' Calculating bond-length distibution',/,                 &
     &        '    Allowed range : ',F6.2,' A to ',F6.2,                &
     &        '  A / File : ',A12,' (',I4,' pts)',/)                    
!5000 FORMAT (F8.3,I12) 
      END SUBROUTINE chem_blen_cluster              
!*****7*****************************************************************
      SUBROUTINE chem_bang (iianz, jjanz, kkanz, werte, wwerte, uwerte, &
      maxw)                                                             
!+                                                                      
!     Calculate bond angle distribution within crystal                  
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
use discus_output_save_mod
      USE atom_name 
      USE atom_env_mod 
      USE chem_mod 
      USE do_find_mod
      USE metric_mod
      USE modify_mod
      USE modify_func_mod
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw, iianz, jjanz, kkanz 
      REAL(KIND=PREC_DP) :: werte (maxw), wwerte (maxw), uwerte (maxw) 
!                                                                       
      INTEGER i, j, k, l, is, js, ibin 
      INTEGER ba_anz (0:maxscat, 0:maxscat) 
      INTEGER ba_env (0:MAX_ATOM_ENV) 
      LOGICAL lspace 
!
      REAL(kind=PREC_DP) , DIMENSION(:), ALLOCATABLE :: xwrt
      REAL(kind=PREC_DP) , DIMENSION(:), ALLOCATABLE :: ywrt
      INTEGER                         :: all_status
!                                                                       
!     von der relativen Reihenfolge der beiden Statements haengt es ab, 
!     ob der zweite Nachbar gefunden wird oder nicht !?!?               
!                                                                       
      REAL(kind=PREC_DP) :: u (3), v (3), w (3), angle 
      REAL(kind=PREC_DP) :: ba_pos (3, MAX_ATOM_ENV) 
                                                                        
      REAL(kind=PREC_DP) :: ba_min (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_max (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_sum (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_s2 (0:maxscat, 0:maxscat) 
      REAL(kind=PREC_DP) :: ba_ave, ba_sig 
!                                                                       
      CHARACTER(9) at_name_i 
      CHARACTER(9) at_name_j 
!                                                                       
      lspace = .true. 
!                                                                       
!------ write output                                                    
!                                                                       
      WRITE (output_io, 500) chem_bang_cut, chem_blen_cut, chem_fname,  &
      chem_bin                                                          
!                                                                       
!------ Reset arrays                                                    
!                                                                       
      DO i = 0, cr_nscat 
      DO j = 0, cr_nscat 
      ba_sum (i, j) = 0.0 
      ba_s2 (i, j) = 0.0 
      ba_anz (i, j) = 0 
      ba_min (i, j) = 1000.0 
      ba_max (i, j) = - 1000.0 
      ENDDO 
      ENDDO 
      DO i = 1, CHEM_MAX_BIN 
      chem_hist (i) = 0 
      ENDDO 
!                                                                       
!------ Start the calculation                                           
!                                                                       
      DO i = 1, cr_natoms 
      IF (atom_allowed (i, werte, iianz, maxw) ) then 
         DO j = 1, 3 
         u (j) = cr_pos (j, i) 
         ENDDO 
         CALL do_find_env (jjanz, wwerte, maxw, u, chem_blen_cut (1),   &
         chem_blen_cut (2), chem_quick, chem_period)                    
         IF (ier_num.ne.0) return 
!                                                                       
!     ----copy environment of first atom type                           
!                                                                       
         ba_env (0) = atom_env (0) 
         DO k = 1, atom_env (0) 
         ba_env (k) = atom_env (k) 
         DO j = 1, 3 
         ba_pos (j, k) = atom_pos (j, k) 
         ENDDO 
         ENDDO 
         CALL do_find_env (kkanz, uwerte, maxw, u, chem_blen_cut (1),   &
         chem_blen_cut (2), chem_quick, chem_period)                    
         IF (ier_num.ne.0) return 
!                                                                       
!     ----For all atoms in the environment, except atom i itself        
!                                                                       
         DO k = 1, ba_env (0) 
         IF (ba_env (k) .ne.i) then 
            is = cr_iscat (ba_env (k) ) 
            DO j = 1, 3 
            v (j) = ba_pos (j, k) 
            ENDDO 
!                                                                       
!     --------For all second neighbor atoms in the environment,         
!             except atom i and first neighbor atom                     
!                                                                       
            DO l = 1, atom_env (0) 
            IF (atom_env (l) .ne.i.and.atom_env (l) .ne.ba_env (k) )    &
            then                                                        
               js = cr_iscat (atom_env (l) ) 
               DO j = 1, 3 
               w (j) = atom_pos (j, l) 
               ENDDO 
               angle = do_bang (lspace, v, u, w) 
               ba_sum (is, js) = ba_sum (is, js) + angle 
               ba_s2 (is, js) = ba_s2 (is, js) + angle**2 
               ba_anz (is, js) = ba_anz (is, js) + 1 
               ba_min (is, js) = min (ba_min (is, js), angle) 
               ba_max (is, js) = max (ba_max (is, js), angle) 
                                                                        
               ibin = int ( (angle-chem_bang_cut (1) ) * chem_bin /     &
               (chem_bang_cut (2) - chem_bang_cut (1) ) ) + 1           
               IF (0.lt.ibin.and.ibin.le.CHEM_MAX_BIN) then 
                  chem_hist (ibin) = chem_hist (ibin) + 1 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
!------ write output                                                    
!                                                                       
      DO i = 0, cr_nscat 
      DO j = i, cr_nscat 
      IF (ba_anz (i, j) .ne.0.or.ba_anz (j, i) .ne.0) then 
         ba_ave = (ba_sum (i, j) + ba_sum (j, i) ) / (ba_anz (i, j)     &
         + ba_anz (j, i) )                                              
         ba_sig = (ba_s2 (i, j) + ba_s2 (j, i) ) / (ba_anz (i, j)       &
         + ba_anz (j, i) )                                              
         ba_sig = ba_sig - (ba_ave**2) 
         IF (ba_sig.gt.0.0) then 
            ba_sig = sqrt (ba_sig) 
         ELSE 
            ba_sig = 0.0 
         ENDIF 
         at_name_i = at_name (i) 
         at_name_j = at_name (j) 
         WRITE (output_io, 1000) at_name_i, at_name_j, ba_ave, ba_sig,  &
         min (ba_min (i, j), ba_min (j, i) ), max (ba_max (i, j),       &
         ba_max (j, i) )                                                
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
!------ write histogramm                                                
!                                                                       
      ALLOCATE(xwrt(1:chem_bin), stat=all_status)
      ALLOCATE(ywrt(1:chem_bin), stat=all_status)
      DO i = 1, chem_bin
         xwrt(i) = chem_bang_cut(1) + (chem_bang_cut(2) -chem_bang_cut(1)) * (i-1)/chem_bin
      ENDDO
      ywrt(1:chem_bin) = real(chem_hist(1:chem_bin), kind=PREC_DP)
      CALL output_save_file_1d( chem_fname, chem_bin, xwrt, ywrt, 0 )
      DEALLOCATE(xwrt, stat=all_status)
      DEALLOCATE(ywrt, stat=all_status)
!     OPEN (unit = 43, file = chem_fname, status = 'unknown', &
!           IOSTAT=ios, IOMSG=message) 
!     IF(ios/=0) THEN
!        ier_num = -2
!        ier_typ = ER_IO
!        ier_msg(1)(1:60) = chem_fname(1:60)
!        ier_msg(3) = message(1:80)
!        RETURN
!     ENDIF
!     DO i = 1, chem_bin 
!     WRITE (43, 5000) chem_bang_cut (1) + (chem_bang_cut (2) -         &
!     chem_bang_cut (1) ) * (i - 1) / chem_bin, chem_hist (i)           
!     ENDDO 
!     CLOSE (43) 
!                                                                       
  500 FORMAT     (' Calculating bond-angle distibution',/,              &
     &        '    Allowed range : ',F6.2,'   to ',F6.2,' Degrees',/,   &
     &        '    Allowed length: ',F6.2,' A to ',F6.2,                &
     &        '  A / File : ',A12,' (',I4,' pts)',/)                    
 1000 FORMAT     ('    ',A9,'- ',A9,': a = ',F7.3,' +- ',F7.3,' Deg ',  &
     &                   '(Min = ',F7.3,', Max = ',F7.3,')')            
!5000 FORMAT     (F8.3,I12) 
      END SUBROUTINE chem_bang                      
!*****7*****************************************************************
      SUBROUTINE chem_mole (lout) 
!+                                                                      
!     Show information about molecules/rel. amounts within crystal      
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
USE atom_env_mod
      USE chem_mod 
      USE molecule_mod 
      USE errlist_mod 
      USE lib_f90_allocate_mod
      USE param_mod 
      USE prompt_mod 
use precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP) :: proz 
      INTEGER nmole (0:MOLE_MAX_TYPE) 
      INTEGER i 
      INTEGER :: n_res 
      LOGICAL lout 
!                                                                       
!     Error condition                                                   
!                                                                       
      IF (mole_num_mole.eq.0) then 
         ier_typ = ER_CHEM 
         ier_num = - 20 
         RETURN 
      ENDIF 
!                                                                       
!------ reset counters, ...                                             
!                                                                       
      DO i = 1, mole_num_type 
      nmole (i) = 0 
      ENDDO 
!                                                                       
!------ get relative amounts of molecule types                          
!                                                                       
      DO i = 1, mole_num_mole 
      nmole (mole_type (i) ) = nmole (mole_type (i) ) + 1 
      ENDDO 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) then 
         WRITE (output_io, 1000) (cr_icc (i), i = 1, 3) 
         WRITE (output_io, 1100) mole_num_mole, mole_num_unit,          &
         mole_num_type                                                  
      ENDIF 
!                                                                       
      IF(mole_num_type                   > MAXPAR_RES) THEN
         n_res = MAX(mole_num_type, MAXPAR_RES, MAX_ATOM_ENV)
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
      res_para (0) = REAL(mole_num_type) 
      DO i = 1, mole_num_type 
      proz = REAL(nmole (i) ) / mole_num_mole 
      IF (lout) then 
         WRITE (output_io, 1200) i, proz, nmole (i) 
      ENDIF 
      IF (i <= MAXPAR_RES) then 
         res_para (i + 1) = proz 
      ELSE 
         ier_typ = ER_CHEM 
         ier_num = - 2 
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT (' Size of the crystal (unit cells)  : ',2(I4,' x '),I4) 
 1100 FORMAT (' Total number of molecules         : ',I9,/              &
     &                   ' Number of molecules per unit cell : ',I9,/   &
     &                   ' Number of different molecules     : ',I9,/)  
 1200 FORMAT     ('    Mol. type : ',I4,5X,' rel. abundance : ',F5.3,   &
     &                   '  (',I9,' molecules)')                        
      END SUBROUTINE chem_mole                      
!
!*******************************************************************************
!
SUBROUTINE chem_trans(zeile,lp)
!
USE crystal_mod
USE chem_mod
USE celltoindex_mod
USE errlist_mod
USE get_params_mod
USE ber_params_mod
USE param_mod
USE precision_mod
USE prompt_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER                  :: MAXW=6
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))),DIMENSION(MAXW) :: cpara
INTEGER            ,DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) ,DIMENSION(MAXW) :: werte
INTEGER                             :: ianz    ! loop index
INTEGER               :: i    ! loop index
INTEGER               :: ia   ! Atom index
INTEGER, DIMENSION(3) :: ic   ! Cell number
INTEGER               :: is   ! site number
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: O_BOUND   = 1
INTEGER, PARAMETER :: O_OUT     = 2
CHARACTER(LEN=   8), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'boundary', 'output'   /
DATA loname /  8,          6         /
opara  =  (/ 'crystal', 'screen '     /)   ! Always provide fresh default values
lopara =  (/  7       ,  6            /)
owerte =  (/  0.0     ,  0.0          /)
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF(ier_NUM/=0) RETURN
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_NUM/=0) RETURN
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num == 0 .AND. cr_natoms /= 0) THEN 
   IF (ianz == 1) THEN 
      ia = NINT(werte (1) ) 
      IF (ia > 0 .AND. ia <= cr_natoms) THEN 
         CALL indextocell (ia, ic, is) 
         IF(opara(O_OUT)=='screen') THEN
            WRITE (output_io, 3000) ia, ic, is 
         ENDIF
         res_para(0  ) = 5
         res_para(1  ) = ia
         res_para(2:4) = ic
         res_para(5  ) = is
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
  ELSEIF (ianz == 4) THEN 
      is = NINT(werte(4) ) 
      DO i = 1, 3 
         ic (i) = NINT(werte(i)) 
      ENDDO 
      IF(opara(O_BOUND)=='periodic') THEN
         DO i = 1, 3 
            IF(chem_period(i)) THEN
               ic(i) = MOD(ic(i) + cr_icc(i) -1, cr_icc(i)) + 1
            ENDIF 
         ENDDO 
      ENDIF
      IF(ic(1) > 0 .AND. ic(1) <= cr_icc(1) .AND.  &
         ic(2) > 0 .AND. ic(2) <= cr_icc(2) .AND.  &
         ic(3) > 0 .AND. ic(3) <= cr_icc(3) .AND.  &
         is >  0   .AND. is  <= cr_ncatoms         ) THEN               
         CALL celltoindex (ic, is, ia) 
         IF(opara(O_OUT)=='screen') THEN
            WRITE (output_io, 3000) ia, ic, is 
         ENDIF
         res_para(0  ) = 5
         res_para(1  ) = ia
         res_para(2:4) = ic
         res_para(5  ) = is
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
!                                                                       
 3000 FORMAT    (' Atomindex ',I6,' : Unitcell ',3(I4,1X),              &
     &                  ' / site ',I2)                                  
END SUBROUTINE chem_trans
!
!*******************************************************************************
!
SUBROUTINE chem_apply_period(isel,lupdate_conn)
!
!  Applies periodic boundary conditions to the atom position of 
!  atom number isel.
!  If the connectivity mode is choosen, then the offsets to
!  and from the neighbors are updated.
!
USE crystal_mod
USE chem_mod
USE molecule_mod
USE conn_mod
!
use precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: isel          ! selected atom number
LOGICAL, INTENT(IN) :: lupdate_conn  ! update the connectivity, if true
!
REAL(kind=PREC_DP)   , PARAMETER  :: DELTA = 0.1
!
INTEGER               :: i, j, k, imole
LOGICAL               :: lshift
REAL(kind=PREC_DP)   , DIMENSION(3) :: shift
!
lshift   = .FALSE.
DO i=1,3
   shift(i) = 0.0
   IF(chem_period(i)) THEN
   IF(    cr_pos(i,isel) < NINT(cr_dim(i,1))          -DELTA) THEN
      shift(i) =  cr_icc(i)
      lshift   = .TRUE.
   ELSEIF(cr_pos(i,isel) > NINT(cr_dim(i,1))+cr_icc(i)+DELTA) THEN
      shift(i) = -cr_icc(i)
      lshift   = .TRUE.
   ENDIF
   ENDIF
ENDDO
!
IF(lshift) THEN   ! Periodic boundary conditions were applied
   IF(cr_mole(isel)==0) THEN    ! Single atom
      cr_pos(:,isel) = cr_pos(:,isel) + shift(:)
      IF(lupdate_conn) THEN
         CALL conn_update(isel, shift)
      ENDIF
   ELSE
      imole = cr_mole(isel)
      DO j=1, mole_len(imole)
         k = mole_cont(mole_off(imole)+j)
         cr_pos(:,k) = cr_pos(:,k) + shift(:)
         IF(lupdate_conn) THEN
            CALL conn_update(k, shift)
         ENDIF
      ENDDO
   ENDIF
!  IF(lupdate_conn) THEN
!     CALL conn_update(isel, shift)
!  ENDIF
ENDIF
!
END SUBROUTINE chem_apply_period
!
!*****7***************************************************************  
!
SUBROUTINE chem_reset

USE discus_allocate_appl_mod
USE chem_mod
!
IMPLICIT NONE
!
CALL alloc_chem_correlation( 1     )
CALL alloc_chem_ang ( 1,  1        )
CALL alloc_chem_aver( 1,  1        )
CALL alloc_chem_disp( 1,  1        )
CALL alloc_chem_env ( 1,  1,  1    )
CALL alloc_chem_vec ( 1,  1        )
CALL alloc_chem_ran ( 1,  1 ,  1   )
CALL alloc_chem_con ( 1,  1        )
CALL alloc_chem_dir (     1        )
CALL alloc_chem_dist(     1        )
CALL alloc_chem_hist(     1        )
!
CHEM_MAXAT_CELL   = 1
CHEM_MAX_AVE_ATOM = 1
CHEM_MAX_VEC      = 1
CHEM_MAX_ANG      = 1
CHEM_MAX_RAN      = 0
CHEM_MAX_CON      = 1
CHEM_MAX_ENV      = 1
!
chem_ndef      = 0
chem_fname     = 'blen.xy'
chem_bin       = 601
chem_ncor      =   0
chem_blen_cut(:)  = (/1.5,  7.5/)
chem_bang_cut(:)  = (/0.0,180.0/)
chem_quick     = .true.
chem_cluster   = .false.
chem_period(:) = .true.
chem_purge     = .FALSE.
chem_sel_atom  = .true.
!
! Vector definitions
IF(ALLOCATED(chem_cvec))    chem_cvec(:,:)    = 0        !  (5,CHEM_MAX_VEC)
IF(ALLOCATED(chem_cvec))    chem_cvec(1,:)    = -9999    !  (5,CHEM_MAX_VEC)
IF(ALLOCATED(chem_use_vec)) chem_use_vec(:,:) = 1     !  (CHEM_MAX_VEC,CHEM_MAX_COR)
!
! Range definitions
IF(ALLOCATED(chem_use_ran))     chem_use_ran(:,:)    = 1     !  (CHEM_MAX_RAN,CHEM_MAX_COR)
IF(ALLOCATED(chem_cran_cent))   chem_cran_cent(:,:)  = -9999   !  (0:CHEM_MAX_ATOM,CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_neig))   chem_cran_neig(:,:)  = 0   !  (0:CHEM_MAX_ATOM,CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_nuvw))   chem_cran_nuvw(:)    = 0   !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_nshort)) chem_cran_nshort(:)  = 0 !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_uvw))    chem_cran_uvw(:,:,:) = 0.0    !  (3,48,CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_sig))    chem_cran_sig(:)     = 0.0    !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_wsig))   chem_cran_wsig(:)    = 0.0   !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_rmax))   chem_cran_rmax(:)    = 0.0   !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_rmin))   chem_cran_rmin(:)    = 0.0   !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_cang))   chem_cran_cang(:)    = .FALSE. !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_lsym))   chem_cran_lsym(:)    = .FALSE. !  (CHEM_MAX_RAN)
IF(ALLOCATED(chem_cran_short))  chem_cran_short(:)   = .FALSE. !  (CHEM_MAX_RAN)
!
! Angle definitions
IF(ALLOCATED(chem_cwin))    chem_cwin(:,:)  = -9999        !  (9,CHEM_MAX_ANG)
IF(ALLOCATED(chem_use_win)) chem_use_win(:,:) = 0     !  (CHEM_MAX_ANG,CHEM_MAX_COR)
!
! Connectivity definitions
IF(ALLOCATED(chem_ccon))    chem_ccon(:,:) = -9999        !  (2,CHEM_MAX_ANG)
IF(ALLOCATED(chem_cname))   chem_cname(:) = ' '  !  MMC Correlation name
IF(ALLOCATED(chem_cname_l)) chem_cname_l(:) = 1     !  (2,CHEM_MAX_ANG)
IF(ALLOCATED(chem_use_con)) chem_use_con(:,:) = 0     !  (CHEM_MAX_CON,CHEM_MAX_COR)
!
! Environment definitions
!
chem_hist(:) = 0 !  (CHEM_MAX_BIN)
!
! Average crystal structure
IF(ALLOCATED(chem_ave_n))     chem_ave_n(:)       = 0       !  (MAXAT_CELL)
IF(ALLOCATED(chem_ave_iscat)) chem_ave_iscat(:,:) = 0   !  (MAXAT_CELL,CHEM_MAX_ATOM)
IF(ALLOCATED(chem_ave_pos))  chem_ave_pos(:,:)   = 0.0     !  (3,MAXAT_CELL)
IF(ALLOCATED(chem_ave_sig))  chem_ave_sig(:,:)   = 0.0     !  (3,MAXAT_CELL)
IF(ALLOCATED(chem_ave_bese)) chem_ave_bese(:,:)  = 0.0    !  (MAXAT_CELL,CHEM_MAX_ATOM)
chem_run_aver     = .TRUE.
chem_run_aver_ind = .TRUE.
!
! Displacement correlations
IF(ALLOCATED(chem_disp_ave)) chem_disp_ave(:,:,:) = 0.0    !  (CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
IF(ALLOCATED(chem_disp_sig)) chem_disp_sig(:,:,:) = 0.0    !  (CHEM_MAX_COR,0:DEF_MAXSCAT,0:DEF_MAXSCAT)
!
! Environment interactions
IF(ALLOCATED(chem_use_env))  chem_use_env(:,:) = -9999     !  (CHEM_MAX_ENV,CHEM_MAX_COR)
IF(ALLOCATED(chem_cenv))     chem_cenv(:,:) = -9999        !  (0:MAX_ATOM_ENV,CHEM_MAX_ENV)
IF(ALLOCATED(chem_env_neig)) chem_env_neig(:) = 0    !  (CHEM_MAX_ENV)
IF(ALLOCATED(chem_rmax_env)) chem_rmax_env(:) = 3.5  !  (CHEM_MAX_ENV)
IF(ALLOCATED(chem_rmin_env)) chem_rmin_env(:) = 0.1  !  (CHEM_MAX_ENV)
!
chem_nran(:) = 0 !  (CHEM_MAX_COR)
chem_nvec(:) = 0 !  (CHEM_MAX_COR)
chem_ncon(:) = 0 !  (CHEM_MAX_COR)
chem_nwin(:) = 0 !  (CHEM_MAX_COR)
chem_ctyp(:) = 0 !  (CHEM_MAX_COR)
chem_nnei(:) = 0 !  (CHEM_MAX_COR)
chem_nenv(:) = 0 !  (CHEM_MAX_COR)
chem_neig(:,:,:) = 0.0        !  (3,48,CHEM_MAX_COR)
chem_dir(:,:,:) = -9999.0     !  (3,2,CHEM_MAX_COR)
chem_rmax(:)     = 0.0        !  (CHEM_MAX_COR)
chem_rmin(:)     = 0.0        !  (CHEM_MAX_COR)
chem_freq_sigma(:)     = 0.0  !  (CHEM_MAX_COR)
chem_wink_sigma(:)     = 0.0  !  (CHEM_MAX_COR)
chem_cang(:)       = .FALSE.  !  (CHEM_MAX_COR)
chem_ldall(:)      = .TRUE.   !  (CHEM_MAX_COR)
!
!
END SUBROUTINE chem_reset
!
!*****7***************************************************************  
!
END MODULE chem_menu
