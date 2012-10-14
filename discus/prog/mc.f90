!*****7*****************************************************************
!                                                                       
      SUBROUTINE mc 
!-                                                                      
!     This sublevel includes all commands and functions for the         
!     Monte-Carlo simulations in DISCUS.                                
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE mc_mod 
      USE modify_mod
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'macro.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'param.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(40) cdummy 
      CHARACTER(1024) line, zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER lp, length 
      INTEGER indxg, lbef, ianz 
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!
      ier_num = +4
      ier_typ = ER_APPL
      ier_msg(1) = 'MC Monte-Carlo menu is no longer valid'
      ier_msg(2) = 'Please move to the mmc Monte-Carlo menu'
      ier_msg(3) = 'MMC has many new and improved features'
      RETURN
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL no_error 
      prom = prompt (1:len_str (prompt) ) //'/mc' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
            CALL do_math (line, indxg, length) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            CALL file_kdo (line (2:length), length - 1) 
!                                                                       
!------ continues a macro 'continue'                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
            CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
            CALL do_eval (zeile, lp) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 11 
               CALL do_hel ('discus mc '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
            cdummy = ' ' 
            IF (zeile.ne.' ') then 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     Waiting for user input                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
            CALL do_input (zeile, lp) 
!                                                                       
!------ ------------------------------------------------------------    
!------ Here start MC specific commands                                 
!------ ------------------------------------------------------------    
!                                                                       
!------ command 'run'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
            CALL mc_run 
!                                                                       
!------ command 'save'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.eq.0) then 
               WRITE (output_io, 1500) cpara (1) (1:len_str (cpara (1) )&
               )                                                        
               CALL save_struc (cpara (1), lpara (1) ) 
            ENDIF 
!                                                                       
!------ command 'set'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
            CALL mc_set (zeile, lp) 
!                                                                       
!------ command 'show'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            CALL mc_show 
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
            IF (lmakro) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
      GOTO 10 
!                                                                       
 9999 CONTINUE 
!                                                                       
 1500 FORMAT   (' Saving structure to file : ',a) 
!                                                                       
      END SUBROUTINE mc                             
!*****7*****************************************************************
      SUBROUTINE mc_show 
!+                                                                      
!     Show parameters of MC section                                     
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(1024) zeile 
      CHARACTER(9) at_name, at_name_i, at_name_j 
      INTEGER i, j, k 
!                                                                       
!------ Information about MC running mode etc ..                        
!                                                                       
      IF (mo_mode.eq.rmc_mode_shift) then 
         WRITE (output_io, 1100) 'move atoms' 
      ELSEIF (mo_mode.eq.rmc_mode_swchem) then 
         WRITE (output_io, 1100) 'switch atoms' 
      ELSEIF (mo_mode.eq.rmc_mode_swdisp) then 
         WRITE (output_io, 1100) 'switch displacements' 
      ELSEIF (mo_mode.eq.rmc_mode_external) then 
         WRITE (output_io, 1100) 'defined in user supplied subroutine' 
      ELSE 
         WRITE (output_io, 1100) '??????????' 
      ENDIF 
!                                                                       
      IF (mo_sel_atom) then 
         WRITE (output_io, 1105) 'atoms' 
      ELSE 
         WRITE (output_io, 1105) 'molecules' 
      ENDIF 
!                                                                       
      IF (mo_energy.eq.mc_none) then 
         WRITE (output_io, 1110) 'not defined' 
      ELSEIF (mo_energy.eq.mc_occ) then 
         WRITE (output_io, 1110) 'correlations occupational' 
      ELSEIF (mo_energy.eq.mc_disp) then 
         WRITE (output_io, 1110) 'correlations displacements' 
      ELSEIF (mo_energy.eq.mc_spring) then 
         WRITE (output_io, 1110) 'distortions (Hookes law)' 
      ELSEIF (mo_energy.eq.mc_angle) then 
         WRITE (output_io, 1110) 'angular distortions' 
      ELSE 
         WRITE (output_io, 1110) '??????????' 
      ENDIF 
!                                                                       
      IF (mo_local.eq.rmc_local_all) then 
         WRITE (output_io, 1120) 'all' 
      ELSEIF (mo_local.eq.rmc_local_loc) then 
         WRITE (output_io, 1120) 'only local (+-1 unit cell)' 
      ELSEIF (mo_local.eq.rmc_local_locsite) then 
         WRITE (output_io, 1120) 'only local, same site' 
      ELSEIF (mo_local.eq.rmc_local_site) then 
         WRITE (output_io, 1120) 'all, same site' 
      ELSE 
         WRITE (output_io, 1120) '??????????' 
      ENDIF 
!                                                                       
      WRITE (output_io, 1250) mo_cyc 
      WRITE (output_io, 1300) mo_feed 
      WRITE (output_io, 1400) mo_kt 
!                                                                       
!------ Information about target values, constants etc ..               
!                                                                       
      WRITE (output_io, 2000) mo_const (0) 
      DO i = 1, chem_ncor 
      WRITE (output_io, 2100) i, mo_const (i), mo_cfac (i), mo_ach_corr &
      (i), mo_target_corr (i)                                           
      ENDDO 
!                                                                       
!------ Information about defined correlations etc ..                   
!                                                                       
      WRITE (output_io, 3000) 
      zeile = 'corr' 
      CALL chem_show (zeile) 
!                                                                       
!------ Information about defined desired distances                     
!                                                                       
      IF (mo_mode.eq.MC_SPRING) then 
         WRITE (output_io, 4000) 
!                                                                       
         IF (mo_sel_atom) then 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            DO k = 1, chem_ncor 
            IF (mo_disp (k, i, j) .ne.0.0) write (output_io, 4100)      &
            at_name_i, at_name_j, k, mo_disp (k, i, j)                  
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            DO k = 1, chem_ncor 
            IF (mo_disp (k, i, j) .ne.0.0) write (output_io, 4200) i, j,&
            k, mo_disp (k, i, j)                                        
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Maxmove settings in 'shift' mode                                
!                                                                       
      IF (mo_mode.eq.rmc_mode_shift.and. (                              &
      mo_energy.eq.mc_disp.or.mo_energy.eq.mc_spring) ) then            
         WRITE (output_io, 4300) 
         IF (mo_sel_atom) then 
            DO i = 1, cr_nscat 
            at_name_i = at_name (i) 
            WRITE (output_io, 4310) at_name_i, (mo_maxmove (j, i),      &
            j = 1, 3)                                                   
            ENDDO 
         ELSE 
            DO i = 1, mole_num_type 
            WRITE (output_io, 4320) i, (mo_maxmove (j, i), j = 1, 3) 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
 1100 FORMAT     (  ' MC running mode                : ',a) 
 1105 FORMAT     (  '   Operation mode for MC        : ',a) 
 1110 FORMAT     (  '   Energy calculation mode      : ',a) 
 1120 FORMAT     (  '   Valid neighbours             : ',a) 
 1250 FORMAT     (  '   Max. number of MC cycles     : ',i8) 
 1300 FORMAT     (  '   Feedback/update intervall    : ',i8) 
 1400 FORMAT     (  '   Temperature [kT]             : ',f8.4) 
 2000 FORMAT     (/,' MC settings (C(0)=',f10.3,')  : ',//,             &
     &        '    #',4x,'interaction',2x,'factor',5x,'curr. corr.',    &
     &             2x,'target corr.',/,3x,55('-'))                      
 2100 FORMAT     (3x,i3,4x,f10.3,2x,f6.3,7x,f6.3,7x,f6.3) 
 3000 FORMAT     (/,' Correlation definitions        : ',/) 
 4000 FORMAT     (/,' Desired distortions for SPRING : ',/) 
 4100 FORMAT     (  '   Atoms : ',a9,' - ',a9,'  neig. #',i3,           &
     &                     '   distance : ',f7.3,' A')                  
 4200 FORMAT     (  '   Molecule types : ',i4,' - ',i4,'  neig. #',i3,  &
     &                     '   distance : ',f7.3,' A')                  
 4300 FORMAT     (/,' Sigmas for MC shifts (l.u.)    : ') 
 4310 FORMAT     (  '                 Atom ',a9,' : ',3(F7.3,1X)) 
 4320 FORMAT     (  '        Molecule type ',i9,' : ',3(F7.3,1X)) 
      END SUBROUTINE mc_show                        
!*****7*****************************************************************
      SUBROUTINE mc_set (zeile, lp) 
!+                                                                      
!     sets parameters for MC section                                    
!-                                                                      
      USE config_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      REAL uerte (maxw) 
      REAL verte (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, iianz, jjanz, kkanz, is, js, ic, i, j 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.2) then 
            CALL do_cap (cpara (1) ) 
!                                                                       
!------ --- 'set atom': setting of atoms to be taken for MC             
!                                                                       
            IF (cpara (1) (1:2) .eq.'AT'.or.cpara (1) (1:3) .eq.'MOL')  &
            then                                                        
               IF (ianz.eq.3) then 
                  mo_sel_atom = (cpara (1) (1:2) .eq.'AT') 
                  mo_atom (1) = cpara (2) (1:lpara (2) ) 
                  mo_atom (2) = cpara (3) (1:lpara (3) ) 
               ELSEIF (ianz.eq.4) then 
                  mo_sel_atom = (cpara (1) (1:2) .eq.'AT') 
                  mo_atom (1) = cpara (2) (1:lpara (2) ) 
                  mo_atom (2) = cpara (3) (1:lpara (3) ) 
                  mo_atom (3) = cpara (4) (1:lpara (4) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- 'set angle' : setting of neighbouring angles                
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'AN') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_angle (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set envir' : setting of neighbouring environment           
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'ENV') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_envir (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set const': setting interactions and feedback factor       
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'CON') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               IF (ianz.eq.2.or.ianz.eq.3) then 
                  ic = nint (werte (1) ) 
                  IF (ic.ge.0.and.ic.le.chem_ncor) then 
                     mo_const (ic) = werte (2) 
                     IF (ianz.eq.3) mo_cfac (ic) = werte (3) 
                  ELSE 
                     ier_num = - 14 
                     ier_typ = ER_CHEM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- 'set cyc' : setting number of MC moves                      
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'CY') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_cyc = nint (werte (1) ) 
!                                                                       
!------ --- 'set energy' : setting of energy expression                 
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'ENE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
!                                                                       
               IF (ianz.ge.1) then 
                  CALL do_cap (cpara (1) ) 
                  IF (cpara (1) (1:2) .eq.'CO') then 
                     mo_energy = MC_OCC 
                  ELSEIF (cpara (1) (1:2) .eq.'CD') then 
                     mo_energy = MC_DISP 
                  ELSEIF (cpara (1) (1:2) .eq.'DI') then 
                     mo_energy = MC_SPRING 
                  ELSEIF (cpara (1) (1:2) .eq.'AN') then 
                     mo_energy = MC_ANGLE 
                  ELSE 
                     ier_num = - 1 
                     ier_typ = ER_MMC 
                     RETURN 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- 'set feed' : setting display/feedback intervall             
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'FE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_feed = nint (werte (1) ) 
!                                                                       
!------ --- 'set mode': sets operation mode for MC                      
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'MOD') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL rmc_set_mode (mo_mode, mo_local, ianz, cpara, lpara,&
               maxw)                                                    
!                                                                       
!------ --- 'set move': sets maxmove for shift MC mode                  
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'MOV') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL rmc_set_move (mo_maxmove, mo_sel_atom, ianz, cpara, &
               werte, lpara, maxw)                                      
!                                                                       
!------ --- 'set neig': setting correlation determination method        
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'NE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_neig (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set target' : setting of target correlations               
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'TA') then 
               IF (ianz.ge.3) then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (1, cpara, lpara, werte, maxw) 
                  ic = nint (werte (1) ) 
                  IF (ic.gt.0.and.ic.le.chem_ncor) then 
                     IF (ianz.eq.2) then 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (1, cpara, lpara, werte, maxw) 
                        mo_target_corr (ic) = werte (1) 
                     ELSEIF (ianz.eq.4) then 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (1, cpara, lpara, werte, maxw) 
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                        CALL mc_set_disp (ic, is, js, werte (1) ) 
                        ENDDO 
                        ENDDO 
                     ELSEIF (ianz.eq.5) then 
!                                                                       
!     ------------Three atom names are listed on the traget instruction 
!                 Required for angular correlations                     
!                                                                       
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        kkanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (kkanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (1, cpara, lpara, werte, maxw) 
                        DO i = 1, jjanz 
                        DO j = 1, kkanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                        CALL mc_set_disp (ic, is, js, werte (1) ) 
                        ENDDO 
                        ENDDO 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 14 
                     ier_typ = ER_CHEM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- 'set temp' : setting kT for MC simulation                   
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'TE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_kt = werte (1) 
!                                                                       
!------ --- 'set vector' : setting of neighbouring vectors              
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'VE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_vec (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- no valid subcommand                                         
!                                                                       
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE mc_set                         
!*****7*****************************************************************
      SUBROUTINE mc_set_disp (ic, is, js, dist) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
      USE config_mod 
      USE crystal_mod 
      USE mc_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER ic, is, js, ii, jj 
      REAL dist 
!                                                                       
                                                                        
      IF (is.ne. - 1.and.js.ne. - 1) then 
         mo_disp (ic, is, js) = dist 
         mo_disp (ic, js, is) = dist 
      ELSEIF (is.eq. - 1.and.js.ne. - 1) then 
         DO ii = 0, cr_nscat 
         mo_disp (ic, ii, js) = dist 
         mo_disp (ic, js, ii) = dist 
         ENDDO 
      ELSEIF (is.ne. - 1.and.js.eq. - 1) then 
         DO ii = 0, cr_nscat 
         mo_disp (ic, ii, is) = dist 
         mo_disp (ic, is, ii) = dist 
         ENDDO 
      ELSE 
         DO ii = 0, cr_nscat 
         DO jj = 0, cr_nscat 
         mo_disp (ic, ii, jj) = dist 
         mo_disp (ic, jj, ii) = dist 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE mc_set_disp                    
!*****7*****************************************************************
      SUBROUTINE mc_run 
!+                                                                      
!     This is the main MC routine                                       
!-                                                                      
      USE config_mod 
      USE mc_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      CHARACTER(1) a 
!                                                                       
      IF (mo_cyc.le.0) then 
         ier_num = - 4 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
      IF (mo_feed.le.0) then 
         ier_num = - 5 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ Correlation occupational mode                                   
!                                                                       
      IF (mo_energy.eq.MC_OCC) then 
         IF (mo_sel_atom) then 
            CALL mc_run_occ 
         ELSE 
            CALL mc_run_occ_mol 
         ENDIF 
!                                                                       
!------ Correlation displacement mode                                   
!                                                                       
      ELSEIF (mo_energy.eq.MC_DISP) then 
         IF (mo_sel_atom) then 
            CALL mc_run_dis 
         ELSE 
            CALL mc_run_dis_mol 
         ENDIF 
!                                                                       
!------ Spring mode                                                     
!                                                                       
      ELSEIF (mo_energy.eq.MC_SPRING) then 
         IF (mo_sel_atom) then 
            CALL mc_run_spring 
         ELSE 
            CALL mc_run_spring_mol 
         ENDIF 
!                                                                       
!------ Angular mode                                                    
!                                                                       
      ELSEIF (mo_energy.eq.MC_ANGLE) then 
         IF (mo_sel_atom) then 
            CALL mc_run_angle 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 1 
         ier_typ = ER_MMC 
      ENDIF 
!                                                                       
      CALL update_cr_dim 
!                                                                       
      END SUBROUTINE mc_run                         
!*****7*****************************************************************
      SUBROUTINE mc_run_occ 
!+                                                                      
!     This is the MC routine for occupational correlations              
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE modify_func_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(1024) cpara (MAXSCAT) 
      INTEGER lpara (MAXSCAT) 
      REAL werte (MAXSCAT), wwerte (MAXSCAT) 
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER isel (chem_max_atom) 
      INTEGER latom (2) 
      INTEGER iscat, is1, is2, iz1 (3), iz2 (3) 
      INTEGER i, natoms, ianz, janz, l 
      INTEGER lbeg (3) 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_occ 
!     LOGICAL atom_allowed 
      INTEGER len_str 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      maxw = MAXSCAT
      IF (mo_mode.ne.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
      lbeg (1) = - 1 
      lbeg (2) = 0 
      lbeg (3) = 0 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), cr_natoms 
!                                                                       
!------ Setup atom arrays                                               
!                                                                       
      ianz = 1 
      cpara (1) = mo_atom (1) 
      lpara (1) = latom (1) 
!                                                                       
      CALL get_iscat (ianz, cpara, lpara, werte, maxw, .false.) 
!                                                                       
      janz = 1 
      cpara (1) = mo_atom (2) 
      lpara (1) = latom (2) 
!                                                                       
      CALL get_iscat (janz, cpara, lpara, wwerte, maxw, .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Obtain starting interactions                                    
!                                                                       
      l = 2 
      CALL chem_corr_occ (l, mo_atom, latom, 2, .false., lbeg) 
      DO i = 1, chem_ncor 
      mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (i) &
      - mo_ach_corr (i) ) / 2.)                                         
      ENDDO 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
      natoms = 2 
      CALL rmc_select (mo_local, natoms, isel, iz1, iz2, is1, is2) 
      laccept = atom_allowed (isel (1), werte, ianz, maxw)              &
      .and.atom_allowed (isel (2), wwerte, janz, maxw) .and.cr_iscat (  &
      isel (1) ) .ne.cr_iscat (isel (2) )                               
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy                                                 
!                                                                       
         e_old = mc_energy_occ (natoms, isel, ianz, werte, maxw) 
!                                                                       
         iscat = cr_iscat (isel (2) ) 
         cr_iscat (isel (2) ) = cr_iscat (isel (1) ) 
         cr_iscat (isel (1) ) = iscat 
!                                                                       
         e_new = mc_energy_occ (natoms, isel, ianz, werte, maxw) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
            iscat = cr_iscat (isel (2) ) 
            cr_iscat (isel (2) ) = cr_iscat (isel (1) ) 
            cr_iscat (isel (1) ) = iscat 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         l = 2 
         CALL chem_corr_occ (l, mo_atom, latom, 2, .false., lbeg) 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         DO i = 1, chem_ncor 
         WRITE (output_io, 2100) i, mo_const (i), mo_ach_corr (i),      &
         mo_target_corr (i)                                             
         mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (&
         i) - mo_ach_corr (i) ) / 2.)                                   
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      WRITE (output_io, 3000) 
      l = 2 
      CALL chem_corr_occ (l, mo_atom, latom, 2, .false., lbeg) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      DO i = 1, chem_ncor 
      WRITE (output_io, 2100) i, mo_const (i),  mo_ach_corr (i),  mo_tar&
     &get_corr (i)                                                      
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' atoms')                             
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2100 FORMAT ('    Neig. ',i3,'  const: ',f10.3,'   achieved: ',f6.3,   &
     &                     '  target: ',f6.3)                           
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_occ                     
!*****7*****************************************************************
      SUBROUTINE mc_run_occ_mol 
!+                                                                      
!     This is the MC routine for occupational correlations              
!     Molecules version ..                                              
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      REAL mol (2) 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER i, imol (2), mtyp (2), amol (2), latom (2), l 
      INTEGER lbeg (3) 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_occ_mol, ran1 
      INTEGER len_str 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      IF (mo_mode.ne.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
      lbeg (1) = - 1 
      lbeg (2) = 0 
      lbeg (3) = 0 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), mole_num_mole 
!                                                                       
!------ Get molecule type numbers                                       
!                                                                       
      CALL ber_params (2, mo_atom, latom, mol, 2) 
      IF (ier_num.ne.0) return 
!                                                                       
      amol (1) = nint (mol (1) ) 
      amol (2) = nint (mol (2) ) 
!                                                                       
!------ Obtain starting interactions                                    
!                                                                       
      l = 2 
      CALL chem_corr_occ_mol (l, mo_atom, latom, 2, .false., lbeg) 
      DO i = 1, chem_ncor 
      mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (i) &
      - mo_ach_corr (i) ) / 2.)                                         
      ENDDO 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
      imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
      imol (2) = int (ran1 (idum) * mole_num_mole) + 1 
      mtyp (1) = mole_type (imol (1) ) 
      mtyp (2) = mole_type (imol (2) ) 
      laccept = (mtyp (1) .eq.amol (1) .or.mtyp (1) .eq.amol (2) )      &
      .and. (mtyp (2) .eq.amol (1) .or.mtyp (2) .eq.amol (2) ) .and. (  &
      mtyp (1) .ne.mtyp (2) )                                           
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy                                                 
!                                                                       
         e_old = mc_energy_occ_mol (2, imol, amol) 
         CALL do_swap_mole (imol (1), imol (2), .true.) 
         e_new = mc_energy_occ_mol (2, imol, amol) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
            CALL do_swap_mole (imol (1), imol (2), .true.) 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         l = 2 
         CALL chem_corr_occ_mol (l, mo_atom, latom, 2, .false., lbeg) 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         DO i = 1, chem_ncor 
         WRITE (output_io, 2100) i, mo_const (i), mo_ach_corr (i),      &
         mo_target_corr (i)                                             
         mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (&
         i) - mo_ach_corr (i) ) / 2.)                                   
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      WRITE (output_io, 3000) 
      l = 2 
      CALL chem_corr_occ_mol (l, mo_atom, latom, 2, .false., lbeg) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      DO i = 1, chem_ncor 
      WRITE (output_io, 2100) i, mo_const (i),  mo_ach_corr (i),  mo_tar&
     &get_corr (i)                                                      
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' molecules')                         
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2100 FORMAT ('    Neig. ',i3,'  const: ',f10.3,'   achieved: ',f6.3,   &
     &                     '  target: ',f6.3)                           
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_occ_mol                 
!*****7*****************************************************************
      SUBROUTINE mc_run_dis 
!+                                                                      
!     This is the MC routine for displacement correlations,             
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE rmc_mod 
      USE modify_mod
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(1024) cpara (MAXSCAT) 
      INTEGER lpara (MAXSCAT) 
      REAL werte (MAXSCAT), wwerte (MAXSCAT) 
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      REAL disp1, disp2 
      REAL idir (3), jdir (3) 
      REAL disp (3, rmc_max_atom) 
      REAL pos (3, rmc_max_atom) 
      REAL rdi (chem_max_cor) 
      REAL rdj (chem_max_cor) 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER isel (chem_max_atom) 
      INTEGER latom (2) 
      INTEGER lbeg (3) 
      INTEGER is1, is2, iz1 (3), iz2 (3) 
      INTEGER i, j, k, natoms, ianz, janz, ia, l 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_dis, ran1, gasdev, skalpro 
      INTEGER len_str 
!     LOGICAL atom_allowed 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      maxw = MAXSCAT
      IF (mo_mode.eq.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
      lbeg (1) = - 1 
      lbeg (2) = 0 
      lbeg (3) = 0 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), cr_natoms 
!                                                                       
!------ Setup atom arrays                                               
!                                                                       
      ianz = 1 
      cpara (1) = mo_atom (1) 
      lpara (1) = latom (1) 
!                                                                       
      CALL get_iscat (ianz, cpara, lpara, werte, maxw, .false.) 
!     call del_params(1,ianz,cpara,lpara,maxw)                          
!                                                                       
      janz = 1 
      cpara (1) = mo_atom (2) 
      lpara (1) = latom (2) 
!                                                                       
      CALL get_iscat (janz, cpara, lpara, wwerte, maxw, .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Normalise correlation directions                                
!                                                                       
      CALL chem_aver (.false.) 
      DO k = 1, chem_ncor 
      DO i = 1, 3 
      idir (i) = chem_dir (i, 1, k) 
      jdir (i) = chem_dir (i, 2, k) 
      ENDDO 
!                                                                       
      rdi (k) = skalpro (idir, idir, cr_gten) 
      rdj (k) = skalpro (jdir, jdir, cr_gten) 
      rdi (k) = sqrt (rdi (k) ) 
      rdj (k) = sqrt (rdj (k) ) 
      ENDDO 
!                                                                       
!------ Obtain starting interactions                                    
!                                                                       
      l = 2 
      CALL chem_corr_dis (l, mo_atom, latom, 2, .false., lbeg) 
      DO i = 1, chem_ncor 
      mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (i) &
      - mo_ach_corr (i) ) / 2.)                                         
      ENDDO 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
!------ - selecting sites                                               
!                                                                       
      IF (mo_mode.eq.rmc_mode_shift) then 
         natoms = 1 
   10    CONTINUE 
         isel (1) = int (ran1 (idum) * cr_natoms) + 1 
         IF (isel (1) .gt.cr_natoms.or.isel (1) .lt.1) goto 10 
         laccept = atom_allowed (isel (1), werte, ianz, maxw) 
         IF (laccept) then 
            DO i = 1, 3 
            disp (i, 1) = gasdev (mo_maxmove (i, cr_iscat (isel (1) ) ) &
            )                                                           
            ENDDO 
         ENDIF 
!                                                                       
      ELSEIF (mo_mode.eq.rmc_mode_swdisp) then 
         natoms = 2 
         CALL rmc_select (mo_local, natoms, isel, iz1, iz2, is1, is2) 
         laccept = atom_allowed (isel (1), werte, ianz, maxw)           &
         .and.atom_allowed (isel (2), wwerte, janz, maxw)               
         IF (laccept) then 
            DO j = 1, 3 
            disp1 = cr_pos (j, isel (1) ) - chem_ave_pos (j, is1)       &
            - float (iz1 (j) - 1) - cr_dim0 (j, 1)                      
            disp2 = cr_pos (j, isel (2) ) - chem_ave_pos (j, is2)       &
            - float (iz2 (j) - 1) - cr_dim0 (j, 1)                      
            disp (j, 1) = - disp1 + disp2 
            disp (j, 2) = - disp2 + disp1 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!----- -- Try move                                                      
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy (energy calculates DELTA E !!)                  
!                                                                       
         e_old = 0.0 
!                                                                       
         DO i = 1, 3 
         DO ia = 1, natoms 
         pos (i, ia) = cr_pos (i, isel (ia) ) 
         cr_pos (i, isel (ia) ) = cr_pos (i, isel (ia) ) + disp (i, ia) 
         ENDDO 
         ENDDO 
!                                                                       
         e_new = mc_energy_dis (natoms, isel, disp, rdi, rdj) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
            DO i = 1, 3 
            DO ia = 1, natoms 
            cr_pos (i, isel (ia) ) = pos (i, ia) 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         l = 2 
         CALL chem_corr_dis (l, mo_atom, latom, 2, .false., lbeg) 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         DO i = 1, chem_ncor 
         WRITE (output_io, 2100) i, mo_const (i), mo_ach_corr (i),      &
         mo_target_corr (i)                                             
         mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (&
         i) - mo_ach_corr (i) ) / 2.)                                   
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      WRITE (output_io, 3000) 
      l = 2 
      CALL chem_corr_dis (l, mo_atom, latom, 2, .false., lbeg) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      DO i = 1, chem_ncor 
      WRITE (output_io, 2100) i, mo_const (i),  mo_ach_corr (i),  mo_tar&
     &get_corr (i)                                                      
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' atoms')                             
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2100 FORMAT ('    Neig. ',i3,'  const: ',f10.3,'   achieved: ',f6.3,   &
     &                     '  target: ',f6.3)                           
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_dis                     
!*****7*****************************************************************
      SUBROUTINE mc_run_dis_mol 
!+                                                                      
!     This is the MC routine for displacement correlations,             
!     Molecule mode ..                                                  
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      REAL disp1, disp2 
      REAL idir (3), jdir (3) 
      REAL disp (3, rmc_max_atom) 
      REAL pos (3, rmc_max_atom) 
      REAL rdi (chem_max_cor) 
      REAL rdj (chem_max_cor) 
      REAL mol (2) 
      INTEGER imol (2), amol (2), mtyp (2), nmol 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER latom (2) 
      INTEGER is1, is2, iz1 (3), iz2 (3) 
      INTEGER i, j, k, ia, l 
      INTEGER is, js, ii, jj, im 
      INTEGER lbeg (3) 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_dis_mol, ran1, gasdev, skalpro 
      INTEGER len_str 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      IF (mo_mode.eq.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
      lbeg (1) = - 1.0 
      lbeg (2) = 0.0 
      lbeg (3) = 0.0 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), mole_num_mole 
!                                                                       
!------ Get molecule type numbers                                       
!                                                                       
      CALL ber_params (2, mo_atom, latom, mol, 2) 
      IF (ier_num.ne.0) return 
!                                                                       
      amol (1) = nint (mol (1) ) 
      amol (2) = nint (mol (2) ) 
!                                                                       
!------ Normalise correlation directions                                
!                                                                       
      DO k = 1, chem_ncor 
      DO i = 1, 3 
      idir (i) = chem_dir (i, 1, k) 
      jdir (i) = chem_dir (i, 2, k) 
      ENDDO 
!                                                                       
      rdi (k) = skalpro (idir, idir, cr_gten) 
      rdj (k) = skalpro (jdir, jdir, cr_gten) 
      rdi (k) = sqrt (rdi (k) ) 
      rdj (k) = sqrt (rdj (k) ) 
      ENDDO 
!                                                                       
!------ Obtain starting interactions                                    
!                                                                       
      l = 2 
      CALL chem_corr_dis_mol (l, mo_atom, latom, 2, .false., lbeg) 
      DO i = 1, chem_ncor 
      mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (i) &
      - mo_ach_corr (i) ) / 2.)                                         
      ENDDO 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
!------ - selecting sites                                               
!                                                                       
      IF (mo_mode.eq.rmc_mode_shift) then 
         nmol = 1 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         mtyp (1) = mole_type (imol (1) ) 
         laccept = (mtyp (1) .eq.amol (1) .or.mtyp (1) .eq.amol (2) ) 
         IF (laccept) then 
            DO i = 1, 3 
            disp (i, 1) = gasdev (mo_maxmove (i, mtyp (1) ) ) 
            ENDDO 
         ENDIF 
!                                                                       
      ELSEIF (mo_mode.eq.rmc_mode_swdisp) then 
         nmol = 2 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         imol (2) = int (ran1 (idum) * mole_num_mole) + 1 
         mtyp (1) = mole_type (imol (1) ) 
         mtyp (2) = mole_type (imol (2) ) 
         is = mole_cont (mole_off (imol (1) ) ) + 1 
         js = mole_cont (mole_off (imol (2) ) ) + 1 
         laccept = (mtyp (1) .eq.amol (1) .or.mtyp (1) .eq.amol (2) )   &
         .and. (mtyp (2) .eq.amol (1) .or.mtyp (2) .eq.amol (2) )       
!                                                                       
         IF (laccept) then 
            CALL indextocell (is, iz1, is1) 
            CALL indextocell (js, iz2, is2) 
            DO j = 1, 3 
            disp1 = cr_pos (j, is) - chem_ave_pos (j, is1) - float (iz1 &
            (j) - 1) - cr_dim0 (j, 1)                                   
            disp2 = cr_pos (j, js) - chem_ave_pos (j, is2) - float (iz2 &
            (j) - 1) - cr_dim0 (j, 1)                                   
            disp (j, 1) = - disp1 + disp2 
            disp (j, 2) = - disp2 + disp1 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!----- -- Try move                                                      
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy (energy calculates DELTA E !!)                  
!                                                                       
         e_old = 0.0 
!                                                                       
         jj = 0 
!                                                                       
         DO ia = 1, nmol 
         DO im = 1, mole_len (imol (ia) ) 
         ii = mole_cont (mole_off (imol (ia) ) + im) 
         jj = jj + 1 
         DO i = 1, 3 
         pos (i, jj) = cr_pos (i, ii) 
         cr_pos (i, ii) = cr_pos (i, ii) + disp (i, ia) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
         e_new = mc_energy_dis_mol (nmol, imol, disp, rdi, rdj) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
            jj = 0 
!                                                                       
            DO ia = 1, nmol 
            DO im = 1, mole_len (imol (ia) ) 
            ii = mole_cont (mole_off (imol (ia) ) + im) 
            jj = jj + 1 
            DO i = 1, 3 
            cr_pos (i, ii) = pos (i, jj) 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         l = 2 
         CALL chem_corr_dis_mol (l, mo_atom, latom, 2, .false., lbeg) 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         DO i = 1, chem_ncor 
         WRITE (output_io, 2100) i, mo_const (i), mo_ach_corr (i),      &
         mo_target_corr (i)                                             
         mo_const (i) = mo_const (i) - mo_cfac (i) * ( (mo_target_corr (&
         i) - mo_ach_corr (i) ) / 2.)                                   
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      WRITE (output_io, 3000) 
      l = 2 
      CALL chem_corr_dis (l, mo_atom, latom, 2, .false., lbeg) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      DO i = 1, chem_ncor 
      WRITE (output_io, 2100) i, mo_const (i),  mo_ach_corr (i),  mo_tar&
     &get_corr (i)                                                      
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' molecules')                         
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2100 FORMAT ('    Neig. ',i3,'  const: ',f10.3,'   achieved: ',f6.3,   &
     &                     '  target: ',f6.3)                           
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_dis_mol                 
!*****7*****************************************************************
      SUBROUTINE mc_run_spring 
!+                                                                      
!     This is the MC routine for distortions (Hook's law)               
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE modify_func_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA =  20 ! A command requires at leaset these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: wwerte
!
      CHARACTER(9) at_name_i, at_name_j, at_name 
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      REAL disp1, disp2 
      REAL disp (3, rmc_max_atom) 
      REAL pos (3, rmc_max_atom) 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER isel (chem_max_atom) 
      INTEGER latom (2) 
      INTEGER is1, is2, iz1 (3), iz2 (3) 
      INTEGER i, j, k, natoms, ianz, janz, ia, l 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_spr, ran1, gasdev 
      INTEGER len_str 
!     LOGICAL atom_allowed 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      maxw = MAX(min_PARA,MAXSCAT+1)
!
      IF (mo_mode.eq.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), cr_natoms 
!                                                                       
!------ Setup atom arrays                                               
!                                                                       
      ianz = 1 
      janz = 1 
      cpara (1) = mo_atom (1) 
      cpara (2) = mo_atom (2) 
      lpara (1) = latom (1) 
      lpara (2) = latom (2) 
!                                                                       
      CALL get_iscat (ianz, cpara, lpara, werte, maxw, .false.) 
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL get_iscat (janz, cpara, lpara, wwerte, maxw, .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
      l = 2 
      cpara (1) = 'all' 
      cpara (2) = 'all' 
      lpara (1) = 3 
      lpara (2) = 3 
!                                                                       
      DO i = 1, chem_ncor 
      IF (mo_const (i) .eq.0.0) mo_const (i) = mo_cfac (i) * 1.0 
      ENDDO 
!                                                                       
      CALL chem_aver (.false.) 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
!------ - selecting sites                                               
!                                                                       
      IF (mo_mode.eq.rmc_mode_shift) then 
         natoms = 1 
   10    CONTINUE 
         isel (1) = int (ran1 (idum) * cr_natoms) + 1 
         IF (isel (1) .gt.cr_natoms.or.isel (1) .lt.1) goto 10 
         laccept = atom_allowed (isel (1), werte, ianz, maxw)           &
         .or.atom_allowed (isel (1), wwerte, ianz, maxw)                
         IF (laccept) then 
            DO i = 1, 3 
            disp (i, 1) = gasdev (mo_maxmove (i, cr_iscat (isel (1) ) ) &
            )                                                           
            ENDDO 
         ENDIF 
!                                                                       
      ELSEIF (mo_mode.eq.rmc_mode_swdisp) then 
         natoms = 2 
         CALL rmc_select (mo_local, natoms, isel, iz1, iz2, is1, is2) 
         laccept = atom_allowed (isel (1), werte, ianz, maxw)           &
         .and.atom_allowed (isel (2), wwerte, janz, maxw)               
         IF (laccept) then 
            DO j = 1, 3 
            disp1 = cr_pos (j, isel (1) ) - chem_ave_pos (j, is1)       &
            - float (iz1 (j) - 1) - cr_dim0 (j, 1)                      
            disp2 = cr_pos (j, isel (2) ) - chem_ave_pos (j, is2)       &
            - float (iz2 (j) - 1) - cr_dim0 (j, 1)                      
            disp (j, 1) = - disp1 + disp2 
            disp (j, 2) = - disp2 + disp1 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!----- -- Try move                                                      
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy                                                 
!                                                                       
         e_old = mc_energy_spr (natoms, isel) 
!                                                                       
         DO i = 1, 3 
         DO ia = 1, natoms 
         pos (i, ia) = cr_pos (i, isel (ia) ) 
         cr_pos (i, isel (ia) ) = cr_pos (i, isel (ia) ) + disp (i, ia) 
         ENDDO 
         ENDDO 
!                                                                       
         e_new = mc_energy_spr (natoms, isel) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
            DO i = 1, 3 
            DO ia = 1, natoms 
            cr_pos (i, isel (ia) ) = pos (i, ia) 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         CALL chem_disp (l, cpara, lpara, werte, maxw, .false.) 
         DO i = 0, cr_nscat 
         DO j = i, cr_nscat 
         at_name_i = at_name (i) 
         at_name_j = at_name (j) 
         DO k = 1, chem_ncor 
         IF (mo_disp (k, i, j) .ne.0.0) then 
            WRITE (output_io, 2300) k, at_name_i, at_name_j, mo_disp (k,&
            i, j), chem_disp_ave (k, i, j), chem_disp_sig (k, i, j)     
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      WRITE (output_io, 3000) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      CALL chem_disp (l, cpara, lpara, werte, maxw, .false.) 
      DO i = 0, cr_nscat 
      DO j = i, cr_nscat 
      at_name_i = at_name (i) 
      at_name_j = at_name (j) 
      DO k = 1, chem_ncor 
      IF (mo_disp (k, i, j) .ne.0.0) then 
         WRITE (output_io, 2300) k, at_name_i, at_name_j, mo_disp (k, i,&
         j), chem_disp_ave (k, i, j), chem_disp_sig (k, i, j)           
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' atoms')                             
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2300 FORMAT     ('    Neig. ',i3,' : ',a9,'- ',a9,'target: ',          &
     &                   f7.3,' ach.: ',f7.3,' +- ',f7.3)               
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_spring                  
!*****7*****************************************************************
      SUBROUTINE mc_run_spring_mol 
!+                                                                      
!     This is the MC routine for distortions (Hook's law)               
!     Molecules version ..                                              
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      REAL disp1, disp2 
      REAL disp (3, rmc_max_atom) 
      REAL pos (3, rmc_max_atom) 
      REAL mol (2) 
      REAL werte (2) 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER amol (2), imol (2), mtyp (2), nmol 
      INTEGER latom (2) 
      INTEGER is1, is2, iz1 (3), iz2 (3) 
      INTEGER i, j, k, is, js, im, ia, ii, jj, l 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_spr_mol, ran1, gasdev 
      INTEGER len_str 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      IF (mo_mode.eq.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      l = 2 
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), mole_num_mole 
!                                                                       
!------ Get molecule type numbers                                       
!                                                                       
      CALL ber_params (2, mo_atom, latom, mol, 2) 
      IF (ier_num.ne.0) return 
!                                                                       
      amol (1) = nint (mol (1) ) 
      amol (2) = nint (mol (2) ) 
!                                                                       
      DO i = 1, chem_ncor 
      IF (mo_const (i) .eq.0.0) mo_const (i) = mo_cfac (i) * 1.0 
      ENDDO 
!                                                                       
      CALL chem_aver (.false.) 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
!------ - selecting sites                                               
!                                                                       
      IF (mo_mode.eq.rmc_mode_shift) then 
         nmol = 1 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         mtyp (1) = mole_type (imol (1) ) 
         laccept = (mtyp (1) .eq.amol (1) .or.mtyp (1) .eq.amol (2) ) 
         IF (laccept) then 
            DO i = 1, 3 
            disp (i, 1) = gasdev (mo_maxmove (i, mtyp (1) ) ) 
            ENDDO 
         ENDIF 
!                                                                       
      ELSEIF (mo_mode.eq.rmc_mode_swdisp) then 
         nmol = 2 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         imol (2) = int (ran1 (idum) * mole_num_mole) + 1 
         mtyp (1) = mole_type (imol (1) ) 
         mtyp (2) = mole_type (imol (2) ) 
         is = mole_cont (mole_off (imol (1) ) ) + 1 
         js = mole_cont (mole_off (imol (2) ) ) + 1 
         laccept = (mtyp (1) .eq.amol (1) .or.mtyp (1) .eq.amol (2) )   &
         .and. (mtyp (2) .eq.amol (1) .or.mtyp (2) .eq.amol (2) )       
!                                                                       
         IF (laccept) then 
            CALL indextocell (is, iz1, is1) 
            CALL indextocell (js, iz2, is2) 
            DO j = 1, 3 
            disp1 = cr_pos (j, is) - chem_ave_pos (j, is1) - float (iz1 &
            (j) - 1) - cr_dim0 (j, 1)                                   
            disp2 = cr_pos (j, js) - chem_ave_pos (j, is2) - float (iz2 &
            (j) - 1) - cr_dim0 (j, 1)                                   
            disp (j, 1) = - disp1 + disp2 
            disp (j, 2) = - disp2 + disp1 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!----- -- Try move                                                      
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy                                                 
!                                                                       
         e_old = mc_energy_spr_mol (nmol, imol) 
!                                                                       
         jj = 0 
!                                                                       
         DO ia = 1, nmol 
         DO im = 1, mole_len (imol (ia) ) 
         ii = mole_cont (mole_off (imol (ia) ) + im) 
         jj = jj + 1 
         DO i = 1, 3 
         pos (i, jj) = cr_pos (i, ii) 
         cr_pos (i, ii) = cr_pos (i, ii) + disp (i, ia) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
         e_new = mc_energy_spr_mol (nmol, imol) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
            jj = 0 
!                                                                       
            DO ia = 1, nmol 
            DO im = 1, mole_len (imol (ia) ) 
            ii = mole_cont (mole_off (imol (ia) ) + im) 
            jj = jj + 1 
            DO i = 1, 3 
            cr_pos (i, ii) = pos (i, jj) 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         CALL chem_disp_mol (l, mo_atom, latom, werte, 2, .false.) 
         DO i = 1, mole_num_type 
         DO j = i, mole_num_type 
         DO k = 1, chem_ncor 
         IF (mo_disp (k, i, j) .ne.0.0) then 
      WRITE (output_io, 2300) k, i, j, mo_disp (k, i, j),  chem_disp_ave&
     & (k, i, j),  chem_disp_sig (k, i, j)                              
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      WRITE (output_io, 3000) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      CALL chem_disp_mol (l, mo_atom, latom, werte, 2, .false.) 
      DO i = 1, mole_num_type 
      DO j = i, mole_num_type 
      DO k = 1, chem_ncor 
      IF (mo_disp (k, i, j) .ne.0.0) then 
      WRITE (output_io, 2300) k, i, j, mo_disp (k, i, j),  chem_disp_ave&
     & (k, i, j),  chem_disp_sig (k, i, j)                              
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' molecules')                         
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2300 FORMAT     ('    Neig. ',i3,' : Mol',i4,'  - Mol',i4,'  target: ',&
     &                   f7.3,' ach.: ',f7.3,' +- ',f7.3)               
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_spring_mol              
!*****7*****************************************************************
      SUBROUTINE mc_run_angle 
!+                                                                      
!     This is the MC routine for angular distortions using displacements
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE modify_func_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(1024) cpara (MAXSCAT) 
      CHARACTER(9) at_name_i, at_name_j, at_name 
      REAL werte (MAXSCAT), wwerte (MAXSCAT), wwwerte (MAXSCAT) 
      REAL verte (MAXSCAT) 
      REAL e_old, e_new 
      REAL start, zeit, seknds 
      REAL disp1, disp2 
      REAL disp (3, 0:CHEM_MAX_NEIG) 
      REAL posz (3) 
      REAL pos (3, 0:CHEM_MAX_NEIG, CHEM_MAX_COR) 
      REAL patom (3, 0:CHEM_MAX_NEIG) 
      INTEGER iatom (0:CHEM_MAX_NEIG) 
      INTEGER jatom (CHEM_MAX_NEIG, CHEM_MAX_COR) 
      INTEGER natom_ic (CHEM_MAX_COR) 
      INTEGER lpara (MAXSCAT) 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER isel (CHEM_MAX_ATOM) 
      INTEGER iselz 
      INTEGER latom (3) 
      INTEGER ic, is1, is2, iz1 (3), iz2 (3) 
      INTEGER i, j, k, natoms, ianz, janz, ia, l, kanz 
      INTEGER natom 
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
!                                                                       
      REAL mc_energy_angle, ran1, gasdev 
      INTEGER len_str 
!     LOGICAL atom_allowed 
!                                                                       
!------ check for valid input and settings                              
!                                                                       
      maxw = MAXSCAT
      IF (mo_mode.eq.rmc_mode_swchem) then 
         ier_num = - 3 
         ier_typ = ER_MMC 
         RETURN 
      ENDIF 
!                                                                       
!------ reset some counters                                             
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      done = .true. 
!                                                                       
      latom (1) = len_str (mo_atom (1) ) 
      latom (2) = len_str (mo_atom (2) ) 
      latom (3) = len_str (mo_atom (3) ) 
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), cr_natoms 
!                                                                       
!------ Setup atom arrays                                               
!                                                                       
      ianz = 1 
      cpara (1) = mo_atom (1) 
      lpara (1) = latom (1) 
      CALL get_iscat (ianz, cpara, lpara, werte, maxw, .false.) 
      janz = 1 
      cpara (1) = mo_atom (2) 
      lpara (1) = latom (2) 
      CALL get_iscat (janz, cpara, lpara, wwerte, maxw, .false.) 
      kanz = 1 
      cpara (1) = mo_atom (3) 
      lpara (1) = latom (3) 
      CALL get_iscat (kanz, cpara, lpara, wwwerte, maxw, .false.) 
                                                                        
!                                                                       
!     call del_params(1,ianz,cpara,lpara,maxw)                          
!     call del_params(1,ianz,cpara,lpara,maxw)                          
      IF (ier_num.ne.0) return 
!                                                                       
      l = 3 
      cpara (1) = mo_atom (1) 
      cpara (2) = mo_atom (2) 
      cpara (3) = mo_atom (3) 
      lpara (1) = latom (1) 
      lpara (2) = latom (2) 
      lpara (3) = latom (3) 
!                                                                       
      DO i = 1, chem_ncor 
      IF (mo_const (i) .eq.0.0) mo_const (i) = mo_cfac (i) * 1.0 
      ENDDO 
!                                                                       
!     call chem_aver(.false.)                                           
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
      l = 3 
      cpara (1) = mo_atom (1) 
      cpara (2) = mo_atom (2) 
      cpara (3) = mo_atom (3) 
      lpara (1) = latom (1) 
      lpara (2) = latom (2) 
      lpara (3) = latom (3) 
!                                                                       
!                                                                       
!------ - selecting sites                                               
!                                                                       
      natoms = 1 
   10 CONTINUE 
      isel (1) = int (ran1 (idum) * cr_natoms) + 1 
      iselz = isel (1) 
      IF (isel (1) .gt.cr_natoms.or.isel (1) .lt.1) goto 10 
      laccept = atom_allowed (iselz, werte, ianz, maxw)                 &
      .or.atom_allowed (iselz, wwerte, ianz, maxw) .or.atom_allowed (   &
      iselz, wwwerte, ianz, maxw)                                       
!       if (laccept) then                                               
!         do i=1,3                                                      
!                 disp(i,1) = gasdev(mo_maxmove(i,cr_iscat(isel(1))))   
!         ENDDO                                                         
!       endif                                                           
!                                                                       
!----- -- Try move                                                      
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!------ --- Calc energy                                                 
!                                                                       
         e_old = mc_energy_angle (natoms, isel) 
!                                                                       
!     ----Modify the central atom                                       
!                                                                       
         DO i = 1, 3 
         disp (i, 0) = gasdev (mo_maxmove (i, cr_iscat (iselz) ) ) 
         posz (i) = cr_pos (i, iselz) 
         cr_pos (i, iselz) = cr_pos (i, iselz) + disp (i, 0) 
         ENDDO 
!                                                                       
         e_new = mc_energy_angle (natoms, isel) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
         IF (.not.laccept) then 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
            DO i = 1, 3 
            cr_pos (i, iselz) = posz (i) 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      loop = (itry.lt.mo_cyc) 
!                                                                       
      IF (igen.gt.1000 * mo_feed.and.itry.eq.0) then 
         ier_num = - 2 
         ier_typ = ER_MMC 
         loop = .false. 
      ENDIF 
!                                                                       
!-------  Feedback ?                                                    
!                                                                       
      IF (mod (itry, mo_feed) .eq.0.and..not.done.and.loop) then 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
         CALL chem_angle (l, cpara, lpara, verte, wwerte, wwwerte, maxw,&
         .false.)                                                       
         DO i = 0, cr_nscat 
         DO j = i, cr_nscat 
         at_name_i = at_name (i) 
         at_name_j = at_name (j) 
         DO k = 1, chem_ncor 
         IF (mo_disp (k, i, j) .ne.0.0) then 
            WRITE (output_io, 2300) k, at_name_i, at_name_j, mo_disp (k,&
            i, j), chem_disp_ave (k, i, j), chem_disp_sig (k, i, j)     
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
         done = .true. 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
      l = 3 
      cpara (1) = mo_atom (1) 
      cpara (2) = mo_atom (2) 
      cpara (3) = mo_atom (3) 
      lpara (1) = latom (1) 
      lpara (2) = latom (2) 
      lpara (3) = latom (3) 
!                                                                       
      WRITE (output_io, 3000) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      CALL chem_angle (l, cpara, lpara, verte, wwerte, wwwerte, maxw,   &
      .false.)                                                          
      DO i = 0, cr_nscat 
      DO j = i, cr_nscat 
      at_name_i = at_name (i) 
      at_name_j = at_name (j) 
      DO k = 1, chem_ncor 
      IF (mo_disp (k, i, j) .ne.0.0) then 
         WRITE (output_io, 2300) k, at_name_i, at_name_j, mo_disp (k, i,&
         j), chem_disp_ave (k, i, j), chem_disp_sig (k, i, j)           
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
 1000 FORMAT   (' Running MC simulation ...',//                         &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' atoms')                             
 2000 FORMAT     (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,    &
     &                     ' / ',I7,'  MC moves ')                      
 2300 FORMAT     ('    Neig. ',i3,' : ',a9,'- ',a9,'target: ',          &
     &                   f7.3,' ach.: ',f7.3,' +- ',f7.3)               
 3000 FORMAT     (/,' --- Final configuration ---') 
 4000 FORMAT     (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/ &
     &                     ' Time/cycle   : ',F9.3,' sec',/)            
!                                                                       
      END SUBROUTINE mc_run_angle                   
!*****7*****************************************************************
      SUBROUTINE mc_test (iacc_good, iacc_bad, e_new, e_old, laccept) 
!+                                                                      
!     Tests performed MC move                                           
!-                                                                      
      USE config_mod 
      USE mc_mod 
      IMPLICIT none 
!                                                                       
      include'random.inc' 
!                                                                       
      INTEGER iacc_good, iacc_bad 
      REAL e_new, e_old, e_del 
      LOGICAL laccept 
!                                                                       
      REAL ran1 
!                                                                       
      e_del = e_new - e_old 
      IF (e_del.lt.0) then 
         laccept = .true. 
      ELSE 
         IF (mo_kt.lt.1.0e-10) then 
            laccept = .false. 
         ELSE 
            e_del = exp ( - e_del / mo_kt) 
            e_del = e_del / (1 + e_del) 
            laccept = (e_del.gt.ran1 (idum) ) 
         ENDIF 
      ENDIF 
!                                                                       
      IF (laccept) then 
         IF (e_del.lt.0.0) then 
            iacc_good = iacc_good+1 
         ELSE 
            iacc_bad = iacc_bad+1 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE mc_test                        
!*****7*****************************************************************
      REAL function mc_energy_occ (natoms, isel, ianz, werte, maxw) 
!+                                                                      
!     Calculates the energy for occupational disorder according to      
!                                                                       
!       E = SUM si * (J1*s(i-1) + J2*s(i+1) + ...) with si = +-1.       
!                                                                       
!     Note there is no need for H term because concentration is         
!     fixed using SWCHEM mode.                                          
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_func_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom) 
      INTEGER natoms, ianz, maxw 
      REAL werte (maxw) 
!                                                                       
      INTEGER iatom (0:maxatom), natom 
      INTEGER ic, in, ia, ival1, ival2 
      REAL patom (3, 0:maxatom) 
!                                                                       
!     LOGICAL atom_allowed 
!                                                                       
      mc_energy_occ = 0.0 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, natoms 
      ival1 = - 1 
      IF (atom_allowed (isel (ia), werte, ianz, maxw) ) ival1 = 1 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
      CALL chem_neighbour (isel (ia), ic, iatom, patom, natom, maxatom) 
      IF (natom.ne.0) then 
         DO in = 1, natom 
         ival2 = - 1 
         IF (atom_allowed (iatom (in), werte, ianz, maxw) ) ival2 = 1 
         mc_energy_occ = mc_energy_occ + mo_const (ic) * ival1 * ival2 
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION mc_energy_occ                    
!*****7*****************************************************************
      REAL function mc_energy_occ_mol (ianz, imol, amol) 
!+                                                                      
!     Calculates the energy for occupational disorder for               
!     molecules.                                                        
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxmol 
      PARAMETER (maxmol = chem_max_neig) 
!                                                                       
      INTEGER ianz 
      INTEGER imol (ianz) 
      INTEGER amol (ianz) 
!                                                                       
      INTEGER ineig (0:maxmol), nneig 
      INTEGER ic, in, ia, ival1, ival2 
!                                                                       
      mc_energy_occ_mol = 0.0 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, ianz 
      ival1 = - 1 
      IF (mole_type (imol (ia) ) .eq.amol (1) ) ival1 = 1 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
      CALL chem_neighbour_mol (imol (ia), ic, ineig, nneig, maxmol) 
      IF (nneig.ne.0) then 
         DO in = 1, nneig 
         ival2 = - 1 
         IF (mole_type (ineig (in) ) .eq.amol (1) ) ival2 = 1 
         mc_energy_occ_mol = mc_energy_occ_mol + mo_const (ic) * ival1 *&
         ival2                                                          
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION mc_energy_occ_mol                
!*****7*****************************************************************
      REAL function mc_energy_dis (natoms, isel, disp, rdi, rdj) 
!+                                                                      
!     Calculates the DELTA energy for distortions according to          
!                                                                       
!       DE = SUM k*disp*x                                               
!                                                                       
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom), natoms 
      REAL disp (3, rmc_max_atom) 
      REAL rdi (chem_max_cor) 
      REAL rdj (chem_max_cor) 
!                                                                       
      INTEGER iatom (0:maxatom), natom 
      INTEGER cell (3), site 
      INTEGER i, ic, in, ia 
      REAL patom (3, 0:maxatom) 
      REAL idir (3), jdir (3) 
      REAL delta, dx, u (3), v (3) 
!                                                                       
      REAL skalpro 
!                                                                       
      mc_energy_dis = 0.0 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, natoms 
!                                                                       
      CALL indextocell (isel (ia), cell, site) 
      DO i = 1, 3 
      v (i) = cr_pos (i, isel (ia) ) - chem_ave_pos (i, site) - float ( &
      cell (i) - 1) - cr_dim0 (i, 1)                                    
      v (i) = v (i) + disp (i, ia) 
      ENDDO 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
      IF (chem_ldall (ic) ) then 
         DO i = 1, 3 
         jdir (i) = v (i) 
         ENDDO 
         rdj (ic) = skalpro (jdir, jdir, cr_gten) 
         IF (rdj (ic) .gt.0.0) then 
            rdj (ic) = sqrt (rdj (ic) ) 
         ELSE 
            rdj (ic) = 1.0 
         ENDIF 
         delta = 1.0 
      ELSE 
         DO i = 1, 3 
         idir (i) = chem_dir (i, 1, ic) 
         jdir (i) = chem_dir (i, 2, ic) 
         ENDDO 
         delta = skalpro (v, idir, cr_gten) / rdi (ic) 
      ENDIF 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
      CALL chem_neighbour (isel (ia), ic, iatom, patom, natom, maxatom) 
      IF (natom.ne.0) then 
         DO in = 1, natom 
         CALL indextocell (iatom (in), cell, site) 
         DO i = 1, 3 
         u (i) = cr_pos (i, iatom (in) ) - chem_ave_pos (i, site)       &
         - float (cell (i) - 1) - cr_dim0 (i, 1)                        
         ENDDO 
         dx = skalpro (u, jdir, cr_gten) / rdj (ic) 
!                                                                       
         mc_energy_dis = mc_energy_dis + mo_const (ic) * delta * dx 
!                                                                       
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION mc_energy_dis                    
!*****7*****************************************************************
      REAL function mc_energy_dis_mol (nmol, imol, disp, rdi, rdj) 
!+                                                                      
!     Calculates the DELTA energy for distortions according to          
!                                                                       
!       DE = SUM k*disp*x                                               
!                                                                       
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxmol 
      PARAMETER (maxmol = chem_max_neig) 
!                                                                       
      INTEGER nmol, imol (2) 
      REAL disp (3, rmc_max_atom) 
      REAL rdi (chem_max_cor) 
      REAL rdj (chem_max_cor) 
!                                                                       
      INTEGER ineig (0:maxmol), nneig 
      INTEGER cell (3), site 
      INTEGER iatom, jatom 
      INTEGER i, ic, in, ia 
      REAL idir (3), jdir (3) 
      REAL delta, dx, u (3), v (3) 
!                                                                       
      REAL skalpro 
!                                                                       
      mc_energy_dis_mol = 0.0 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, nmol 
      iatom = mole_cont (mole_off (imol (ia) ) + 1) 
!                                                                       
      CALL indextocell (iatom, cell, site) 
      DO i = 1, 3 
      v (i) = cr_pos (i, iatom) - chem_ave_pos (i, site) - float (cell (&
      i) - 1) - cr_dim0 (i, 1)                                          
      v (i) = v (i) + disp (i, ia) 
      ENDDO 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
      IF (chem_ldall (ic) ) then 
         DO i = 1, 3 
         jdir (i) = v (i) 
         ENDDO 
         rdj (ic) = skalpro (jdir, jdir, cr_gten) 
         IF (rdj (ic) .gt.0.0) then 
            rdj (ic) = sqrt (rdj (ic) ) 
         ELSE 
            rdj (ic) = 1.0 
         ENDIF 
         delta = 1.0 
      ELSE 
         DO i = 1, 3 
         idir (i) = chem_dir (i, 1, ic) 
         jdir (i) = chem_dir (i, 2, ic) 
         ENDDO 
         delta = skalpro (v, idir, cr_gten) / rdi (ic) 
      ENDIF 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
      CALL chem_neighbour_mol (imol (ia), ic, ineig, nneig, maxmol) 
      IF (nneig.ne.0) then 
         DO in = 1, nneig 
         jatom = mole_cont (mole_off (ineig (in) ) + 1) 
         CALL indextocell (jatom, cell, site) 
         DO i = 1, 3 
         u (i) = cr_pos (i, iatom) - chem_ave_pos (i, site) - float (   &
         cell (i) - 1) - cr_dim0 (i, 1)                                 
         ENDDO 
         dx = skalpro (u, jdir, cr_gten) / rdj (ic) 
!                                                                       
         mc_energy_dis_mol = mc_energy_dis_mol + mo_const (ic) * delta *&
         dx                                                             
!                                                                       
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION mc_energy_dis_mol                
!*****7*****************************************************************
      REAL function mc_energy_spr (natoms, isel) 
!+                                                                      
!     Calculates the energy for distortions according to                
!                                                                       
!       E = SUM k*(d-d0)**2                                             
!                                                                       
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom), natoms 
!                                                                       
      INTEGER iatom (0:maxatom), natom 
      INTEGER cell (3), site 
      INTEGER i, is, js, ic, in, ia 
      REAL patom (3, 0:maxatom) 
      REAL d, u (3), v (3) 
!                                                                       
      REAL do_blen 
!                                                                       
      mc_energy_spr = 0.0 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, natoms 
!                                                                       
!------ - Restoring force to lattice (displ. from average position)     
!                                                                       
      CALL indextocell (isel (ia), cell, site) 
      DO i = 1, 3 
      u (i) = cr_pos (i, isel (ia) ) - chem_ave_pos (i, site) - float ( &
      cell (i) - 1) - cr_dim0 (i, 1)                                    
      ENDDO 
      d = do_blen (.true., u, u) 
!                                                                       
      mc_energy_spr = mc_energy_spr + mo_const (0) * d**2 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
      IF (chem_ctyp (ic) .eq.CHEM_VEC) then 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
         CALL chem_neighbour (isel (ia), ic, iatom, patom, natom,       &
         maxatom)                                                       
         IF (natom.ne.0) then 
            DO in = 1, natom 
            is = cr_iscat (isel (ia) ) 
            js = cr_iscat (iatom (in) ) 
            IF (mo_disp (ic, is, js) .ne.0.0) then 
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               v (i) = patom (i, in) 
               ENDDO 
               d = do_blen (.true., u, v) 
!                                                                       
               mc_energy_spr = mc_energy_spr + mo_const (ic) * (d-      &
               mo_disp (ic, is, js) ) **2                               
!                                                                       
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION mc_energy_spr                    
!*****7*****************************************************************
      REAL function mc_energy_spr_mol (nmol, imol) 
!+                                                                      
!     Calculates the energy for distortions for molecules ..            
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxmol 
      PARAMETER (maxmol = chem_max_neig) 
!                                                                       
      INTEGER imol (2), nmol 
!                                                                       
      INTEGER ineig (0:maxmol), nneig 
      INTEGER iatom, jatom 
      INTEGER cell (3), site 
      INTEGER i, is, js, ic, in, ia 
      REAL d, u (3), v (3) 
!                                                                       
      REAL do_blen 
!                                                                       
      mc_energy_spr_mol = 0.0 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, nmol 
      iatom = mole_cont (mole_off (imol (ia) ) + 1) 
!                                                                       
!------ - Restoring force to lattice (displ. from average position)     
!                                                                       
      CALL indextocell (iatom, cell, site) 
      DO i = 1, 3 
      u (i) = cr_pos (i, iatom) - chem_ave_pos (i, site) - float (cell (&
      i) - 1) - cr_dim0 (i, 1)                                          
      ENDDO 
      d = do_blen (.true., u, u) 
!                                                                       
      mc_energy_spr_mol = mc_energy_spr_mol + mo_const (0) * d**2 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
      IF (chem_ctyp (ic) .eq.CHEM_VEC) then 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
         CALL chem_neighbour_mol (imol (ia), ic, ineig, nneig, maxmol) 
         IF (nneig.ne.0) then 
            DO in = 1, nneig 
            is = mole_type (imol (ia) ) 
            js = mole_type (ineig (in) ) 
            jatom = mole_cont (mole_off (ineig (in) ) + 1) 
            IF (mo_disp (ic, is, js) .ne.0.0) then 
               DO i = 1, 3 
               u (i) = cr_pos (i, iatom) 
               v (i) = cr_pos (i, jatom) 
               ENDDO 
               d = do_blen (.true., u, v) 
!                                                                       
               mc_energy_spr_mol = mc_energy_spr_mol + mo_const (ic)    &
               * (d-mo_disp (ic, is, js) ) **2                          
!                                                                       
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION mc_energy_spr_mol                
!*****7*****************************************************************
      REAL function mc_energy_angle (natoms, isel) 
!+                                                                      
!     Calculates the energy for angular deviations according to         
!                                                                       
!       E = SUM k*(a-a0)**2                                             
!                                                                       
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE modify_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (CHEM_MAX_ATOM), natoms 
!                                                                       
      INTEGER iatom (0:maxatom), natom 
      INTEGER cell (3), site 
      INTEGER i, is, js, ic, in, ia, io 
      INTEGER ii 
      LOGICAL lnoneig 
      REAL patom (3, 0:maxatom) 
      REAL a, u (3), v (3), w (3) 
!                                                                       
      REAL do_bang 
!                                                                       
      mc_energy_angle = 0.0 
      lnoneig = .true. 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, natoms 
!                                                                       
!------ - Position of zentral atom                                      
!                                                                       
      CALL indextocell (isel (ia), cell, site) 
!                                                                       
!       mc_energy_angle = mc_energy_angle + mo_const(0)*d**2            
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
      IF (chem_ctyp (ic) .eq.CHEM_ANG) then 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
         CALL chem_neighbour (isel (ia), ic, iatom, patom, natom,       &
         maxatom)                                                       
         IF (natom.ne.0) then 
            DO ii = 1, natom, 3 
            DO i = 1, 3 
            u (i) = patom (i, ii) 
            ENDDO 
            in = ii + 1 
            is = cr_iscat (iatom (in) ) 
            js = cr_iscat (iatom (in + 1) ) 
            IF (mo_disp (ic, is, js) .ne.0.0) then 
               DO i = 1, 3 
               v (i) = patom (i, in) 
               w (i) = patom (i, in + 1) 
               ENDDO 
               a = do_bang (.true., v, u, w) 
!                                                                       
               mc_energy_angle = mc_energy_angle+mo_const (ic) *        &
               (a - mo_disp (ic, is, js) ) **2                          
!                                                                       
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
!DBG                                                                    
      IF (lnoneig) then 
         WRITE ( * , * ) ' No neighbour found ', natom 
      WRITE ( * ,  * ) ' Zentral Atom  : ', isel (ia) , cr_iscat (ia) , &
     &u                                                                 
      ENDIF 
!                                                                       
      END FUNCTION mc_energy_angle                  
!*****7*****************************************************************
      SUBROUTINE errlist_mc 
!-                                                                      
!     Displays error Messages for the error type MC                     
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = - 5, IO = 0) 
!                                                                       
      CHARACTER(41) ERROR (IU:IO) 
!                                                                       
      DATA ERROR / 'Number of feedback intervalls is zero      ', 'Numbe&
     &r of MC cycles is zero                ', 'Invalid mode selected fo&
     &r COCC MC run      ', 'No valid move after 1000 cycles            &
     &', 'Invalid or no energy type selected         ', ' ' /           
                                                            ! -5        
                                                            ! -4        
                                                            ! -3        
                                                            ! -2        
                                                            ! -1        
                  !  0                                                  
!                                                                       
      CALL disp_error ('MC  ', error, iu, io) 
      END SUBROUTINE errlist_mc                     
