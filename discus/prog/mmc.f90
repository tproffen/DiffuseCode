MODULE mmc_menu
!
CONTAINS
!*****7**************************************************************** 
!                                                                       
SUBROUTINE mmc 
!-                                                                      
!     This sublevel includes all commands and functions for the         
!     Monte-Carlo simulations in DISCUS.                                
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE discus_allocate_appl_mod
USE chem_mod
USE discus_init_mod
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE chem_symm_mod
!
USE doact_mod 
USE errlist_mod 
USE learn_mod 
USE class_macro_internal
USE param_mod 
USE prompt_mod 
IMPLICIT none 
!                                                                       
       
!                                                                       
INTEGER, PARAMETER :: MIN_PARA =  20 ! A command requires at least these no of parameters
INTEGER            :: maxw           ! Array size for cpara, lpara, werte
!                                                                       
CHARACTER(LEN=5)    :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=40)   :: cdummy 
CHARACTER(LEN=1024) :: line, zeile 
INTEGER             :: lp, length 
INTEGER             :: indxg, lbef
LOGICAL, PARAMETER  :: lold = .false. 
!
INTEGER             :: n_corr ! dummy for allocation
INTEGER             :: n_scat ! dummy for allocation
!                                                                       
INTEGER, EXTERNAL   :: len_str 
LOGICAL, EXTERNAL   :: str_comp 
!                                                                       
maxw = MAX(MIN_PARA,MAXSCAT+1)
!
! Basic allocation
!
n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
! call alloc_chem ! NEEDS WORK
call alloc_mmc ( n_corr, MC_N_ENERGY, n_scat )
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/mmc' 
!
   10 CONTINUE 
!                                                                       
      CALL no_error 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
             line == char(13) .or. line(1:1) == '!'  ) GOTO 10
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0.and.                                      &
          .not. (str_comp (befehl, 'echo', 2, lbef, 4) ) .and. &
          .not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and. &
          .not. (str_comp (befehl, 'help', 2, lbef, 4) .or.    &
                 str_comp (befehl, '?   ', 2, lbef, 4) ) ) then                                              
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
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or. &
              str_comp (befehl, '?   ', 1, lbef, 4) ) then                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 11 
               CALL do_hel ('discus mmc '//zeile, lp) 
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
!------ command 'rese'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'rese', 2, lbef, 3) ) then 
            CALL mmc_init 
!                                                                       
!------ command 'run'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
            CALL mmc_run_multi 
!                                                                       
!------ command 'save'                                                  
!                                                                       
!        ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
!           CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!           CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
!           IF (ier_num.eq.0) then 
!              WRITE (output_io, 1500) cpara (1)(1:len_str (cpara (1) ))
!              CALL save_struc (cpara (1), lpara (1) ) 
!           ENDIF 
!                                                                       
!------ command 'sel' selecting/deselecting atoms                       
!                                                                       
         ELSEIF (str_comp (befehl, 'sele', 3, lbef, 4) .or.   &
                 str_comp (befehl, 'dese', 2, lbef, 4) ) then
!                                                                       
            CALL atom_select (zeile, lp, 0, MMC_MAX_SCAT, mmc_latom, &
            mmc_sel_atom, lold, str_comp (befehl,  &
            'sele', 3, lbef, 4) )                                       
!                                                                       
!------ command 'set'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) then 
            CALL mmc_set (zeile, lp) 
!                                                                       
!------ command 'show'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            CALL mmc_show 
!
!------ command 'symmetry'
!
         ELSEIF (str_comp (befehl, 'apply_symm', 2, lbef, 10) ) then 
            CALL chem_symm(zeile, lp)
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
               IF(sprompt /= prompt) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in mmc menu'
                  prompt_status = PROMPT_ON 
                  RETURN
               ELSE
                  CALL macro_close 
                  prompt_status = PROMPT_ON 
               ENDIF 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               prompt_status = PROMPT_ON 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
      GOTO 10 
!                                                                       
 9999 CONTINUE 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE mmc                            
!*****7*****************************************************************
      SUBROUTINE mmc_show 
!+                                                                      
!     Show parameters of MMC section                                    
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE chem_mod 
      USE chem_menu
      USE mc_mod 
      USE mmc_mod 
      USE molecule_mod 
      USE rmc_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(1024) zeile 
!      CHARACTER(LEN=24), DIMENSION(0:MC_N_ENERGY) :: c_energy = & !(0:MC_N_ENERGY) 
!         (/'                        ', &
!           'Occupation correlation  ', &
!           'Displacement correlation', &
!           'Displacement energy     ', &
!           'Angular energy          ', &
!           'Vector energy           ', &
!           'Bond length energy      ', &
!           'Lennard Jones Potential ', &
!           'Buckingham    Potential ', &
!           'Repulsive  Potential    '/)
      CHARACTER (LEN=20), DIMENSION(MC_N_MOVE)    :: c_move = & !  (MC_N_MOVE) 
         (/ 'switch chemistry    ', &
            'switch displacement ', &
            'shift atom          ', &
            'inverse displacement' /)
      CHARACTER (LEN=20), DIMENSION(4)            :: c_site = & !(4) 
         (/ 'all                 ', &
            'local +- 1 unit cell', &
            'local, same site    ', &
            'all,   same site    ' /)
      CHARACTER(9) at_name_i, at_name_j 
      INTEGER i, j, k
!                                                                       
!                                                                       
      IF (mo_sel_atom) then 
         WRITE (output_io, 1105) 'atoms' 
      ELSE 
         WRITE (output_io, 1105) 'molecules' 
      ENDIF 
!                                                                       
      WRITE (output_io, 1250) mo_cyc 
      WRITE (output_io, 1300) mo_feed 
      WRITE (output_io, 1400) mo_kt 
!                                                                       
!     Information about the moves                                       
!                                                                       
      WRITE (output_io, 5000) 
      DO i = 1, MC_N_MOVE 
         WRITE (output_io, 5100) c_move (i), mmc_move_prob (i), &
                                 c_site (mmc_local (i) )
      ENDDO 
!                                                                       
!------ Information about defined correlations etc ..                   
!                                                                       
      WRITE (output_io, 3000) 
      zeile = 'corr' 
      CALL chem_show (zeile) 
!                                                                       
!------ Information about defined potential energies                    
!                                                                       
      IF (mmc_cor_energy (0, MC_OCC) ) then 
         WRITE (output_io, 2100) 
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_OCC, i, j) .ne.0.0) then 
               WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_OCC, i, j),                       &
               mmc_depth (k, MC_OCC,i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_OCC, i, j) .ne.0.0) then 
               WRITE (output_io, 4200) i, j, k,         &
                     mmc_target_corr (k,MC_OCC, i, j),  &
                     mmc_depth (k, MC_OCC, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (mmc_cor_energy (0, MC_DISP) ) then 
         WRITE (output_io, 2200) 
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_DISP, i, j) .ne.0.0) then 
               WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_DISP, i, j),                      &
               mmc_depth (k,MC_DISP, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_DISP, i, j) .ne.0.0) then 
               WRITE (output_io, 4200) i, j, k,                   &
               mmc_target_corr (k,MC_DISP, i, j),                 &
               mmc_depth (k, MC_DISP, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (mmc_cor_energy (0, MC_SPRING) ) then 
         WRITE (output_io, 2300) 
!                                                                       
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_SPRING, i, j) .ne.0.0) then 
               WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_SPRING, i, j),                    &
               mmc_depth (k,MC_SPRING, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_SPRING, i, j) .ne.0.0) then 
               WRITE (output_io, 4200) i, j, k,                         &
               mmc_target_corr (k,MC_SPRING, i, j),                     &
               mmc_depth (k, MC_SPRING, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Information about defined desired angles                        
!                                                                       
      IF (mmc_cor_energy (0, MC_ANGLE) ) then 
         WRITE (output_io, 2010) 
!                                                                       
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_ANGLE, i, j) .ne.0.0) then 
               WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_ANGLE, i, j), mmc_depth (k,       &
               MC_ANGLE, i, j)                                          
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_ANGLE, i, j) .ne.0.0) then 
               WRITE (output_io, 4200) i, j, k, mmc_target_corr (k,     &
               MC_ANGLE, i, j), mmc_depth (k, MC_ANGLE, i, j)           
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Information about defined desired vector differences            
!                                                                       
!     IF (mmc_cor_energy (0, MC_VECTOR) ) then 
!        WRITE (output_io, 2020) 
!                                                                       
!        IF (mo_sel_atom) then 
!           DO k = 1, chem_ncor 
!           DO i = 0, cr_nscat 
!           DO j = i, cr_nscat 
!           at_name_i = at_name (i) 
!           at_name_j = at_name (j) 
!           DO l = 1, mmc_nvec (k, i, j) 
!           WRITE (output_io, 4120) at_name_i, at_name_j, k,  &
!                  (mmc_vec ( m, l, k, i, j), m = 1, 4)
!           ENDDO 
!           ENDDO 
!           ENDDO 
!           ENDDO 
!        ELSE 
!           DO k = 1, chem_ncor 
!           DO i = 1, mole_num_type 
!           DO j = i, mole_num_type 
!           DO l = 1, mmc_nvec (k, i, j) 
!           WRITE (output_io, 4220) i, j, k,               &
!                  (mmc_vec (m, l, k, i, j), m = 1, 4)
!           ENDDO 
!           ENDDO 
!           ENDDO 
!           ENDDO 
!        ENDIF 
!     ENDIF 
!                                                                       
!------ Information about defined desired Lennard Jones Potentiale      
!                                                                       
      IF (mmc_cor_energy (0, MC_LENNARD) ) then 
         WRITE (output_io, 2700) 
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_LENNARD, i, j) .ne.0.0) then 
               WRITE (output_io, 7700) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_LENNARD, i, j),                   &
               mmc_depth (k, MC_LENNARD, i, j),                         &
               mmc_len_a (k, i, j), mmc_len_b (k, i, j),                &
               mmc_len_m (k, i, j), mmc_len_n (k, i, j)             
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_LENNARD, i, j) .ne.0.0) then 
               WRITE (output_io, 7700) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_LENNARD, i, j),                   &
               mmc_depth (k,MC_LENNARD, i, j),                          &
               mmc_len_a (k, i, j), mmc_len_b (k, i,  j),               &
               mmc_len_m (k, i, j), mmc_len_n (k, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Information about defined desired Repulsive Potentiale      
!                                                                       
      IF (mmc_cor_energy (0, MC_REPULSIVE) ) then 
         WRITE (output_io, 2750) 
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_REPULSIVE, i, j) .ne.0.0) then 
               WRITE (output_io, 7750) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_REPULSIVE, i, j),                 &
               mmc_depth (k, MC_REPULSIVE, i, j),                       &
               mmc_rep_a (k,i,j), mmc_rep_b (k,i,j),                    &
               mmc_rep_c (k,i,j), mmc_rep_m (k, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_REPULSIVE, i, j) .ne.0.0) then 
               WRITE (output_io, 7750) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_REPULSIVE, i, j),                 &
               mmc_depth (k,MC_REPULSIVE, i, j),                        &
               mmc_rep_a (k,i,j), mmc_rep_b (k,i,j),                    &
               mmc_rep_c (k,i,j), mmc_rep_m (k, i, j)
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Information about defined desired Buckingham Potentiale         
!                                                                       
      IF (mmc_cor_energy (0, MC_BUCKING) ) then 
         WRITE (output_io, 2800) 
         IF (mo_sel_atom) then 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_BUCKING, i, j) .ne.0.0) then 
               WRITE (output_io, 7800) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_BUCKING, i, j), mmc_depth (k,     &
               MC_BUCKING, i, j), mmc_buck_a (k, i, j), mmc_buck_rho (k,&
               i, j), mmc_buck_b (k, i, j)                              
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ELSE 
            DO k = 1, chem_ncor 
            DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
            IF (mmc_target_corr (k, MC_BUCKING, i, j) .ne.0.0) then 
               WRITE (output_io, 7800) at_name_i, at_name_j, k,         &
               mmc_target_corr (k, MC_BUCKING, i, j), mmc_depth (k,     &
               MC_BUCKING, i, j), mmc_buck_a (k, i, j), mmc_buck_rho (k,&
               i, j), mmc_buck_b (k, i, j)                              
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Maxmove settings in 'shift' mode                                
!                                                                       
      IF (mmc_move_prob (3) .gt.0.0) then 
         WRITE (output_io, 6000) 
         IF (mo_sel_atom) then 
            DO i = 0, cr_nscat 
               at_name_i = at_name (i) 
               IF(mo_maxmove(4,i)==0.0) THEN
                  WRITE(output_io, 6010) at_name_i,(mo_maxmove(j,i), j = 1, 3)
               ELSE
                  WRITE(output_io, 6015) at_name_i,(mo_maxmove(j,i), j = 1, 4)
               ENDIF
            ENDDO 
         ELSE 
            DO i = 1, mole_num_type 
               IF(mo_maxmove(4,i)==0.0) THEN
                  WRITE (output_io, 6020) i, (mo_maxmove (j, i), j = 1, 3) 
               ELSE
                  WRITE (output_io, 6025) i, (mo_maxmove (j, i), j = 1, 3) 
               ENDIF
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
 1105 FORMAT (  '   Operation mode for MC        : ',a) 
 1250 FORMAT (  '   Max. number of MC cycles     : ',i8) 
 1300 FORMAT (  '   Feedback/update intervall    : ',i8) 
 1400 FORMAT (  '   Temperature [kT]             : ',f8.4) 
 3000 FORMAT (/,' Correlation definitions        : ',/) 
 2100 FORMAT (/,' Desired correlations for Chemical Occupancy : ',/,/,  &
     &         12x,'Pairs',10x,'neigh. #',3x,'correl. ',7x,'depth')     
 2200 FORMAT (/,' Desired correlations for Displacement : ',/,/,        &
     &         12x,'Pairs',10x,'neigh. #',3x,'correl. ',7x,'depth')     
 2300 FORMAT (/,' Desired distortions for SPRING : ',/,/,               &
     &         12x,'Pairs',10x,'neigh. #',3x,'distance',7x,'depth')     
 2010 FORMAT (/,' Desired distortions for ANGLE  : ',/) 
 2700 FORMAT (/,' Desired distance for LENNARD : ',/,/,                 &
     &         12x,'Pairs',10x,'neigh. #',3x,'distance/A',7x,           &
     &             'depth/B',7x,'m',7x,'n')                             
 2750 FORMAT (/,' Desired distance for REPULSIVE : ',/,/,               &
     &         12x,'Pairs',10x,'neigh. #',3x,'distance/A',7x,           &
     &             'depth/B',7x,'short/C',7x,'m'       )                             
 2800 FORMAT (/,' Desired distortions for BUCKING : ',/,/,              &
     &         12x,'Pairs',10x,'neigh. #',3x,'distance/A',7x,           &
     &             'depth/Rho',6x,'B')                                  
 7300 FORMAT (  5x,a8,' - ',a8,5x,i3,5x,f7.3,5x,g18.8e2) 
 7700 FORMAT (  5x,a8,' - ',a8,5x,i3,5x,f7.3,5x,g18.8e2,/,              &
     &         30x,g18.8e2,g18.8e2,f7.3,2x,f7.3)                        
 7750 FORMAT (  5x,a8,' - ',a8,5x,i3,5x,f7.3,5x,g18.8e2,/,              &
     &         30x,g18.8e2,g18.8e2,g18.8e2,f7.3,2x     )                        
 7800 FORMAT (  5x,a8,' - ',a8,5x,i3,5x,f7.3,5x,g18.8e2,/,              &
     &         30x,g18.8e2,g18.8e2,g18.8e3)                             
 4200 FORMAT (  '   Molecule types : ',i4,' - ',i4,'  neig. #',i3,      &
     &          '   distance : ',f7.3,' A',' depth ',g18.8e2)           
 6000 FORMAT (/,' Sigmas for MC shifts (l.u.)    : ') 
 6010 FORMAT (  '                 Atom ',a9,' : ',3(F9.5,1X)) 
 6015 FORMAT (  '   Atom Vector; shift ',a9,' : ',3(F9.5,1X),'; ',F9.5) 
 6020 FORMAT (  '        Molecule type ',i9,' : ',3(F9.5,1X)) 
 6025 FORMAT (  '   Mol type; Vec; shft',i9,' : ',3(F9.5,1X),'; ',F9.5) 
 5000 FORMAT (/' Operation modes for MMC '/                             &
     &        '   Mode                  Probability',                   &
     &        '  two-atom correlation'/3x,55('-'))                      
 5100 FORMAT (3x,a20,2x,f7.5,6x,a20) 
      END SUBROUTINE mmc_show                       
!*****7*****************************************************************
      SUBROUTINE mmc_set (zeile, lp) 
!+                                                                      
!     sets parameters for MMC section                                   
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE chem_menu
      USE mc_mod 
      USE mmc_mod 
      USE modify_mod
      USE rmc_sup_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
!                                                                       
      CHARACTER ( LEN=* ), INTENT(INOUT)    :: zeile 
      INTEGER            , INTENT(INOUT)    :: lp 
!
      CHARACTER ( LEN=1024)                 :: line
      INTEGER                               :: length
      INTEGER                               :: ianz1
      INTEGER                               :: ianz2
!                                                                       
      CHARACTER (LEN=1024), DIMENSION(MAXW) ::  cpara !(maxw) 
      CHARACTER (LEN=1024), DIMENSION(MAXW) ::  cpara1 !(maxw) 
      CHARACTER (LEN=1024), DIMENSION(MAXW) ::  cpara2 !(maxw) 
      REAL uerte (maxw) 
      REAL verte (maxw) 
      REAL werte (maxw) 
      REAL   , DIMENSION(MAXW) :: werte1 (maxw) 
      REAL   , DIMENSION(MAXW) :: werte2 (maxw) 
      REAL a, b 
      INTEGER lpara (maxw) 
      INTEGER, DIMENSION(MAXW) :: lpara1 ! (maxw) 
      INTEGER, DIMENSION(MAXW) :: lpara2 ! (maxw) 
      INTEGER ianz, iianz, jjanz, kkanz, is, js, ls, ic, i, j 
      INTEGER is_start, is_end 
      INTEGER                :: n_corr ! Dummy for allocation
      INTEGER                :: n_scat ! Dummy for allocation
      INTEGER                :: n_angles ! Dummy for allocation
!                                                                       
!     INTEGER angles2index 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      n_corr = 0
      n_scat = 0
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.2) then 
            CALL do_cap (cpara (1) ) 
!
!------ --- 'set allowed', set which atom are allowed in a move
!
            IF (cpara (1) (1:2) .eq.'AL') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL get_iscat (ianz, cpara, lpara, werte,     &
                        maxw, .false.)                                  
               IF ( werte(1).eq.-1) THEN
                  mmc_allowed = .true.
               ELSE
                  DO i=1,ianz
                     mmc_allowed(nint(werte(i))) = .true.
                  ENDDO
               ENDIF
!
!------ --- 'set disallowed', set which atom are allowed in a move
!
            ELSEIF (cpara (1) (1:2) .eq.'DI') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL get_iscat (ianz, cpara, lpara, werte,     &
                        maxw, .false.)                                  
               IF ( werte(1).eq.-1) THEN
                  mmc_allowed = .true.
               ELSE
                  DO i=1,ianz
                     mmc_allowed(nint(werte(i))) = .false.
                  ENDDO
               ENDIF
!                                                                       
!------ --- 'set angle' : setting of neighbouring angles                
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'AN') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_angle (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set conn' : setting of neighbouring connectivity           
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'CON') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_con (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set envir' : setting of neighbouring environment           
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'ENV') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_envir (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set cyc' : setting number of MC moves                      
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'CY') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_cyc = nint (werte (1) ) 
!                                                                       
!------ --- 'set feed' : setting display/feedback intervall             
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'FE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_feed = nint (werte (1) ) 
!                                                                       
!------ --- 'set limited': sets limited selction range for atoms        
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'LIM') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (cpara (1) (1:4) .eq.'cell') then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) return 
                  DO i = 1, 3 
                  mmc_l_center (i) = int( werte (i) )
                  mmc_l_extend (i) = int( werte (i + 3) )
                  ENDDO 
                  mmc_l_limited = .true. 
                  mmc_l_type = MMC_L_CELLS 
               ELSEIF (cpara (1) (1:4) .eq.'atom') then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) return 
                  mmc_l_lower = int( werte (1) )
                  mmc_l_upper = int( werte (2) )
                  mmc_l_limited = .true. 
                  mmc_l_type = MMC_L_ATOMS 
               ELSEIF (cpara (1) (1:3) .eq.'OFF') then 
                  mmc_l_limited = .false. 
               ENDIF 
!                                                                       
!------ --- 'set mode': sets operation mode for MMC                     
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'MOD') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               j = 1 
               CALL ber_params (j, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL mmc_set_mode (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set move': sets maxmove for shift MMC mode                 
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'MOV') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL rmc_set_move (mo_maxmove, mo_sel_atom, ianz, cpara, &
               werte, lpara, maxw, 4)                                      
!                                                                       
!------ --- 'set neig': setting correlation determination method        
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'NE') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_neig (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set range' : setting of neighbouring ranges                
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'RA') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_ranges (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set target' : setting of target correlations               
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'TA') then 
               IF (ianz.ge.3) then 
                  is = - 2 
                  js = - 2 
                  ls = - 2 
                  werte  = 0.0
                  uerte  = 0.0
                  verte  = 0.0
                  werte1 = 0.0
                  werte2 = 0.0
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (1, cpara, lpara, werte, maxw) 
                  ic = nint (werte (1) ) 
                  IF( ic > CHEM_MAX_COR .or. ic > MMC_MAX_CORR ) THEN
!
!                      Basic allocation
!
                     n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
                     n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
                     ! call alloc_chem ! NEEDS WORK
                     call alloc_mmc ( n_corr, MC_N_ENERGY, n_scat )
                  ENDIF
                  IF (ic.gt.0.and.ic.le.chem_ncor) then 
                     IF (str_comp (cpara (2) , 'corr', 2, lpara (2) , 4)) then                                             
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
!                       Get atom types in the two allowed groups
!                       Allowed parameters are: Au, Cu            ! OLd style 
!                                               (Au), (Cu)        ! New style  :binary correlation
!                                               (Au,Pt), (Cu,Zn)  ! New style  :quaternary correlation
                        IF ( cpara(1)(1:1)=='(' .and. cpara(1)(lpara(1):lpara(1))==')') THEN
                           line   = cpara(1)(2:lpara(1)-1)
                           length = lpara(1)-2
                           CALL get_params (line, ianz1, cpara1, lpara1, maxw, length) 
                           CALL get_iscat (ianz1, cpara1, lpara1, werte1, maxw, .false.)                                  
                        ELSEIF ( cpara(1)(1:1)/='(' .and. cpara(1)(lpara(1):lpara(1))/=')') THEN
                           ianz1 = 1 
                           CALL get_iscat (ianz1, cpara, lpara, werte1, maxw, .false.)                                  
                        ELSE
                           ier_num = -6
                           ier_typ = ER_COMM
                           RETURN
                        ENDIF
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF ( cpara(1)(1:1)=='(' .and. cpara(1)(lpara(1):lpara(1))==')') THEN
                           line   = cpara(1)(2:lpara(1)-1)
                           length = lpara(1)-2
                           CALL get_params (line, ianz2, cpara2, lpara2, maxw, length) 
                           CALL get_iscat (ianz2, cpara2, lpara2, werte2, maxw, .false.)                                  
                        ELSEIF ( cpara(1)(1:1)/='(' .and. cpara(1)(lpara(1):lpara(1))/=')') THEN
                           ianz2 = 1 
                           CALL get_iscat (ianz2, cpara, lpara, werte2, maxw, .false.)                                  
                        ELSE
                           ier_num = -6
                           ier_typ = ER_COMM
                           RETURN
                        ENDIF
!
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF (cpara (ianz) (1:2) .eq.'CO') then 
                           mmc_cfac (ic, MC_OCC) = 1.0 
                           ianz = ianz - 1 
                        ELSEIF (cpara (ianz) (1:2) .eq.'EN') then 
                           mmc_cfac (ic, MC_OCC) = 0.0 
                           ianz = ianz - 1 
                        ENDIF 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)
                        IF (mmc_cfac (ic, MC_OCC) ==1.0) then 
                        CALL mmc_set_disp_occ (ic, MC_OCC, ianz1, ianz2, &
                             MAXW, werte1, werte2, werte(1) , -werte(1) )                                   
                           mmc_depth (ic, MC_OCC, 0, 0) = -werte (1) 
                        ELSEIF (mmc_cfac (ic, MC_OCC) ==0.0) then 
                        CALL mmc_set_disp_occ (ic, MC_OCC, ianz1, ianz2, &
                             MAXW, werte1, werte2, werte(1) , werte(2) )                                   
                           mmc_depth (ic, MC_OCC, 0, 0) = werte (2) 
                        ENDIF
                        mmc_cor_energy (ic, MC_OCC) = .true. 
                        mmc_cor_energy (0, MC_OCC) = .true. 
                     ELSEIF (str_comp (cpara (2) , 'cd', 2, lpara (2) , &
                     2) ) then                                          
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (2, cpara, lpara, werte, maxw) 
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           mmc_allowed(is) = .true. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .true. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_DISP, is, js, werte ( &
                        1), werte (2) )                                 
                        ENDDO 
                        ENDDO 
                        chem_ldall (ic) = .true. 
                        IF (ianz.gt.2) then 
                           CALL del_params (2, ianz, cpara, lpara, maxw) 
                           IF (cpara (1) (1:3) .eq.'all') then 
                              chem_ldall (ic) = .true. 
                           ELSE 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.ne.0) return 
!                                                                       
                              IF (ianz.eq.3) then 
                                 chem_ldall (ic) = .false. 
                                 chem_dir (1, 1, ic) = werte (1) 
                                 chem_dir (2, 1, ic) = werte (2) 
                                 chem_dir (3, 1, ic) = werte (3) 
                                 chem_dir (1, 2, ic) = werte (1) 
                                 chem_dir (2, 2, ic) = werte (2) 
                                 chem_dir (3, 2, ic) = werte (3) 
                              ELSEIF (ianz.eq.6) then 
                                 chem_ldall (ic) = .false. 
                                 chem_dir (1, 1, ic) = werte (1) 
                                 chem_dir (2, 1, ic) = werte (2) 
                                 chem_dir (3, 1, ic) = werte (3) 
                                 chem_dir (1, 2, ic) = werte (4) 
                                 chem_dir (2, 2, ic) = werte (5) 
                                 chem_dir (3, 2, ic) = werte (6) 
                              ELSE 
                                 ier_num = - 6 
                                 ier_typ = ER_COMM 
                              ENDIF 
                           ENDIF 
                        ENDIF 
                        mmc_cor_energy (ic, MC_DISP) = .true. 
                        mmc_cor_energy (0, MC_DISP) = .true. 
                     ELSEIF (str_comp (cpara (2) , 'spring', 2, lpara ( &
                     2) , 6) ) then                                     
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (2, cpara, lpara, werte, maxw) 
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           mmc_allowed(is) = .true. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .true. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_SPRING, is, js, werte &
                        (1), werte (2) )                                
                        ENDDO 
                        ENDDO 
                        mmc_cor_energy (ic, MC_SPRING) = .true. 
                        mmc_cor_energy (0, MC_SPRING) = .true. 
                     ELSEIF (str_comp (cpara (2) , 'angle', 2, lpara (2)&
                     , 5) ) then                                        
!                                                                       
!     ------------Three atom names are listed on the target instruction 
!                 Required for angular correlations                     
!                                                                       
                        IF (mmc_n_angles == 0          .or.     &
                            mmc_n_angles >= MMC_MAX_ANGLES) then 
                           n_angles = max(mmc_n_angles+20, int(MMC_MAX_ANGLES*1.025))
                           CALL alloc_mmc_angle (CHEM_MAX_COR,n_angles)
                        ENDIF
                        IF (mmc_n_angles.lt.MMC_MAX_ANGLES) then 
                           CALL del_params (2, ianz, cpara, lpara, maxw) 
                           mmc_n_angles = mmc_n_angles + 1 
                           iianz = 1 
                           jjanz = 1 
                           kkanz = 1 
                           CALL get_iscat (iianz, cpara, lpara, uerte,  &
                           maxw, .false.)                               
                           IF (uerte (1) .eq. - 1) then 
                              is_start = 0 
                              is_end = cr_nscat 
                              mmc_allowed = .true.           ! all atoms are allowed in mmc moves
                           ELSE 
                              is_start = int( uerte (1) )
                              is_end = int( uerte (1) )
                              mmc_allowed(is_start) = .true. ! this atom is allowed in mmc moves
                           ENDIF 
                           is = int( uerte (1) )
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL get_iscat (jjanz, cpara, lpara, uerte,  &
                           maxw, .false.)                               
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL get_iscat (kkanz, cpara, lpara, verte,  &
                           maxw, .false.)                               
                           js = min1 (uerte (1), verte (1) ) 
                           ls = max1 (uerte (1), verte (1) ) 
                           mmc_allowed(js) = .true. ! this atom is allowed in mmc moves
                           mmc_allowed(ls) = .true. ! this atom is allowed in mmc moves
                           i = angles2index (ic, mmc_n_angles, is, js,  &
                           ls, MAXSCAT)                 
                           mmc_angles (mmc_n_angles) = i 
                           CALL index2angles (i, ic, mmc_n_angles, is,  &
                           js, ls, MAXSCAT)             
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (2, cpara, lpara, werte,     &
                           maxw)                                        
                           mmc_target_angl (mmc_n_angles) = werte (1) 
                           mmc_target_corr (ic, MC_ANGLE, js, ls)       &
                           = werte (1)                                  
                           mmc_target_corr (ic, MC_ANGLE, ls, js)       &
                           = werte (1)                                  
                           mmc_depth_angl (mmc_n_angles) = werte (2) 
                           mmc_depth (ic, MC_ANGLE, js, ls) = werte (2) 
                           mmc_depth (ic, MC_ANGLE, ls, js) = werte (2) 
                           mmc_cor_energy (ic, MC_ANGLE) = .true. 
                           mmc_cor_energy (0, MC_ANGLE) = .true. 
                        ENDIF 
!DBG_VEC For later development                                          
!                    ELSEIF (str_comp (cpara (2) , 'vector', 2, lpara ( &
!                    2) , 6) ) then                                     
!                       CALL del_params (2, ianz, cpara, lpara, maxw) 
!                       iianz = 1 
!                       jjanz = 1 
!                       CALL get_iscat (iianz, cpara, lpara, uerte,     &
!                       maxw, .false.)                                  
!                       CALL del_params (1, ianz, cpara, lpara, maxw) 
!                       CALL get_iscat (jjanz, cpara, lpara, verte,     &
!                       maxw, .false.)                                  
!                       CALL del_params (1, ianz, cpara, lpara, maxw) 
!                       CALL ber_params (4, cpara, lpara, werte, maxw) 
!                       DO i = 1, iianz 
!                       DO j = 1, jjanz 
!                       is = nint (uerte (i) ) 
!                       js = nint (verte (j) ) 
!                       CALL mmc_set_vec (ic, is, js, werte, maxw) 
!                       ENDDO 
!                       ENDDO 
!                       mmc_cor_energy (ic, MC_VECTOR) = .true. 
!                       mmc_cor_energy (0, MC_VECTOR) = .true. 
                     ELSEIF (str_comp (cpara (2) , 'lennard', 2, lpara (&
                     2) , 6) ) then                                     
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ianz.eq.3) then 
                           werte (4) = 12.0 
                           werte (5) = 6.0 
                        ENDIF 
                        IF(ic > MMC_LENN_CORR .or.  & ! Allocate Lennard
                           ic > CHEM_MAX_COR  .or.  &
                           MAXSCAT > MMC_LENN_SCAT   ) THEN
                           n_corr = MAX(n_corr, CHEM_MAX_COR, &
                                                MMC_LENN_CORR)
                           n_scat = MAX(MAXSCAT,MMC_LENN_SCAT)
                           call alloc_mmc_lenn (n_corr, n_scat)
                           IF(ier_num /= 0) THEN
                              RETURN
                           ENDIF
                        ENDIF
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           mmc_allowed(is) = .true. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .true. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_LENNARD, is, js,      &
                        ABS(werte (1)), werte (2) )                          
                        a = - ABS(werte(2)) * werte(4) / (werte(4)       &
                        - werte (3) ) * werte (1) **werte (3)           
                        b = - ABS(werte(2)) * werte(3) / (werte(4)       &
                        - werte (3) ) * werte (1) **werte (4)           
                        CALL mmc_set_lenn (ic, is, js, a, b,&
                        werte (3), werte (4) )                          
                        ENDDO 
                        ENDDO 
                        mmc_cor_energy (ic, MC_LENNARD) = .true. 
                        mmc_cor_energy (0, MC_LENNARD) = .true. 
                     ELSEIF (str_comp (cpara (2) , 'repulsive', 2, &
                                       lpara (2) , 9)        ) then                                     
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF (ianz >  0) then 
                          CALL ber_params (ianz, cpara, lpara, werte,   &
                                           maxw)
                        ENDIF
                        IF (ianz == 0) then 
                           werte (1) =  0.0 
                           werte (2) =  1.0 
                           werte (3) =  0.0 
                           werte (4) =  1.0 
                        ELSEIF (ianz == 1) then 
                           werte (2) =  1.0 
                           werte (3) =  0.0 
                           werte (4) =  1.0 
                        ELSEIF (ianz == 2) then 
                           werte (3) =  0.0 
                           werte (4) =  1.0 
                        ELSEIF (ianz == 3) then 
                           werte (4) =  1.0 
                        ENDIF 
                        IF(ic > MMC_REP_CORR .or.  & ! Allocate Repulsive
                           ic > CHEM_MAX_COR  .or.  &
                           MAXSCAT > MMC_REP_SCAT   ) THEN
                           n_corr = MAX(n_corr, CHEM_MAX_COR, &
                                                MMC_REP_CORR)
                           n_scat = MAX(MAXSCAT,MMC_REP_SCAT)
                           call alloc_mmc_rep (n_corr, n_scat)
                           IF(ier_num /= 0) THEN
                              RETURN
                           ENDIF
                        ENDIF
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           mmc_allowed(is) = .true. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .true. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_REPULSIVE, is, js,    &
                        100.0    , ABS(werte (1)) )                          
                        CALL mmc_set_rep  (ic, is, js,    &
                        ABS(werte (1)), werte(2) , werte(3), werte(4) )                          
                        ENDDO 
                        ENDDO 
                        mmc_cor_energy (ic, MC_REPULSIVE) = .true. 
                        mmc_cor_energy (0, MC_REPULSIVE) = .true. 
                     ELSEIF (str_comp (cpara (2) , 'bucking', 2,        &
                                       lpara (2) , 6)           )THEN
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .false.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF(ic > MMC_BUCK_CORR .or.  & ! Allocate Buckingham
                           ic > CHEM_MAX_COR  .or.  &
                           MAXSCAT > MMC_BUCK_SCAT   ) THEN
                           n_corr = MAX(n_corr, CHEM_MAX_COR, &
                                                MMC_BUCK_CORR)
                           n_scat = MAX(MAXSCAT,MMC_BUCK_SCAT)
                           call alloc_mmc_buck (n_corr, n_scat)
                           IF(ier_num /= 0) THEN
                              RETURN
                           ENDIF
                        ENDIF
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           mmc_allowed(is) = .true. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .true. ! this atom is allowed in mmc moves
                        mmc_buck_a (ic, is, js) = werte (2) 
                        mmc_buck_rho (ic, is, js) = werte (3) 
                        mmc_buck_b (ic, is, js) = werte (4) 
                        mmc_buck_a (ic, js, is) = werte (2) 
                        mmc_buck_rho (ic, js, is) = werte (3) 
                        mmc_buck_b (ic, js, is) = werte (4) 
                        CALL find_bucking (werte, MAXW) 
                        mmc_buck_rmin (ic, js, is) = werte (7) 
                        mmc_buck_atmin (ic, js, is) = werte (8) 
                        mmc_buck_rmin (ic, is, js) = werte (7) 
                        mmc_buck_atmin (ic, is, js) = werte (8) 
                        CALL mmc_set_disp (ic, MC_BUCKING, is, js,      &
                        werte (5), abs (werte (6) ) )                   
                        ENDDO 
                        ENDDO 
                        mmc_cor_energy (ic, MC_BUCKING) = .true. 
                        mmc_cor_energy (0, MC_BUCKING) = .true. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
!                                                                       
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
!------ --- 'set fixed' : setting fixed atom coordinate ranges          
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'FI') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (str_comp (cpara (1) , 'OFF', 2, lpara (1) , 3) )     &
               then                                                     
                  mmc_l_constrains = .false. 
               ELSEIF (str_comp (cpara (1) , 'X', 1, lpara (1) , 1) )   &
               then                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_min (1) = werte (1) 
                  mmc_c_max (1) = werte (2) 
                  mmc_constrain_type = MMC_C_XYZ 
               ELSEIF (str_comp (cpara (1) , 'Y', 1, lpara (1) , 1) )   &
               then                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_min (2) = werte (1) 
                  mmc_c_max (2) = werte (2) 
                  mmc_constrain_type = MMC_C_XYZ 
               ELSEIF (str_comp (cpara (1) , 'Z', 1, lpara (1) , 1) )   &
               then                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_min (3) = werte (1) 
                  mmc_c_max (3) = werte (2) 
                  mmc_constrain_type = MMC_C_XYZ 
               ELSEIF (str_comp (cpara (1) , 'RAD', 1, lpara (1) , 3) ) &
               then                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_rad = werte (1) 
                  mmc_constrain_type = MMC_C_RADIUS 
               ENDIF 
!                                                                       
!------ --- no valid subcommand                                         
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
      END SUBROUTINE mmc_set                        
!*****7*****************************************************************
      SUBROUTINE mmc_set_disp (ic, ie, is, js, dist, depth) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE mmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER ic, ie, is, js, ii, jj 
      REAL dist 
      REAL depth 
!                                                                       
                                                                        
      IF (is.ne. - 1.and.js.ne. - 1) then 
         mmc_target_corr (ic, ie, is, js) = dist 
         mmc_target_corr (ic, ie, js, is) = dist 
         mmc_depth (ic, ie, is, js) = depth 
         mmc_depth (ic, ie, js, is) = depth 
         mmc_pair (ic, ie, is, js) = -1     
         mmc_pair (ic, ie, js, is) = -1     
      ELSEIF (is.eq. - 1.and.js.ne. - 1) then 
         DO ii = 0, cr_nscat 
         mmc_target_corr (ic, ie, ii, js) = dist 
         mmc_target_corr (ic, ie, js, ii) = dist 
         mmc_depth (ic, ie, ii, js) = depth 
         mmc_depth (ic, ie, js, ii) = depth 
         mmc_pair (ic, ie, ii, js) = -1     
         mmc_pair (ic, ie, js, ii) = -1     
         ENDDO 
      ELSEIF (is.ne. - 1.and.js.eq. - 1) then 
         DO ii = 0, cr_nscat 
         mmc_target_corr (ic, ie, ii, is) = dist 
         mmc_target_corr (ic, ie, is, ii) = dist 
         mmc_depth (ic, ie, ii, is) = depth 
         mmc_depth (ic, ie, is, ii) = depth 
         mmc_pair (ic, ie, ii, is) = -1     
         mmc_pair (ic, ie, is, ii) = -1     
         ENDDO 
      ELSE 
         DO ii = 0, cr_nscat 
         DO jj = 0, cr_nscat 
         mmc_target_corr (ic, ie, ii, jj) = dist 
         mmc_target_corr (ic, ie, jj, ii) = dist 
         mmc_depth (ic, ie, ii, jj) = depth 
         mmc_depth (ic, ie, jj, ii) = depth 
         mmc_pair (ic, ie, ii, jj) = -1     
         mmc_pair (ic, ie, jj, ii) = -1     
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE mmc_set_disp                   
!*****7*****************************************************************
      SUBROUTINE mmc_set_disp_occ (ic, ie, ianz1, ianz2, &
                             MAXW, werte1, werte2, corr, depth )
!
      USE discus_config_mod 
      USE crystal_mod 
      USE mmc_mod 
      IMPLICIT NONE
!
      INTEGER                 , INTENT(IN) :: ic     ! Correlation number
      INTEGER                 , INTENT(IN) :: ie     ! Energy number == MC_OCC
      INTEGER                 , INTENT(IN) :: ianz1  ! No of atom types in first group
      INTEGER                 , INTENT(IN) :: ianz2  ! No of atom types in second group
      INTEGER                 , INTENT(IN) :: MAXW   ! Array Dimension 
      REAL   , DIMENSION(MAXW), INTENT(IN) :: werte1 ! Actual atom types group1
      REAL   , DIMENSION(MAXW), INTENT(IN) :: werte2 ! Actual atom types group1
      REAL   ,                  INTENT(IN) :: corr   ! Desired correlation
      REAL   ,                  INTENT(IN) :: depth  ! Energy Depth
! 
      INTEGER                              :: is, js ! Dummy atom types
      INTEGER                              :: i, j   ! Loop indices
!
      DO i=1,ianz1              ! Set "equal" pairs first group
         is = NINT(werte1(i))
         mmc_allowed(is) = .true.
         DO j=1,ianz1
            js = NINT(werte1(j))
         mmc_target_corr (ic, ie, is, js) = corr 
         mmc_depth       (ic, ie, is, js) = depth 
         mmc_pair        (ic, ie, is, js) = +1     ! These pairs contribute positively to energy
         END DO
      END DO
      DO i=1,ianz2              ! Set "equal" pairs second group
         is = NINT(werte2(i))
         mmc_allowed(is) = .true.
         DO j=1,ianz1
            js = NINT(werte2(j))
         mmc_target_corr (ic, ie, is, js) = corr 
         mmc_depth       (ic, ie, is, js) = depth 
         mmc_pair        (ic, ie, is, js) = +2     ! These pairs contribute positively to energy
         END DO
      END DO
      DO i=1,ianz1              ! Set "opposite" pairs
         is = NINT(werte1(i))
         DO j=1,ianz2
            js = NINT(werte2(j))
            mmc_target_corr (ic, ie, is, js) = corr 
            mmc_target_corr (ic, ie, js, is) = corr 
            mmc_depth       (ic, ie, is, js) = depth 
            mmc_depth       (ic, ie, js, is) = depth 
            mmc_pair        (ic, ie, is, js) = -1     ! These pairs contribute negatively to energy
            mmc_pair        (ic, ie, js, is) = -2     ! Second ==> first group
         END DO
      END DO
!
      END SUBROUTINE mmc_set_disp_occ
!*****7*****************************************************************
      SUBROUTINE mmc_set_lenn (ic, is, js, a, b, m, n) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE mmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER ic, is, js, ii, jj 
      REAL a, b, m, n 
!                                                                       
                                                                        
      IF (is.ne. - 1.and.js.ne. - 1) then 
         mmc_len_a (ic, is, js) = a 
         mmc_len_b (ic, is, js) = b 
         mmc_len_m (ic, is, js) = m 
         mmc_len_n (ic, is, js) = n 
         mmc_len_a (ic, js, is) = a 
         mmc_len_b (ic, js, is) = b 
         mmc_len_m (ic, js, is) = m 
         mmc_len_n (ic, js, is) = n 
      ELSEIF (is.eq. - 1.and.js.ne. - 1) then 
         DO ii = 0, cr_nscat 
         mmc_len_a (ic, ii, js) = a 
         mmc_len_b (ic, ii, js) = b 
         mmc_len_m (ic, ii, js) = m 
         mmc_len_n (ic, ii, js) = n 
         mmc_len_a (ic, js, ii) = a 
         mmc_len_b (ic, js, ii) = b 
         mmc_len_m (ic, js, ii) = m 
         mmc_len_n (ic, js, ii) = n 
         ENDDO 
      ELSEIF (is.ne. - 1.and.js.eq. - 1) then 
         DO ii = 0, cr_nscat 
         mmc_len_a (ic, ii, is) = a 
         mmc_len_b (ic, ii, is) = b 
         mmc_len_m (ic, ii, is) = m 
         mmc_len_n (ic, ii, is) = n 
         mmc_len_a (ic, is, ii) = a 
         mmc_len_b (ic, is, ii) = b 
         mmc_len_m (ic, is, ii) = m 
         mmc_len_n (ic, is, ii) = n 
         ENDDO 
      ELSE 
         DO ii = 0, cr_nscat 
         DO jj = 0, cr_nscat 
         mmc_len_a (ic, ii, jj) = a 
         mmc_len_b (ic, ii, jj) = b 
         mmc_len_m (ic, ii, jj) = m 
         mmc_len_n (ic, ii, jj) = n 
         mmc_len_a (ic, jj, ii) = a 
         mmc_len_b (ic, jj, ii) = b 
         mmc_len_m (ic, jj, ii) = m 
         mmc_len_n (ic, jj, ii) = n 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE mmc_set_lenn                   
!*****7*****************************************************************
      SUBROUTINE mmc_set_rep (ic, is, js, a, b, c, m) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE mmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER ic, is, js, ii, jj 
      REAL a,  b, c, m
!                                                                       
!     Define mimimum energy at infinity distance
!
      mmc_rep_low = MIN(mmc_rep_low, ABS(a))
                                                                        
      IF (is.ne. - 1.and.js.ne. - 1) then 
         mmc_rep_a (ic, is, js) = a 
         mmc_rep_b (ic, is, js) = b 
         mmc_rep_c (ic, is, js) = c 
         mmc_rep_m (ic, is, js) = m 
         mmc_rep_a (ic, js, is) = a 
         mmc_rep_b (ic, js, is) = b 
         mmc_rep_c (ic, js, is) = c 
         mmc_rep_m (ic, js, is) = m 
      ELSEIF (is.eq. - 1.and.js.ne. - 1) then 
         DO ii = 0, cr_nscat 
         mmc_rep_a (ic, ii, js) = a 
         mmc_rep_b (ic, ii, js) = b 
         mmc_rep_c (ic, ii, js) = c 
         mmc_rep_m (ic, ii, js) = m 
         mmc_rep_a (ic, js, ii) = a 
         mmc_rep_b (ic, js, ii) = b 
         mmc_rep_c (ic, js, ii) = c 
         mmc_rep_m (ic, js, ii) = m 
         ENDDO 
      ELSEIF (is.ne. - 1.and.js.eq. - 1) then 
         DO ii = 0, cr_nscat 
         mmc_rep_a (ic, ii, is) = a 
         mmc_rep_b (ic, ii, is) = b 
         mmc_rep_c (ic, ii, is) = c 
         mmc_rep_m (ic, ii, is) = m 
         mmc_rep_a (ic, is, ii) = a 
         mmc_rep_b (ic, is, ii) = b 
         mmc_rep_c (ic, is, ii) = c 
         mmc_rep_m (ic, is, ii) = m 
         ENDDO 
      ELSE 
         DO ii = 0, cr_nscat 
         DO jj = 0, cr_nscat 
         mmc_rep_a (ic, ii, jj) = a 
         mmc_rep_b (ic, ii, jj) = b 
         mmc_rep_c (ic, ii, jj) = c 
         mmc_rep_m (ic, ii, jj) = m 
         mmc_rep_a (ic, jj, ii) = a 
         mmc_rep_b (ic, jj, ii) = b 
         mmc_rep_c (ic, jj, ii) = c 
         mmc_rep_m (ic, jj, ii) = m 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE mmc_set_rep                   
!*****7*****************************************************************
!     SUBROUTINE mmc_set_vec (ic, is, js, werte, maxw) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
!     USE discus_config_mod 
!     USE crystal_mod 
!     USE mmc_mod 
!     USE errlist_mod 
!     IMPLICIT none 
!                                                                       
!      
!                                                                       
!     INTEGER maxw 
!                                                                       
!     INTEGER ic, is, js, ii, jj 
!     INTEGER i, nv 
!     REAL werte (maxw) 
!                                                                       
!                                                                       
!     IF (is.ne. - 1.and.js.ne. - 1) then 
!        mmc_nvec (ic, js, is) = mmc_nvec (ic, js, is) + 1 
!        mmc_nvec (ic, is, js) = mmc_nvec (ic, js, is) 
!        nv = mmc_nvec (ic, is, js) 
!        DO i = 1, 4 
!        mmc_vec (i, nv, ic, is, js) = werte (i) 
!        mmc_vec (i, nv, ic, js, is) = werte (i) 
!        ENDDO 
!     ELSEIF (is.eq. - 1.and.js.ne. - 1) then 
!        DO ii = 0, cr_nscat 
!        mmc_nvec (ic, js, ii) = mmc_nvec (ic, js, ii) + 1 
!        mmc_nvec (ic, ii, js) = mmc_nvec (ic, js, is) 
!        nv = mmc_nvec (ic, ii, js) 
!        DO i = 1, 4 
!        mmc_vec (i, nv, ic, ii, js) = werte (i) 
!        mmc_vec (i, nv, ic, js, ii) = werte (i) 
!        ENDDO 
!        ENDDO 
!     ELSEIF (is.ne. - 1.and.js.eq. - 1) then 
!        DO ii = 0, cr_nscat 
!        mmc_nvec (ic, is, ii) = mmc_nvec (ic, is, ii) + 1 
!        mmc_nvec (ic, ii, is) = mmc_nvec (ic, is, ii) 
!        nv = mmc_nvec (ic, ii, is) 
!        DO i = 1, 4 
!        mmc_vec (i, nv, ic, ii, is) = werte (i) 
!        mmc_vec (i, nv, ic, is, ii) = werte (i) 
!        ENDDO 
!        ENDDO 
!     ELSE 
!        DO ii = 0, cr_nscat 
!        DO jj = ii, cr_nscat 
!        mmc_nvec (ic, jj, ii) = mmc_nvec (ic, jj, ii) + 1 
!        mmc_nvec (ic, ii, jj) = mmc_nvec (ic, jj, ii) 
!        nv = mmc_nvec (ic, ii, jj) 
!        DO i = 1, 4 
!        mmc_vec (i, nv, ic, ii, jj) = werte (i) 
!        mmc_vec (i, nv, ic, jj, ii) = werte (i) 
!        ENDDO 
!        ENDDO 
!        ENDDO 
!     ENDIF 
!                                                                       
!     END SUBROUTINE mmc_set_vec                    
!*****7*****************************************************************
      SUBROUTINE mmc_set_mode (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     Sets MMC    mode                                                  
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE rmc_mod 
      USE mmc_mod 
      USE modify_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER :: ianz, imode=MC_MOVE_NONE, i 
      INTEGER is 
      REAL werte (maxw) 
      REAL sump 
!                                                                       
      IF (ianz.ge.1) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3) .eq.'SHI') then 
            imode = MC_MOVE_DISP 
         ELSEIF (cpara (1) (1:3) .eq.'SWD') then 
            imode = MC_MOVE_SWDISP 
         ELSEIF (cpara (1) (1:3) .eq.'SWC') then 
            imode = MC_MOVE_SWCHEM 
         ELSEIF (cpara (1) (1:3) .eq.'INV') then 
            imode = MC_MOVE_INVDISP 
         ELSE 
            ier_typ = ER_RMC 
            ier_num = - 9 
         ENDIF 
!                                                                       
         IF (ianz.eq.1) then 
            mmc_local (imode) = rmc_local_all 
         ELSEIF (ianz.ge.2) then 
            CALL do_cap (cpara (2) ) 
            IF (cpara (2) (1:1) .eq.'A') then 
               mmc_local (imode) = rmc_local_all 
            ELSEIF (cpara (2) (1:1) .eq.'L') then 
               mmc_local (imode) = rmc_local_loc 
            ELSEIF (cpara (2) (1:2) .eq.'SL') then 
               mmc_local (imode) = rmc_local_locsite 
            ELSEIF (cpara (2) (1:2) .eq.'SI') then 
               mmc_local (imode) = rmc_local_site 
            ELSE 
               ier_typ = ER_RMC 
               ier_num = - 9 
            ENDIF 
            IF (ianz.gt.2) then 
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL get_iscat (ianz, cpara, lpara, werte, maxw, .false.) 
               IF (ier_num.ne.0) return 
!                                                                       
               IF (ier_num.eq.0) then 
                  IF (werte (1) .eq. - 1) then 
                     DO i = 0, cr_nscat 
                     mmc_allowed (i) = .true. 
                     ENDDO 
                  ELSE 
                     DO i = 1, ianz 
                     is = nint (werte (i) ) 
                     IF (is.ge.0.and.is.le.cr_nscat) then 
                        mmc_allowed (is) = .true. 
                     ELSE 
                        ier_num = - 27 
                        ier_typ = ER_APPL 
                     ENDIF 
                     ENDDO 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
!     --Set probabilities for the different moves                       
!                                                                       
         mmc_move_prob (imode) = werte (1) 
         sump = 0.0 
         DO i = 1, MC_N_MOVE 
         sump = sump + mmc_move_prob (i) 
         ENDDO 
         IF (sump.gt.0) then 
            mmc_move_cprob (1) = mmc_move_prob (1) / sump 
            DO i = 2, MC_N_MOVE 
            mmc_move_cprob (i) = mmc_move_cprob (i - 1) + mmc_move_prob &
            (i) / sump                                                  
            ENDDO 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE mmc_set_mode                   
!*****7*****************************************************************
      SUBROUTINE mmc_run_multi 
!+                                                                      
!     This is the MC routine for multiple energy calculations           
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_menu
      USE celltoindex_mod
      USE atom_name
      USE metric_mod
      USE mc_mod 
      USE rmc_menu
      USE rmc_mod 
      USE mmc_mod 
      USE modify_func_mod
      USE random_mod
      USE rmc_sup_mod
!
      USE debug_mod 
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = DEF_MAXSCAT) 
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = CHEM_MAX_NEIG) 
!                                                                       
      CHARACTER(24) c_energy (0:MC_N_ENERGY) 
      REAL werte (maxw), wwerte (maxw), wwwerte (maxw) 
      REAL start, zeit, seknds 
      REAL disp1, disp2 
      REAL disp (3, 0:CHEM_MAX_NEIG, 2) 
      REAL idir (3), jdir (3) 
      REAL rdi (CHEM_MAX_COR) 
      REAL rdj (CHEM_MAX_COR) 
      REAL delta 
      REAL :: posz (3) = 0.0
      REAL :: posz2 (3) = 0.0
      REAL :: rel_cycl    ! how far are we in the desired number of cycles
      REAL patom (3, 0:CHEM_MAX_NEIG, CHEM_MAX_CENT) 
      REAL :: rrrr 
      INTEGER iatom (0:CHEM_MAX_NEIG, CHEM_MAX_CENT) 
      INTEGER igen, itry, iacc_good, iacc_bad 
      INTEGER isel (CHEM_MAX_ATOM) 
      INTEGER :: iselz=0, iselz2=0
      INTEGER lbeg (3) 
      INTEGER ic, is (2), iz1 (3), iz2 (3), ie 
      INTEGER iz (2, 3) 
      INTEGER i, j, natoms, ia
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER ncent 
      INTEGER icent 
      INTEGER iscat 
      INTEGER :: NALLOWED   ! Current size mmc_allowed
      INTEGER zh, zm, zs 
      LOGICAL loop, laccept, done 
      LOGICAL valid_e 
      LOGICAL :: valid_all = .false.
      LOGICAL lout 
!                                                                       
      REAL v (3) 
      REAL z 
!                                                                       
      REAL e_old (0:MC_N_ENERGY) 
      REAL e_new (0:MC_N_ENERGY) 
!                                                                       
!      REAL mmc_energy_angle 
!      REAL mmc_energy_occ 
!      REAL mmc_energy_dis 
!      REAL mmc_energy_spr 
!      REAL mmc_energy_vec 
!      REAL mmc_energy_len 
!      REAL mmc_energy_buck 
!      REAL mmc_energy_rep 
      REAL ran1, gasdev 
!      INTEGER len_str 
!     LOGICAL atom_allowed 
!     LOGICAL check_select_status 
!     REAL do_blen 
!     REAL skalpro 
!                                                                       
      DATA c_energy /                    &
           '                        ',   &
           'Occupation correlation  ',   &
           'Displacement correlation',   &
           'Displacement energy     ',   &
           'Angular energy          ',   &
           'Vector energy           ',   &
           'Bond length energy      ',   &
           'Lennard Jones Potential ',   &
           'Buckingham Potential    ',   &
           'Repulsive     Potential ' /                                                              !                                                                       
      DO i = 0, MC_N_ENERGY 
      n_e_av_p (i) = 0 
      n_e_av_m (i) = 0 
      n_e_av_z (i) = 0 
      e_aver_p (i) = 0.0 
      e_aver_m (i) = 0.0 
      ENDDO 
      DO i = 1, maxw 
      werte (i) = 0.0 
      wwerte (i) = 0.0 
      wwwerte (i) = 0.0 
      ENDDO 
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
!     Normalize the correlation directions                              
!                                                                       
      DO ic = 1, chem_ncor 
      IF (mmc_cor_energy (0, MC_DISP) ) then 
         DO i = 1, 3 
         idir (i) = chem_dir (i, 1, ic) 
         jdir (i) = chem_dir (i, 2, ic) 
         ENDDO 
         rdi (ic) = skalpro (idir, idir, cr_gten) 
         rdj (ic) = skalpro (jdir, jdir, cr_gten) 
         rdi (ic) = sqrt (rdi (ic) ) 
         rdj (ic) = sqrt (rdj (ic) ) 
      ENDIF 
      ENDDO 
!                                                                       
!     Initialize the different energies                                 
!                                                                       
      lout = .false. 
      CALL mmc_correlations (lout, 0.0) 
      IF (ier_num.ne.0) return 
!
!     current size of mmc_allowed
!
      NALLOWED = UBOUND(mmc_allowed,1)
!                                                                       
      IF (mo_cyc.eq.0) loop = .false. 
!                                                                       
!------ Main MC loop here                                               
!                                                                       
      start = seknds (0.0) 
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
!                                                                       
!------ - selecting sites                                               
!                                                                       
      natoms = 1 
   10 CONTINUE 
!                                                                       
!     --- Choose a move at random                                       
!                                                                       
      z = ran1 (idum) 
      i = 1 
      DO while (z.gt.mmc_move_cprob (i) .and.i.lt.MC_N_MOVE) 
      i = i + 1 
      ENDDO 
      mmc_move = i 
      mo_local = mmc_local (i) 
!                                                                       
      IF (mmc_move.eq.MC_MOVE_DISP.or.mmc_move.eq.MC_MOVE_INVDISP) then 
         IF (mmc_l_limited) then 
            CALL mmc_limit_selection (isel, natoms) 
         ELSE 
            isel (1) = int (ran1 (idum) * cr_natoms) + 1 
         ENDIF 
         iselz = isel (1) 
         IF (isel (1) .gt.cr_natoms.or.isel (1) .lt.1) goto 10 
      laccept = mmc_allowed (cr_iscat (isel (1) ) ) .and.check_select_st&
     &atus (.true., cr_prop (isel (1) ),  cr_sel_prop)                  
         IF (cr_ncatoms.gt.0) then 
            CALL indextocell (isel (1), iz1, is (1) ) 
         ELSE 
            DO i = 1, 3 
               iz1 (i) = 1 
            ENDDO 
            is(1) = 1
         ENDIF 
         DO i = 1, 3 
         iz (1, i) = iz1 (i) 
         ENDDO 
      ELSEIF (mmc_move.eq.MC_MOVE_SWDISP) then 
         natoms = 2 
         CALL rmc_select (mo_local, isel, iz1, iz2, is (1), is (2) , &
                          NALLOWED, mmc_allowed)                                                       
         iselz = isel (1) 
         iselz2 = isel (2) 
         DO i = 1, 3 
         iz (1, i) = iz1 (i) 
         iz (2, i) = iz2 (i) 
         ENDDO 
         laccept = mmc_allowed (cr_iscat (isel (1) ) ) .and.mmc_allowed &
         (cr_iscat (isel (2) ) ) .and.check_select_status (.true.,      &
         cr_prop (isel (1) ), cr_sel_prop) .and.check_select_status (   &
         .true., cr_prop (isel (2) ), cr_sel_prop)                      
      ELSEIF (mmc_move.eq.MC_MOVE_SWCHEM) then 
         natoms = 2 
         CALL rmc_select (mo_local, isel, iz1, iz2, is (1), is (2) , &
                          NALLOWED, mmc_allowed)                                                       
!RBN_REP
!isel(1) = 34
!isel(2) = 37
! call indextocell (isel(1), iz1,is(1))
! call indextocell (isel(2), iz1,is(1))
!iz1(1)  =  2
!iz1(2)  =  2
!iz1(3)  =  1
!iz2(1)  =  2
!iz2(2)  =  5
!iz2(3)  =  1
!is (1)  =  1
!is (2)  =  1
         iselz = isel (1) 
         iselz2 = isel (2) 
         DO i = 1, 3 
         iz (1, i) = iz1 (i) 
         iz (2, i) = iz2 (i) 
         ENDDO 
         laccept = cr_iscat (isel (1) ) .ne.cr_iscat (isel (2) ) .and.                  &
                 ( mmc_allowed (cr_iscat (isel (1) ) ) .and.                            &
                   mmc_allowed (cr_iscat (isel (2) ) )      )    .and.                  &
                   check_select_status (.true., cr_prop (isel (1) ), cr_sel_prop) .and. &
                   check_select_status (.true., cr_prop (isel (2) ), cr_sel_prop)                              
      ENDIF 
!                                                                       
!-----      ----Check whether geometrical constrains apply              
!                                                                       
      IF (laccept) then 
         CALL check_geometry_multi (isel, natoms, laccept) 
      ENDIF 
!                                                                       
!----- -- Try move                                                      
!                                                                       
      IF (laccept) then 
         itry = itry + 1 
         done = .false. 
!                                                                       
!     ----Calculate old energy                                          
!                                                                       
         DO ic = 1, chem_ncor 
         DO ie = 0, MC_N_ENERGY 
         e_old (ie) = 0.0 
         ENDDO 
         ENDDO 
!                                                                       
!     ----Loop over all modified atoms                                  
!                                                                       
         DO ia = 1, natoms 
!                                                                       
!     ------Set the assumption of at least one propper energy to FALSE  
!           Will be set to TRUE if at least one energy calculation      
!               is fine                                                 
!                                                                       
         valid_all = .false. 
!                                                                       
!     ------Loop over all defined neighbour interactions                
!                                                                       
         DO ic = 1, chem_ncor 
         CALL chem_neighbour_multi (isel (ia), ic, iatom, patom, natom, &
         ncent, maxatom)                                                
         DO icent = 1, ncent 
         IF (natom (icent) .gt.0) then 
!                                                                       
!------ ------- Calc old energy for those energies that are affected by 
!               the move                                                
!                                                                       
            IF (mmc_cor_energy (ic, MC_OCC) ) then 
!                                                                       
!     ------- Occupation Correlation, only if both atoms are different  
!     ------- and the move is switch chemistry                          
!                                                                       
               IF (mmc_move.eq.MC_MOVE_SWCHEM) then 
                  IF (cr_iscat (isel (1) ) .ne.cr_iscat (isel (2) ) )   &
                  then                                                  
                     e_old (MC_OCC) = e_old (MC_OCC) + mmc_energy_occ ( &
                     isel, ia, ic, iatom, icent, natom, valid_e) 
                     valid_all = valid_all.or.valid_e 
                  ENDIF 
               ENDIF 
            ENDIF 
!                                                                       
!     ------- Displacement correlation                                  
!                                                                       
            IF (mmc_cor_energy (ic, MC_DISP) ) then 
               e_old (MC_DISP) = 0.0 
               valid_all = .true. 
            ENDIF 
!                                                                       
!     ------- Displacement (Hooke's law) '                              
!                                                                       
            IF (mmc_cor_energy (ic, MC_SPRING) ) then 
               e_old (MC_SPRING) = e_old (MC_SPRING) + mmc_energy_spr ( &
               isel, ia, ic, iatom, patom, icent, natom, valid_e)       
               valid_all = valid_all.or.valid_e 
            ENDIF 
!                                                                       
!     ------- Angular energy                                            
!                                                                       
            IF (mmc_cor_energy (ic, MC_ANGLE) ) then 
               e_old (MC_ANGLE) = e_old (MC_ANGLE) + mmc_energy_angle ( &
               ic, iatom, patom, icent, natom, valid_e)           
               valid_all = valid_all.or.valid_e 
            ENDIF 
!                                                                       
!     ------- Vector energy, i.e. directional bond length               
!                                                                       
!           IF (mmc_cor_energy (ic, MC_VECTOR) ) then 
!              e_old (MC_VECTOR) = e_old (MC_VECTOR) + mmc_energy_vec ( &
!              isel, ia, ic, iatom, patom, icent, natom, valid_e)       
!              valid_all = valid_all.or.valid_e 
!           ENDIF 
!                                                                       
!     ------- Bond length energy, non - directional                     
!                                                                       
            IF (mmc_cor_energy (ic, MC_BLEN) ) then 
               CONTINUE 
            ENDIF 
!                                                                       
!     ------- Lennard-Jones Potential                                   
!                                                                       
            IF (mmc_cor_energy (ic, MC_LENNARD) ) then 
               e_old (MC_LENNARD) = e_old (MC_LENNARD) + mmc_energy_len &
               (isel, ia, ic, iatom, patom, icent, natom, valid_e)
               valid_all = valid_all.or.valid_e 
            ENDIF 
!
!     ------- Repulsive     Potential                                   
!
            IF (mmc_cor_energy (ic, MC_REPULSIVE) ) then 
               e_old (MC_REPULSIVE) = e_old (MC_REPULSIVE) +            &
                     mmc_energy_rep (isel, ia, ic, iatom, patom, icent, &
                                     natom, valid_e)
               valid_all = valid_all.or.valid_e 
            ENDIF 
!                                                                       
!     ------- Buckingham    Potential                                   
!                                                                       
            IF (mmc_cor_energy (ic, MC_BUCKING) ) then 
               e_old (MC_BUCKING) = e_old (MC_BUCKING) +                &
               mmc_energy_buck (isel, ia, ic, iatom, patom, icent,      &
               natom, valid_e)                                          
               valid_all = valid_all.or.valid_e 
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!     ----If a proper energy was found, modify the atom                 
!                                                                       
         IF (valid_all) then 
!                                                                       
!     ----Perform the modifications of the atoms                        
!                                                                       
            IF (mmc_move.eq.MC_MOVE_DISP) then 
!                                                                       
!     ------Displace atoms from current positions                       
!                                                                       
!                                                                       
!     ------Modify the central atom                                     
!                                                                       
!               CALL indextocell (isel (1), iz1, is (1) ) 
               IF(mo_maxmove(4, cr_iscat(iselz) )==0.0) THEN
                  DO i = 1, 3 
                     disp(i, 0, 1) = gasdev(mo_maxmove(i, cr_iscat(iselz) ))
                     posz (i) = cr_pos (i, iselz) 
                     cr_pos (i, iselz) = cr_pos (i, iselz) + disp (i, 0, 1) 
                  ENDDO 
               ELSE
!     -- Move along a vector direction
                  rrrr = gasdev(mo_maxmove(4, cr_iscat(iselz) ))
                  DO i = 1, 3 
                     disp(i, 0, 1) = rrrr* (mo_maxmove(i, cr_iscat(iselz) ))
                     posz (i) = cr_pos (i, iselz) 
                     cr_pos (i, iselz) = cr_pos (i, iselz) + disp (i, 0, 1) 
                   ENDDO 
               ENDIF
            ELSEIF (mmc_move.eq.MC_MOVE_SWDISP) then 
!                                                                       
!-----      ------Switch displacement of two selected atoms and their   
!             neighbours                                                
!                                                                       
               CALL indextocell (isel (1), iz1, is (1) ) 
               CALL indextocell (isel (2), iz2, is (2) ) 
!                                                                       
!     ------Modify the central atom                                     
!                                                                       
               DO j = 1, 3 
               disp1 = cr_pos (j, isel (1) ) - chem_ave_pos (j, is (1) )&
               - float (iz1 (j) - 1) - cr_dim0 (j, 1)                   
               disp2 = cr_pos (j, isel (2) ) - chem_ave_pos (j, is (2) )&
               - float (iz2 (j) - 1) - cr_dim0 (j, 1)                   
               disp (j, 0, 1) = - disp1 + disp2 
               disp (j, 0, 2) = - disp2 + disp1 
               posz (j) = cr_pos (j, isel (1) ) 
               cr_pos (j, isel (1) ) = cr_pos (j, isel (1) ) + disp (j, &
               0, 1)                                                    
               posz2 (j) = cr_pos (j, isel (2) ) 
               cr_pos (j, isel (2) ) = cr_pos (j, isel (2) ) + disp (j, &
               0, 2)                                                    
               ENDDO 
            ELSEIF (mmc_move.eq.MC_MOVE_INVDISP) then 
!                                                                       
!-----      ------Switch displacement of a selected atom                
!                                                                       
               CALL indextocell (isel (1), iz1, is (1) ) 
               DO j = 1, 3 
               disp1 = cr_pos (j, isel (1) ) - chem_ave_pos (j, is (1) )&
               - float (iz1 (j) - 1) - cr_dim0 (j, 1)                   
               disp (j, 0, 1) = - 2 * disp1 
               posz (j) = cr_pos (j, isel (1) ) 
               cr_pos (j, isel (1) ) = cr_pos (j, isel (1) ) + disp (j, &
               0, 1)                                                    
               ENDDO 
            ELSEIF (mmc_move.eq.MC_MOVE_SWCHEM) then 
!                                                                       
!     ------Switch the Chemistry of two selected atoms                  
!                                                                       
               iscat = cr_iscat (isel (2) ) 
               cr_iscat (isel (2) ) = cr_iscat (isel (1) ) 
               cr_iscat (isel (1) ) = iscat 
               iscat = cr_prop (isel (2) ) 
               cr_prop (isel (2) ) = cr_prop (isel (1) ) 
               cr_prop (isel (1) ) = iscat 
!                                                                       
!-----      ------End of Modification of atoms according to different   
!             moves                                                     
!                                                                       
            ENDIF 
!                                                                       
!                                                                       
!     ----Calculate new energy                                          
!                                                                       
            DO ie = 0, MC_N_ENERGY 
            e_new (ie) = 0.0 
            ENDDO 
!                                                                       
!     ------Set the assumption of at least one propper energy to FALSE  
!           Will be set to TRUE if at least one energy calculation      
!             is fine                                                   
!                                                                       
            valid_all = .false. 
!                                                                       
!     ----Loop over all modified atoms                                  
!                                                                       
            DO ia = 1, natoms 
            DO i = 1, 3 
            v (i) = cr_pos (i, isel (ia) ) - chem_ave_pos (i, is (ia) ) &
            - float (iz (ia, i) - 1) - cr_dim0 (i, 1)                   
            v (i) = v (i) - disp (i, 0, ia) 
            ENDDO 
!                                                                       
!     ------Loop over all defined neighbour interactions                
!                                                                       
            DO ic = 1, chem_ncor 
            CALL chem_neighbour_multi (isel (ia), ic, iatom, patom,     &
            natom, ncent, maxatom)                                      
            DO icent = 1, ncent 
            IF (natom (icent) .gt.0) then 
!                                                                       
!------ ------- Calc new energy for those energies that are affected    
!               by the move                                             
!                                                                       
               IF (mmc_cor_energy (ic, MC_OCC) ) then 
!                                                                       
!     ------- Occupation Correlation                                    
!                                                                       
                  valid_all = valid_all.or.valid_e 
                  e_new (MC_OCC) = e_new (MC_OCC) + mmc_energy_occ (    &
                  isel, ia, ic, iatom, icent, natom, valid_e)    
                  valid_all = valid_all.or.valid_e 
               ENDIF 
!                                                                       
!     ------- Displacement correlation                                  
!                                                                       
               IF (mmc_cor_energy (ic, MC_DISP) ) then 
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
                  e_new (MC_DISP) = e_new (MC_DISP) + mmc_energy_dis (  &
                  isel, ia, ic, iatom, icent, natom, jdir, delta,&
                  rdj, valid_e)                                         
                  loop = .false. 
               ENDIF 
!                                                                       
!     ------- Displacement (Hooke's law) '                              
!                                                                       
               IF (mmc_cor_energy (ic, MC_SPRING) ) then 
                  e_new (MC_SPRING) = e_new (MC_SPRING) +               &
                  mmc_energy_spr (isel, ia, ic, iatom, patom, icent,    &
                  natom, valid_e)                                       
                  valid_all = valid_all.or.valid_e 
               ENDIF 
!                                                                       
!     ------- Angular energy                                            
!                                                                       
               IF (mmc_cor_energy (ic, MC_ANGLE) ) then 
                  e_new (MC_ANGLE) = e_new (MC_ANGLE) +                 &
                  mmc_energy_angle (ic, iatom, patom, icent,      &
                  natom, valid_e)                                       
                  valid_all = valid_all.or.valid_e 
               ENDIF 
!                                                                       
!     ------- Vector energy, i.e. directional bond length               
!                                                                       
!              IF (mmc_cor_energy (ic, MC_VECTOR) ) then 
!                 e_new (MC_VECTOR) = e_new (MC_VECTOR) +               &
!                 mmc_energy_vec (isel, ia, ic, iatom, patom, icent,    &
!                 natom, valid_e)                                       
!                 valid_all = valid_all.or.valid_e 
!              ENDIF 
!                                                                       
!     ------- Bond length energy, non - directional                     
!                                                                       
               IF (mmc_cor_energy (ic, MC_BLEN) ) then 
                  CONTINUE 
               ENDIF 
!                                                                       
!     ------- Lennard-Jones Potential                                   
!                                                                       
                                                                        
               IF (mmc_cor_energy (ic, MC_LENNARD) ) then 
                  e_new (MC_LENNARD) = e_new (MC_LENNARD) +             &
                  mmc_energy_len (isel, ia, ic, iatom, patom, icent,    &
                  natom, valid_e)
                  valid_all = valid_all.or.valid_e 
               ENDIF 
!
!     ------- Repulsive     Potential                                   
!
               IF (mmc_cor_energy (ic, MC_REPULSIVE) ) then 
                  e_new (MC_REPULSIVE) = e_new (MC_REPULSIVE) +         &
                     mmc_energy_rep (isel, ia, ic, iatom, patom, icent, &
                                     natom, valid_e)
                  valid_all = valid_all.or.valid_e 
               ENDIF 
!                                                                       
!     ------- Buckingham Potential                                      
!                                                                       
               IF (mmc_cor_energy (ic, MC_BUCKING) ) then 
                  e_new (MC_BUCKING) = e_new (MC_BUCKING) +             &
                  mmc_energy_buck (isel, ia, ic, iatom, patom, icent,   &
                  natom, valid_e)                                       
                  valid_all = valid_all.or.valid_e 
               ENDIF 
            ENDIF 
            ENDDO 
            ENDDO 
            ENDDO 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------The comparison of the energies is done only if a proper     
!           new energy was found. Otherwise the move is automatically   
!           rejected. This might happen, if an atom has moved outside   
!           the allowed sphere of influence.                            
!                                                                       
            IF (valid_all) then 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
               CALL mmc_test_multi (iacc_good, iacc_bad, e_new, e_old,  &
               laccept)                                                 
            ELSE 
               laccept = .false. 
            ENDIF 
!                                                                       
!     ----The move was not accepted, move atoms back to old places      
!                                                                       
            IF (.not.laccept) then 
               IF (mmc_move.ne.MC_MOVE_SWCHEM) then 
                  DO i = 1, 3 
                  cr_pos (i, iselz) = posz (i) 
                  ENDDO 
                  IF (mmc_move.eq.MC_MOVE_SWDISP) then 
                     DO i = 1, 3 
                     cr_pos (i, iselz2) = posz2 (i) 
                     ENDDO 
                  ENDIF 
               ELSE 
!                                                                       
!     ------Switch the Chemistry of two selected atoms                  
!                                                                       
                  iscat = cr_iscat (isel (2) ) 
                  cr_iscat (isel (2) ) = cr_iscat (isel (1) ) 
                  cr_iscat (isel (1) ) = iscat 
                  iscat = cr_prop (isel (2) ) 
                  cr_prop (isel (2) ) = cr_prop (isel (1) ) 
                  cr_prop (isel (1) ) = iscat 
               ENDIF 
            ELSE     ! Move accepted check periodic bounday conditions
               IF(.NOT.chem_quick .AND.                                 &
                  (chem_period(1).OR.chem_period(2).OR.chem_period(3))) THEN
!                 normal and periodic mode
                  CALL chem_apply_period(iselz, .TRUE.)
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
!     --End of modification of atoms if a proper "old" energy was found 
!                                                                       
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
         done = .true. 
         WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
!                                                                       
!     ----New mmc_correlations for all energies                         
!                                                                       
         lout = .true. 
         rel_cycl = float(itry)/float(mo_cyc)
         CALL mmc_correlations (lout, rel_cycl)
                                                                        
!        IF (mmc_cor_energy (0, MC_VECTOR) ) then 
!                                                                       
!     ------VECTOR  energy was selected, give feedback                  
!                                                                       
!           l = 2 
!           cpara (1) = 'all' 
!           cpara (2) = 'all' 
!           lpara (1) = 3 
!           lpara (2) = 3 
!           CALL chem_vector_multi (l, cpara, lpara, verte, maxw,       &
!           .false.)                                                    
!           DO k = 1, chem_ncor 
!           IF (mmc_cor_energy (k, MC_VECTOR) ) then 
!              DO i = 0, cr_nscat 
!              DO j = i, cr_nscat 
!              at_name_i = at_name (i) 
!              at_name_j = at_name (j) 
!              DO l = 1, mmc_nvec (k, i, j) 
!     WRITE (output_io, 2400) mmc_vec (1, l, k, i, j),  chem_vect_ave (1&
!    &, l, k, i, j),  chem_vect_sig (1, l, k, i, j)                     
!     WRITE (output_io, 2410) k, at_name_i, at_name_j, mmc_vec (2, l, k,&
!    & i, j),  chem_vect_ave (2, l, k, i, j),  chem_vect_sig (2, l, k, i&
!    &, j)                                                              
!     WRITE (output_io, 2420) mmc_vec (3, l, k, i, j),  chem_vect_ave (3&
!    &, l, k, i, j),  chem_vect_sig (3, l, k, i, j)                     
!              ENDDO 
!              ENDDO 
!              ENDDO 
!           ENDIF 
!           ENDDO 
!        ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ Loop finished                                                   
!                                                                       
!                                                                       
      WRITE (output_io, 3000) 
      WRITE (output_io, 2000) igen, itry, iacc_good, iacc_bad 
      lout = .true. 
      CALL mmc_correlations (lout, 1.0) 
!     IF (mmc_cor_energy (0, MC_VECTOR) ) then 
!                                                                       
!     --VECTOR  energy was selected, give feedback                      
!                                                                       
!        l = 2 
!        cpara (1) = 'all' 
!        cpara (2) = 'all' 
!        lpara (1) = 3 
!        lpara (2) = 3 
!        CALL chem_vector_multi (l, cpara, lpara, verte, maxw, .false.) 
!        DO i = 0, cr_nscat 
!        DO j = i, cr_nscat 
!        at_name_i = at_name (i) 
!        at_name_j = at_name (j) 
!        DO k = 1, chem_ncor 
!        IF (mmc_cor_energy (k, MC_VECTOR) ) then 
!           DO l = 1, mmc_nvec (k, i, j) 
!     WRITE (output_io, 2400) mmc_vec (1, l, k, i, j),  chem_vect_ave (1&
!    &, l, k, i, j),  chem_vect_sig (1, l, k, i, j)                     
!           WRITE (output_io, 2410) k, at_name_i, at_name_j, mmc_vec (2,&
!           l, k, i, j), chem_vect_ave (2, l, k, i, j), chem_vect_sig ( &
!           2, l, k, i, j)                                              
!     WRITE (output_io, 2420) mmc_vec (3, l, k, i, j),  chem_vect_ave (3&
!    &, l, k, i, j),  chem_vect_sig (3, l, k, i, j)                     
!           ENDDO 
!        ENDIF 
!        ENDDO 
!        ENDDO 
!        ENDDO 
!     ENDIF 
!                                                                       
!     Give average energy changes for the different energy terms        
!                                                                       
      WRITE ( output_io, 5000) 
      DO i = 1, MC_N_ENERGY 
      IF (n_e_av_p (i) .gt.0) then 
         e_aver_p (i) = e_aver_p (i) / float (n_e_av_p (i) ) 
      ENDIF 
      IF (n_e_av_m (i) .gt.0) then 
         e_aver_m (i) = e_aver_m (i) / float (n_e_av_m (i) ) 
      ENDIF 
      WRITE ( output_io, 5010) c_energy (i), n_e_av_m (i), e_aver_m (i),        &
      n_e_av_z (i), n_e_av_p (i), e_aver_p (i)                          
      n_e_av_p (i) = 0 
      n_e_av_m (i) = 0 
      n_e_av_z (i) = 0 
      e_aver_p (i) = 0.0 
      e_aver_m (i) = 0.0 
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
 2000 FORMAT (/,' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I7,        &
     &          ' / ',I7,'  MC moves ')                                 
 3000 FORMAT (/,' --- Final multiple energy configuration ---') 
 4000 FORMAT (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/     &
     &          ' Time/cycle   : ',F9.3,' sec',/)                       
 5000 FORMAT (/,' --- Average energy changes ---',//,                   &
     &        ' Energy type',14x,'number',17x,'number',2x,'number',/,   &
     &        29x,'E < 0',5x,' <E>',11x,'E = 0',5x,'E > 0',4x,' <E>')   
 5010 FORMAT(a24, 2x,i8,f15.6,2x,i8,2x,i8,f15.6) 
!                                                                       
      END SUBROUTINE mmc_run_multi                  
!*****7*****************************************************************
      SUBROUTINE mmc_test_multi (iacc_good, iacc_bad, e_new, e_old,     &
      laccept)                                                          
!+                                                                      
!     Tests performed MC move                                           
!-                                                                      
      USE discus_config_mod 
      USE mc_mod 
      USE mmc_mod 
      USE random_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iacc_good, iacc_bad 
      INTEGER i 
!                                                                       
      REAL e_old (0:MC_N_ENERGY) 
      REAL e_new (0:MC_N_ENERGY) 
      REAL e_del 
      REAL e_delta 
      LOGICAL laccept 
!                                                                       
      REAL ran1 
!                                                                       
      e_del = 0 
      DO i = 1, MC_N_ENERGY 
      IF (mmc_cor_energy (0, i) ) then 
         e_delta = (e_new (i) - e_old (i) ) 
         e_del = e_del + e_delta 
         IF (e_delta.gt.0) then 
            e_aver_p (i) = e_aver_p (i) + e_delta 
            n_e_av_p (i) = n_e_av_p (i) + 1 
         ELSEIF (e_delta.lt.0) then 
            e_aver_m (i) = e_aver_m (i) + e_delta 
            n_e_av_m (i) = n_e_av_m (i) + 1 
         ELSE 
            n_e_av_z (i) = n_e_av_z (i) + 1 
         ENDIF 
      ENDIF 
      ENDDO 
      IF (e_del.lt.0) then 
         laccept = .true. 
!     ELSEIF(e_del.eq.0) then                                           
!       laccept = .true.                                                
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
      END SUBROUTINE mmc_test_multi                 
!*****7*****************************************************************
      REAL function mmc_energy_occ (isel, ia, ic, iatom, icent,  &
      natom, valid_e)                                                   
!+                                                                      
!     Calculates the energy for chemical disorder                       
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE mc_mod 
      USE mmc_mod 
      USE modify_mod
      USE modify_func_mod
      USE param_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
!                                                                       
      INTEGER isel (chem_max_atom) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER is, js, ic, in, ia
      INTEGER in_a, in_e 
      INTEGER ncalc 
      INTEGER ival1
!                                                                       
!     LOGICAL check_select_status 
!                                                                       
      mmc_energy_occ = 0.0 
      ncalc = 0 
      valid_e = .false. 
!                                                                       
      IF (chem_ctyp(ic) == CHEM_VEC    .or. &
          chem_ctyp(ic) == CHEM_ENVIR  .or. &
          chem_ctyp(ic) == CHEM_RANGE  .or. &
          chem_ctyp(ic) == CHEM_DIST   .or. &
          chem_ctyp(ic) == CHEM_CON        )   then                                               
!                                                                       
         IF (natom (icent) .ne.0) then 
            IF (isel (ia) .eq.iatom (0, icent) ) then 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
               DO in = in_a, in_e 
                  IF (check_select_status (.true., cr_prop (iatom (in,  &
                                           icent) ), cr_sel_prop) ) then                         
                     ival1 = 0 
                     js = cr_iscat (iatom (in, icent) ) 
                     ival1 = sign(1,mmc_pair(ic, MC_OCC,is,js))
                     mmc_energy_occ = mmc_energy_occ +                  &
                                      mmc_depth (ic,MC_OCC, 0, 0) * ival1
                  ENDIF 
               ENDDO 
!              DO jjs = 0, cr_nscat 
!              IF (mmc_pair (ic, MC_OCC, is, jjs) == -1 ) then 
!                 DO in = in_a, in_e 
!                 IF (check_select_status (.true., cr_prop (iatom (in,  &
!                 icent) ), cr_sel_prop) ) then                         
!                    ival1 = 0 
!                    js = cr_iscat (iatom (in, icent) ) 
!                    IF (is.eq.js) then 
!                       ival1 = 1 
!DBG      ELSE                                                          
!                    ELSEIF (mmc_pair (ic, MC_OCC, is, js) == -1 ) then 
!                       ival1 = - 1 
!                    ENDIF 
!write(*,*) ' Neigbour     : ', iatom(in,icent), cr_iscat(iatom(in,icent)), mmc_pair (ic, MC_OCC, is, js), ival1, &
!                     mmc_depth (ic,   MC_OCC, 0, 0) * ival1                              
!                    mmc_energy_occ = mmc_energy_occ + mmc_depth (ic,   &
!                    MC_OCC, 0, 0) * ival1                              
!                 ENDIF 
!                 ENDDO 
!              ENDIF 
!              ENDDO 
            ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
               in_a = 0 
               in_e = 0 
               is = cr_iscat (isel (ia) ) 
               in = 0
                  IF (check_select_status (.true., cr_prop (iatom (in,  &
                                           icent) ), cr_sel_prop) ) then                         
                     ival1 = 0 
                     js = cr_iscat (iatom (in, icent) ) 
                     ival1 = mmc_pair(ic, MC_OCC,is,js)
                     mmc_energy_occ = mmc_energy_occ +                  &
                                      mmc_depth (ic,MC_OCC, 0, 0) * ival1
                  ENDIF 
!              
!              DO jjs = 0, cr_nscat 
!              IF (mmc_pair (ic, MC_OCC, is, jjs) == -1 ) then 
!                 in = 0 
!                 IF (check_select_status (.true., cr_prop (iatom (in,  &
!                 icent) ), cr_sel_prop) ) then                         
!                    js = cr_iscat (iatom (in, icent) ) 
!                    ival1 = 0 
!                    IF (is.eq.js) then 
!                       ival1 = 1 
!                    ELSEIF (mmc_pair (ic, MC_OCC, is, js) == -1 ) then 
!                       ival1 = - 1 
!                    ENDIF 
!                    mmc_energy_occ = mmc_energy_occ + mmc_depth (ic,   &
!                    MC_OCC, 0, 0) * ival1                              
!                                                                       
!                    ncalc = ncalc + 1 
!                                                                       
!                 ENDIF 
!              ENDIF 
!              ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
      valid_e = .true. 
!                                                                       
      END FUNCTION mmc_energy_occ                   
!*****7*****************************************************************
      REAL function mmc_energy_occ_mol (ianz, imol, amol, valid_e) 
!+                                                                      
!     Calculates the energy for occupational disorder for               
!     molecules.                                                        
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_menu
      USE mc_mod 
      USE mmc_mod 
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
      LOGICAL valid_e 
!                                                                       
      INTEGER ineig (0:maxmol), nneig 
      INTEGER ic, in, ia, ival1, ival2 
!                                                                       
      mmc_energy_occ_mol = 0.0 
      valid_e = .true. 
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
         mmc_energy_occ_mol = mmc_energy_occ_mol + mmc_const (ic,       &
         MC_OCC) * ival1 * ival2                                        
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (.not.valid_e) then 
         mmc_energy_occ_mol = 0.0 
      ENDIF 
!                                                                       
      END FUNCTION mmc_energy_occ_mol               
!*****7*****************************************************************
      REAL function mmc_energy_dis (isel, ia, ic, iatom, icent,  &
      natom, jdir, delta, rdj, valid_e)                                 
!+                                                                      
!     Calculates the energy for chemical disorder                       
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE celltoindex_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
!     USE modify_mod
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
!                                                                       
      INTEGER isel (chem_max_atom) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER i, is, js, ic, in, ia, jjs 
      INTEGER in_a, in_e 
      INTEGER cell (3), site 
      REAL u (3)
!                                                                       
      REAL jdir (3) 
      REAL rdj (CHEM_MAX_COR)
      REAL delta 
      REAL dx 
!                                                                       
!     LOGICAL check_select_status 
!     REAL skalpro 
!                                                                       
      mmc_energy_dis = 0.0 
      valid_e = .false. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .or.  &
          chem_ctyp (ic) == CHEM_ENVIR .or.  &
          chem_ctyp (ic) == CHEM_RANGE .or.  &
          chem_ctyp (ic) == CHEM_DIST  .or.  &
          chem_ctyp (ic) == CHEM_CON        ) then                                               
!                                                                       
         IF (natom (icent) .ne.0) then 
            IF (isel (ia) .eq.iatom (0, icent) ) then 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
               DO jjs = 0, cr_nscat 
               IF (mmc_pair (ic, MC_DISP, is, jjs) == -1 ) then 
                  DO in = in_a, in_e 
                  js = cr_iscat (iatom (in, icent) ) 
                  IF (is.eq.js.or.mmc_pair (ic, MC_DISP, is, js) == -1 ) then 
                     IF (check_select_status (.true., cr_prop (iatom (  &
                     in, icent) ), cr_sel_prop) ) then                  
                        CALL indextocell (iatom (in, icent), cell, site) 
                        DO i = 1, 3 
                        u (i) = cr_pos (i, iatom (in, icent) ) -        &
                        chem_ave_pos (i, site) - float (cell (i)        &
                        - 1) - cr_dim0 (i, 1)                           
                        ENDDO 
                        dx = skalpro (u, jdir, cr_gten) / rdj (ic) 
!                                                                       
                        mmc_energy_dis = mmc_energy_dis + mmc_depth (ic,&
                        MC_DISP, 0, 0) * delta * dx                     
                        valid_e = .true. 
                     ENDIF 
                  ENDIF 
                  ENDDO 
               ENDIF 
               ENDDO 
            ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
               in_a = 0 
               in_e = 0 
               is = cr_iscat (isel (ia) ) 
               DO jjs = 0, cr_nscat 
               IF (mmc_pair (ic, MC_DISP, is, jjs) == -1 ) then 
                  in = 0 
                  js = cr_iscat (iatom (in, icent) ) 
                  IF (is.eq.js.or.mmc_pair (ic, MC_DISP, is, js) == -1 ) then 
                     IF (check_select_status (.true., cr_prop (iatom (  &
                     in, icent) ), cr_sel_prop) ) then                  
                        CALL indextocell (iatom (in, icent), cell, site) 
                        DO i = 1, 3 
                        u (i) = cr_pos (i, iatom (in, icent) ) -        &
                        chem_ave_pos (i, site) - float (cell (i)        &
                        - 1) - cr_dim0 (i, 1)                           
                        ENDDO 
                        dx = skalpro (u, jdir, cr_gten) / rdj (ic) 
!                                                                       
                        mmc_energy_dis = mmc_energy_dis + mmc_depth (ic,&
                        MC_DISP, 0, 0) * delta * dx                     
                        valid_e = .true. 
!                                                                       
                     ENDIF 
                  ENDIF 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION mmc_energy_dis                   
!*****7*****************************************************************
      REAL function mmc_energy_spr (isel, ia, ic, iatom, patom, icent,  &
      natom, valid_e)                                                   
!+                                                                      
!     Calculates the energy for distortions according to                
!                                                                       
!       E = SUM k*(d-d0)**2                                             
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER i, is, js, ic, in, ia 
      INTEGER in_a, in_e 
      INTEGER ncalc 
      REAL patom (3, 0:maxatom, CHEM_MAX_CENT) 
      REAL d, u (3), v (3) 
!                                                                       
!     LOGICAL check_select_status 
!     REAL do_blen 
!                                                                       
      mmc_energy_spr = 0.0 
      ncalc = 0 
      valid_e = .false. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .or.  &
          chem_ctyp (ic) == CHEM_ENVIR .or.  &
          chem_ctyp (ic) == CHEM_RANGE .or.  &
          chem_ctyp (ic) == CHEM_DIST  .or.  &
          chem_ctyp (ic) == CHEM_CON        ) then                                               
!                                                                       
         IF (natom (icent) .ne.0) then 
            IF (isel (ia) .eq.iatom (0, icent) ) then 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               DO i = 1, 3 
               u (i) = patom (i, 0, icent) 
               ENDDO 
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (iatom (in, icent) ) 
            IF (mmc_target_corr (ic, MC_SPRING, is, js) .ne.0.0) then 
               IF (check_select_status (.true., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) then                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.true., u, v) 
!                                                                       
                  mmc_energy_spr = mmc_energy_spr + mmc_depth (ic,      &
                  MC_SPRING, is, js) * (d-mmc_target_corr (ic,          &
                  MC_SPRING, is, js) ) **2                              
                  ncalc = ncalc + 1 
!                                                                       
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (ncalc.gt.0) then 
         mmc_energy_spr = mmc_energy_spr / float (ncalc) 
         valid_e = .true. 
      ELSE 
         mmc_energy_spr = 0.0 
         valid_e = .false. 
      ENDIF 
!                                                                       
      END FUNCTION mmc_energy_spr                   
!*****7*****************************************************************
      REAL function mmc_energy_spr_mol (nmol, imol, valid_e) 
!+                                                                      
!     Calculates the energy for distortions for molecules ..            
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_menu
      USE celltoindex_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE molecule_mod 
      USE rmc_mod 
!     USE modify_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxmol 
      PARAMETER (maxmol = chem_max_neig) 
!                                                                       
      INTEGER imol (2), nmol 
      LOGICAL valid_e 
!                                                                       
      INTEGER ineig (0:maxmol), nneig 
      INTEGER iatom, jatom 
      INTEGER cell (3), site 
      INTEGER i, is, js, ic, in, ia 
      REAL d, u (3), v (3) 
!                                                                       
!     REAL do_blen 
!                                                                       
      mmc_energy_spr_mol = 0.0 
      valid_e = .false. 
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
      mmc_energy_spr_mol = mmc_energy_spr_mol + mmc_const (0, MC_SPRING)&
      * d**2                                                            
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
      IF (chem_ctyp (ic) .eq.CHEM_VEC.or.chem_ctyp (ic) .eq.CHEM_ENVIR) &
      then                                                              
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
         CALL chem_neighbour_mol (imol (ia), ic, ineig, nneig, maxmol) 
         IF (nneig.ne.0) then 
            DO in = 1, nneig 
            is = mole_type (imol (ia) ) 
            js = mole_type (ineig (in) ) 
            jatom = mole_cont (mole_off (ineig (in) ) + 1) 
            IF (mmc_target_corr (ic, MC_SPRING, is, js) .ne.0.0) then 
               DO i = 1, 3 
               u (i) = cr_pos (i, iatom) 
               v (i) = cr_pos (i, jatom) 
               ENDDO 
               d = do_blen (.true., u, v) 
!                                                                       
               mmc_energy_spr_mol = mmc_energy_spr_mol + mmc_const (ic, &
               MC_SPRING) * (d-mmc_target_corr (ic, MC_SPRING, is, js) )&
               **2                                                      
!                                                                       
               valid_e = .true. 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      IF (.not.valid_e) then 
         mmc_energy_spr_mol = 0.0 
      ENDIF 
!                                                                       
      END FUNCTION mmc_energy_spr_mol               
!*****7*****************************************************************
      REAL function mmc_energy_len (isel, ia, ic, iatom, patom, icent,  &
      natom, valid_e)
!+                                                                      
!     Calculates the energy for distortions according to a Lennard-Jones
!     potential                                                         
!                                                                       
!       E = SUM A/d**12 - b/d**6                                        
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER i, is, js, ic, in, ia 
      INTEGER in_a, in_e 
      INTEGER ncalc 
      REAL patom (3, 0:maxatom, CHEM_MAX_CENT) 
      REAL d, u (3), v (3) 
!                                                                       
!     LOGICAL check_select_status 
!     REAL do_blen 
!                                                                       
      mmc_energy_len = 0.0 
      ncalc = 0 
      valid_e = .false. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .or.  &
          chem_ctyp (ic) == CHEM_ENVIR .or.  &
          chem_ctyp (ic) == CHEM_RANGE .or.  &
          chem_ctyp (ic) == CHEM_DIST  .or.  &
          chem_ctyp (ic) == CHEM_CON        ) then                                               
!                                                                       
         IF (natom (icent) .ne.0) then 
            IF (isel (ia) .eq.iatom (0, icent) ) then 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               DO i = 1, 3 
               u (i) = patom (i, 0, icent) 
               ENDDO 
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (iatom (in, icent) ) 
            IF (mmc_target_corr (ic, MC_LENNARD, is, js) .ne.0.0) then 
               IF (check_select_status (.true., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) then                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.true., u, v) 
!                                                                       
                  mmc_energy_len = mmc_energy_len + mmc_depth (ic,      &
                  MC_LENNARD, is, js) * (mmc_len_a (ic, is, js) / d**   &
                  mmc_len_m (ic, is, js) - mmc_len_b (ic, is, js)       &
                  / d**mmc_len_n (ic, is, js) )                         
                  ncalc = ncalc + 1 
!                                                                       
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (ncalc.gt.0) then 
         mmc_energy_len = mmc_energy_len / float (ncalc) 
         valid_e = .true. 
      ELSE 
         mmc_energy_len = 0.0 
         valid_e = .false. 
      ENDIF 
!                                                                       
      END FUNCTION mmc_energy_len                   
!*****7*****************************************************************
      REAL function mmc_energy_rep (isel, ia, ic, iatom, patom, icent,  &
      natom, valid_e)
!+                                                                      
!     Calculates the energy for distortions according to a Repulsive
!     potential                                                         
!                                                                       
!       E = 1/d**n + depth
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER i, is, js, ic, in, ia 
      INTEGER in_a, in_e 
      INTEGER ncalc 
      REAL patom (3, 0:maxatom, CHEM_MAX_CENT) 
      REAL d, u (3), v (3) 
!                                                                       
!     LOGICAL check_select_status 
!     REAL do_blen 
!                                                                       
      mmc_energy_rep = 0.0 
      ncalc = 0 
      valid_e = .false. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .or.  &
          chem_ctyp (ic) == CHEM_ENVIR .or.  &
          chem_ctyp (ic) == CHEM_RANGE .or.  &
          chem_ctyp (ic) == CHEM_DIST  .or.  &
          chem_ctyp (ic) == CHEM_CON        ) then                                               
!                                                                       
         IF (natom (icent) .ne.0) then 
            IF (isel (ia) .eq.iatom (0, icent) ) then 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               DO i = 1, 3 
               u (i) = patom (i, 0, icent) 
               ENDDO 
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (iatom (in, icent) ) 
            IF (mmc_target_corr (ic, MC_REPULSIVE, is, js) .ne.0.0) then 
               IF (check_select_status (.true., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) then                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.true., u, v) 
!                                                                       
                  IF(d.gt.mmc_rep_c (ic, is,js)) THEN
                  mmc_energy_rep = mmc_energy_rep +                     &
                     (-1.)*(ABS(mmc_rep_a (ic, is,js))) +               &
                     1./((d-mmc_rep_c (ic, is,js))/                     &
                          mmc_rep_b (ic, is,js)    )                    &
                                                **mmc_rep_m (ic, is, js)
                  ELSE
                  mmc_energy_rep = mmc_energy_rep +                     &
                      (-1.)*(ABS(mmc_rep_a (ic, is,js))) +              &
                     1./(0.000001                  )                    &
                                                **mmc_rep_m (ic, is, js)
                  ENDIF
                  ncalc = ncalc + 1 
!                                                                       
               ENDIF 
            ELSE
               mmc_energy_rep = mmc_energy_rep - ABS(mmc_rep_low)
               ncalc = ncalc + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (ncalc.gt.0) then 
         mmc_energy_rep = mmc_energy_rep / float (ncalc) 
         valid_e = .true. 
      ELSE 
         mmc_energy_rep =    -1.*ABS(mmc_rep_low)
         valid_e = .true. 
      ENDIF 
!write(*,*) ' Ncalc, energy_rep',ncalc, mmc_energy_rep
!                                                                       
      END FUNCTION mmc_energy_rep
!*****7*****************************************************************
      REAL function mmc_energy_buck (isel, ia, ic, iatom, patom, icent, &
      natom, valid_e)                                                   
!+                                                                      
!     Calculates the energy for distortions according to a Buckingham   
!     potential                                                         
!                                                                       
!       E = SUM A * exp(-d/rho) - B/r**6                                
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER isel (chem_max_atom) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER i, is, js, ic, in, ia 
      INTEGER in_a, in_e 
      INTEGER ncalc 
      REAL patom (3, 0:maxatom, CHEM_MAX_CENT) 
      REAL d, u (3), v (3) 
!                                                                       
!     LOGICAL check_select_status 
!     REAL do_blen 
!                                                                       
      mmc_energy_buck = 0.0 
      ncalc = 0 
      valid_e = .false. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .or.  &
          chem_ctyp (ic) == CHEM_ENVIR .or.  &
          chem_ctyp (ic) == CHEM_RANGE .or.  &
          chem_ctyp (ic) == CHEM_DIST  .or.  &
          chem_ctyp (ic) == CHEM_CON        ) then                                               
!                                                                       
         IF (natom (icent) .ne.0) then 
            IF (isel (ia) .eq.iatom (0, icent) ) then 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               DO i = 1, 3 
               u (i) = patom (i, 0, icent) 
               ENDDO 
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (iatom (in, icent) ) 
            IF (mmc_target_corr (ic, MC_BUCKING, is, js) .ne.0.0) then 
               IF (check_select_status (.true., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) then                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.true., u, v) 
!                                                                       
                  IF (d.gt.mmc_buck_rmin (ic, is, js) ) then 
                     mmc_energy_buck = mmc_energy_buck + mmc_depth (ic, &
                     MC_BUCKING, is, js) * (mmc_buck_a (ic, is, js)     &
                     * exp ( - d / mmc_buck_rho (ic, is, js) ) -        &
                     mmc_buck_b (ic, is, js) / d**6)                    
                  ELSE 
                     mmc_energy_buck = mmc_energy_buck + mmc_buck_atmin &
                     (ic, is, js)                                       
                  ENDIF 
                  ncalc = ncalc + 1 
!                                                                       
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (ncalc.gt.0) then 
         mmc_energy_buck = mmc_energy_buck / float (ncalc) 
         valid_e = .true. 
      ELSE 
         mmc_energy_buck = 0.0 
         valid_e = .false. 
      ENDIF 
!                                                                       
      END FUNCTION mmc_energy_buck                  
!*****7*****************************************************************
      REAL function mmc_energy_angle (ic, iatom, patom, icent,    &
      natom, valid_e)                                                   
!+                                                                      
!     Calculates the energy for angular deviations according to         
!                                                                       
!       E = SUM k*(a-a0)**2                                             
!                                                                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
      USE modify_func_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = chem_max_neig) 
!                                                                       
      INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
      INTEGER icent 
      LOGICAL valid_e 
!                                                                       
      INTEGER i, is, js, ic 
      INTEGER ii, jj 
      LOGICAL lnoneig 
      REAL patom (3, 0:maxatom, CHEM_MAX_CENT) 
      REAL a, b, u (3), v (3), w (3) 
!                                                                       
!     LOGICAL check_select_status 
!     REAL do_bang 
!                                                                       
      lnoneig = .true. 
      mmc_energy_angle = 0.0 
      valid_e = .false. 
!                                                                       
!                                                                       
      IF (chem_ctyp (ic) == CHEM_ANG   .or.  &
          chem_ctyp (ic) == CHEM_VEC   .or.  &
          chem_ctyp (ic) == CHEM_CON   .or.  &
          chem_ctyp (ic) == CHEM_ENVIR      ) then                                               
         IF (natom (icent) .ne.0) then 
            DO i = 1, 3 
            u (i) = patom (i, 0, icent) 
            ENDDO 
            DO ii = 1, natom (icent) - 1 
            DO jj = ii + 1, natom (icent) 
            lnoneig = .false. 
            is = cr_iscat (iatom (ii, icent) ) 
            js = cr_iscat (iatom (jj, icent) ) 
            IF (mmc_target_corr (ic, MC_ANGLE, is, js) .ne.0.0) then 
               IF (check_select_status (.true., cr_prop (iatom (ii,     &
               icent) ), cr_sel_prop) ) then                            
                  IF (check_select_status (.true., cr_prop (iatom (jj,  &
                  icent) ), cr_sel_prop) ) then                         
                     DO i = 1, 3 
                     v (i) = patom (i, ii, icent) 
                     w (i) = patom (i, jj, icent) 
                     ENDDO 
                     a = do_bang (.true., v, u, w) 
                     b = mmc_target_corr (ic, MC_ANGLE, is, js) 
!                                                                       
                     IF (b.le.90) then 
                        IF (a.gt.b) then 
                           mmc_energy_angle = mmc_energy_angle+         &
                           mmc_depth (ic, MC_ANGLE, is, js) * (mod (a + &
                           b / 2., b) - b / 2.) **2                     
                        ELSE 
                           mmc_energy_angle = mmc_energy_angle+         &
                           mmc_depth (ic, MC_ANGLE, is, js) * (a - b) **&
                           2                                            
                        ENDIF 
                     ELSE 
                        mmc_energy_angle = mmc_energy_angle+mmc_depth ( &
                        ic, MC_ANGLE, is, js) * (a - b) **2             
                     ENDIF 
!                                                                       
                     valid_e = .true. 
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (lnoneig) then 
         valid_e = .false. 
      ELSE 
      ENDIF 
!                                                                       
      IF (.not.valid_e) then 
         mmc_energy_angle = 0.0 
      ENDIF 
      END FUNCTION mmc_energy_angle                 
!*****7*****************************************************************
!     REAL function mmc_energy_vec (isel, ia, ic, iatom, patom, icent,  &
!     natom, valid_e)                                                   
!+                                                                      
!     Calculates the energy for vector distributions according to       
!                                                                       
!       E = SUM k*|(v-v0)|**2                                           
!                                                                       
!-                                                                      
!     USE discus_config_mod 
!     USE crystal_mod 
!     USE chem_mod 
!     USE mc_mod 
!     USE mmc_mod 
!     IMPLICIT none 
!                                                                       
       
!     USE debug_mod 
!                                                                       
!     INTEGER maxatom 
!     PARAMETER (maxatom = chem_max_neig) 
!                                                                       
!     INTEGER isel (chem_max_atom) 
!                                                                       
!     INTEGER iatom (0:maxatom, CHEM_MAX_CENT) 
!     INTEGER natom (CHEM_MAX_CENT) 
!     INTEGER icent 
!     LOGICAL valid_e 
!     INTEGER i, is, js, ic, in, ia, nv 
!     INTEGER in_a, in_e 
!     LOGICAL found 
!     REAL patom (3, 0:maxatom, CHEM_MAX_CENT) 
!     REAL an, dd, d, u (3), v (3), w (3) 
!     REAL ddmin 
!     REAL null (3) 
!                                                                       
!     REAL do_blen 
!     REAL do_bang 
!                                                                       
!     DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
!     dbg = .false. 
!     IF (dbg) then 
!     WRITE ( * ,  * ) ' ================energy_vec', '=================&
!    &==========='                                                      
!     ENDIF 
!     mmc_energy_vec = 0.0 
!     valid_e = .false. 
!     ddmin = 1.1 * chem_rmax_env (ic) 
!     found = .false. 
!                                                                       
!     IF (chem_ctyp (ic) .eq.CHEM_VEC.or.chem_ctyp (ic) .eq.CHEM_ENVIR) &
!     then                                                              
!                                                                       
!DBG                                                                    
!        IF (dbg) then 
!           WRITE ( * , * ) ' selected atom,typ ', isel (ia) , cr_iscat &
!           (isel (ia) )                                                
!     WRITE ( * ,  * ) ' selected pos      ',  (cr_pos (i, isel (ia) ) ,&
!    & i = 1, 3)                                                        
!     WRITE ( * ,  * ) ' icent,natom       ', icent, natom (icent) 
!     WRITE ( * ,  * ) ' iatom             ',  (iatom (i, icent) , i = 0&
!    &, natom (icent) )                                                 
!           DO nv = 0, natom (icent) 
!     WRITE ( * ,  * ) ' patom             ',  (patom (i, nv, icent) , i&
!    & = 1, 3)                                                          
!           ENDDO 
!        ENDIF 
!        IF (natom (icent) .ne.0) then 
!           IF (isel (ia) .eq.iatom (0, icent) ) then 
!              IF (dbg) then 
!                 WRITE ( * , * ) ' selected = central' 
!              ENDIF 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
!              DO i = 1, 3 
!              u (i) = patom (i, 0, icent) 
!              ENDDO 
!              in_a = 1 
!              in_e = natom (icent) 
!              is = cr_iscat (iatom (0, icent) ) 
!           ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
!              IF (dbg) then 
!                 WRITE ( * , * ) ' selected = neighbour' 
!              ENDIF 
!              DO i = 1, 3 
!              u (i) = cr_pos (i, isel (ia) ) 
!              ENDDO 
!              in_a = 0 
!              in_e = 0 
!              is = cr_iscat (isel (ia) ) 
!           ENDIF 
!           IF (dbg) then 
!     WRITE ( * ,  * ) ' central           ', u 
!           ENDIF 
!           DO in = in_a, in_e 
!           IF (dbg) then 
!     WRITE ( * ,  * ) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' 
!           ENDIF 
!           js = cr_iscat (iatom (in, icent) ) 
!           DO nv = 1, mmc_nvec (ic, is, js) 
!           IF (dbg) then 
!     WRITE ( * ,  * ) ' is,js ,mmc_vec    ', is, js,  (mmc_vec (i, nv, &
!    &ic, is, js) , i = 1, 3)                                           
!           ENDIF 
!           DO i = 1, 3 
!           v (i) = patom (i, in, icent) - u (i) 
!           w (i) = mmc_vec (i, nv, ic, is, js) 
!           ENDDO 
!           IF (dbg) then 
!     WRITE ( * ,  * ) ' nachbar           ',  (patom (i, in, icent) , i&
!    & = 1, 3)                                                          
!           ENDIF 
!           an = do_bang (.true., v, null, w) 
!           IF (dbg) then 
!     WRITE ( * ,  * ) ' Winkel            ', an 
!     WRITE ( * ,  * ) ' Vektor Zen->Nei   ',  (v (i) , i = 1, 3) 
!              WRITE ( * , * ) ' Differenz zu soll ', (w (i) - v (i) ,  &
!              i = 1, 3)                                                
!           ENDIF 
!           IF (an.lt.mmc_vec (4, nv, ic, is, js) ) then 
!              dd = do_blen (.true., null, v) 
!              IF (dbg) then 
!     WRITE ( * ,  * ) ' Laenge            ', dd 
!              ENDIF 
!              IF (dd.lt.ddmin) then 
!                 ddmin = dd 
!                 d = do_blen (.true., v, w) 
!DBG          mmc_energy_vec = mmc_energy_vec +                         
!DBG     &                                mmc_const(ic,MC_VECTOR) * d**2
!                 found = .true. 
!                 IF (dbg) then 
!     WRITE ( * ,  * ) ' distance           ', d 
!                 ENDIF 
!              ENDIF 
!           ELSEIF (180. - an.lt.mmc_vec (4, nv, ic, is, js) ) then 
!              DO i = 1, 3 
!              v (i) = - v (i) 
!              ENDDO 
!              dd = do_blen (.true., null, v) 
!              IF (dbg) then 
!     WRITE ( * ,  * ) ' Vektor Zen->Nei   ',  (v (i) , i = 1, 3) 
!                 WRITE ( * , * ) ' Differenz zu soll ', (w (i) - v (i) &
!                 , i = 1, 3)                                           
!     WRITE ( * ,  * ) ' Laenge            ', dd 
!              ENDIF 
!              IF (dd.lt.ddmin) then 
!                 ddmin = dd 
!                 d = do_blen (.true., v, w) 
!DBG          mmc_energy_vec = mmc_energy_vec +                         
!DBG     &                                mmc_const(ic,MC_VECTOR) * d**2
!                 found = .true. 
!                 IF (dbg) then 
!     WRITE ( * ,  * ) ' distance           ', d 
!                 ENDIF 
!              ENDIF 
!           ENDIF 
!           WRITE (88, * ) ' nv ', nv 
!           ENDDO 
!           mmc_energy_vec = mmc_energy_vec + 1.0 * d**2 
!DBG     &                                mmc_const(ic,MC_VECTOR) * d**2
!           ENDDO 
!        ENDIF 
!        IF (found) then 
!DBG          mmc_energy_vec = mmc_const(ic,MC_VECTOR) * d**2           
!           valid_e = .true. 
!           CONTINUE 
!        ELSE 
!           mmc_energy_vec = 0.0 
!           valid_e = .false. 
!        ENDIF 
!        IF (dbg) then 
!           WRITE ( * , * ) 'energy_vec ', mmc_energy_vec 
!        ENDIF 
!     ENDIF 
!     dbg = .false. 
!                                                                       
!     IF (.not.valid_e) then 
!        mmc_energy_vec = 0.0 
!     ENDIF 
!     END FUNCTION mmc_energy_vec                   
!*****7*****************************************************************
      SUBROUTINE mmc_correlations (lout, rel_cycl) 
!-                                                                      
!     Determines the achieved correlations                              
!                                                                       
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_menu
      USE chem_aver_mod
      USE celltoindex_mod
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
!
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
       
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL , INTENT(IN) :: lout     ! Flag for output yes/no
      REAL    , INTENT(IN) :: rel_cycl ! Relative progress along cycles
! 
!                                                                       
      CHARACTER(30) energy_name (0:MC_N_ENERGY) 
!                                                                       
      INTEGER maxatom 
      PARAMETER (maxatom = CHEM_MAX_NEIG) 
!                                                                       
      INTEGER ic, je, ic_a 
      INTEGER is, js, ls 
      INTEGER iis, jjs, lls, iic, kk 
      INTEGER i, j, k, l 
      INTEGER icent 
!                                                                       
      LOGICAL searching 
      LOGICAL   :: lfirst = .true.  ! Flag to write output only at first instance
!                                                                       
      INTEGER ncent 
      REAL patom (3, 0:CHEM_MAX_NEIG, CHEM_MAX_CENT) 
      INTEGER iatom (0:CHEM_MAX_NEIG, CHEM_MAX_CENT) 
      INTEGER natom (CHEM_MAX_CENT) 
!                                                                       
      INTEGER bl_anz (0:DEF_maxscat, 0:DEF_maxscat) 
      REAL bl_sum (0:DEF_maxscat, 0:DEF_maxscat) 
      REAL bl_s2 (0:DEF_maxscat, 0:DEF_maxscat) 
      REAL u (3), v (3), d (3) 
      REAL dist 
!                                                                       
      INTEGER pneig (0:DEF_MAXSCAT, 0:DEF_MAXSCAT) 
      INTEGER pair11, pair12, pair21, pair22 
      INTEGER nneigh 
      REAL :: prob11=0.0, prob12, prob22 
      REAL :: thet = 0.0
      REAL :: damp = 1.0
!                                                                       
      REAL wi, wis 
!                                                                       
      REAL ba_sum (CHEM_MAX_COR * MMC_MAX_ANGLES) 
      REAL ba_s2 (CHEM_MAX_COR * MMC_MAX_ANGLES) 
      INTEGER ba_anz (CHEM_MAX_COR * MMC_MAX_ANGLES) 
!                                                                       
      INTEGER icc (3), jcc (3) 
      REAL idir (3), jdir (3), disi (3), disj (3) 
      REAL :: rdi=1.0, rdj=1.0, dpi=1.0, dpj
      INTEGER xnn (0:maxscat, 0:maxscat) 
      REAL xij (0:maxscat, 0:maxscat) 
      REAL xi2 (0:maxscat, 0:maxscat) 
      REAL xj2 (0:maxscat, 0:maxscat) 
!                                                                       
!     REAL do_blen 
!     REAL do_bang 
!     INTEGER angles2index 
!     REAL skalpro 
!                                                                       
      DATA energy_name / 'none', 'Chemical correlation    ', 'Displaceme&
     &nt correlation', 'Distance correlation (Hooke)', 'Angular      cor&
     &relation', 'Vector       correlation', 'Distance     correlation',&
     & 'Lennard Jones potential ', 'Buckingham    potential ' ,         &
       'Repulsive     potential '/         
!
      damp = 0.01 + 0.99*exp(-4.0*rel_cycl)
!                                                                       
!------ Write title line                                                
!                                                                       
      IF (lout) then 
         WRITE (output_io, 410) 
      ENDIF 
!                                                                       
!     Get the average structure for the distance energies               
!                                                                       
      IF (mmc_cor_energy (0, MC_DISP)    .or.mmc_cor_energy (0, MC_SPRING) &
      .or.mmc_cor_energy (0, MC_LENNARD) .or.mmc_cor_energy (0, MC_BUCKING)&
      .or.mmc_cor_energy (0,MC_REPULSIVE) ) then                                                
         CALL chem_aver (.false., .true.) 
      ENDIF 
!                                                                       
!     Reset all achieved correlations                                   
!                                                                       
!
      DO ic = 1, CHEM_MAX_COR 
         DO je = 1, MC_N_ENERGY 
            DO is = - 1, MAXSCAT 
               DO js = - 1, MAXSCAT 
                  mmc_ach_corr (ic, je, is, js) = 0.0 
                  mmc_ach_sigm (ic, je, is, js) = 0.0 
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
      DO i = 1, CHEM_MAX_COR * MMC_MAX_ANGLES 
         mmc_ach_angl (i) = 0.0 
         mmc_ang_sigm (i) = 0.0 
      ENDDO 
!                                                                       
!     Loop over all correlations                                        
!                                                                       
main_corr: DO ic = 1, chem_ncor 
         IF (lout) then 
            WRITE (output_io, * ) 
         ENDIF 
         DO is = 0, MAXSCAT 
            DO js = 0, MAXSCAT 
               pneig (is, js) = 0 
               bl_sum (is, js) = 0.0 
               bl_s2 (is, js) = 0.0 
               bl_anz (is, js) = 0 
               xnn (is, js) = 0 
               xij (is, js) = 0 
               xi2 (is, js) = 0 
               xj2 (is, js) = 0 
            ENDDO 
         ENDDO 
         DO i = 1, CHEM_MAX_COR * MMC_MAX_ANGLES 
            ba_sum (i) = 0.0 
            ba_s2 (i) = 0.0 
            ba_anz (i) = 0 
         ENDDO 
         IF (mmc_cor_energy (0, MC_DISP) ) then 
            DO i = 1, 3 
            idir (i) = chem_dir (i, 1, ic) 
            jdir (i) = chem_dir (i, 2, ic) 
            ENDDO 
!                                                                       
!------ calculate correlations                                          
!                                                                       
            rdi = skalpro (idir, idir, cr_gten) 
            rdj = skalpro (jdir, jdir, cr_gten) 
            IF (rdi.gt.0.0) rdi = sqrt (rdi) 
            IF (rdj.gt.0.0) rdj = sqrt (rdj) 
         ENDIF 
!                                                                       
!     -- Loop over all atoms                                            
!                                                                       
main_atoms:         DO i = 1, cr_natoms 
            is = cr_iscat (i) 
            CALL chem_neighbour_multi (i, ic, iatom, patom, natom, ncent,     &
            maxatom)                                                          
            IF (ier_num.ne.0) return 
!                                                                       
!------ ---- In case of Displacement correlation, calculate             
!            displacement of central atom                                                 
!                                                                       
is_mc_disp: IF (mmc_cor_energy (ic, MC_DISP) ) then 
               CALL indextocell (i, icc, is) 
               DO j = 1, 3 
                  disi (j) = cr_pos (j, i) - chem_ave_pos (j, i) - &
                             float (icc (j) - 1) - cr_dim0 (j, 1)                                          
               ENDDO 
!                                                                       
               IF (chem_ldall (ic) ) then 
                  DO j = 1, 3 
                     jdir (j) = disi (j) 
                  ENDDO 
                  rdj = skalpro (jdir, jdir, cr_gten) 
                  IF (rdj.gt.0.0) then 
                     rdj = sqrt (rdj) 
                  ELSE 
                     rdj = 1.0 
                  ENDIF 
                  dpi = 1.0 
               ELSE 
                  dpi = skalpro (disi, idir, cr_gten) / rdi 
               ENDIF 
            ENDIF is_mc_disp
!                                                                       
!     ---- Loop over all centers                                        
!                                                                       
main_cent: DO icent = 1, ncent 
!                                                                       
!     ------- Since the loop is over all atoms, do central atoms only   
!                                                                       
is_cent:      IF (i.eq.iatom (0, icent) ) then 
is_energy:       IF (mmc_cor_energy (ic, MC_OCC)        .or. &
                     mmc_cor_energy (ic, MC_DISP)       .or. &
                     mmc_cor_energy (ic, MC_SPRING)     .or. &
                     mmc_cor_energy (ic, MC_LENNARD)    .or. &
                     mmc_cor_energy (ic, MC_REPULSIVE)  .or. &
                     mmc_cor_energy (ic, MC_BUCKING)         ) then    
!                                                                       
!     ---------- Loop over all neighbours                               
!                                                                       
            DO j = 1, natom (icent) 
            js = cr_iscat (iatom (j, icent) ) 
            DO k = 1, 3 
            u (k) = patom (k, 0, icent) 
            ENDDO 
!                                                                       
!     --------- Accumulate values for all Energies                      
!                                                                       
            IF (mmc_cor_energy (ic, MC_OCC) ) then 
!                                                                       
!     ----------- Chemical correlation, add number of atom pairs        
!                                                                       
               pneig (is, js) = pneig (is, js) + 1 
            ENDIF 
            IF (mmc_cor_energy (ic, MC_SPRING)    .or.     &
                mmc_cor_energy (ic, MC_LENNARD)   .or.     &
                mmc_cor_energy (ic, MC_REPULSIVE) .or.     &
                mmc_cor_energy (ic, MC_BUCKING)     ) then      
               DO k = 1, 3 
               v (k) = patom (k, j, icent) 
               d (k) = v (k) - u (k) 
               ENDDO 
               dist = do_blen (.true., u, v) 
               js   = cr_iscat (iatom (j, icent) ) 
               bl_sum (is, js) = bl_sum (is, js) + dist 
               bl_s2  (is, js) = bl_s2 (is, js) + dist**2 
               bl_anz (is, js) = bl_anz (is, js) + 1 
               pneig (is, js) = pneig (is, js) + 1 
            ENDIF 
            IF (mmc_cor_energy (ic, MC_DISP) ) then 
               CALL indextocell (iatom (j, icent), jcc, js) 
               DO k = 1, 3 
               disj (k) = cr_pos (k, iatom (j, icent) ) - chem_ave_pos (&
               k, js) - float (jcc (k) - 1) - cr_dim0 (k, 1)            
               ENDDO 
               dpj = skalpro (disj, jdir, cr_gten) / rdj 
               xij (is, js) = xij (is, js) + dpi * dpj 
               xi2 (is, js) = xi2 (is, js) + dpi**2 
               xj2 (is, js) = xj2 (is, js) + dpj**2 
               xnn (is, js) = xnn (is, js) + 1 
            ENDIF 
            ENDDO 
!                     j ! Loop over all neighbours                      
         ENDIF is_energy
         IF (mmc_cor_energy (ic, MC_ANGLE) ) then 
!                                                                       
!     ---------- Angular Correlations                                   
!                                                                       
            is = cr_iscat (i) 
            DO k = 1, 3 
            u (k) = patom (k, 0, icent) 
            ENDDO 
!                                                                       
!     ---------- Double loop over all neighbours                        
!                                                                       
            DO j = 1, natom (icent) - 1 
            js = cr_iscat (iatom (j, icent) ) 
            DO k = 1, 3 
            v (k) = patom (k, j, icent) 
            ENDDO 
            DO l = j + 1, natom (icent) 
            ls = cr_iscat (iatom (l, icent) ) 
!                                                                       
!     -------------- Find proper entry in correlation table             
!                                                                       
            k = 0 
            ic_a = 0 
            searching = .true. 
            DO while (searching.and.k.le.mmc_n_angles) 
            k = k + 1 
            CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls,  &
            MAXSCAT)
            IF (iic.eq.ic) then 
               IF (iis.eq.is.or.iis.eq. - 1) then 
                  IF (jjs.eq. - 1.or.jjs.eq.min (js, ls) ) then 
                     IF (lls.eq. - 1.or.lls.eq.max (js, ls) ) then 
                        searching = .false. 
                        ic_a = k 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
!                         ! Find proper entry in correlation table      
            IF (.not.searching) then 
               DO k = 1, 3 
               d (k) = patom (k, l, icent) 
               ENDDO 
               wi = do_bang (.true., v, u, d) 
               wis = mmc_target_angl (ic_a) 
               IF (wis.le.90.) then 
                  IF (wi.gt.1.5 * wis) then 
                     wi = mod (wi + wis / 2., wis) + wis / 2. 
                  ENDIF 
               ENDIF 
               ba_sum (ic_a) = ba_sum (ic_a) + wi 
               ba_s2 (ic_a) = ba_s2 (ic_a) + wi**2 
               ba_anz (ic_a) = ba_anz (ic_a) + 1 
            ENDIF 
            ENDDO 
            ENDDO 
!                     j ! Double loop over neighbours                   
               ENDIF 
               ENDIF is_cent    !       ! center atoms only                                   
            ENDDO main_cent     ! icent ! Loop over centers                               
         ENDDO main_atoms       ! i     ! Loop over all atoms                                   
!                                                                       
!------ -- Summ up all energies, write output                           
!                                                                       
!     ----- Chemical correlation                                        
!                                                                       
         pair11 = 0
         pair12 = 0
         pair21 = 0
         pair22 = 0
         DO is = 0, cr_nscat 
            DO js =  0, cr_nscat 
               IF     (mmc_pair (ic, MC_OCC, is, js) == -1 ) then 
                  pair12 = pair12 + pneig (is,js)
               ELSEIF (mmc_pair (ic, MC_OCC, is, js) == -2 ) then 
                  pair21 = pair21 + pneig (is,js)
               ELSEIF (mmc_pair (ic, MC_OCC, is, js) == +1 ) then 
                  pair11 = pair11 + pneig (is,js)
               ELSEIF (mmc_pair (ic, MC_OCC, is, js) == +2 ) then 
                  pair22 = pair22 + pneig (is,js)
               ENDIF
            ENDDO
         ENDDO
         je = MC_OCC 
!                                                                       
!                                                                       
         nneigh = pair11 + pair12 + pair21 + pair22 
         IF (nneigh.gt.0.) then 
            prob11 =  pair11           / float (nneigh) 
            prob12 = (pair12 + pair21) / float (nneigh) 
            prob22 =  pair22           / float (nneigh) 
            thet = 0.5 * (2.0 * pair11 + pair12 + pair21) / float(nneigh)                                                     
         ENDIF 
         lfirst = .true.
corr_pair: DO is = 0, cr_nscat 
            DO js = is, cr_nscat 
               IF     (mmc_pair (ic, MC_OCC, is, js) /=  0 ) then 
                 IF (thet.ne.0.0.and.thet.ne.1.0) then 
                    mmc_ach_corr (ic, je, is, js) = (prob11 - thet**2) /&
                                                    (thet * (1 - thet) )
                    mmc_ach_corr (ic, je, js, is) = (prob11 - thet**2) /&
                                                    (thet * (1 - thet) )
!               Feedback mechanism                                      
                    mmc_depth (ic, MC_OCC, 0, 0) = mmc_depth (ic, MC_OCC, 0, 0) - &
                    mmc_cfac (ic, MC_OCC) * (mmc_target_corr (ic, MC_OCC, is,js)- &
                                             mmc_ach_corr (ic, MC_OCC, is, js) ) / 2. &
                    *ABS(mmc_target_corr (ic, MC_OCC, is, js)) &
                    * damp
                 ELSE 
                    mmc_ach_corr (ic, je, is, js) = 0.0 
                    mmc_ach_corr (ic, je, js, is) = 0.0 
                 ENDIF 
!                                                                       
                 IF (lout .and. mmc_pair(ic,MC_OCC,is,js) < 0 .and. lfirst) THEN
                    lfirst = .false.
                    WRITE (output_io, 3100) ic, cr_at_lis (is), cr_at_lis (js),         &
                        mmc_target_corr (ic, je, is, js),                               &
                        mmc_ach_corr (ic, je, is, js),                                  &
                        mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic,je, is, js),&
                        nneigh
                 ENDIF 
!                                                                       
               ENDIF 
            ENDDO 
         ENDDO corr_pair
!                                                                       
!     ----- correlation of displacements                                
!                                                                       
disp_pair: DO is = 0, cr_nscat 
         DO js = is, cr_nscat 
         IF (mmc_pair (ic, MC_DISP, is, js) == -1 ) then 
         je = MC_DISP 
         xnn (is, js) = xnn (is, js) + xnn (js, is) 
         xij (is, js) = xij (is, js) + xij (js, is) 
         xi2 (is, js) = xi2 (is, js) + xi2 (js, is) 
         xj2 (is, js) = xj2 (is, js) + xj2 (js, is) 
         IF (xnn (is, js) .ne.0) then 
            xij (is, js) = xij (is, js) / float (xnn (is, js) ) 
            xi2 (is, js) = xi2 (is, js) / float (xnn (is, js) ) 
            xj2 (is, js) = xj2 (is, js) / float (xnn (is, js) ) 
!                                                                       
            IF (xi2 (is, js) .ne.0.and.xj2 (is, js) .ne.0.0) then 
               mmc_ach_corr (ic, je, is, js) = xij (is, js) / sqrt (xi2 &
               (is, js) * xj2 (is, js) )                                
               mmc_ach_corr (ic, je, js, is) = mmc_ach_corr (ic, je, is,&
               js)                                                      
!               Feedback mechanism                                      
               mmc_depth (ic, MC_DISP, 0, 0) = mmc_depth (ic, MC_DISP, 0,0) - &
               mmc_cfac  (ic, MC_DISP) * (mmc_target_corr (ic, MC_DISP, is, js) - &
               mmc_ach_corr (ic, MC_DISP, is, js) ) / 2. &
                    *ABS(mmc_target_corr (ic, MC_DISP, is, js)) &
                    * damp
            ELSE 
               mmc_ach_corr (ic, je, is, js) = 0.0 
               mmc_ach_corr (ic, je, js, is) = 0.0 
            ENDIF 
         ELSE 
            mmc_ach_corr (ic, je, is, js) = 0.0 
            mmc_ach_corr (ic, je, js, is) = 0.0 
         ENDIF 
!                                                                       
         IF (lout) then 
            WRITE (output_io, 3200) ic, cr_at_lis (is), cr_at_lis (js), &
            mmc_target_corr (ic, je, is, js), mmc_ach_corr (ic, je, is, &
            js), mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic,   &
            je, is, js), nneigh                                         
         ENDIF 
         ENDIF 
         ENDDO 
         ENDDO disp_pair
!                                                                       
!                                                                       
!     -- Loop over all atom pairs to do Hooke potential                 
!                                                                       
spri_pair: DO is = 0, cr_nscat 
         DO js = is, cr_nscat 
         IF (mmc_pair (ic, MC_SPRING, is, js) == -1 ) then 
         je = MC_SPRING 
!                                                                       
!     ----- Spring                                                      
!                                                                       
         IF (bl_anz (is, js) .ne.0.or.bl_anz (js, is) .ne.0) then 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum ( &
            js, is) ) / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2 (is, js) + bl_s2 (js,&
            is) ) / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is,  &
            js) - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js) .gt.0) then 
               mmc_ach_sigm (ic, je, is, js) = sqrt (mmc_ach_sigm (ic,  &
               je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
            IF (lout) then 
      WRITE (output_io, 3300) ic, cr_at_lis (is),  cr_at_lis (js),  mmc_&
     &target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  mmc&
     &_ach_sigm (ic, je, is, js),  mmc_target_corr (ic, je, is, js)  - m&
     &mc_ach_corr (ic, je, is, js),  bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
         ENDIF 
         ENDDO 
         ENDDO  spri_pair
!                                                                       
!     -- Loop over all defined angle correlations                       
!                                                                       
angl_pair: IF (mmc_cor_energy (ic, MC_ANGLE) ) then 
         je = MC_ANGLE 
         DO k = 1, mmc_n_angles 
         CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls,     &
         MAXSCAT)
         IF (ba_anz (k) .ne.0) then 
            mmc_ach_angl (k) = (ba_sum (k) ) / (ba_anz (k) ) 
            mmc_ang_sigm (k) = (ba_s2 (k) ) / (ba_anz (k) ) 
            mmc_ang_sigm (k) = (mmc_ang_sigm (k) - (mmc_ach_angl (k) ** &
            2) )                                                        
            IF (mmc_ang_sigm (k) .gt.0) then 
               mmc_ang_sigm (k) = sqrt (mmc_ang_sigm (k) ) 
            ELSE 
               mmc_ang_sigm (k) = 0.0 
            ENDIF 
            mmc_ach_corr (ic, je, jjs, lls) = mmc_ach_angl (k) 
            mmc_ach_sigm (ic, je, jjs, lls) = mmc_ang_sigm (k) 
            IF (lout) then 
               IF (iic.eq.ic) then 
      WRITE (output_io, 3400) ic, cr_at_lis (iis),  cr_at_lis (jjs),  cr&
     &_at_lis (lls),  mmc_target_angl (k),  mmc_ach_angl (k),  mmc_ang_s&
     &igm (k),  mmc_target_angl (k)  - mmc_ach_angl (k),  ba_anz (k)  + &
     &ba_anz (k)                                                        
                  IF (iis.eq. - 1) then 
                     IF (jjs.eq. - 1) then 
                        IF (lls.eq. - 1) then 
                           searching = .false. 
                           ic_a = k 
                        ELSE 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSE 
            IF (lout) then 
               IF (iic.eq.ic) then 
      WRITE (output_io, 3410) ic, cr_at_lis (iis),  cr_at_lis (jjs),  cr&
     &_at_lis (lls),  mmc_target_angl (k),  0                           
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF angl_pair
!                                                                       
!     -- Loop over all atom pairs to do Lennard Jones potential         
!                                                                       
lenn_pair: DO is = 0, cr_nscat 
      DO js = is, cr_nscat 
      IF (mmc_pair (ic, MC_LENNARD, is, js) == -1 ) then 
         je = MC_LENNARD 
!                                                                       
!     ----- Lennard                                                     
!                                                                       
         IF (bl_anz (is, js) .ne.0.or.bl_anz (js, is) .ne.0) then 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum ( &
            js, is) ) / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2 (is, js) + bl_s2 (js,&
            is) ) / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is,  &
            js) - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js) .gt.0) then 
               mmc_ach_sigm (ic, je, is, js) = sqrt (mmc_ach_sigm (ic,  &
               je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
            IF (lout) then 
      WRITE (output_io, 3700) ic, cr_at_lis (is),  cr_at_lis (js),  mmc_&
     &target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  mmc&
     &_ach_sigm (ic, je, is, js),  mmc_target_corr (ic, je, is, js)  - m&
     &mc_ach_corr (ic, je, is, js),  bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO lenn_pair
!                                                                       
!     -- Loop over all atom pairs to do Repulsive     potential         
!                                                                       
repu_pair: DO is = 0, cr_nscat 
      DO js = is, cr_nscat 
      IF (mmc_pair (ic, MC_REPULSIVE, is, js) == -1 ) then 
         je = MC_REPULSIVE 
!                                                                       
!     ----- REPULSIVE                                                     
!                                                                       
         IF (bl_anz (is, js) .ne.0.or.bl_anz (js, is) .ne.0) then 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum (js, is) ) &
                                          / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2  (is, js) + bl_s2  (js, is) ) &
                                          / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is, js)       &
                                          - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js) .gt.0) then 
               mmc_ach_sigm (ic, je, is, js) = sqrt(mmc_ach_sigm(ic, je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
            IF (lout) then 
               WRITE (output_io, 3900) ic, cr_at_lis (is),  cr_at_lis (js),       &
               mmc_target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  &
               mmc_ach_sigm (ic, je, is, js),                                     &
               mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js), &
               bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
         ENDDO repu_pair      
!                                                                       
!     -- Loop over all atom pairs to do Buckingham potential            
!                                                                       
buck_pair: DO is = 0, cr_nscat 
            DO js = is, cr_nscat 
               IF (mmc_pair (ic, MC_BUCKING, is, js) == -1 ) then 
                  je = MC_BUCKING 
!                                                                       
!     ----- Buckingham                                                  
!                                                                       
                  IF (bl_anz (is, js) .ne.0.or.bl_anz (js, is) .ne.0) then 
                     mmc_ach_corr (ic, je, is, js) =              &
                           (bl_sum (is, js) + bl_sum (js, is) ) / &
                           (bl_anz (is, js) + bl_anz (js, is) )            
                     mmc_ach_sigm (ic, je, is, js) =              &
                           (bl_s2  (is, js) + bl_s2  (js, is) ) / &
                           (bl_anz (is, js) + bl_anz (js, is) )                
                     mmc_ach_sigm (ic, je, is, js) =              &
                           (mmc_ach_sigm (ic, je, is, js) -       &
                           (mmc_ach_corr (ic, je, is, js) **2) )                 
                     IF (mmc_ach_sigm (ic, je, is, js) .gt.0) then 
                        mmc_ach_sigm (ic, je, is, js) =           &
                            sqrt (mmc_ach_sigm (ic, je, is, js) )                                            
                     ELSE 
                        mmc_ach_sigm (ic, je, is, js) = 0.0 
                     ENDIF 
                     IF (lout) then 
                        WRITE (output_io, 2100) ic, cr_at_lis (is), cr_at_lis (  &
                        js), mmc_ach_corr (ic, je, is, js), mmc_ach_sigm (ic, je,&
                        is, js), bl_anz (is, js) + bl_anz (js, is)               
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDDO 
         ENDDO buck_pair
      ENDDO main_corr
!                                                                       
  410 FORMAT ( 45x,'Correlations/',/                                    &
     &   ' Neig.- Energy-',7x,'Atoms',11x,'Target',2x,'Distance/',4x,   &
     &    'Sigma',5x,'Diff',6x,'Number',/                               &
     &    ' Def.   Type    central  Neighbors',13x,'Angle'              &
     &   ,25x,'of pairs')                                               
 2100 FORMAT (1x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3100 FORMAT (1x,i3,3x,'Occupancy',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,i8)
 3200 FORMAT (1x,i3,3x,'Disp.Cor.',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,i8)                                           
 3300 FORMAT (1x,i3,3x,'Hooke    ',a5,3x,a5,      8x,4(f7.3,3x)         &
     &             ,i8)                                                 
 3400 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,4(f7.3,3x)         &
     &             ,i8)                                                 
 3410 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,  f7.3,33x         &
     &             ,i8)                                                 
 3700 FORMAT (1x,i3,3x,'Lennard  ',a5,3x,a5,      8x,4(f7.3,3x)         &
     &             ,i8)                                                 
 3900 FORMAT (1x,i3,3x,'Repuls.  ',a5,3x,a5,      8x,4(f7.3,3x)         &
     &             ,i8)                                                 
!                                                                       
      END SUBROUTINE mmc_correlations               
!*****7*****************************************************************
      SUBROUTINE check_geometry_multi (isel, natoms, laccept) 
!-                                                                      
!     Check whether geometrical constrains apply to the movement of     
!     the selected atoms                                                
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE metric_mod
      USE mc_mod 
      USE mmc_mod 
!
      USE debug_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER isel (CHEM_MAX_ATOM) 
      INTEGER natoms 
      LOGICAL laccept 
!                                                                       
      INTEGER i, j 
      LOGICAL lspace 
      PARAMETER (lspace = .true.) 
      REAL d 
      REAL u (3), null (3) 
!                                                                       
!     REAL do_blen 
!                                                                       
      DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
      laccept = .true. 
!                                                                       
      IF (mmc_l_constrains) then 
         IF (mmc_constrain_type.eq.MMC_C_XYZ) then 
            DO i = 1, natoms 
            DO j = 1, 3 
            IF (cr_pos (j, isel (i) ) .lt.mmc_c_min (j) .or.mmc_c_max ( &
            j) .lt.cr_pos (j, isel (i) ) ) then                         
               laccept = .false. 
               GOTO 9000 
            ENDIF 
            ENDDO 
            ENDDO 
         ELSEIF (mmc_constrain_type.eq.MMC_C_RADIUS) then 
            DO i = 1, natoms 
            DO j = 1, 3 
            u (j) = cr_pos (j, isel (i) ) 
            ENDDO 
            d = do_blen (lspace, u, null) 
            IF (mmc_c_rad.lt.d) then 
               laccept = .false. 
               GOTO 9000 
            ENDIF 
            ENDDO 
         ENDIF 
      ELSE 
         laccept = .true. 
      ENDIF 
!                                                                       
 9000 CONTINUE 
!                                                                       
      END SUBROUTINE check_geometry_multi           
!*****7*****************************************************************
      SUBROUTINE mmc_limit_selection (isel, natoms) 
!-                                                                      
!     Selects an atom from a limited subset                             
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE mmc_mod 
!     USE modify_mod
      USE errlist_mod 
      USE random_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER isel (CHEM_MAX_ATOM) 
      INTEGER natoms 
!                                                                       
      INTEGER i, j 
      INTEGER icell (3), isite, nsel 
      INTEGER ntrial 
!                                                                       
      REAL ran1 
!                                                                       
      IF (mmc_l_type.eq.MMC_L_CELLS) then 
         DO i = 1, natoms 
         ntrial = mmc_l_extend (1) * mmc_l_extend (2) * mmc_l_extend (3)&
         * cr_ncatoms                                                   
         nsel = int( ran1 (idum) * ntrial + 1 )
         CALL mmc_indextocell (nsel, icell, isite, mmc_l_extend) 
         DO j = 1, 3 
         icell (j) = icell (j) + mmc_l_center (j) - 1 
         ENDDO 
         CALL celltoindex (icell, isite, nsel) 
         isel (i) = nsel 
         ENDDO 
      ELSEIF (mmc_l_type.eq.MMC_L_ATOMS) then 
         DO i = 1, natoms 
         isel (i) = int( mmc_l_lower + ran1 (idum) * (mmc_l_upper -          &
         mmc_l_lower + 1) )
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE mmc_limit_selection            
!*****7*****************************************************************
      SUBROUTINE mmc_indextocell (iatom, icell, isite, icc) 
!-                                                                      
!       calculates in which unit cell on which site the atom <ia> is    
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER iatom, isite, icell (3), icc (3) 
      INTEGER ia 
!                                                                       
      ia = iatom - 1 
!                                                                       
      icell (3) = int (ia / icc (1) / icc (2) / cr_ncatoms) + 1 
      ia = ia - (icell (3) - 1) * icc (1) * icc (2) * cr_ncatoms 
      icell (2) = int (ia / icc (1) / cr_ncatoms) + 1 
      ia = ia - (icell (2) - 1) * icc (1) * cr_ncatoms 
      icell (1) = int (ia / cr_ncatoms) + 1 
      isite = ia - (icell (1) - 1) * cr_ncatoms + 1 
!                                                                       
      END SUBROUTINE mmc_indextocell                
!*****7*****************************************************************
      INTEGER function angles2index (ic, nr, is, js, ls, MAXSCAT)
!-                                                                      
!     Calculates a unique number for the angle triplet                  
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER ic 
      INTEGER nr 
      INTEGER is 
      INTEGER js 
      INTEGER ls 
      INTEGER MAXSCAT 
      INTEGER, PARAMETER :: X_ANGLES =  200  ! This needs work, may not be unique!
      INTEGER, PARAMETER :: X_SCAT   =   50  ! This needs work, may not be unique!
!                                                                       
      INTEGER na 
!                                                                       
      na = MAXSCAT + 2 
      na = X_SCAT
!                                                                       
!     angles2index = (ic - 1) * MMC_MAX_ANGLES * na * na * na +  &
      angles2index = (ic - 1) *       X_ANGLES * na * na * na +  &
                     (nr - 1) * na * na * na +                   &
                     (is + 1) * na * na + (js + 1) * na + ls + 2      
      END FUNCTION angles2index                     
!*****7*****************************************************************
      SUBROUTINE index2angles (in, ic, nr, is, js, ls, MAXSCAT)
!-                                                                      
!     Calculates the correlation number, and the angle triplet from     
!     a unique number that was determined by angles2index               
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: X_ANGLES =  200  ! This needs work, may not be unique!
      INTEGER, PARAMETER :: X_SCAT   =   50  ! This needs work, may not be unique!
!
      INTEGER in 
      INTEGER ic 
      INTEGER nr 
      INTEGER is 
      INTEGER js 
      INTEGER ls 
      INTEGER MAXSCAT 
!                                                                       
      INTEGER i 
      INTEGER na 
!                                                                       
      i = in 
      na = MAXSCAT + 2 
      na = X_SCAT
!                                                                       
      ic = i / (      X_ANGLES * na * na * na) + 1 
      i = i - (ic - 1) *       X_ANGLES * na * na * na 
      nr = i / (na * na * na) + 1 
      i = i - (nr - 1) * na * na * na 
      is = i / (na * na) - 1 
      i = i - (is + 1) * na * na 
      js = i / (na) - 1 
      i = i - (js + 1) * na 
      ls = i - 2 
!                                                                       
      END SUBROUTINE index2angles                   
!*****7**************************************************************** 
      SUBROUTINE find_bucking (werte, MAXW) 
!-                                                                      
!     Finds the minimum of a Buckingham function                        
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER MAXW 
      REAL werte (MAXW) 
!                                                                       
      REAL last 
      REAL next 
      REAL step 
      REAL old 
      REAL new 
      REAL minstep 
!                                                                       
!     REAL buckingham 
!                                                                       
      minstep = 0.0001 
!                                                                       
      last = 10.0 
      step = - 0.5 
      old = buckingham (werte (2), werte (3), werte (4), last) 
      next = last + step 
!                                                                       
      DO while (minstep.lt.abs (step) ) 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      DO while ( (new - old) .lt.0) 
      old = new 
      last = next 
      next = last + step 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      ENDDO 
      old = new 
      last = next 
      step = - 0.5 * step 
      next = last + step 
      ENDDO 
!                                                                       
      werte (5) = last + 2. * step 
      werte (6) = buckingham (werte (2), werte (3), werte (4), werte (5)&
      )                                                                 
!                                                                       
      last = werte (5) 
      old = werte (6) 
      step = werte (5) / 100. 
      next = last - step 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      DO while ( (new - old) .gt.0) 
      old = new 
      last = next 
      next = last - step 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      ENDDO 
!                                                                       
      werte (7) = last 
      werte (8) = old 
!                                                                       
      RETURN 
      END SUBROUTINE find_bucking                   
!*****7**************************************************************** 
      REAL function buckingham (a, rho, b, x) 
!-                                                                      
!     calculates the Buckingham potential at point x                    
!+                                                                      
      IMPLICIT none 
!                                                                       
      REAL a 
      REAL rho 
      REAL b 
      REAL x 
!                                                                       
      buckingham = a * exp ( - x / rho) - b / x**6 
!                                                                       
      RETURN 
      END FUNCTION buckingham                       
END MODULE mmc_menu
