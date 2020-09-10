MODULE mmc_menu
!
CONTAINS
!
!*****7**************************************************************** 
!                                                                       
SUBROUTINE mmc 
!-                                                                      
!     This sublevel includes all commands and functions for the         
!     Monte-Carlo simulations in DISCUS.                                
!+                                                                      
USE crystal_mod 
USE discus_allocate_appl_mod
USE chem_mod
USE discus_init_mod
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE chem_symm_mod
!
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
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
USE param_mod 
USE precision_mod
USE prompt_mod 
USE str_comp_mod
USE sup_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MIN_PARA =  20 ! A command requires at least these no of parameters
INTEGER            :: maxw           ! Array size for cpara, lpara, werte
!                                                                       
CHARACTER(LEN=5)    :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=40)   :: cdummy 
CHARACTER(LEN=PREC_STRING) :: line, zeile 
INTEGER             :: lp, length 
INTEGER             :: indxg, lbef
LOGICAL, PARAMETER  :: lold = .FALSE. 
!
INTEGER             :: n_corr = 1 ! dummy for allocation
INTEGER             :: n_scat = 1 ! dummy for allocation
INTEGER             :: n_site = 1 ! dummy for allocation
!                                                                       
!                                                                       
maxw = MAX(MIN_PARA,MAXSCAT+1)
!
! Basic allocation
!
n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
n_site = MAX(MAXSCAT, MMC_MAX_SITE)
! call alloc_chem ! NEEDS WORK
CALL alloc_mmc ( n_corr, MC_N_ENERGY, n_scat, n_site )
CALL alloc_mmc_move(n_corr, n_scat)
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/mmc' 
!
   10 CONTINUE 
!                                                                       
      CALL no_error 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num == 0) THEN 
         IF (line (1:1)  == ' '.OR.line (1:1)  == '#' .OR.   & 
             line == char(13) .OR. line(1:1) == '!'  ) THEN
            IF(linteractive .OR. lmakro) THEN
               GOTO 10
            ELSE
               RETURN
            ENDIF
         ENDIF
!                                                                       
!------ search for "="                                                  
!                                                                       
indxg = index (line, '=') 
IF (indxg /= 0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
            CALL do_math (line, indxg, length) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1)  == '@') THEN 
            line(1:length-1) = line(2:length)
            length = length - 1
            CALL file_kdo(line, length)
!                                                                       
!------ continues a macro 'continue'                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN 
            CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
            CALL do_eval (zeile, lp,.TRUE.) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) THEN 
            GOTO 9999 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
              str_comp (befehl, '?   ', 1, lbef, 4) ) THEN                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 11 
               CALL do_hel ('discus mmc '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
            cdummy = ' ' 
            IF (zeile /= ' ') THEN 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     Waiting for user input                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
            CALL do_input (zeile, lp) 
!                                                                       
!------ ------------------------------------------------------------    
!------ Here start MC specific commands                                 
!------ ------------------------------------------------------------    
!                                                                       
!------ command 'rese'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'rese', 2, lbef, 3) ) THEN 
            CALL mmc_init 
!                                                                       
!------ command 'run'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) THEN 
            CALL mmc_run_multi (.TRUE.)
!                                                                       
!------ command 'save'                                                  
!                                                                       
!        ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) THEN 
!           CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!           CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
!           IF (ier_num == 0) THEN 
!              WRITE (output_io, 1500) cpara (1)(1:len_str (cpara (1) ))
!              CALL save_struc (cpara (1), lpara (1) ) 
!           ENDIF 
!                                                                       
!------ command 'sel' selecting/deselecting atoms                       
!                                                                       
         ELSEIF (str_comp (befehl, 'sele', 3, lbef, 4) .OR.   &
                 str_comp (befehl, 'dese', 2, lbef, 4) ) THEN
!                                                                       
            CALL atom_select (zeile, lp, 0, MMC_MAX_SCAT, mmc_latom, &
            mmc_lsite, 0, MMC_MAX_SITE,            &
            mmc_sel_atom, lold, str_comp (befehl,  &
            'sele', 3, lbef, 4) )                                       
!                                                                       
!------ command 'set'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) THEN 
            CALL mmc_set (zeile, lp) 
!                                                                       
!------ command 'show'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
            CALL mmc_show 
!                                                                       
!------ command 'grow'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'grow', 2, lbef, 4) ) THEN 
            CALL mmc_grow 
!
!------ command 'symmetry'
!
         ELSEIF (str_comp (befehl, 'apply_symm', 2, lbef, 10) ) THEN 
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
      IF (ier_num /= 0) THEN 
         CALL errlist 
         IF (ier_sta /= ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in mmc menu'
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
      IF(linteractive .OR. lmakro) THEN
         GOTO 10
      ELSE
         RETURN
      ENDIF
!                                                                       
 9999 CONTINUE 
!
prompt = orig_prompt
!                                                                       
END SUBROUTINE mmc                            
!
!*****7*****************************************************************
!
SUBROUTINE mmc_show 
!+                                                                      
!     Show parameters of MMC section                                    
!-                                                                      
USE crystal_mod 
USE atom_name
USE chem_mod 
USE chem_menu
USE discus_allocate_appl_mod
USE mc_mod 
USE mmc_mod 
USE molecule_mod 
USE rmc_mod 
!
USE errlist_mod 
USE precision_mod
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=PREC_STRING) :: zeile 
CHARACTER (LEN=20), DIMENSION(MC_N_MOVE)    :: c_move = & !  (MC_N_MOVE) 
   (/ 'switch chemistry    ', &
      'switch displacement ', &
      'shift atom          ', &
      'inverse displacement' /)
CHARACTER (LEN=20), DIMENSION(5)            :: c_site = & !(4) 
   (/ 'all                 ', &
      'local +- 1 unit cell', &
      'local, same site    ', &
      'all,   same site    ', &
      'connectivity        ' /)
CHARACTER(9) at_name_i, at_name_j 
INTEGER :: i, j, k
INTEGER :: ie        ! Energy for current correlation 
INTEGER :: n_corr
INTEGER :: n_scat
INTEGER :: n_site
!                                                                       
n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
n_site = MAX(MAXSCAT, MMC_MAX_SITE)
! call alloc_chem ! NEEDS WORK
call alloc_mmc ( n_corr, MC_N_ENERGY, n_scat, n_site )
CALL alloc_mmc_move(n_corr, n_scat)
!                                                                       
IF (mo_sel_atom) THEN 
   WRITE (output_io, 1105) 'atoms' 
ELSE 
   WRITE (output_io, 1105) 'molecules' 
ENDIF 
!                                                                       
WRITE (output_io, 1250) mo_cyc 
WRITE (output_io, 1300) mo_feed 
WRITE (output_io, 1350) mmc_no_valid 
WRITE (output_io, 1400) mo_kt 
IF(mmc_h_stop) THEN
   WRITE(output_io, '(3x,a)'     ) 'MMC will stop at convergence'
   WRITE(output_io, '(3x,a,f8.4,a)') '   Max rel. difference       : ', mmc_h_conv_m, &
                                   ' at last cycle'
   WRITE(output_io, '(3x,a,f8.4,a,i4,a)') '   Max change rel. difference: ', mmc_h_conv_c, &
                                   ' observed over ',MMC_H_NNNN ,' feedbacks'
   WRITE(output_io, '(3x,a,f8.4,a,i4,a)') '   Av. change rel. difference: ', mmc_h_conv_a, &
                                   ' averaged over ',MMC_H_NNNN ,' feedbacks'
ELSE
   WRITE(output_io, '(3x,a)'     ) 'MMC will stop at full cycles'
ENDIF
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
IF (mmc_cor_energy (0, MC_OCC) ) THEN
   WRITE (output_io, 2100) 
   IF (mo_sel_atom) THEN 
      DO k = 1, chem_ncor 
         DO ie = 1, 1
            DO i = 0, cr_nscat 
               DO j = i, cr_nscat 
                  at_name_i = at_name (i) 
                  at_name_j = at_name (j) 
                  IF (mmc_pair        (k, ie, i, j) /=  0.0) THEN 
                     WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
                     mmc_target_corr (k, MC_OCC, i, j),                       &
                     mmc_depth (k, MC_OCC,i, j)
                  ENDIF 
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
   ELSE 
      DO k = 1, chem_ncor 
         DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
               IF (mmc_pair        (k, MC_OCC, i, j) /=  0.0) THEN 
                  WRITE (output_io, 4200) i, j, k,         &
                        mmc_target_corr (k,MC_OCC, i, j),  &
                        mmc_depth (k, MC_OCC, i, j)
                  ENDIF 
               ENDDO 
            ENDDO 
         ENDDO 
      ENDIF 
   ENDIF 
IF (mmc_cor_energy (0, MC_UNI)) THEN
   WRITE (output_io, 2100) 
   IF (mo_sel_atom) THEN 
      DO k = 1, chem_ncor 
         DO i = 0, cr_nscat 
            DO j = 0, cr_nscat 
               at_name_i = at_name (i) 
               at_name_j = at_name (j) 
!           IF (mmc_target_corr (k, MC_OCC, i, j)  /= 0.0) THEN 
               IF (mmc_pair        (k, MC_UNI, i, j) /=  0.0) THEN 
                  WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
                  mmc_target_corr (k, MC_UNI, i, j),                       &
                  mmc_depth (k, MC_UNI,i, j)
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 
   ELSE 
      DO k = 1, chem_ncor 
         DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
!           IF (mmc_target_corr (k, MC_OCC, i, j)  /= 0.0) THEN 
               IF (mmc_pair        (k, MC_UNI, i, j) /=  0.0) THEN 
                  WRITE (output_io, 4200) i, j, k,         &
                        mmc_target_corr (k,MC_UNI, i, j),  &
                        mmc_depth (k, MC_UNI, i, j)
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 
   ENDIF 
ENDIF 
!
IF (mmc_cor_energy (0, MC_DISP) ) THEN 
   WRITE (output_io, 2200) 
   IF (mo_sel_atom) THEN 
      DO k = 1, chem_ncor 
         DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
               at_name_i = at_name (i) 
               at_name_j = at_name (j) 
               IF(mmc_pair        (k, MC_DISP, i, j) /=  0.0) THEN 
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
               IF(mmc_pair        (k, MC_DISP, i, j) /=  0.0) THEN 
                  WRITE (output_io, 4200) i, j, k,                   &
                  mmc_target_corr (k,MC_DISP, i, j),                 &
                  mmc_depth (k, MC_DISP, i, j)
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 
   ENDIF 
ENDIF 
!
IF(mmc_cor_energy(0, MC_COORDNUM) ) THEN 
   WRITE (output_io, 2400) 
!                                                                       
   IF (mo_sel_atom) THEN 
      DO k = 1, chem_ncor 
         DO i = 0, cr_nscat 
            DO j = 0, cr_nscat 
               at_name_i = at_name (i) 
               at_name_j = at_name (j) 
               IF (mmc_target_corr (k, MC_COORDNUM, i, j)  /= 0.0) THEN 
                  WRITE (output_io, 7300) at_name_i, at_name_j, k,         &
                  mmc_target_corr (k, MC_COORDNUM, i, j),                    &
                  mmc_depth (k,MC_COORDNUM, i, j)
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 
   ELSE 
      DO k = 1, chem_ncor 
         DO i = 1, mole_num_type 
            DO j = i, mole_num_type 
               IF (mmc_target_corr (k, MC_COORDNUM, i, j)  /= 0.0) THEN 
                  WRITE (output_io, 4200) i, j, k,                         &
                  mmc_target_corr (k,MC_COORDNUM, i, j),                     &
                  mmc_depth (k, MC_COORDNUM, i, j)
               ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 
   ENDIF 
ENDIF 
!
      IF (mmc_cor_energy (0, MC_SPRING) ) THEN 
         WRITE (output_io, 2300) 
!                                                                       
         IF (mo_sel_atom) THEN 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_SPRING, i, j)  /= 0.0) THEN 
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
            IF (mmc_target_corr (k, MC_SPRING, i, j)  /= 0.0) THEN 
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
      IF (mmc_cor_energy (0, MC_ANGLE) ) THEN 
         WRITE (output_io, 2010) 
!                                                                       
         IF (mo_sel_atom) THEN 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_ANGLE, i, j)  /= 0.0) THEN 
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
            IF (mmc_target_corr (k, MC_ANGLE, i, j)  /= 0.0) THEN 
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
!     IF (mmc_cor_energy (0, MC_VECTOR) ) THEN 
!        WRITE (output_io, 2020) 
!                                                                       
!        IF (mo_sel_atom) THEN 
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
      IF (mmc_cor_energy (0, MC_LENNARD) ) THEN 
         WRITE (output_io, 2700) 
         IF (mo_sel_atom) THEN 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_LENNARD, i, j)  /= 0.0) THEN 
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
            IF (mmc_target_corr (k, MC_LENNARD, i, j)  /= 0.0) THEN 
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
      IF (mmc_cor_energy (0, MC_REPULSIVE) ) THEN 
         WRITE (output_io, 2750) 
         IF (mo_sel_atom) THEN 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_REPULSIVE, i, j)  /= 0.0) THEN 
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
            IF (mmc_target_corr (k, MC_REPULSIVE, i, j)  /= 0.0) THEN 
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
      IF (mmc_cor_energy (0, MC_BUCKING) ) THEN 
         WRITE (output_io, 2800) 
         IF (mo_sel_atom) THEN 
            DO k = 1, chem_ncor 
            DO i = 0, cr_nscat 
            DO j = i, cr_nscat 
            at_name_i = at_name (i) 
            at_name_j = at_name (j) 
            IF (mmc_target_corr (k, MC_BUCKING, i, j)  /= 0.0) THEN 
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
            IF (mmc_target_corr (k, MC_BUCKING, i, j)  /= 0.0) THEN 
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
      IF (mmc_move_prob (3)  > 0.0) THEN 
         WRITE (output_io, 6000) 
         IF (mo_sel_atom) THEN 
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
 1350 FORMAT (  '   Maximum no. non_valid cycles : ',i8) 
 1400 FORMAT (  '   Temperature [kT]             : ',f8.4) 
 3000 FORMAT (/,' Correlation definitions        : ',/) 
 2100 FORMAT (/,' Desired correlations for Chemical Occupancy : ',/,/,  &
     &         12x,'Pairs',10x,'neigh. #',3x,'correl. ',7x,'depth')     
 2200 FORMAT (/,' Desired correlations for Displacement : ',/,/,        &
     &         12x,'Pairs',10x,'neigh. #',3x,'correl. ',7x,'depth')     
 2400 FORMAT (/,' Desired Rcoordination numbers  : ',/,/,               &
     &         12x,'Pairs',10x,'neigh. #',3x,'coordno.',7x,'depth')     
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
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set (zeile, lp) 
!+                                                                      
!     sets parameters for MMC section                                   
!-                                                                      
USE discus_allocate_appl_mod
USE crystal_mod 
USE chem_mod 
USE chem_menu
USE get_iscat_mod
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE rmc_sup_mod
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile 
INTEGER         , INTENT(INOUT) :: lp 
!                                                                       
INTEGER, PARAMETER :: maxw = 200 
!
CHARACTER(LEN=PREC_STRING)                 :: line
      INTEGER                               :: length
      INTEGER                               :: ianz1
      INTEGER                               :: ianz2
!                                                                       
      CHARACTER (LEN=PREC_STRING), DIMENSION(MAXW) ::  cpara !(maxw) 
      CHARACTER (LEN=PREC_STRING), DIMENSION(MAXW) ::  cpara1 !(maxw) 
      CHARACTER (LEN=PREC_STRING), DIMENSION(MAXW) ::  cpara2 !(maxw) 
      REAL(KIND=PREC_DP) :: uerte (maxw) 
      REAL(KIND=PREC_DP) :: verte (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte1 (maxw) 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte2 (maxw) 
      REAL(KIND=PREC_DP) :: a, b 
      INTEGER lpara (maxw) 
      INTEGER, DIMENSION(MAXW) :: lpara1 ! (maxw) 
      INTEGER, DIMENSION(MAXW) :: lpara2 ! (maxw) 
      INTEGER ianz, iianz, jjanz, kkanz, is, js, ls, ic, i, j 
      INTEGER is_start, is_end 
      INTEGER                :: n_corr ! Dummy for allocation
      INTEGER                :: n_scat ! Dummy for allocation
      INTEGER                :: n_site ! Dummy for allocation
      INTEGER                :: n_angles ! Dummy for allocation
LOGICAL :: is_corr ! Current target is correlation energy
!                                                                       
!     INTEGER angles2index 
!                                                                       
!                                                                       
      n_corr = 0
      n_scat = 0
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num == 0) THEN 
   IF (ianz >= 2) THEN 
      CALL do_cap (cpara (1) ) 
!
!------ --- 'set allowed', set which atom are allowed in a move
!
      IF(cpara(1)(1:2)  == 'AL') THEN 
         CALL del_params(1, ianz, cpara, lpara, maxw) 
         CALL get_iscat(ianz, cpara, lpara, werte,     &
                        maxw, .FALSE.)                                  
         IF (NINT(werte(1)) == -1) THEN
                  mmc_allowed = .TRUE.
         ELSE
                  DO i=1,ianz
                     mmc_allowed(nint(werte(i))) = .TRUE.
                  ENDDO
         ENDIF
!
!------ --- 'set disallowed', set which atom are allowed in a move
!
      ELSEIF (cpara (1) (1:2)  == 'DI') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL get_iscat (ianz, cpara, lpara, werte,     &
                        maxw, .FALSE.)                                  
               IF (NINT(werte(1)) == -1) THEN
                  mmc_allowed = .TRUE.
               ELSE
                  DO i=1,ianz
                     mmc_allowed(nint(werte(i))) = .FALSE.
                  ENDDO
               ENDIF
!                                                                       
!------ --- 'set angle' : setting of neighbouring angles                
!                                                                       
      ELSEIF (cpara (1) (1:2)  == 'AN') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_angle (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set conn' : setting of neighbouring connectivity           
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'CON') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_con (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set envir' : setting of neighbouring environment           
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'ENV') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_envir (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set cyc' : setting number of MC moves                      
!                                                                       
      ELSEIF (cpara (1) (1:2)  == 'CY') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_cyc = NINT(werte (1), PREC_INT_LARGE ) 
!                                                                       
!------ --- 'set feed' : setting display/feedback intervall             
!                                                                       
      ELSEIF(cpara (1) (1:2)  == 'FE') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_feed = nint (werte (1) ) 
      ELSEIF(str_comp(cpara(1), 'FINISH', 3, lpara(1), 6)) THEN
         CALL del_params(1, ianz, cpara, lpara, maxw)
         CALL get_finish(ianz, cpara, lpara, werte, MAXW)
!                                                                       
!------ --- 'set limited': sets limited selction range for atoms        
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'LIM') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (cpara (1) (1:4)  == 'cell') THEN 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num /= 0) return 
                  DO i = 1, 3 
                  mmc_l_center (i) = int( werte (i) )
                  mmc_l_extend (i) = int( werte (i + 3) )
                  ENDDO 
                  mmc_l_limited = .TRUE. 
                  mmc_l_type = MMC_L_CELLS 
               ELSEIF (cpara (1) (1:4)  == 'atom') THEN 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num /= 0) return 
                  mmc_l_lower = int( werte (1) )
                  mmc_l_upper = int( werte (2) )
                  mmc_l_limited = .TRUE. 
                  mmc_l_type = MMC_L_ATOMS 
               ELSEIF (cpara (1) (1:3)  == 'OFF') THEN 
                  mmc_l_limited = .FALSE. 
               ENDIF 
!                                                                       
!------ --- 'set mode': sets operation mode for MMC                     
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'MOD') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num /= 0) return 
               j = 1 
               CALL ber_params (j, cpara, lpara, werte, maxw) 
               IF (ier_num /= 0) return 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num /= 0) return 
               CALL mmc_set_mode (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set move': sets maxmove for shift MMC mode                 
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'MOV') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num /= 0) return 
               CALL rmc_set_move (mo_maxmove, mo_sel_atom, ianz, cpara, &
               werte, lpara, maxw, 4)                                      
!                                                                       
!------ --- 'set neig': setting correlation determination method        
!                                                                       
      ELSEIF (cpara (1) (1:2)  == 'NE') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_neig (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set range' : setting of neighbouring ranges                
!                                                                       
      ELSEIF (cpara (1) (1:2)  == 'RA') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_ranges (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set target' : setting of target correlations               
!                                                                       
      ELSEIF (cpara (1) (1:2)  == 'TA') THEN 
         IF (ianz >= 3) THEN 
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
            IF( ic > CHEM_MAX_COR .OR. ic > MMC_MAX_CORR ) THEN
!
!                      Basic allocation
!
               n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
               n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
               n_site = MAX(MAXSCAT, MMC_MAX_SITE)
               ! call alloc_chem ! NEEDS WORK
               CALL alloc_mmc ( n_corr, MC_N_ENERGY, n_scat, n_site )
               CALL alloc_mmc_move(n_corr, n_scat)
            ENDIF
            IF (ic > 0.AND.ic <= chem_ncor) THEN 
               IF (str_comp (cpara (2) , 'corr', 2, lpara (2) , 4) .OR.  &
                   str_comp (cpara (2) , 'unid', 2, lpara (2) , 4)) THEN                                             
                  is_corr= str_comp (cpara (2) , 'corr', 2, lpara (2) , 4)
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
!                       Get atom types in the two allowed groups
!                       Allowed parameters are: Au, Cu            ! OLd style 
!                                               (Au), (Cu)        ! New style  :binary correlation
!                                               (Au,Pt), (Cu,Zn)  ! New style  :quaternary correlation
                  IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
                     line   = cpara(1)(2:lpara(1)-1)
                     length = lpara(1)-2
                     CALL get_params (line, ianz1, cpara1, lpara1, maxw, length) 
                     CALL get_iscat (ianz1, cpara1, lpara1, werte1, maxw, .FALSE.)                                  
                  ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
                     ianz1 = 1 
                     CALL get_iscat (ianz1, cpara, lpara, werte1, maxw, .FALSE.)                                  
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN
                  ENDIF
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
                     line   = cpara(1)(2:lpara(1)-1)
                     length = lpara(1)-2
                     CALL get_params (line, ianz2, cpara2, lpara2, maxw, length) 
                     CALL get_iscat (ianz2, cpara2, lpara2, werte2, maxw, .FALSE.)                                  
                  ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
                     ianz2 = 1 
                     CALL get_iscat (ianz2, cpara, lpara, werte2, maxw, .FALSE.)                                  
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN
                  ENDIF
!
                  IF(is_corr                                    ) THEN                                             
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     IF (cpara (ianz) (1:2)  == 'CO') THEN 
                        mmc_cfac (ic, MC_OCC) =  1.0 
                        ianz = ianz - 1 
                     ELSEIF (cpara (ianz) (1:2)  == 'EN') THEN 
                        mmc_cfac (ic, MC_OCC) = 0.0 
                        ianz = ianz - 1 
                     ENDIF 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw)
                     IF (mmc_cfac (ic, MC_OCC) > 0.0) THEN 
                        CALL mmc_set_disp_occ (ic, MC_OCC, ianz1, ianz2, &
                             MAXW, werte1, werte2, werte(1) , -1.00*werte(1) )                                   
                           mmc_depth (ic, MC_OCC, 0, 0) = -1.00*werte (1) 
                     ELSEIF (mmc_cfac (ic, MC_OCC) ==0.0) THEN 
                     CALL mmc_set_disp_occ (ic, MC_OCC, ianz1, ianz2, &
                             MAXW, werte1, werte2, werte(1) , werte(2) )                                   
                           mmc_depth (ic, MC_OCC, 0, 0) = werte (2) 
                     ENDIF
                     mmc_cor_energy (ic, MC_OCC) = .TRUE. 
                     mmc_cor_energy (0, MC_OCC) = .TRUE. 
                  ELSE
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     IF (cpara (ianz) (1:2)  == 'CO') THEN 
                        mmc_cfac (ic, MC_UNI) =  1.0 
                        ianz = ianz - 1 
                     ELSEIF (cpara (ianz) (1:2)  == 'EN') THEN 
                        mmc_cfac (ic, MC_UNI) = 0.0 
                        ianz = ianz - 1 
                     ENDIF 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw)
                     IF (mmc_cfac (ic, MC_UNI) > 0.0) THEN 
                        CALL mmc_set_unid_occ (ic, MC_UNI, ianz1, ianz2, &
                             MAXW, werte1, werte2, werte(1) , -0.50*werte(1) )                                   
                             mmc_depth (ic, MC_UNI, 0, 0) = -0.50*werte (1) 
                     ELSEIF (mmc_cfac (ic, MC_UNI) ==0.0) THEN 
                        CALL mmc_set_unid_occ (ic, MC_UNI, ianz1, ianz2, &
                             MAXW, werte1, werte2, werte(1) , werte(2) )                                   
                             mmc_depth (ic, MC_UNI, 0, 0) = werte (2) 
                     ENDIF
                     mmc_cor_energy (ic, MC_UNI) = .TRUE. 
                     mmc_cor_energy (0, MC_UNI) = .TRUE. 
                  ENDIF
               ELSEIF (str_comp (cpara (2) , 'cd', 2, lpara (2) , &
                     2) ) THEN                                          
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF (cpara (ianz) (1:2)  == 'CO') THEN 
                           mmc_cfac (ic, MC_DISP) =  1.00 
                           ianz = ianz - 1 
                        ELSEIF (cpara (ianz) (1:2)  == 'EN') THEN 
                           mmc_cfac (ic, MC_DISP) = 0.0 
                           ianz = ianz - 1 
                        ENDIF 
                        CALL ber_params (2, cpara, lpara, werte, maxw) 
                        IF (mmc_cfac (ic, MC_DISP) /= 0.0) THEN 
                           werte(2) = -500.*werte(1)
                        ENDIF
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           IF(is<0 .OR. is>cr_nscat .OR. &
                              js<0 .OR. js>cr_nscat)  THEN
                              ier_num = -97
                              ier_typ = ER_APPL
                              WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
                              WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
                              ier_msg(3) = 'Check parameter list for wrong/missing atoms'
                              RETURN
                           ENDIF
                           mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_DISP, is, js, werte ( &
                        1), werte (2) )                                 
                        ENDDO 
                        ENDDO 
                        chem_ldall (ic) = .TRUE. 
                        IF (ianz > 2) THEN 
                           CALL del_params (2, ianz, cpara, lpara, maxw) 
                           IF (cpara (1) (1:3)  == 'all') THEN 
                              chem_ldall (ic) = .TRUE. 
                           ELSE 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num /= 0) return 
!                                                                       
                              IF (ianz == 3) THEN 
                                 chem_ldall (ic) = .FALSE. 
                                 chem_dir (1, 1, ic) = werte (1) 
                                 chem_dir (2, 1, ic) = werte (2) 
                                 chem_dir (3, 1, ic) = werte (3) 
                                 chem_dir (1, 2, ic) = werte (1) 
                                 chem_dir (2, 2, ic) = werte (2) 
                                 chem_dir (3, 2, ic) = werte (3) 
                              ELSEIF (ianz == 6) THEN 
                                 chem_ldall (ic) = .FALSE. 
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
                        mmc_cor_energy (ic, MC_DISP) = .TRUE. 
                        mmc_cor_energy (0, MC_DISP) = .TRUE. 
                     ELSEIF(str_comp(cpara(2) , 'cn', 2, lpara(2), 2) ) THEN   !Coordination number
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
!                       Get atom types in the two allowed groups
!                       Allowed parameters are: Au, Cu            ! OLd style 
!                                               (Au), (Cu)        ! New style  :binary correlation
!                                               (Au,Pt), (Cu,Zn)  ! New style  :quaternary correlation
                        IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
                           line   = cpara(1)(2:lpara(1)-1)
                           length = lpara(1)-2
                           CALL get_params (line, ianz1, cpara1, lpara1, maxw, length) 
                           CALL get_iscat (ianz1, cpara1, lpara1, werte1, maxw, .FALSE.)                                  
                        ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
                           ianz1 = 1 
                           CALL get_iscat (ianz1, cpara, lpara, werte1, maxw, .FALSE.)                                  
                        ELSE
                           ier_num = -6
                           ier_typ = ER_COMM
                           RETURN
                        ENDIF
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
                           line   = cpara(1)(2:lpara(1)-1)
                           length = lpara(1)-2
                           CALL get_params (line, ianz2, cpara2, lpara2, maxw, length) 
                           CALL get_iscat (ianz2, cpara2, lpara2, werte2, maxw, .FALSE.)                                  
                        ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
                           ianz2 = 1 
                           CALL get_iscat (ianz2, cpara, lpara, werte2, maxw, .FALSE.)                                  
                        ELSE
                           ier_num = -6
                           ier_typ = ER_COMM
                           RETURN
                        ENDIF
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (2, cpara, lpara, werte, maxw) 
!                       DO i = 1, iianz 
!                       DO j = 1, jjanz 
!                       is = nint (uerte (i) ) 
!                       js = nint (verte (j) ) 
!                          IF(is<0 .OR. is>cr_nscat .OR. &
!                             js<0 .OR. js>cr_nscat)  THEN
!                             ier_num = -97
!                             ier_typ = ER_APPL
!                             WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
!                             WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
!                             ier_msg(3) = 'Check parameter list for wrong/missing atoms'
!                             RETURN
!                          ENDIF
!                          mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
!                          mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
                           CALL mmc_set_cn(ic, MC_COORDNUM, ianz1, ianz2, MAXW, werte1, werte2, werte(1), werte(2)) 
!                       ENDDO 
!                       ENDDO 
                        mmc_cor_energy (ic, MC_COORDNUM) = .TRUE. 
                        mmc_cor_energy (0, MC_COORDNUM) = .TRUE. 
                        mmc_cfac (ic, MC_COORDNUM) = 1.0 
                     ELSEIF (str_comp (cpara (2) , 'spring', 2, lpara ( &
                     2) , 6) ) THEN                                     
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (2, cpara, lpara, werte, maxw) 
                        DO i = 1, iianz 
                        DO j = 1, jjanz 
                        is = nint (uerte (i) ) 
                        js = nint (verte (j) ) 
                           IF(is<0 .OR. is>cr_nscat .OR. &
                              js<0 .OR. js>cr_nscat)  THEN
                              ier_num = -97
                              ier_typ = ER_APPL
                              WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
                              WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
                              ier_msg(3) = 'Check parameter list for wrong/missing atoms'
                              RETURN
                           ENDIF
                           mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_SPRING, is, js, werte &
                        (1), werte (2) )                                
                        ENDDO 
                        ENDDO 
                        mmc_cor_energy (ic, MC_SPRING) = .TRUE. 
                        mmc_cor_energy (0, MC_SPRING) = .TRUE. 
                     ELSEIF (str_comp (cpara (2) , 'angle', 2, lpara (2)&
                     , 5) ) THEN                                        
!                                                                       
!     ------------Three atom names are listed on the target instruction 
!                 Required for angular correlations                     
!                                                                       
                        IF (mmc_n_angles == 0          .OR.     &
                            mmc_n_angles >= MMC_MAX_ANGLES) THEN 
                           n_angles = max(mmc_n_angles+20, int(MMC_MAX_ANGLES*1.025))
                           CALL alloc_mmc_angle (CHEM_MAX_COR,n_angles)
                        ENDIF
                        IF (mmc_n_angles < MMC_MAX_ANGLES) THEN 
                           CALL del_params (2, ianz, cpara, lpara, maxw) 
                           mmc_n_angles = mmc_n_angles + 1 
                           iianz = 1 
                           jjanz = 1 
                           kkanz = 1 
                           CALL get_iscat (iianz, cpara, lpara, uerte,  &
                           maxw, .FALSE.)                               
                           IF (uerte (1)  ==  - 1) THEN 
                              is_start = 0 
                              is_end = cr_nscat 
                              mmc_allowed = .TRUE.           ! all atoms are allowed in mmc moves
                           ELSE 
                              is_start = int( uerte (1) )
                              is_end = int( uerte (1) )
                              mmc_allowed(is_start) = .TRUE. ! this atom is allowed in mmc moves
                           ENDIF 
                           is = int( uerte (1) )
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL get_iscat (jjanz, cpara, lpara, uerte,  &
                           maxw, .FALSE.)                               
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL get_iscat (kkanz, cpara, lpara, verte,  &
                           maxw, .FALSE.)                               
                           js = min (uerte (1), verte (1) ) 
                           ls = max (uerte (1), verte (1) ) 
                           mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
                           mmc_allowed(ls) = .TRUE. ! this atom is allowed in mmc moves
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
                           mmc_cor_energy (ic, MC_ANGLE) = .TRUE. 
                           mmc_cor_energy (0, MC_ANGLE) = .TRUE. 
                        ENDIF 
!DBG_VEC For later development                                          
!                    ELSEIF (str_comp (cpara (2) , 'vector', 2, lpara ( &
!                    2) , 6) ) THEN                                     
!                       CALL del_params (2, ianz, cpara, lpara, maxw) 
!                       iianz = 1 
!                       jjanz = 1 
!                       CALL get_iscat (iianz, cpara, lpara, uerte,     &
!                       maxw, .FALSE.)                                  
!                       CALL del_params (1, ianz, cpara, lpara, maxw) 
!                       CALL get_iscat (jjanz, cpara, lpara, verte,     &
!                       maxw, .FALSE.)                                  
!                       CALL del_params (1, ianz, cpara, lpara, maxw) 
!                       CALL ber_params (4, cpara, lpara, werte, maxw) 
!                       DO i = 1, iianz 
!                       DO j = 1, jjanz 
!                       is = nint (uerte (i) ) 
!                       js = nint (verte (j) ) 
!                       CALL mmc_set_vec (ic, is, js, werte, maxw) 
!                       ENDDO 
!                       ENDDO 
!                       mmc_cor_energy (ic, MC_VECTOR) = .TRUE. 
!                       mmc_cor_energy (0, MC_VECTOR) = .TRUE. 
                     ELSEIF (str_comp (cpara (2) , 'lennard', 2, lpara (&
                     2) , 6) ) THEN                                     
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ianz == 3) THEN 
                           werte (4) = 12.0 
                           werte (5) = 6.0 
                        ENDIF 
                        IF(ic > MMC_LENN_CORR .OR.  & ! Allocate Lennard
                           ic > CHEM_MAX_COR  .OR.  &
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
                           IF(is<0 .OR. is>cr_nscat .OR. &
                              js<0 .OR. js>cr_nscat)  THEN
                              ier_num = -97
                              ier_typ = ER_APPL
                              WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
                              WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
                              ier_msg(3) = 'Check parameter list for wrong/missing atoms'
                              RETURN
                           ENDIF
                           mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
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
                        mmc_cor_energy (ic, MC_LENNARD) = .TRUE. 
                        mmc_cor_energy (0, MC_LENNARD) = .TRUE. 
                     ELSEIF (str_comp (cpara (2) , 'repulsive', 2, &
                                       lpara (2) , 9)        ) THEN                                     
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF (ianz >  0) THEN 
                          CALL ber_params (ianz, cpara, lpara, werte,   &
                                           maxw)
                        ENDIF
                        IF (ianz == 0) THEN 
                           werte (1) =  0.0 
                           werte (2) =  1.0 
                           werte (3) =  0.0 
                           werte (4) =  1.0 
                        ELSEIF (ianz == 1) THEN 
                           werte (2) =  1.0 
                           werte (3) =  0.0 
                           werte (4) =  1.0 
                        ELSEIF (ianz == 2) THEN 
                           werte (3) =  0.0 
                           werte (4) =  1.0 
                        ELSEIF (ianz == 3) THEN 
                           werte (4) =  1.0 
                        ENDIF 
                        IF(ic > MMC_REP_CORR .OR.  & ! Allocate Repulsive
                           ic > CHEM_MAX_COR  .OR.  &
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
                           IF(is<0 .OR. is>cr_nscat .OR. &
                              js<0 .OR. js>cr_nscat)  THEN
                              ier_num = -97
                              ier_typ = ER_APPL
                              WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
                              WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
                              ier_msg(3) = 'Check parameter list for wrong/missing atoms'
                              RETURN
                           ENDIF
                           mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
                        CALL mmc_set_disp (ic, MC_REPULSIVE, is, js,    &
                        100.0D0    , ABS(werte (1)) )                          
                        CALL mmc_set_rep  (ic, is, js,    &
                        ABS(werte (1)), werte(2) , werte(3), werte(4) )                          
                        ENDDO 
                        ENDDO 
                        mmc_cor_energy (ic, MC_REPULSIVE) = .TRUE. 
                        mmc_cor_energy (0, MC_REPULSIVE) = .TRUE. 
                     ELSEIF (str_comp (cpara (2) , 'bucking', 2,        &
                                       lpara (2) , 6)           )THEN
                        CALL del_params (2, ianz, cpara, lpara, maxw) 
                        iianz = 1 
                        jjanz = 1 
                        CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        maxw, .FALSE.)                                  
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF(ic > MMC_BUCK_CORR .OR.  & ! Allocate Buckingham
                           ic > CHEM_MAX_COR  .OR.  &
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
                           IF(is<0 .OR. is>cr_nscat .OR. &
                              js<0 .OR. js>cr_nscat)  THEN
                              ier_num = -97
                              ier_typ = ER_APPL
                              WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
                              WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
                              ier_msg(3) = 'Check parameter list for wrong/missing atoms'
                              RETURN
                           ENDIF
                           mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
                           mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
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
                        mmc_cor_energy (ic, MC_BUCKING) = .TRUE. 
                        mmc_cor_energy (0, MC_BUCKING) = .TRUE. 
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
               IF(ier_num==0) mmc_h_number = mmc_h_number + 1
!                                                                       
!------ --- 'set temp' : setting kT for MC simulation                   
!                                                                       
            ELSEIF (cpara (1) (1:2)  == 'TE') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_kt = werte (1) 
!                                                                       
!------ --- 'set valid' : setting number of cycleas after which
!           no valid move is found
!                                                                       
            ELSEIF (cpara (1) (1:2)  == 'VA') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mmc_no_valid = NINT(werte(1) )
!                                                                       
!------ --- 'set vector' : setting of neighbouring vectors              
!                                                                       
            ELSEIF (cpara (1) (1:2)  == 'VE') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL chem_set_vec (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!------ --- 'set fixed' : setting fixed atom coordinate ranges          
!                                                                       
            ELSEIF (cpara (1) (1:3)  == 'FIX') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (str_comp (cpara (1) , 'OFF', 2, lpara (1) , 3) )     &
               THEN                                                     
                  mmc_l_constrains = .FALSE. 
               ELSEIF (str_comp (cpara (1) , 'X', 1, lpara (1) , 1) )   &
               THEN                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_min (1) = werte (1) 
                  mmc_c_max (1) = werte (2) 
                  mmc_constrain_type = MMC_C_XYZ 
               ELSEIF (str_comp (cpara (1) , 'Y', 1, lpara (1) , 1) )   &
               THEN                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_min (2) = werte (1) 
                  mmc_c_max (2) = werte (2) 
                  mmc_constrain_type = MMC_C_XYZ 
               ELSEIF (str_comp (cpara (1) , 'Z', 1, lpara (1) , 1) )   &
               THEN                                                     
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  mmc_c_min (3) = werte (1) 
                  mmc_c_max (3) = werte (2) 
                  mmc_constrain_type = MMC_C_XYZ 
               ELSEIF (str_comp (cpara (1) , 'RAD', 1, lpara (1) , 3) ) &
               THEN                                                     
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
!
!*****7*****************************************************************
!
SUBROUTINE get_finish(ianz, cpara, lpara, werte, MAXW)
!
USE mmc_mod
!
USE errlist_mod
USE precision_mod
USE str_comp_mod
USE take_param_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(INOUT) :: ianz
INTEGER, INTENT(IN)    :: MAXW
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(INOUT) :: werte
!
INTEGER, PARAMETER :: NOPTIONAL = 5
INTEGER, PARAMETER :: O_FEED    = 1                  ! Number of feedbacks to average
INTEGER, PARAMETER :: O_DIFF    = 2                  ! Maximum difference
INTEGER, PARAMETER :: O_CHANGE  = 3                  ! Maximum change in differences
INTEGER, PARAMETER :: O_AVER    = 4                  ! Maximum change in differences
INTEGER, PARAMETER :: O_STOP    = 5                  ! How to stop cycles/convergence
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte    ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 4 ! Number of values to calculate
!
DATA oname  / 'feed  ', 'diff  ', 'change', 'aver  ', 'stop  ' /   ! mmc_set capitalizes only the first parameter
DATA loname /  4      ,  4      ,  6      ,  4      ,  4       /
!
opara  = (/'3     ', '0.0900', '0.0500', '0.0010', 'cycles' /)
lopara = (/ 1      ,  6      ,  6      ,  6      , 6        /)
owerte = (/ 3.     ,  0.0900 ,  0.0500 ,  0.001  , 0.0      /)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
!write(*,*) ' LPRESENT ', lpresent, owerte
IF(ier_num==0) THEN
   IF(opara(O_STOP) == 'converge') THEN
      mmc_h_stop = .TRUE.
   ELSEIF(opara(O_STOP) == 'cycles') THEN
      mmc_h_stop = .FALSE.
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'Optional stop must be stop:converge or '
      ier_msg(2) = '                      stop:cycles '
      RETURN
   ENDIF
   MMC_H_NNNN   = NINT(owerte(O_FEED))    ! Number of feedbacks to average
   mmc_h_conv_m = owerte(O_DIFF)          ! Largest difference (target-achieved)/target
   mmc_h_conv_c = owerte(O_CHANGE)        ! Largest change in difference (target-achieved)/target
   mmc_h_conv_a = owerte(O_AVER  )        ! Average change in difference (target-achieved)/target
ENDIF
!
END SUBROUTINE get_finish
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_disp (ic, ie, is, js, dist, depth) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
USE crystal_mod 
USE mmc_mod 
USE errlist_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!
INTEGER, INTENT(IN) :: ic         ! Correlation number
INTEGER, INTENT(IN) :: ie         ! Energy nmber
INTEGER, INTENT(IN) :: is         ! scattering type atom i
INTEGER, INTENT(IN) :: js         ! Scattering type atom j
REAL(KIND=PREC_DP), INTENT(IN) :: dist 
REAL(KIND=PREC_DP), INTENT(IN) :: depth 
!                                                                       
INTEGER :: ii, jj 
!                                                                       
                                                                        
IF (is /=  - 1.AND.js /=  - 1) THEN 
   mmc_target_corr (ic, ie, is, js) = dist 
   mmc_target_corr (ic, ie, js, is) = dist 
   mmc_depth (ic, ie, is, js) = depth 
   mmc_depth (ic, ie, js, is) = depth 
   mmc_pair (ic, ie, is, js) = -1     
   mmc_pair (ic, ie, js, is) = -1     
ELSEIF (is ==  - 1.AND.js /=  - 1) THEN 
   DO ii = 0, cr_nscat 
      mmc_target_corr (ic, ie, ii, js) = dist 
      mmc_target_corr (ic, ie, js, ii) = dist 
      mmc_depth (ic, ie, ii, js) = depth 
      mmc_depth (ic, ie, js, ii) = depth 
      mmc_pair (ic, ie, ii, js) = -1     
      mmc_pair (ic, ie, js, ii) = -1     
   ENDDO 
ELSEIF (is /=  - 1.AND.js ==  - 1) THEN 
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
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_cn(ic, ie, ianz1, ianz2, MAXW, werte1, werte2, coord, depth) 
!+                                                                      
!     Set desired coordination numbers                                  
!-                                                                      
!                                                                       
USE crystal_mod 
USE mmc_mod 
USE errlist_mod 
USE precision_mod
!
IMPLICIT none 
!

INTEGER           , INTENT(IN) :: ic     ! Correlation number
INTEGER           , INTENT(IN) :: ie     ! Energy number == MC_COORDNUM
INTEGER           , INTENT(IN) :: ianz1  ! No of atom types in first group
INTEGER           , INTENT(IN) :: ianz2  ! No of atom types in second group
INTEGER           , INTENT(IN) :: MAXW   ! Array Dimension 
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN) :: werte1 ! Actual atom types group1
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN) :: werte2 ! Actual atom types group1
REAL(KIND=PREC_DP),                  INTENT(IN) :: coord   ! Desired correlation
REAL(KIND=PREC_DP),                  INTENT(IN) :: depth  ! Energy Depth
!
INTEGER :: is, js,      ii, jj 
!                                                                       
DO ii = 1, ianz1
   is = NINT(werte1(ii))
   mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
   DO jj = 1, ianz2
      js = NINT(werte2(jj))
      mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
      mmc_target_corr (ic, ie, is, js) = coord 
      mmc_depth (ic, ie, is, js) = depth 
      mmc_pair (ic, ie, is, js) = -1     
   ENDDO
ENDDO
!IF (is /=  - 1.AND.js /=  - 1) THEN 
!   mmc_target_corr (ic, ie, is, js) = dist 
!   mmc_depth (ic, ie, is, js) = depth 
!   mmc_pair (ic, ie, is, js) = -1     
!ELSEIF (is ==  - 1.AND.js /=  - 1) THEN 
!   DO ii = 0, cr_nscat 
!      mmc_target_corr (ic, ie, ii, js) = dist 
!!      mmc_depth (ic, ie, ii, js) = depth 
!      mmc_pair (ic, ie, ii, js) = -1     
!      ENDDO 
!ELSEIF (is /=  - 1.AND.js ==  - 1) THEN 
!   DO ii = 0, cr_nscat 
!      mmc_target_corr (ic, ie, is, ii) = dist 
!      mmc_depth (ic, ie, is, ii) = depth 
!      mmc_pair (ic, ie, is, ii) = -1     
!   ENDDO 
!ELSE 
!   DO ii = 0, cr_nscat 
!      DO jj = 0, cr_nscat 
!         mmc_target_corr (ic, ie, ii, jj) = dist 
!         mmc_depth (ic, ie, ii, jj) = depth 
!         mmc_pair (ic, ie, ii, jj) = -1     
!      ENDDO 
!   ENDDO 
!ENDIF 
!                                                                       
END SUBROUTINE mmc_set_cn                   
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_disp_occ (ic, ie, ianz1, ianz2, &
                             MAXW, werte1, werte2, corr, depth )
!
USE crystal_mod 
USE mmc_mod 
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                 , INTENT(IN) :: ic     ! Correlation number
INTEGER                 , INTENT(IN) :: ie     ! Energy number == MC_OCC
INTEGER                 , INTENT(IN) :: ianz1  ! No of atom types in first group
INTEGER                 , INTENT(IN) :: ianz2  ! No of atom types in second group
INTEGER                 , INTENT(IN) :: MAXW   ! Array Dimension 
REAL(KIND=PREC_DP)   , DIMENSION(MAXW), INTENT(IN) :: werte1 ! Actual atom types group1
REAL(KIND=PREC_DP)   , DIMENSION(MAXW), INTENT(IN) :: werte2 ! Actual atom types group1
REAL(KIND=PREC_DP)   ,                  INTENT(IN) :: corr   ! Desired correlation
REAL(KIND=PREC_DP)   ,                  INTENT(IN) :: depth  ! Energy Depth
! 
INTEGER                              :: is, js ! Dummy atom types
INTEGER                              :: i, j   ! Loop indices
!
      DO i=1,ianz1              ! Set "equal" pairs first group
         is = NINT(werte1(i))
         mmc_allowed(is) = .TRUE.
         DO j=1,ianz1
            js = NINT(werte1(j))
         mmc_target_corr (ic, ie, is, js) = corr 
         mmc_depth       (ic, ie, is, js) = depth 
         mmc_pair        (ic, ie, is, js) = +1     ! These pairs contribute positively to energy
         END DO
      END DO
      DO i=1,ianz2              ! Set "equal" pairs second group
         is = NINT(werte2(i))
         mmc_allowed(is) = .TRUE.
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
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_unid_occ (ic, ie, ianz1, ianz2, &
                             MAXW, werte1, werte2, corr, depth )
!
! Chemical correlation, unidirectional case.
! Atoms in group one are at center, atoms in group two are at neighbor position
!
USE crystal_mod 
USE mmc_mod 
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                 , INTENT(IN) :: ic     ! Correlation number
INTEGER                 , INTENT(IN) :: ie     ! Energy number == MC_OCC
INTEGER                 , INTENT(IN) :: ianz1  ! No of atom types in first group
INTEGER                 , INTENT(IN) :: ianz2  ! No of atom types in second group
INTEGER                 , INTENT(IN) :: MAXW   ! Array Dimension 
REAL(KIND=PREC_DP)   , DIMENSION(MAXW), INTENT(IN) :: werte1 ! Actual atom types group1
REAL(KIND=PREC_DP)   , DIMENSION(MAXW), INTENT(IN) :: werte2 ! Actual atom types group1
REAL(KIND=PREC_DP)   ,                  INTENT(IN) :: corr   ! Desired correlation
REAL(KIND=PREC_DP)   ,                  INTENT(IN) :: depth  ! Energy Depth
! 
INTEGER                              :: is, js ! Dummy atom types
INTEGER                              :: i, j   ! Loop indices
!
DO i=1,ianz1              ! Set "equal" pairs first group
   is = NINT(werte1(i))
   mmc_allowed(is) = .TRUE.
   DO j=1,ianz1
      js = NINT(werte1(j))
   mmc_target_corr (ic, ie, is, js) = corr 
   mmc_depth       (ic, ie, is, js) = depth 
   mmc_pair        (ic, ie, is, js) = +1     ! These pairs contribute positively to energy
   END DO
END DO
DO i=1,ianz2              ! Set "equal" pairs second group
   is = NINT(werte2(i))
   mmc_allowed(is) = .TRUE.
   DO j=1,ianz1
      js = NINT(werte2(j))
      mmc_target_corr (ic, ie, is, js) = corr 
      mmc_depth       (ic, ie, is, js) = depth 
      mmc_pair        (ic, ie, is, js) = +2     ! These pairs contribute positively to energy
   END DO
END DO
!
DO i=1,ianz1              ! Set "opposite" pairs
   is = NINT(werte1(i))
   DO j=1,ianz2
      js = NINT(werte2(j))
      mmc_target_corr (ic, ie, is, js) = corr 
      mmc_target_corr (ic, ie, js, is) = corr 
      mmc_depth       (ic, ie, is, js) = depth 
      mmc_depth       (ic, ie, js, is) = depth 
      mmc_pair        (ic, ie, is, js) = -1     ! These pairs contribute negatively to energy
      mmc_pair        (ic, ie, js, is) = +2     ! Second ==> first group
   END DO
END DO
!
!
END SUBROUTINE mmc_set_unid_occ
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_lenn (ic, is, js, a, b, m, n) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
USE crystal_mod 
USE mmc_mod 
USE errlist_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: is
INTEGER, INTENT(IN) :: js
REAL(KIND=PREC_DP), INTENT(IN) :: a, b, m, n 
!
INTEGER :: ii, jj 
!                                                                       
                                                                        
      IF (is /=  - 1.AND.js /=  - 1) THEN 
         mmc_len_a (ic, is, js) = a 
         mmc_len_b (ic, is, js) = b 
         mmc_len_m (ic, is, js) = m 
         mmc_len_n (ic, is, js) = n 
         mmc_len_a (ic, js, is) = a 
         mmc_len_b (ic, js, is) = b 
         mmc_len_m (ic, js, is) = m 
         mmc_len_n (ic, js, is) = n 
      ELSEIF (is ==  - 1.AND.js /=  - 1) THEN 
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
      ELSEIF (is /=  - 1.AND.js ==  - 1) THEN 
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
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_rep (ic, is, js, a, b, c, m) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
USE crystal_mod 
USE mmc_mod 
USE errlist_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: is
INTEGER, INTENT(IN) :: js
REAL(KIND=PREC_DP), INTENT(IN) :: a, b, c, m 
!
INTEGER :: ii, jj 
!                                                                       
!     Define mimimum energy at infinity distance
!
      mmc_rep_low = MIN(mmc_rep_low, REAL(ABS(a)))
                                                                        
      IF (is /=  - 1.AND.js /=  - 1) THEN 
         mmc_rep_a (ic, is, js) = a 
         mmc_rep_b (ic, is, js) = b 
         mmc_rep_c (ic, is, js) = c 
         mmc_rep_m (ic, is, js) = m 
         mmc_rep_a (ic, js, is) = a 
         mmc_rep_b (ic, js, is) = b 
         mmc_rep_c (ic, js, is) = c 
         mmc_rep_m (ic, js, is) = m 
      ELSEIF (is ==  - 1.AND.js /=  - 1) THEN 
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
      ELSEIF (is /=  - 1.AND.js ==  - 1) THEN 
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
!
!*****7*****************************************************************
!     SUBROUTINE mmc_set_vec (ic, is, js, werte, maxw) 
!+                                                                      
!     Set desired displacements                                         
!-                                                                      
!                                                                       
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
!     IF (is /=  - 1.AND.js /=  - 1) THEN 
!        mmc_nvec (ic, js, is) = mmc_nvec (ic, js, is) + 1 
!        mmc_nvec (ic, is, js) = mmc_nvec (ic, js, is) 
!        nv = mmc_nvec (ic, is, js) 
!        DO i = 1, 4 
!        mmc_vec (i, nv, ic, is, js) = werte (i) 
!        mmc_vec (i, nv, ic, js, is) = werte (i) 
!        ENDDO 
!     ELSEIF (is ==  - 1.AND.js /=  - 1) THEN 
!        DO ii = 0, cr_nscat 
!        mmc_nvec (ic, js, ii) = mmc_nvec (ic, js, ii) + 1 
!        mmc_nvec (ic, ii, js) = mmc_nvec (ic, js, is) 
!        nv = mmc_nvec (ic, ii, js) 
!        DO i = 1, 4 
!        mmc_vec (i, nv, ic, ii, js) = werte (i) 
!        mmc_vec (i, nv, ic, js, ii) = werte (i) 
!        ENDDO 
!        ENDDO 
!     ELSEIF (is /=  - 1.AND.js ==  - 1) THEN 
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
!
SUBROUTINE mmc_set_mode (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     Sets MMC    mode                                                  
!-                                                                      
USE crystal_mod 
USE get_iscat_mod
USE rmc_mod 
USE mmc_mod 
USE modify_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) ::  MAXW 
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: cpara !(MAXW) 
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara !(MAXW) 
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(INOUT) :: werte !(MAXW) 
!                                                                       
INTEGER :: ianz, imode=MC_MOVE_NONE, i 
INTEGER :: is 
REAL :: sump 
!                                                                       
      IF (ianz >= 1) THEN 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3)  == 'SHI') THEN 
            imode = MC_MOVE_DISP 
         ELSEIF (cpara (1) (1:3)  == 'SWD') THEN 
            imode = MC_MOVE_SWDISP 
         ELSEIF (cpara (1) (1:3)  == 'SWC') THEN 
            imode = MC_MOVE_SWCHEM 
         ELSEIF (cpara (1) (1:3)  == 'INV') THEN 
            imode = MC_MOVE_INVDISP 
         ELSE 
            ier_typ = ER_RMC 
            ier_num = - 9 
         ENDIF 
!                                                                       
         IF (ianz == 1) THEN 
            mmc_local (imode) = rmc_local_all 
         ELSEIF (ianz >= 2) THEN 
            CALL do_cap (cpara (2) ) 
            IF (cpara (2) (1:1)  == 'A') THEN 
               mmc_local (imode) = rmc_local_all 
            ELSEIF (cpara (2) (1:1)  == 'L') THEN 
               mmc_local (imode) = rmc_local_loc 
            ELSEIF (cpara (2) (1:2)  == 'SL') THEN 
               mmc_local (imode) = rmc_local_locsite 
            ELSEIF (cpara (2) (1:2)  == 'SI') THEN 
               mmc_local (imode) = rmc_local_site 
            ELSEIF (cpara (2) (1:2)  == 'CO') THEN 
               mmc_local (imode) = rmc_local_conn 
            ELSE 
               ier_typ = ER_RMC 
               ier_num = - 9 
            ENDIF 
            IF (ianz > 2) THEN 
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL get_iscat (ianz, cpara, lpara, werte, maxw, .FALSE.) 
               IF (ier_num /= 0) return 
!                                                                       
               IF (ier_num == 0) THEN 
                  IF(NINT(werte(1)) == -1) THEN 
                     DO i = 0, cr_nscat 
                     mmc_allowed (i) = .TRUE. 
                     ENDDO 
                  ELSE 
                     DO i = 1, ianz 
                     is = nint (werte (i) ) 
                     IF (is >= 0.AND.is <= cr_nscat) THEN 
                        mmc_allowed (is) = .TRUE. 
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
         IF (sump > 0) THEN 
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
!
!*****7*****************************************************************
!
SUBROUTINE mmc_run_multi (lout_feed)
!+                                                                      
!     This is the MC routine for multiple energy calculations           
!-                                                                      
USE discus_allocate_appl_mod
USE crystal_mod 
USE chem_mod 
USE chem_neig_multi_mod
USE chem_menu
!
USE celltoindex_mod
USE atom_env_mod
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
!$ USE omp_lib
USE parallel_mod
!
USE debug_mod 
USE errlist_mod 
USE lib_random_func
USE param_mod 
USE precision_mod
USE prompt_mod 
USE support_mod
!                                                                       
IMPLICIT none 
!
LOGICAL, INTENT(IN) :: lout_feed     ! Write output upon feedback?
!                                                                       
CHARACTER(LEN=24) :: c_energy (0:MC_N_ENERGY) 
!
REAL :: start, zeit
REAL(KIND=PREC_SP), DIMENSION(:,:,:), ALLOCATABLE :: disp !(3, 0:MAX_ATOM_ENV, 2) 
REAL, DIMENSION(3) :: idir, jdir
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: rdi ! (CHEM_MAX_COR) 
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: rdj ! (CHEM_MAX_COR) 
REAL :: rel_cycl    ! how far are we in the desired number of cycles
REAL(KIND=PREC_SP), DIMENSION(:,:,:), ALLOCATABLE :: patom ! (3, 0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER(KIND=PREC_INT_LARGE) :: itry, igen, imodulus
INTEGER           , DIMENSION(  :,:), ALLOCATABLE :: iatom ! (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: natom ! (MMC_MAX_CENT) 
INTEGER :: iacc_good, iacc_neut, iacc_bad 
INTEGER :: isel(MMC_MAX_ATOM) 
INTEGER :: lbeg (3) 
INTEGER :: ic
INTEGER :: nocc
INTEGER :: i, natoms
INTEGER :: ncent 
INTEGER :: NALLOWED   ! Current size mmc_allowed
INTEGER :: zh, zm, zs 
INTEGER :: n_angles   ! Dummy for angle allocation
LOGICAL, DIMENSION(3) :: old_chem_period
LOGICAL :: lserial    ! serial calculation if TRUE
LOGICAL :: loop, laccept, done 
LOGICAL :: lout, lfinished
!                                                                       
REAL, DIMENSION(0:MC_N_ENERGY) :: e_old
REAL, DIMENSION(0:MC_N_ENERGY) :: e_new
!
INTEGER :: tid
INTEGER(KIND=PREC_INT_LARGE) :: nthreads
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
     'Repulsive     Potential ',   &
     'Coordination Number     ',   &
     'Unidirectional Corr     ' /
!
old_chem_period = chem_period
DO i=1,3
  IF(cr_icc(i)==1) chem_period(i) = .FALSE.
ENDDO
!CALL alloc_mmc      (CHEM_MAX_COR, MC_N_ENERGY, MAXSCAT, 
IF (mmc_n_angles == 0          .OR.     &
    mmc_n_angles >= MMC_MAX_ANGLES) THEN 
   n_angles = max(mmc_n_angles+20, int(MMC_MAX_ANGLES*1.025))
   CALL alloc_mmc_angle (CHEM_MAX_COR,n_angles)
ENDIF
IF(CHEM_MAX_COR>UBOUND(mmc_buck_a,1) .OR. MAXSCAT>UBOUND(mmc_buck_a,2)) THEN
   CALL alloc_mmc_buck( CHEM_MAX_COR, MAXSCAT )
ENDIF
IF(CHEM_MAX_COR>UBOUND(mmc_len_a ,1) .OR. MAXSCAT>UBOUND(mmc_len_a ,2)) THEN
   CALL alloc_mmc_lenn( CHEM_MAX_COR, MAXSCAT )
ENDIF
IF(CHEM_MAX_COR>UBOUND(mmc_rep_a ,1) .OR. MAXSCAT>UBOUND(mmc_rep_a ,2)) THEN
   CALL alloc_mmc_rep ( CHEM_MAX_COR, MAXSCAT )
ENDIF
!
ALLOCATE(disp (3, 0:MAX_ATOM_ENV, 2))
ALLOCATE(patom(3, 0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(rdi(CHEM_MAX_COR))
ALLOCATE(rdj(CHEM_MAX_COR))
ALLOCATE(iatom(   0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(natom(                    MMC_MAX_CENT))
IF(ALLOCATED(mmc_h_diff) ) DEALLOCATE(mmc_h_diff)
IF(ALLOCATED(mmc_h_targ)) DEALLOCATE(mmc_h_targ)
IF(ALLOCATED(mmc_h_aver)) DEALLOCATE(mmc_h_aver)
IF(ALLOCATED(mmc_h_maxd)) DEALLOCATE(mmc_h_maxd)
ALLOCATE(mmc_h_diff(mmc_h_number+10,0:MMC_H_NNNN-1))
ALLOCATE(mmc_h_targ(mmc_h_number+10))
ALLOCATE(mmc_h_aver(mmc_h_number+10))
ALLOCATE(mmc_h_maxd(mmc_h_number+10, 0:MMC_H_NNNN-1))
mmc_h_diff= 0.0
mmc_h_targ= 0.0
mmc_h_aver= 1.0E20
mmc_h_maxd= 1.0E20
!ALLOCATE(disp(3,0:MMC_MAX_ATOM_L,2))
!
!DO i = 0, MC_N_ENERGY 
n_e_av_p = 0    !(i) = 0 
n_e_av_m = 0    !(i) = 0 
n_e_av_z = 0    !(i) = 0 
e_aver_p = 0.0  !(i) = 0.0 
e_aver_m = 0.0  !(i) = 0.0 
!ENDDO 
!                                                                       
!------ reset some counters                                             
!                                                                       
igen = 0 
itry = 0 
iacc_good = 0 
iacc_neut = 0 
iacc_bad = 0 
loop = .TRUE. 
done = .TRUE. 
lbeg (1) = - 1 
lbeg (2) = 0 
lbeg (3) = 0 
!                                                                       
!     Normalize the correlation directions                              
!                                                                       
nocc = 0
DO ic = 1, chem_ncor 
   IF (mmc_cor_energy (0, MC_DISP) ) THEN 
      DO i = 1, 3 
         idir (i) = chem_dir (i, 1, ic) 
         jdir (i) = chem_dir (i, 2, ic) 
      ENDDO 
      rdi (ic) = skalpro (idir, idir, cr_gten) 
      rdj (ic) = skalpro (jdir, jdir, cr_gten) 
      rdi (ic) = sqrt (rdi (ic) ) 
      rdj (ic) = sqrt (rdj (ic) ) 
   ELSEIF(mmc_cor_energy(ic, MC_OCC)) THEN
      nocc = nocc + 1
   ENDIF 
ENDDO 
IF(nocc>1) THEN
   DO ic = 1, chem_ncor
      IF(mmc_cor_energy(ic, MC_OCC)) THEN
         mmc_cfac(ic,MC_OCC) = 2.5/REAL(nocc)
      ENDIF
   ENDDO
ENDIF
!read(*,*) ic
!                                                                       
!     Initialize the different energies                                 
!                                                                       
lout = .FALSE. 
lfinished = .FALSE.
CALL mmc_correlations (lout, 0.0, done, lfinished) 
IF(ier_num /= 0) RETURN 
!
!     current size of mmc_allowed
!
NALLOWED = UBOUND(mmc_allowed,1)
!                                                                       
loop = .TRUE.
IF(mo_cyc == 0) loop = .FALSE. 
!
lserial = .FALSE.
check_conn: DO IC = 1, chem_ncor         ! Check is we have a connectivity
   IF(chem_ctyp(ic) == CHEM_CON) THEN    ! if so, we need serial computation
      lserial = .TRUE.                   ! Needs serious debugging...
      EXIT check_conn
   ENDIF
ENDDO check_conn
tid = 0
nthreads = 1
IF(.NOT.lserial .AND. par_omp_use) THEN
!$OMP PARALLEL PRIVATE(tid)
!$   tid = OMP_GET_THREAD_NUM()
!$   IF (tid == 0) THEN
!$      IF(par_omp_maxthreads == -1) THEN
!$         nthreads = OMP_GET_NUM_THREADS()
!$      ELSE
!$         nthreads = MAX(1,MIN(par_omp_maxthreads, OMP_GET_NUM_THREADS()))
!$      ENDIF
!$   END IF
!$OMP END PARALLEL
ENDIF
imodulus=MAX(1_PREC_INT_LARGE, mo_feed/nthreads)
!                                                                       
!------ Main MC loop here                                               
!                                                                       
start = seknds(0.0) 
!
IF(ier_ctrlc) THEN
      DEALLOCATE(patom)
      DEALLOCATE(disp)
      DEALLOCATE(rdi)
      DEALLOCATE(rdj)
      DEALLOCATE(iatom)
      DEALLOCATE(natom)
   ier_num = -14
   ier_typ = ER_COMM
   chem_period = old_chem_period
   RETURN
ENDIF
IF(ier_num/=0) THEN
   IF(ALLOCATED(patom)) DEALLOCATE(patom)
   IF(ALLOCATED(disp )) DEALLOCATE(disp )
   IF(ALLOCATED(rdi  )) DEALLOCATE(rdi  )
   IF(ALLOCATED(rdj  )) DEALLOCATE(rdj  )
   IF(ALLOCATED(iatom)) DEALLOCATE(iatom)
   IF(ALLOCATED(natom)) DEALLOCATE(natom)
   chem_period = old_chem_period
   RETURN      ! An error occured or CTRL-C
ENDIF
done = .FALSE. 
!
!!! !$OMP PARALLEL PRIVATE(tid, isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
!!!$OMP                  disp, kgen, ktry, e_old, e_new)                     &
IF(nthreads > 1) THEN
   !$OMP PARALLEL PRIVATE(tid, isel, natoms, &
   !$OMP                  iatom, patom, natom,ncent, laccept,                 &
   !$OMP                  disp,             e_old, e_new)                     &
   !$OMP          SHARED(done)
   !$   tid = OMP_GET_THREAD_NUM()
   !$OMP DO SCHEDULE(DYNAMIC, mo_cyc/nthreads/32)
   parallel_loop: DO itry=1, mo_cyc     ! Do mmc in parallel
      IF(done) CYCLE parallel_loop    ! Quickly cycle to end if an error occuredd
      CALL    mmc_run_loop(tid, nthreads, igen, itry, natoms, &
                           iatom, patom, natom,ncent, laccept,                       &
                           rdi, rdj,         e_old, e_new, done, loop,               &
                           iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                           lout_feed, imodulus,                                      &
                           NALLOWED, MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
!
   ENDDO parallel_loop
   !$OMP END DO NOWAIT
   !$ IF(tid==0) igen = igen*nthreads
   !$OMP END PARALLEL
   itry = igen
ELSE     ! Use nonparallel code
   serial_loop: DO itry=1, mo_cyc     ! Do mmc serially
      CALL    mmc_run_loop(tid, nthreads, igen, itry, &
                           natoms, &
                           iatom, patom, natom,ncent, laccept,                       &
                           rdi, rdj,         e_old, e_new, done, loop,               &
                           iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                           lout_feed, imodulus,                                      &
                           NALLOWED, MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
      IF(ier_num/=0 .OR. done) EXIT serial_loop
   ENDDO serial_loop
ENDIF
!
IF(ier_num/=0) THEN
   IF(ALLOCATED(patom)) DEALLOCATE(patom)
   IF(ALLOCATED(disp )) DEALLOCATE(disp )
   IF(ALLOCATED(rdi  )) DEALLOCATE(rdi  )
   IF(ALLOCATED(rdj  )) DEALLOCATE(rdj  )
   IF(ALLOCATED(iatom)) DEALLOCATE(iatom)
   IF(ALLOCATED(natom)) DEALLOCATE(natom)
   chem_period = old_chem_period
   RETURN      ! An error occured or CTRL-C
ENDIF
!                                                                       
!------ Loop finished                                                   
!                                                                       
IF(lout_feed) THEN
   WRITE (output_io, 3000) 
   WRITE (output_io, 2000) igen, itry, iacc_good, iacc_neut, iacc_bad 
ENDIF
!     lout = .TRUE. 
lfinished = .TRUE.
CALL mmc_correlations (lout_feed, rel_cycl, done, lfinished) 
!                                                                       
!     Give average energy changes for the different energy terms        
!                                                                       
IF(lout_feed) WRITE ( output_io, 5000) 
DO i = 1, MC_N_ENERGY 
   IF(n_e_av_p(i) > 0) THEN 
      e_aver_p(i) = e_aver_p(i) / REAL(n_e_av_p(i) ) 
   ENDIF 
   IF(n_e_av_m(i) > 0) THEN 
      e_aver_m(i) = e_aver_m(i) / REAL(n_e_av_m(i) ) 
   ENDIF 
   IF(lout_feed) WRITE(output_io, 5010) c_energy(i), n_e_av_m(i), e_aver_m(i),        &
                                        n_e_av_z(i), n_e_av_p(i), e_aver_p(i)
   n_e_av_p (i) = 0 
   n_e_av_m (i) = 0 
   n_e_av_z (i) = 0 
   e_aver_p (i) = 0.0 
   e_aver_m (i) = 0.0 
ENDDO 
!                                                                       
!------ Write timing results                                            
!                                                                       
IF(lout_feed) THEN
   zeit = seknds (start) 
   zh = int (zeit / 3600.) 
   zm = int ( (zeit - zh * 3600.) / 60.) 
   zs = int (zeit - zh * 3600 - zm * 60.) 
   WRITE (output_io, 4000) zh, zm, zs, zeit / mo_cyc 
ENDIF
!
DEALLOCATE(patom)
DEALLOCATE(disp )
DEALLOCATE(rdi)
DEALLOCATE(rdj)
DEALLOCATE(iatom)
DEALLOCATE(natom)
!
chem_period = old_chem_period
!                                                                       
 2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')                                 
 3000 FORMAT (/,' --- Final multiple energy configuration ---') 
 4000 FORMAT (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/     &
     &          ' Time/cycle   : ',F9.3,' sec',/)                       
 5000 FORMAT (/,' --- Average energy changes ---',//,                   &
     &        ' Energy type',14x,'number',17x,'number',2x,'number',/,   &
     &        29x,'E < 0',5x,' <E>',11x,'E = 0',5x,'E > 0',4x,' <E>')   
 5010 FORMAT(a24, 2x,i8,f15.6,2x,i8,2x,i8,f15.6) 
!                                                                       
END SUBROUTINE mmc_run_multi                  
!
!*******************************************************************************
!
SUBROUTINE mmc_run_loop(tid, nthreads, igen, itry, &
                        natoms, &
                        iatom, patom, natom,ncent, laccept,                       &
                        rdi, rdj,         e_old, e_new, done, loop,               &
                        iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                        lout_feed, imodulus,                                      &
                        NALLOWED, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!
USE crystal_mod
USE chem_mod
USE chem_menu
USE mc_mod
USE mmc_mod
!
USE errlist_mod
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                           , INTENT(IN) :: NALLOWED
INTEGER                           , INTENT(IN) :: MAX_ATOM_ENV_L 
INTEGER                           , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                           , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER                           , INTENT(IN)  :: tid
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)  :: nthreads
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(INOUT)  :: igen
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)     :: itry
INTEGER                           , INTENT(OUT) :: natoms
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
INTEGER                                                 , INTENT(INOUT) :: ncent
LOGICAL                           , INTENT(OUT) :: laccept
REAL(KIND=PREC_SP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdi ! (CHEM_MAX_COR) 
REAL(KIND=PREC_SP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdj ! (CHEM_MAX_COR) 
REAL   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_old
REAL   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_new
LOGICAL                           , INTENT(INOUT) :: done
LOGICAL                           , INTENT(OUT) :: loop
INTEGER                           , INTENT(INOUT)  :: iacc_good
INTEGER                           , INTENT(INOUT)  :: iacc_neut
INTEGER                           , INTENT(INOUT)  :: iacc_bad
REAL                              , INTENT(INOUT)  :: rel_cycl
LOGICAL                           , INTENT(IN ) :: lout_feed
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)  :: imodulus
!
INTEGER, DIMENSION(MMC_MAX_ATOM_L)              :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)                           :: is
INTEGER, DIMENSION(2, 3)                        :: iz
INTEGER, DIMENSION(3)                           :: iz1
INTEGER, DIMENSION(3)                           :: iz2
INTEGER                                         :: iselz
INTEGER                                         :: iselz2
LOGICAL :: valid_all
REAL   , DIMENSION(3, 0:nthreads-1)                     :: posz !(3) = 0.0
REAL   , DIMENSION(3, 0:nthreads-1)                     :: posz2 !(3) = 0.0
REAL(KIND=PREC_SP), DIMENSION(3,0:MAX_ATOM_ENV_l,2) :: disp
   IF(tid==0) igen = igen + 1
!
   IF(done) RETURN                 ! Quickly cycle to end if an error occuredd
!
!  -- Choose move and atoms
!
   CALL mmc_select_atoms(isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                         laccept, loop, NALLOWED, MMC_MAX_ATOM_L)
   IF(ier_num/=0) THEN             ! Error, cycle to end of loop
      done = .TRUE.
      RETURN
   ENDIF
!
!-- Try move                                                      
!                                                                       
!write(*,*) ' ISEL ', isel(1), i, ' : ', isel(2), j , string(1:4)
!write(*,*) ' iscat', cr_iscat(isel(1)), cr_iscat(i), ' : ', cr_iscat(isel(2)), cr_iscat(j)
!write(*,*) ' ISEL ', isel(1), ' : ', isel(2)
!write(*,*) ' iscat', cr_iscat(isel(1)),cr_iscat(isel(2))
!read(*,*) i
!                                                                       
!--Calculate old energy                                          
!                                                                       
   e_old      = 0.0    ! e_old(:)
!
   valid_all = .FALSE.
   CALL mmc_energies(isel, is, iz, natoms, iatom, patom, natom, ncent, &
                     rdi, rdj, valid_all, e_old, CHEM_MAX_COR,         &
                     MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   IF(ier_num/=0) THEN             ! Error, cycle to end of loop
      done = .TRUE.
      RETURN
   ENDIF
   IF(valid_all) THEN 
!                                                                       
!-----Perform the modifications of the atoms                        
!                                                                       
      CALL mmc_modify(isel, posz(:,tid), posz2(:,tid), disp, MMC_MAX_ATOM_L)
      IF(ier_num/=0) THEN             ! Error, cycle to end of loop
         done = .TRUE.
         RETURN
      ENDIF
!                                                                       
!-----Calculate new energy                                          
!                                                                       
      e_new = 0.0   !      e_new (ie) = 0.0 
!                                                                       
!-----Set the assumption of at least one propper energy to FALSE  
!     Will be set to TRUE if at least one energy calculation is fine
!                                                                       
      valid_all = .FALSE. 
      CALL mmc_energies(isel, is, iz, natoms, iatom, patom, natom, ncent, &
                        rdi, rdj, valid_all, e_new, CHEM_MAX_COR,         &
                        MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
      IF(ier_num/=0) THEN             ! Error, cycle to end of loop
         done = .TRUE.
         RETURN
      ENDIF
!                                                                       
!     ------The comparison of the energies is done only if a proper     
!           new energy was found. Otherwise the move is automatically   
!           rejected. This might happen, if an atom has moved outside   
!           the allowed sphere of influence.                            
!                                                                       
      IF (valid_all) THEN 
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
         CALL mmc_test_multi (iacc_good, iacc_neut, iacc_bad,  &
              e_new, e_old,  laccept)                                                 
      ELSE 
         laccept = .FALSE. 
      ENDIF 
      IF(ier_num/=0) THEN             ! Error, cycle to end of loop
         done = .TRUE.
         RETURN
      ENDIF
!                                                                       
!     ----The move was not accepted, move atoms back to old places      
!                                                                       
      IF (.NOT.laccept) THEN 
         CALL mmc_unmodify(isel, posz(:,tid), posz2(:,tid), MMC_MAX_ATOM_L)
      ELSE     ! Move accepted check periodic bounday conditions
         IF(.NOT.chem_quick .AND.                                 &
            (chem_period(1).OR.chem_period(2).OR.chem_period(3))) THEN
!                 normal and periodic mode
            CALL chem_apply_period(iselz, .TRUE.)
         ENDIF 
      ENDIF 
      IF(ier_num/=0) THEN             ! Error, cycle to end of loop
         done = .TRUE.
         RETURN
      ENDIF
   ENDIF 
!                                                                       
!     --End of modification of atoms if a proper "old" energy was found 
!                                                                       
   IF(tid==0) THEN
      IF(MOD(igen, imodulus)==0) THEN ! .AND. .NOT.done .AND. loop) THEN
!         done = .TRUE. 
         IF(lout_feed) WRITE (output_io, 2000) igen*nthreads, itry*nthreads, &
             iacc_good, iacc_neut, iacc_bad 
!                                                                       
!     ----New mmc_correlations for all energies                         
!                                                                       
         rel_cycl = REAL(igen)/REAL(mo_cyc)*REAL(NTHREADS)
         CALL mmc_correlations (lout_feed, rel_cycl, done, .FALSE.)
      ENDIF 
      IF(igen> mo_cyc/nthreads) THEN
         done = .TRUE.
      ENDIF
   ENDIF 
!
 2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')                                 
!
END SUBROUTINE mmc_run_loop
!
!*****7*****************************************************************
!
SUBROUTINE mmc_select_atoms(isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                            laccept, loop, NALLOWED, MMC_MAX_ATOM_L)
!                                                                       
!     --- Choose a move at random                                       
!                                                                       
USE crystal_mod
USE chem_mod
USE celltoindex_mod
USE metric_mod
USE mc_mod
USE mmc_mod
USE modify_func_mod
USE rmc_sup_mod
!
USE errlist_mod
USE lib_random_func
!
!SAVE
!
IMPLICIT NONE
!
INTEGER                           , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER                           , INTENT(IN)  :: NALLOWED   ! Current size mmc_allowed
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(OUT) :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)             , INTENT(OUT) :: is
INTEGER                           , INTENT(OUT) :: iselz
INTEGER                           , INTENT(OUT) :: iselz2
INTEGER                           , INTENT(OUT) :: natoms
INTEGER, DIMENSION(2, 3)          , INTENT(OUT) :: iz
INTEGER, DIMENSION(3)             , INTENT(OUT) :: iz1
INTEGER, DIMENSION(3)             , INTENT(OUT) :: iz2
LOGICAL                           , INTENT(OUT) :: laccept
LOGICAL                           , INTENT(OUT) :: loop
!
INTEGER  :: i, j
REAL :: z
REAL, DIMENSION(3) ::  v, u
REAL :: r1
!
laccept = .TRUE. 
j = 0
main: DO
   j = j + 1
   IF(j>10000000) THEN   ! Failure did not find a valid pair
      isel = 0
            ier_num = -22
            ier_typ = ER_RMC
            ier_msg(1) = 'RMC did not find a valid pair'
            ier_msg(2) = 'Check composition and properties'
      RETURN
   ENDIF
!                                                                       
!                                                                       
!------ - selecting sites                                               
!                                                                       
   natoms = 1 
   CALL RANDOM_NUMBER(z)
   i = 1 
   DO WHILE (z >= mmc_move_cprob (i) .AND. i <  MC_N_MOVE) 
      i = i + 1 
   ENDDO 
      mmc_move = i 
      mo_local = mmc_local (i) 
!                                                                       
      isel(:) = 0
      IF (mmc_move == MC_MOVE_DISP.OR.mmc_move == MC_MOVE_INVDISP) THEN 
         IF (mmc_l_limited) THEN 
            CALL mmc_limit_selection (isel, natoms, MMC_MAX_ATOM) 
         ELSE 
            CALL RANDOM_NUMBER(r1)
            isel (1) = INT(r1          * cr_natoms) + 1 
!           isel (1) = INT(ran1 (idum) * cr_natoms) + 1 
         ENDIF 
         isel(1) = MAX(1,MIN(isel(1),cr_natoms))
         iselz = isel (1) 
!         IF (isel (1)  > cr_natoms.OR.isel (1)  < 1) CYCLE inner
            laccept = mmc_allowed(cr_iscat(isel(1) ) ) .AND. &
                      check_select_status(isel(1), .TRUE., cr_prop (isel (1) ),  cr_sel_prop)                  
         IF (cr_ncatoms > 0) THEN 
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
      ELSEIF (mmc_move == MC_MOVE_SWDISP) THEN 
         natoms = 2 
         CALL rmc_select (mo_local, isel, iz1, iz2, is (1), is (2) , &
                          NALLOWED, MMC_MAX_ATOM, mmc_allowed)                                                       
         IF(isel(2)==0) THEN
!           DEALLOCATE(patom)
!           DEALLOCATE(disp )
!           DEALLOCATE(rdi )
!           DEALLOCATE(rdj )
!           DEALLOCATE(iatom)
!           DEALLOCATE(natom)
            ier_num = -22
            ier_typ = ER_RMC
            ier_msg(1) = 'RMC did not find a valid pair'
            ier_msg(2) = 'Check composition and properties'
            RETURN
         ENDIF
         iselz = isel (1) 
         iselz2 = isel (2) 
         DO i = 1, 3 
         iz (1, i) = iz1 (i) 
         iz (2, i) = iz2 (i) 
         ENDDO 
         DO i = 1, 3 
            v(i) = cr_pos(i, isel( 1)) - chem_ave_pos(i, is( 1) ) &
                   -REAL(iz( 1, i) - 1) - cr_dim0(i, 1)
            u(i) = cr_pos(i, isel( 2)) - chem_ave_pos(i, is( 2) ) &
                   -REAL(iz( 2, i) - 1) - cr_dim0(i, 1)
         ENDDO 
         laccept = mmc_allowed(cr_iscat(isel(1) ) ) .AND.&
                   mmc_allowed(cr_iscat(isel(2) ) ) .AND.&
                   check_select_status(isel(1),.TRUE., cr_prop(isel(1) ), cr_sel_prop) .AND.&
                   check_select_status(isel(2),.TRUE., cr_prop(isel(2) ), cr_sel_prop) .AND.&
                   skalpro(u, v, cr_gten) < 0.0
      ELSEIF (mmc_move == MC_MOVE_SWCHEM) THEN 
         natoms = 2 
         CALL rmc_select (mo_local, isel, iz1, iz2, is (1), is (2) , &
                          NALLOWED, MMC_MAX_ATOM, mmc_allowed)                                                       
         IF(isel(2)==0) THEN
!           DEALLOCATE(patom)
!           DEALLOCATE(disp )
!           DEALLOCATE(rdi )
!           DEALLOCATE(rdj )
!           DEALLOCATE(iatom)
!           DEALLOCATE(natom)
            ier_num = -22
            ier_typ = ER_RMC
            ier_msg(1) = 'RMC did not find a valid pair'
            ier_msg(2) = 'Check composition and properties'
            RETURN
         ENDIF
         iselz = isel (1) 
         iselz2 = isel (2) 
         DO i = 1, 3 
         iz (1, i) = iz1 (i) 
         iz (2, i) = iz2 (i) 
         ENDDO 
         laccept = cr_iscat (isel (1) )  /= cr_iscat (isel (2) ) .AND.                  &
                 ( mmc_allowed (cr_iscat (isel (1) ) ) .AND.                            &
                   mmc_allowed (cr_iscat (isel (2) ) )      )    .AND.                  &
                   check_select_status (isel(1), .TRUE., cr_prop (isel (1) ), cr_sel_prop) .AND. &
                   check_select_status (isel(2), .TRUE., cr_prop (isel (2) ), cr_sel_prop)                              
      ENDIF 
!                                                                       
!-----      ----Check whether geometrical constrains apply              
!                                                                       
      IF (laccept) THEN 
         CALL check_geometry_multi (isel, natoms, laccept, MMC_MAX_ATOM) 
      ENDIF 
      IF(laccept) EXIT main
ENDDO main
!
END SUBROUTINE mmc_select_atoms
!
!*****7*****************************************************************
!
SUBROUTINE mmc_energies(isel, is, iz, natoms, iatom, patom, natom, ncent, &
                        rdi, rdj, valid_all, e_cur, CHEM_MAX_COR_L,       &
                        MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!-
!
!+
!
USE crystal_mod
USE chem_mod
USE chem_neig_multi_mod
USE mmc_mod
USE  mc_mod
USE metric_mod
!
use errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                                                 , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                                                 , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                                                 , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER                                                 , INTENT(IN) :: CHEM_MAX_COR_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L)                      , INTENT(IN) :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)                                   , INTENT(IN) :: is
INTEGER                                                 , INTENT(IN) :: natoms
INTEGER, DIMENSION(2, 3)                                , INTENT(IN) :: iz
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
INTEGER                                                 , INTENT(INOUT) :: ncent
REAL(KIND=PREC_SP), DIMENSION(CHEM_MAX_COR_L)           , INTENT(INOUT) :: rdi
REAL(KIND=PREC_SP), DIMENSION(CHEM_MAX_COR_L)           , INTENT(INOUT) :: rdj
LOGICAL                                                 , INTENT(OUT) :: valid_all
REAL   , DIMENSION(0:MC_N_ENERGY)                         , INTENT(INOUT) :: e_cur
!
INTEGER :: i      ! Dummy loop variable
INTEGER :: ia     ! Dummy loop variable for modified atoms
INTEGER :: ic     ! Dummy loop variable for correlations
INTEGER :: icent  ! Dummy loop variable for centers
REAL    :: delta
REAL, DIMENSION(3) :: v
REAL, DIMENSION(3) :: idir, jdir
LOGICAL            :: valid_e
!
!                                                                       
!     ----Loop over all modified atoms                                  
!                                                                       
loop_natoms:DO ia = 1, natoms 
!                                                                       
!     ------Set the assumption of at least one proper energy to FALSE  
!           Will be set to TRUE if at least one energy calculation      
!               is fine                                                 
!                                                                       
   valid_all = .FALSE. 
!                                                                       
!     ------Loop over all defined neighbour interactions                
!                                                                       
!write(*,*) ' OLD ENERGIES', e_old(MC_UNI), cr_iscat(isel(1)), cr_iscat(isel(2))
   loop_corr: DO ic = 1, chem_ncor 
      CALL chem_neighbour_multi(isel(ia), ic, iatom, patom, natom, &
                                ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)                                                
!QWwrite(*,*) ' central ', isel(ia), cr_iscat(isel(ia)), ic
!QWwrite(*,*) ' neigb   ', iatom(0:natom(1),1), cr_iscat(iatom(0:natom(1),1)), ' ::', ncent

      loop_cent: DO icent = 1, ncent 
         IF (natom (icent)  > 0) THEN 
!write(*,*) ' AT ATOM ', ia, isel(ia), ic, iatom(0:natom(1),icent)
!                                                                       
!------ ------- Calc old energy for those energies that are affected by 
!               the move                                                
!                                                                       
            IF (mmc_cor_energy (ic, MC_OCC) ) THEN 
!                                                                       
!     ------- Occupation Correlation, only if both atoms are different  
!     ------- and the move is switch chemistry                          
!                                                                       
               IF (mmc_move == MC_MOVE_SWCHEM) THEN 
                  IF (cr_iscat (isel (1) )  /= cr_iscat (isel (2) ) ) THEN
                     e_cur (MC_OCC) = e_cur (MC_OCC) + mmc_energy_occ ( &
                     isel, ia, ic, iatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L) 
                     valid_all = valid_all.OR.valid_e 
                  ENDIF 
               ENDIF 
            ENDIF 
!                                                                       
!     ------- Unidirectional Occupation Correlation, only if both atoms are different  
!     ------- and the move is switch chemistry                          
!                                                                       
            IF (mmc_cor_energy (ic, MC_UNI) ) THEN 
               IF (mmc_move == MC_MOVE_SWCHEM) THEN 
                  IF (cr_iscat (isel (1) )  /= cr_iscat (isel (2) ) ) THEN                                                  
                     e_cur (MC_UNI) = e_cur (MC_UNI) + mmc_energy_uni ( &
                     isel, ia, ic, iatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L) 
                     valid_all = valid_all.OR.valid_e 
                  ENDIF 
               ENDIF 
            ENDIF 
!
!     ------- Coordination number
!
            IF(mmc_cor_energy(ic, MC_COORDNUM) ) THEN 
               e_cur(MC_COORDNUM) = e_cur(MC_COORDNUM) + mmc_energy_cn ( &
               isel, ia, ic, iatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)       
               valid_all = valid_all.OR.valid_e 
            ENDIF 
!                                                                       
!     ------- Displacement correlation                                  
!                                                                       
            IF (mmc_cor_energy (ic, MC_DISP) ) THEN 
            DO i = 1, 3 
            v (i) = cr_pos (i, isel (ia) ) - chem_ave_pos (i, is (ia) ) &
            - REAL(iz (ia, i) - 1) - cr_dim0 (i, 1)                   
!              u(i) = v(i)
!           v (i) = v (i) - disp (i, 0, ia) 
            ENDDO 
                  IF (chem_ldall (ic) ) THEN 
                     DO i = 1, 3 
                     jdir (i) = v (i) 
                     ENDDO 
                     rdj (ic) = skalpro (jdir, jdir, cr_gten) 
                     IF (rdj (ic)  > 0.0) THEN 
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
                  e_cur (MC_DISP) = e_cur (MC_DISP) + mmc_energy_dis (  &
                  isel, ia, ic, iatom, icent, natom, jdir, delta,&
                  rdj, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                          
                  valid_all = valid_all.OR.valid_e 
               valid_all = .TRUE. 
            ENDIF 
!                                                                       
!     ------- Displacement (Hooke's law) '                              
!                                                                       
            IF (mmc_cor_energy (ic, MC_SPRING) ) THEN 
               e_cur (MC_SPRING) = e_cur (MC_SPRING) + mmc_energy_spr ( &
               isel, ia, ic, iatom, patom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
               valid_all = valid_all.OR.valid_e 
            ENDIF 
!                                                                       
!     ------- Angular energy                                            
!                                                                       
            IF (mmc_cor_energy (ic, MC_ANGLE) ) THEN 
               e_cur (MC_ANGLE) = e_cur (MC_ANGLE) + mmc_energy_angle ( &
               ic, iatom, patom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)           
               valid_all = valid_all.OR.valid_e 
            ENDIF 
!                                                                       
!     ------- Vector energy, i.e. directional bond length               
!                                                                       
!           IF (mmc_cor_energy (ic, MC_VECTOR) ) THEN 
!              e_cur (MC_VECTOR) = e_cur (MC_VECTOR) + mmc_energy_vec ( &
!              isel, ia, ic, iatom, patom, icent, natom, valid_e)       
!              valid_all = valid_all.OR.valid_e 
!           ENDIF 
!                                                                       
!     ------- Bond length energy, non - directional                     
!                                                                       
            IF (mmc_cor_energy (ic, MC_BLEN) ) THEN 
               CONTINUE 
            ENDIF 
!                                                                       
!     ------- Lennard-Jones Potential                                   
!                                                                       
            IF (mmc_cor_energy (ic, MC_LENNARD) ) THEN 
               e_cur (MC_LENNARD) = e_cur (MC_LENNARD) + mmc_energy_len &
               (isel, ia, ic, iatom, patom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
               valid_all = valid_all.OR.valid_e 
            ENDIF 
!
!     ------- Repulsive     Potential                                   
!
            IF (mmc_cor_energy (ic, MC_REPULSIVE) ) THEN 
               e_cur (MC_REPULSIVE) = e_cur (MC_REPULSIVE) +            &
                     mmc_energy_rep (isel, ia, ic, iatom, patom, icent, &
                                     natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
               valid_all = valid_all.OR.valid_e 
            ENDIF 
!                                                                       
!     ------- Buckingham    Potential                                   
!                                                                       
            IF (mmc_cor_energy (ic, MC_BUCKING) ) THEN 
               e_cur (MC_BUCKING) = e_cur (MC_BUCKING) +                &
               mmc_energy_buck (isel, ia, ic, iatom, patom, icent,      &
               natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
               valid_all = valid_all.OR.valid_e 
            ENDIF 
         ENDIF 
      ENDDO loop_cent
   ENDDO loop_corr
ENDDO loop_natoms
!
END SUBROUTINE mmc_energies
!
!*******************************************************************************
!
SUBROUTINE mmc_modify(isel, posz, posz2, disp, MMC_MAX_ATOM_L)
!                                                                       
!     ----Perform the modifications of the atoms                        
!                                                                       
USE crystal_mod
USE chem_mod
USE celltoindex_mod
USE mc_mod
USE mmc_mod
!
USE lib_random_func
!
IMPLICIT NONE
!
!SAVE
!
INTEGER                                                 , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(IN) :: isel !(chem_max_atom) 
REAL   , DIMENSION(3)             , INTENT(OUT) :: posz
REAL   , DIMENSION(3)             , INTENT(OUT) :: posz2
REAL(KIND=PREC_SP), DIMENSION(3,0:MMC_MAX_ATOM_L,2), INTENT(INOUT) :: disp
!
INTEGER :: i, j                 ! Dummy loop indices
INTEGER :: iscat
INTEGER, DIMENSION(3) :: iz1
INTEGER, DIMENSION(3) :: iz2
INTEGER, DIMENSION(2) :: is
REAL :: disp1, disp2 
REAL :: rrrr
!
!
IF (mmc_move == MC_MOVE_DISP) THEN 
!                                                                       
!     ------Displace atoms from current positions                       
!                                                                       
!                                                                       
!     ------Modify the central atom                                     
!                                                                       
   IF(mo_maxmove(4, cr_iscat(isel(1)))==0.0) THEN
      DO i = 1, 3 
         disp(i, 0, 1)      = gasdev(DBLE(mo_maxmove(i, cr_iscat(isel(1)) )))
         posz(i)            = cr_pos(i,isel(1)) 
         cr_pos(i, isel(1)) = cr_pos(i, isel(1)) + disp(i, 0, 1) 
      ENDDO 
   ELSE
!     -- Move along a vector direction
      rrrr = gasdev(DBLE(mo_maxmove(4, cr_iscat(isel(1)) )))
      DO i = 1, 3 
         disp(i, 0, 1) = rrrr* (mo_maxmove(i, cr_iscat(isel(1)) ))
         posz (i) = cr_pos (i, isel(1)) 
         cr_pos (i, isel(1)) = cr_pos (i, isel(1)) + disp (i, 0, 1) 
      ENDDO 
   ENDIF
ELSEIF (mmc_move == MC_MOVE_SWDISP) THEN 
!                                                                       
!--Switch displacement of two selected atoms and their neighbours
!                                                                       
   CALL indextocell(isel(1), iz1, is(1)) 
   CALL indextocell(isel(2), iz2, is(2)) 
!                                                                       
!     ------Modify the central atom                                     
!                                                                       
   DO j = 1, 3 
      disp1 = cr_pos(j, isel(1)) - chem_ave_pos(j, is(1))&
               - REAL(iz1 (j) - 1) - cr_dim0 (j, 1)                   
      disp2 = cr_pos(j, isel(2)) - chem_ave_pos(j, is(2))&
               - REAL(iz2 (j) - 1) - cr_dim0 (j, 1)                   
      disp(j, 0, 1) = - disp1 + disp2 
      disp(j, 0, 2) = - disp2 + disp1 
      posz(j)            = cr_pos(j, isel(1)) 
      cr_pos(j, isel(1)) = cr_pos(j, isel(1)) + disp(j, 0, 1)
      posz2(j)           = cr_pos(j, isel(2)) 
      cr_pos(j, isel(2)) = cr_pos(j, isel(2)) + disp(j, 0, 2)
   ENDDO 
ELSEIF (mmc_move == MC_MOVE_INVDISP) THEN 
!                                                                       
!-----      ------Switch displacement of a selected atom                
!                                                                       
   CALL indextocell (isel (1), iz1, is (1) ) 
   DO j = 1, 3 
      disp1 = cr_pos(j, isel(1)) - chem_ave_pos(j, is(1))&
               - REAL(iz1 (j) - 1) - cr_dim0(j, 1)                   
      disp(j, 0, 1) = - 2 * disp1 
      posz(j) = cr_pos(j, isel(1) ) 
      cr_pos(j, isel(1)) = cr_pos(j, isel(1)) + disp(j, 0, 1)
   ENDDO 
ELSEIF (mmc_move == MC_MOVE_SWCHEM) THEN 
!                                                                       
!     ------Switch the Chemistry of two selected atoms                  
!                                                                       
   iscat             = cr_iscat(isel(2)) 
   cr_iscat(isel(2)) = cr_iscat(isel(1)) 
   cr_iscat(isel(1)) = iscat 
   iscat            = cr_prop(isel(2)) 
   cr_prop(isel(2)) = cr_prop(isel(1)) 
   cr_prop(isel(1)) = iscat 
!                                                                       
!--End of Modification of atoms according to different moves
!                                                                       
ENDIF 
!
!DEALLOCATE(disp)
!
END SUBROUTINE mmc_modify
!
!*******************************************************************************
!
SUBROUTINE mmc_unmodify(isel, posz, posz2, MMC_MAX_ATOM_L)
!                                                                       
!     ----Perform the modifications of the atoms                        
!                                                                       
USE crystal_mod
USE chem_mod
USE celltoindex_mod
USE mc_mod
USE mmc_mod
!
USE lib_random_func
!
IMPLICIT NONE
!
INTEGER                           , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(IN) :: isel !(chem_max_atom) 
REAL   , DIMENSION(3)             , INTENT(IN) :: posz
REAL   , DIMENSION(3)             , INTENT(IN) :: posz2
!
INTEGER :: i                 ! Dummy loop indices
INTEGER :: iscat
!
IF (mmc_move == MC_MOVE_DISP) THEN 
!                                                                       
!     ------Displace atoms from current positions                       
!                                                                       
   DO i = 1, 3 
      cr_pos(i, isel(1)) = posz(i)
   ENDDO 
ELSEIF (mmc_move == MC_MOVE_SWDISP) THEN 
   DO i = 1, 3 
      cr_pos(i, isel(1)) = posz(i)
   ENDDO 
   DO i = 1, 3 
      cr_pos(i, isel(2)) = posz2(i)
   ENDDO 
ELSEIF (mmc_move == MC_MOVE_INVDISP) THEN 
   DO i = 1, 3 
      cr_pos(i, isel(1)) = posz(i)
   ENDDO 
ELSEIF (mmc_move == MC_MOVE_SWCHEM) THEN 
!                                                                       
!     ------Switch the Chemistry of two selected atoms                  
!                                                                       
   iscat             = cr_iscat(isel(2)) 
   cr_iscat(isel(2)) = cr_iscat(isel(1)) 
   cr_iscat(isel(1)) = iscat 
   iscat            = cr_prop(isel(2)) 
   cr_prop(isel(2)) = cr_prop(isel(1)) 
   cr_prop(isel(1)) = iscat 
!                                                                       
!--End of Modification of atoms according to different moves
!                                                                       
ENDIF 
!
END SUBROUTINE mmc_unmodify
!
!*****7*****************************************************************
!
SUBROUTINE mmc_test_multi (iacc_good, iacc_neut, iacc_bad, &
           e_new, e_old, laccept)                                                          
!+                                                                      
!     Tests performed MC move                                           
!-                                                                      
USE mc_mod 
USE mmc_mod 
USE lib_random_func
USE random_mod
!                                                                       
IMPLICIT none 
!SAVE
!                                                                       
INTEGER, INTENT(INOUT) :: iacc_good, iacc_neut, iacc_bad 
REAL,    INTENT(IN)    :: e_old (0:MC_N_ENERGY) 
REAL,    INTENT(IN)    :: e_new (0:MC_N_ENERGY) 
LOGICAL, INTENT(OUT)   :: laccept 
!                                                                       
INTEGER :: i 
REAL    :: e_del 
REAL    :: e_ran 
REAL    :: e_delta 
REAL    :: r1
!                                                                       
!                                                                       
e_del = 0.0
DO i = 1, MC_N_ENERGY 
   IF (mmc_cor_energy(0, i) ) THEN 
      e_delta = NINT((e_new (i) - e_old (i) ) *1.0E4)*1.0E-4
      e_del = e_del + e_delta 
!     IF (e_delta > 0) THEN 
!        e_aver_p (i) = e_aver_p (i) + e_delta 
!        n_e_av_p (i) = n_e_av_p (i) + 1 
!     ELSEIF (e_delta < 0) THEN 
!        e_aver_m (i) = e_aver_m (i) + e_delta 
!        n_e_av_m (i) = n_e_av_m (i) + 1 
!     ELSE 
!        n_e_av_z (i) = n_e_av_z (i) + 1 
!     ENDIF 
   ENDIF 
ENDDO 
IF (e_del <  0) THEN 
   laccept = .TRUE. 
!     ELSEIF(e_del == 0) THEN                                           
!       laccept = .TRUE.                                                
ELSE 
   IF (mo_kt <  1.0e-10) THEN 
      laccept = .FALSE. 
   ELSE 
      e_ran = exp ( - e_del / mo_kt) 
      e_ran = e_ran / (1 + e_ran) 
      CALL RANDOM_NUMBER(r1)
      laccept = (e_ran > r1          ) 
!     laccept = (e_ran > ran1 (idum) ) 
   ENDIF 
ENDIF 
!                                                                       
IF (laccept) THEN 
   IF (e_del <  0.0) THEN 
      iacc_good = iacc_good+1 
!write(*,*) ' ACCEPTED ?', laccept, ' Good'
   ELSEIF(e_del==0) THEN
      iacc_neut = iacc_neut + 1
!write(*,*) ' ACCEPTED ?', laccept, ' Neutral'
   ELSE 
      iacc_bad = iacc_bad+1 
!write(*,*) ' ACCEPTED ?', laccept, ' Bad '
   ENDIF 
   DO i = 1, MC_N_ENERGY 
      IF (mmc_cor_energy(0, i) ) THEN 
         e_delta = NINT((e_new (i) - e_old (i) ) *1.0E4)*1.0E-4
         IF (e_delta > 0.0) THEN 
            e_aver_p (i) = e_aver_p (i) + e_delta 
            n_e_av_p (i) = n_e_av_p (i) + 1 
         ELSEIF (e_delta < 0.0) THEN 
            e_aver_m (i) = e_aver_m (i) + e_delta 
            n_e_av_m (i) = n_e_av_m (i) + 1 
         ELSE 
            n_e_av_z (i) = n_e_av_z (i) + 1 
         ENDIF 
      ENDIF 
   ENDDO
ENDIF 
!                                                                       
END SUBROUTINE mmc_test_multi                 
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_occ(isel, ia, ic, iatom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
!+                                                                      
!     Calculates the energy for chemical disorder                       
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE modify_func_mod
USE param_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(IN) :: isel
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ic
!                                                                       
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: iatom
INTEGER,                                                 INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L)                    , INTENT(IN) :: natom
LOGICAL,                                                 INTENT(OUT) :: valid_e 
!                                                                       
INTEGER :: is, js, ind
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
INTEGER :: ival1
!                                                                       
mmc_energy_occ = 0.0 
ncalc = 0 
valid_e = .FALSE. 
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
   IF (natom (icent)  /= 0) THEN 
      IF (isel (ia)  == iatom (0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         in_a = 1 
         in_e = natom (icent) 
         is = cr_iscat (iatom (0, icent) ) 
         DO ind = in_a, in_e 
            IF (check_select_status (iatom (ind, icent),  &
                                     .TRUE., cr_prop (iatom (ind,  &
                                     icent) ), cr_sel_prop) ) THEN                         
               ival1 = 0 
               js = cr_iscat (iatom (ind, icent) ) 
               ival1 = sign(1,mmc_pair(ic, MC_OCC,is,js))
               mmc_energy_occ = mmc_energy_occ +                  &
                                mmc_depth (ic,MC_OCC, is, js) * ival1
            ENDIF 
         ENDDO 
      ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
         in_a = 0 
         in_e = 0 
         is = cr_iscat (isel (ia) ) 
         ind = 0
         IF (check_select_status (iatom (ind, icent),  &
                                  .TRUE., cr_prop (iatom (ind,  &
                                  icent) ), cr_sel_prop) ) THEN                         
            ival1 = 0 
            js = cr_iscat (iatom (ind, icent) ) 
            ival1 = mmc_pair(ic, MC_OCC,is,js)
            mmc_energy_occ = mmc_energy_occ +                  &
                             mmc_depth (ic,MC_OCC, is, js) * ival1
         ENDIF 
      ENDIF 
   ENDIF 
ENDIF 
valid_e = .TRUE. 
!                                                                       
END FUNCTION mmc_energy_occ                   
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_uni(isel, ia, ic, iatom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!+                                                                      
!     Calculates the energy for unidirectional chemical disorder                       
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE modify_func_mod
USE param_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MAX_ATOM_ENV_L), INTENT(IN) :: isel
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ic
!                                                                       
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: iatom
INTEGER,                                                 INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L)                    , INTENT(IN) :: natom
LOGICAL,                                                 INTENT(OUT) :: valid_e 
!                                                                       
!                                                                       
INTEGER :: is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
INTEGER :: ival1
!                                                                       
mmc_energy_uni = 0.0 
ncalc   = 0 
valid_e = .FALSE. 
in = 1
is = cr_iscat (iatom (0, icent) )
js = cr_iscat (iatom (in, icent) )
!write(*,*) ' IN ENERGY UNI   ', iatom (0, icent), iatom (in, icent), is,js, mmc_pair(ic, MC_UNI, is,js), isel(ia)
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
   IF(natom(icent) /= 0) THEN 
      IF(isel(ia) == iatom (0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         in_a = 1 
         in_e = natom (icent) 
         is = cr_iscat (iatom (0, icent) ) 
         DO in = in_a, in_e 
            IF (check_select_status (iatom (in, icent),  &
                                     .TRUE., cr_prop (iatom (in,  &
                                     icent) ), cr_sel_prop) ) THEN                         
               IF(mmc_pair(ic, MC_UNI,is,js) /= 0) THEN
                  ival1 = 0 
                  js = cr_iscat (iatom (in, icent) ) 
                  ival1 = sign(1,mmc_pair(ic, MC_UNI,is,js))
                  mmc_energy_uni = mmc_energy_uni +                  &
                                   mmc_depth (ic,MC_UNI,  0,  0) * ival1
!write(*,*) ' Selected central',  is, js, ival1, mmc_depth (ic,MC_UNI, is, js) * ival1
               ENDIF 
            ENDIF 
         ENDDO 
      ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
         in_a = 0 
         in_e = 0 
         is = cr_iscat (isel (ia) ) 
         in = 0
            IF (check_select_status (iatom (in, icent),  &
                                     .TRUE., cr_prop (iatom (in,  &
                                     icent) ), cr_sel_prop) ) THEN                         
                  ival1 = 0 
                  js = cr_iscat (iatom (in, icent) ) 
               IF(mmc_pair(ic, MC_UNI,js,is) /= 0) THEN
                  ival1 = mmc_pair(ic, MC_UNI,js,is)
                  mmc_energy_uni = mmc_energy_uni +                  &
                                   mmc_depth (ic,MC_UNI, 0, 0) * ival1
!write(*,*) ' Selected neighbo', is, js, ival1, mmc_depth (ic,MC_UNI,js,is) * ival1
!write(*,*) ' Unidirec neighbo', js, is, ival1, mmc_depth (ic,MC_UNI,js,is) * ival1
               ENDIF 
            ENDIF 
!              
      ENDIF 
   ENDIF 
ENDIF 
valid_e = .TRUE. 
!write(*,*) ' CALCULATED ENERGY ', mmc_energy_uni
!read(*,*) is
!                                                                       
END FUNCTION mmc_energy_uni                   
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_occ_mol(ianz, imol, amol, valid_e, MAXMOL) 
!+                                                                      
!     Calculates the energy for occupational disorder for               
!     molecules.                                                        
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE chem_menu
USE mc_mod 
USE mmc_mod 
USE molecule_mod 
USE rmc_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ianz 
INTEGER, INTENT(IN) :: MAXMOL
INTEGER, DIMENSION(ianz), INTENT(IN) ::  imol !(ianz) 
INTEGER, DIMENSION(ianz), INTENT(IN) ::  amol !(ianz) 
LOGICAL                 , INTENT(OUT) :: valid_e 
!                                                                       
INTEGER, DIMENSION(0:MAXMOL) :: ineig
INTEGER :: nneig 
INTEGER :: ic, in, ia, ival1, ival2 
!                                                                       
      mmc_energy_occ_mol = 0.0 
      valid_e = .TRUE. 
!                                                                       
!------ Loop over all modified atoms                                    
!                                                                       
      DO ia = 1, ianz 
      ival1 = - 1 
      IF (mole_type (imol (ia) )  == amol (1) ) ival1 = 1 
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
      CALL chem_neighbour_mol (imol (ia), ic, ineig, nneig, maxmol) 
      IF (nneig /= 0) THEN 
         DO in = 1, nneig 
         ival2 = - 1 
         IF (mole_type (ineig (in) )  == amol (1) ) ival2 = 1 
         mmc_energy_occ_mol = mmc_energy_occ_mol + mmc_const (ic,       &
         MC_OCC) * ival1 * ival2                                        
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (.NOT.valid_e) THEN 
         mmc_energy_occ_mol = 0.0 
      ENDIF 
!                                                                       
END FUNCTION mmc_energy_occ_mol               
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_cn(isel, ia, ic, iatom, icent, natom, valid_e,         &
                            MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
!+                                                                      
!     Calculates the energy for coordination number
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE modify_func_mod
USE param_mod 
!                                                                       
IMPLICIT none 
!                                                                       
!
INTEGER                          , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                          , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                          , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM), INTENT(IN) :: isel
INTEGER                          , INTENT(IN)  ::  ia
INTEGER                          , INTENT(IN)  ::  ic
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT), INTENT(IN) :: iatom
INTEGER                                             , INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L), INTENT(IN)                  :: natom
LOGICAL, INTENT(OUT) :: valid_e 
!                                                                       
INTEGER, DIMENSION(:), ALLOCATABLE           :: n_neig   ! Number of neighbors of type js
INTEGER :: is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
!
mmc_energy_cn = 0.0
ncalc = 0 
valid_e = .FALSE. 
IF(ALLOCATED(n_neig)) DEALLOCATE(n_neig)
ALLOCATE(n_neig(0:CR_NSCAT))
n_neig(:) = 0
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
   IF(natom (icent) /= 0) THEN 
      IF(isel(ia) == iatom(0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         in_a = 1 
         in_e = natom(icent) 
         is   = cr_iscat(iatom (0, icent) ) 
         js   = -1                            ! Flag negative to indicate no neighbor
         DO in = in_a, in_e 
            IF(check_select_status(iatom(in, icent),  .TRUE.,          &
                                   cr_prop (iatom (in, icent)), cr_sel_prop)) THEN
               js    = cr_iscat(iatom (in, icent) ) 
               IF(mmc_pair(ic, MC_COORDNUM, is, js)==-1) THEN
                  n_neig(js) = n_neig(js) + 1
               ENDIF
            ENDIF
         ENDDO 
         DO js = 0, cr_nscat
            IF(mmc_pair(ic, MC_COORDNUM, is, js) == -1) THEN
               mmc_energy_cn = mmc_energy_cn +                                 &
                               mmc_depth (ic,MC_COORDNUM, is, js) *            &
                               (n_neig(js) - mmc_target_corr(ic, MC_COORDNUM, is, js))**2
            ENDIF 
         ENDDO 
      ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
!        in_a = 0 
!        in_e = 0 
!        is = cr_iscat (isel (ia) ) 
!        in = 0
!        IF (check_select_status (iatom (in, icent),  &
!                                 .TRUE., cr_prop (iatom (in,  &
!                                 icent) ), cr_sel_prop) ) THEN                         
!           ival1 = 0 
!           js = cr_iscat (iatom (in, icent) ) 
!           ival1 = mmc_pair(ic, MC_OCC,is,js)
!           mmc_energy_cn = mmc_energy_cn +                  &
!                            mmc_depth (ic,MC_OCC, 0, 0) * ival1
!        ENDIF 
!              
      ENDIF 
   ENDIF 
ENDIF 
valid_e = .TRUE. 
DEALLOCATE(n_neig)
!                                                                       
END FUNCTION mmc_energy_cn                   
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_dis (isel, ia, ic, iatom, icent,  &
      natom, jdir, delta, rdj, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L,&
      MMC_MAX_ATOM_L)                                 
!+                                                                      
!     Calculates the energy for chemical disorder                       
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE celltoindex_mod
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
!
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER                          , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                          , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                          , INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L) , INTENT(IN) :: isel
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ic
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: iatom
INTEGER, INTENT(IN) :: icent 
INTEGER , DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL , INTENT(OUT) ::valid_e 
REAL(KIND=PREC_SP), DIMENSION(3), INTENT(IN) :: jdir
REAL(KIND=PREC_SP),               INTENT(IN) ::  delta 
REAL(KIND=PREC_SP), DIMENSION(CHEM_MAX_COR), INTENT(IN) :: rdj
!                                                                       
INTEGER :: i, is, js, in, jjs 
INTEGER :: in_a, in_e 
INTEGER :: cell (3), site 
REAL :: u (3)
!                                                                       
REAL :: dx 
!                                                                       
mmc_energy_dis = 0.0 
valid_e = .FALSE. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .OR.  &
          chem_ctyp (ic) == CHEM_ENVIR .OR.  &
          chem_ctyp (ic) == CHEM_RANGE .OR.  &
          chem_ctyp (ic) == CHEM_DIST  .OR.  &
          chem_ctyp (ic) == CHEM_CON        ) THEN                                               
!                                                                       
         IF (natom (icent)  /= 0) THEN 
            IF (isel (ia)  == iatom (0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
               DO jjs = 0, cr_nscat 
               IF (mmc_pair (ic, MC_DISP, is, jjs) == -1 ) THEN 
                  DO in = in_a, in_e 
                  js = cr_iscat (iatom (in, icent) ) 
                  IF (is == js.OR.mmc_pair (ic, MC_DISP, is, js) == -1 ) THEN 
                     IF (check_select_status (iatom (in, icent),  &
                                           .TRUE., cr_prop (iatom (  &
                     in, icent) ), cr_sel_prop) ) THEN                  
                        CALL indextocell (iatom (in, icent), cell, site) 
                        DO i = 1, 3 
                        u (i) = cr_pos (i, iatom (in, icent) ) -        &
                        chem_ave_pos (i, site) - REAL(cell (i)        &
                        - 1) - cr_dim0 (i, 1)                           
                        ENDDO 
                        dx = skalpro (u, jdir, cr_gten) / rdj (ic) 
!                                                                       
                        mmc_energy_dis = mmc_energy_dis + mmc_depth (ic,&
                        MC_DISP, 0, 0) * delta * dx                     
                        valid_e = .TRUE. 
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
               IF (mmc_pair (ic, MC_DISP, is, jjs) == -1 ) THEN 
                  in = 0 
                  js = cr_iscat (iatom (in, icent) ) 
                  IF (is == js.OR.mmc_pair (ic, MC_DISP, is, js) == -1 ) THEN 
                     IF (check_select_status (iatom (in, icent),  &
                                           .TRUE., cr_prop (iatom (  &
                     in, icent) ), cr_sel_prop) ) THEN                  
                        CALL indextocell (iatom (in, icent), cell, site) 
                        DO i = 1, 3 
                        u (i) = cr_pos (i, iatom (in, icent) ) -        &
                        chem_ave_pos (i, site) - REAL(cell (i)        &
                        - 1) - cr_dim0 (i, 1)                           
                        ENDDO 
                        dx = skalpro (u, jdir, cr_gten) / rdj (ic) 
!                                                                       
                        mmc_energy_dis = mmc_energy_dis + mmc_depth (ic,&
                        MC_DISP, 0, 0) * delta * dx                     
                        valid_e = .TRUE. 
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
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_spr (isel, ia, ic, iatom, patom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                  
!+                                                                      
!     Calculates the energy for distortions according to                
!                                                                       
!       E = SUM k*(d-d0)**2                                             
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(IN) :: isel
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ic
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: patom
INTEGER                                                   , INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT)                         , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL    :: d, u (3), v (3) 
!                                                                       
!                                                                       
      mmc_energy_spr = 0.0 
      ncalc = 0 
      valid_e = .FALSE. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .OR.  &
          chem_ctyp (ic) == CHEM_ENVIR .OR.  &
          chem_ctyp (ic) == CHEM_RANGE .OR.  &
          chem_ctyp (ic) == CHEM_DIST  .OR.  &
          chem_ctyp (ic) == CHEM_CON        ) THEN                                               
!                                                                       
         IF (natom (icent)  /= 0) THEN 
            IF (isel (ia)  == iatom (0, icent) ) THEN 
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
            IF (mmc_target_corr (ic, MC_SPRING, is, js)  /= 0.0) THEN 
               IF (check_select_status (iatom (in, icent),  &
                                           .TRUE., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) THEN                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.TRUE., u, v) 
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
      IF (ncalc > 0) THEN 
         mmc_energy_spr = mmc_energy_spr / REAL(ncalc) 
         valid_e = .TRUE. 
      ELSE 
         mmc_energy_spr = 0.0 
         valid_e = .FALSE. 
      ENDIF 
!                                                                       
END FUNCTION mmc_energy_spr                   
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_spr_mol (nmol, imol, valid_e, MAX_ATOM_ENV_L) 
!+                                                                      
!     Calculates the energy for distortions for molecules ..            
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE chem_menu
USE celltoindex_mod
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE molecule_mod 
USE rmc_mod 
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: nmol
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, DIMENSION(2), INTENT(IN) :: imol
LOGICAL, INTENT(OUT) :: valid_e 
!                                                                       
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L) :: ineig 
INTEGER :: nneig
INTEGER ::  iatom, jatom 
INTEGER ::  cell (3), site 
INTEGER  :: i, is, js, ic, in, ia 
REAL  :: d, u (3), v (3) 
!                                                                       
      mmc_energy_spr_mol = 0.0 
      valid_e = .FALSE. 
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
      u (i) = cr_pos (i, iatom) - chem_ave_pos (i, site) - REAL(cell (&
      i) - 1) - cr_dim0 (i, 1)                                          
      ENDDO 
      d = do_blen (.TRUE., u, u) 
!                                                                       
      mmc_energy_spr_mol = mmc_energy_spr_mol + mmc_const (0, MC_SPRING)&
      * d**2                                                            
!                                                                       
!------ - Loop over all defined interactions                            
!                                                                       
      DO ic = 1, chem_ncor 
      IF (chem_ctyp (ic)  == CHEM_VEC.OR.chem_ctyp (ic)  == CHEM_ENVIR) &
      THEN                                                              
!                                                                       
!------ --- Get the neighbours                                          
!                                                                       
         CALL chem_neighbour_mol (imol (ia), ic, ineig, nneig, MAX_ATOM_ENV_L) 
         IF (nneig /= 0) THEN 
            DO in = 1, nneig 
            is = mole_type (imol (ia) ) 
            js = mole_type (ineig (in) ) 
            jatom = mole_cont (mole_off (ineig (in) ) + 1) 
            IF (mmc_target_corr (ic, MC_SPRING, is, js)  /= 0.0) THEN 
               DO i = 1, 3 
               u (i) = cr_pos (i, iatom) 
               v (i) = cr_pos (i, jatom) 
               ENDDO 
               d = do_blen (.TRUE., u, v) 
!                                                                       
               mmc_energy_spr_mol = mmc_energy_spr_mol + mmc_const (ic, &
               MC_SPRING) * (d-mmc_target_corr (ic, MC_SPRING, is, js) )&
               **2                                                      
!                                                                       
               valid_e = .TRUE. 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      IF (.NOT.valid_e) THEN 
         mmc_energy_spr_mol = 0.0 
      ENDIF 
!                                                                       
END FUNCTION mmc_energy_spr_mol               
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_len (isel, ia, ic, iatom, patom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!+                                                                      
!     Calculates the energy for distortions according to a Lennard-Jones
!     potential                                                         
!                                                                       
!       E = SUM A/d**12 - b/d**6                                        
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(IN) ::  isel !(chem_max_atom) 
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: ia
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) :: valid_e 
!                                                                       
      INTEGER :: i, is, js, in
      INTEGER :: in_a, in_e 
      INTEGER :: ncalc 
      REAL :: d, u (3), v (3) 
!                                                                       
!QWwrite(*,*) ' DIMENSIONS  ', MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L
      mmc_energy_len = 0.0 
      ncalc = 0 
      valid_e = .FALSE. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .OR.  &
          chem_ctyp (ic) == CHEM_ENVIR .OR.  &
          chem_ctyp (ic) == CHEM_RANGE .OR.  &
          chem_ctyp (ic) == CHEM_DIST  .OR.  &
          chem_ctyp (ic) == CHEM_CON        ) THEN                                               
!                                                                       
         IF (natom (icent)  /= 0) THEN 
            IF (isel (ia)  == iatom (0, icent) ) THEN 
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
            IF (mmc_target_corr (ic, MC_LENNARD, is, js)  /= 0.0) THEN 
               IF (check_select_status (iatom (in, icent),  &
                                           .TRUE., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) THEN                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.TRUE., u, v) 
!                                                                       
!QWwrite(*,*) ' ENERGY ', is, js, d,  mmc_target_corr (ic, MC_LENNARD, is, js),&
!QW                  d-mmc_target_corr (ic, MC_LENNARD, is, js), &
!QW                                                    mmc_depth (ic,      &
!QW                  MC_LENNARD, is, js) * (mmc_len_a (ic, is, js) / d**   &
!QW                  mmc_len_m (ic, is, js) - mmc_len_b (ic, is, js)       &
!QW                  / d**mmc_len_n (ic, is, js) )
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
      IF (ncalc > 0) THEN 
         mmc_energy_len = mmc_energy_len / REAL(ncalc) 
         valid_e = .TRUE. 
      ELSE 
         mmc_energy_len = 0.0 
         valid_e = .FALSE. 
      ENDIF 
!                                                                       
END FUNCTION mmc_energy_len                   
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_rep (isel, ia, ic, iatom, patom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!+                                                                      
!     Calculates the energy for distortions according to a Repulsive
!     potential                                                         
!                                                                       
!       E = 1/d**n + depth
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MAX_ATOM_ENV_L) , INTENT(IN) :: isel
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ic
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL  :: d, u (3), v (3) 
!                                                                       
      mmc_energy_rep = 0.0 
      ncalc = 0 
      valid_e = .FALSE. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .OR.  &
          chem_ctyp (ic) == CHEM_ENVIR .OR.  &
          chem_ctyp (ic) == CHEM_RANGE .OR.  &
          chem_ctyp (ic) == CHEM_DIST  .OR.  &
          chem_ctyp (ic) == CHEM_CON        ) THEN                                               
!                                                                       
         IF (natom (icent)  /= 0) THEN 
            IF (isel (ia)  == iatom (0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
               DO i = 1, 3 
               u (i) = patom (i, 0, icent) 
               ENDDO 
               in_a = 1 
               in_e = natom (icent) 
               is = cr_iscat (iatom (0, icent) ) 
!!write(*,*) 'CENT ', isel(ia), cr_iscat(isel(ia)), (cr_iscat(iatom(i, icent)),i=1, natom(icent))
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
            IF (mmc_target_corr (ic, MC_REPULSIVE, is, js)  /= 0.0) THEN 
               IF (check_select_status (iatom (in, icent),  &
                                           .TRUE., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) THEN                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.TRUE., u, v) 
!                                                                       
                  IF(d > mmc_rep_c (ic, is,js)) THEN
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
      IF (ncalc > 0) THEN 
         mmc_energy_rep = mmc_energy_rep / REAL(ncalc) 
         valid_e = .TRUE. 
      ELSE 
         mmc_energy_rep =    -1.*ABS(mmc_rep_low)
         valid_e = .TRUE. 
      ENDIF 
!write(*,*) ' Ncalc, energy_rep',ncalc, mmc_energy_rep
!                                                                       
END FUNCTION mmc_energy_rep
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_buck (isel, ia, ic, iatom, patom, icent, &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!+                                                                      
!     Calculates the energy for distortions according to a Buckingham   
!     potential                                                         
!                                                                       
!       E = SUM A * exp(-d/rho) - B/r**6                                
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(IN) :: isel ! (chem_max_atom) 
INTEGER, INTENT(IN) :: ia
INTEGER, INTENT(IN) :: ic
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL :: d, u (3), v (3) 
!                                                                       
      mmc_energy_buck = 0.0 
      ncalc = 0 
      valid_e = .FALSE. 
!                                                                       
      IF (chem_ctyp (ic) == CHEM_VEC   .OR.  &
          chem_ctyp (ic) == CHEM_ENVIR .OR.  &
          chem_ctyp (ic) == CHEM_RANGE .OR.  &
          chem_ctyp (ic) == CHEM_DIST  .OR.  &
          chem_ctyp (ic) == CHEM_CON        ) THEN                                               
!                                                                       
         IF (natom (icent)  /= 0) THEN 
            IF (isel (ia)  == iatom (0, icent) ) THEN 
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
            IF (mmc_target_corr (ic, MC_BUCKING, is, js)  /= 0.0) THEN 
               IF (check_select_status (iatom (in, icent),  &
                                           .TRUE., cr_prop (iatom (in,     &
               icent) ), cr_sel_prop) ) THEN                            
                  DO i = 1, 3 
                  v (i) = patom (i, in, icent) 
                  ENDDO 
                  d = do_blen (.TRUE., u, v) 
!                                                                       
                  IF (d > mmc_buck_rmin (ic, is, js) ) THEN 
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
      IF (ncalc > 0) THEN 
         mmc_energy_buck = mmc_energy_buck / REAL(ncalc) 
         valid_e = .TRUE. 
      ELSE 
         mmc_energy_buck = 0.0 
         valid_e = .FALSE. 
      ENDIF 
!                                                                       
END FUNCTION mmc_energy_buck                  
!
!*****7*****************************************************************
!
REAL FUNCTION mmc_energy_angle (ic, iatom, patom, icent,    &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!+                                                                      
!     Calculates the energy for angular deviations according to         
!                                                                       
!       E = SUM k*(a-a0)**2                                             
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
USE rmc_mod 
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, INTENT(IN) :: ic
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: iatom
REAL   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js
INTEGER :: ii, jj 
LOGICAL :: lnoneig 
REAL :: a, b, u (3), v (3), w (3) 
!                                                                       
      lnoneig = .TRUE. 
      mmc_energy_angle = 0.0 
      valid_e = .FALSE. 
!                                                                       
!                                                                       
      IF (chem_ctyp (ic) == CHEM_ANG   .OR.  &
          chem_ctyp (ic) == CHEM_VEC   .OR.  &
          chem_ctyp (ic) == CHEM_CON   .OR.  &
          chem_ctyp (ic) == CHEM_ENVIR      ) THEN                                               
         IF (natom (icent)  /= 0) THEN 
            DO i = 1, 3 
            u (i) = patom (i, 0, icent) 
            ENDDO 
            DO ii = 1, natom (icent) - 1 
            DO jj = ii + 1, natom (icent) 
            lnoneig = .FALSE. 
            is = cr_iscat (iatom (ii, icent) ) 
            js = cr_iscat (iatom (jj, icent) ) 
            IF (mmc_target_corr (ic, MC_ANGLE, is, js)  /= 0.0) THEN 
               IF (check_select_status (iatom (ii, icent),  &
                                           .TRUE., cr_prop (iatom (ii,     &
               icent) ), cr_sel_prop) ) THEN                            
                  IF (check_select_status (jj, .TRUE., cr_prop (iatom (jj,  &
                  icent) ), cr_sel_prop) ) THEN                         
                     DO i = 1, 3 
                     v (i) = patom (i, ii, icent) 
                     w (i) = patom (i, jj, icent) 
                     ENDDO 
                     a = do_bang (.TRUE., v, u, w) 
                     b = mmc_target_corr (ic, MC_ANGLE, is, js) 
!                                                                       
                     IF (b <= 90) THEN 
                        IF (a > b) THEN 
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
                     valid_e = .TRUE. 
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      IF (lnoneig) THEN 
         valid_e = .FALSE. 
      ELSE 
      ENDIF 
!                                                                       
      IF (.NOT.valid_e) THEN 
         mmc_energy_angle = 0.0 
      ENDIF 
      END FUNCTION mmc_energy_angle                 
!
!*****7*****************************************************************
!     REAL FUNCTION mmc_energy_vec (isel, ia, ic, iatom, patom, icent,  &
!     natom, valid_e)                                                   
!+                                                                      
!     Calculates the energy for vector distributions according to       
!                                                                       
!       E = SUM k*|(v-v0)|**2                                           
!                                                                       
!-                                                                      
!     USE crystal_mod 
!     USE chem_mod 
!     USE mc_mod 
!     USE mmc_mod 
!     IMPLICIT none 
!                                                                       
       
!     USE debug_mod 
!                                                                       
!                                                                       
!     INTEGER isel (chem_max_atom) 
!                                                                       
!     INTEGER iatom (0:maxatom, MMC_MAX_CENT) 
!     INTEGER natom (MMC_MAX_CENT) 
!     INTEGER icent 
!     LOGICAL valid_e 
!     INTEGER i, is, js, ic, in, ia, nv 
!     INTEGER in_a, in_e 
!     LOGICAL found 
!     REAL patom (3, 0:maxatom, MMC_MAX_CENT) 
!     REAL an, dd, d, u (3), v (3), w (3) 
!     REAL ddmin 
!     REAL null (3) 
!                                                                       
!     REAL do_blen 
!     REAL do_bang 
!                                                                       
!     DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
!     dbg = .FALSE. 
!     IF (dbg) THEN 
!     WRITE ( * ,  * ) ' ================energy_vec', '=================&
!    &==========='                                                      
!     ENDIF 
!     mmc_energy_vec = 0.0 
!     valid_e = .FALSE. 
!     ddmin = 1.1 * chem_rmax_env (ic) 
!     found = .FALSE. 
!                                                                       
!     IF (chem_ctyp (ic)  == CHEM_VEC.OR.chem_ctyp (ic)  == CHEM_ENVIR) &
!     THEN                                                              
!                                                                       
!DBG                                                                    
!        IF (dbg) THEN 
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
!        IF (natom (icent)  /= 0) THEN 
!           IF (isel (ia)  == iatom (0, icent) ) THEN 
!              IF (dbg) THEN 
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
!              IF (dbg) THEN 
!                 WRITE ( * , * ) ' selected = neighbour' 
!              ENDIF 
!              DO i = 1, 3 
!              u (i) = cr_pos (i, isel (ia) ) 
!              ENDDO 
!              in_a = 0 
!              in_e = 0 
!              is = cr_iscat (isel (ia) ) 
!           ENDIF 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) ' central           ', u 
!           ENDIF 
!           DO in = in_a, in_e 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' 
!           ENDIF 
!           js = cr_iscat (iatom (in, icent) ) 
!           DO nv = 1, mmc_nvec (ic, is, js) 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) ' is,js ,mmc_vec    ', is, js,  (mmc_vec (i, nv, &
!    &ic, is, js) , i = 1, 3)                                           
!           ENDIF 
!           DO i = 1, 3 
!           v (i) = patom (i, in, icent) - u (i) 
!           w (i) = mmc_vec (i, nv, ic, is, js) 
!           ENDDO 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) ' nachbar           ',  (patom (i, in, icent) , i&
!    & = 1, 3)                                                          
!           ENDIF 
!           an = do_bang (.TRUE., v, null, w) 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) ' Winkel            ', an 
!     WRITE ( * ,  * ) ' Vektor Zen->Nei   ',  (v (i) , i = 1, 3) 
!              WRITE ( * , * ) ' Differenz zu soll ', (w (i) - v (i) ,  &
!              i = 1, 3)                                                
!           ENDIF 
!           IF (an < mmc_vec (4, nv, ic, is, js) ) THEN 
!              dd = do_blen (.TRUE., null, v) 
!              IF (dbg) THEN 
!     WRITE ( * ,  * ) ' Laenge            ', dd 
!              ENDIF 
!              IF (dd < ddmin) THEN 
!                 ddmin = dd 
!                 d = do_blen (.TRUE., v, w) 
!DBG          mmc_energy_vec = mmc_energy_vec +                         
!DBG     &                                mmc_const(ic,MC_VECTOR) * d**2
!                 found = .TRUE. 
!                 IF (dbg) THEN 
!     WRITE ( * ,  * ) ' distance           ', d 
!                 ENDIF 
!              ENDIF 
!           ELSEIF (180. - an < mmc_vec (4, nv, ic, is, js) ) THEN 
!              DO i = 1, 3 
!              v (i) = - v (i) 
!              ENDDO 
!              dd = do_blen (.TRUE., null, v) 
!              IF (dbg) THEN 
!     WRITE ( * ,  * ) ' Vektor Zen->Nei   ',  (v (i) , i = 1, 3) 
!                 WRITE ( * , * ) ' Differenz zu soll ', (w (i) - v (i) &
!                 , i = 1, 3)                                           
!     WRITE ( * ,  * ) ' Laenge            ', dd 
!              ENDIF 
!              IF (dd < ddmin) THEN 
!                 ddmin = dd 
!                 d = do_blen (.TRUE., v, w) 
!DBG          mmc_energy_vec = mmc_energy_vec +                         
!DBG     &                                mmc_const(ic,MC_VECTOR) * d**2
!                 found = .TRUE. 
!                 IF (dbg) THEN 
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
!        IF (found) THEN 
!DBG          mmc_energy_vec = mmc_const(ic,MC_VECTOR) * d**2           
!           valid_e = .TRUE. 
!           CONTINUE 
!        ELSE 
!           mmc_energy_vec = 0.0 
!           valid_e = .FALSE. 
!        ENDIF 
!        IF (dbg) THEN 
!           WRITE ( * , * ) 'energy_vec ', mmc_energy_vec 
!        ENDIF 
!     ENDIF 
!     dbg = .FALSE. 
!                                                                       
!     IF (.NOT.valid_e) THEN 
!        mmc_energy_vec = 0.0 
!     ENDIF 
!     END FUNCTION mmc_energy_vec                   
!*****7*****************************************************************
!
SUBROUTINE mmc_correlations (lout, rel_cycl, done, lfinished) 
!-                                                                      
!     Determines the achieved correlations                              
!                                                                       
!+                                                                      
USE crystal_mod 
USE chem_mod 
USE chem_menu
USE chem_aver_mod
USE chem_neig_multi_mod
USE atom_env_mod
USE celltoindex_mod
USE metric_mod
USE mc_mod 
USE mmc_mod 
!
USE debug_mod 
USE errlist_mod 
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
LOGICAL , INTENT(IN) :: lout     ! Flag for output yes/no
REAL    , INTENT(IN) :: rel_cycl ! Relative progress along cycles
LOGICAL , INTENT(INOUT) :: done     ! MMC is converged/ stagnates
LOGICAL , INTENT(IN)    :: lfinished ! MMC is finished
! 
!                                                                       
CHARACTER(LEN=30) :: energy_name (0:MC_N_ENERGY) 
!                                                                       
INTEGER :: ic, je, ic_a 
INTEGER :: is, js, ls 
INTEGER :: iis, jjs, lls, iic, kk 
INTEGER :: i, j, k, l 
INTEGER :: icent 
!                                                                       
LOGICAL :: searching 
!LOGICAL   :: lfirst = .TRUE.  ! Flag to write output only at first instance
!                                                                       
INTEGER :: ncent 
REAL(KIND=PREC_SP), DIMENSION(:,:,:), ALLOCATABLE :: patom ! (3, 0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(  :,:), ALLOCATABLE :: iatom ! (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: natom ! ( MMC_MAX_CENT) 
!                                                                       
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: bl_anz ! (0:DEF_maxscat, 0:DEF_maxscat) 
REAL(KIND=PREC_SP), DIMENSION(:,:), ALLOCATABLE :: bl_sum ! (0:DEF_maxscat, 0:DEF_maxscat) 
REAL(KIND=PREC_SP), DIMENSION(:,:), ALLOCATABLE :: bl_s2  ! (0:DEF_maxscat, 0:DEF_maxscat) 
REAL :: u (3), v (3), d (3) 
REAL :: dist 
!REAL :: divisor
!
INTEGER :: n_cn
INTEGER, DIMENSION(:), ALLOCATABLE :: ncentral
INTEGER, DIMENSION(:,:), ALLOCATABLE :: p_cn
!                                                                       
INTEGER, DIMENSION(:,:), ALLOCATABLE :: pneig ! (0:DEF_MAXSCAT, 0:DEF_MAXSCAT) 
!INTEGER pair11, pair12, pair21, pair22 
INTEGER nneigh 
!REAL :: prob11=0.0, prob12, prob22 
!REAL :: thet = 0.0
REAL :: damp = 1.0
!                                                                       
REAL :: wi, wis 
REAL :: divider 
!                                                                       
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: ba_sum ! (CHEM_MAX_COR * MMC_MAX_ANGLES) 
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: ba_s2  ! (CHEM_MAX_COR * MMC_MAX_ANGLES) 
INTEGER           , DIMENSION(:), ALLOCATABLE :: ba_anz ! (CHEM_MAX_COR * MMC_MAX_ANGLES) 
!                                                                       
INTEGER :: icc (3), jcc (3) 
REAL :: idir (3), jdir (3), disi (3), disj (3) 
REAL :: rdi=1.0, rdj=1.0, dpi=1.0, dpj
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: xnn !(0:maxscat, 0:maxscat) 
REAL(KIND=PREC_SP), DIMENSION(:,:), ALLOCATABLE :: xij !(0:maxscat, 0:maxscat) 
REAL(KIND=PREC_SP), DIMENSION(:,:), ALLOCATABLE :: xi2 !(0:maxscat, 0:maxscat) 
REAL(KIND=PREC_SP), DIMENSION(:,:), ALLOCATABLE :: xj2 !(0:maxscat, 0:maxscat) 
!
!                                                                       
      DATA energy_name / 'none', 'Chemical correlation    ', 'Displaceme&
     &nt correlation', 'Distance correlation (Hooke)', 'Angular      cor&
     &relation', 'Vector       correlation', 'Distance     correlation',&
     & 'Lennard Jones potential ', 'Buckingham    potential ' ,         &
       'Repulsive     potential ', 'Coordination number     ', &
       'Unidirectional corr     '/         
!
      damp = 0.01 + 0.99*exp(-4.0*rel_cycl)
!
ALLOCATE(patom(3, 0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(iatom(0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(natom(MMC_MAX_CENT))
ALLOCATE( xnn (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( xij (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( xi2 (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( xj2 (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( ba_sum (CHEM_MAX_COR * MMC_MAX_ANGLES) )
ALLOCATE( ba_s2  (CHEM_MAX_COR * MMC_MAX_ANGLES) )
ALLOCATE( ba_anz (CHEM_MAX_COR * MMC_MAX_ANGLES) )
ALLOCATE( bl_anz (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( bl_sum (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( bl_s2 (0:MAXSCAT, 0:MAXSCAT) )
ALLOCATE( p_cn(0:MAXSCAT, 0:MAXSCAT))
ALLOCATE( pneig(0:MAXSCAT, 0:MAXSCAT) )
!
IF(.NOT.ALLOCATED(mmc_h_diff) ) THEN
   ALLOCATE(mmc_h_diff(1:mmc_h_number + 10,0:MMC_H_NNNN-1))
   ALLOCATE(mmc_h_targ(mmc_h_number+10))
   ALLOCATE(mmc_h_aver(mmc_h_number+10))
   ALLOCATE(mmc_h_maxd(mmc_h_number+10,0:MMC_H_NNNN-1))
   mmc_h_diff  = 0.0
   mmc_h_targ  = 0.0
   mmc_h_aver  = 1.0E20
   mmc_h_maxd  = 1.0E20
   mmc_h_index =  -1
   mmc_h_ncycl =  0
   mmc_h_aver  =  0.0
   mmc_h_maxd  =  1.0E20
ENDIF
mmc_h_ctarg = 0          ! Start with no targets in history
mmc_h_index = MOD(mmc_h_index+1,MMC_H_NNNN)   ! increment current index
!                                                                       
!------ Write title line                                                
!                                                                       
      IF (lout) THEN 
         WRITE (output_io, 410) 
      ENDIF 
!                                                                       
!     Get the average structure for the distance energies               
!                                                                       
!     IF (mmc_cor_energy (0, MC_DISP)    .OR.mmc_cor_energy (0, MC_SPRING) &
!     .OR.mmc_cor_energy (0, MC_LENNARD) .OR.mmc_cor_energy (0, MC_BUCKING)&
!     .OR.mmc_cor_energy (0,MC_REPULSIVE) ) THEN                                                
      IF (mmc_cor_energy (0, MC_DISP)                                      &
         .OR. mmc_move_prob(MC_MOVE_SWDISP) > 0                            &
         .OR. mmc_move_prob(MC_MOVE_INVDISP) > 0                           &
         ) THEN
         CALL chem_aver (.FALSE., .TRUE.) 
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
ALLOCATE(ncentral(0:MAXSCAT))
!                                                                       
!     Loop over all correlations                                        
!                                                                       
main_corr: DO ic = 1, chem_ncor 
   ncentral(:) = 0
   IF (lout) THEN 
      WRITE (output_io, * ) 
   ENDIF 
   DO is = 0, MAXSCAT 
      DO js = 0, MAXSCAT 
         pneig (is, js) = 0 
         p_cn  (is, js) = 0 
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
   IF (mmc_cor_energy (0, MC_DISP) ) THEN 
      DO i = 1, 3 
         idir (i) = chem_dir (i, 1, ic) 
         jdir (i) = chem_dir (i, 2, ic) 
      ENDDO 
!                                                                       
!------ calculate correlations                                          
!                                                                       
      rdi = skalpro (idir, idir, cr_gten) 
      rdj = skalpro (jdir, jdir, cr_gten) 
      IF (rdi > 0.0) rdi = sqrt (rdi) 
      IF (rdj > 0.0) rdj = sqrt (rdj) 
   ENDIF 
!                                                                       
!     -- Loop over all atoms                                            
!                                                                       
main_atoms: DO i = 1, cr_natoms 
      is = cr_iscat (i) 
      ncentral(is) = ncentral(is) + 1
      CALL chem_neighbour_multi (i, ic, iatom, patom, natom, ncent, MAX_ATOM_ENV, MMC_MAX_CENT)
      IF (ier_num /= 0) THEN
         RETURN 
      ENDIF
!                                                                       
!------ ---- In case of Displacement correlation, calculate             
!            displacement of central atom                                                 
!                                                                       
is_mc_disp: IF (mmc_cor_energy (ic, MC_DISP) ) THEN 
         CALL indextocell (i, icc, is) 
         DO j = 1, 3 
            disi (j) = cr_pos (j, i) - chem_ave_pos (j, is) - &
                       REAL(icc (j) - 1) - cr_dim0 (j, 1)                                          
         ENDDO 
!                                                                       
      IF (chem_ldall (ic) ) THEN 
         DO j = 1, 3 
            jdir (j) = disi (j) 
         ENDDO 
         rdj = skalpro (jdir, jdir, cr_gten) 
         IF (rdj > 0.0) THEN 
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
is_cent:  IF (i == iatom (0, icent) ) THEN 
is_energy:  IF (mmc_cor_energy (ic, MC_OCC)        .OR. &
                mmc_cor_energy (ic, MC_UNI)        .OR. &
                mmc_cor_energy (ic, MC_DISP)       .OR. &
                mmc_cor_energy (ic, MC_SPRING)     .OR. &
                mmc_cor_energy (ic, MC_LENNARD)    .OR. &
                mmc_cor_energy (ic, MC_REPULSIVE)  .OR. &
                mmc_cor_energy (ic, MC_COORDNUM )  .OR. &
                mmc_cor_energy (ic, MC_BUCKING)         ) THEN    
!                                                                       
!     ---------- Loop over all neighbours                               
!                                                                       
loop_neig:  DO j = 1, natom (icent) 
               js = cr_iscat (iatom (j, icent) ) 
               DO k = 1, 3 
                  u (k) = patom (k, 0, icent) 
               ENDDO 
!                                                                       
!     --------- Accumulate values for all Energies                      
!                                                                       
               IF (mmc_cor_energy (ic, MC_OCC)  ) THEN
!                                                                       
!     ----------- Chemical correlation, add number of atom pairs        
!                                                                       
                  pneig (is, js) = pneig (is, js) + 1 
               ENDIF 
               IF (mmc_cor_energy (ic, MC_UNI) ) THEN 
!                                                                       
!     ----------- Chemical correlation, add number of atom pairs        
!                                                                       
                  IF(mmc_pair(ic, MC_UNI,is,js) /=0) THEN
                  pneig (is, js) = pneig (is, js) + 1 
                  ENDIF
               ENDIF 
!
!
!--- Coordination number
!
               IF(mmc_cor_energy(ic, MC_COORDNUM)) THEN
                  IF(mmc_pair(ic,MC_COORDNUM, is,js)==-1) THEN
                     p_cn (is, js) = p_cn (is, js) + 1 
                  ENDIF
               ENDIF 
               IF (mmc_cor_energy (ic, MC_SPRING)    .OR.     &
                   mmc_cor_energy (ic, MC_LENNARD)   .OR.     &
                   mmc_cor_energy (ic, MC_REPULSIVE) .OR.     &
                   mmc_cor_energy (ic, MC_BUCKING)     ) THEN      
                  DO k = 1, 3 
                     v (k) = patom (k, j, icent) 
                     d (k) = v (k) - u (k) 
                  ENDDO 
!if(ic==3) THEN
!write(*,*) ' should add for atom pair is,js', i,is,js
!endif
                  dist = do_blen (.TRUE., u, v) 
                  js   = cr_iscat (iatom (j, icent) ) 
                  bl_sum (is, js) = bl_sum (is, js) + dist 
                  bl_s2  (is, js) = bl_s2 (is, js) + dist**2 
                  bl_anz (is, js) = bl_anz (is, js) + 1 
                  pneig (is, js) = pneig (is, js) + 1 
               ENDIF 
               IF (mmc_cor_energy (ic, MC_DISP) ) THEN 
                  CALL indextocell (iatom (j, icent), jcc, js) 
                  DO k = 1, 3 
                     disj (k) = cr_pos (k, iatom (j, icent) ) - chem_ave_pos (&
                     k, js) - REAL(jcc (k) - 1) - cr_dim0 (k, 1)            
                  ENDDO 
                  dpj = skalpro (disj, jdir, cr_gten) / rdj 
                  xij (is, js) = xij (is, js) + dpi * dpj 
                  xi2 (is, js) = xi2 (is, js) + dpi**2 
                  xj2 (is, js) = xj2 (is, js) + dpj**2 
                  xnn (is, js) = xnn (is, js) + 1 
               ENDIF 
            ENDDO loop_neig
!                     j ! Loop over all neighbours                      
         ENDIF is_energy
         IF (mmc_cor_energy (ic, MC_ANGLE) ) THEN 
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
                  searching = .TRUE. 
                  DO WHILE (searching.AND.k <= mmc_n_angles) 
                     k = k + 1 
                     CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls,  &
                                          MAXSCAT)
                     IF (iic == ic) THEN 
                        IF (iis == is.OR.iis ==  - 1) THEN 
                           IF (jjs ==  - 1.OR.jjs == min (js, ls) ) THEN 
                              IF (lls ==  - 1.OR.lls == max (js, ls) ) THEN 
                                 searching = .FALSE. 
                                 ic_a = k 
                              ENDIF 
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDDO 
!                         ! Find proper entry in correlation table      
                  IF (.NOT.searching) THEN 
                     DO k = 1, 3 
                        d (k) = patom (k, l, icent) 
                     ENDDO 
                     wi = do_bang (.TRUE., v, u, d) 
                     wis = mmc_target_angl (ic_a) 
                     IF (wis <= 90.) THEN 
                        IF (wi > 1.5 * wis) THEN 
                           wi = mod (wi + wis / 2., wis) + wis / 2. 
                        ENDIF 
                  ENDIF 
                  ba_sum (ic_a) = ba_sum (ic_a) + wi 
                  ba_s2 (ic_a) = ba_s2 (ic_a) + wi**2 
                  ba_anz (ic_a) = ba_anz (ic_a) + 1 
               ENDIF 
            ENDDO 
         ENDDO 
!                      ! j     ! Double loop over neighbours                   
      ENDIF 
      ENDIF is_cent    !       ! center atoms only                                   
   ENDDO main_cent     ! icent ! Loop over centers                               
ENDDO main_atoms       ! i     ! Loop over all atoms                                   
!                                                                       
!------ -- Summ up all energies, write output                           
!                                                                       
!     ----- Chemical correlation                                        
!                                                                       
CALL mmc_correlations_occ(ic, pneig, rel_cycl, damp, lout, MAXSCAT)
!                                                                       
!     ----- Unidirectional Chemical correlation                                        
!                                                                       
CALL mmc_correlations_uni(ic, pneig, rel_cycl, damp, lout, MAXSCAT)
!
!  Coordination number
!
   je = MC_COORDNUM
   n_cn = 0
   cn_pair: DO is = 0, cr_nscat 
      DO js = 0 , cr_nscat 
         IF(mmc_pair(ic, MC_COORDNUM, is, js) /=  0 ) THEN 
            n_cn = n_cn + p_cn(is, js)
         ENDIF
      ENDDO
   ENDDO cn_pair
!write(*,*) ' p_cn 1 ', p_cn(1,1:cr_nscat), n_cn, ncentral(1)
!write(*,*) ' p_cn 2 ', p_cn(2,1:cr_nscat), n_cn, ncentral(2)
!write(*,*) ' p_cn 3 ', p_cn(3,1:cr_nscat), n_cn, ncentral(3)
   cn_out: DO is = 0, cr_nscat 
      DO js = 0 , cr_nscat 
         IF(mmc_pair(ic, MC_COORDNUM, is, js) /=  0 ) THEN 
            mmc_ach_corr(ic, je, is, js) = REAL(n_cn        )/REAL(ncentral(is))
!           Feedback mechanism                                      
!           mmc_depth (ic, MC_COORDNUM, 0, 0) = mmc_depth (ic, MC_COORDNUM, 0, 0) - &
!           mmc_cfac (ic, MC_COORDNUM) * (mmc_target_corr (ic, MC_COORDNUM, is,js)- &
!                                            mmc_ach_corr (ic, MC_COORDNUM, is, js) ) / 2. &
!           *ABS(mmc_target_corr (ic, MC_OCC, is, js)) &
!           * damp
            mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
            mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            IF(lout) THEN
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE(output_io, 3150) ic, cr_at_lis(is), cr_at_lis(js),         &
                  mmc_target_corr(ic, je, is, js),                              &
                  mmc_ach_corr(ic, je, is, js),                                 &
                  mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js),&
                 (mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js))/divider,&
                  ncentral(is)
            ENDIF
         ENDIF
      ENDDO
   ENDDO cn_out
!                                                                       
!     ----- correlation of displacements                                
!                                                                       
disp_pair: DO is = 0, cr_nscat 
         DO js = is, cr_nscat 
         IF (mmc_pair (ic, MC_DISP, is, js) == -1 ) THEN 
         je = MC_DISP 
         xnn (is, js) = xnn (is, js) + xnn (js, is) 
         xij (is, js) = xij (is, js) + xij (js, is) 
         xi2 (is, js) = xi2 (is, js) + xi2 (js, is) 
         xj2 (is, js) = xj2 (is, js) + xj2 (js, is) 
         IF (xnn (is, js)  /= 0) THEN 
            xij (is, js) = xij (is, js) / REAL(xnn (is, js) ) 
            xi2 (is, js) = xi2 (is, js) / REAL(xnn (is, js) ) 
            xj2 (is, js) = xj2 (is, js) / REAL(xnn (is, js) ) 
!                                                                       
            IF (xi2 (is, js)  /= 0.AND.xj2 (is, js)  /= 0.0) THEN 
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
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
!                                                                       
         IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
            WRITE (output_io, 3200) ic, cr_at_lis (is), cr_at_lis (js), &
               mmc_target_corr (ic, je, is, js), mmc_ach_corr(ic, je, is,js), &
               mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic, je, is, js), &
               (mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js))/divider,&
               nneigh                                         
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
         IF (mmc_pair (ic, MC_SPRING, is, js) == -1 ) THEN 
         je = MC_SPRING 
!                                                                       
!     ----- Spring                                                      
!                                                                       
         IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum ( &
            js, is) ) / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2 (is, js) + bl_s2 (js,&
            is) ) / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is,  &
            js) - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
               mmc_ach_sigm (ic, je, is, js) = sqrt (mmc_ach_sigm (ic,  &
               je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
!                                                                       
            IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3300) ic, cr_at_lis (is),  cr_at_lis (js),   &
                  mmc_target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  &
                  mmc_ach_sigm (ic, je, is, js),                                     &
                  mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js), &
                 (mmc_target_corr(ic, je, is, js) - mmc_ach_corr(ic,je, is, js))/divider,&
                  bl_anz (is, js)  + bl_anz (js, is) 
            ENDIF 
         ENDIF 
         ENDIF 
         ENDDO 
         ENDDO  spri_pair
!                                                                       
!     -- Loop over all defined angle correlations                       
!                                                                       
angl_pair: IF (mmc_cor_energy (ic, MC_ANGLE) ) THEN 
         je = MC_ANGLE 
         DO k = 1, mmc_n_angles 
         CALL index2angles (mmc_angles (k), iic, kk, iis, jjs, lls,     &
         MAXSCAT)
         IF (ba_anz (k)  /= 0) THEN 
            mmc_ach_angl (k) = (ba_sum (k) ) / (ba_anz (k) ) 
            mmc_ang_sigm (k) = (ba_s2 (k) ) / (ba_anz (k) ) 
            mmc_ang_sigm (k) = (mmc_ang_sigm (k) - (mmc_ach_angl (k) ** &
            2) )                                                        
            IF (mmc_ang_sigm (k)  > 0) THEN 
               mmc_ang_sigm (k) = sqrt (mmc_ang_sigm (k) ) 
            ELSE 
               mmc_ang_sigm (k) = 0.0 
            ENDIF 
            mmc_ach_corr (ic, je, jjs, lls) = mmc_ach_angl (k) 
            mmc_ach_sigm (ic, je, jjs, lls) = mmc_ang_sigm (k) 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_angl (k               ) - mmc_ach_angl(k               )
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_angl(k               )
!                                                                       
            IF (lout) THEN 
               IF (iic == ic) THEN 
               IF(mmc_target_angl(k               )/=0.0) THEN
                  divider = mmc_target_angl(k)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3400) ic, cr_at_lis (iis),  cr_at_lis (jjs), &
                  cr_at_lis (lls),  mmc_target_angl (k),  mmc_ach_angl (k),   &
                  mmc_ang_sigm (k),                                           &
                  mmc_target_angl (k)  - mmc_ach_angl (k),                    &
                 (mmc_target_angl (k)  - mmc_ach_angl (k))/divider,           &
                  ba_anz (k)  + ba_anz (k)                                                        
                  IF (iis ==  - 1) THEN 
                     IF (jjs ==  - 1) THEN 
                        IF (lls ==  - 1) THEN 
                           searching = .FALSE. 
                           ic_a = k 
                        ELSE 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSE 
            IF (lout) THEN 
               IF (iic == ic) THEN 
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
      IF (mmc_pair (ic, MC_LENNARD, is, js) == -1 ) THEN 
         je = MC_LENNARD 
!                                                                       
!     ----- Lennard                                                     
!                                                                       
         IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum ( &
            js, is) ) / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2 (is, js) + bl_s2 (js,&
            is) ) / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is,  &
            js) - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
               mmc_ach_sigm (ic, je, is, js) = sqrt (mmc_ach_sigm (ic,  &
               je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3700) ic, cr_at_lis (is),  cr_at_lis (js), &
                   mmc_target_corr(ic, je, is, js),  mmc_ach_corr(ic, je, is, js), &
                   mmc_ach_sigm(ic, je, is, js),                                   &
                   mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic, je, is, js),&
                  (mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic, je, is, js))/divider,&
                   bl_anz (is, js)  + bl_anz (js, is) 
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
      IF (mmc_pair (ic, MC_REPULSIVE, is, js) == -1 ) THEN 
!write(*,*) is,js,ic, mmc_pair (ic, MC_REPULSIVE, is, js), bl_anz (is, js), bl_anz (js, is)
         je = MC_REPULSIVE 
!                                                                       
!     ----- REPULSIVE                                                     
!                                                                       
         IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
            mmc_ach_corr (ic, je, is, js) = (bl_sum (is, js) + bl_sum (js, is) ) &
                                          / (bl_anz (is, js) + bl_anz (js, is) )            
            mmc_ach_sigm (ic, je, is, js) = (bl_s2  (is, js) + bl_s2  (js, is) ) &
                                          / (bl_anz (is, js) + bl_anz (js, is) )                
            mmc_ach_sigm (ic, je, is, js) = (mmc_ach_sigm (ic, je, is, js)       &
                                          - (mmc_ach_corr (ic, je, is, js) **2) )                 
            IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
               mmc_ach_sigm (ic, je, is, js) = sqrt(mmc_ach_sigm(ic, je, is, js) )                                            
            ELSE 
               mmc_ach_sigm (ic, je, is, js) = 0.0 
            ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js)  - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
            IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
               WRITE (output_io, 3900) ic, cr_at_lis (is),  cr_at_lis (js),       &
               mmc_target_corr (ic, je, is, js),  mmc_ach_corr (ic, je, is, js),  &
               mmc_ach_sigm (ic, je, is, js),                                     &
               mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js), &
              (mmc_target_corr (ic, je, is, js)  - mmc_ach_corr (ic, je, is, js))/divider, &
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
         IF (mmc_pair (ic, MC_BUCKING, is, js) == -1 ) THEN 
            je = MC_BUCKING 
!                                                                       
!     ----- Buckingham                                                  
!                                                                       
            IF (bl_anz (is, js)  /= 0.OR.bl_anz (js, is)  /= 0) THEN 
               mmc_ach_corr (ic, je, is, js) =              &
                     (bl_sum (is, js) + bl_sum (js, is) ) / &
                     (bl_anz (is, js) + bl_anz (js, is) )            
               mmc_ach_sigm (ic, je, is, js) =              &
                     (bl_s2  (is, js) + bl_s2  (js, is) ) / &
                     (bl_anz (is, js) + bl_anz (js, is) )                
               mmc_ach_sigm (ic, je, is, js) =              &
                     (mmc_ach_sigm (ic, je, is, js) -       &
                     (mmc_ach_corr (ic, je, is, js) **2) )                 
               IF (mmc_ach_sigm (ic, je, is, js)  > 0) THEN 
                  mmc_ach_sigm (ic, je, is, js) =           &
                      sqrt (mmc_ach_sigm (ic, je, is, js) )                                            
               ELSE 
                  mmc_ach_sigm (ic, je, is, js) = 0.0 
               ENDIF 
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr (ic, je, is, js) - mmc_ach_corr(ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
               IF (lout) THEN 
               IF(mmc_target_corr(ic, je, is, js)/=0.0) THEN
                  divider = mmc_target_corr(ic, je, is, js)
               ELSE
                  divider = 1.0
               ENDIF
                  WRITE (output_io, 2100) ic, cr_at_lis(is), cr_at_lis(js),  &
                     mmc_ach_corr(ic, je, is, js), mmc_ach_sigm(ic, je,is, js), &
                     bl_anz (is, js) + bl_anz (js, is)               
               ENDIF 
            ENDIF 
         ENDIF 
      ENDDO 
   ENDDO buck_pair
ENDDO main_corr
!
DEALLOCATE(ncentral)
DEALLOCATE(patom)
DEALLOCATE(iatom)
DEALLOCATE(natom)
DEALLOCATE(xnn)
DEALLOCATE(xij)
DEALLOCATE(xi2)
DEALLOCATE(xj2)
DEALLOCATE(ba_sum)
DEALLOCATE(ba_s2)
DEALLOCATE(ba_anz)
DEALLOCATE( bl_anz)
DEALLOCATE( bl_sum)
DEALLOCATE( bl_s2)
DEALLOCATE( p_cn)
DEALLOCATE( pneig)
!
!  Check for convergence
!
done = .FALSE.
IF(mmc_h_stop) THEN     ! Apply convergence criteria
   IF(.NOT.lfinished .AND.mmc_h_ctarg>0) THEN
   DO i=1, mmc_h_number
      j = MOD(mmc_h_index+2,MMC_H_NNNN)     ! Increment Feedback number
      mmc_h_aver(i) = mmc_h_diff(i,0)       ! initialize average over MMC_H_NNNN feedbacks
      DO is = 1, MMC_H_NNNN-1
         mmc_h_aver(i) = mmc_h_aver(i) + mmc_h_diff(i,is)
      ENDDO
      mmc_h_aver(i) = mmc_h_aver(i)/MMC_H_NNNN
      mmc_h_maxd(i,mmc_h_index) = (mmc_h_diff(i,mmc_h_index) - mmc_h_diff(i,j))
      IF(mmc_h_targ(i) /= 0.0) THEN         ! Normalize relative to target value
         mmc_h_aver(i) = ABS(mmc_h_aver(i)/mmc_h_targ(i))
         mmc_h_maxd(i,mmc_h_index) = (mmc_h_maxd(i,mmc_h_index)/mmc_h_targ(i))
      ENDIF
   ENDDO
!mmc_h_ctarg
   ENDIF
   done = MAXVAL(    mmc_h_aver(1:mmc_h_ctarg))    < mmc_h_conv_m .AND.     &
          MAXVAL(ABS(mmc_h_maxd(1:mmc_h_ctarg,:))) < mmc_h_conv_c .AND.     &
          ABS(SUM(mmc_h_maxd(1:mmc_h_ctarg,0:MMC_H_NNNN-1)))/               &
             (REAL(mmc_h_ctarg*MMC_H_NNNN))   < mmc_h_conv_a
!if(mmc_h_ctarg>0) then
!WRITE(output_io,'(a,10F8.4)') ' Convergence criterium ', mmc_h_aver(1:mmc_h_ctarg), mmc_h_conv_m
!DO is = 0, MMC_H_NNNN-1
!WRITE(output_io,'(a,10F8.4)') ' Stagnant    changes   ', (mmc_h_maxd(1:mmc_h_ctarg,is))
!enddo
!WRITE(output_io,'(a,10F8.4)') ' Stagnant    largest   ', MAXVAL(ABS(mmc_h_maxd(1:mmc_h_ctarg,:))), &
! mmc_h_conv_c
!WRITE(output_io,'(a,10F8.4)') ' Stagnant    average   ', &
!SUM(mmc_h_maxd(1:mmc_h_ctarg,0:MMC_H_NNNN-1))/(REAL(mmc_h_ctarg*MMC_H_NNNN)), mmc_h_conv_a
!write(output_io,* ) ' DONE, finished ', done,lfinished, rel_cycl, mmc_h_ctarg
!endif
   IF(lout .AND. lfinished) THEN
      WRITE(output_io,*)
      IF(done) THEN
         WRITE(output_io,'(a,f6.0,a,i3,a)') ' Convergence reached after ',  &
            rel_cycl*100., ' % of user cycles', &
            MMC_H_NNNN, ' Feedbacks averaged'
      ELSEIF(rel_cycl==1.0) THEN
         WRITE(output_io,*) 'Maximum Cycles reached '
      ENDIF
      IF(done .OR. rel_cycl==1.0) THEN
      WRITE(output_io,'(a,f8.4,a, f8.4)') ' Largest relative difference ',    &
            MAXVAL(mmc_h_aver(1:mmc_h_number)), ' < ',                  &
            mmc_h_conv_m
      WRITE(output_io,'(a,f8.4,a, f8.4)') ' Stagnant: largest changes   ',    &
            MAXVAL(ABS(mmc_h_maxd(1:mmc_h_number,:))), ' < ',           &
            mmc_h_conv_c
      WRITE(output_io,'(a,f8.4,a, f8.4)') ' Stagnant: average changes   ',    &
            ABS(SUM(mmc_h_maxd(1:mmc_h_ctarg,0:MMC_H_NNNN-1)))/         &
            (REAL(mmc_h_ctarg*MMC_H_NNNN)), ' < ', mmc_h_conv_a
      WRITE(output_io,*)
      ENDIF
   ENDIF
ENDIF
!                                                                       
  410 FORMAT ( 45x,'Correlations/',/                                    &
     &   ' Neig.- Energy-',7x,'Atoms',11x,'Target',2x,'Distance/',4x,   &
     &    'Sigma',5x,'Diff',6x,'Diff/',5x,'Number',/                               &
     &    ' Def.   Type    central  Neighbors',13x,'Angle'              &
     &   ,27x,'Target  of pairs')                                               
 2100 FORMAT (1x,i3,3x,a9,3x,a9,5x,f7.3,3x,f7.3,3x,i8) 
 3100 FORMAT (1x,i3,3x,'Occupancy',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,i8)
 3150 FORMAT (1x,i3,3x,'Coord.No.',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
 3200 FORMAT (1x,i3,3x,'Disp.Cor.',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)                                           
 3300 FORMAT (1x,i3,3x,'Hooke    ',a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
 3400 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,5(f7.3,3x)         &
     &             ,i8)                                                 
 3410 FORMAT (1x,i3,3x,'Angle    ',a5,3x,a5,2x,a5,1x,  f7.3,23x         &
     &             ,f7.3,3x,i8)                                                 
 3700 FORMAT (1x,i3,3x,'Lennard  ',a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
 3900 FORMAT (1x,i3,3x,'Repuls.  ',a5,3x,a5,      8x,5(f7.3,3x)         &
     &             ,i8)                                                 
!                                                                       
END SUBROUTINE mmc_correlations               
!
!*******************************************************************************
!
SUBROUTINE mmc_correlations_occ(ic, pneig, rel_cycl, damp, lout, MAXSCAT_L)
!
!     ----- Chemical correlation                                        
!
USE crystal_mod 
USE mc_mod 
USE mmc_mod 
!
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: MAXSCAT_L
INTEGER, DIMENSION(0:MAXSCAT_L, 0:MAXSCAT_L) , INTENT(IN) :: pneig
REAL   , INTENT(IN) :: rel_cycl ! Relative progress along cycles
REAL   , INTENT(IN) :: damp
LOGICAL, INTENT(IN) :: lout
!
INTEGER :: pair11
INTEGER :: pair12
INTEGER :: pair21
INTEGER :: pair22
INTEGER :: is, js, je
INTEGER :: nneigh
LOGICAL :: lfirst
      REAL :: prob11=0.0, prob12, prob22 
REAL(PREC_SP) :: thet
REAL(PREC_SP) :: divisor
!
pair11 = 0
pair12 = 0
pair21 = 0
pair22 = 0
thet   = 0.0
DO is = 0, cr_nscat 
   DO js =  0, cr_nscat 
      IF     (mmc_pair (ic, MC_OCC, is, js) == -1 ) THEN 
         pair12 = pair12 + pneig (is,js)
      ELSEIF (mmc_pair (ic, MC_OCC, is, js) == -2 ) THEN 
         pair21 = pair21 + pneig (is,js)
      ELSEIF (mmc_pair (ic, MC_OCC, is, js) == +1 ) THEN 
         pair11 = pair11 + pneig (is,js)
      ELSEIF (mmc_pair (ic, MC_OCC, is, js) == +2 ) THEN 
         pair22 = pair22 + pneig (is,js)
      ENDIF
   ENDDO
ENDDO
je = MC_OCC 
!                                                                       
nneigh = pair11 + pair12 + pair21 + pair22 
IF (nneigh > 0.) THEN 
   prob11 =  pair11           / REAL(nneigh) 
   prob12 = (pair12 + pair21) / REAL(nneigh) 
   prob22 =  pair22           / REAL(nneigh) 
   thet = 0.5 * (2.0 * pair11 + pair12 + pair21) / REAL(nneigh)                                                     
!           thet = 0.5 * ((pair22 + pair11) + pair12 + pair21) / REAL(nneigh)                                                     
ENDIF 
lfirst = .TRUE.
corr_pair: DO is = 0, cr_nscat 
   DO js = is, cr_nscat 
      IF     (mmc_pair (ic, MC_OCC, is, js) /=  0 ) THEN 
        IF (thet /= 0.0.AND.thet /= 1.0) THEN 
           mmc_ach_corr (ic, je, is, js) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
           mmc_ach_corr (ic, je, js, is) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
        ELSE 
           IF((pair11>0 .OR. pair22>0) .AND. (pair12==0 .AND. pair21==0)) THEN
              mmc_ach_corr (ic, je, is, js) = 1.0 
              mmc_ach_corr (ic, je, js, is) = 1.0 
           ELSEIF((pair12>0 .OR. pair12>0) .AND. (pair11==0 .AND. pair22==0)) THEN
              mmc_ach_corr (ic, je, is, js) =-1.0 
              mmc_ach_corr (ic, je, js, is) =-1.0 
           ELSE 
              mmc_ach_corr (ic, je, is, js) = 0.0 
              mmc_ach_corr (ic, je, js, is) = 0.0 
           ENDIF 
        ENDIF 
!               Feedback mechanism                                      
        IF(mmc_target_corr (ic, MC_OCC, is, js) /= 0.0) THEN
           divisor = ABS(mmc_target_corr (ic, MC_OCC, is, js))
        ELSE
           divisor = 1.0
        ENDIF
        IF(rel_cycl>0.0) THEN
           mmc_depth (ic, MC_OCC, is, js) = mmc_depth      (ic, MC_OCC, is, js) -      &
                    mmc_cfac(ic, MC_OCC) * (mmc_target_corr(ic, MC_OCC, is, js) -      &
                                            mmc_ach_corr   (ic, MC_OCC, is, js) ) / 1. &
                   *ABS(mmc_target_corr(ic, MC_OCC, is, js)) * damp
           mmc_depth(ic, MC_OCC, js, is) = mmc_depth (ic, MC_OCC, is, js)
           mmc_depth(ic, MC_OCC,  0,  0) = mmc_depth (ic, MC_OCC, is, js)
!
!write(*,'(2(i3), 1x,6(f8.4,1x))') is, js,mmc_depth (ic, MC_OCC, is, js), & 
!mmc_depth (ic, MC_OCC, js, is), mmc_depth (ic, MC_OCC, 0,0), &
! (mmc_target_corr(ic, MC_OCC, is, js)-                               &
!  mmc_ach_corr   (ic, MC_OCC, is, js) ) / 2., damp,                  &
!-mmc_cfac (ic, MC_OCC) * (mmc_target_corr(ic, MC_OCC, is, js)- &
!                          mmc_ach_corr   (ic, MC_OCC, is, js) ) / 1. &
!        *ABS(mmc_target_corr (ic, MC_OCC, is, js)) &
!        * damp
!
        ENDIF
!                                                                       
        IF (lout .AND. mmc_pair(ic,MC_OCC,is,js) < 0 .AND. lfirst) THEN
           mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
           mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - &
                                                  mmc_ach_corr   (ic, je, is, js)
           mmc_h_targ(mmc_h_ctarg)              = mmc_target_corr(ic, je, is, js)
           lfirst = .FALSE.
           WRITE (output_io, 3100) ic, cr_at_lis (is), cr_at_lis (js),         &
               mmc_target_corr(ic, je, is, js),                                &
               mmc_ach_corr   (ic, je, is, js),                                &
               mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic,je, is, js), &
              (mmc_target_corr(ic, je, is, js) - mmc_ach_corr (ic,je, is, js))/divisor, &
               nneigh
!                                                                     
        ENDIF 
!                                                                       
      ENDIF 
   ENDDO 
ENDDO corr_pair
!
 3100 FORMAT (1x,i3,3x,'Occupancy',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
!
END SUBROUTINE mmc_correlations_occ
!
!*******************************************************************************
!
SUBROUTINE mmc_correlations_uni(ic, pneig, rel_cycl, damp, lout, MAXSCAT_L)
!
!     ----- Chemical correlation                                        
!
USE crystal_mod 
USE mc_mod 
USE mmc_mod 
!
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ic
INTEGER, INTENT(IN) :: MAXSCAT_L
INTEGER, DIMENSION(0:MAXSCAT_L, 0:MAXSCAT_L) , INTENT(IN) :: pneig
REAL   , INTENT(IN) :: rel_cycl ! Relative progress along cycles
REAL   , INTENT(IN) :: damp
LOGICAL, INTENT(IN) :: lout
!
INTEGER :: pair11
INTEGER :: pair12
INTEGER :: pair21
INTEGER :: pair22
INTEGER :: is, js, je
INTEGER :: nneigh
LOGICAL :: lfirst
      REAL :: prob11=0.0, prob12, prob22 
REAL(PREC_SP) :: thet
REAL(PREC_SP) :: divisor
!
pair11 = 0
pair12 = 0
pair21 = 0
pair22 = 0
thet   = 0.0
DO is = 0, cr_nscat 
   DO js =  0, cr_nscat 
      IF     (mmc_pair (ic, MC_UNI, is, js) == -1 ) THEN 
         pair12 = pair12 + pneig (is,js)
      ELSEIF (mmc_pair (ic, MC_UNI, is, js) == -2 ) THEN 
         pair21 = pair21 + pneig (is,js)
      ELSEIF (mmc_pair (ic, MC_UNI, is, js) == +1 ) THEN 
         pair11 = pair11 + pneig (is,js)
      ELSEIF (mmc_pair (ic, MC_UNI, is, js) == +2 ) THEN 
         pair22 = pair22 + pneig (is,js)
      ENDIF
   ENDDO
ENDDO
je = MC_UNI 
!                                                                       
nneigh = pair11 + pair12 + pair21 + pair22 
IF (nneigh > 0.) THEN 
   prob11 =  pair11           / REAL(nneigh) 
   prob12 = (pair12 + pair21) / REAL(nneigh) 
   prob22 =  pair22           / REAL(nneigh) 
   thet = 0.5 * (2.0 * pair11 + pair12 + pair21) / REAL(nneigh)                                                     
!           thet = 0.5 * ((pair22 + pair11) + pair12 + pair21) / REAL(nneigh)                                                     
ENDIF 
lfirst = .TRUE.
corr_pair: DO is = 0, cr_nscat 
   DO js = 0 , cr_nscat 
      IF     (mmc_pair (ic, MC_UNI, is, js) /=  0 ) THEN 
        IF (thet /= 0.0.AND.thet /= 1.0) THEN 
           mmc_ach_corr (ic, je, is, js) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
           mmc_ach_corr (ic, je, js, is) = (prob11 - thet**2) /&
                                           (thet * (1 - thet) )
        ELSE 
           IF((pair11>0 .OR. pair22>0) .AND. (pair12==0 .AND. pair21==0)) THEN
              mmc_ach_corr (ic, je, is, js) = 1.0 
              mmc_ach_corr (ic, je, js, is) = 1.0 
           ELSEIF((pair12>0 .OR. pair12>0) .AND. (pair11==0 .AND. pair22==0)) THEN
              mmc_ach_corr (ic, je, is, js) =-1.0 
              mmc_ach_corr (ic, je, js, is) =-1.0 
           ELSE 
              mmc_ach_corr (ic, je, is, js) = 0.0 
              mmc_ach_corr (ic, je, js, is) = 0.0 
           ENDIF 
        ENDIF 
!               Feedback mechanism                                      
        IF(mmc_target_corr (ic, MC_UNI, is, js) /= 0.0) THEN
           divisor = ABS(mmc_target_corr (ic, MC_UNI, is, js))
        ELSE
           divisor = 1.0
        ENDIF
        IF(rel_cycl>0.0) THEN
           mmc_depth (ic, MC_UNI, is, js) = mmc_depth      (ic, MC_UNI, is, js)    - &
                   mmc_cfac (ic, MC_UNI) * (mmc_target_corr(ic, MC_UNI, is, js)-     &
                                            mmc_ach_corr   (ic, MC_UNI, is, js) )/2. &
                  *ABS(mmc_target_corr (ic, MC_UNI, is, js)) * damp
           mmc_depth(ic, MC_UNI, js, is) = mmc_depth (ic, MC_UNI, is, js)
           mmc_depth(ic, MC_UNI,  0,  0) = mmc_depth (ic, MC_UNI, is, js)
        ENDIF
!                                                                       
        IF (lout .AND. mmc_pair(ic,MC_UNI,is,js) < 0 .AND. lfirst) THEN
         mmc_h_ctarg = mmc_h_ctarg + 1          ! Increment targets in history
         mmc_h_diff(mmc_h_ctarg, mmc_h_index) = mmc_target_corr(ic, je, is, js) - &
                                                mmc_ach_corr   (ic, je, is, js)
            mmc_h_targ(mmc_h_ctarg)           = mmc_target_corr(ic, je, is, js)
           lfirst = .FALSE.
           WRITE (output_io, 3100) ic, cr_at_lis (is), cr_at_lis (js),         &
               mmc_target_corr (ic, je, is, js),                               &
               mmc_ach_corr (ic, je, is, js),                                  &
               mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic,je, is, js),&
              (mmc_target_corr (ic, je, is, js) - mmc_ach_corr (ic,je, is, js))/divisor,&
               nneigh
        ENDIF 
!                                                                       
      ENDIF 
   ENDDO 
ENDDO corr_pair
!
 3100 FORMAT (1x,i3,3x,'Unidirect',a5,3x,a5,      8x,2(f7.3,3x),        &
     &        10x,f7.3,3x,f7.3,3x,i8)
!
END SUBROUTINE mmc_correlations_uni
!
!*****7*****************************************************************
!
SUBROUTINE check_geometry_multi (isel, natoms, laccept, MMC_MAX_ATOM_L) 
!-                                                                      
!     Check whether geometrical constrains apply to the movement of     
!     the selected atoms                                                
!+                                                                      
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
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, INTENT(IN) :: isel (MMC_MAX_ATOM_L) 
INTEGER, INTENT(IN) :: natoms 
LOGICAL, INTENT(OUT) :: laccept 
!                                                                       
INTEGER :: i, j 
LOGICAL,PARAMETER ::lspace = .TRUE.
REAL               :: d 
REAL, DIMENSION(3) :: u (3), null (3) 
!                                                                       
DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
laccept = .TRUE. 
!                                                                       
IF (mmc_l_constrains) THEN 
   IF(mmc_constrain_type == MMC_C_XYZ) THEN 
      loop_1: DO i = 1, natoms 
         DO j = 1, 3 
            IF(cr_pos(j, isel(i))  <  mmc_c_min(j) .OR.  &
               mmc_c_max(j) < cr_pos(j, isel(i))) THEN                         
               laccept = .FALSE. 
               EXIT loop_1
            ENDIF 
         ENDDO 
      ENDDO  loop_1
   ELSEIF(mmc_constrain_type == MMC_C_RADIUS) THEN 
      loop_2: DO i = 1, natoms 
         DO j = 1, 3 
            u(j) = cr_pos(j, isel(i)) 
         ENDDO 
         d = do_blen (lspace, u, null) 
         IF (mmc_c_rad <  d) THEN 
            laccept = .FALSE. 
            EXIT loop_2
         ENDIF 
      ENDDO  loop_2
   ENDIF 
ELSE 
   laccept = .TRUE. 
ENDIF 
!                                                                       
 9000 CONTINUE 
!                                                                       
END SUBROUTINE check_geometry_multi           
!
!*****7*****************************************************************
!
SUBROUTINE mmc_limit_selection(isel, natoms, MMC_MAX_ATOM_L) 
!-                                                                      
!     Selects an atom from a limited subset                             
!+                                                                      
USE crystal_mod 
USE celltoindex_mod
USE mmc_mod 
!USE errlist_mod 
USE lib_random_func
USE random_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER                           , INTENT(IN)  :: MMC_MAX_ATOM_L
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(OUT) :: isel 
INTEGER                           , INTENT(IN)  :: natoms 
!                                                                       
INTEGER :: i, j 
INTEGER, DIMENSION(3) :: icell  !(3)
INTEGER :: isite, nsel 
INTEGER :: ntrial 
REAL :: r1
!                                                                       
!                                                                       
IF(mmc_l_type == MMC_L_CELLS) THEN 
   DO i = 1, natoms 
      ntrial = mmc_l_extend(1)*mmc_l_extend(2)*mmc_l_extend(3)*cr_ncatoms
      CALL RANDOM_NUMBER(r1)
      nsel = INT(r1         * ntrial + 1 )
!     nsel = INT(ran1(idum) * ntrial + 1 )
      CALL mmc_indextocell(nsel, icell, isite, mmc_l_extend) 
      DO j = 1, 3 
         icell(j) = icell(j) + mmc_l_center(j) - 1 
      ENDDO 
      CALL celltoindex(icell, isite, nsel) 
      isel(i) = nsel 
   ENDDO 
ELSEIF(mmc_l_type == MMC_L_ATOMS) THEN 
   DO i = 1, natoms 
      CALL RANDOM_NUMBER(r1)
      isel(i) = INT(mmc_l_lower + r1         * (mmc_l_upper - mmc_l_lower + 1) )
!     isel(i) = INT(mmc_l_lower + ran1(idum) * (mmc_l_upper - mmc_l_lower + 1) )
   ENDDO 
ENDIF 
!                                                                       
END SUBROUTINE mmc_limit_selection            
!
!*****7*****************************************************************
!
SUBROUTINE mmc_indextocell (iatom, icell, isite, icc) 
!-                                                                      
!       calculates in which unit cell on which site the atom <ia> is    
!+                                                                      
USE crystal_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iatom
INTEGER, INTENT(OUT) :: isite
INTEGER, DIMENSION(3), INTENT(OUT) :: icell(3)
INTEGER, DIMENSION(3), INTENT(IN ) :: icc
INTEGER :: ia 
!                                                                       
ia = iatom - 1 
!                                                                       
icell(3) = INT(ia / icc(1) / icc(2) / cr_ncatoms) + 1 
ia       = ia - (icell(3) - 1) * icc(1) * icc(2) * cr_ncatoms 
icell(2) = INT(ia / icc (1) / cr_ncatoms) + 1 
ia       = ia - (icell(2) - 1) * icc(1) * cr_ncatoms 
icell(1) = INT(ia / cr_ncatoms) + 1 
isite    = ia - (icell(1) - 1) * cr_ncatoms + 1 
!                                                                       
END SUBROUTINE mmc_indextocell                
!
!*****7*****************************************************************
!
INTEGER FUNCTION angles2index (ic, nr, is, js, ls, MAXSCAT)
!-                                                                      
!     Calculates a unique number for the angle triplet                  
!+                                                                      
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ic 
INTEGER, INTENT(IN) :: nr 
INTEGER, INTENT(IN) :: is 
INTEGER, INTENT(IN) :: js 
INTEGER, INTENT(IN) :: ls 
INTEGER, INTENT(IN) :: MAXSCAT 
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
!
END FUNCTION angles2index                     
!
!*****7*****************************************************************
!
SUBROUTINE index2angles (ind, ic, nr, is, js, ls, MAXSCAT)
!-                                                                      
!     Calculates the correlation number, and the angle triplet from     
!     a unique number that was determined by angles2index               
!+                                                                      
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: X_ANGLES =  200  ! This needs work, may not be unique!
INTEGER, PARAMETER :: X_SCAT   =   50  ! This needs work, may not be unique!
!
INTEGER, INTENT(IN)  :: ind 
INTEGER, INTENT(OUT) :: ic 
INTEGER, INTENT(OUT) :: nr 
INTEGER, INTENT(OUT) :: is 
INTEGER, INTENT(OUT) :: js 
INTEGER, INTENT(OUT) :: ls 
INTEGER, INTENT(IN)  :: MAXSCAT 
!                                                                       
INTEGER i 
INTEGER na 
!                                                                       
i = ind 
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
!
!*****7**************************************************************** 
!
SUBROUTINE find_bucking (werte, MAXW) 
!-                                                                      
!     Finds the minimum of a Buckingham function                        
!+                                                                      
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAXW 
REAL(KIND=PREC_DP), INTENT(INOUT) :: werte (MAXW) 
!                                                                       
REAL(KIND=PREC_DP) :: last 
REAL(KIND=PREC_DP) :: next 
REAL(KIND=PREC_DP) :: step 
REAL(KIND=PREC_DP) :: old 
REAL(KIND=PREC_DP) :: new 
REAL(KIND=PREC_DP) :: minstep 
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
      DO WHILE (minstep < abs (step) ) 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      DO WHILE ( (new - old)  < 0) 
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
      werte (6) = buckingham (werte (2), werte (3), werte (4), werte (5))
!                                                                       
      last = werte (5) 
      old = werte (6) 
      step = werte (5) / 100. 
      next = last - step 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      DO WHILE ( (new - old)  > 0) 
      old = new 
      last = next 
      next = last - step 
      new = buckingham (werte (2), werte (3), werte (4), next) 
      ENDDO 
!                                                                       
      werte (7) = last 
      werte (8) = old 
!                                                                       
END SUBROUTINE find_bucking                   
!
!*****7**************************************************************** 
!
REAL FUNCTION buckingham(a, rho, b, x) 
!-                                                                      
!     calculates the Buckingham potential at point x                    
!+                                                                      
USE precision_mod
      IMPLICIT none 
!                                                                       
REAL(KIND=PREC_DP), INTENT(IN) :: a 
REAL(KIND=PREC_DP), INTENT(IN) :: rho 
REAL(KIND=PREC_DP), INTENT(IN) :: b 
REAL(KIND=PREC_DP), INTENT(IN) :: x 
!                                                                       
buckingham = a * exp ( - x / rho) - b / x**6 
!                                                                       
END FUNCTION buckingham                       
!
!*****7*************************************************************************
!
SUBROUTINE mmc_grow
!
! "grows" a correlation  
!
!  Initial quick and dirty method Special solution for PEROVSKITE
!
!
USE crystal_mod
USE atom_env_mod
USE chem_mod
USE do_find_mod
USE mc_mod
!
USE precision_mod
USE lib_random_func
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: MAXW = 6
!
INTEGER :: i, j, k                        ! Dummy loop indices
INTEGER :: ianz                           ! Number of atom types
INTEGER :: iatom
INTEGER :: istart                         ! Start index to minimize bias
INTEGER :: is_max                         ! Index of most neighbors
INTEGER :: nn_max                         ! Number of most neighbors
INTEGER, DIMENSION(6) :: nneig            ! Number of neighbors of type j
REAL(KIND=prec_SP), DIMENSION(3) :: x     ! Atom position
REAL(KIND=prec_DP), DIMENSION(3) :: werte ! Atom types allowed
REAL               :: rmin  ! minimum distance for find_env
REAL               :: rmax  ! maximum distance for find_env
REAL :: rel_cycl    ! how far are we in the desired number of cycles
LOGICAL :: lout_feed, done
REAL :: r1
!
done = .FALSE.
rmin = 0.0
rmax = cr_a0(1)*SQRT(3.0) + 0.1
werte(1) = -1.0D0
ianz     = 1
lout_feed = .TRUE.
DO i=1, mo_cyc
   rel_cycl = REAL(i)/REAL(mo_cyc)
   CALL RANDOM_NUMBER(r1)
   iatom = INT(r1     *cr_natoms) + 1    ! randomly choose an atom
!  iatom = INT(ran1(0)*cr_natoms) + 1    ! randomly choose an atom
   x = cr_pos(:, iatom)
   CALL do_find_env (ianz, werte, MAXW, x, rmin, rmax, chem_quick, chem_period)
   nneig = 0
   DO j=1,atom_env(0)                   ! Loop over all neighbors
      nneig(cr_iscat(atom_env(j))) = nneig(cr_iscat(atom_env(j))) + 1
   ENDDO
   is_max = 0
   nn_max = 0
   CALL RANDOM_NUMBER(r1)
   istart = INT(r1     *6.0)
!  istart = INT(ran1(0)*6.0)
   DO j=0,5
      k = MOD(istart+j,6) + 1           ! Current index in neighors
      IF(nneig(k) > nn_max) THEN
         nn_max = nneig(k)
         is_max = k
      ENDIF
   ENDDO
   cr_iscat(iatom) = is_max
   IF(MOD(INT(i,PREC_INT_LARGE), mo_feed)==0) THEN
      CALL mmc_correlations (lout_feed, rel_cycl, done, .FALSE.)
   ENDIF
ENDDO
!
END SUBROUTINE mmc_grow
!
!*****7*************************************************************************
!
END MODULE mmc_menu
