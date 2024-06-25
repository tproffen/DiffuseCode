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
use mmc_basic_mod
USE mc_mod 
USE mmc_mod 
USE mmc_mole
USE modify_mod
USE chem_symm_mod
use molecule_mod
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
CHARACTER(LEN=10)   :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=40)   :: cdummy 
CHARACTER(LEN=PREC_STRING) :: line, zeile 
INTEGER             :: lp, length 
INTEGER             :: indxg, lbef
!
INTEGER             :: n_corr = 1 ! dummy for allocation
INTEGER             :: n_scat = 1 ! dummy for allocation
INTEGER             :: n_site = 1 ! dummy for allocation
INTEGER             :: n_mole = 1 ! dummy for allocation
logical :: lout
logical :: lfeed
logical :: lfinished
logical :: done
real(kind=PREC_DP), dimension(2)                    :: maxdev =(/0.0, 0.0/)
real(kind=PREC_DP) :: rel_cycl
!                                                                       
!                                                                       
maxw = MAX(MIN_PARA,MAXSCAT+1)
!
! Basic allocation
!
n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
n_scat = MAX(MAXSCAT, MMC_MAX_SCAT, 3)    ! 3 is needed for 'group'
n_site = MAX(MAXSCAT, MMC_MAX_SITE)
n_mole = MOLE_MAX_TYPE
! call alloc_chem ! NEEDS WORK
CALL alloc_mmc ( n_corr, MC_N_ENERGY, n_scat, n_site )
CALL alloc_mmc_move(n_corr, n_scat, n_mole)
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
IF (indxg /= 0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
              .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ',   2, lbef, 4) )    &
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
         ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) THEN 
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
         ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) THEN 
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
         ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) THEN 
            CALL mmc_init 
!                                                                       
!------ command 'run'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) THEN 
            if(mmc_style .eqv. MMC_IS_ATOM) then
               CALL mmc_run_multi (.TRUE.)
            elseif(mmc_style .eqv. MMC_IS_MOLE) then
               call mmc_run_multi_mole (.TRUE.)
            endif
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
!SELECT  ELSEIF (str_comp (befehl, 'select', 3, lbef, 6) .OR.   &
!SELECT          str_comp (befehl, 'deselect', 2, lbef, 8) ) THEN
!SELECT                                                                 
!SELECT      CALL atom_select (zeile, lp, 0, MMC_MAX_SCAT, mmc_latom, &
!SELECT      mmc_lsite, 0, MMC_MAX_SITE,            &
!SELECT      mmc_sel_atom, lold, str_comp (befehl,  &
!SELECT      'sele', 3, lbef, 4) )                                       
!                                                                       
!------ command 'calc'                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'calc', 2, lbef, 4) ) THEN 
            if(zeile== ' ' .or.                    & ! No parameter
               str_comp(zeile, 'corr', 4, lp, 4)) then
               lout      = .false.
               rel_cycl  = 1.0D0
               done      = .true.
               lfinished = .true.
               CALL mmc_correlations (lout, rel_cycl, done, lfinished, lfeed, maxdev)
            else
               ier_num = -6
               ier_typ = -6
               ier_msg(1) = 'Parameter must be ''corr'' or absent '
            endif
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
      'inverse displacement', &
      'rotate molecule     ', &
      'switch neighbors    '  &
   /)
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
INTEGER :: n_mole = 1 ! dummy for allocation
!                                                                       
n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
n_scat = MAX(MAXSCAT, MMC_MAX_SCAT, 3)
n_site = MAX(MAXSCAT, MMC_MAX_SITE)
n_mole = MOLE_MAX_TYPE
! call alloc_chem ! NEEDS WORK
call alloc_mmc ( n_corr, MC_N_ENERGY, n_scat, n_site )
CALL alloc_mmc_move(n_corr, n_scat, n_mole)
!                                                                       
IF (mo_sel_atom) THEN 
   WRITE (output_io, 1105) 'atoms' 
ELSE 
   WRITE (output_io, 1105) 'molecules' 
ENDIF 
!                                                                       
WRITE (output_io, 1250) mo_cyc 
if(mmc_feed_auto) then
  write(output_io,'(a)') '   Feedback/update intervall adjusted automatically'
else
  WRITE (output_io, 1300) mo_feed 
endif
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
if( mo_sel_atom) then
   WRITE (output_io, 5000) 
   DO i = 1, MC_N_MOVE - 1
      WRITE (output_io, 5100) c_move (i), mmc_move_prob (i), &
                              c_site (mmc_local (i) )
   ENDDO 
else
  write(output_io, 5001)
   DO i = 1, MC_N_MOVE 
      WRITE (output_io, 5100) c_move (i), mmc_move_prob (i), &
                              c_site (mmc_local (i) )
   ENDDO 
endif
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
!DO i = 0, cr_nscat 
!!write(*,*) 'MMC_PAIR ',i,' : ', mmc_pair(k, MC_OCC, i,0:cr_nscat)
!write(fff,'(a,i1,a,i1,a)') '(a,i2,a,', cr_nscat+1, 'i4,', cr_nscat+1,  'f7.2)'
!write(*,fff                 ) 'MMC_PAIR ',i,' : ', mmc_pair(k, MC_OCC, i,0:cr_nscat),&
!                                                   mmc_depth(k,MC_OCC, i,0:cr_nscat)
!enddo
         DO ie = 1, 1
            DO i = 0, cr_nscat 
               DO j = i+1, cr_nscat 
                  at_name_i = at_name (i) 
                  at_name_j = at_name (j) 
                  IF (mmc_pair        (k, ie, i, j) <   0.0) THEN 
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
               IF (mmc_pair        (k, MC_OCC, i, j) <   0.0) THEN 
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
!DO i = 0, cr_nscat 
!write(fff,'(a,i1,a,i1,a)') '(a,i2,a,', cr_nscat+1, 'i4,', cr_nscat+1,  'f7.2)'
!write(*,fff                 ) 'MMC_PAIR ',i,' : ', mmc_pair(k, MC_UNI, i,0:cr_nscat),&
!                                                   mmc_depth(k,MC_UNI, i,0:cr_nscat)
!enddo
!read(*,*) i
         DO i = 0, cr_nscat 
            DO j = 0, cr_nscat 
               at_name_i = at_name (i) 
               at_name_j = at_name (j) 
!           IF (mmc_target_corr (k, MC_OCC, i, j)  /= 0.0) THEN 
               IF (mmc_pair        (k, MC_UNI, i, j) <   0.0) THEN 
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
               IF (mmc_pair        (k, MC_UNI, i, j) <   0.0) THEN 
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
               IF (mmc_depth (k, MC_COORDNUM, i, j)  /= 0.0) THEN 
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
!        IF (mo_sel_atom) THEN 
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
!        ELSE 
!           DO k = 1, chem_ncor 
!           DO i = 1, mole_num_type 
!           DO j = i, mole_num_type 
!           IF (mmc_target_corr (k, MC_LENNARD, i, j)  /= 0.0) THEN 
!              WRITE (output_io, 7700) at_name_i, at_name_j, k,         &
!              mmc_target_corr (k, MC_LENNARD, i, j),                   &
!              mmc_depth (k,MC_LENNARD, i, j),                          &
!              mmc_len_a (k, i, j), mmc_len_b (k, i,  j),               &
!              mmc_len_m (k, i, j), mmc_len_n (k, i, j)
!           ENDIF 
!           ENDDO 
!           ENDDO 
!           ENDDO 
!        ENDIF 
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
               IF(mo_maxmove_mole(4,i)==0.0) THEN
                  WRITE (output_io, 6020) i, (mo_maxmove_mole (j, i), j = 1, 3) 
               ELSE
                  WRITE (output_io, 6025) i, (mo_maxmove_mole (j, i), j = 1, 3) 
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
 2400 FORMAT (/,' Desired coordination numbers   : ',/,/,               &
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
 5000 FORMAT (/' Operation modes for MMC; atoms selected '/             &
     &        '   Mode                  Probability',                   &
     &        '  two-atom correlation'/3x,55('-'))                      
 5001 FORMAT (/' Operation modes for MMC; molecules selected '/         &
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
use mmc_basic_mod
USE modify_mod
use molecule_mod
USE rmc_sup_mod
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
USE string_convert_mod
use take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile 
INTEGER         , INTENT(INOUT) :: lp 
!                                                                       
INTEGER :: maxw = 200 
!
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(:), allocatable :: cpara !(maxw) 
REAL(KIND=PREC_DP)        , dimension(:), allocatable :: werte ! (maxw) 
integer                   , dimension(:), allocatable :: lpara! (maxw) 
INTEGER :: ianz  !, iianz, jjanz, kkanz, is, js, ls, ic, i, j 
integer :: i,j
INTEGER                :: n_corr ! Dummy for allocation
INTEGER                :: n_scat ! Dummy for allocation
!
INTEGER, PARAMETER :: NOPTIONAL = 4
INTEGER, PARAMETER :: O_TYPE    = 1
INTEGER, PARAMETER :: O_FEED    = 2
INTEGER, PARAMETER :: O_FINAL   = 3
INTEGER, PARAMETER :: O_SITE    = 4
CHARACTER(LEN=   5), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'type ', 'feed ', 'final', 'site'   /
DATA loname /  4     ,  4     ,   5    ,  4       /
opara  =  (/ 'atoms', 'on   ', 'on   ', '[0]  '   /)   ! Always provide fresh default values
lopara =  (/  5     ,  3     ,  3     ,  3  /)
owerte =  (/  0.0   ,  1.0   ,  1.0   ,  0.0  /)
!
MAXW = max(200,cr_nscat+10)
allocate(cpara(MAXW))
allocate(lpara(MAXW))
allocate(werte(MAXW))
!
!                                                                       
!     INTEGER angles2index 
!                                                                       
!                                                                       
n_corr = 0
n_scat = 0
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num == 0) THEN 
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   CALL do_cap(cpara(1))
   IF (ianz >= 2 .or. str_comp(cpara(1), 'OUTPUT', 2, lpara(1), 6)) THEN 
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
      ELSEIF (cpara(1)(1:3) == 'RBN') THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         if(cpara(1)=='classic') then
            mmc_algo = MMC_CLASSIC
         elseif(cpara(1)=='growth') then
            mmc_algo = MMC_GROWTH
         endif
         cpara(1) = '0'
         lpara(1) = 1
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         mmc_g_rate = werte(2)
         mmc_g_bad  = werte(3)
         mmc_g_neut = werte(4)
!
!------ --- 'set disallowed', set which atom are allowed in a move
!
      ELSEIF (cpara (1) (1:2)  == 'DI') THEN 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL get_iscat (ianz, cpara, lpara, werte,     &
                        maxw, .FALSE.)                                  
               IF (NINT(werte(1)) == -1) THEN
                  mmc_allowed = .false.
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
!     ELSEIF (cpara (1) (1:2)  == 'CY') THEN 
      elseif(str_comp(cpara(1), 'CYCLE', 2, lpara(1), 5)) then
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               mo_cyc = NINT(werte (1), PREC_INT_LARGE ) 
!                                                                       
!------ --- 'set feed' : setting display/feedback intervall             
!                                                                       
!     ELSEIF(cpara(1)(1:2)  == 'FE') THEN 
      elseif(str_comp(cpara(1), 'FEED', 2, lpara(1), 4)) then
         CALL del_params(1, ianz, cpara, lpara, maxw) 
         if(str_comp(cpara(1), 'auto', 3, lpara(1), 4)) then
            mmc_feed_auto = .TRUE.
         else
            CALL ber_params(ianz, cpara, lpara, werte, maxw) 
            mo_feed = nint (werte (1) ) 
         endif
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
               CALL mmc_set_mode (ianz, cpara, lpara, werte, maxw, opara(O_SITE), lopara(O_SITE)) 
!                                                                       
!------ --- 'set move': sets maxmove for shift MMC mode                 
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'MOV') THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (ier_num /= 0) return 
         IF(opara(O_TYPE) == 'atoms') THEN
            mmc_style   = MMC_IS_ATOM
            mo_sel_atom = .true.
            CALL rmc_set_move (mo_maxmove, mo_sel_atom, ianz, cpara, &
                 werte, lpara, maxw, 4, MAXSCAT)                                      
         ELSEIF(opara(O_TYPE) == 'mole') THEN
            mmc_style   = MMC_IS_MOLE
            mo_sel_atom = .false.
            CALL rmc_set_move(mo_maxmove_mole, .FALSE.,  &
                 ianz, cpara, werte, lpara, maxw, 4, MOLE_MAX_TYPE)
         ENDIF
!
!------ --- Output log
!
      elseif(str_comp(cpara(1), 'OUTPUT', 2, lpara(1), 6)) then
         if(str_comp(opara(O_FEED), 'on', 2, lopara(O_FEED), 2)) then
            mmc_out_feed = .true.
         elseif(str_comp(opara(O_FEED), 'off', 2, lopara(O_FEED), 3)) then
            mmc_out_feed = .false.
         else
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = 'Feedback parameter must be ''on'' or ''off'''
            return
         endif
!
         if(str_comp(opara(O_FINAL), 'on', 2, lopara(O_FINAL), 2)) then
            mmc_out_final = .true.
         elseif(str_comp(opara(O_FINAL), 'off', 2, lopara(O_FINAL), 3)) then
            mmc_out_final = .false.
         else
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = 'Feedback parameter must be ''on'' or ''off'''
            return
         endif
!                                                                       
!------ --- 'set rota': sets maxmove for molecule rotation MMC mode                 
!                                                                       
      ELSEIF (cpara (1) (1:3)  == 'ROT') THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF(opara(O_TYPE) == 'mole') THEN
            mmc_style   = MMC_IS_MOLE
            mo_sel_atom = .false.
            CALL rmc_set_move(mo_maxrota_mole, .FALSE.,  &
                 ianz, cpara, werte, lpara, maxw, 4, MOLE_MAX_TYPE)
         ENDIF
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
!------ --- 'set style' : setting of style atom/molecules               
!                                                                       
      ELSEIF(cpara(1)(1:2) == 'ST') THEN 
         if(cpara(2)(1:4) == 'atom') then
            mmc_style   = MMC_IS_ATOM
            mo_sel_atom = .true.
         elseif(cpara(2)(1:4) == 'mole') then
            mmc_style   = MMC_IS_MOLE
            mo_sel_atom = .false.
         endif
!                                                                       
!------ --- 'set target' : setting of target correlations               
!                                                                       
      ELSEIF (cpara (1) (1:2)  == 'TA') THEN 
         call mmc_set_target(ianz, cpara, lpara, MAXW)
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
               ELSEIF (str_comp (cpara (1) , 'RADIUS', 1, lpara (1) , 6) ) &
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
deallocate(cpara)
deallocate(lpara)
deallocate(werte)
!                                                                       
END SUBROUTINE mmc_set                        
!
!*******************************************************************************
!
subroutine mmc_set_target(ianz, cpara, lpara, MAXW)
!
! Handle 'set target' instruction
!+
use chem_mod
use crystal_mod
use discus_allocate_appl_mod
use mmc_mod
use mmc_basic_mod
use molecule_mod
!
use ber_params_mod
use get_params_mod
use get_iscat_mod
use errlist_mod
use precision_mod
use str_comp_mod
!
implicit none
!
integer                          , intent(inout) :: ianz
integer                          , intent(in   ) :: MAXW
character(len=*), dimension(MAXW), intent(inout) :: cpara
integer         , dimension(MAXW), intent(inout) :: lpara
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING), dimension(:), allocatable :: cpara1
character(len=PREC_STRING), dimension(:), allocatable :: cpara2
integer                   , dimension(:), allocatable :: lpara1
integer                   , dimension(:), allocatable :: lpara2
integer :: ianz1, ianz2, iianz, jjanz, kkanz
integer :: i, j         ! Loop index
integer :: is, js, ls   ! Atoms on site is, js, ls
integer :: ic           ! Correlation number
integer :: length       ! A string length
integer :: n_corr
integer :: n_scat
integer :: n_site
integer :: n_mole
integer :: n_angles
integer :: is_start, is_end
logical :: is_corr
real(kind=PREC_DP) :: a, b 
real(kind=PREC_DP) :: winit
real(kind=PREC_DP), dimension(MAXW)           :: werte
real(kind=PREC_DP), dimension(:), allocatable :: uerte
real(kind=PREC_DP), dimension(:), allocatable :: verte
real(kind=PREC_DP), dimension(:), allocatable :: werte1
real(kind=PREC_DP), dimension(:), allocatable :: werte2
!
if(ianz<3) then
   ier_num = -6
   ier_typ = ER_COMM
   return
endif
!        IF (ianz >= 3) THEN 
n_corr = 0
n_scat = 0
!
is = - 2 
js = - 2 
ls = - 2 
CALL del_params (1, ianz, cpara, lpara, maxw) 
CALL ber_params (1, cpara, lpara, werte, maxw) 
ic = nint (werte (1) ) 
!
IF( ic > CHEM_MAX_COR .OR. ic > MMC_MAX_CORR ) THEN
!
!  Basic allocation
!
   n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
   n_scat = MAX(MAXSCAT, MMC_MAX_SCAT, 3)
   n_site = MAX(MAXSCAT, MMC_MAX_SITE)
   n_mole = MOLE_MAX_TYPE
   ! call alloc_chem ! NEEDS WORK
   CALL alloc_mmc ( n_corr, MC_N_ENERGY, n_scat, n_site )
   CALL alloc_mmc_move(n_corr, n_scat, n_mole)
ENDIF
!
if(ic<=0 .or. ic>chem_ncor) then
   ier_num = - 14 
   ier_typ = ER_CHEM 
   return
endif
!
allocate(cpara1(MAXW))
allocate(cpara2(MAXW))
allocate(lpara1(MAXW))
allocate(lpara2(MAXW))
allocate(werte1(MAXW))
allocate(werte2(MAXW))
allocate(uerte (MAXW))
allocate(verte (MAXW))
werte  = 0.0D0
uerte  = 0.0D0
verte  = 0.0D0
werte1 = 0.0D0
werte2 = 0.0D0
!
!           IF (ic > 0.AND.ic <= chem_ncor) THEN 
cond_type: IF(str_comp(cpara(2), 'corr', 2, lpara(2), 4) .OR.  &
   str_comp(cpara(2), 'unid', 2, lpara(2), 4)) THEN
!
   is_corr= str_comp (cpara (2) , 'corr', 2, lpara (2) , 4)
   CALL del_params (2, ianz, cpara, lpara, maxw) 
!     Get atom typeDA(ianz, cpara, lpara, MAXW)s in the two allowed groups
!     Allowed parameters are: Au, Cu            ! OLd style 
!                             (Au), (Cu)        ! New style  :binary correlation
!                             (Au,Pt), (Cu,Zn)  ! New style  :quaternary correlation
   IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
      line   = cpara(1)(2:lpara(1)-1)
      length = lpara(1)-2
      CALL get_params (line, ianz1, cpara1, lpara1, MAXW, length) 
      CALL get_iscat (ianz1, cpara1, lpara1, werte1, MAXW, .FALSE.)                                  
   ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
      ianz1 = 1 
      CALL get_iscat(ianz1, cpara, lpara, werte1, MAXW, .FALSE.)                                  
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      exit cond_type
   ENDIF
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
      line   = cpara(1)(2:lpara(1)-1)
      length = lpara(1)-2
      CALL get_params (line, ianz2, cpara2, lpara2, MAXW, length) 
      CALL get_iscat (ianz2, cpara2, lpara2, werte2, MAXW, .FALSE.)                                  
   ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
      ianz2 = 1 
      CALL get_iscat (ianz2, cpara, lpara, werte2, MAXW, .FALSE.)                                  
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      exit cond_type
   ENDIF
!
   IF(is_corr                                    ) THEN                                             
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      IF (cpara (ianz) (1:2)  == 'CO') THEN 
         mmc_cfac (ic, MC_OCC) =  1.0 
         mmc_lfeed(ic, MC_OCC) =  .true.
         ianz = ianz - 1 
      ELSEIF (cpara (ianz) (1:2)  == 'EN') THEN 
         mmc_cfac (ic, MC_OCC) = 0.0 
         mmc_lfeed(ic, MC_OCC) =  .false.
         ianz = ianz - 1 
      ENDIF 
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF (mmc_cfac (ic, MC_OCC) > 0.0) THEN 
         if(abs(werte(1))<0.1) then
            winit = 0.0D0
         else
         winit = -2.0D0*werte(1) - 3.0D0*werte(1)**5
         winit = -1.0D0*werte(1)  - 3.0D0*werte(1)**5
         endif
         CALL mmc_set_disp_occ (ic, MC_OCC, ianz1, ianz2, &
                          MAXW, werte1, werte2, werte(1) , winit)
                        mmc_depth (ic, MC_OCC, 0, 0) = winit
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
         mmc_lfeed(ic, MC_UNI) =  .true.
         ianz = ianz - 1 
      ELSEIF (cpara (ianz) (1:2)  == 'EN') THEN 
         mmc_cfac (ic, MC_UNI) = 0.0 
         mmc_lfeed(ic, MC_UNI) =  .false.
         ianz = ianz - 1 
      ENDIF 
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF (mmc_cfac (ic, MC_UNI) > 0.0) THEN 
         winit = -2.0D0*werte(1) - 3.0D0*werte(1)**5
         CALL mmc_set_unid_occ (ic, MC_UNI, ianz1, ianz2, &
                       MAXW, werte1, werte2, werte(1) , winit)
                       mmc_depth (ic, MC_UNI, 0, 0) = winit
      ELSEIF (mmc_cfac (ic, MC_UNI) ==0.0) THEN 
         CALL mmc_set_unid_occ (ic, MC_UNI, ianz1, ianz2, &
                       MAXW, werte1, werte2, werte(1) , werte(2) )                                   
                       mmc_depth (ic, MC_UNI, 0, 0) = werte (2) 
      ENDIF
      mmc_cor_energy (ic, MC_UNI) = .TRUE. 
      mmc_cor_energy (0, MC_UNI) = .TRUE. 
   ENDIF
ELSEIF(str_comp(cpara(2), 'group', 2, lpara(2) , 5)) THEN  cond_type ! Group wise correlations
   CALL set_target_group(MAXW, ianz, cpara, lpara, werte, ic)
ELSEIF(str_comp(cpara(2), 'cd', 2, lpara(2), 2) ) THEN     cond_type 
   CALL del_params (2, ianz, cpara, lpara, maxw) 
   iianz = 1 
   jjanz = 1 
   CALL get_iscat(iianz, cpara, lpara, uerte, MAXW, .FALSE.)
   CALL del_params(1, ianz, cpara, lpara, maxw) 
   CALL get_iscat (jjanz, cpara, lpara, verte, MAXW, .FALSE.)
   CALL del_params(1, ianz, cpara, lpara, maxw) 
   IF (cpara (ianz) (1:2)  == 'CO') THEN 
      mmc_cfac (ic, MC_DISP) =  1.00 
      mmc_lfeed(ic, MC_DISP) =  .true.
      ianz = ianz - 1 
   ELSEIF (cpara (ianz) (1:2)  == 'EN') THEN 
      mmc_cfac (ic, MC_DISP) = 0.0 
      mmc_lfeed(ic, MC_DISP) =  .false.
      ianz = ianz - 1 
   ENDIF 
   CALL ber_params (2, cpara, lpara, werte, maxw) 
   IF (mmc_cfac (ic, MC_DISP) /= 0.0) THEN 
      werte(2) = -500.*werte(1)
   ENDIF
   DO i = 1, iianz 
      DO j = 1, jjanz 
         is = nint(uerte (i) ) 
         js = nint(verte (j) ) 
            IF(is<0 .OR. is>cr_nscat .OR. &
               js<0 .OR. js>cr_nscat)  THEN
               ier_num = -97
               ier_typ = ER_APPL
               WRITE(ier_msg(1),'(a,i10)') 'Atom type ',is
               WRITE(ier_msg(2),'(a,i10)') 'Atom type ',js
               ier_msg(3) = 'Check parameter list for wrong/missing atoms'
               exit cond_type
            ENDIF
            mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
            mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
         CALL mmc_set_disp (ic, MC_DISP, is, js, werte ( &
                        1), werte (2) )                                 
      ENDDO 
   ENDDO 
!
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
ELSEIF(str_comp(cpara(2) , 'cn', 2, lpara(2), 2) ) THEN  cond_type  !Coordination number
   CALL del_params (2, ianz, cpara, lpara, maxw) 
!  Get atom types in the two allowed groups
!  Allowed parameters are: Au, Cu            ! OLd style 
!                          (Au), (Cu)        ! New style  :binary correlation
!                          (Au,Pt), (Cu,Zn)  ! New style  :quaternary correlation
   IF ( cpara(1)(1:1)=='(' .AND. cpara(1)(lpara(1):lpara(1))==')') THEN
      line   = cpara(1)(2:lpara(1)-1)
      length = lpara(1)-2
      CALL get_params (line, ianz1, cpara1, lpara1, MAXW, length) 
      CALL get_iscat (ianz1, cpara1, lpara1, werte1, MAXW, .FALSE.)                                  
   ELSEIF ( cpara(1)(1:1)/='(' .AND. cpara(1)(lpara(1):lpara(1))/=')') THEN
      ianz1 = 1 
      CALL get_iscat (ianz1, cpara, lpara, werte1, MAXW, .FALSE.)                                  
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      exit cond_type
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
      exit cond_type
   ENDIF
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   if(ianz==1) werte(2) = 1.0_PREC_DP
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
ELSEIF(str_comp(cpara(2), 'spring', 2, lpara(2), 6) ) THEN cond_type
   CALL del_params (2, ianz, cpara, lpara, maxw) 
   iianz = 1 
   jjanz = 1 
   CALL get_iscat (iianz, cpara, lpara, uerte,     &
                        MAXW, .FALSE.)                                  
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL get_iscat (jjanz, cpara, lpara, verte,     &
                        MAXW, .FALSE.)                                  
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL ber_params (2, cpara, lpara, werte, maxw) 
!
   if(iianz==1 .and. nint(uerte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         uerte(i+1) = real(i,kind=PREC_DP)
      enddo
      iianz = cr_nscat+1
   endif
   if(jjanz==1 .and. nint(verte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         verte(i+1) = real(i,kind=PREC_DP)
      enddo
      jjanz = cr_nscat+1
   endif
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
            exit cond_type
         ENDIF
         mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
         mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
         CALL mmc_set_disp (ic, MC_SPRING, is, js, werte &
                        (1), werte (2) )                                
      ENDDO 
   ENDDO 
   mmc_cor_energy (ic, MC_SPRING) = .TRUE. 
   mmc_cor_energy (0, MC_SPRING) = .TRUE. 
ELSEIF(str_comp(cpara(2) , 'angle', 2, lpara(2), 5) ) THEN cond_type
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
      CALL get_iscat(iianz, cpara, lpara, uerte, maxw, .FALSE.)
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
      CALL get_iscat (jjanz, cpara, lpara, uerte, MAXW, .FALSE.)                               
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL get_iscat (kkanz, cpara, lpara, verte, MAXW, .FALSE.)                               
      js = min (uerte (1), verte (1) ) 
      ls = max (uerte (1), verte (1) ) 
      mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
      mmc_allowed(ls) = .TRUE. ! this atom is allowed in mmc moves
      i = angles2index (ic, mmc_n_angles, is, js, ls, MAXSCAT)                 
      mmc_angles (mmc_n_angles) = i 
      CALL index2angles (i, ic, mmc_n_angles, is,  js, ls, MAXSCAT)             
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL ber_params (2, cpara, lpara, werte, maxw)                                        
      mmc_target_angl (mmc_n_angles) = werte (1) 
      mmc_target_corr (ic, MC_ANGLE, js, ls) = werte (1)                                  
      mmc_target_corr (ic, MC_ANGLE, ls, js) = werte (1)                                  
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
ELSEIF(str_comp(cpara(2), 'lennard', 2, lpara(2), 6) ) THEN cond_type
   CALL del_params (2, ianz, cpara, lpara, maxw) 
   iianz = 1 
   jjanz = 1 
   CALL get_iscat(iianz, cpara, lpara, uerte, MAXW, .FALSE.)                                  
   CALL del_params(1, ianz, cpara, lpara, maxw) 
   CALL get_iscat(jjanz, cpara, lpara, verte, MAXW, .FALSE.)                                  
   CALL del_params(1, ianz, cpara, lpara, maxw) 
   CALL ber_params(ianz, cpara, lpara, werte, maxw)                                           
!
   if(iianz==1 .and. nint(uerte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         uerte(i+1) = real(i,kind=PREC_DP)
      enddo
      iianz = cr_nscat+1
   endif
   if(jjanz==1 .and. nint(verte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         verte(i+1) = real(i,kind=PREC_DP)
      enddo
      jjanz = cr_nscat+1
   endif
!
   IF (ianz == 3) THEN 
      werte(4) = 12.0D0
      werte(5) = 6.0D0
   ENDIF 
   IF(ic > MMC_LENN_CORR .OR.  & ! Allocate Lennard
      ic > CHEM_MAX_COR  .OR.  &
      MAXSCAT > MMC_LENN_SCAT   ) THEN
      n_corr = MAX(n_corr, CHEM_MAX_COR, &
                           MMC_LENN_CORR)
      n_scat = MAX(MAXSCAT,MMC_LENN_SCAT)
      call alloc_mmc_lenn (n_corr, n_scat)
      IF(ier_num /= 0) THEN
         exit cond_type
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
               exit cond_type
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
ELSEIF(str_comp(cpara(2), 'repulsive', 2, lpara(2), 9)        ) THEN cond_type
   CALL del_params (2, ianz, cpara, lpara, maxw) 
   iianz = 1 
   jjanz = 1 
   CALL get_iscat (iianz, cpara, lpara, uerte, MAXW, .FALSE.)                                  
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL get_iscat (jjanz, cpara, lpara, verte, MAXW, .FALSE.)                                  
   CALL del_params (1, ianz, cpara, lpara, maxw) 
!
   if(iianz==1 .and. nint(uerte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         uerte(i+1) = real(i,kind=PREC_DP)
      enddo
      iianz = cr_nscat+1
   endif
   if(jjanz==1 .and. nint(verte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         verte(i+1) = real(i,kind=PREC_DP)
      enddo
      jjanz = cr_nscat+1
   endif
!
   IF (ianz >  0) THEN 
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
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
         exit cond_type
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
            exit cond_type
         ENDIF
         mmc_allowed(is) = .TRUE. ! this atom is allowed in mmc moves
         mmc_allowed(js) = .TRUE. ! this atom is allowed in mmc moves
         CALL mmc_set_disp (ic, MC_REPULSIVE, is, js, 100.0D0, ABS(werte(1)) )                          
         CALL mmc_set_rep(ic, is, js, ABS(werte(1)), werte(2) , werte(3), werte(4))
      ENDDO 
   ENDDO 
   mmc_cor_energy (ic, MC_REPULSIVE) = .TRUE. 
   mmc_cor_energy (0, MC_REPULSIVE) = .TRUE. 
ELSEIF(str_comp(cpara(2), 'bucking', 2, lpara(2), 7) ) THEN cond_type
   CALL del_params (2, ianz, cpara, lpara, maxw) 
   iianz = 1 
   jjanz = 1 
   CALL get_iscat (iianz, cpara, lpara, uerte, MAXW, .FALSE.)                                  
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL get_iscat (jjanz, cpara, lpara, verte, MAXW, .FALSE.)                                  
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   CALL ber_params (ianz, cpara, lpara, werte, maxw)                                           
!
   if(iianz==1 .and. nint(uerte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         uerte(i+1) = real(i,kind=PREC_DP)
      enddo
      iianz = cr_nscat+1
   endif
   if(jjanz==1 .and. nint(verte(1))==-1) then   ! 'all' atom typs selected
      do i=0,cr_nscat
         verte(i+1) = real(i,kind=PREC_DP)
      enddo
      jjanz = cr_nscat+1
   endif
!
   IF(ic > MMC_BUCK_CORR .OR.  & ! Allocate Buckingham
      ic > CHEM_MAX_COR  .OR.  &
      MAXSCAT > MMC_BUCK_SCAT   ) THEN
      n_corr = MAX(n_corr, CHEM_MAX_COR, MMC_BUCK_CORR)
      n_scat = MAX(MAXSCAT,MMC_BUCK_SCAT)
      call alloc_mmc_buck (n_corr, n_scat)
      IF(ier_num /= 0) THEN
         exit cond_type
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
            exit cond_type
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
         CALL mmc_set_disp(ic, MC_BUCKING, is, js, werte(5), abs(werte(6) ) )
      ENDDO 
   ENDDO 
   mmc_cor_energy (ic, MC_BUCKING) = .TRUE. 
   mmc_cor_energy (0, MC_BUCKING) = .TRUE. 
ELSE  cond_type
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF  cond_type
!                                                                       
!                 ELSE 
!                    ier_num = - 14 
!                    ier_typ = ER_CHEM 
!                 ENDIF 
!              ELSE 
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!              ENDIF 
!
if(ier_num==0) mmc_h_number = mmc_h_number + 1
!
deallocate(cpara1)
deallocate(cpara2)
deallocate(lpara1)
deallocate(lpara2)
deallocate(werte1)
deallocate(werte2)
deallocate(uerte )
deallocate(verte )
!
end subroutine mmc_set_target
!
!*****7*****************************************************************
!
SUBROUTINE set_target_group(MAXW, ianz, cpara, lpara, werte, ic)
!-
!  Set the target values for the "group" correlations 
!  set target, group, 1, (a,b,c), (d,e,f), corr, ener, "CORR" | "ENER"
!
USE get_iscat_mod
USE mc_mod
USE mmc_mod
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)    :: MAXW
INTEGER, INTENT(INOUT) :: ianz
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(INOUT) :: werte
INTEGER, INTENT(IN)    :: ic
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara1
INTEGER                   , DIMENSION(MAXW) :: lpara1
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte1
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte2
!
INTEGER :: ianz1
INTEGER :: ianz2
!
!
cpara1 = ' '
lpara1 =  0
cpara1(1) = cpara(3)
lpara1(1) = lpara(3)
ianz1     = 1
CALL get_iscat(ianz1, cpara1, lpara1, werte1, MAXW, .FALSE.)
!
!
cpara1 = ' '
lpara1 =  0
cpara1(1) = cpara(4)
lpara1(1) = lpara(4)
ianz2     = 1
CALL get_iscat(ianz2, cpara1, lpara1, werte2, MAXW, .FALSE.)
!
!
CALL del_params (4, ianz, cpara, lpara, maxw) 
!
!
IF (cpara(ianz)(1:2)  == 'CO') THEN 
   mmc_cfac(ic, MC_GROUP) =  1.0 
   mmc_lfeed(ic, MC_GROUP) =  .true.
   ianz = ianz - 1 
ELSEIF(cpara(ianz)(1:2)  == 'EN') THEN 
   mmc_cfac(ic, MC_GROUP) = 0.0 
   mmc_lfeed(ic, MC_GROUP) =  .false.
   ianz = ianz - 1 
ENDIF 
!
CALL ber_params(ianz, cpara, lpara, werte, maxw)
!
IF(mmc_cfac(ic, MC_GROUP) > 0.0) THEN 
   CALL mmc_set_group_occ (ic, MC_GROUP, ianz1, ianz2, &
        MAXW, werte1, werte2, werte(1) , -0.50*werte(1) )                                   
ELSEIF(mmc_cfac(ic, MC_GROUP) ==0.0) THEN 
   CALL mmc_set_group_occ (ic, MC_GROUP, ianz1, ianz2, &
        MAXW, werte1, werte2, werte(1) , werte(2) )                                   
ENDIF
mmc_cor_energy(ic, MC_GROUP) = .TRUE. 
mmc_cor_energy(0,  MC_GROUP) = .TRUE. 
!
END SUBROUTINE set_target_group
!
!*****7*****************************************************************
!
SUBROUTINE get_finish(ianz, cpara, lpara, werte, MAXW)
!
USE mmc_mod
!
USE errlist_mod
USE precision_mod
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
INTEGER, PARAMETER :: NOPTIONAL = 7
INTEGER, PARAMETER :: O_FEED    = 1                  ! Number of feedbacks to average
INTEGER, PARAMETER :: O_DIFF    = 2                  ! Maximum difference
INTEGER, PARAMETER :: O_RDIFF   = 3                  ! Maximum difference
INTEGER, PARAMETER :: O_CHANGE  = 4                  ! Maximum change in differences
INTEGER, PARAMETER :: O_AVER    = 5                  ! Maximum change in differences
INTEGER, PARAMETER :: O_STOP    = 6                  ! How to stop cycles/convergence
INTEGER, PARAMETER :: O_LOG     = 7                  ! Screen log ?
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte    ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 5 ! Number of values to calculate
!
DATA oname  / 'feed  ', 'diff  ', 'rdiff ', 'change', 'aver  ', 'stop  ' , 'log   ' /   ! mmc_set capitalizes only the first parameter
DATA loname /  4      ,  4      ,  5      ,  6      ,  4      ,  4       ,  3       /
!
opara  = (/'3     ', '99999.', '99999.', '99999.', '99999.', 'cycles' ,'none  '/)
lopara = (/ 1      ,  6      ,  6      ,  6      ,  6      , 6        , 4      /)
owerte = (/ 3.     ,  99999. ,  99999. ,  99999. ,  99999. , 0.0      , 0.0    /)
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
   IF(opara(O_LOG) == 'screen') THEN
      mmc_h_log = .TRUE.
   ELSEIF(opara(O_LOG) == 'none') THEN
      mmc_h_log = .FALSE.
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'Optional log must be log:none or '
      ier_msg(2) = '                     log:screen '
      RETURN
   ENDIF
   MMC_H_NNNN   = NINT(owerte(O_FEED))    ! Number of feedbacks to average
   mmc_h_conv_m = owerte(O_DIFF)          ! Largest difference (target-achieved)
   mmc_h_conv_r = owerte(O_RDIFF)         ! Largest difference (target-achieved)/target
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
      mmc_depth (ic, ie, is, js) = depth * coord
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
         DO j=1,ianz2
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
!write(*,*) ' First group ', werte1(1:ianz1)
!write(*,*) ' 2nd   group ', werte2(1:ianz2)
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
   DO j=1,ianz2
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
      mmc_pair        (ic, ie, js, is) = -2     ! Second ==> first group
   END DO
END DO
!
      mmc_target_corr (ic, ie, 0 , 0 ) = corr 
!
END SUBROUTINE mmc_set_unid_occ
!
!*****7*****************************************************************
!
SUBROUTINE mmc_set_group_occ(ic, ie, ianz1, ianz2, &
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
mmc_left (ic, ie, :) = 0
mmc_right(ic, ie, :) = 0
DO i=1,ianz1              ! Set "equal" pairs first group
   is = NINT(werte1(i))
   mmc_allowed(is) = .TRUE.
   DO j=1,ianz1
      js = NINT(werte1(j))
      mmc_target_corr(ic, ie, is, js) = corr 
      mmc_depth      (ic, ie, is, js) = depth 
      mmc_pair       (ic, ie, is, js) = +1     ! These pairs contribute positively to energy
   END DO
   mmc_left(ic, ie, is    ) = +1     ! These pairs contribute positively to energy
END DO
!
DO i=1,ianz2              ! Set "equal" pairs second group
   is = NINT(werte2(i))
   mmc_allowed(is) = .TRUE.
   DO j=1,ianz2
      js = NINT(werte2(j))
      mmc_target_corr(ic, ie, is, js) = corr 
      mmc_depth      (ic, ie, is, js) = depth 
      mmc_pair       (ic, ie, is, js) = +2     ! These pairs contribute positively to energy
   END DO
   mmc_right(ic, ie, is    ) = +2     ! These pairs contribute positively to energy
END DO
!
DO i=1,ianz2              ! Set "opposite" pairs second group
   is = NINT(werte2(i))
   mmc_allowed(is) = .TRUE.
   DO j=1,ianz1
      js = NINT(werte1(j))
      mmc_target_corr(ic, ie, is, js) = corr 
      mmc_depth      (ic, ie, is, js) = depth 
      mmc_pair       (ic, ie, is, js) = -1     ! These pairs contribute positively to energy
   END DO
END DO
!
DO i=1,ianz1              ! Set "opposite" pairs first group
   is = NINT(werte1(i))
   mmc_allowed(is) = .TRUE.
   DO j=1,ianz2
      js = NINT(werte2(j))
      mmc_target_corr(ic, ie, is, js) = corr 
      mmc_depth      (ic, ie, is, js) = depth 
      mmc_pair       (ic, ie, is, js) = -1     ! These pairs contribute positively to energy
   END DO
END DO
!
!write(*,*)  ' W1 ', werte1(1:ianz1),' ENERGY ', ie
!write(*,*)  ' W2 ', werte2(1:ianz2),' ENERGY ', ie
!write(*,'(a, 12i4)') ' MMC_LEFT  ', mmc_left (ic, ie,:)
!write(*,'(a, 12i4)') ' MMC_RIGHT ', mmc_right(ic, ie,:)
!write(*,'(a, f8.2)') ' TARGET IS ', corr
!do is=0, cr_nscat
!write(*,'(a,12f7.2)') ' TARGET ',mmc_target_corr(ic, ie, is, :)
!enddo
!write(*,'(a, f8.2)') ' DEPTH  IS ', depth
!do is=0, cr_nscat
!write(*,'(a,12f7.2)') ' Depth  ',mmc_depth(ic, ie, is, :)
!enddo
!do is=0, cr_nscat
!write(*,'(a,12i4  )') ' Pair   ',mmc_pair(ic, ie, is, :)
!enddo
!read (*,*) is
!
mmc_target_corr(ic, ie, 0 , 0 ) = corr 
mmc_depth_def(ic) = (depth)
!write(*,*) ' mmc_depth_def ', ic, mmc_depth_def(ic)
!read (*,*) is
!
END SUBROUTINE mmc_set_group_occ
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
      mmc_rep_low = MIN(mmc_rep_low, ABS(a))
                                                                        
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
SUBROUTINE mmc_set_mode (ianz, cpara, lpara, werte, maxw, oopara, llopara) 
!+                                                                      
!     Sets MMC    mode                                                  
!-                                                                      
USE crystal_mod 
USE get_iscat_mod
USE rmc_mod 
USE mmc_mod 
!
USE modify_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE string_convert_mod
use take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) ::  MAXW 
CHARACTER(LEN=*)  , DIMENSION(MAXW), INTENT(INOUT) :: cpara !(MAXW) 
INTEGER           , DIMENSION(MAXW), INTENT(INOUT) :: lpara !(MAXW) 
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(INOUT) :: werte !(MAXW) 
character(len=PREC_STRING) :: oopara    ! Single string for get_optional_multi
integer                    :: llopara
!                                                                       
INTEGER :: ianz, imode=MC_MOVE_NONE, i , j
INTEGER :: is 
REAL(kind=PREC_DP) :: sump 
!
integer :: MAXWW
real(kind=PREC_DP) , dimension(:), allocatable :: wwerte   ! Calculated values
!
IF (ianz >= 1) THEN 
   CALL do_cap (cpara (1) ) 
   IF (cpara (1) (1:3)  == 'SHI') THEN 
      imode = MC_MOVE_DISP 
   ELSEIF (cpara (1) (1:3)  == 'SWD') THEN 
      imode = MC_MOVE_SWDISP 
   ELSEIF (cpara (1) (1:3)  == 'SWC') THEN 
      imode = MC_MOVE_SWCHEM 
   ELSEIF (cpara (1) (1:3)  == 'SWN') THEN 
!
! in addition interpret optional site:[] string
      imode = MC_MOVE_SWNEIG 
      MAXWW = cr_ncatoms
      allocate(wwerte(MAXWW))
      call get_optional_multi(MAXWw, oopara, llopara, wwerte, ianz)
      if(allocated(mmc_lsite)) deallocate(mmc_lsite)
      allocate(mmc_lsite(1:cr_ncatoms))
      mmc_lsite = .FALSE.
      do i=1, ianz
         j = nint(wwerte(i))
         if(j>0 .and. j<=cr_ncatoms) then
            mmc_lsite(j) = .TRUE.
         endif
      enddo
      deallocate(wwerte)
   ELSEIF (cpara (1) (1:3)  == 'INV') THEN 
      imode = MC_MOVE_INVDISP 
   ELSEIF (cpara (1) (1:3)  == 'ROT') THEN 
      imode = MC_MOVE_ROTATE  
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
SUBROUTINE mmc_run_multi (lout_feed_in)
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
use mmc_basic_mod
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
LOGICAL, INTENT(IN) :: lout_feed_in     ! Write output upon feedback?
!                                                                       
CHARACTER(LEN=24) :: c_energy (0:MC_N_ENERGY) 
!
REAL(kind=PREC_DP) :: start, zeit
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: disp !(3, 0:MAX_ATOM_ENV, 2) 
REAL(kind=PREC_DP), DIMENSION(3) :: idir, jdir
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: rdi ! (CHEM_MAX_COR) 
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: rdj ! (CHEM_MAX_COR) 
REAL(kind=PREC_DP) :: rel_cycl    ! how far are we in the desired number of cycles
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: patom ! Cooridnates for neighbors (3, 0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER(KIND=PREC_INT_LARGE) :: itry, igen, imodulus
INTEGER(KIND=PREC_INT_LARGE) :: iitry, iigen
INTEGER           , DIMENSION(  :,:), ALLOCATABLE :: iatom ! Indizes of neighbors      (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
LOGICAL           , DIMENSION(  :,:), ALLOCATABLE :: tatom ! True is atom is neighbor to atom iatom(0,*) (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: natom ! Number of neighbors       (MMC_MAX_CENT) 
INTEGER :: iacc_good, iacc_neut, iacc_bad 
INTEGER :: isel(MMC_MAX_ATOM) 
!INTEGER :: lbeg (3) 
INTEGER :: ic
INTEGER :: nocc
INTEGER :: i, natoms
INTEGER :: ncent 
INTEGER :: NALLOWED   ! Current size mmc_allowed
INTEGER :: zh, zm, zs 
LOGICAL, DIMENSION(3) :: old_chem_period
LOGICAL :: lserial    ! serial calculation if TRUE
logical :: lout_feed  ! Output on/off as combination of system and user settings
LOGICAL :: lfeed      ! Perform feedback algorithm
LOGICAL :: loop, laccept, done 
LOGICAL :: lout, lfinished
logical :: lmodulus
!                                                                       
REAL(kind=PREC_DP), DIMENSION(0:MC_N_ENERGY) :: e_old
REAL(kind=PREC_DP), DIMENSION(0:MC_N_ENERGY) :: e_new
!
INTEGER :: tid
INTEGER(KIND=PREC_INT_LARGE) :: nthreads, nnthreads
real(kind=PREC_DP), dimension(2)                    :: maxdev =(/0.0, 0.0/)
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
     'Unidirectional Corr     ',   &
     'Groupwise correlation   ' /
!
call mmc_initial(old_chem_period, itry, igen, iacc_good, iacc_neut, iacc_bad, & 
           done, loop)                 ! Perform initialization
!
!
ALLOCATE(disp (3, 0:MAX_ATOM_ENV, 2))
ALLOCATE(patom(3, 0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(rdi(CHEM_MAX_COR))
ALLOCATE(rdj(CHEM_MAX_COR))
ALLOCATE(iatom(   0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(tatom(   0:MAX_ATOM_ENV, MMC_MAX_CENT))
ALLOCATE(natom(                    MMC_MAX_CENT))
if(mmc_cor_energy (0, MC_DISP)) then
   call alloc_chem_dir(chem_ncor)
endif
!
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
lfeed = .FALSE.
CALL mmc_correlations (lout, 0.0D0, done, lfinished, lfeed, maxdev) 
lfeed = .TRUE.
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
check_conn: DO IC = 1, chem_ncor         ! Check is we have a connectivity, or Environment
   IF(chem_ctyp(ic) == CHEM_CON .or. chem_ctyp(ic) == CHEM_ENVIR) THEN    ! if so, we need serial computation
      lserial = .TRUE.                   ! Needs serious debugging...
      EXIT check_conn
   ENDIF
ENDDO check_conn
!                                                                       
!------ Loop will be started                                            
!                                                                       
lout_feed = lout_feed_in .and. mmc_out_feed    ! Combine system and user feedback settings
lfinished = .false.
done      = .false.
CALL mmc_correlations (.false., rel_cycl, done, lfinished, lfeed, maxdev) 
mmc_ini_corr  = mmc_ach_corr 
mmc_ini_sigm  = mmc_ach_sigm 
mmc_ini_angl  = mmc_ach_angl 
mmc_ini_sang  = mmc_ang_sigm 
mmc_ini_pairs = mmc_ach_pairs
if(lout_feed) then
   call mmc_correlation_write
endif
!write(output_io, '(a)') '---------------------------------------------'
!
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
!
if(mmc_feed_auto) then            ! Feedback cycles parameters set automatically
   imodulus=MAX(1_PREC_INT_LARGE, min(int(cr_natoms,PREC_INT_LARGE),mo_cyc/50/nthreads))
else                              ! Feedback cycles set by user
   imodulus=MAX(1_PREC_INT_LARGE, mo_feed/nthreads)
endif
!write(*,*) ' FEED ', min(cr_natoms,mo_cyc/50/nthreads), cr_natoms, mo_cyc/50/nthreads
!read(*,*) itry
!if(abs(maxdev(1))> 0.2) then
! imodulus = imodulus * 4
!endif
lmodulus = .true.
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
      DEALLOCATE(tatom)
      DEALLOCATE(natom)
      call alloc_mmc_pid  (1, 1, 1, 1)
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
   IF(ALLOCATED(tatom)) DEALLOCATE(tatom)
   IF(ALLOCATED(natom)) DEALLOCATE(natom)
   call alloc_mmc_pid  (1, 1, 1, 1)
   chem_period = old_chem_period
   RETURN      ! An error occured or CTRL-C
ENDIF
done = .FALSE. 
mmc_h_index = -1
mmc_m_index = -1
mmc_h_nfeed =  0
!write(*,*) ' *********************************'
!write(*,*) ' *********************************'
!write(*,*) ' *********************************'
!
IF(mmc_algo .EQV. MMC_CLASSIC) THEN           !Use the classical MMC algorithm
!
IF(nthreads > 1) THEN
   !$OMP PARALLEL PRIVATE(tid, isel, natoms, &
   !$OMP                  iatom, patom, tatom, natom,ncent, laccept,                 &
   !$OMP                  disp,             e_old, e_new)                            &
   !$OMP          SHARED(done)
   !$   tid = OMP_GET_THREAD_NUM()
   !$OMP DO SCHEDULE(DYNAMIC, mo_cyc/nthreads/32)
   parallel_loop: DO itry=1, mo_cyc     ! Do mmc in parallel
      IF(done) CYCLE parallel_loop    ! Quickly cycle to end if an error occuredd
      CALL    mmc_run_loop(tid, nthreads, igen, itry, natoms, &
                           iatom, patom, tatom, natom,ncent, laccept,                &
                           rdi, rdj,         e_old, e_new, done, loop,               &
                           iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                           lout_feed, lfeed, imodulus, lmodulus,                     &
                           NALLOWED, MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
!
   ENDDO parallel_loop
   !$OMP END DO NOWAIT
   !$ IF(tid==0) igen = igen*nthreads
   !$OMP END PARALLEL
   itry = igen
!
   tid = 0
   iigen = 0
   iitry = 1
   nnthreads = 1
!
   if(mmc_feed_auto) then            ! Feedback cycles parameters set automatically
      imodulus=MAX(1_PREC_INT_LARGE, min(int(cr_natoms,PREC_INT_LARGE),mo_cyc/50/nnthreads))
   else                              ! Feedback cycles set by user
      imodulus=MAX(1_PREC_INT_LARGE, mo_feed/nnthreads)
   endif
   done = .false.
   i = max(1_PREC_INT_LARGE, min(nthreads/4*cr_natoms, 2*mo_feed, mo_cyc))
   final_loop: DO iitry=1, i! nthreads/4*cr_natoms !2*mo_feed !nthreads**2*cr_natoms    ! Do mmc serially
      CALL    mmc_run_loop(tid,nnthreads, iigen, iitry, natoms, &
                           iatom, patom, tatom, natom,ncent, laccept,                &
                           rdi, rdj,         e_old, e_new, done, loop,               &
                           iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                           .false.  , lfeed, imodulus, lmodulus,                     &
                           NALLOWED, MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
!     IF(ier_num/=0 .OR. done) EXIT  final_loop
      IF(ier_num/=0          ) EXIT  final_loop
   ENDDO  final_loop
ELSE     ! Use nonparallel code
   serial_loop: DO itry=1, mo_cyc     ! Do mmc serially
      CALL    mmc_run_loop(tid, nthreads, igen, itry, &
                           natoms, &
                           iatom, patom, tatom, natom,ncent, laccept,                &
                           rdi, rdj,         e_old, e_new, done, loop,               &
                           iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                           lout_feed, lfeed, imodulus, lmodulus,                     &
                           NALLOWED, MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
      IF(ier_num/=0 .OR. done) EXIT serial_loop
   ENDDO serial_loop
ENDIF
ELSEIF(mmc_algo .EQV. MMC_GROWTH) THEN           !Use the growth MMC algorithm
      serial_grow: DO itry=1, mo_cyc     ! Do mmc serially
         CALL    mmc_run_grow(tid, nthreads, igen, itry, &
                              natoms, &
                              iatom, patom, tatom, natom,ncent, laccept,                &
                              rdi, rdj,         e_old, e_new, done, loop,               &
                              iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                              lout_feed, lfeed, imodulus,                               &
                              NALLOWED, MAX_ATOM_ENV, MMC_MAX_CENT, MMC_MAX_ATOM)
         IF(ier_num/=0 .OR. done) EXIT serial_grow
      ENDDO serial_grow
ENDIF
!
IF(ier_num/=0) THEN
   IF(ALLOCATED(patom)) DEALLOCATE(patom)
   IF(ALLOCATED(disp )) DEALLOCATE(disp )
   IF(ALLOCATED(rdi  )) DEALLOCATE(rdi  )
   IF(ALLOCATED(rdj  )) DEALLOCATE(rdj  )
   IF(ALLOCATED(iatom)) DEALLOCATE(iatom)
   IF(ALLOCATED(tatom)) DEALLOCATE(tatom)
   IF(ALLOCATED(natom)) DEALLOCATE(natom)
   call alloc_mmc_pid  (1, 1, 1, 1)
   chem_period = old_chem_period
   RETURN      ! An error occured or CTRL-C
ENDIF
!
lout_feed = lout_feed_in .and. mmc_out_final  ! Combine system and user final output setting
!
if(lout_feed) then
   write(output_io,*)
   write(output_io,'(a)')  ' ---------------------------------------------'
   write(output_io,'(a)')  ' --- multiple energy simulation is finished---'
   write(output_io,'(a)')  ' ---------------------------------------------'
   write(output_io,*)
   call mmc_correlation_write
endif
!                                                                       
!------ Loop finished                                                   
!                                                                       
IF(lout_feed) THEN
   WRITE (output_io, 3000)
   WRITE (output_io, 2000) igen, itry, iacc_good, iacc_neut, iacc_bad 
ENDIF
!     lout = .TRUE. 
lfinished = .TRUE.
lfeed     = .FALSE.   ! no feedback algorithm
done      = .TRUE.
CALL mmc_correlations (lout_feed, rel_cycl, done, lfinished, lfeed, maxdev) 
!                                                                       
!     Give average energy changes for the different energy terms        
!                                                                       
IF(lout_feed) WRITE ( output_io, 5000) 
DO i = 1, MC_N_ENERGY 
   if(n_e_av_p(i)>0 .or. n_e_av_z(i)>0 .or. n_e_av_m(i)>0) then
   IF(n_e_av_p(i) > 0) THEN 
      e_aver_p(i) = e_aver_p(i) / REAL(n_e_av_p(i) ) 
   ENDIF 
   IF(n_e_av_m(i) > 0) THEN 
      e_aver_m(i) = e_aver_m(i) / REAL(n_e_av_m(i) ) 
   ENDIF 
   IF(lout_feed) then
      WRITE(output_io, 5010) c_energy(i), n_e_av_m(i), e_aver_m(i),             & 
                                        n_e_av_z(i), n_e_av_p(i), e_aver_p(i)
   endif
   n_e_av_p (i) = 0 
   n_e_av_m (i) = 0 
   n_e_av_z (i) = 0 
   e_aver_p (i) = 0.0 
   e_aver_m (i) = 0.0 
endif
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
DEALLOCATE(tatom)
DEALLOCATE(natom)
call alloc_mmc_pid  (1, 1, 1, 1)
!
chem_period = old_chem_period
!                                                                       
 2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ')                                 
 2500 FORMAT (/,' --- Initial multiple energy configuration ---') 
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
                        iatom, patom, tatom, natom,ncent, laccept,                &
                        rdi, rdj,         e_old, e_new, done, loop,               &
                        iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                        lout_feed, lfeed, imodulus, lmodulus,                     &
                        NALLOWED, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!
USE crystal_mod
USE chem_mod
USE chem_menu
USE mc_mod
USE mmc_mod
use mmc_basic_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER                                                 , INTENT(INOUT) :: ncent
LOGICAL                           , INTENT(OUT) :: laccept
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdi ! (CHEM_MAX_COR) 
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdj ! (CHEM_MAX_COR) 
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_old
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_new
LOGICAL                           , INTENT(INOUT) :: done
LOGICAL                           , INTENT(in ) :: loop
INTEGER                           , INTENT(INOUT)  :: iacc_good
INTEGER                           , INTENT(INOUT)  :: iacc_neut
INTEGER                           , INTENT(INOUT)  :: iacc_bad
REAL(kind=PREC_DP)                              , INTENT(INOUT)  :: rel_cycl
LOGICAL                           , INTENT(IN ) :: lout_feed
LOGICAL                           , INTENT(IN ) :: lfeed
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(INout)  :: imodulus
LOGICAL                           , INTENT(INout ) :: lmodulus
!
INTEGER, DIMENSION(MMC_MAX_ATOM_L)              :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)                           :: is
INTEGER, DIMENSION(2, 3)                        :: iz
INTEGER, DIMENSION(3)                           :: iz1
INTEGER, DIMENSION(3)                           :: iz2
INTEGER                                         :: iselz
INTEGER                                         :: iselz2
LOGICAL :: valid_all
REAL(kind=PREC_DP)   , DIMENSION(3, 0:nthreads-1)                     :: posz !(3) = 0.0
REAL(kind=PREC_DP)   , DIMENSION(3, 0:nthreads-1)                     :: posz2 !(3) = 0.0
REAL(KIND=PREC_DP), DIMENSION(3,0:MAX_ATOM_ENV_l,2) :: disp
real(kind=PREC_DP), dimension(2)                    :: maxdev =(/0.0, 0.0/)
!integer :: i
!
IF(tid==0) then
   igen = igen + 1
!  if(rel_cycl>0.1 .and. lmodulus .and. maxdev(1)<0.1) then
!     imodulus = imodulus/5
!     lmodulus = .false.
!  endif
endif
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
!write(*,*) ' ISEL ', isel(1), ' : ', isel(2), ' | ', isel(3:4)
!if(mod(isel(1),2)==1) then
!write(*,*) ' First Atom at site 1 ' , isel(1), isel(1)+1, ' ISCAT ', cr_iscat(1,isel(1)), cr_iscat(1,isel(1)+1) 
!else
!write(*,*) ' First Atom at site 2 ' , isel(1)-1, isel(1), ' ISCAT ', cr_iscat(1,isel(1)-1), cr_iscat(1,isel(1)) 
!endif
!if(mod(isel(2),2)==1) then
!write(*,*) ' Scnd  Atom at site 1 ' , isel(2), isel(2)+1, ' ISCAT ', cr_iscat(1,isel(2)), cr_iscat(1,isel(2)+1) 
!else
!write(*,*) ' Scnd  Atom at site 2 ' , isel(2)-1, isel(2), ' ISCAT ', cr_iscat(1,isel(2)-1), cr_iscat(1,isel(2)) 
!endif
!write(*,*) '  mmc_pair 1 ', mmc_pair(1,MC_COORDNUM, cr_iscat(1,isel(1)),:)
!write(*,*) '  mmc_pair 2 ', mmc_pair(1,MC_COORDNUM, cr_iscat(1,isel(1)),:)
!write(*,*) '  MMC_PAIR11 ', mmc_pair(:,MC_COORDNUM, 1                  ,:)
!write(*,*) '  MMC_PAIR 1 ', any(mmc_pair(:,MC_COORDNUM, 1                  ,:)==-1)
!write(*,*) '  MMC_PAIR 2 ', any(mmc_pair(:,MC_COORDNUM, 2                  ,:)==-1)
!write(*,*) '  MMC_PAIR 3 ', any(mmc_pair(:,MC_COORDNUM, 3                  ,:)==-1)
!write(*,*) '  MMC_PAIR 4 ', any(mmc_pair(:,MC_COORDNUM, 4                  ,:)==-1)
!write(*,*)
!write(*,*) '  MMC_PAIR13 ', mmc_pair(:,MC_COORDNUM, 3                  ,:)
!write(,*,*) ' iscat', cr_iscat(1,isel(1)),cr_iscat(1,isel(2))
!write(*,*)
!read(*,*) i
!                                                                       
!--Calculate old energy                                          
!                                                                       
e_old      = 0.0    ! e_old(:)
!
valid_all = .FALSE.
!write(*,*)
!write(*,*) ' CALL ENERGIES OLD'
CALL mmc_energies(isel, is, iz, natoms, iatom, patom, tatom, natom, ncent, &
                  rdi, rdj, valid_all, e_old, CHEM_MAX_COR,         &
                  MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L, .true.)
!write(*,*) ' Old Energy ', e_old(MC_COORDNUM)
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
!write(*,*) ' CALL ENERGIES NEW'
   CALL mmc_energies(isel, is, iz, natoms, iatom, patom, tatom, natom, ncent, &
                     rdi, rdj, valid_all, e_new, CHEM_MAX_COR,         &
                     MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L, .false.)
!write(*,*) ' New Energy ', e_new(MC_COORDNUM)
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
!               normal and periodic mode
         CALL chem_apply_period(iselz, .TRUE.)
      ENDIF 
   ENDIF 
!if(e_new(MC_COORDNUM)-e_old(MC_COORDNUM)>0.0d0) then
!write(*,*) 'New - Old ', e_new(MC_COORDNUM)-e_old(MC_COORDNUM), laccept, &
!iacc_good, iacc_neut, iacc_bad
!!read(*,*) kk
!endif
   IF(ier_num/=0) THEN             ! Error, cycle to end of loop
      done = .TRUE.
      RETURN
   ENDIF
ENDIF 
!                                                                       
!     --End of modification of atoms if a proper "old" energy was found 
!                                                                       
!$OMP CRITICAL
   IF(tid==0) THEN
      IF(MOD(igen, imodulus)==0) THEN ! .AND. .NOT.done .AND. loop) THEN
!         done = .TRUE. 
         IF(lout_feed) WRITE (output_io, 2000) igen*nthreads, itry*nthreads, &
             iacc_good, iacc_neut, iacc_bad , rel_cycl
!                                                                       
!     ----New mmc_correlations for all energies                         
!                                                                       
      rel_cycl = REAL(igen)/REAL(mo_cyc)*REAL(NTHREADS)
      maxdev = 0.0
      CALL mmc_correlations (lout_feed, rel_cycl, done, .FALSE., lfeed, maxdev)
      if(mmc_feed_auto) then
         if(maxdev(1)<0.1) then
!           imodulus = max(min(imodulus-1, nint(imodulus*0.999)),nint(cr_natoms*0.05))
            imodulus = max(min(imodulus-1, nint(imodulus*0.950,PREC_INT_LARGE)), &
                           nint(cr_natoms*0.50, PREC_INT_LARGE))
         else
!           imodulus = max(imodulus+1,nint(imodulus*1.001))
            imodulus=MAX(1_PREC_INT_LARGE, min(int(cr_natoms,PREC_INT_LARGE),mo_cyc/50/nthreads))
         endif
      endif
!      imodulus = max(1000, nint(imodulus*0.90))
   ENDIF 
   IF(igen> mo_cyc/nthreads) THEN
      done = .TRUE.
   ENDIF
ENDIF 
!$OMP END CRITICAL
!
 2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,        &
     &          ' / ',I8,' / ',I8,'  MC moves ', f6.4)                                 
!
END SUBROUTINE mmc_run_loop
!
!*******************************************************************************
!
SUBROUTINE mmc_run_grow(tid, nthreads, igen, itry, &
                        natoms, &
                        iatom, patom, tatom, natom,ncent, laccept,                       &
                        rdi, rdj,         e_old, e_new, done, loop,               &
                        iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                        lout_feed, lfeed, imodulus,                               &
                        NALLOWED, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!
!-
!  The procedure starts at a random atom and performs a loop over all following atoms
!  The current atom is tested with respect to all energies. If a lower/equal energy
!  state exists, a second atom is chosen and tested.
!+
!
USE celltoindex_mod
USE crystal_mod
USE chem_mod
USE chem_menu
USE chem_neig_multi_mod
USE mc_mod
USE mmc_mod
use mmc_basic_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
INTEGER                                                 , INTENT(INOUT) :: ncent
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
LOGICAL                           , INTENT(OUT) :: laccept
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdi ! (CHEM_MAX_COR) 
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdj ! (CHEM_MAX_COR) 
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_old
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_new
LOGICAL                           , INTENT(INOUT) :: done
LOGICAL                           , INTENT(in ) :: loop
INTEGER                           , INTENT(INOUT)  :: iacc_good
INTEGER                           , INTENT(INOUT)  :: iacc_neut
INTEGER                           , INTENT(INOUT)  :: iacc_bad
REAL(kind=PREC_DP)                              , INTENT(INOUT)  :: rel_cycl
LOGICAL                           , INTENT(IN ) :: lout_feed
LOGICAL                           , INTENT(IN ) :: lfeed
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)  :: imodulus
!
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L,2) :: iatom_l      !Neighbors for isel(1:2)
INTEGER, DIMENSION(                  MMC_MAX_CENT_L,2) :: natom_l      ! number of neighbors
INTEGER, DIMENSION(                                 2) :: ncent_l      ! number of centers
!
!LOGICAL, PARAMETER :: LOUT=.FALSE.
!
!INTEGER :: MAXTEST
INTEGER, DIMENSION(MMC_MAX_ATOM_L)              :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)                           :: is
INTEGER, DIMENSION(2, 3)                        :: iz
INTEGER, DIMENSION(3)                           :: iz1
INTEGER, DIMENSION(3)                           :: iz2
INTEGER                                         :: iselz
INTEGER                                         :: iselz2
!
INTEGER, DIMENSION(2) :: iposit ! Original atom number at selected position 1,2
INTEGER, DIMENSION(2) :: iorig  ! Original atom type at selected position 1,2
INTEGER :: i !, j, k
INTEGER :: imode = 1  ! Select atom that could be improved most
!INTEGER :: latom
!INTEGER :: itest
INTEGER :: ichoose    ! Which atom within isel to choose
!INTEGER :: istart     ! Start for loop to search for new atom type
!INTEGER :: iend       ! End   for loop to search for new atom type
!INTEGER :: jmin       ! Atom type with minimum energy at selected atom position
INTEGER :: ic                             ! Loop index correlations
!INTEGER :: ie                             ! Loop index energy type
!INTEGER :: jc                             ! Current correlation to improve
!INTEGER :: js                             ! Loop index atom type 
!INTEGER :: ks                             ! Loop index atom type 
!INTEGER :: nt                             ! Number of targets for chosen energy
!INTEGER :: it                             ! Loop index over targets for chosen energy
INTEGER :: wic, wie                       ! Index of Correlation, energy, atom types at worst deviation
INTEGER :: icent                          ! Current center 
INTEGER :: iaccept                        ! DACounter if we want toaccept correlation changes
INTEGER, DIMENSION(2) :: wjks             ! Index of Correlation, energy, atom types at worst deviation
integer, dimension(0:2) :: isnei
LOGICAL :: valid_e
REAL(kind=PREC_DP) :: r1
REAL(kind=PREC_DP), DIMENSION(2,0:CHEM_MAX_COR) :: en_old    ! Old energy at position 1,2
REAL(kind=PREC_DP), DIMENSION(2,0:CHEM_MAX_COR) :: en_new    ! New energy at position 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  :: cold  ! old correlation at atom 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  :: cnew  ! old correlation at atom 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  ::ccold  ! old correlation at atom 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  ::ccnew  ! old correlation at atom 1,2
real(kind=PREC_DP), dimension(2)                    :: maxdev =(/0.0, 0.0/)
!integer, dimension(0:3, 0:3, 11) :: rbn_pneig = 0
!
REAL(kind=PREC_DP) :: damp = 1.0
!
damp = 0.01 + 0.99*exp(-4.0*rel_cycl)
!
mmc_move = MC_MOVE_SWCHEM      ! Needs Work
IF(tid==0) igen = igen + 1
!
IF(done) RETURN                 ! Quickly cycle to end if an error occuredd
!
!ALLOCATE(energies_old(0:cr_nscat,0:MC_N_ENERGY))    ! Should probably be done in main mmc_run
!
!                               ! Select correlation that needs work
CALL select_correlation(wic, wie, wjks)
!
natoms  = 1
cold    = 0.0
cnew    = 0.0
icent   = 1
ichoose = 1
isel = 0
!write(*,*) ' SELECT ISEL(1) ', wic, wie, wjks, ichoose, imode
CALL mmc_select_grow(wic, wie, wjks, ichoose, imode,                  &
                     isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                         iatom, patom, tatom, natom, ncent, &
                         laccept, loop,  cold, cnew, NALLOWED, &
                         MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
iatom_l(:,:,1) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,1) = natom
ncent_l(    1) = ncent
!write(*,*)
ichoose = 2
!write(*,*) ' SELECT ISEL(2) '
CALL mmc_select_grow(wic, wie, wjks, ichoose, imode,             &
                     isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                         iatom, patom, tatom, natom, ncent, &
                         laccept, loop,  cold, cnew, NALLOWED, &
                         MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
iatom_l(:,:,2) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,2) = natom
ncent_l(    2) = ncent
iposit    = isel
iorig(1)  = cr_iscat(1,isel(1))
iorig(2)  = cr_iscat(1,isel(2))
!write(*,*) ' OLD ic ', wic
!write(*,*) ' CHOICE ', isel(1), isel(2)
!write(*,*) ' Type   ', iorig
!write(*,*) ' Neig   ',          iatom_l(1:natom_l(1,1),1,1) , '|',         iatom_l(1:natom_l(1,2),1,2)
!write(*,*) ' Ntypes ', cr_iscat(1,iatom_l(1:natom_l(1,1),1,1)), '|',cr_iscat(1,iatom_l(1:natom_l(1,2),1,2))
!write(*,*) ' corr   ', cold, ' : ', SUM(cold)
!
!                                              CALCULATE OLD ENERGY
ccold = 0.0
en_old = 0.0
cold = 0.0
!WRITE(*,*) ' '
!WRITE(*,*) ' OLD STATUS', chem_ncor
!write(*,*) ' CHOICE   ',          isel(1) , '|',          isel(2)
!write(*,*) ' Type     ', cr_iscat(1,isel(1)), '|', cr_iscat(1,isel(2))
DO ic=1,chem_ncor                              ! Loop over all targets
   IF(mmc_cor_energy(ic,wie)) THEN
i = 1
!write(*,*) ' CALLING CORRELATION FOR ISEL(1) '
CALL chem_neighbour_multi(isel(i) , ic, iatom, patom, tatom, natom, &
                          ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
iatom_l(:,:,1) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,1) = natom
ncent_l(    1) = ncent
IF(mmc_cor_energy ( ic, MC_OCC)) THEN
   en_old(i,ic) =    mmc_occ_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cold(i)=cold(i) + en_old(i,ic)
ELSEIF(mmc_cor_energy ( ic, MC_UNI)) THEN
   en_old(i,ic) =    mmc_uni_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cold(i)=cold(i) + en_old(i,ic)
ENDIF
!
isnei(0) = iatom(1, 1)
isnei(1) = iatom(0, 1)
isnei(2) = iatom(2, 1)
!write(*,*) ' latom', latom - iinc
!write(*,*) ' ATOM ', isel(i),          ' : ',          isnei(0:2)          
!write(*,*) ' TYPES', cr_iscat(1,isel(i)),' : ', cr_iscat(1,isnei(0:2)          ), cnew(i)
i = 2
!write(*,*) ' CALLING CORRELATION FOR ISEL(2) '
CALL chem_neighbour_multi(isel(i) , ic, iatom, patom, tatom, natom, &
                          ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
iatom_l(:,:,2) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,2) = natom
ncent_l(    2) = ncent
IF(mmc_cor_energy ( ic, MC_OCC)) THEN
   en_old(i,ic) =    mmc_occ_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cold(i)=cold(i) + en_old(i,ic)
ELSEIF(mmc_cor_energy ( ic, MC_UNI)) THEN
   en_old(i,ic) =    mmc_uni_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cold(i)=cold(i) + en_old(i,ic)
ENDIF
!write(*,*) ' Neig     ',          iatom_l(1:natom_l(1,1),1,1) , '|',         iatom_l(1:natom_l(1,2),1,2)
!write(*,*) ' Ntypes   ', cr_iscat(1,iatom_l(1:natom_l(1,1),1,1)), '|',cr_iscat(1,iatom_l(1:natom_l(1,2),1,2))
!write(*,*) ' ic, corr ', ic, en_old(i,ic), '| ', cold, ' : ', SUM(cold)
   ENDIF
ENDDO
!
cr_iscat(1,isel(1)) = iorig(2)              ! Swap chemistry
cr_iscat(1,isel(2)) = iorig(1)
!
ccnew = 0.0
en_new = 0.0
cnew = 0.0
!WRITE(*,*) ' NEW STATUS'
!write(*,*) ' CHOICE   ',          isel(1) , '|',          isel(2)
!write(*,*) ' Type     ', cr_iscat(1,isel(1)), '|', cr_iscat(1,isel(2))
DO ic=1,chem_ncor                              ! Loop over all targets
   IF(mmc_cor_energy(ic,wie)) THEN
i = 1
!write(*,*) ' CALLING CORRELATION FOR ISEL(1) '
CALL chem_neighbour_multi(isel(i) , ic, iatom, patom, tatom, natom, &
                          ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
iatom_l(:,:,1) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,1) = natom
ncent_l(    1) = ncent
IF(mmc_cor_energy (wic, MC_OCC)) THEN
   en_new(i,ic) =    mmc_occ_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cnew(i)=cnew(i) + en_new(i,ic)
ELSEIF(mmc_cor_energy (wic, MC_UNI)) THEN
   en_new(i,ic) =    mmc_uni_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cnew(i)=cnew(i) + en_new(i,ic)
ENDIF
!
isnei(0) = iatom(1, 1)
isnei(1) = iatom(0, 1)
isnei(2) = iatom(2, 1)
!write(*,*) ' latom', latom - iinc
!write(*,*) ' ATOM ', isel(i),          ' : ',          isnei(0:2)          
!write(*,*) ' TYPES', cr_iscat(1,isel(i)),' : ', cr_iscat(1,isnei(0:2)          ), cnew(i)
i = 2
!write(*,*) ' CALLING CORRELATION FOR ISEL(2) '
CALL chem_neighbour_multi(isel(i) , ic, iatom, patom, tatom, natom, &
                          ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
iatom_l(:,:,2) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,2) = natom
ncent_l(    2) = ncent
IF(mmc_cor_energy (wic, MC_OCC)) THEN
   en_new(i,ic) =    mmc_occ_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cnew(i)=cnew(i) + ccnew(i)
ELSEIF(mmc_cor_energy (wic, MC_UNI)) THEN
   en_new(i,ic) =    mmc_uni_correl(isel, i, ic, iatom, icent,  natom, valid_e, &
                          MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
   cnew(i)=cnew(i) + en_new(i,ic)
ENDIF
!write(*,*) ' Neig     ',          iatom_l(1:natom_l(1,1),1,1) , '|',         iatom_l(1:natom_l(1,2),1,2)
!write(*,*) ' Ntypes   ', cr_iscat(1,iatom_l(1:natom_l(1,1),1,1)), '|',cr_iscat(1, iatom_l(1:natom_l(1,2),1,2))
!write(*,*) ' ic, corr ', ic, en_new(i,ic), '| ', cnew, ' : ', SUM(cnew)
   ENDIF
ENDDO
!
isnei(0) = iatom(1, 1)
isnei(1) = iatom(0, 1)
isnei(2) = iatom(2, 1)
!write(*,*) ' ATOM ', isel(i),          ' : ',          isnei(0:1)          
!write(*,*) ' TYPES', cr_iscat(1,isel(i)),' : ', cr_iscat(1,isnei(0:1)          ), cnew(i)
!write(*,*) ' NEW ic ', wic
!write(*,*) ' CHOICE ', isel(1), isel(2)
!write(*,*) ' Type   ', cr_iscat(1,isel(1)), cr_iscat(1,isel(2))
!write(*,*) ' Neig   ',          iatom_l(1:natom_l(1,1),1,1) , '|',         iatom_l(1:natom_l(1,2),1,2)
!write(*,*) ' Ntypes ', cr_iscat(1,iatom_l(1:natom_l(1,1),1,1)), '|',cr_iscat(1,iatom_l(1:natom_l(1,2),1,2))
!write(*,*) ' corr   ', cnew, ' : ', SUM(cnew)
!
laccept = .FALSE.
iaccept = 0
DO ic=1, chem_ncor
!IF(    mmc_target_corr( ic, wie, wjks(1), wjks(2))-           &
!       mmc_ach_corr   ( ic, wie, wjks(1), wjks(2))  <0) THEN       ! Target < Achieved
IF(    mmc_target_corr( ic, wie, 0      , 0      )-           &
       mmc_ach_corr   ( ic, wie, 0      , 0      )  <0) THEN       ! Target < Achieved
!  IF(SUM(cnew)<SUM(cold)) THEN                                    ! Correlation closer to target, accept
   IF(SUM(en_new(:,ic))<SUM(en_old(:,ic))) THEN                                    ! Correlation closer to target, accept
      iaccept = iaccept - 1
!     laccept = .TRUE.
!     iacc_good = iacc_good + 1
!write(*,*) ' GOOD: ACCEPTED T < A, CNEW < COLD ', SUM(en_new(:,ic)), SUM(en_old(:,ic)), ic, &
!mmc_target_corr( ic, wie, 0      , 0      ), mmc_ach_corr   ( ic, wie, 0      , 0      )
   ELSEIF(SUM(en_new(:,ic))==SUM(en_old(:,ic))) THEN                                    ! Correlation closer to target, accept
      iaccept = iaccept + 0
!     CALL RANDOM_NUMBER(r1)                    ! 
!     IF(r1<mmc_g_neut) THEN 
!        laccept = .TRUE.
!        iacc_neut = iacc_neut + 1
!write(*,*) ' NEUT: ACCEPTED T < A, CNEW = COLD ', SUM(en_new(:,ic)), SUM(en_old(:,ic)), ic, &
!mmc_target_corr( ic, wie, 0      , 0      ), mmc_ach_corr   ( ic, wie, 0      , 0      )
!     ENDIF
   ELSEIF(SUM(en_new(:,ic))>SUM(en_old(:,ic))) THEN                                    ! Correlation closer to target, accept
      iaccept = iaccept + 1
!     CALL RANDOM_NUMBER(r1)                    ! 
!     IF(r1<mmc_g_bad) THEN                          ! NEEDS WORK :: Boltzman !!
!        laccept = .TRUE.
!        iacc_bad  = iacc_bad  + 1
!write(*,*) ' BAD : ACCEPTED T < A, CNEW > COLD ', SUM(en_new(:,ic)), SUM(en_old(:,ic)), ic, &
!mmc_target_corr( ic, wie, 0      , 0      ), mmc_ach_corr   ( ic, wie, 0      , 0      )
!     ENDIF
   ENDIF
!ELSEIF(mmc_target_corr(wic, wie, wjks(1), wjks(2))-           &
!       mmc_ach_corr   (wic, wie, wjks(1), wjks(2))  >0) THEN       ! Target > Achieved
ELSEIF(mmc_target_corr(wic, wie, 0      , 0      )-           &
       mmc_ach_corr   (wic, wie, 0      , 0      )  >0) THEN       ! Target > Achieved
   IF(SUM(en_new(:,ic))>SUM(en_old(:,ic))) THEN                                    ! Correlation closer to target, accept
      iaccept = iaccept - 1
!     laccept = .TRUE.
!     iacc_good = iacc_good + 1
!write(*,*) ' GOOD: ACCEPTED T > A, CNEW > COLD ', SUM(en_new(:,ic)), SUM(en_old(:,ic)), ic, &
!mmc_target_corr( ic, wie, 0      , 0      ), mmc_ach_corr   ( ic, wie, 0      , 0      )
   ELSEIF(SUM(cnew)==SUM(en_old(:,ic))) THEN                                    ! Correlation closer to target, accept
      iaccept = iaccept + 0
!     CALL RANDOM_NUMBER(r1)                    ! 
!     IF(r1<0.25) THEN 
!        laccept = .TRUE.
!        iacc_neut = iacc_neut + 1
!write(*,*) ' NEUT: ACCEPTED T > A, CNEW = COLD ', SUM(en_new(:,ic)), SUM(en_old(:,ic)), ic, &
!mmc_target_corr( ic, wie, 0      , 0      ), mmc_ach_corr   ( ic, wie, 0      , 0      )
!     ENDIF
   ELSEIF(SUM(en_new(:,ic))<SUM(en_old(:,ic))) THEN                                    ! Correlation closer to target, accept
      iaccept = iaccept + 1
!     CALL RANDOM_NUMBER(r1)                    ! 
!     IF(r1<0.05) THEN                          ! NEEDS WORK :: Boltzman !!
!        laccept = .TRUE.
!        iacc_bad  = iacc_bad  + 1
!write(*,*) ' BAD : ACCEPTED T > A, CNEW < COLD ', SUM(en_new(:,ic)), SUM(en_old(:,ic)), ic, &
!mmc_target_corr( ic, wie, 0      , 0      ), mmc_ach_corr   ( ic, wie, 0      , 0      )
!     ENDIF
   ENDIF
ENDIF
ENDDO
!
laccept = .FALSE.
IF(    iaccept < 0) THEN       ! Accept move , most correlations closer to target
   laccept = .TRUE.
   iacc_good = iacc_good + 1
!   write(*,*) ' GOOD: ACCEPTED ', iaccept
ELSEIF(iaccept == 0) THEN  ! Accept move with neutral probability
   CALL RANDOM_NUMBER(r1)                    ! 
   IF(r1<mmc_g_neut) THEN 
      laccept = .TRUE.
      iacc_neut = iacc_neut + 1
!   write(*,*) ' NEUT: ACCEPTED ', iaccept
   ELSE
      laccept = .FALSE.
!   write(*,*) ' NEUT: REJECTED ', iaccept
   ENDIF
ELSEIF(iaccept >  0) THEN  ! Accept move with bad     probability
   CALL RANDOM_NUMBER(r1)                    ! 
   IF(r1<mmc_g_bad ) THEN 
      laccept = .TRUE.
      iacc_bad  = iacc_bad  + 1
!   write(*,*) ' BAD : ACCEPTED ', iaccept
   ELSE
      laccept = .FALSE.
!   write(*,*) ' BAD : REJECTED ', iaccept
   ENDIF
ENDIF
!
IF(laccept) THEN                                ! Move accepted update correlation
!P!   DO ic=1, chem_ncor                           ! Update pneig for all correlatiosn with OCC or UNI energy
!P!      IF(mmc_cor_energy(ic, mc_OCC) .OR. mmc_cor_energy(ic, mc_UNI)) THEN    ! to include side effects of current change
!P!do js=0,cr_nscat
!P!write(*,*) ' ACCEPTED:  PAIRS OLD ', mmc_pneig(js,0:cr_nscat, ic)
!P!enddo
!P!write(*,*) ' k=1 ncent ', ncent_l(1), '|', natom_l(:,1), ' IC ', ic
!P!write(*,*) ' k=2 ncent ', ncent_l(2), '|', natom_l(:,2), ' IC ', ic
!P!   DO k=1, 2
!P!      js = cr_iscat(1,isel(k))
!P!      DO i=1, ncent_l(k)
!P!         DO j=1, natom_l(i,k)
!P!            ks = cr_iscat(1,iatom_l(j,i,k))
!P!write(*,*) 'k, i, ks ', k,i,ks, ' iorig(k), js, ks ',iorig(k), js, ks, 'O ', mmc_pneig(iorig(k),ks, ic), mmc_pneig(js      ,ks, ic)
!P!!           IF(mmc_pair(wic, wie, iorig(k),ks)/=0) mmc_pneig(iorig(k),ks, wic) = mmc_pneig(iorig(k),ks, wic) - 1
!P!!           IF(mmc_pair(wic, wie, js      ,ks)/=0) mmc_pneig(js      ,ks, wic) = mmc_pneig(js      ,ks, wic) + 1
!P!                                                   mmc_pneig(iorig(k),ks, ic) = mmc_pneig(iorig(k),ks, ic) - 1
!P!                                                   mmc_pneig(js      ,ks, ic) = mmc_pneig(js      ,ks, ic) + 1
!P!            i = iatom_l(j,i,k)
!P!            CALL chem_neighbour_multi(i , wic, iatom, patom, natom, &
!P!                          ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
!P!            DO i=1, natom(1)   ! CENTERS !!!!!
!P!               IF(iatom(i,1)==isel(k)) Found current selected atom as neighbor to a neighbor
!P!            ENDDO
!P!write(*,*) 'k, i, ks ', k,i,ks, ' iorig(k), js, ks ',iorig(k), js, ks, 'N ', mmc_pneig(iorig(k),ks, ic), mmc_pneig(js      ,ks, ic)
!P!         ENDDO
!P!      ENDDO
!P!   ENDDO
!P!do js=0,cr_nscat
!P!write(*,*) ' ACCEPTED:  PAIRS NEW ', mmc_pneig(js,0:cr_nscat, ic)
!P!enddo
!P!rbn_pneig = mmc_pneig
!P!if(minval(mmc_pneig(1:3,1:3,ic)) <0) then
!P!read(*,*) js
!P!endif
!P!!  IF(mmc_cor_energy(ic, MC_OCC)) THEN         ! Update achieved correlations for SRO occupational correlations
!P!!     CALL mmc_correlations_occ(ic, mmc_pneig, rel_cycl, damp, LOUT, MAXSCAT, CHEM_MAX_COR)
!P!!  ENDIF
!P!write(*,*) ' MOVE ACCEPTED IC WIE',  wic, wie
!P!   IF(mmc_cor_energy(ic, MC_UNI)) THEN         ! Update achieved correlations for Unidirectional energy
!P!!Q      CALL mmc_correlations_uni(ic, mmc_pneig, rel_cycl, damp, LOUT, MAXSCAT, CHEM_MAX_COR)
!P!   ENDIF
!P!      ENDIF
!P!   ENDDO
!P!      CALL mmc_correlations (.true.   , rel_cycl, done, .FALSE.)
!P!if(MAXVAL(ABS(rbn_pneig -mmc_pneig))>2) THEN
!P!do ic=1,chem_ncor
!P!write(*,*) ' CORRELATION ', ic
!P!do js=0,cr_nscat
!P!write(*,*) ' mmc, rbn ', mmc_pneig(js,0:cr_nscat, ic), ' | ',rbn_pneig(js,0:cr_nscat, ic)
!P!enddo
!P!enddo
!P!
!P!read(*,*) js
!P!ENDIF
ELSE                                            ! Move not accepted restore atom types
   cr_iscat(1,isel(1)) = iorig(1)
   cr_iscat(1,isel(2)) = iorig(2)
!write(*,*) ' NOT   ACCEPTED      ; CNEW ? COLD ', SUM(cnew), SUM(cold)
!write(*,*) ' CHOICE ', isel(1), isel(2)
!write(*,*) ' Type   ', cr_iscat(1,isel(1)), cr_iscat(1,isel(2))
ENDIF
!                                                                       
!------ --- Test and accept/reject move                                 
!                                                                       
!e_old = en_old(1,:) + en_old(2,:)
!e_new = en_new(1,:) + en_new(2,:)
!         CALL mmc_test_multi (iacc_good, iacc_neut, iacc_bad,  &
!              e_new, e_old,  laccept)                                                 
!IF(.NOT.laccept) THEN
!   cr_iscat(1,iposit(1)) = iorig(1)
!   cr_iscat(1,iposit(2)) = iorig(2)
!ENDIF
!write(*,*) ' FINAL    ', SUM(e_old), sum(e_new), laccept, iacc_good, iacc_neut, iacc_bad
!write(*,*) ' Atoms pos 1', iposit(1)-1          , iposit(1)          , iposit(1)+1 
!write(*,*) ' Types pos 1', cr_iscat(1,iposit(1)-1), cr_iscat(1,iposit(1)), cr_iscat(1,iposit(1)+1)
!write(*,*) ' Atoms pos 2', iposit(2)-1          , iposit(2)          , iposit(2)+1 
!write(*,*) ' Types pos 2', cr_iscat(1,iposit(2)-1), cr_iscat(1,iposit(2)), cr_iscat(1,iposit(2)+1)
!write(*,*)
!
!DEALLOCATE(energies_old)                    ! Should probably be done in main mmc_run
!
!     --End of modification of atoms if a proper "old" energy was found 
!
!$OMP CRITICAL
IF(tid==0) THEN
   IF(MOD(igen, imodulus)==0) THEN ! .AND. .NOT.done .AND. loop) THEN
!         done = .TRUE. 
      IF(lout_feed) WRITE (output_io, 2000) igen*nthreads, itry*nthreads, &
          iacc_good, iacc_neut, iacc_bad 
!                                                                       
!     ----New mmc_correlations for all energies                         
!                                                                       
      rel_cycl = REAL(igen)/REAL(mo_cyc)*REAL(NTHREADS)
      CALL mmc_correlations (lout_feed, rel_cycl, done, .FALSE., lfeed, maxdev)
   ENDIF 
   IF(igen> mo_cyc/nthreads) THEN
      done = .TRUE.
   ENDIF
ENDIF 
!$OMP END CRITICAL
!
2000 FORMAT (/,' Gen: ',I10,' try: ',I10,' acc: (g/n/b): ',I8,        &
             ' / ',I8,' / ',I8,'  MC moves ')                                 
!
END SUBROUTINE mmc_run_grow
!
!*******************************************************************************
!
SUBROUTINE mmc_run_oswald(tid, nthreads, igen, itry, &
                        natoms, &
                        iatom, patom, tatom, natom,ncent, laccept,                       &
                        rdi, rdj,         e_old, e_new, done, loop,               &
                        iacc_good, iacc_neut, iacc_bad, rel_cycl,                 &
                        lout_feed, lfeed, imodulus,                               &
                        NALLOWED, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!
!-
!  The procedure starts at a random atom that cannot be improved any further.
!  The neighbors are searched (in turn) until one of these neighbors could be 
!  improved. A complementary atom is chosen that can be improved as well.
!+
!
USE celltoindex_mod
USE crystal_mod
USE chem_mod
USE chem_menu
USE chem_neig_multi_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER                                                 , INTENT(INOUT) :: ncent
LOGICAL                           , INTENT(in ) :: laccept
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdi ! (CHEM_MAX_COR) 
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(INOUT) :: rdj ! (CHEM_MAX_COR) 
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_old
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY) , INTENT(INOUT) :: e_new
LOGICAL                           , INTENT(INOUT) :: done
LOGICAL                           , INTENT(in ) :: loop
INTEGER                           , INTENT(INOUT)  :: iacc_good
INTEGER                           , INTENT(INOUT)  :: iacc_neut
INTEGER                           , INTENT(INOUT)  :: iacc_bad
REAL(kind=PREC_DP)                              , INTENT(INOUT)  :: rel_cycl
LOGICAL                           , INTENT(IN ) :: lout_feed
LOGICAL                           , INTENT(IN ) :: lfeed
INTEGER(KIND=PREC_INT_LARGE)      , INTENT(IN)  :: imodulus
!
INTEGER, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L,2) :: iatom_l      !Neighbors for isel(1:2)
INTEGER, DIMENSION(                  MMC_MAX_CENT_L,2) :: natom_l      ! number of neighbors
INTEGER, DIMENSION(                                 2) :: ncent_l      ! number of centers
!
INTEGER, DIMENSION(MMC_MAX_ATOM_L)              :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)                           :: is
INTEGER, DIMENSION(2, 3)                        :: iz
INTEGER, DIMENSION(3)                           :: iz1
INTEGER, DIMENSION(3)                           :: iz2
INTEGER                                         :: iselz
INTEGER                                         :: iselz2
!
INTEGER :: imode = 2  ! Select atom that cannot be improved any further
INTEGER :: ichoose    ! Which atom within isel to choose
INTEGER :: wic, wie                       ! Index of Correlation, energy, atom types at worst deviation
INTEGER :: icent                          ! Current center 
INTEGER, DIMENSION(2) :: wjks             ! Index of Correlation, energy, atom types at worst deviation
REAL(KIND=PREC_DP), DIMENSION(2)  :: cold  ! old correlation at atom 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  :: cnew  ! old correlation at atom 1,2
!
REAL(kind=PREC_DP) :: damp = 1.0
!
!write(*,*) ' IN MMC OSWALD '
damp = 0.01 + 0.99*exp(-4.0*rel_cycl)
!
mmc_move = MC_MOVE_SWCHEM      ! Needs Work
IF(tid==0) igen = igen + 1
!
IF(done) RETURN                 ! Quickly cycle to end if an error occuredd
!
!ALLOCATE(energies_old(0:cr_nscat,0:MC_N_ENERGY))    ! Should probably be done in main mmc_run
!
!                               ! Select correlation that needs work
CALL select_correlation(wic, wie, wjks)
!
natoms  = 1
cold    = 0.0
cnew    = 0.0
icent   = 1
ichoose = 1
isel = 0
!write(*,*) ' SELECT ISEL(1) '
CALL mmc_select_grow(wic, wie, wjks, ichoose, imode,                  &
                     isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                         iatom, patom, tatom, natom, ncent, &
                         laccept, loop,  cold, cnew, NALLOWED, &
                         MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!write(*,*) ' BACK FROM mmc_select_grow ', isel(1), ier_num, ier_typ
!write(*,*) ' NCENT, Neighbors ', ncent, '|',natom
!write(*,*) ' IATOM            ', iatom(0:natom(1),1)
iatom_l(:,:,1) = iatom                         ! Copy neighbors of isel(1)
natom_l(  :,1) = natom
ncent_l(    1) = ncent
!write(*,*) ' OLD ic ', wic
!write(*,*) ' CHOICE ', isel(1), isel(2)
!write(*,*) ' Type   ', iorig
!write(*,*) ' Neig   ',          iatom_l(1:natom_l(1,1),1,1) , '|',         iatom_l(1:natom_l(1,2),1,2)
!write(*,*) ' Ntypes ', cr_iscat(1,iatom_l(1:natom_l(1,1),1,1)), '|',cr_iscat(1,iatom_l(1:natom_l(1,2),1,2))
!write(*,*) ' corr   ', cold, ' : ', SUM(cold)
!write(*,*) ' ENTER NUMBER '
!read(*,*) ic
!
END SUBROUTINE mmc_run_oswald
!
!*****7*****************************************************************
!
SUBROUTINE select_correlation(wic, wie, wjks)
!
USE crystal_mod
USE chem_mod
USE mc_mod
USE mmc_mod
!
USE precision_mod
!
IMPLICIT NONE
!
INTEGER              , INTENT(OUT) :: wic    ! Selected coordination
INTEGER              , INTENT(OUT) :: wie    ! Selected energy
INTEGER, DIMENSION(2), INTENT(OUT) :: wjks   ! Selected atom types
!
INTEGER :: ic                             ! Loop index correlations
INTEGER :: ie                             ! Loop index energy type
!INTEGER :: jc                             ! Current correlation to improve
INTEGER :: js                             ! Loop index atom type 
INTEGER :: ks                             ! Loop index atom type 
INTEGER :: nt                             ! Number of targets for chosen energy
INTEGER :: it                             ! Loop index over targets for chosen energy
!
REAL(KIND=PREC_DP) :: deviation           ! Max relative deviation (target-achieved)/target
REAL(KIND=PREC_DP) :: dev                 ! Relative deviation
REAL(KIND=PREC_DP) :: r1                  ! Random numberation
!
wic  = 0
wie  = 0
wjks = 0
!ie = MC_UNI    ! NEEDS WORK
deviation = 0.0
dev       = 0.0
!write(*,*) ' IC, cr_nscat, ie ', chem_ncor, ie, cr_nscat
!write(*,*) ' mmc_pair ', mmc_pair(1,ie,:,:)
!write(*,*) ' TARGET   ', mmc_target_corr(1,ie,:,:)
!write(*,*) ' ACHIEVED ', mmc_ach_corr   (1,ie,:,:)
DO ic=1,chem_ncor
   DO ie=1, MC_N_ENERGY
      IF(mmc_cor_energy(ic, ie)) THEN    ! A target is set for this energy
   DO js=0,cr_nscat-1
      DO ks=js+1,cr_nscat
         IF(mmc_pair(ic,ie,js,ks)/=0) THEN
            IF(mmc_target_corr(ic, ie, js, ks)/=0.0) THEN
              dev = ABS((mmc_target_corr(ic, ie, js, ks)-      &
                         mmc_ach_corr (ic, ie, js, ks)   )/    &
                         mmc_target_corr(ic, ie, js, ks)    )
            ELSE
               dev = ABS(mmc_target_corr(ic, ie, js, ks))
            ENDIF
         ENDIF
         IF(dev>deviation) THEN
            deviation = dev
            wic = ic
            wie = ie
            wjks(1) = js
            wjks(2) = ks
         ENDIF
      ENDDO
   ENDDO
      ENDIF
   ENDDO
ENDDO
CALL RANDOM_NUMBER(r1)                    ! 
wic = int(r1*chem_ncor) + 1               ! Choose a correlation at random
nt = 0
DO it = 1, MC_N_ENERGY
   IF(mmc_cor_energy(wic,it)) THEN
      nt = nt + 1                         ! Accumulate number of targets for this correlation
      wie = it                            ! Make sure wie is set to at least one energy
   ENDIF
ENDDO
CALL RANDOM_NUMBER(r1)                    ! 
ie = INT(r1*nt) + 1                       ! Choose one of the NT energies
nt = 0
find_e: DO it = 1, MC_N_ENERGY            ! Loop to find the nt's energy for this correlations
   IF(mmc_cor_energy(wic,it)) nt = nt + 1
   IF(nt==ie) THEN
      wie = it                            ! Replace wie with current index
      EXIT find_e
   ENDIF
ENDDO find_e
!
DO js=0,cr_nscat
   IF(mmc_pair(wic, wie, js,js)==1) THEN
      wjks(1) = js
   ELSEIF(mmc_pair(wic, wie, js,js)==2) then
      wjks(2) = js
   ENDIF
ENDDO
!IF(wic==1) THEN
!   wjks(1) = 1
!   wjks(2) = 2
!ELSEIF(wic==2) THEN
!   wjks(1) = 2
!   wjks(2) = 3
!ELSEIF(wic==3) THEN
!   wjks(1) = 3
!   wjks(2) = 1
!ENDIF
!write(*,*) ' WORST CASE ', wic, wie, wjks, deviation
!
END SUBROUTINE select_correlation
!
!*****7*****************************************************************
!
SUBROUTINE mmc_select_grow(wic, wie, wjks, ichoose, imode,                &
                           isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                         iatom, patom, tatom, natom, ncent, &
                         laccept, loop, cold, cnew, NALLOWED,             &
                         MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!-
!  Choose a move that will improve an achieved correlation towards
!  the target value.
!  All other correlations around the central atom must in sum
!  either improve or be neutral as well
!+
!
USE crystal_mod
USE chem_mod
USE chem_neig_multi_mod
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
INTEGER                           , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                           , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                           , INTENT(IN)  :: NALLOWED   ! Current size mmc_allowed
INTEGER                           , INTENT(IN)  :: wic             ! Index of Correlation, energy, atom types at worst deviation
INTEGER                           , INTENT(IN)  :: wie             ! Index of Correlation, energy, atom types at worst deviation
INTEGER, DIMENSION(2)             , INTENT(IN)  :: wjks            ! Index of Correlation, energy, atom types at worst deviation
INTEGER                           , INTENT(IN)  :: ichoose         ! Which atom to choose within isel
INTEGER                           , INTENT(IN)  :: imode           ! How to chose an atom 
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(OUT) :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)             , INTENT(in ) :: is
INTEGER                           , INTENT(in ) :: iselz
INTEGER                           , INTENT(in ) :: iselz2
INTEGER                           , INTENT(in ) :: natoms
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
INTEGER                                                 , INTENT(INOUT) :: ncent
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER, DIMENSION(2, 3)          , INTENT(in ) :: iz
INTEGER, DIMENSION(3)             , INTENT(in ) :: iz1
INTEGER, DIMENSION(3)             , INTENT(in ) :: iz2
LOGICAL                           , INTENT(in ) :: laccept
LOGICAL                           , INTENT(in ) :: loop
REAL(KIND=PREC_DP), DIMENSION(2)  , INTENT(INOUT) :: cold  ! old correlation at atom 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  , INTENT(INOUT) :: cnew  ! old correlation at atom 1,2
!
INTEGER, PARAMETER :: ISANY     = 0
INTEGER, PARAMETER :: ISWORST   = 1
INTEGER, PARAMETER :: ISPERFECT = 2
!
REAL(kind=PREC_DP) :: r1
INTEGER :: ngrand                         ! End   atom number for random search
INTEGER :: istart                         ! Start atom number for random search
INTEGER :: icounter                       ! End   atom number for random search
INTEGER :: inumber                        ! End   atom number for random search
INTEGER :: iinc                           ! Increment         for random search
INTEGER :: latom                          ! Current atom number
!INTEGER :: ic                             ! Loop index correlations
!INTEGER :: ie                             ! Loop index energy type
!INTEGER :: jc                             ! Current correlation to improve
!INTEGER :: js                             ! Loop index atom type 
!INTEGER :: ks                             ! Loop iin x atom type 
!INTEGER :: i                              ! Dummy loop index
INTEGER :: ia                             ! Dummy loop index
INTEGER :: icent                          ! Dummy loop index
integer, dimension(0:2) :: isnei
!
LOGICAL :: valid_e
!
ngrand = 0
grand: DO
   ngrand = ngrand + 1
!
   iinc = 1
   inumber = cr_icc(1)
!
   CALL RANDOM_NUMBER(r1)                    ! 
!
IF(r1<1./2.) THEN                         ! Increment along x-axis
   IF(cr_icc(1) > 1) THEN
      CALL RANDOM_NUMBER(r1)
      iinc = NINT(SIGN(1.0D0, r1-0.5D0))*cr_ncatoms  ! Randomly choose increment -1 or +1 * atoms per unit cell
      inumber = cr_icc(1)
   ELSE
      CYCLE grand   
   ENDIF
ELSEIF(r1<2./2+1) THEN                      ! Increment along y-axis
   IF(cr_icc(2) > 1) THEN
      CALL RANDOM_NUMBER(r1)
      iinc = NINT(SIGN(1.0D0, r1-0.5D0))*cr_ncatoms*cr_icc(1)  ! Randomly choose increment -1 or +1 * atoms per unit cell
      inumber = cr_icc(2)
   ELSE
      CYCLE grand   
   ENDIF
!ELSE                                      ! Increment along z-axis
!   CYCLE grand   
ENDIF
!iinc   = NINT(SIGN(1.0, r1-0.5))          ! Randomly choose increment -1 or +1
!CALL RANDOM_NUMBER(r1)                    ! 
istart = INT(r1*cr_natoms) + 1            ! Initalize loop over all atoms at random start
icounter = 0
!iend   = istart + iinc*(cr_natoms-1)      ! End after full list of atoms
latom = istart
!write(*,*) ' LOOP       ', istart, inumber, iinc, ngrand
search_ini: DO                            ! Loop until we find an atom that could be improved
   isel(ichoose) = MOD(latom-1+cr_natoms, cr_natoms) + 1
!write(*,*) ' latom isel(1), wjs ', latom, isel(ichoose), wjks(ichoose)
   latom = latom + iinc
   icounter = icounter + 1
   IF(icounter==inumber) THEN
      EXIT search_ini   
   ENDIF
   IF(cr_iscat(1,isel(ichoose))==wjks(ichoose)) THEN          ! Found atom of proper type
      CALL chem_neighbour_multi(isel(ichoose) ,wic, iatom, patom, tatom, natom, &
                                ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
      ia = ichoose
      icent = 1
      IF(mmc_cor_energy (wic, MC_OCC)) THEN
         cold(ichoose)=mmc_occ_correl(isel, ia,wic, iatom, icent,     &
                                      natom, valid_e, MAX_ATOM_ENV_L, &
                                      MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
      ELSEIF(mmc_cor_energy (wic, MC_UNI)) THEN
         cold(ichoose)=mmc_uni_correl(isel, ia,wic, iatom, icent,     &
                                      natom, valid_e, MAX_ATOM_ENV_L, &
                                      MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
      ENDIF
      isnei(0) = iatom(2, 1)
      isnei(1) = iatom(0, 1)
      isnei(2) = iatom(1, 1)
!write(*,*) ' latom', latom - iinc
!write(*,*) ' ATOM ', isel(ichoose),          ' : ',          isnei(0:2)          
!write(*,*) ' TYPES', cr_iscat(1,isel(ichoose)),' : ', cr_iscat(1,isnei(0:2)          ), cold(ichoose)
!write(*,*) ' Targ, ach, diff ', mmc_target_corr(wic, wie, wjks(1), wjks(2)), &
!mmc_ach_corr   (wic, wie, wjks(1), wjks(2)), &
!mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2)) , ' T ? A '
      IF(mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2))<0) THEN
!write(*,*) ' Targ, ach, diff ', mmc_target_corr(wic, wie, wjks(1), wjks(2)), &
!mmc_ach_corr   (wic, wie, wjks(1), wjks(2)), &
!mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2)) , ' T < A '
         IF(cold(ichoose)>mmc_target_corr(wic, wie, wjks(1), wjks(2))) THEN
!write(*,*) '1 IMPROVEMENT POSSIBLE if cnew becomes smaller            T < corr'
            IF(imode==ISWORST .OR. IMODE==ISANY) THEN 
               EXIT search_ini
            ENDIF
         ELSEIF(cold(ichoose)<mmc_target_corr(wic, wie, wjks(1), wjks(2))) THEN
!write(*,*) '2 IMPROVEMENT POSSIBLE if cnew becomes even smaller       T > corr'
            IF(IMODE==ISANY) THEN 
               CALL RANDOM_NUMBER(r1)                    ! 
               IF(r1<mmc_g_rate) EXIT search_ini
            ENDIF
         ELSE
            IF(IMODE==ISPERFECT) THEN 
!write(*,*) '3 ALREADY PERFECT '
               EXIT search_ini
            ENDIF
         ENDIF
      ELSEIF(mmc_target_corr(wic, wie, wjks(1), wjks((2)))-mmc_ach_corr (wic, wie, wjks(1), wjks(2))>0) THEN
!write(*,*) ' Targ, ach, diff ', mmc_target_corr(wic, wie, wjks(1), wjks(2)), &
!mmc_ach_corr   (wic, wie, wjks(1), wjks(2)), &
!mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2)) , ' T > A '
         IF(cold(ichoose)>mmc_target_corr(wic, wie, wjks(1),wjks(2))) THEN
!write(*,*) '4 IMPROVEMENT POSSIBLE if cnew becomes even larger        T < corr'
         ELSEIF(cold(ichoose)<mmc_target_corr(wic, wie, wjks(1), wjks(2))) THEN
!write(*,*) '5 IMPROVEMENT POSSIBLE if cnew becomes larger             T > corr'
exit search_ini
         ELSE
!write(*,*) '6 ALREADY PERFECT '
         ENDIF
      ELSE
!write(*,*) '7 ALREADY PERFECT '
      ENDIF
!exit search_ini
   ENDIF
ENDDO search_ini
   IF(isel(ichoose) /=0) EXIT grand
   IF(ngrand>10   ) THEN
!write(*,*) ' NO MOVE FOUND ', iinc, istart, inumber
      ier_num=-2
      ier_typ = ER_MMC
      RETURN
   ENDIF
ENDDO grand
!if(isel(ichoose) == 0) then
!  write(*,*) ' WRONG MOVE B' , iinc, istart, inumber, ngrand
!  write(*,*) ' WRONG MOVE B' , isel, ' : ' , ichoose, latom
!endif
!
END SUBROUTINE mmc_select_grow
!
!*****7*****************************************************************
!
SUBROUTINE mmc_select_emin(wic, wie, wjks, ichoose, imode,                &
                           isel, is, iz, iz1, iz2, iselz, iselz2, natoms, &
                         iatom, patom, tatom, natom, ncent, &
                         laccept, loop, cold, cnew, NALLOWED,             &
                         MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)
!-
!  Choose an atom that is at best 
!  the target value.
!+
!
USE crystal_mod
USE chem_mod
USE chem_neig_multi_mod
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
INTEGER                           , INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER                           , INTENT(IN) :: MMC_MAX_CENT_L
INTEGER                           , INTENT(IN)  :: NALLOWED   ! Current size mmc_allowed
INTEGER                           , INTENT(IN)  :: wic             ! Index of Correlation, energy, atom types at worst deviation
INTEGER                           , INTENT(IN)  :: wie             ! Index of Correlation, energy, atom types at worst deviation
INTEGER, DIMENSION(2)             , INTENT(IN)  :: wjks            ! Index of Correlation, energy, atom types at worst deviation
INTEGER                           , INTENT(IN)  :: ichoose         ! Which atom to choose within isel
INTEGER                           , INTENT(IN)  :: imode           ! How to chose an atom 
INTEGER, DIMENSION(MMC_MAX_ATOM_L), INTENT(OUT) :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)             , INTENT(in ) :: is
INTEGER                           , INTENT(in ) :: iselz
INTEGER                           , INTENT(in ) :: iselz2
INTEGER                           , INTENT(in ) :: natoms
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
INTEGER                                                 , INTENT(INOUT) :: ncent
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER, DIMENSION(2, 3)          , INTENT(in ) :: iz
INTEGER, DIMENSION(3)             , INTENT(in ) :: iz1
INTEGER, DIMENSION(3)             , INTENT(in ) :: iz2
LOGICAL                           , INTENT(in ) :: laccept
LOGICAL                           , INTENT(in ) :: loop
REAL(KIND=PREC_DP), DIMENSION(2)  , INTENT(INOUT) :: cold  ! old correlation at atom 1,2
REAL(KIND=PREC_DP), DIMENSION(2)  , INTENT(INOUT) :: cnew  ! old correlation at atom 1,2
!
INTEGER, PARAMETER :: ISANY     = 0
INTEGER, PARAMETER :: ISWORST   = 1
INTEGER, PARAMETER :: ISPERFECT = 2
!
REAL(kind=PREC_DP) :: r1
INTEGER :: ngrand                         ! End   atom number for random search
INTEGER :: istart                         ! Start atom number for random search
INTEGER :: icounter                       ! End   atom number for random search
INTEGER :: inumber                        ! End   atom number for random search
INTEGER :: iinc                           ! Increment         for random search
INTEGER :: latom                          ! Current atom number
INTEGER :: ic                             ! Loop index correlations
!INTEGER :: ie                             ! Loop index energy type
!INTEGER :: jc                             ! Current correlation to improve
!INTEGER :: js                             ! Loop index atom type 
!INTEGER :: ks                             ! Loop index atom type 
!INTEGER :: i                              ! Dummy loop index
INTEGER :: ia                             ! Dummy loop index
INTEGER :: icent                          ! Dummy loop index
integer, dimension(0:2) :: isnei
!
LOGICAL :: valid_e
!
ngrand = 0
grand: DO
!
   iinc = 1
   inumber = cr_icc(1)
!
   ngrand = ngrand + 1
   CALL RANDOM_NUMBER(r1)                    ! 
IF(r1<1./2.) THEN                         ! Increment along x-axis
   IF(cr_icc(1) > 1) THEN
      CALL RANDOM_NUMBER(r1)
      iinc = NINT(SIGN(1.0D0, r1-0.5D0))*cr_ncatoms  ! Randomly choose increment -1 or +1 * atoms per unit cell
      inumber = cr_icc(1)
   ELSE
      CYCLE grand   
   ENDIF
ELSEIF(r1<2./2+1) THEN                      ! Increment along y-axis
   IF(cr_icc(2) > 1) THEN
      CALL RANDOM_NUMBER(r1)
      iinc = NINT(SIGN(1.0D0, r1-0.5D0))*cr_ncatoms*cr_icc(1)  ! Randomly choose increment -1 or +1 * atoms per unit cell
      inumber = cr_icc(2)
   ELSE
      CYCLE grand   
   ENDIF
!ELSE                                      ! Increment along z-axis
!   CYCLE grand   
ENDIF
!iinc   = NINT(SIGN(1.0, r1-0.5))          ! Randomly choose increment -1 or +1
!CALL RANDOM_NUMBER(r1)                    ! 
istart = INT(r1*cr_natoms) + 1            ! Initalize loop over all atoms at random start
icounter = 0
!iend   = istart + iinc*(cr_natoms-1)      ! End after full list of atoms
latom = istart
!write(*,*) ' LOOP       ', istart, inumber, iinc, ngrand
search_ini: DO                            ! Loop until we find an atom that could be improved
   isel(ichoose) = MOD(latom-1+cr_natoms, cr_natoms) + 1
!write(*,*) ' latom isel(1), wjs ', latom, isel(ichoose), wjks(ichoose)
   latom = latom + iinc
   icounter = icounter + 1
   IF(icounter==inumber) THEN
      EXIT search_ini   
   ENDIF
   IF(cr_iscat(1,isel(ichoose))==wjks(ichoose)) THEN          ! Found atom of proper type
      CALL chem_neighbour_multi(isel(ichoose) , ic, iatom, patom, tatom, natom, &
                                ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)
      ia = ichoose
      icent = 1
      IF(mmc_cor_energy (ic, MC_OCC)) THEN
         cold(ichoose)=mmc_occ_correl(isel, ia, ic, iatom, icent,     &
                                      natom, valid_e, MAX_ATOM_ENV_L, &
                                      MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
      ELSEIF(mmc_cor_energy (ic, MC_UNI)) THEN
         cold(ichoose)=mmc_uni_correl(isel, ia, ic, iatom, icent,     &
                                      natom, valid_e, MAX_ATOM_ENV_L, &
                                      MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
      ENDIF
      isnei(0) = iatom(2, 1)
      isnei(1) = iatom(0, 1)
      isnei(2) = iatom(1, 1)
!write(*,*) ' latom', latom - iinc
!write(*,*) ' ATOM ', isel(ichoose),          ' : ',          isnei(0:2)          
!write(*,*) ' TYPES', cr_iscat(1,isel(ichoose)),' : ', cr_iscat(1,isnei(0:2)          ), cold(ichoose)
!write(*,*) ' Targ, ach, diff ', mmc_target_corr(wic, wie, wjks(1), wjks(2)), &
!mmc_ach_corr   (wic, wie, wjks(1), wjks(2)), &
!mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2)) , ' T ? A '
      IF(mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2))<0) THEN
!write(*,*) ' Targ, ach, diff ', mmc_target_corr(wic, wie, wjks(1), wjks(2)), &
!mmc_ach_corr   (wic, wie, wjks(1), wjks(2)), &
!mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2)) , ' T < A '
         IF(cold(ichoose)>mmc_target_corr(wic, wie, wjks(1), wjks(2))) THEN
!write(*,*) '1 IMPROVEMENT POSSIBLE if cnew becomes smaller            T < corr'
            IF(imode==ISWORST .OR. IMODE==ISANY) THEN 
               EXIT search_ini
            ENDIF
         ELSEIF(cold(ichoose)<mmc_target_corr(wic, wie, wjks(1), wjks(2))) THEN
!write(*,*) '2 IMPROVEMENT POSSIBLE if cnew becomes even smaller       T > corr'
            IF(IMODE==ISANY) THEN 
               CALL RANDOM_NUMBER(r1)                    ! 
               IF(r1<mmc_g_rate) EXIT search_ini
            ENDIF
         ELSE
            IF(IMODE==ISPERFECT) THEN 
!write(*,*) '3 ALREADY PERFECT '
               EXIT search_ini
            ENDIF
         ENDIF
      ELSEIF(mmc_target_corr(wic, wie, wjks(1), wjks((2)))-mmc_ach_corr (wic, wie, wjks(1), wjks(2))>0) THEN
!write(*,*) ' Targ, ach, diff ', mmc_target_corr(wic, wie, wjks(1), wjks(2)), &
!mmc_ach_corr   (wic, wie, wjks(1), wjks(2)), &
!mmc_target_corr(wic, wie, wjks(1), wjks(2))-mmc_ach_corr (wic, wie, wjks(1), wjks(2)) , ' T > A '
         IF(cold(ichoose)>mmc_target_corr(wic, wie, wjks(1),wjks(2))) THEN
!write(*,*) '4 IMPROVEMENT POSSIBLE if cnew becomes even larger        T < corr'
         ELSEIF(cold(ichoose)<mmc_target_corr(wic, wie, wjks(1), wjks(2))) THEN
!write(*,*) '5 IMPROVEMENT POSSIBLE if cnew becomes larger             T > corr'
exit search_ini
         ELSE
!write(*,*) '6 ALREADY PERFECT '
         ENDIF
      ELSE
!write(*,*) '7 ALREADY PERFECT '
      ENDIF
!exit search_ini
   ENDIF
ENDDO search_ini
   IF(isel(ichoose) /=0) EXIT grand
   IF(ngrand>10   ) THEN
!write(*,*) ' NO MOVE FOUND ', iinc, istart, inumber
      ier_num=-2
      ier_typ = ER_MMC
      RETURN
   ENDIF
ENDDO grand
!if(isel(ichoose) == 0) then
!  write(*,*) ' WRONG MOVE B' , iinc, istart, inumber, ngrand
!  write(*,*) ' WRONG MOVE B' , isel, ' : ' , ichoose, latom
!endif
!
END SUBROUTINE mmc_select_emin
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
use precision_mod
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
LOGICAL                           , INTENT(in ) :: loop
!
INTEGER  :: i, j
REAL(kind=PREC_DP) :: z
REAL(kind=PREC_DP), DIMENSION(3) ::  v, u
REAL(kind=PREC_DP) :: r1
!
!write(*,*) ' MMC_SEL_ATOM ', MMC_MAX_ATOM_L, mmc_move
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
            laccept = mmc_allowed(cr_iscat(1,isel(1) ) ) .AND. &
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
         laccept = mmc_allowed(cr_iscat(1,isel(1) ) ) .AND.&
                   mmc_allowed(cr_iscat(1,isel(2) ) ) .AND.&
                   check_select_status(isel(1),.TRUE., cr_prop(isel(1) ), cr_sel_prop) .AND.&
                   check_select_status(isel(2),.TRUE., cr_prop(isel(2) ), cr_sel_prop) .AND.&
                   skalpro(u, v, cr_gten) < 0.0
      ELSEIF (mmc_move == MC_MOVE_SWCHEM) then
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
         laccept = cr_iscat (1,isel (1))  /= cr_iscat (1,isel (2) ) .AND.                  &
                 ( mmc_allowed (cr_iscat (1,isel (1) ) ) .AND.                            &
                   mmc_allowed (cr_iscat (1,isel (2) ) )      )    .AND.                  &
                   check_select_status (isel(1), .TRUE., cr_prop (isel (1) ), cr_sel_prop) .AND. &
                   check_select_status (isel(2), .TRUE., cr_prop (isel (2) ), cr_sel_prop) 
      ELSEIF (mmc_move == MC_MOVE_SWNEIG) then
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
         laccept = cr_iscat (1,isel (1))  /= cr_iscat (1,isel (2) ) .AND.                  &
                 ( mmc_allowed (cr_iscat (1,isel (1) ) ) .AND.                            &
                   mmc_allowed (cr_iscat (1,isel (2) ) )      )    .AND.                  &
                   check_select_status (isel(1), .TRUE., cr_prop (isel (1) ), cr_sel_prop) .AND. &
                   check_select_status (isel(2), .TRUE., cr_prop (isel (2) ), cr_sel_prop) .and. &
         any(mmc_pair(:,MC_COORDNUM, cr_iscat(1, isel(1)),:)==-1) .and.  &
         any(mmc_pair(:,MC_COORDNUM, cr_iscat(1, isel(2)),:)==-1)
      ENDIF 
!                                                                       
!-----      ----Check whether geometrical constrains apply              
!                                                                       
      IF (laccept) THEN 
         CALL check_geometry_multi (isel, natoms, laccept, MMC_MAX_ATOM) 
      ENDIF 
      IF(laccept) EXIT main
ENDDO main
!write(*,*) ' ISEL 1 ', isel(1), cr_iscat(1,isel(1)), any(mmc_pair(:,MC_COORDNUM, cr_iscat(1, isel(1)),:)==-1)
!write(*,*) ' ISEL 2 ', isel(2), cr_iscat(1,isel(2)), any(mmc_pair(:,MC_COORDNUM, cr_iscat(1, isel(2)),:)==-1)
!write(*,'(a, 2i7, a, 2i3, a, 2i7)') ' SELECTED PAIR ', isel(1:2), ' | ', cr_iscat(1,isel(1:2)), ' || ', isel(3:4)
!
END SUBROUTINE mmc_select_atoms
!
!*****7*****************************************************************
!
subroutine mmc_select_site(isite, nsel)
!-
!  Select a site from the allowed list in a random unit cell number
!+
!
use crystal_mod
use celltoindex_mod
!
use precision_mod
!
implicit none
!
integer                  , intent(in)  :: isite
integer                  , intent(out) :: nsel

real(kind=prec_SP) :: r1   ! random number
integer :: i         ! Dummy indices
integer, dimension(3) :: icell ! Unit cell number
!
do i=1, 3
   call random_number(r1)
   icell(i) = int(cr_icc(i)*r1) + 1
enddo
!
call celltoindex(icell, isite, nsel) 
!
end subroutine mmc_select_site
!
!*****7*****************************************************************
!
SUBROUTINE mmc_energies(isel, is, iz, natoms, iatom, patom, tatom, natom, ncent, &
                        rdi, rdj, valid_all, e_cur, CHEM_MAX_COR_L,       &
                        MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L, is_old)
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
INTEGER, DIMENSION(MMC_MAX_ATOM_L)                      , INTENT(inout) :: isel !(chem_max_atom) 
INTEGER, DIMENSION(2)                                   , INTENT(IN) :: is
INTEGER                                                 , INTENT(IN) :: natoms
INTEGER, DIMENSION(2, 3)                                , INTENT(IN) :: iz
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: iatom
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: patom
INTEGER, DIMENSION(     MMC_MAX_CENT_L)                 , INTENT(INOUT) :: natom
LOGICAL, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(INOUT) :: tatom
INTEGER                                                 , INTENT(INOUT) :: ncent
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR_L)           , INTENT(INOUT) :: rdi
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR_L)           , INTENT(INOUT) :: rdj
LOGICAL                                                 , INTENT(OUT) :: valid_all
REAL(kind=PREC_DP)   , DIMENSION(0:MC_N_ENERGY)                         , INTENT(INOUT) :: e_cur
logical                                                 , intent(in)    :: is_old
!
INTEGER :: i      ! Dummy loop variable
INTEGER :: ia     ! Dummy loop variable for modified atoms
INTEGER :: ic     ! Dummy loop variable for correlations
INTEGER :: icent  ! Dummy loop variable for centers
REAL(kind=PREC_DP)    :: delta
REAL(kind=PREC_DP)    :: x             ! random number
REAL(kind=PREC_DP), DIMENSION(3) :: v
REAL(kind=PREC_DP), DIMENSION(3) :: idir, jdir
LOGICAL            :: valid_e
!
!write(*,*) '=========='
!write(*,'(a, 2i7, 3x,2i3, a, i3)') ' MMC_ENERGIES  ', isel(1:2), cr_iscat(1,isel(1:2)), ' | ', natoms
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
!write(*,*) ' OLD ENERGIES', e_cur(MC_DISP), isel(1:2), ' SCAT: ',cr_iscat(1,isel(1)), cr_iscat(1,isel(2))
   loop_corr: DO ic = 1, chem_ncor 
      CALL chem_neighbour_multi(isel(ia), ic, iatom, patom, tatom, natom, &
                                ncent, MAX_ATOM_ENV_L, MMC_MAX_CENT_L)                                                
!write(*,*) ' central ', isel(ia), cr_iscat(1,isel(ia)), ic, 'NATOM : ', natom(1)
!write(*,*) ' neigb   ', iatom(0:natom(1),1), ' : ',cr_iscat(1,iatom(0:natom(1),1)), ' ::', ncent
!write(*,*) ' type is ', tatom(0:natom(1),1)

      loop_cent: DO icent = 1, ncent 
         IF (natom (icent)  > 0) THEN 
!write(*,*) ' AT ATOM ', ia, isel(ia), ic, iatom(0:natom(1),icent), mmc_cor_energy (ic, MC_OCC), mmc_move == MC_MOVE_SWCHEM
!read(*,*) i
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
!                 IF (cr_iscat (1,isel (1) )  /= cr_iscat (1,isel (2) ) ) THEN
                     e_cur (MC_OCC) = e_cur (MC_OCC) + mmc_energy_occ ( &
                     isel, ia, ic, iatom, tatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L) 
                     valid_all = valid_all.OR.valid_e 
!                 ENDIF 
               ENDIF 
            ENDIF 
!                                                                       
!     ------- Unidirectional Occupation Correlation, only if both atoms are different  
!     ------- and the move is switch chemistry                          
!                                                                       
!write(*,*)
!write(*,*) ' TEST FOR UNI ', mmc_cor_energy (ic, MC_UNI), mmc_move == MC_MOVE_SWCHEM
            IF (mmc_cor_energy (ic, MC_UNI) ) THEN 
               IF (mmc_move == MC_MOVE_SWCHEM) THEN 
!                 IF (cr_iscat (1,isel (1) )  /= cr_iscat (1,isel (2) ) ) THEN                                                  
                     e_cur (MC_UNI) = e_cur (MC_UNI) + mmc_energy_uni ( &
                     isel, ia, ic, iatom, tatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L) 
                     valid_all = valid_all.OR.valid_e 
!                 ENDIF 
               ENDIF 
            ENDIF 
!                                                                       
!     ------- Occupation Correlation, for group wise correlations
!     ------- and the move is switch chemistry                          
            IF(mmc_cor_energy(ic, MC_GROUP) ) THEN 
!                                                                       
               IF(mmc_move == MC_MOVE_SWCHEM) THEN 
!                 IF (cr_iscat (1,isel (1) )  /= cr_iscat (1,isel (2) ) ) THEN
                     e_cur (MC_GROUP) = e_cur (MC_GROUP) + mmc_energy_group ( &
                     isel, ia, ic, iatom, tatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L) 
                     valid_all = valid_all.OR.valid_e 
!                 ENDIF 
               ENDIF 
            ENDIF 
!
!     ------- Coordination number
!
            IF(mmc_cor_energy(ic, MC_COORDNUM) ) THEN 
               e_cur(MC_COORDNUM) = e_cur(MC_COORDNUM) + mmc_energy_cn ( &
               isel, ia, ic, iatom, icent, natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)       
               valid_all = valid_all.OR.valid_e 
               call random_number(x)
               if(is_old) then
               isel(ia+2) = iatom(int(natom(icent)*x) + 1, icent)
!write(*,'(a, i3, i6, i3, i3, a, 7i7,a, i7)') ' ENERGIES    ', ia, isel(ia), icent, natom(icent), ' | ', &
! iatom(0:6, icent),  ' > ',isel(ia+2)
               endif
            ENDIF 
!                                                                       
!     ------- Displacement correlation                                  
!                                                                       
            IF(mmc_cor_energy(ic, MC_DISP) ) THEN 
               idir = 0.0_PREC_DP
               DO i = 1, 3 
                 v(i) = cr_pos(i, isel(ia) ) - chem_ave_pos(i, is(ia)) &
                        -REAL(iz(ia, i) - 1, kind=PREC_DP) - cr_dim0(i, 1)                   
!              u(i) = v(i)
!           v (i) = v (i) - disp (i, 0, ia) 
               ENDDO 
               IF(chem_ldall (ic) ) THEN 
                  DO i = 1, 3 
                     jdir(i) = v(i) 
                  ENDDO 
                  rdj(ic) = skalpro (jdir, jdir, cr_gten) 
                  IF(rdj(ic)  > 0.0) THEN 
                     rdj(ic) = sqrt(rdj(ic) ) 
                  ELSE 
                     rdj(ic) = 1.0 
                  ENDIF 
                  delta = 1.0 
               ELSE 
                  DO i = 1, 3 
                     idir(i) = chem_dir(i, 1, ic) 
                     jdir(i) = chem_dir(i, 2, ic) 
                  ENDDO 
                  delta = skalpro(v, idir, cr_gten) / rdi(ic) 
               ENDIF 
!do i=1, 3
!write(*,'(11f10.4)') cr_pos(i, isel(ia) ), chem_ave_pos(i, is(ia)),  &
!              REAL(iz(ia, i) - 1, kind=PREC_DP), cr_dim0(i, 1), v(i) , &
!              idir, jdir
!enddo
!write(*,*) ' DELTA ', delta, chem_ldall(ic)
               e_cur(MC_DISP) = e_cur(MC_DISP) + mmc_energy_dis (  &
               isel, ia, ic, iatom, icent, natom, jdir, delta,&
               rdj, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                          
               valid_all = valid_all.OR.valid_e 
               valid_all = .TRUE. 
!write(*,*) ' CURR  ', e_cur(MC_DISP)
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
USE atom_env_mod
USE chem_mod
USE chem_neig_multi_mod
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
INTEGER                                            , INTENT(IN)    :: MMC_MAX_ATOM_L
INTEGER           , DIMENSION(MMC_MAX_ATOM_L)      , INTENT(inout) :: isel !(chem_max_atom) 
REAL(kind=PREC_DP), DIMENSION(3)                   , INTENT(OUT)   :: posz
REAL(kind=PREC_DP), DIMENSION(3)                   , INTENT(OUT)   :: posz2
REAL(KIND=PREC_DP), DIMENSION(3,0:MMC_MAX_ATOM_L,2), INTENT(INOUT) :: disp
!
INTEGER :: i, j, k              ! Dummy loop indices
INTEGER :: iscat
INTEGER, DIMENSION(3) :: iz1
INTEGER, DIMENSION(3) :: iz2
INTEGER, DIMENSION(2) :: is
REAL(kind=PREC_DP) :: disp1, disp2 
REAL(kind=PREC_DP) :: rrrr
!
integer :: ic, jc
INTEGER :: ncent
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: patom ! (3, 0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(  :,:), ALLOCATABLE :: iatom ! (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
LOGICAL           , DIMENSION(  :,:), ALLOCATABLE :: tatom ! (0:MAX_ATOM_ENV, MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: natom ! ( MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: pair1 ! ( MMC_MAX_CENT) 
INTEGER           , DIMENSION(    :), ALLOCATABLE :: pair2 ! ( MMC_MAX_CENT) 
integer :: npair1
integer :: npair2
REAL(KIND=PREC_DP) :: r1
integer :: istart1
integer :: istart2
integer, dimension(:,:), allocatable :: loc_site
integer :: nloc
integer :: nsel
!
!
IF (mmc_move == MC_MOVE_DISP) THEN 
!                                                                       
!     ------Displace atoms from current positions                       
!                                                                       
!                                                                       
!     ------Modify the central atom                                     
!                                                                       
   IF(mo_maxmove(4, cr_iscat(1,isel(1)))==0.0) THEN
      DO i = 1, 3 
         disp(i, 0, 1)      = gasdev(DBLE(mo_maxmove(i, cr_iscat(1,isel(1)) )))
         posz(i)            = cr_pos(i,isel(1)) 
         cr_pos(i, isel(1)) = cr_pos(i, isel(1)) + disp(i, 0, 1) 
      ENDDO 
   ELSE
!     -- Move along a vector direction
      rrrr = gasdev(DBLE(mo_maxmove(4, cr_iscat(1,isel(1)) )))
      DO i = 1, 3 
         disp(i, 0, 1) = rrrr* (mo_maxmove(i, cr_iscat(1,isel(1)) ))
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
   iscat             = cr_iscat(1,isel(2)) 
   cr_iscat(1,isel(2)) = cr_iscat(1,isel(1)) 
   cr_iscat(1,isel(1)) = iscat 
   iscat            = cr_prop(isel(2)) 
   cr_prop(isel(2)) = cr_prop(isel(1)) 
   cr_prop(isel(1)) = iscat 
ELSEIF (mmc_move == MC_MOVE_SWNEIG) THEN 
!                                                                       
!     ------Switch the Chemistry of two selected neighbor atoms                  
!                                                                       
!  allocate space for neighbors
   ALLOCATE(patom(3, 0:MAX_ATOM_ENV, MMC_MAX_CENT))
   ALLOCATE(iatom(0:MAX_ATOM_ENV, MMC_MAX_CENT))
   ALLOCATE(tatom(0:MAX_ATOM_ENV, MMC_MAX_CENT))
   ALLOCATE(natom(MMC_MAX_CENT))
!======================================
!
   j = 0
   nloc = 0
   do i=1, ubound(mmc_lsite,1)
      if(mmc_lsite(i)) j=j+1         ! Count number of additional sites
   enddo
   if(j>0)  then
      nloc = j                       !Total number of additional sites
      allocate(loc_site(1:nloc,2))   ! Room for additional sites
      loc_site = 0
      do k =1, 2                     ! Choose additional for isel(1) and isel(2)
         j = 0
         do i=1, ubound(mmc_lsite,1)
            if(mmc_lsite(i)) then        ! This site is allowed, pick an atom
               call mmc_select_site(i, nsel)
               j = j+1
               loc_site(j,k) = nsel
            endif
         enddo
      enddo
   else
      nloc = 0                      ! No aditional sites requested
      allocate(loc_site(1:1,2))
      loc_site = 0
   endif
!======================================
!write(*,*) ' SELECT NEIGHBORS ***********************'
!write(*,*) ' MMC_LSITE ', ubound(mmc_lsite,1), mmc_lsite
!write(*,*) 'LOC_SITE ', loc_site
   isel(3:4) = 0                          ! Default to no atoms switched
   j = isel(1) 
   jc = 1   !WORK 
   ic = 0   !WORK 
   loop_cn_ncor1:do jc=1, CHEM_MAX_COR    ! Find neighbors for ISEL(1)
      if(any(mmc_pair(jc,MC_COORDNUM, cr_iscat(1, isel(1)), :)==-1)) then
         ic = jc               ! This target used Coordination number
         exit loop_cn_ncor1
      endif
   enddo loop_cn_ncor1
!   write(*,*) ' ISEL 1 ', isel(1), cr_iscat(1, isel(1)), ic
   
   if(ic> 0) then                         ! Found a target
      CALL chem_neighbour_multi (j, ic, iatom, patom, tatom, natom, ncent, MAX_ATOM_ENV, MMC_MAX_CENT)
   else
      ncent = 1
      natom(ncent) = 0
   endif
   cond_pair1:if(natom(ncent)+nloc > 0) then         ! Found any atoms to switch
      allocate(pair1(0:natom(ncent)+nloc-1))   ! All neighbors plus additional sites
!write(*,*) ' PAIR 1 ', lbound(pair1), ubound(pair1)
      pair1(0:natom(ncent)-1) = iatom(1:natom(ncent), ncent)
      if(nloc>0) pair1(natom(ncent):natom(ncent)+nloc-1) = loc_site(:,1)
      npair1 = natom(ncent)+nloc
!write(*,'(a,i7,a,9i7)') ' Neighbors isel1     ', iatom(0,1) , ' | ', pair1!(0:npair1-1)
!
      j = isel(2) 
      jc = 1   !WORK 
      ic = 0   !WORK 
      loop_cn_ncor2:do jc=1, CHEM_MAX_COR    ! Find neighbors for ISEL(1)
         if(any(mmc_pair(jc,MC_COORDNUM, cr_iscat(1, isel(2)), :)==-1)) then
            ic = jc
            exit loop_cn_ncor2
         endif
      enddo loop_cn_ncor2
!  write(*,*) ' ISEL 1 ', isel(2), cr_iscat(1, isel(2)), ic
      if(ic> 0) then                         ! Found a target
         CALL chem_neighbour_multi (j, ic, iatom, patom, tatom, natom, ncent, MAX_ATOM_ENV, MMC_MAX_CENT)
      else
         ncent = 1
         natom(ncent) = 0
      endif
      cond_pair2:if(natom(ncent)+nloc > 0) then         ! Found any atoms to switch
         allocate(pair2(0:natom(ncent)+nloc-1))
         pair2(0:natom(ncent)-1) = iatom(1:natom(ncent), ncent)
         if(nloc>0) pair2(natom(ncent):natom(ncent)+nloc-1) = loc_site(:,2)
         npair2 = natom(ncent)+nloc
!write(*,'(a,i7,a,9i7)') ' Neighbors isel2     ', iatom(0,1) , ' | ', pair2!(0:npair2-1)
!
!write(*,'(a,4i7)') '   modify PRIOR MOVE ', isel(1:4)
!write(*,'(a,4i7)') '          atom types ', cr_iscat(1,isel(1:4))
   
         call random_number(r1)
         istart1 = int(6*r1)
         call random_number(r1)
         istart2 = int(6*r1)
         loop_search: do i=0,npair1+nloc-1     ! Assign a random pair with different chemistry
            ic = mod(istart1 + i, npair1)
            do j=0,npair2+nloc-1
               jc = mod(istart2 + j, npair2)
               if(cr_iscat(1, pair1(ic)) /= cr_iscat(1, pair2(jc))) exit loop_search
            enddo
         enddo loop_search
         isel(3) = pair1(ic)
         isel(4) = pair2(jc)
!
         iscat             = cr_iscat(1,isel(4)) 
         cr_iscat(1,isel(4)) = cr_iscat(1,isel(3)) 
         cr_iscat(1,isel(3)) = iscat 
         iscat            = cr_prop(isel(4)) 
         cr_prop(isel(4)) = cr_prop(isel(3)) 
         cr_prop(isel(3)) = iscat 
         deallocate(pair2)
      endif cond_pair2
      deallocate(pair1)
   endif cond_pair1
! write(*,'(a,4i7)') '   modify POST  MOVE ', isel(1:4)
! write(*,'(a,4i7)') '          atom types ', cr_iscat(1,isel(1:4))
   deallocate(patom)
   deallocate(iatom)
   deallocate(tatom)
   deallocate(natom)
!  deallocate(lsite)
   deallocate(loc_site)
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
REAL(kind=PREC_DP)   , DIMENSION(3)             , INTENT(IN) :: posz
REAL(kind=PREC_DP)   , DIMENSION(3)             , INTENT(IN) :: posz2
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
   iscat             = cr_iscat(1,isel(2)) 
   cr_iscat(1,isel(2)) = cr_iscat(1,isel(1)) 
   cr_iscat(1,isel(1)) = iscat 
   iscat            = cr_prop(isel(2)) 
   cr_prop(isel(2)) = cr_prop(isel(1)) 
   cr_prop(isel(1)) = iscat 
ELSEIF (mmc_move == MC_MOVE_SWNEIG) THEN 
!                                                                       
!     ------Switch the Chemistry of two selected atoms                  
!                                                                       
!write(*,'(a,4i7)') ' UNMODIFY PRIOR MOVE ', isel(1:4)
!write(*,'(a,4i7)') '          atom types ', cr_iscat(1,isel(1:4))
   if(isel(3)/=0 .and. isel(4)/=0) then
      iscat             = cr_iscat(1,isel(4)) 
      cr_iscat(1,isel(4)) = cr_iscat(1,isel(3)) 
      cr_iscat(1,isel(3)) = iscat 
      iscat            = cr_prop(isel(4)) 
      cr_prop(isel(4)) = cr_prop(isel(3)) 
      cr_prop(isel(3)) = iscat 
!write(*,'(a,4i7)') ' UNMODIFY POST  MOVE ', isel(1:4)
!write(*,'(a,4i7)') '          atom types ', cr_iscat(1,isel(1:4))
   endif
!                                                                       
!--End of Modification of atoms according to different moves
!                                                                       
ENDIF 
!
END SUBROUTINE mmc_unmodify
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_occ_correl(isel, ia, ic, iatom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
!+                                                                      
!     Calculates the local correlation for chemical disorder                       
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE modify_func_mod
USE param_mod 
use precision_mod
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
mmc_occ_correl = 0.0 
ncalc = 0 
valid_e = .FALSE. 
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
   IF(natom(icent)  /= 0) THEN 
      IF(isel(ia) == iatom(0, icent)) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         in_a = 1 
         in_e = natom (icent) 
         is = cr_iscat (1,iatom (0, icent) ) 
         DO ind = in_a, in_e 
               js = cr_iscat (1,iatom (ind, icent) ) 
            IF(mmc_pair(ic, MC_OCC,is,js)/=0) THEN
            IF(check_select_status(iatom (ind, icent),  &
                                     .TRUE., cr_prop (iatom (ind,  &
                                     icent) ), cr_sel_prop) ) THEN                         
               ival1 = 0 
               js = cr_iscat (1, iatom (ind, icent) ) 
               ival1 = SIGN(1,mmc_pair(ic, MC_OCC,is,js))
               mmc_occ_correl = mmc_occ_correl + ival1
            ENDIF 
            ELSE
               mmc_occ_correl = mmc_occ_correl - 1
            ENDIF 
         ENDDO 
      ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
         in_a = 0 
         in_e = 0 
         is = cr_iscat (1,isel (ia) ) 
         ind = 0
            js = cr_iscat (1,iatom (ind, icent) ) 
            IF(mmc_pair(ic, MC_OCC,is,js)/=0) THEN
         IF (check_select_status (iatom (ind, icent),  &
                                  .TRUE., cr_prop (iatom (ind,  &
                                  icent) ), cr_sel_prop) ) THEN                         
            ival1 = 0 
            js = cr_iscat (1,iatom (ind, icent) ) 
            ival1 = SIGN(1,mmc_pair(ic, MC_OCC,is,js))
            mmc_occ_correl = mmc_occ_correl + ival1
         ENDIF 
            ELSE
               mmc_occ_correl = mmc_occ_correl - 1
         ENDIF 
      ENDIF 
   ENDIF 
ENDIF 
valid_e = .TRUE. 
!                                                                       
END FUNCTION mmc_occ_correl                   
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_uni_correl(isel, ia, ic, iatom, icent,  &
      natom, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L)                                                   
!+                                                                      
!     Calculates the local correlation for chemical disorder                       
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE mc_mod 
USE mmc_mod 
USE modify_mod
USE modify_func_mod
USE param_mod 
use precision_mod
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
mmc_uni_correl = 0.0 
ncalc = 0 
valid_e = .FALSE. 
is = cr_iscat(1,iatom(0,icent))
js = cr_iscat(1,iatom(1,icent))
!write(*,*) ' IN MMC_CORREL_UNI         ', iatom(0:natom(icent),icent), is, js
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
   IF(natom(icent)  /= 0) THEN 
      IF(isel(ia) == iatom(0, icent)) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         in_a = 1 
         in_e = natom (icent) 
         is = cr_iscat (1,iatom (0, icent) ) 
         IF(ABS(mmc_pair(ic, MC_UNI,is,is))==1) THEN               ! "is" is of "left==starting" atom type
!write(*,*) ' Selected is central atom type ic, ia, is            ', &
! ic, ia, isel(ia), cr_iscat(1,isel(ia)) 
         DO ind = in_a, in_e 
            js = cr_iscat (1,iatom (ind, icent) ) 
            IF(ABS(mmc_pair(ic, MC_UNI,is,js))==1) THEN            ! "is,is" or "is,js"  with js=="ending"
            IF(check_select_status(iatom (ind, icent),  &
                                     .TRUE., cr_prop (iatom (ind,  &
                                     icent) ), cr_sel_prop) ) THEN                         
               ival1 = SIGN(1,mmc_pair(ic, MC_UNI,is,js))
               mmc_uni_correl = mmc_uni_correl + ival1
!write(*,*) ' Pair is of type  PAIR is, js, ival1,CORR, iatoms    ', &
! is, js, ival1, mmc_uni_correl, iatom(0,icent), iatom(ind,icent)
            ENDIF 
            ELSE
               mmc_uni_correl = mmc_uni_correl + 1
!write(*,*) ' PAIR is of type  OTHER                              ', &
!  is, js, 1, mmc_uni_correl
            ENDIF 
         ENDDO 
         ELSE
!write(*,*) ' Selected is NOT central atom type                   ', &
! ic, ia, isel(ia), cr_iscat(1,isel(ia))
!write(*,*) ' PAIR is of type  OTHER                              ', 0.0 
         ENDIF
      ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
!write(*,*) ' Selected is neigh   atom type ', ic, ia, cr_iscat(1,ia)
         in_a = 0 
         in_e = 0 
         is = cr_iscat (1,isel (ia) ) 
         ind = 0
            js = cr_iscat (1,iatom (ind, icent) ) 
            IF(mmc_pair(ic, MC_UNI,is,js)/=0) THEN
         IF (check_select_status (iatom (ind, icent),  &
                                  .TRUE., cr_prop (iatom (ind,  &
                                  icent) ), cr_sel_prop) ) THEN                         
            ival1 = 0 
            js = cr_iscat (1,iatom (ind, icent) ) 
            ival1 = SIGN(1,mmc_pair(ic, MC_UNI,is,js))
            mmc_uni_correl = mmc_uni_correl + ival1
         ENDIF 
            ELSE
               mmc_uni_correl = mmc_uni_correl + 1
         ENDIF 
      ENDIF 
   ENDIF 
ENDIF 
valid_e = .TRUE. 
!                                                                       
END FUNCTION mmc_uni_correl                   
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_energy_occ(isel, ia, ic, iatom, tatom, icent,  &
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
use precision_mod
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
LOGICAL, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: tatom
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
   IF(natom(icent) /= 0) THEN              ! We do have neighbors
      in_a = 1 
      in_e = natom (icent) 
      DO ind = in_a, in_e 
         IF(tatom(ind, icent)) THEN  ! Selected atom is central
            is = cr_iscat (1,iatom (0, icent) )
            js = cr_iscat (1,iatom (ind, icent) ) 
!              ival1 = 0 
!           ival1 = sign(1,mmc_pair(ic, MC_OCC,is,js))
         ELSE                        ! Selected atom is a neighbor
            js = cr_iscat (1, iatom (0, icent) )
            is = cr_iscat (1, iatom (ind, icent) ) 
!           ival1 = mmc_pair(ic, MC_OCC,is,js)
         ENDIF
!write(*,*) 'OCC ', tatom(ind, icent), is, js, mmc_depth(ic,MC_OCC, is, js)
         IF(mmc_pair(ic, MC_OCC,is,js)/=0) THEN
            IF (check_select_status (iatom (ind, icent),  &
                                     .TRUE., cr_prop (iatom (ind,  &
                                     icent) ), cr_sel_prop) ) THEN                         
!              ival1 = 0 
               ival1 = sign(1,mmc_pair(ic, MC_OCC,is,js))
               mmc_energy_occ = mmc_energy_occ +                  &
                                mmc_depth (ic,MC_OCC, is, js) * ival1
            ENDIF 
         ELSE
            mmc_energy_occ = mmc_energy_occ - 1
         ENDIF 
!IF(tatom(ind, icent)) THEN  ! Selected atom is central
!write(*,'(2(a,i6, i3), f7.2)') ' MMC_OCC Central ', iatom (0, icent), is, ' : ', iatom (ind, icent), js, &
!mmc_depth (ic,MC_OCC, is, js) * ival1
!else
!write(*,'(2(a,i6, i3), f7.2)') ' MMC_OCC NEIG    ', iatom (0, icent), js, ' : ', iatom (ind, icent), is, &
!mmc_depth (ic,MC_OCC, is, js) * ival1
!endif
      ENDDO
   ENDIF 
ENDIF 
!!      IF (isel (ia)  == iatom (0, icent) ) THEN 
!!!                                                                       
!!!     ----The selected atom is the central atom, check all atoms        
!!!                                                                       
!!         in_a = 1 
!!         in_e = natom (icent) 
!!         is = cr_iscat (1,iatom (0, icent) ) 
!!         DO ind = in_a, in_e 
!!            js = cr_iscat (1,iatom (ind, icent) ) 
!!            IF(mmc_pair(ic, MC_OCC,is,js)/=0) THEN
!!            IF (check_select_status (iatom (ind, icent),  &
!!                                     .TRUE., cr_prop (iatom (ind,  &
!!                                     icent) ), cr_sel_prop) ) THEN                         
!!               ival1 = 0 
!!               ival1 = sign(1,mmc_pair(ic, MC_OCC,is,js))
!!               mmc_energy_occ = mmc_energy_occ +                  &
!!                                mmc_depth (ic,MC_OCC, is, js) * ival1
!!            ENDIF 
!!            ELSE
!!               mmc_energy_occ = mmc_energy_occ - 1
!!            ENDIF 
!!         ENDDO 
!!      ELSE 
!!!                                                                       
!!!     The selected atom is a neighbour, use this atom only              
!!!                                                                       
!!         in_a = 0 
!!         in_e = 0 
!!         is = cr_iscat (1,isel (ia) ) 
!!         ind = 0
!!            js = cr_iscat (1,iatom (ind, icent) ) 
!!            IF(mmc_pair(ic, MC_OCC,is,js)/=0) THEN
!!         IF (check_select_status (iatom (ind, icent),  &
!!                                  .TRUE., cr_prop (iatom (ind,  &
!!                                  icent) ), cr_sel_prop) ) THEN                         
!!            ival1 = 0 
!!            ival1 = mmc_pair(ic, MC_OCC,is,js)
!!            mmc_energy_occ = mmc_energy_occ +                  &
!!                             mmc_depth (ic,MC_OCC, is, js) * ival1
!!         ENDIF 
!!            ELSE
!!               mmc_energy_occ = mmc_energy_occ - 1
!!            ENDIF 
!!      ENDIF 
!!   ENDIF 
!!ENDIF 
valid_e = .TRUE. 
!                                                                       
END FUNCTION mmc_energy_occ                   
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_energy_uni(isel, ia, ic, iatom, tatom, icent,  &
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
use precision_mod
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
LOGICAL, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: tatom
INTEGER                                              , INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L)                   , INTENT(IN) :: natom
LOGICAL                                              , INTENT(OUT) :: valid_e 
!                                                                       
!                                                                       
INTEGER :: is, js, ind
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
INTEGER :: ival1
!                                                                       
mmc_energy_uni = 0.0 
ncalc   = 0 
valid_e = .FALSE. 
!                                                                       
!QIF (chem_ctyp(ic) == CHEM_VEC    .OR. &
!Q    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
!Q    chem_ctyp(ic) == CHEM_RANGE  .OR. &
!Q    chem_ctyp(ic) == CHEM_DIST   .OR. &
!Q    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!Q!                                                                       
!Q   IF(natom(icent) /= 0) THEN              ! We do have neighbors
!Q      in_a = 1 
!Q      in_e = natom (icent) 
!Q      DO ind = in_a, in_e 
!Q         IF(tatom(ind, icent)) THEN  ! Selected atom is central
!Q            is = cr_iscat (1,iatom (0, icent) )
!Q            js = cr_iscat (1,iatom (ind, icent) ) 
!Q!              ival1 = 0 
!Q!           ival1 = sign(1,mmc_pair(ic, MC_UNI,is,js))
!Q         ELSE                        ! Selected atom is a neighbor
!Q            js = cr_iscat (1,iatom (0, icent) )
!Q            is = cr_iscat (1,iatom (ind, icent) ) 
!Q!           ival1 = mmc_pair(ic, MC_UNI,is,js)
!Q         ENDIF
!Q         IF(mmc_pair(ic, MC_UNI,is,js)/=0) THEN
!Q            IF (check_select_status (iatom (ind, icent),  &
!Q                                     .TRUE., cr_prop (iatom (ind,  &
!Q                                     icent) ), cr_sel_prop) ) THEN                         
!Q!              ival1 = 0 
!Q               ival1 = sign(1,mmc_pair(ic, MC_UNI,is,js))
!Q               mmc_energy_uni = mmc_energy_uni +                  &
!Q                                mmc_depth (ic,MC_UNI, is, js) * ival1
!Q            ENDIF 
!Q         ELSE
!Q            mmc_energy_uni = mmc_energy_uni - 1
!Q         ENDIF 
!QIF(tatom(ind, icent)) THEN  ! Selected atom is central
!!QQwrite(*,'(2(a,i6, i3), f7.2)') ' MMC_UNI Central ', iatom (0, icent), is, ' : ', iatom (ind, icent), js, &
!Qmmc_depth (ic,MC_UNI, is, js) * ival1
!Qelse
!Qwrite(*,'(2(a,i6, i3), f7.2)') ' MMC_UNI NEIG    ', iatom (0, icent), js, ' : ', iatom (ind, icent), is, &
!Qmmc_depth (ic,MC_UNI, is, js) * ival1
!Qendif
!Q      ENDDO
!Q   ENDIF 
!QENDIF 
!QRETURN
!
ind = 1
is = cr_iscat(1,iatom(0,   icent) )
js = cr_iscat(1,iatom(ind, icent) )
!write(*,'(a, 6i5,a,3i5)') ' IN ENERGY UNI   ', ic,iatom (0, icent), iatom (ind, icent), &
!is,js, mmc_pair(ic, MC_UNI, is,js), ' : ',isel(ia), ia, natom(icent)
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
!                                                                       
   IF(natom(icent) /= 0) THEN              ! We do have neighbors
      in_a = 1 
      in_e = natom (icent) 
      DO ind = in_a, in_e 
         IF(tatom(ind, icent)) THEN  ! Selected atom is central
            is = cr_iscat (1,iatom (0, icent) )
            js = cr_iscat (1,iatom (ind, icent) ) 
!write(*,'(a,4i5)') 'CENTRAL        ', is, js,    (mmc_pair(ic, MC_UNI,is,is)),  &
!                                                 (mmc_pair(ic, MC_UNI,is,js))
            IF(is/=js                            ) THEN               ! Different atom types
            IF(ABS(mmc_pair(ic, MC_UNI,is,is))==1 .AND.    &          ! "is" is of "left==starting" atom type
               ABS(mmc_pair(ic, MC_UNI,is,js))==1) THEN               ! "is,is" or "is,js"  with js=="starting"
               IF(check_select_status(iatom(ind, icent),  .TRUE. ,  &
                                      cr_prop(iatom(ind, icent)) ,  &
                                                     cr_sel_prop) ) THEN                         
                  ival1 = sign(1,mmc_pair(ic, MC_UNI,is,js))
                  mmc_energy_uni = mmc_energy_uni +                  &
                                   mmc_depth (ic,MC_UNI, is, js) * ival1
!write(*,*) ' C: 1 Different atoms, ',(is),' ==> ',(js), &
!'  deltaE ', mmc_depth (ic,MC_UNI, is, js) * ival1, &
!iatom (0  , icent), &
!iatom (ind, icent), &
!' '
               ENDIF
            ELSEIF(ABS(mmc_pair(ic, MC_UNI,is,is))==2 .AND.    &          ! "is" is of "right==ending" atom type
                   ABS(mmc_pair(ic, MC_UNI,is,js))==2) THEN               ! "is,is" or "is,js"  with is=="ending"
               IF(check_select_status(iatom(ind, icent),  .TRUE. ,  &
                                      cr_prop(iatom(ind, icent)) ,  &
                                                     cr_sel_prop) ) THEN                         
                  mmc_energy_uni = mmc_energy_uni +                  &
                               ABS(mmc_depth (ic,MC_UNI, is, js))
!write(*,*) ' C: 2 Different atoms, ',(is),' ==> ',(js), &
!'  deltaE ', ABS(mmc_depth (ic,MC_UNI, is, js)), &
!iatom (0  , icent), &
!iatom (ind, icent), &
!' '
               ENDIF
            ENDIF
            ELSE                                                       ! Identical atom types
                  ival1 =  SIGN(1,mmc_pair(ic, MC_UNI,is,js))
                  mmc_energy_uni = mmc_energy_uni +                  &
                                  (mmc_depth (ic,MC_UNI, is, js)) * ival1
!write(*,*) ' C: 3 Equal     atoms, ',(is),' ==> ',(is), &
!'  deltaE ', (mmc_depth (ic,MC_UNI, is, js)) * ival1,&
!iatom (0  , icent), &
!iatom (ind, icent), &
!' '
            ENDIF
         ELSE                        ! Selected atom is a neighbor
            js = cr_iscat (1,iatom (0, icent) )                         ! Scattering type at end of vector
            is = cr_iscat (1,iatom (ind, icent) )                       ! Scattering type at start of vector
!write(*,'(a,4i5)') 'Neighb  ', is, js, ABS(mmc_pair(ic, MC_UNI,js,js)),  &
!                                       ABS(mmc_pair(ic, MC_UNI,js,is))
            IF(is/=js                            ) THEN               ! Different atom types
            IF(ABS(mmc_pair(ic, MC_UNI,is,is))==2 .AND.    &          ! "js" is of "right==ending" atom type
               ABS(mmc_pair(ic, MC_UNI,is,js))==2) THEN               ! "is,is" or "is,js"  with is=="ending"
!              ival1 = 0 
               IF(check_select_status(iatom(ind, icent),  .TRUE. ,  &
                                      cr_prop(iatom(ind, icent)) ,  &
                                                     cr_sel_prop) ) THEN                         
                  ival1 = sign(1,mmc_pair(ic, MC_UNI,is,js))
                  mmc_energy_uni = mmc_energy_uni +                  &
                               ABS(mmc_depth (ic,MC_UNI, is, js))
!write(*,*) ' N: 4 Different atoms, ',(is),' ==> ',(js), &
!'  deltaE ', ABS(mmc_depth (ic,MC_UNI, is, js)), &
!iatom (ind, icent), &
!iatom (0  , icent), &
!' '
               ENDIF
            ELSEIF(ABS(mmc_pair(ic, MC_UNI,is,is))==1 .AND.    &          ! "is" is of "left==starting" atom type
                   ABS(mmc_pair(ic, MC_UNI,is,js))==1) THEN               ! "is,is" or "is,js"  with is=="starting"
               IF(check_select_status(iatom(ind, icent),  .TRUE. ,  &
                                      cr_prop(iatom(ind, icent)) ,  &
                                                     cr_sel_prop) ) THEN                         
                  ival1 = sign(1,mmc_pair(ic, MC_UNI,is,js))
                  mmc_energy_uni = mmc_energy_uni +                  &
                                  (mmc_depth (ic,MC_UNI, is, js))*ival1
!write(*,*) ' N: 5 Different atoms, ',(is),' ==> ',(js), &
!'  deltaE ',    (mmc_depth (ic,MC_UNI, is, js))*ival1, &
!iatom (ind, icent), &
!iatom (0  , icent), &
!' '
               ENDIF
            ENDIF
            ELSE                                                       ! Identical atom types
                  ival1 =  SIGN(1,mmc_pair(ic, MC_UNI,is,js))
                  mmc_energy_uni = mmc_energy_uni +                  &
                                  (mmc_depth (ic,MC_UNI, is, js)) * ival1
!write(*,*) ' N: 6 Equal     atoms, ',(is),' ==> ',(is), &
!'  deltaE ', (mmc_depth (ic,MC_UNI, is, js)) * ival1,&
!iatom (ind, icent), &
!iatom (0  , icent), &
!' '
            ENDIF
         ENDIF
!        IF(mmc_pair(ic, MC_OCC,is,js)/=0) THEN
!           IF (check_select_status (iatom (ind, icent),  &
!                                    .TRUE., cr_prop (iatom (ind,  &
!                                    icent) ), cr_sel_prop) ) THEN                         
!              ival1 = 0 
!              ival1 = sign(1,mmc_pair(ic, MC_OCC,is,js))
!              mmc_energy_occ = mmc_energy_occ +                  &
!                               mmc_depth (ic,MC_OCC, is, js) * ival1
!           ENDIF 
!        ELSE
!           mmc_energy_occ = mmc_energy_occ - 1
!        ENDIF 
      ENDDO
   ENDIF 
ENDIF 
!
!!!   IF(natom(icent) /= 0) THEN 
!!!      IF(isel(ia) == iatom (0, icent) ) THEN 
!!!!                                                                       
!!!!     ----The selected atom is the central atom, check all atoms        
!!!!                                                                       
!!!         in_a = 1 
!!!         in_e = natom (icent) 
!!!         is = cr_iscat (1,iatom (0, icent) ) 
!!!         IF(ABS(mmc_pair(ic, MC_UNI,is,is))==1) THEN               ! "is" is of "left==starting" atom type
!!!!write(*,*) ' Selected is central atom type '
!!!         DO in = in_a, in_e 
!!!            js = cr_iscat (1,iatom (in, icent) ) 
!!!            IF(ABS(mmc_pair(ic, MC_UNI,is,js))==1) THEN            ! "is,is" or "is,js"  with js=="ending"
!!!               IF (check_select_status (iatom (in, icent),  &
!!!                                        .TRUE., cr_prop (iatom (in,  &
!!!                                        icent) ), cr_sel_prop) ) THEN                         
!!!                  ival1 = SIGN(1,mmc_pair(ic, MC_UNI,is,js))
!!!                  mmc_energy_uni = mmc_energy_uni +                  &
!!!                                   mmc_depth (ic,MC_UNI,  0,  0) * ival1
!!!!write(*,*) ' Pair is of type  PAIR                               ', &
!!!! is, js, ival1, mmc_depth (ic,MC_UNI, is, js) * ival1
!!!               ENDIF 
!!!            ELSE
!!!               mmc_energy_uni = mmc_energy_uni +                  &
!!!                                ABS(mmc_depth(ic,MC_UNI,  0,  0))
!!!!write(*,*) ' PAIR is of type  OTHER                              ', &
!!!!  is, js, 1    , ABS(mmc_depth(ic,MC_UNI,  0,  0))
!!!            ENDIF 
!!!         ENDDO 
!!!         ELSE
!!!!write(*,*) ' Selected is of other type                           '
!!!         ENDIF
!!!      ELSE 
!!!!                                                                       
!!!!     The selected atom is a neighbour, use this atom only              
!!!!                                                                       
!!!         in_a = 0 
!!!         in_e = 0 
!!!         is = cr_iscat (1,isel (ia) ) 
!!!         in = 0
!!!                  js = cr_iscat (1,iatom (in, icent) ) 
!!!            IF(mmc_pair(ic, MC_UNI,is,js)/=0) THEN
!!!            IF (check_select_status (iatom (in, icent),  &
!!!                                     .TRUE., cr_prop (iatom (in,  &
!!!                                     icent) ), cr_sel_prop) ) THEN                         
!!!                  ival1 = 0 
!!!                  js = cr_iscat (1,iatom (in, icent) ) 
!!!               IF(mmc_pair(ic, MC_UNI,js,is) /= 0) THEN
!!!                  ival1 = mmc_pair(ic, MC_UNI,js,is)
!!!                  mmc_energy_uni = mmc_energy_uni +                  &
!!!                                   mmc_depth (ic,MC_UNI, 0, 0) * ival1
!!!!write(*,*) ' Selected neighbo', is, js, ival1, mmc_depth (ic,MC_UNI,js,is) * ival1
!!!!write(*,*) ' Unidirec neighbo                                    ', &
!!!! js, is, ival1, mmc_depth (ic,MC_UNI,js,is) * ival1
!!!               ENDIF 
!!!            ENDIF 
!!!            ELSE
!!!                  ival1 = sign(1,mmc_pair(ic, MC_UNI,is,is))
!!!                  mmc_energy_uni = mmc_energy_uni +                  &
!!!                                   mmc_depth (ic,MC_UNI,  0,  0) * ival1
!!!!write(*,*) ' Unidirec other                                      ', &
!!!! js, is, ival1, mmc_depth (ic,MC_UNI, 0, 0) * ival1
!!!            ENDIF 
!!!!              
!!!      ENDIF 
!!!   ENDIF 
!!!ENDIF 
valid_e = .TRUE. 
!write(*,*) ' CALCULATED ENERGY ', mmc_energy_uni
!read(*,*) is
!                                                                       
END FUNCTION mmc_energy_uni                   
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_energy_group(isel, ia, ic, iatom, tatom, icent,  &
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
use precision_mod
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
LOGICAL, DIMENSION(0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: tatom
INTEGER                                              , INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L)                   , INTENT(IN) :: natom
LOGICAL                                              , INTENT(OUT) :: valid_e 
!                                                                       
!                                                                       
INTEGER :: is, js, ind
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
!                                                                       
mmc_energy_group = 0.0 
ncalc   = 0 
valid_e = .FALSE. 
!                                                                       
IF (chem_ctyp(ic) == CHEM_VEC    .OR. &
    chem_ctyp(ic) == CHEM_ENVIR  .OR. &
    chem_ctyp(ic) == CHEM_RANGE  .OR. &
    chem_ctyp(ic) == CHEM_DIST   .OR. &
    chem_ctyp(ic) == CHEM_CON        )   THEN                                               
!                                                                       
   IF(natom(icent) /= 0) THEN              ! We do have neighbors
      in_a = 1 
      in_e = natom (icent) 
      DO ind = in_a, in_e 
         IF(tatom(ind, icent)) THEN  ! Selected atom is central
            is = cr_iscat (1,iatom (0, icent) )
            js = cr_iscat (1,iatom (ind, icent) ) 
         ELSE                        ! Selected atom is a neighbor
            js = cr_iscat (1,iatom (0, icent) )
            is = cr_iscat (1,iatom (ind, icent) ) 
         ENDIF
         IF(mmc_left(ic, MC_GROUP,is)/=0) THEN          ! starting atom is in left  group
            IF(mmc_right(ic, MC_GROUP,js)/=0) THEN      ! ending   atom is in right group (A) => (B)
               IF(check_select_status(iatom (ind, icent),  &
                                      .TRUE., cr_prop (iatom (ind,  &
                                      icent) ), cr_sel_prop) ) THEN                         
                  mmc_energy_group = mmc_energy_group - mmc_depth   (ic, MC_GROUP, is, js) !mmc_depth_def(ic) !(ic,MC_GROUP, is, js)
!write(*,*) ' GROUP ENREGY ', is, js, ic, mmc_depth   (ic, MC_GROUP, is, js) !mmc_depth_def(ic)
!IF(tatom(ind, icent)) THEN  ! Selected atom is central
!write(*,'(2(a,i6, i3), f10.4)') ' MMC_GROUP Central ', iatom (0, icent), is, ' : ', iatom (ind, icent), js, &
!-mmc_depth_def(ic) ! (ic,MC_GROUP, is, js)
!else
!write(*,'(2(a,i6, i3), f10.4)') ' MMC_GROUP NEIG    ', iatom (0, icent), js, ' : ', iatom (ind, icent), is, &
!-mmc_depth_def(ic) ! (ic,MC_GROUP, is, js) 
!endif
               ENDIF 
            ENDIF 
         ENDIF 
      ENDDO
   ENDIF 
ENDIF 
valid_e = .TRUE. 
!                                                                       
END FUNCTION mmc_energy_group                   
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_energy_occ_mol(ianz, imol, amol, valid_e, MAXMOL) 
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
use precision_mod
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_cn(isel, ia, ic, iatom, icent, natom, valid_e,         &
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
         is   = cr_iscat(1,iatom (0, icent) ) 
         js   = -1                            ! Flag negative to indicate no neighbor
         DO in = in_a, in_e 
            IF(check_select_status(iatom(in, icent),  .TRUE.,          &
                                   cr_prop (iatom (in, icent)), cr_sel_prop)) THEN
               js    = cr_iscat(1,iatom (in, icent) ) 
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
!        is = cr_iscat (1,isel (ia) ) 
!        in = 0
!        IF (check_select_status (iatom (in, icent),  &
!                                 .TRUE., cr_prop (iatom (in,  &
!                                 icent) ), cr_sel_prop) ) THEN                         
!           ival1 = 0 
!           js = cr_iscat (1,iatom (in, icent) ) 
!           ival1 = mmc_pair(ic, MC_OCC,is,js)
!           mmc_energy_cn = mmc_energy_cn +                  &
!                            mmc_depth (ic,MC_OCC, 0, 0) * ival1
!        ENDIF 
!              
      ENDIF 
   ENDIF 
ENDIF 
valid_e = .TRUE. 
!write(*,'(a, i3, i6, i3, i3, a, 7i7,a,5i2, f10.2)') &
!' energy_cn   ', ia, isel(ia), icent, natom(icent), ' | ', &
! cr_iscat(1,iatom(0:6, icent)),  ' > ', n_neig(0:4), mmc_energy_cn
DEALLOCATE(n_neig)
!                                                                       
END FUNCTION mmc_energy_cn                   
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION mmc_energy_dis (isel, ia, ic, iatom, icent,  &
      natom, jdir, delta, rdj, valid_e, MAX_ATOM_ENV_L, MMC_MAX_CENT_L,&
      MMC_MAX_ATOM_L)                                 
!+                                                                      
!     Calculates the energy for displacement disorder                       
!                                                                       
!-                                                                      
USE crystal_mod 
USE chem_mod 
USE celltoindex_mod
USE metric_mod
USE mc_mod 
USE mmc_mod 
USE modify_func_mod
use precision_mod
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
REAL(KIND=PREC_DP), DIMENSION(3), INTENT(IN) :: jdir
REAL(KIND=PREC_DP),               INTENT(IN) ::  delta 
REAL(KIND=PREC_DP), DIMENSION(CHEM_MAX_COR), INTENT(IN) :: rdj
!                                                                       
INTEGER :: i, is, js, in, jjs 
INTEGER :: in_a, in_e 
INTEGER :: cell (3), site 
REAL(kind=PREC_DP) :: u (3)
!                                                                       
REAL(kind=PREC_DP) :: dx 
!                                                                       
mmc_energy_dis = 0.0 
valid_e = .FALSE. 
!write(*,*) ' DISPLACEMENT CENTRAL ? ', isel(ia), isel(ia)==iatom(0, icent)
!                                                                       
IF(chem_ctyp(ic) == CHEM_VEC   .OR.  &
   chem_ctyp(ic) == CHEM_ENVIR .OR.  &
   chem_ctyp(ic) == CHEM_RANGE .OR.  &
   chem_ctyp(ic) == CHEM_DIST  .OR.  &
   chem_ctyp(ic) == CHEM_CON        ) THEN                                               
!                                                                       
   IF(natom(icent)  /= 0) THEN 
      IF(isel(ia)  == iatom(0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         in_a = 1 
         in_e = natom(icent) 
         is = cr_iscat(1,iatom(0, icent) ) 
         loop_scat1:DO jjs = 0, cr_nscat 
            IF(mmc_pair(ic, MC_DISP, is, jjs) == -1 ) THEN 
               DO in = in_a, in_e 
                  js = cr_iscat(1,iatom(in, icent) ) 
                  IF(is == js.OR.mmc_pair(ic, MC_DISP, is, js) == -1 ) THEN 
                     IF(check_select_status(iatom(in, icent), .TRUE.,              &
                                            cr_prop(iatom(in, icent)), cr_sel_prop &
                                           ) ) THEN
                        CALL indextocell(iatom(in, icent), cell, site) 
                        DO i = 1, 3 
                           u(i) = cr_pos(i, iatom(in, icent) )    - &
                                  chem_ave_pos(i, site)          - &
                                  REAL(cell(i) - 1, kind=PREC_DP) - &
                                  cr_dim0(i, 1)                           
                        ENDDO 
                        dx = skalpro(u, jdir, cr_gten) / rdj(ic) 
!                                                                       
!do i=1, 3
!write(*,'(11f10.4)') cr_pos(i, iatom(in, icent)), chem_ave_pos(i, site  ),  &
!              REAL(cell(i) - 1, kind=PREC_DP), cr_dim0(i, 1), u(i) , &
!              0.0,0.0,0.0, jdir
!enddo
!write(*,*) 'DX ', dx, mmc_depth(ic, MC_DISP, 0,0) *delta*dx
                        mmc_energy_dis = mmc_energy_dis + mmc_depth(ic,&
                        MC_DISP, 0, 0) * delta * dx                     
                        valid_e = .TRUE. 
                     ENDIF 
                  ENDIF 
               ENDDO 
            ENDIF 
         ENDDO loop_scat1
      ELSE 
!                                                                       
!     The selected atom is a neighbour, use this atom only              
!                                                                       
           in_a = 0 
           in_e = 0 
           is = cr_iscat(1,isel(ia) ) 
           DO jjs = 0, cr_nscat 
               IF(mmc_pair(ic, MC_DISP, is, jjs) == -1 ) THEN 
                  in = 0 
                  js = cr_iscat(1, iatom(in, icent) ) 
                  IF(is == js.OR.mmc_pair(ic, MC_DISP, is, js) == -1 ) THEN 
                     IF(check_select_status(iatom(in, icent),           &
                                            .TRUE.,                     &
                                            cr_prop(iatom(in, icent) ), &
                                            cr_sel_prop) ) THEN
                     CALL indextocell(iatom(in, icent), cell, site) 
                     DO i = 1, 3 
                        u(i) = cr_pos(i, iatom(in, icent) )    -   &
                               chem_ave_pos(i, site)           -   &
                               REAL(cell(i) - 1, kind=PREC_DP) -   &
                               cr_dim0(i, 1)                           
                     ENDDO 
                     dx = skalpro(u, jdir, cr_gten) / rdj(ic) 
!                                                                       
                     mmc_energy_dis = mmc_energy_dis + mmc_depth(ic,&
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_spr (isel, ia, ic, iatom, patom, icent,  &
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
use precision_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: patom
INTEGER                                                   , INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT)                         , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL(kind=PREC_DP)    :: d, u (3), v (3) 
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
               is = cr_iscat (1,iatom (0, icent) ) 
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (1,isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (1,iatom (in, icent) ) 
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_spr_mol (nmol, imol, valid_e, MAX_ATOM_ENV_L) 
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
use precision_mod
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
REAL(kind=PREC_DP)  :: d, u (3), v (3) 
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_len (isel, ia, ic, iatom, patom, icent,  &
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
use precision_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L) , INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) :: valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL(kind=PREC_DP) :: d, u (3), v (3) 
!                                                                       
!QWwrite(*,*) ' DIMENSIONS  ', MAX_ATOM_ENV_L, MMC_MAX_CENT_L, MMC_MAX_ATOM_L
mmc_energy_len = 0.0 
ncalc = 0 
valid_e = .FALSE. 
!                                                                       
IF(chem_ctyp(ic) == CHEM_VEC   .OR.  &
    chem_ctyp(ic) == CHEM_ENVIR .OR.  &
    chem_ctyp(ic) == CHEM_RANGE .OR.  &
    chem_ctyp(ic) == CHEM_DIST  .OR.  &
    chem_ctyp(ic) == CHEM_CON        ) THEN                                               
!                                                                       
   IF(natom(icent)  /= 0) THEN 
      IF(isel(ia)  == iatom(0, icent) ) THEN 
!                                                                       
!     ----The selected atom is the central atom, check all atoms        
!                                                                       
         DO i = 1, 3 
               u(i) = patom(i, 0, icent) 
         ENDDO 
         in_a = 1 
         in_e = natom(icent) 
         is = cr_iscat(1,iatom(0, icent) ) 
      ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
         DO i = 1, 3 
            u(i) = cr_pos(i, isel(ia) ) 
         ENDDO 
         in_a = 0 
         in_e = 0 
         is = cr_iscat(1,isel(ia) ) 
      ENDIF 
      DO in = in_a, in_e 
         js = cr_iscat(1,iatom(in, icent) ) 
         IF(mmc_target_corr(ic, MC_LENNARD, is, js)  /= 0.0) THEN 
            IF(check_select_status(iatom(in, icent),  &
                                        .TRUE., cr_prop(iatom(in,     &
            icent) ), cr_sel_prop) ) THEN                            
            DO i = 1, 3 
               v(i) = patom(i, in, icent) 
            ENDDO 
            d = do_blen(.TRUE., u, v) 
!                                                                       
!QWwrite(*,*) ' ENERGY ', is, js, d,  mmc_target_corr(ic, MC_LENNARD, is, js),&
!QW                  d-mmc_target_corr(ic, MC_LENNARD, is, js), &
!QW                                                    mmc_depth(ic,      &
!QW                  MC_LENNARD, is, js) *(mmc_len_a(ic, is, js) / d**   &
!QW                  mmc_len_m(ic, is, js) - mmc_len_b(ic, is, js)       &
!QW                  / d**mmc_len_n(ic, is, js) )
!DEPTH         mmc_energy_len = mmc_energy_len +                                   &
!DEPTH                          mmc_depth(ic, MC_LENNARD, is, js) *                &
!DEPTH                         (mmc_len_a(ic, is, js) / d**mmc_len_m(ic, is, js) - &
!DEPTH                          mmc_len_b(ic, is, js) / d**mmc_len_n(ic, is, js)  )                         
               mmc_energy_len = mmc_energy_len +                                   &
                               (mmc_len_a(ic, is, js) / d**mmc_len_m(ic, is, js) - &
                                mmc_len_b(ic, is, js) / d**mmc_len_n(ic, is, js)  )                         
               ncalc = ncalc + 1 
!                                                                       
            ENDIF 
         ENDIF 
      ENDDO 
   ENDIF 
ENDIF 
!
IF(ncalc > 0) THEN 
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_rep (isel, ia, ic, iatom, patom, icent,  &
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
use precision_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL(kind=PREC_DP)  :: d, u (3), v (3) 
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
               is = cr_iscat (1,iatom (0, icent) ) 
!!write(*,*) 'CENT ', isel(ia), cr_iscat(1,isel(ia)), (cr_iscat(1,iatom(i, icent)),i=1, natom(icent))
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (1,isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (1,iatom (in, icent) ) 
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_buck (isel, ia, ic, iatom, patom, icent, &
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
use precision_mod
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
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js, in
INTEGER :: in_a, in_e 
INTEGER :: ncalc 
REAL(kind=PREC_DP) :: d, u (3), v (3) 
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
               is = cr_iscat (1,iatom (0, icent) ) 
            ELSE 
!                                                                       
!     The selcted atom is a neighbour, use this atom only               
!                                                                       
               DO i = 1, 3 
               u (i) = cr_pos (i, isel (ia) ) 
               ENDDO 
               in_a = 0 
               in_e = 0 
               is = cr_iscat (1,isel (ia) ) 
            ENDIF 
            DO in = in_a, in_e 
            js = cr_iscat (1,iatom (in, icent) ) 
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
REAL(kind=PREC_DP) FUNCTION mmc_energy_angle (ic, iatom, patom, icent,    &
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
use precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: MAX_ATOM_ENV_L
INTEGER, INTENT(IN) :: MMC_MAX_CENT_L
INTEGER, INTENT(IN) :: MMC_MAX_ATOM_L
INTEGER, INTENT(IN) :: ic
INTEGER, DIMENSION(   0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: iatom
REAL(kind=PREC_DP)   , DIMENSION(3, 0:MAX_ATOM_ENV_L, MMC_MAX_CENT_L), INTENT(IN) :: patom
INTEGER, INTENT(IN) :: icent 
INTEGER, DIMENSION(MMC_MAX_CENT_L) , INTENT(IN) :: natom
LOGICAL, INTENT(OUT) ::  valid_e 
!                                                                       
INTEGER :: i, is, js
INTEGER :: ii, jj 
LOGICAL :: lnoneig 
REAL(kind=PREC_DP) :: a, b, u (3), v (3), w (3) 
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
            is = cr_iscat (1,iatom (ii, icent) ) 
            js = cr_iscat (1,iatom (jj, icent) ) 
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
!           (1,isel (ia) )                                                
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
!              is = cr_iscat (1,iatom (0, icent) ) 
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
!              is = cr_iscat (1,isel (ia) ) 
!           ENDIF 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) ' central           ', u 
!           ENDIF 
!           DO in = in_a, in_e 
!           IF (dbg) THEN 
!     WRITE ( * ,  * ) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' 
!           ENDIF 
!           js = cr_iscat (1,iatom (in, icent) ) 
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
use precision_mod
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
REAL(kind=PREC_DP)               :: d 
REAL(kind=PREC_DP), DIMENSION(3) :: u (3), nullv(3) 
!                                                                       
DATA nullv/ 0.0, 0.0, 0.0 / 
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
         d = do_blen (lspace, u, nullv) 
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
REAL(kind=PREC_DP) :: r1
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
REAL(kind=PREC_DP) FUNCTION buckingham(a, rho, b, x) 
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
use mmc_basic_mod
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
REAL(KIND=prec_DP), DIMENSION(3) :: x     ! Atom position
REAL(KIND=prec_DP), DIMENSION(3) :: werte ! Atom types allowed
REAL(kind=PREC_DP) :: rmin  ! minimum distance for find_env
REAL(kind=PREC_DP) :: rmax  ! maximum distance for find_env
REAL(kind=PREC_DP) :: rel_cycl    ! how far are we in the desired number of cycles
LOGICAL :: lout_feed, done
REAL(kind=PREC_DP) :: r1
real(kind=PREC_DP), dimension(2)                    :: maxdev =(/0.0, 0.0/)
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
      nneig(cr_iscat(1,atom_env(j))) = nneig(cr_iscat(1,atom_env(j))) + 1
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
   cr_iscat(1,iatom) = is_max
   IF(MOD(INT(i,PREC_INT_LARGE), mo_feed)==0) THEN
      CALL mmc_correlations (lout_feed, rel_cycl, done, .FALSE., .TRUE., maxdev)
   ENDIF
ENDDO
!
END SUBROUTINE mmc_grow
!
!*****7*************************************************************************
!
!
!*****7*************************************************************************
!
END MODULE mmc_menu
