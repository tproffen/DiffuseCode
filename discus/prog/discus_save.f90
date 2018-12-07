MODULE save_menu
!
CONTAINS
!+                                                                      
!     Menu for the structure saveing part of DISCUS. The user           
!     can define what structural information is saves to the structure  
!     file.                                                             
!     If the Menu is called with one or two parameters, the structure   
!     is written immedeately to the file defined by the first parameter.
!     In this case, the format is taken from the current settings.      
!                                                                       
!*****7*****************************************************************
SUBROUTINE save_struc (string, lcomm) 
!-                                                                      
!     Main menu for generalized transformation operations               
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE modify_mod
      USE discus_save_mod 
      USE discus_show_menu
      USE prop_char_mod
      USE prop_para_func
!
      USE ber_params_mod
      USE build_name_mod
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE get_params_mod
      USE learn_mod 
      USE class_macro_internal 
      USE prompt_mod 
      USE sup_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA =  20 ! A command requires at least these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
!
      CHARACTER ( LEN=*    ) string 
      CHARACTER ( LEN=1024 ) zeile 
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt 
      CHARACTER(1024) line
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz
      INTEGER lcomm, sav_flen 
      LOGICAL lend 
!      REAL, DIMENSION(SAV_MAXSCAT) :: repl ! Dummy variable needed for atom_select
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      DATA sav_flen / 1 / 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
!
!     Interpret parameters used by 'save' command                       
!                                                                       
      CALL get_params (string, ianz, cpara, lpara, maxw, lcomm) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.gt.0) THEN 
!                                                                       
!     ----Parameters were found, write file immedeatly, and return      
!         to calling routine                                            
!                                                                       
            CALL save_check_alloc()
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.eq.0) THEN 
               sav_file = cpara (1) (1:lpara(1))
               sav_flen = lpara (1) 
!                                                                       
               IF (sav_keyword) THEN 
!                                                                       
!-----      ------Write new type of structur file                       
!                                                                       
                  CALL save_keyword (cpara (1) ) 
               ELSE 
!                                                                       
!-----      ------Write old type of structur file                       
!                                                                       
                  CALL save_nokeyword (cpara (1) ) 
               ENDIF 
            ENDIF 
!                                                                       
!     ----File was written, return to caller                            
!                                                                       
            RETURN 
         ENDIF 
      ELSE 
         RETURN 
      ENDIF 
!********************************************************************** 
!                                                                       
!     Beginn of the normal menu structure                               
!                                                                       
!                                                                       
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/save' 
!                                                                       
main: DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) THEN 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') THEN 
                  IF (length.ge.2) THEN 
                     CALL file_kdo (line (2:length), length - 1) 
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
               ELSEIF (str_comp (befehl, 'chem', 2, lbef, 4) ) THEN 
                  CALL show_chem 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
                  CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) THEN                                      
                  line = zeile 
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                     lp = lp + 7 
                     CALL do_hel ('discus '//line, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus save '//line, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
                  IF (zeile.ne.' ') THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 2, lbef, 4) ) THEN 
                  CALL do_input (zeile, lp) 
!
!     ----Reset save menu 'rese'n'                                      
!
               ELSEIF (str_comp (befehl, 'rese', 2, lbef, 4) ) THEN 
                  CALL save_reset
               ELSE
!
!---------SAVE specific commands
!                                                                       
               CALL save_check_alloc()
               IF ( ier_num < 0 ) THEN
                  RETURN
               ENDIF
!                                                                       
!     ----Deselect which atoms are included in the wave 'dese'          
!                                                                       
               IF (str_comp (befehl, 'dese', 1, lbef, 4) ) THEN 
                   CALL atom_select (zeile, lp, 0, SAV_MAXSCAT, sav_latom, &
                   sav_sel_atom, .false., .false.)              
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num.eq.0) THEN 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num.eq.0) THEN 
!                       IF (werte (1) .eq. - 1) THEN 
!                          DO i = 0, cr_nscat 
!                          sav_latom (i) = .false. 
!                          ENDDO 
!                       ELSE 
!                          DO i = 1, ianz 
!                          sav_latom (nint (werte (i) ) ) = .false. 
!                          ENDDO 
!                       ENDIF 
!                    ENDIF 
!                 ENDIF 
!                                                                       
!     define format of output file 'format'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'format', 1, lbef, 6) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) THEN 
                     IF (str_comp (cpara (1) , 'keyword', 1, lpara (1) ,&
                     7) ) THEN                                          
                        sav_keyword = .true. 
                     ELSEIF (str_comp (cpara (1) , 'nokeywo', 1, lpara (&
                     1) , 7) ) THEN                                     
                        sav_keyword = .false. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select range of atoms within crystal to be included 'incl'    
!                                                                       
               ELSEIF (str_comp (befehl, 'incl', 1, lbef, 4) ) THEN 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) THEN 
                     IF (ianz.eq.2) THEN 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) THEN 
                           sav_start = nint (werte (1) ) 
                           sav_end   = nint (werte (2) ) 
                           IF(sav_end<sav_start .or. sav_start<1) THEN
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                              ier_msg(1)='The upper limit must be larger than lower'
                              ier_msg(2)='The lower limit must be larger than zero'
                           ENDIF 
                        ENDIF 
                     ELSEIF (ianz.eq.1) THEN 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) THEN                                       
                           sav_start = 1 
                           sav_end = - 1 
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
!     ----omit parameters from the keyword list  'omit'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'omit', 1, lbef, 4) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) THEN 
                     IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3)) THEN
                        sav_w_gene = .false. 
                        sav_w_ncell= .false. 
                        sav_w_symm = .false. 
                        sav_w_scat = .false. 
                        sav_w_adp  = .false. 
                        sav_w_occ  = .false. 
                        sav_w_obje = .false. 
                        sav_w_mole = .false. 
                        sav_w_doma = .false. 
                        sav_w_prop = .false. 
                     ELSEIF(str_comp(cpara(1), 'gene', 1, lpara(1), 4) ) THEN
                        sav_w_gene = .false. 
                     ELSEIF(str_comp(cpara (1), 'mole', 1, lpara(1), 4) ) THEN
                        sav_w_mole = .false. 
                     ELSEIF(str_comp(cpara (1), 'obje', 1, lpara(1), 4) ) THEN
                        sav_w_obje = .false. 
                     ELSEIF(str_comp(cpara (1), 'doma', 1, lpara(1), 4) ) THEN
                        sav_w_doma = .false. 
                     ELSEIF(str_comp(cpara (1), 'ncell', 1, lpara(1), 5) ) THEN
                        sav_w_ncell = .false. 
                     ELSEIF(str_comp(cpara (1), 'symm', 2, lpara(1), 4) ) THEN
                        sav_w_symm = .false. 
                     ELSEIF(str_comp(cpara (1), 'scat', 2, lpara(1), 4) ) THEN
                        sav_w_scat = .false. 
                     ELSEIF(str_comp(cpara (1), 'adp', 1, lpara(1), 3) ) THEN
                        sav_w_adp = .false. 
                     ELSEIF(str_comp(cpara (1), 'occ', 1, lpara(1), 3) ) THEN
                        sav_w_occ = .false. 
                     ELSEIF(str_comp(cpara (1), 'prop', 1, lpara(1), 4) ) THEN
                        sav_w_prop = .false. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!------ --Handle property settings 'property'                           
!                                                                       
               ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) then
!                                                                       
                  CALL property_select (zeile, lp,sav_sel_prop)
!                                                                       
!     define name of output file 'outfile'                              
!                                                                       
               ELSEIF (str_comp (befehl, 'outfile', 1, lbef, 7) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) THEN 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
                     IF (ier_num.eq.0) THEN 
                        sav_file = cpara (1) (1:lpara(1))
                        sav_flen = lpara (1) 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----run transformation 'run'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) THEN 
                  IF (sav_keyword) THEN 
                     IF(sav_w_scat .EQV. sav_w_adp .AND. &
                        sav_w_scat .EQV. sav_w_occ .AND. &
                        sav_w_adp  .EQV. sav_w_occ )  THEN
                     IF (str_comp (sav_file(1:8),'internal',8,8,8)) THEN
                        CALL save_internal (sav_file) 
                     ELSE 
                        CALL save_keyword (sav_file) 
                     ENDIF 
                     ELSE 
                        ier_num = -159
                        ier_typ = ER_APPL
                        ier_msg(1) = 'Write/Omit differ for SCAT, ADP, OCC'
                        ier_msg(2) = 'Check write/omit statements in save'
                        ier_msg(3) = 'Segments and homogenize'
                     ENDIF 
                  ELSE 
                     CALL save_nokeyword (sav_file) 
                  ENDIF 
!                                                                       
!     ----Select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'sele', 2, lbef, 4) ) THEN 
                   CALL atom_select (zeile, lp, 0, SAV_MAXSCAT, sav_latom, &
                   sav_sel_atom, .false., .true.)               
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num.eq.0) THEN 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num.eq.0) THEN 
!                       IF (werte (1) .eq. - 1) THEN 
!                          DO i = 0, cr_nscat 
!                          sav_latom (i) = .true. 
!                          ENDDO 
!                       ELSE 
!                          DO i = 1, ianz 
!                          sav_latom (nint (werte (i) ) ) = .true. 
!                          ENDDO 
!                       ENDIF 
!                    ENDIF 
!                 ENDIF 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
                  CALL save_show
!                                                                       
!     ----write parameters from the keyword list  'write'               
!                                                                       
               ELSEIF (str_comp (befehl, 'write', 2, lbef, 5) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) THEN 
                     IF(str_comp(cpara(1), 'all', 1, lpara(1), 3) ) THEN                                             
                        sav_w_gene = .true. 
                        sav_w_ncell= .true. 
                        sav_w_symm = .true. 
                        sav_w_scat = .true. 
                        sav_w_adp  = .true. 
                        sav_w_occ  = .true. 
                        sav_w_mole = .true. 
                        sav_w_obje = .true. 
                        sav_w_doma = .true. 
                        sav_w_prop = .true. 
                     ELSEIF(str_comp(cpara (1), 'gene', 1, lpara(1), 4) ) THEN
                        sav_w_gene = .true. 
                     ELSEIF(str_comp(cpara (1), 'mole', 1, lpara(1), 4) ) THEN
                        sav_w_mole = .true. 
                     ELSEIF(str_comp(cpara (1), 'obje', 1, lpara(1), 4) ) THEN
                        sav_w_obje = .true. 
                     ELSEIF(str_comp(cpara (1), 'doma', 1, lpara(1), 4) ) THEN
                        sav_w_doma = .true. 
                     ELSEIF(str_comp(cpara (1), 'ncell', 1, lpara(1), 5) ) THEN
                        sav_w_ncell = .true. 
                     ELSEIF(str_comp(cpara (1), 'symm', 2, lpara(1), 4) ) THEN
                        sav_w_symm = .true. 
                     ELSEIF(str_comp(cpara (1), 'scat', 2, lpara(1), 4) ) THEN
                        sav_w_scat = .true. 
                     ELSEIF(str_comp(cpara (1), 'adp', 1, lpara(1),3) ) THEN
                        sav_w_adp = .true. 
                     ELSEIF(str_comp(cpara (1), 'occ', 1, lpara(1),3) ) THEN
                        sav_w_occ = .true. 
                     ELSEIF(str_comp(cpara (1), 'prop', 1, lpara(1), 4) ) THEN
                        sav_w_prop = .true. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
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
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in save menu'
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
            CYCLE main
         ELSE
            EXIT main
         ENDIF
      ENDDO  main
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE save_struc                     
!
!*****7*****************************************************************
!
      SUBROUTINE save_check_alloc()
!
      USE crystal_mod
      USE discus_allocate_appl_mod
      USE discus_save_mod 
      USE errlist_mod
!
      INTEGER                      :: n_nscat
!
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
         CALL alloc_save (  n_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      END SUBROUTINE save_check_alloc
!*****7*****************************************************************
      SUBROUTINE save_nokeyword (strucfile) 
!+                                                                      
!     This subroutine saves the structure and/or the unit cell          
!     onto a file in the format used until DISCUS 3.0 i.e. without      
!     keywords.                                                         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER ( * ) strucfile 
      INTEGER ist, i, j 
      LOGICAL lread 
!                                                                       
      INTEGER len_str 
!                                                                       
      DATA ist / 7 / 
!                                                                       
      lread = .false. 
      CALL oeffne (ist, strucfile, 'unknown') 
      IF (ier_num.eq.0) THEN 
!                                                                       
!-----      --Write old type of structur file                           
!                                                                       
         j = len_str (cr_name) 
         WRITE (ist, 3) cr_name (1:j) 
         WRITE (ist, 33) cr_spcgr 
         WRITE (ist, 34) (cr_a0 (i), i = 1, 3), (cr_win (i), i = 1, 3) 
         DO i = 1, cr_natoms 
         WRITE (ist, 4) cr_at_lis (cr_iscat (i) ), (cr_pos (j, i),      &
         j = 1, 3), cr_dw (cr_iscat (i) )                               
         ENDDO 
      ENDIF 
      CLOSE (ist) 
!                                                                       
    3 FORMAT (a) 
   33 FORMAT  (a16) 
   34 FORMAT  (6(1x,f10.6)) 
    4 FORMAT (a4,3(2x,f14.6),5x,f8.4) 
      END SUBROUTINE save_nokeyword                 
!********************************************************************** 
      SUBROUTINE save_keyword (strucfile) 
!+                                                                      
!     This subroutine saves the structure and/or the unit cell          
!     onto a file. The format uses keyword description.                 
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE gen_add_mod 
      USE modify_func_mod
      USE molecule_mod 
      USE sym_add_mod 
      USE discus_save_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER ( LEN=* ) , INTENT(in) :: strucfile 
      CHARACTER(31) fform 
      CHARACTER(15) C_MOLE ( - 4:4) 
      INTEGER ist, i, j, k
      INTEGER i_start, i_end 
      INTEGER is, ie 
      INTEGER ::   wr_prop = 1
      INTEGER ::   wr_mole = 0
      INTEGER ::   wr_cont = 0
      LOGICAL lread 
      LOGICAL                            :: lsave
      LOGICAL, DIMENSION(:), ALLOCATABLE :: lwrite ! flag if atom needs write
!                                                                       
      INTEGER len_str 
!                                                                       
!                                                                       
      DATA ist / 67 / 
      DATA C_MOLE / 'domain_fuzzy   ', 'domain_sphere  ', 'domain_cylinder',&
                    'domain_cube    ', 'atoms          ', 'cube           ',&
                    'cylinder       ', 'sphere         ', 'cube           ' /             
!
!     Test if any atom type is selected for write
!
      lsave = .true.
      DO i = 0, cr_nscat 
         lsave = lsave .or. sav_latom(i)
      ENDDO 
      IF(.not. lsave ) THEN
         ier_num = -58
         ier_typ = ER_APPL
         ier_msg(1) = 'No atom types were selected for write'
         ier_msg(2) = 'The output file is not written'
         ier_msg(3) = 'Select atoms, or set error to live'
         RETURN
      ENDIF
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sav_start 
      i_end = sav_end 
      IF (sav_end.eq. - 1) i_end = cr_natoms 
!                                                                       
      lread = .false. 
      CALL oeffne (ist, strucfile, 'unknown') 
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
!                                                                       
!-----      --Write new type of structur file                           
!                                                                       
      j = len_str (cr_name) 
      IF (j.eq.0) THEN 
         cr_name = ' ' 
         j = 1 
      ENDIF 
      WRITE (ist, 3000) cr_name (1:j) 
      IF (spcgr_ianz.eq.0) THEN 
         WRITE (ist, 3010) cr_spcgr 
      ELSEIF (spcgr_ianz.eq.1) THEN 
         WRITE (ist, 3011) cr_spcgr, spcgr_para 
      ENDIF 
      WRITE (ist, 3020) (cr_a0 (i), i = 1, 3), (cr_win (i), i = 1, 3) 
      IF (sav_w_scat) THEN 
         j = (cr_nscat - 1) / 7 
         DO i = 1, j 
         is = (i - 1) * 7 + 1 
         ie = is + 6 
         WRITE (ist, 3110) (cr_at_lis (k), k = is, ie) 
         ENDDO 
         IF (cr_nscat - j * 7.eq.1) THEN 
            WRITE (ist, 3111) cr_at_lis (cr_nscat) 
         ELSEIF (cr_nscat - j * 7.gt.1) THEN 
            WRITE (fform, 7010) cr_nscat - j * 7 - 1 
            WRITE (ist, fform) (cr_at_lis (i), i = j * 7 + 1, cr_nscat) 
         ENDIF 
      ENDIF 
      IF (sav_w_adp) THEN 
         j = (cr_nscat - 1) / 7 
         DO i = 1, j 
         is = (i - 1) * 7 + 1 
         ie = is + 6 
         WRITE (ist, 3120) (cr_dw (k), k = is, ie) 
         ENDDO 
         IF (cr_nscat - j * 7.eq.1) THEN 
            WRITE (ist, 3121) cr_dw (cr_nscat) 
         ELSEIF (cr_nscat - j * 7.gt.1) THEN 
            WRITE (fform, 7020) cr_nscat - j * 7 - 1 
            WRITE (ist, fform) (cr_dw (i), i = j * 7 + 1, cr_nscat) 
         ENDIF 
      ENDIF 
      IF (sav_w_occ) THEN 
         j = (cr_nscat - 1) / 7 
         DO i = 1, j 
         is = (i - 1) * 7 + 1 
         ie = is + 6 
         WRITE (ist, 3220) (cr_occ(k), k = is, ie) 
         ENDDO 
         IF (cr_nscat - j * 7.eq.1) THEN 
            WRITE (ist, 3221) cr_occ(cr_nscat) 
         ELSEIF (cr_nscat - j * 7.gt.1) THEN 
            WRITE (fform, 7030) cr_nscat - j * 7 - 1 
            WRITE (ist, fform) (cr_occ(i), i = j * 7 + 1, cr_nscat) 
         ENDIF 
      ENDIF 
      IF (sav_w_gene) THEN 
         DO k = 1, gen_add_n 
         WRITE (ist, 3021) ( (gen_add (i, j, k), j = 1, 4), i = 1, 3),  &
         gen_add_power (k)                                              
         ENDDO 
      ENDIF 
      IF (sav_w_symm) THEN 
         DO k = 1, sym_add_n 
         WRITE (ist, 3022) ( (sym_add (i, j, k), j = 1, 4), i = 1, 3),  &
         sym_add_power (k)                                              
         ENDDO 
      ENDIF 
      IF (sav_w_ncell) THEN 
         WRITE (ist, 3030) cr_icc, cr_ncatoms , cr_natoms
      ENDIF 
      WRITE (ist, 3900) 
!
      ALLOCATE(lwrite(1:cr_natoms))
      lwrite(:) = .true.
!                                                                       
!     Write content of objects                                          
!                                                                       
      IF (sav_w_obje) THEN 
         DO i = 1, mole_num_mole 
         IF (mole_char (i) .gt.MOLE_ATOM) THEN 
            WRITE (ist, 4000) 'object' 
            WRITE (ist, 4002) 'object', mole_type (i) 
            WRITE (ist, 4100) 'object', c_mole (mole_char (i) ) 
            WRITE (ist, 4200) 'object', mole_dens (i) 
!           DO j = 1, mole_len (i) 
!           k = mole_cont (mole_off (i) + j) 
!           WRITE (ist, 4) cr_at_lis (cr_iscat (k) ), (cr_pos (l, k),   &
!           l = 1, 3), cr_dw (cr_iscat (k) ), cr_prop (k)               
!           lwrite(k) = .false.
!           ENDDO 
            WRITE (ist, 4900) 'object' 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Write content of micro domains                                    
!                                                                       
      IF (sav_w_doma) THEN 
         DO i = 1, mole_num_mole 
         IF (mole_char (i) .lt.MOLE_ATOM) THEN 
            WRITE (ist, 4000) 'domain' 
            WRITE (ist, 4002) 'domain', mole_type (i) 
            WRITE (ist, 4100) 'domain', c_mole (mole_char (i) ) 
            k = len_str (mole_file (i) ) 
            IF (k.gt.0) THEN 
               WRITE (ist, 4300) 'domain', mole_file (i) (1:k) 
            ELSE 
               WRITE (ist, 4300) 'domain' 
            ENDIF 
            WRITE (ist, 4400) 'domain', mole_fuzzy (i) 
!           DO j = 1, mole_len (i) 
!           k = mole_cont (mole_off (i) + j) 
!           WRITE (ist, 4) cr_at_lis (cr_iscat (k) ), (cr_pos (l, k),   &
!           l = 1, 3), cr_dw (cr_iscat (k) ), cr_prop (k)               
!           lwrite(k) = .false.
!           ENDDO 
            WRITE (ist, 4900) 'domain' 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Write content of molecules                                        
!                                                                       
      IF (sav_w_mole) THEN 
         DO i = 1, mole_num_mole 
         IF (mole_char (i) .eq.MOLE_ATOM) THEN 
            WRITE (ist, 4000) 'molecule' 
            WRITE (ist, 4002) 'molecule', mole_type (i) 
            WRITE (ist, 4100) 'molecule', c_mole (mole_char (i) ) 
            WRITE (ist, 4500) 'molecule', mole_biso (mole_type(i))
            WRITE (ist, 4600) 'molecule', mole_clin (mole_type(i))
            WRITE (ist, 4700) 'molecule', mole_cqua (mole_type(i))
!           DO j = 1, mole_len (i) 
!           k = mole_cont (mole_off (i) + j) 
!           WRITE (ist, 4) cr_at_lis (cr_iscat (k) ), (cr_pos (l, k),   &
!           l = 1, 3), cr_dw (cr_iscat (k) ), cr_prop (k)               
!           lwrite(k) = .false.
!           ENDDO 
            WRITE (ist, 4900) 'molecule' 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Write atoms                                                       
!                                                                       
      DO i = i_start, i_end 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
!     IF (sav_latom (cr_iscat (i) ) ) THEN 
      IF (check_select_status (i, sav_latom (cr_iscat (i) ), cr_prop (i),   &
                               sav_sel_prop)     ) THEN
!        IF(lwrite(i)) THEN
         wr_prop = 1
         wr_mole = 0
         wr_cont = 0
         IF(sav_w_prop) wr_prop = cr_prop(i)
         IF (sav_w_mole .OR. sav_w_doma .OR. sav_w_obje) THEN 
            IF(cr_mole(i)/=0) THEN
               wr_mole = cr_mole(i)
               check_mole: DO j = 1, mole_len (cr_mole(i))
                  IF(mole_cont (mole_off(cr_mole(i))+j) == i) THEN
                     wr_cont = j
                     EXIT check_mole
                  ENDIF
               ENDDO check_mole
            ENDIF
         ENDIF
         WRITE (ist, 4) cr_at_lis (cr_iscat (i) ),         &
                        (cr_pos (j, i),j = 1, 3),          &
                        cr_dw (cr_iscat (i) ), wr_prop,    &
                        wr_mole, wr_cont, cr_occ(cr_iscat(i))
!           IF(sav_w_prop) THEN
!              WRITE (ist, 4) cr_at_lis (cr_iscat (i) ),         &
!                             (cr_pos (j, i),j = 1, 3),          &
!                             cr_dw (cr_iscat (i) ), cr_prop (i)               
!           ELSE
!              WRITE (ist, 4) cr_at_lis (cr_iscat (i) ),         &
!                             (cr_pos (j, i),j = 1, 3),          &
!                             cr_dw (cr_iscat (i) ), 1
!           ENDIF 
!        ENDIF 
      ENDIF 
      ENDDO 
!
      CLOSE (ist) 
      DEALLOCATE(lwrite)
!                                                                       
 3000 FORMAT    ('title ',a) 
 3010 FORMAT    ('spcgr ',a16) 
 3011 FORMAT    ('spcgr ',a16,',',i4) 
 3020 FORMAT    ('cell  ',5(f10.6,','),f10.6) 
 3021 FORMAT    ('gene  ',12(f9.6,','),i3) 
 3022 FORMAT    ('symm  ',12(f9.6,','),i3) 
 3030 FORMAT    ('ncell ',3(i8,','),i10,',',i12) 
 3110 FORMAT    (('scat  ', a4,6(',',5x,a4))) 
 3111 FORMAT    (('scat  ', a4)) 
 3120 FORMAT    (('adp   ', f9.6,6(',',5x,f9.6))) 
 3121 FORMAT    (('adp   ', f9.6)) 
 3220 FORMAT    (('occ   ', f9.6,6(',',5x,f9.6))) 
 3221 FORMAT    (('occ   ', f9.6)) 
 3900 FORMAT    ('atoms      x,',14x,'y,',14x,'z,',13x,'Biso,', 4x,'Property,', &
                 2x,'MoleNo,  MoleAt,   Occ')
 4000 FORMAT    (a) 
 4002 FORMAT    (a,' type,',i8) 
 4100 FORMAT    (a,' character,',a15) 
 4200 FORMAT    (a,' density  ,',f12.4) 
 4300 FORMAT    (a,' file     ,',a) 
 4400 FORMAT    (a,' fuzzy    ,',f12.4) 
 4500 FORMAT    (a,' biso     ,',f12.4)
 4600 FORMAT    (a,' clin     ,',f12.8)
 4700 FORMAT    (a,' cqua     ,',f12.8)
 4900 FORMAT    (a,' end') 
    4 FORMAT (a4,3(1x,f14.6,','),4x,f10.6,',',i8, ',', I8, ',', I8,', ', F10.6) 
 7010 FORMAT    ('(''scat  '', a4  ,',i1,'('','',5x,a4  ))') 
 7020 FORMAT    ('(''adp   '', f9.6,',i1,'('','',5x,f9.6))') 
 7030 FORMAT    ('(''occ   '', f9.6,',i1,'('','',5x,f9.6))') 
      END SUBROUTINE save_keyword                   
!********************************************************************** 
      SUBROUTINE save_internal (strucfile) 
!+                                                                      
!     This subroutine saves the structure and/or the unit cell          
!     onto a file. The format uses keyword description.                 
!+                                                                      
!     USE discus_config_mod 
      USE crystal_mod 
!     USE gen_add_mod 
      USE class_internal
!     USE molecule_mod 
!     USE sym_add_mod 
      USE discus_save_mod 
      USE errlist_mod
      IMPLICIT none 
!                                                                       
      CHARACTER ( LEN=* ), INTENT(IN) :: strucfile 
!
      LOGICAL                         :: lsave
      INTEGER                         :: i
      INTEGER                         :: istatus
!
      CALL save_check_alloc()
      IF ( ier_num < 0 ) THEN
         RETURN
      ENDIF
!
!     Test if any atom type is selected for write
!
      lsave = .true.
      DO i = 0, cr_nscat 
         lsave = lsave .or. sav_latom(i)
      ENDDO 
      IF(.not. lsave ) THEN
         ier_num = -58
         ier_typ = ER_APPL
         ier_msg(1) = 'No atom types were selected for write'
         ier_msg(2) = 'The output file is not written'
         ier_msg(3) = 'Select atoms, or set error to live'
         RETURN
      ENDIF
!                                                                       
IF(ASSOCIATED(store_root)) THEN    ! The tree exists already, add a node
!
      ALLOCATE(store_temp, STAT = istatus)        ! Allocate a temporary storage
      NULLIFY(store_temp%before)
      NULLIFY(store_temp%after )
      IF ( istatus /= 0) THEN
         ier_num = -114
         ier_typ = ER_APPL
         ier_msg(1) = 'Temporary storage failed'
         RETURN
      ENDIF
      store_temp%strucfile = strucfile            ! Copy filename
!
!write(*,*) ' ADDING A NODE ' , ASSOCIATED(store_root),' ', ASSOCIATED(store_temp) 
!write(*,*) ' ROOT FILE    ' ,store_root%strucfile(1:len_trim(store_root%strucfile))
!write(*,*) ' STRU FILE    ' ,strucfile(1:len_trim(strucfile))
      CALL store_add_node(store_root, store_temp) ! add this node to storage tree
      ALLOCATE(store_temp%crystal, STAT=istatus)  ! Allocate the crystal at this node
      IF ( istatus /= 0) THEN
         ier_num = -114
         ier_typ = ER_APPL
         ier_msg(1) = 'Could not allocate storage crystal'
         RETURN
      ENDIF
      CALL save_internal_node(store_temp, strucfile)
ELSE
!
      ALLOCATE(store_root, STAT = istatus)        ! Allocate a temporary storage
      NULLIFY(store_root%before)
      NULLIFY(store_root%after )
      IF ( istatus /= 0) THEN
         ier_num = -114
         ier_typ = ER_APPL
         ier_msg(1) = 'Temporary storage failed'
         RETURN
      ENDIF
      store_root%strucfile = strucfile            ! Copy filename
      ALLOCATE(store_root%crystal, STAT=istatus)  ! Allocate the crystal at this node
      IF ( istatus /= 0) THEN
         ier_num = -114
         ier_typ = ER_APPL
         ier_msg(1) = 'Could not allocate storage crystal'
         RETURN
      ENDIF
      CALL save_internal_node(store_root, strucfile)
ENDIF
!
      END SUBROUTINE save_internal
!
!*******************************************************************************
!
SUBROUTINE save_internal_node(ptr, strucfile)
!
      USE discus_config_mod 
      USE crystal_mod 
!     USE gen_add_mod 
      USE class_internal
      USE molecule_mod 
!     USE sym_add_mod 
      USE discus_save_mod 
      USE errlist_mod
      IMPLICIT none 
!
TYPE(internal_storage), POINTER :: ptr
CHARACTER ( LEN=* ), INTENT(IN) :: strucfile 
!
CHARACTER ( LEN=80 ) :: ier_msg_local 
!
!     Allocate sufficient space, even for all headers, and atom type, if they are omitted
!
      CALL ptr%crystal%alloc_arrays(cr_natoms,MAXSCAT, &
           mole_max_mole, mole_max_type, mole_max_atom ) ! Allocate the crystal arrays
!
!     An internal crystal has ALL headers saved, logical flags are used to indicate
!     whether they were supposed to be saved or not.
!     n_latom = UBOUND(sav_latom,1)     ! Make sure we send correct array size
      CALL ptr%crystal%set_crystal_save_flags (sav_w_scat, & 
           sav_w_adp, sav_w_occ, sav_w_gene, sav_w_symm,                &
           sav_w_ncell, sav_w_obje, sav_w_doma, sav_w_mole, sav_w_prop, &
           sav_sel_prop,MAXSCAT,sav_latom)
!
      CALL ptr%crystal%set_crystal_from_standard(strucfile) ! Copy complete crystal
!
      IF(ier_num == -157) THEN
         ier_msg_local = ier_msg(1)
         CALL store_remove_single( strucfile, ier_num)
         ier_num = -157
         ier_typ = ER_APPL
         ier_msg(1) = ier_msg_local
         ier_msg(2) = 'Atom is deseclected but part of a molecule'
         ier_msg(3) = 'Remove this atom type and purge prior to save '
      ENDIF
!
!     CALL store_write_node(store_root)
!
!     NULLIFY(store_temp)
!
END SUBROUTINE save_internal_node
!
!*******************************************************************************
!
      SUBROUTINE save_store_setting
!
      USE crystal_mod 
      USE discus_save_mod 
      USE discus_allocate_appl_mod
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: MIN_PARA = 20
      INTEGER            :: maxw
      INTEGER            :: n_nscat
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
         CALL alloc_save (  n_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      SAV_T_MAXSCAT  = SAV_MAXSCAT
      sav_t_latom(:) = sav_latom(:) 
!
      sav_t_sel_atom = sav_sel_atom
!
      sav_t_file    = sav_file
!
      sav_t_keyword = sav_keyword
!
      sav_t_w_scat  = sav_w_scat
      sav_t_w_adp   = sav_w_adp
      sav_t_w_occ   = sav_w_occ
      sav_t_r_ncell = sav_r_ncell
      sav_t_w_ncell = sav_w_ncell
      sav_t_w_gene  = sav_w_gene
      sav_t_w_symm  = sav_w_symm
      sav_t_w_mole  = sav_w_mole
      sav_t_w_obje  = sav_w_obje
      sav_t_w_doma  = sav_w_doma
      sav_t_w_prop  = sav_w_prop
!
      sav_t_start   = sav_start
      sav_t_end     = sav_end
      sav_t_ncell(:)= sav_ncell(:)
      sav_t_sel_prop(:)=  sav_sel_prop(:)
      sav_t_ncatoms = sav_ncatoms
      sav_t_size_of = sav_size_of
!
      END SUBROUTINE save_store_setting
!
!*******************************************************************************
!
      SUBROUTINE save_restore_setting
!
      USE discus_config_mod
      USE crystal_mod 
      USE discus_allocate_appl_mod
      USE discus_save_mod 
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: MIN_PARA = 20
      INTEGER            :: maxw
      INTEGER            :: n_nscat
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
         CALL alloc_save (  n_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      SAV_MAXSCAT  = SAV_T_MAXSCAT
      sav_latom(:) = sav_t_latom(:) 
!
      sav_sel_atom = sav_t_sel_atom
!
      sav_file    = sav_t_file
!
      sav_keyword = sav_t_keyword
!
      sav_w_scat  = sav_t_w_scat
      sav_w_adp   = sav_t_w_adp
      sav_w_occ   = sav_t_w_occ
      sav_r_ncell = sav_t_r_ncell
      sav_w_ncell = sav_t_w_ncell
      sav_w_gene  = sav_t_w_gene
      sav_w_symm  = sav_t_w_symm
      sav_w_mole  = sav_t_w_mole
      sav_w_obje  = sav_t_w_obje
      sav_w_doma  = sav_t_w_doma
      sav_w_prop  = sav_t_w_prop
!
      sav_start   = sav_t_start
      sav_end     = sav_t_end
      sav_ncell(:)= sav_t_ncell(:)
      sav_sel_prop(:)=  sav_t_sel_prop(:)
      sav_ncatoms = sav_t_ncatoms
      sav_size_of = sav_t_size_of
!
      END SUBROUTINE save_restore_setting
!
!*******************************************************************************
!
      SUBROUTINE save_default_setting
!
      USE discus_config_mod
      USE crystal_mod 
      USE discus_allocate_appl_mod
      USE discus_save_mod 
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: MIN_PARA = 20
      INTEGER            :: maxw
      INTEGER            :: n_nscat
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
         CALL alloc_save (  n_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      SAV_MAXSCAT  = MAXSCAT
      sav_latom(:) = .true.
!
      sav_sel_atom = .true.
!
      sav_file    = 'crystal.stru'
!
      sav_keyword = .true.
!
      sav_w_scat  = .false.
      sav_w_adp   = .false.
      sav_w_occ   = .false.
      sav_r_ncell = .true.
      sav_w_ncell = .true.
      sav_w_gene  = .true.
      sav_w_symm  = .true.
      sav_w_mole  = .true.
      sav_w_obje  = .true.
      sav_w_doma  = .true.
      sav_w_prop  = .true.
!
      sav_start   =  1
      sav_end     = -1
      sav_ncell(:)=  1
      sav_sel_prop(:)=  0
      sav_ncatoms =  1
      sav_size_of = 1
!
      END SUBROUTINE save_default_setting
!
!*******************************************************************************
!
SUBROUTINE save_reset
!
USE discus_allocate_appl_mod
USE discus_save_mod

IMPLICIT NONE
!
CALL alloc_save(1)
!
SAV_T_MAXSCAT = 1
SAV_MAXSCAT  =  1
!
IF(ALLOCATED(sav_latom))   sav_latom(:)   = .TRUE. ! (0:MAXSCAT)
IF(ALLOCATED(sav_t_latom)) sav_t_latom(:) = .TRUE. ! (0:MAXSCAT)
!
sav_t_sel_atom = .TRUE.
sav_sel_atom   = .TRUE.
!
sav_t_file    = 'crystal.stru' 
sav_file      = 'crystal.stru'
!
sav_t_keyword = .TRUE.
sav_keyword   = .TRUE.
!
sav_t_w_scat  = .FALSE.
sav_w_scat    = .FALSE.
sav_t_w_adp   = .FALSE.
sav_w_adp     = .FALSE.
sav_t_w_occ   = .FALSE.
sav_w_occ     = .FALSE.
sav_t_r_ncell = .FALSE.
sav_r_ncell   = .FALSE.
sav_t_w_ncell = .FALSE.
sav_w_ncell   = .FALSE.
sav_t_w_gene  = .TRUE.
sav_w_gene    = .TRUE.
sav_t_w_symm  = .TRUE.
sav_w_symm    = .TRUE.
sav_t_w_mole  = .TRUE.
sav_w_mole    = .TRUE.
sav_t_w_obje  = .TRUE.
sav_w_obje    = .TRUE.
sav_t_w_doma  = .TRUE.
sav_w_doma    = .TRUE.
sav_t_w_prop  = .TRUE.
sav_w_prop    = .TRUE.
!
sav_t_start   =  1
sav_start     =  1
sav_t_end     =  1
sav_end       = -1
sav_t_ncell(:)   =  1 
sav_ncell(:)     =  1
sav_t_sel_prop(:)=  (/0,0/)
sav_sel_prop(:)  =  (/0,0/)
sav_t_ncatoms =  1
sav_ncatoms   =  1
!
END SUBROUTINE save_reset
!
!*******************************************************************************
!
SUBROUTINE save_show
!
USE crystal_mod
USE discus_save_mod
USE prompt_mod 
      USE prop_char_mod
      USE prop_para_func
!
CHARACTER(LEN=32) :: c_property
INTEGER :: i
INTEGER :: length
INTEGER :: sav_flen 
!      
sav_flen = LEN_TRIM(sav_file)
IF (sav_flen.gt.0) THEN 
   WRITE (output_io, 3000) sav_file (1:sav_flen) 
ELSE 
   WRITE (output_io, 3000) ' ' 
ENDIF 
IF (sav_keyword) THEN 
   WRITE (output_io, 3005) 
   IF (sav_w_scat) THEN 
      WRITE (output_io, 3008) 'written' 
   ELSE 
      WRITE (output_io, 3008) 'omitted' 
   ENDIF 
   IF (sav_w_adp) THEN 
      WRITE (output_io, 3009) 'written' 
   ELSE 
      WRITE (output_io, 3009) 'omitted' 
   ENDIF 
   IF (sav_w_occ) THEN 
      WRITE (output_io, 4009) 'written' 
   ELSE 
      WRITE (output_io, 4009) 'omitted' 
   ENDIF 
   IF (sav_w_gene) THEN 
      WRITE (output_io, 3010) 'written' 
   ELSE 
      WRITE (output_io, 3010) 'omitted' 
   ENDIF 
   IF (sav_w_symm) THEN 
      WRITE (output_io, 3020) 'written' 
   ELSE 
      WRITE (output_io, 3020) 'omitted' 
   ENDIF 
   IF (sav_w_ncell) THEN 
      WRITE (output_io, 3030) 'written' 
   ELSE 
      WRITE (output_io, 3030) 'omitted' 
   ENDIF 
   IF (sav_w_mole) THEN 
      WRITE (output_io, 3040) 'written' 
   ELSE 
      WRITE (output_io, 3040) 'omitted' 
   ENDIF 
   IF (sav_w_obje) THEN 
      WRITE (output_io, 3050) 'written' 
   ELSE 
      WRITE (output_io, 3050) 'omitted' 
   ENDIF 
   IF (sav_w_doma) THEN 
      WRITE (output_io, 3060) 'written' 
   ELSE 
      WRITE (output_io, 3060) 'omitted' 
   ENDIF 
   CALL char_prop_2 (c_property,sav_sel_prop (1), sav_sel_prop (0),   &
      length)
   WRITE (output_io, 3131) c_property (1:length)

!                                                                       
   WRITE (output_io, 3090) 
   WRITE (output_io, 3091) 
   DO i = 0, cr_nscat 
      IF (sav_latom (i) ) THEN 
         WRITE (output_io, 3092) i, cr_at_lis (i) 
      ENDIF 
   ENDDO 
   IF (sav_end.eq. - 1) THEN 
      WRITE (output_io, 3080) 
   ELSE 
      WRITE (output_io, 3081) sav_start, sav_end 
   ENDIF 
ELSE 
   WRITE (output_io, 3006) 
ENDIF 
!                                                                       
 3000 FORMAT(20x,'    Data written to structure file '/                 &
     &       20x,' ===================================='//              &
     &       ' Structure file                            : ',a/)        
 3005 FORMAT(' Format of structure file                  :',            &
     &       ' Keyword controlled'/                                     &
     &       ' Title           ',26x,': written'/                       &
     &       ' Space group     ',26x,': written'/                       &
     &       ' Cell constants  ',26x,': written')                       
 3006 FORMAT(' Format of structure file                  :',            &
     &       ' No keywords in structure file'/                          &
     &       ' First line      ',26x,': title'/                         &
     &       ' Second line     ',26x,': space group'/                   &
     &       ' Third line      ',26x,': cell constants'/                &
     &       ' Following lines ',26x,': atoms'/)                        
 3008 FORMAT(' Explicit list of atomic names             : ',a7) 
 3009 FORMAT(' Explicit list of atomic displacement par. : ',a7) 
 4009 FORMAT(' Explicit list of occupancy parameters     : ',a7) 
 3010 FORMAT(' Additional generator matrizes             : ',a7) 
 3020 FORMAT(' Additional symmetry operators             : ',a7) 
 3030 FORMAT(' Number of unit cells, atoms per unit cell : ',a7) 
 3040 FORMAT(' Molecule information: content etc.        : ',a7) 
 3050 FORMAT(' Object   information: content etc.        : ',a7) 
 3060 FORMAT(' Domain   information: content etc.        : ',a7) 
 3080 FORMAT(' Range of atoms from to    : All atoms included') 
 3081 FORMAT(' Range of atoms from to    : ',2(2x,i9)) 
 3131 FORMAT    (/' Atom properties         : ','NMDOEI'/               &
                  '      absent=- ignored=. : ',a)

 3090 FORMAT(' Selected atoms    :') 
 3091 FORMAT('                      type name') 
 3092 FORMAT(20x,2(4x,i2,1x,a4)) 
END SUBROUTINE save_show
!
!*******************************************************************************
!
!
END MODULE save_menu
