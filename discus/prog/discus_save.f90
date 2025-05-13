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
USE lib_errlist_func
USE lib_help
USE lib_do_operating_mod
USE lib_echo
USE lib_length
USE lib_macro_func
      USE class_macro_internal 
!
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      USE sup_mod
use take_param_mod
      IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: string 
INTEGER         , INTENT(INOUT) :: lcomm
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA =  20 ! A command requires at least these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(string))), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
      REAL(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
!
      CHARACTER ( LEN=MAX(PREC_STRING,LEN(string)) ) :: zeile 
      CHARACTER(len=9) :: befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) :: line
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz
integer :: is_style=0    ! style for output file format
      INTEGER sav_flen 
      LOGICAL lend 
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_FORMAT  = 1
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'format' /
data loname /  6       /
!
DATA sav_flen / 1 / 
!
opara  =  (/ 'keyword' /)   ! Always provide fresh default values
lopara =  (/  7        /)
owerte =  (/  0.0      /)
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
loop_main: DO while (.not.lend) 
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
   CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
!
   if(ier_num/=0) cycle loop_main
   if (line == ' '      .or. line(1:1) == '#' .or. &
       line == char(13) .or. line(1:1) == '!'        ) cycle loop_main
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=') 
   cond_equal: IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
      cycle loop_main
   endif cond_equal
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
   cond_general: IF (befehl (1:1) .eq.'@') THEN 
                  IF (length.ge.2) THEN 
                     line(1:length-1) = line(2:length)
                     line(length:length) = ' '
                     length = length - 1
                     CALL file_kdo(line, length)
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
      cycle loop_main
!                                                                       
!     ----list asymmetric unit 'asym'                                   
!                                                                       
   ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) THEN cond_general
                  CALL show_asym 
      cycle loop_main
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
   ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) THEN cond_general
                  CALL macro_continue (zeile, lp) 
      cycle loop_main
!                                                                       
!     ----list atoms present in the crystal 'chem'                      
!                                                                       
   ELSEIF (str_comp (befehl, 'chemistry', 2, lbef, 9) ) THEN cond_general
                  CALL show_chem 
      cycle loop_main
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
   ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN cond_general
                  CALL echo (zeile, lp) 
      cycle loop_main
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
   ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) THEN cond_general
                  CALL do_eval (zeile, lp, .TRUE.) 
      cycle loop_main
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .true. 
      cycle loop_main
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
   ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) THEN                                      cond_general
                  line = zeile 
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                     lp = lp + 7 
                     CALL do_hel ('discus '//line, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus save '//line, lp) 
                  ENDIF 
      cycle loop_main
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
   ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN cond_general
                  IF (zeile.ne.' ') THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
      cycle loop_main
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
   ELSEIF (str_comp (befehl, 'wait', 2, lbef, 4) ) THEN cond_general
                  CALL do_input (zeile, lp) 
      cycle loop_main
!
!     ----Reset save menu 'rese'n'                                      
!
   ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) THEN cond_general
      CALL save_reset
      is_style = 0
      cycle loop_main
   endif cond_general
!  ELSE cond_general
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
               IF (str_comp (befehl, 'deselect', 1, lbef, 8) ) THEN 
                   CALL atom_select (zeile, lp, 0, SAV_MAXSCAT, sav_latom, &
                   sav_lsite, 0, SAV_MAXSITE,                              &
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
      if(ier_num/=0) cycle loop_main
      call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
           oname, loname, opara, lopara, lpresent, owerte)
      if(ier_num/=0) cycle loop_main
!
      if(lpresent(O_FORMAT)) then             ! User provided format:
         if(str_comp(opara(O_FORMAT), 'keyword', 1, lpara(O_FORMAT), 7)) then
            sav_keyword = .true.
            is_style = 0
         elseif(str_comp(opara(O_FORMAT), 'nexus', 3, lopara(O_FORMAT), 5)) then
            sav_keyword = .true.
            is_style = 1
         elseif(str_comp(opara(O_FORMAT), 'nokeyword', 3, lpara(O_FORMAT), 9)) then
            sav_keyword = .false.
            is_style = 2
         endif
      else
         IF(str_comp(cpara(1) , 'keyword', 1, lpara(1), 7) ) THEN                                          
            sav_keyword = .true. 
            is_style = 0
         ELSEIF (str_comp (cpara (1) , 'nexus', 3, lpara(1), 5) ) THEN
            sav_keyword = .true. 
            is_style = 1
         ELSEIF (str_comp (cpara (1) , 'nokeyword', 3, lpara(1), 9) ) THEN
            sav_keyword = .false. 
            is_style = 2
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      endif
!     ENDIF 
!                                                                       
!     ----Select range of atoms within crystal to be included 'incl'    
!                                                                       
   ELSEIF (str_comp (befehl, 'include', 1, lbef, 7) ) THEN 
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
                        sav_w_surf = .false. 
                        sav_w_magn = .false. 
                        sav_w_obje = .false. 
                        sav_w_mole = .false. 
                        sav_w_doma = .false. 
                        sav_w_prop = .false. 
                     ELSEIF(str_comp(cpara(1), 'generator', 1, lpara(1), 9) ) THEN
                        sav_w_gene = .false. 
                     ELSEIF(str_comp(cpara(1), 'molecules', 1, lpara(1), 9) ) THEN
                        sav_w_mole = .false. 
                     ELSEIF(str_comp(cpara(1), 'objects', 1, lpara(1), 7) ) THEN
                        sav_w_obje = .false. 
                     ELSEIF(str_comp(cpara(1), 'domains', 1, lpara(1), 7) ) THEN
                        sav_w_doma = .false. 
                     ELSEIF(str_comp(cpara(1), 'ncell', 1, lpara(1), 5) ) THEN
                        sav_w_ncell = .false. 
                     ELSEIF(str_comp(cpara(1), 'symmetry', 2, lpara(1), 8) ) THEN
                        sav_w_symm = .false. 
                     ELSEIF(str_comp(cpara(1), 'scat', 2, lpara(1), 4) ) THEN
                        sav_w_scat = .false. 
                     ELSEIF(str_comp(cpara(1), 'adp', 1, lpara(1), 3) ) THEN
                        sav_w_adp = .false. 
                     ELSEIF(str_comp(cpara(1), 'occ', 1, lpara(1), 3) ) THEN
                        sav_w_occ = .false. 
                     ELSEIF(str_comp(cpara(1), 'surface', 1, lpara(1), 7) ) THEN
                        sav_w_surf = .false. 
                     ELSEIF(str_comp(cpara(1), 'magnetic', 1, lpara(1), 8) ) THEN
                        sav_w_magn = .false. 
                     ELSEIF(str_comp(cpara(1), 'property', 1, lpara(1), 8) ) THEN
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
            IF(str_comp(sav_file(1:8),'internal',8,8,8)) THEN
               CALL save_internal (sav_file) 
            elseif(is_style==0) then
               CALL save_keyword (sav_file) 
            elseif(is_style==1) then
               CALL save_keyword_nexus(sav_file) 
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
               ELSEIF (str_comp (befehl, 'select', 2, lbef, 6) ) THEN 
                   CALL atom_select (zeile, lp, 0, SAV_MAXSCAT, sav_latom, &
                   sav_lsite, 0, SAV_MAXSITE,                              &
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
                        sav_w_surf = .true. 
                        sav_w_magn = .FALSE.   ! MAGNETIC_WORK
                        sav_w_mole = .true. 
                        sav_w_obje = .true. 
                        sav_w_doma = .true. 
                        sav_w_prop = .true. 
                     ELSEIF(str_comp(cpara (1), 'generator', 1, lpara(1), 9) ) THEN
                        sav_w_gene = .true. 
                     ELSEIF(str_comp(cpara (1), 'molecules', 1, lpara(1), 9) ) THEN
                        sav_w_mole = .true. 
                     ELSEIF(str_comp(cpara (1), 'objects',   1, lpara(1), 7) ) THEN
                        sav_w_obje = .true. 
                     ELSEIF(str_comp(cpara (1), 'domains',   1, lpara(1), 7) ) THEN
                        sav_w_doma = .true. 
                     ELSEIF(str_comp(cpara (1), 'ncell',     1, lpara(1), 5) ) THEN
                        sav_w_ncell = .true. 
                     ELSEIF(str_comp(cpara (1), 'symmetry',  2, lpara(1), 8) ) THEN
                        sav_w_symm = .true. 
                     ELSEIF(str_comp(cpara (1), 'scat',      2, lpara(1), 4) ) THEN
                        sav_w_scat = .true. 
                     ELSEIF(str_comp(cpara (1), 'adp',       1, lpara(1), 3) ) THEN
                        sav_w_adp = .true. 
                     ELSEIF(str_comp(cpara (1), 'occ',       1, lpara(1), 3) ) THEN
                        sav_w_occ = .true. 
                     ELSEIF(str_comp(cpara (1), 'surface',   1, lpara(1), 7) ) THEN
                        sav_w_surf = .true. 
                     ELSEIF(str_comp(cpara (1), 'magnetic',  1, lpara(1), 8) ) THEN
                        sav_w_magn = .FALSE.    ! MAGNETIC_WORK
                     ELSEIF(str_comp(cpara (1), 'property',  1, lpara(1), 8) ) THEN
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
!  ENDIF  cond_general
!  ENDIF cond_equal
   IF(linteractive .OR. lmakro) THEN
      CYCLE loop_main
   ELSE
      EXIT loop_main
   ENDIF
ENDDO  loop_main
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
      INTEGER                      :: n_nscat = 1
      INTEGER                      :: n_nsite = 1
!
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
         n_nsite = MAX(cr_ncatoms, SAV_MAXSCAT, MAXSCAT) 
         CALL alloc_save (  n_nscat, n_nsite )
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
USE lib_length
USE support_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER ( * ) strucfile 
      INTEGER ist, i, j 
      LOGICAL lread 
!                                                                       
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
         WRITE (ist, 4) cr_at_lis (cr_iscat (1,i) ), (cr_pos (j, i),      &
         j = 1, 3), cr_dw (cr_iscat (1,i) )                               
         ENDDO 
      ENDIF 
      CLOSE (ist) 
!                                                                       
    3 FORMAT (a) 
   33 FORMAT  (a16) 
   34 FORMAT  (6(1x,f10.6)) 
    4 FORMAT (a4,3(2x,f14.6),5x,f8.4) 
      END SUBROUTINE save_nokeyword                 
!
!********************************************************************** 
!
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
USE surface_mod
!
use allocate_generic
USE errlist_mod 
USE lib_length
USE support_mod
!
IMPLICIT none 
!
CHARACTER(LEN=*), INTENT(IN) :: strucfile 
!
real(kind=PREC_DP), parameter :: TOL = 2.0D-5
INTEGER, PARAMETER :: ist = 67
CHARACTER(LEN=15), DIMENSION(-4:4) :: C_MOLE !( - 4:4) 
CHARACTER(LEN=1), DIMENSION(0:SURF_MAXTYPE) :: c_surf
character(len=PREC_STRING) :: scat_string
INTEGER :: i, j, k
INTEGER :: i_start, i_end 
integer :: iscat
integer :: ianis
integer :: itype
INTEGER :: is, ie 
INTEGER :: wr_prop = 1
INTEGER :: wr_mole = 0
INTEGER :: wr_cont = 0
INTEGER :: wr_doma_current = 0               ! Current domain number
LOGICAL :: active_domain    =.FALSE.         ! Are we writing a domain
INTEGER, DIMENSION(0:3) :: wr_surf
REAL(KIND=PREC_DP)   , DIMENSION(0:3) :: wr_magn
LOGICAL lread 
LOGICAL                            :: lsave
LOGICAL, DIMENSION(:), ALLOCATABLE :: lwrite ! flag if atom needs write
integer, dimension(:,:), allocatable :: look_anis   ! atoms with unit cell have this ANIS ADP
!                                                                       
!                                                                       
!                                                                       
DATA C_MOLE / 'domain_fuzzy   ', 'domain_sphere  ', 'domain_cylinder',&
              'domain_cube    ', 'atoms          ', 'cube           ',&
                    'cylinder       ', 'sphere         ', 'cube           ' /             
!
DATA c_surf(0:SURF_MAXTYPE) /'_','P', 'S', 'Y', 'E', 'C', 'L', 'T'/
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
   WRITE (ist, 3010) cr_spcgr, cr_set
ELSEIF (spcgr_ianz.eq.1) THEN 
   WRITE (ist, 3011) cr_spcgr, spcgr_para, cr_set
ENDIF 
WRITE (ist, 3020) (cr_a0 (i), i = 1, 3), (cr_win (i), i = 1, 3) 
!
!   Construct ADP lookup table
!
allocate(look_anis(1:cr_nscat, 0:max(cr_ncatoms,cr_nanis)))    ! For all atom types , all ADPs
look_anis = 0
loop_makelook: do j=1, cr_natoms                  ! Find all adps for atom types
   iscat = cr_iscat(1, j)
   ianis = cr_iscat(3, j)
   if(iscat==0) cycle loop_makelook            ! Ignore VOID
   do i=1, look_anis(iscat,0)                  ! Loop over all previous entries
      if(look_anis(iscat,i)==ianis) cycle loop_makelook
   enddo
   look_anis(iscat,0) = look_anis(iscat,0) + 1
   if(look_anis(iscat,0)> ubound(look_anis,2)) then
      k = ubound(look_anis,2) + 10
      call alloc_arr(look_anis, 1, cr_nscat, 0, k, i, 0  )
   endif
   look_anis(iscat,look_anis(iscat,0)) = ianis
enddo loop_makelook
!
!   Write Atom names 
!
IF (sav_w_scat) THEN 
   scat_string    = 'scat  '                   ! Construct "SCAT" line
   itype = 0
   do j=1, cr_nscat
      do i=1,look_anis(j,0)
         itype = itype + 1
         ianis = look_anis(j,i)
         is = 7 + (itype-1)*7
         ie = is + 6
         scat_string(is:ie) = cr_at_lis(j)//',  '
         if(itype==7) then
            write(ist, '(a)') scat_string(1:len_trim(scat_string)-1)
            scat_string    = 'scat '
            itype = 0
         endif
      enddo
   enddo
   if(itype> 0) then
      write(ist, '(a)') scat_string(1:len_trim(scat_string)-1)
   endif
!  j = (cr_nscat - 1) / 7 
!  DO i = 1, j 
!     is = (i - 1) * 7 + 1 
!     ie = is + 6 
!     WRITE (ist, 3110) (cr_at_lis (k), k = is, ie) 
!  ENDDO 
!  IF (cr_nscat - j * 7.eq.1) THEN 
!     WRITE (ist, 3111) cr_at_lis (cr_nscat) 
!  ELSEIF (cr_nscat - j * 7.gt.1) THEN 
!     WRITE (fform, 7010) cr_nscat - j * 7 - 1 
!     WRITE (ist, fform) (cr_at_lis (i), i = j * 7 + 1, cr_nscat) 
!  ENDIF 
ENDIF 
!
IF (sav_w_adp) THEN 
   scat_string    = 'adp   '                   ! Construct "ADP" line
   itype = 0
   do j=1, cr_nscat
      do i=1,look_anis(j,0)
         itype = itype + 1
         ianis = look_anis(j,i)
         is = 6 + (itype-1)*15
         ie = is + 15
         write(scat_string(is:ie),'(f10.6,a1)') cr_dw(j),','
         if(itype==7) then
            write(ist, '(a)') scat_string(1:len_trim(scat_string)-1)
            scat_string    = 'adp  '
            itype = 0
         endif
      enddo
   enddo
   if(itype> 0) then
      write(ist, '(a)') scat_string(1:len_trim(scat_string)-1)
   endif
!
!  j = (cr_nscat - 1) / 7 
!  DO i = 1, j 
!     is = (i - 1) * 7 + 1 
!     ie = is + 6 
!     WRITE (ist, 3120) (cr_dw (k), k = is, ie) 
!  ENDDO 
!  IF (cr_nscat - j * 7.eq.1) THEN 
!     WRITE (ist, 3121) cr_dw (cr_nscat) 
!  ELSEIF (cr_nscat - j * 7.gt.1) THEN 
!     WRITE (fform, 7020) cr_nscat - j * 7 - 1 
!     WRITE (ist, fform) (cr_dw (i), i = j * 7 + 1, cr_nscat) 
!  ENDIF 
!  allocate(look_anis(cr_ncatoms))
!  look_anis = 1
!  do j=1, cr_ncatoms    ! loop over atoms in a unit cell
!     look_anis(cr_iscat(1,j)) = cr_iscat(3,j)
!  enddo
!  do j= 1, cr_nanis
   itype = 0
   do j=1, cr_nscat
      do i=1,look_anis(j,0)
         itype = itype + 1
         ianis = look_anis(j,i)
!     if(abs(cr_prin(4,1,cr_iscat(3,j))-cr_prin(4,2,cr_iscat(3,j)))>TOL  .or.  &
!        abs(cr_prin(4,1,cr_iscat(3,j))-cr_prin(4,3,cr_iscat(3,j)))>TOL      ) then
!   do j= 1, cr_ncatoms
!      if(abs(cr_prin(4,1,look_anis(j,ianis) )-cr_prin(4,2,look_anis(j,ianis) ))>TOL  .or.  &
!         abs(cr_prin(4,1,look_anis(j,ianis) )-cr_prin(4,3,look_anis(j,ianis) ))>TOL      ) then
!write(*,*)  abs(cr_prin(4,1,ianis              )-cr_prin(4,2,ianis              ))>TOL , &
!           abs(cr_prin(4,1,ianis              )-cr_prin(4,3,ianis              ))>TOL  , &
!           maxval(abs(cr_prin(:,:,ianis)))<TOL                                          

       if(abs(cr_prin(4,1,ianis              )-cr_prin(4,2,ianis              ))>TOL  .or.  &
          abs(cr_prin(4,1,ianis              )-cr_prin(4,3,ianis              ))>TOL  .or.  &
          maxval(abs(cr_prin(:,:,ianis)))<TOL                                             ) then
          write(ist, '(a,i2.2,a,5(f10.6,'',''),f10.6,a)') 'anis type:', itype, ', values:[', cr_anis_full(:,ianis),']'
       else
          write(ist, '(a,i2.2,a,f10.6,a)')  'anis type:', itype, ', values:[', cr_prin(4,1,    ianis         ),']'
       endif
      enddo
   enddo
ENDIF 
!
IF (sav_w_occ) THEN 
   scat_string    = 'occ   '                   ! Construct "ADP" line
   itype = 0
   do j=1, cr_nscat
      do i=1,look_anis(j,0)
         itype = itype + 1
         ianis = look_anis(j,i)
         is = 6 + (itype-1)*15
         ie = is + 15
         write(scat_string(is:ie),'(f10.6,a1)') cr_occ(j),','
         if(itype==7) then
            write(ist, '(a)') scat_string(1:len_trim(scat_string)-1)
            scat_string    = 'occ  '
            itype = 0
         endif
      enddo
   enddo
   if(itype> 0) then
      write(ist, '(a)') scat_string(1:len_trim(scat_string)-1)
   endif
!  j = (cr_nscat - 1) / 7 
!  DO i = 1, j 
!     is = (i - 1) * 7 + 1 
!     ie = is + 6 
!     WRITE (ist, 3220) (cr_occ(k), k = is, ie) 
!  ENDDO 
!  IF (cr_nscat - j * 7.eq.1) THEN 
!     WRITE (ist, 3221) cr_occ(cr_nscat) 
!  ELSEIF (cr_nscat - j * 7.gt.1) THEN 
!     WRITE (fform, 7030) cr_nscat - j * 7 - 1 
!     WRITE (ist, fform) (cr_occ(i), i = j * 7 + 1, cr_nscat) 
!  ENDIF 
ENDIF 
!
IF (sav_w_gene) THEN 
   DO k = 1, gen_add_n 
      WRITE (ist, 3021) ( (gen_add (i, j, k), j = 1, 4), i = 1, 3),  &
      gen_add_power (k)                                              
   ENDDO 
ENDIF 
!
IF (sav_w_symm) THEN 
   DO k = 1, sym_add_n 
      WRITE (ist, 3022) ( (sym_add (i, j, k), j = 1, 4), i = 1, 3),  &
      sym_add_power (k)                                              
   ENDDO 
ENDIF 
!
IF (sav_w_ncell) THEN 
   WRITE (ist, 3030) cr_icc, cr_ncatoms , cr_natoms, cr_nscat, &
   mole_num_mole, mole_num_type, mole_num_atom
ENDIF 
!
write(ist, 3800)
WRITE(ist, 3900) 
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
         WRITE (ist, 4900) 'object' 
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
         WRITE (ist, 4900) 'molecule' 
      ENDIF 
   ENDDO 
ENDIF 
!                                                                       
!     Write atoms                                                       
!                                                                       
loop_atoms: DO i = i_start, i_end 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
!     IF (sav_latom (cr_iscat (1,i) ) ) THEN 
   IF (check_select_status (i, sav_latom (cr_iscat (1,i) ), cr_prop (i),   &
                               sav_sel_prop)     ) THEN
!        IF(lwrite(i)) THEN
      wr_prop = 1
      wr_mole = 0
      wr_cont = 0
      wr_surf(:) = 0
      wr_magn(:) = 0.0
      IF(sav_w_prop) wr_prop = cr_prop(i)
      IF (sav_w_mole .OR. sav_w_doma .OR. sav_w_obje) THEN 
         IF(active_domain) THEN
            IF(cr_mole(i)/=wr_doma_current) THEN   ! DOMAIN is finished
               WRITE(ist, '(''domain end'')')
               active_domain = .FALSE.
               wr_doma_current = 0
            ENDIF
         ENDIF
         IF(cr_mole(i)/=0) THEN
            wr_mole = cr_mole(i)
            check_mole: DO j = 1, mole_len (cr_mole(i))
               IF(mole_cont (mole_off(cr_mole(i))+j) == i) THEN
                  wr_cont = j
                  EXIT check_mole
               ENDIF
            ENDDO check_mole
            IF(mole_char(wr_mole)<MOLE_ATOM .AND. .NOT.active_domain .AND. &
               wr_mole/=wr_doma_current                                   ) THEN  ! Need header
               WRITE(ist, 4000) 'domain' 
               WRITE(ist, 4002) 'domain', mole_type (wr_mole) 
               WRITE(ist, 4100) 'domain', c_mole (mole_char (wr_mole) ) 
               k = len_str (mole_file (wr_mole) ) 
               IF(k.gt.0) THEN 
                  WRITE(ist, 4300) 'domain', mole_file(wr_mole)(1:k) 
               ELSE 
                  WRITE(ist, 4300) 'domain' 
               ENDIF 
               WRITE(ist, 4400) 'domain', mole_fuzzy (wr_mole) 
               active_domain = .TRUE.
               wr_doma_current = wr_mole
            ENDIF
         ENDIF
      ENDIF
      IF(sav_w_surf) wr_surf(0:3) = cr_surf(0:3,i)
      IF(sav_w_magn) wr_surf(0:3) = nint(cr_magn(0:3,i))   ! MAGNETIC_WORK
      WRITE (ist, 4) cr_at_lis (cr_iscat (1,i) ),         &
                     (cr_pos (j, i),j = 1, 3),          &
                     cr_dw (cr_iscat (1,i) ), wr_prop,    &
                     wr_mole, wr_cont, cr_occ(cr_iscat(1,i)), &
                     c_surf(wr_surf(0  )), wr_surf(1:3  )!, &
!                    cr_anis_full(:,i)
   ENDIF 
ENDDO loop_atoms
!
IF(active_domain) THEN
   IF(cr_mole(i)/=wr_doma_current) THEN   ! DOMAIN is finished
      WRITE(ist, '(''domain end'')')
      active_domain = .FALSE.
      wr_doma_current = 0
   ENDIF
ENDIF
!
CLOSE(ist) 
DEALLOCATE(lwrite)
!                                                                       
 3000 FORMAT    ('title ',a) 
 3010 FORMAT    ('spcgr ',a16, ', setting:',a3)
 3011 FORMAT    ('spcgr ',a16,',',i4, ', setting:',a3) 
 3020 FORMAT    ('cell  ',5(f10.6,','),f10.6) 
 3021 FORMAT    ('gene  ',12(f9.6,','),i3) 
 3022 FORMAT    ('symm  ',12(f9.6,','),i3) 
 3030 FORMAT    ('ncell ',3(i8,','),i10,5(',',i12)) 
!3110 FORMAT    (('scat  ', a4,6(',',5x,a4))) 
!3111 FORMAT    (('scat  ', a4)) 
!3120 FORMAT    (('adp   ', f9.6,6(',',5x,f9.6))) 
!3121 FORMAT    (('adp   ', f9.6)) 
!3220 FORMAT    (('occ   ', f9.6,6(',',5x,f9.6))) 
!3221 FORMAT    (('occ   ', f9.6)) 
 3800 format    ('format numbers,XYZBPMMOS')
!3900 FORMAT    ('atoms      x,',14x,'y,',14x,'z,',13x,'Biso,', 4x,'Property,', &
!                2x,'MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl',            &
!                ',  U11,       U22,       U33,       U23,       U13,       U12')
 3900 FORMAT    ('atoms      x,',14x,'y,',14x,'z,',13x,'Biso,', 4x,'Property,', &
                 2x,'MoleNo,  MoleAt,   Occ,     St,  Sh,  Sk,  Sl')
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
    4 FORMAT (a4,3(1x,f14.6,','),4x,f10.6,',',i8, ',', I8, ',', &
              I8,', ', F10.6,', ',A1,3(', ',I3):,6(',',f10.6)) 
!7010 FORMAT    ('(''scat  '', a4  ,',i1,'('','',5x,a4  ))') 
!7020 FORMAT    ('(''adp   '', f9.6,',i1,'('','',5x,f9.6))') 
!7030 FORMAT    ('(''occ   '', f9.6,',i1,'('','',5x,f9.6))') 
!
END SUBROUTINE save_keyword                   
!
!********************************************************************** 
!
subroutine save_keyword_nexus (outfile) 
!-
! Interfact to the generiv NeXuS nx_write_structure
!+
!
use celltoindex_mod
use chem_mod
use chem_aver_mod
use crystal_mod
use guess_atoms_mod
use molecule_mod
use wyckoff_mod
!
use element_data_mod
use envir_mod
use errlist_mod
use nx_write_mod
use lib_nx_transfer_mod
use param_mod
use wink_mod
!
implicit none
!
character(len=*), intent(in) :: outfile
!
!character(len=PREC_STRING) :: script_path
character(len=4) :: element
character(len=PREC_string) :: program_version
integer :: i, j, k, l, ia      ! Dummy indices
integer :: ios    ! I/O status
logical :: lout = .false.
logical :: lsite = .true.
!
integer, dimension(3,3) :: unit_cells
real(kind=PREC_DP), dimension(:,:,:), allocatable :: symmetry_mat
character(len=4), dimension(:), allocatable :: types_names
integer, dimension(:), allocatable :: types_ordinal
integer, dimension(:), allocatable :: types_charge
integer, dimension(:), allocatable :: types_isotope
real(kind=PREC_DP), dimension(:), allocatable :: types_occupancy
!
real(kind=PREC_DP), dimension(:,:), allocatable :: atom_pos
integer, dimension(:)  , allocatable :: atom_id
integer, dimension(:)  , allocatable :: atom_type
integer, dimension(:,:), allocatable :: atom_unit_cell
integer, dimension(:)  , allocatable :: atom_site
!
type(anis_adp_type) :: anis_adp
type(molecule_data) :: molecules
type(average_structure) :: average_struc
logical, dimension(6)  :: status_flags
!
program_version = 'DISCUS ' // version_discus(1:len_trim(version_discus))
unit_cells = 0
unit_cells(1,1) = cr_icc(1)
unit_cells(2,2) = cr_icc(2) 
unit_cells(3,3) = cr_icc(3)
allocate(symmetry_mat(3,4,spc_n))
symmetry_mat = spc_mat(1:3, 1:4, 1:spc_n)
!
call guess_atom_all
!
allocate(types_names    (1:cr_nscat))
allocate(types_ordinal  (1:cr_nscat))
allocate(types_charge   (1:cr_nscat))
allocate(types_isotope  (1:cr_nscat))
allocate(types_occupancy(1:cr_nscat))
types_names    = ' '
types_ordinal  = 0
types_charge   = 0
types_isotope  = 0
types_occupancy= 0.0_PREC_DP
!
do i=1, cr_nscat
   types_names(i) = cr_at_lis(i)
   if(cr_scat_equ(i)) then
     element = cr_at_equ(i)
   else
      element = cr_at_lis(i)
   endif
   types_ordinal(i) = get_ordi(element)
   j = len_trim(element)
   if(element(j:j)=='+') then
      read(element(j-1:j-1),'(i1)',iostat=ios) types_charge(i)
   elseif(element(j:j)=='-') then
      read(element(j-1:j-1),'(i1)',iostat=ios) types_charge(i)
      types_charge(i) = -abs(types_charge(i))
   endif
!  types_typeID(i)  = i
   types_occupancy(i) = cr_occ(i)
enddo
!
allocate(atom_pos(3,cr_natoms))
atom_pos = cr_pos(1:3, 1:cr_natoms)
!
allocate(atom_id(1:cr_natoms))
allocate(atom_type(1:cr_natoms))
allocate(atom_unit_cell(3,1:cr_natoms))
allocate(atom_site(1:cr_natoms))
!
do i=1, cr_natoms
   atom_id(i) = i
   atom_type(i) = cr_iscat(1,i)
   call indextocell(i, atom_unit_cell(:,i), atom_site(i))
enddo
!
!  Define status flags
!
status_flags = .false.
if(any(cr_icc>1)   ) status_flags(1) = .true.   ! Superstructure, more than one unit cell
if(any(cr_is_sym>1)) status_flags(2) = .false.  ! Crystal symmetry was likely applied
status_flags(3:5) = chem_period                 ! Copy periodic boundary conditions
status_flags(6)   = cr_is_homo                  ! Copy homogeneity status
!
! Build anisotropic atom structure
!
anis_adp%anis_n_type = cr_nanis
anis_adp%anis_n_atom = cr_natoms
allocate(anis_adp%anis_adp(7, cr_nanis))
anis_adp%anis_adp(1:6,1:cr_nanis) = cr_anis_full(1:6,1:cr_nanis)
do i=1, cr_nanis
   anis_adp%anis_adp(7,i) = (cr_prin(4,   1,i) + cr_prin(4,   2,i) + cr_prin(4,   3,i))/3.0_PREC_DP
enddo
allocate(anis_adp%atom_index(1:cr_natoms))
do i=1, cr_natoms
   anis_adp%atom_index(i) = cr_iscat(3,i)
enddo
!
!  Build molecule structure
!
if(mole_num_mole>0) then    ! Crystal has molecules
!rite(*,*) ' Molecules ', mole_num_mole, mole_num_type, mole_num_atom
!write(*,*) ' Molecules type ', ubound(mole_type)
!write(*,*) ' Molecules char ', ubound(mole_char)
!write(*,*) ' Molecules clin ', ubound(mole_clin)
!write(*,*) ' Molecules cqua ', ubound(mole_cqua)
!write(*,*) ' Molecules biso ', ubound(mole_biso)
   allocate(molecules%mole_int (3, mole_num_mole))
   allocate(molecules%mole_real(3, mole_num_type))
!  allocate(mole_data%mole_char(mole_num_char))
!  allocate(mole_data%mole_clin(mole_num_clin))
!  allocate(mole_data%mole_cqua(mole_num_cqua))
!  allocate(mole_data%mole_ueqv(mole_num_ueqv))
   molecules%number_moles = mole_num_mole
   molecules%number_types = mole_num_type
   molecules%mole_int(1, :) = mole_type(1:mole_num_mole)
   molecules%mole_int(2, :) = mole_char(1:mole_num_mole)
   molecules%mole_int(3, :) = mole_len (1:mole_num_mole)
   molecules%mole_real(1,:) = mole_biso(1:mole_num_type)/8.0_PREC_DP/PI**2
   molecules%mole_real(2,:) = mole_clin(1:mole_num_type)
   molecules%mole_real(3,:) = mole_cqua(1:mole_num_type)
   allocate(molecules%atom_index(1:maxval(mole_len(1:mole_num_mole)),mole_num_mole))
!write(*,*) ' MOLECULES ', ubound(molecules%atom_index), ' >> ', maxval(mole_len(1:mole_num_mole)), mole_num_mole 
!j=1
!write(*,*) ' CONTEN, k)T   ', ubound(mole_cont), ' >> ', mole_cont(mole_off(j)+1: mole_off(j)+mole_len(j))
   do j=1, mole_num_mole
      i = mole_len(j)
!write(*,*) 'MOLE No. ', j, i, mole_len(j), mole_off(j), mole_off(j)+mole_len(j)
      molecules%atom_index(1:i,j) = mole_cont(mole_off(j)+1: mole_off(j)+mole_len(j))
   enddo
else
   molecules%number_moles = 0
endif
!
average_struc%aver_n_atoms = 0  ! Assume no average structure
if(all(chem_period) .or. 1==1) then       ! Crystal has periodic boundary conditions
   lout  = .false.
   lsite = .false.
   call chem_aver(lout, lsite)
   k = 0
   do i = 1, cr_ncatoms
      do j = 1, chem_ave_n(i)
         k = k + 1
      enddo
   enddo
   allocate(average_struc%atom_type(k))
   allocate(average_struc%position(3,k))
   allocate(average_struc%occupancy(k))
   allocate(average_struc%anis_adp(7,k))
   allocate(average_struc%site_number(k))
!
   ia = cr_icc(1) * cr_icc(2) * cr_icc(3)
   k = 0
   do i = 1, cr_ncatoms
      do j = 1, chem_ave_n(i)
         k = k + 1
         average_struc%atom_type(k)   = chem_ave_iscat(i, j)
         average_struc%position(:,k)  = chem_ave_posit(:, i, j)
         average_struc%occupancy(k)   = chem_ave_bese (i, j)/real(ia, kind=PREC_DP)
         l                            = chem_ave_anis(i,j)        ! ADP type
         average_struc%anis_adp(1:6,k) = cr_anis_full(1:6,l)
         average_struc%anis_adp(7  ,k) = (cr_prin(4, 1,l) + cr_prin(4, 2,l) + cr_prin(4, 3,l))/3.0_PREC_DP
         average_struc%site_number(k) = i
      enddo
   enddo
   average_struc%aver_n_atoms = k  ! Set number of atoms in average structure
else
   average_struc%aver_n_atoms = 0  ! Assume no average structure
endif
!
call nx_write_structure(python_script_dir, outfile, program_version, author,          &
           cr_a0, cr_win, cr_gten, cr_spcgr, spcgr_para, cr_set,                &
           spc_n, symmetry_mat, unit_cells,                                     &
           cr_nscat, types_names, types_ordinal, types_charge, types_isotope,   &
           cr_natoms, atom_id, atom_type, atom_pos, atom_unit_cell, atom_site,  &
           status_flags, cr_nanis, ier_num,                                     &
           property_flags = cr_prop,                                            &
           anis_adp       = anis_adp,                                           &
           molecules      = molecules,                                          &
           types_occupancy= types_occupancy,                                    &
           average_struc  = average_struc                                       &
           )
!
deallocate(symmetry_mat)
deallocate(types_names)
deallocate(types_ordinal)
deallocate(types_charge)
deallocate(types_isotope)
deallocate(types_occupancy)
deallocate(atom_id)
deallocate(atom_type)
deallocate(atom_pos )
deallocate(atom_unit_cell)
deallocate(atom_site)
deallocate(anis_adp%atom_index)
deallocate(anis_adp%anis_adp)
if(allocated(average_struc%atom_type)) deallocate(average_struc%atom_type)
if(allocated(average_struc%position)) deallocate(average_struc%position)
if(allocated(average_struc%occupancy)) deallocate(average_struc%occupancy)
if(allocated(average_struc%anis_adp)) deallocate(average_struc%anis_adp)
if(allocated(average_struc%site_number)) deallocate(average_struc%site_number)
!
end subroutine save_keyword_nexus
!
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
subroutine save_internal_node(ptr, strucfile)
!
USE discus_config_mod 
USE crystal_mod 
USE class_internal
USE molecule_mod 
USE discus_save_mod 
USE errlist_mod
!
IMPLICIT none 
!
TYPE(internal_storage), POINTER :: ptr
CHARACTER(LEN=*), INTENT(IN)    :: strucfile 
!
CHARACTER(LEN=80) :: ier_msg_local 
!
!     Allocate sufficient space, even for all headers, and atom type, if they are omitted
!
call ptr%crystal%alloc_arrays(cr_natoms, cr_nscat, cr_ncatoms, cr_nanis, &
     mole_num_mole, mole_num_type, mole_num_atom ) ! Allocate the crystal arrays
!    mole_max_mole, mole_max_type, mole_max_atom ) ! Allocate the crystal arrays
!
!     An internal crystal has ALL headers saved, logical flags are used to indicate
!     whether they were supposed to be saved or not.
!     n_latom = UBOUND(sav_latom,1)     ! Make sure we send correct array size
call ptr%crystal%set_crystal_save_flags (sav_w_scat, & 
     sav_w_adp, sav_w_occ, sav_w_gene, sav_w_symm,                &
     sav_w_ncell, sav_w_obje, sav_w_doma, sav_w_mole, sav_w_prop, &
     sav_sel_prop,cr_nscat,sav_latom)
!          sav_sel_prop,MAXSCAT,sav_latom)
!
call ptr%crystal%set_crystal_from_standard(strucfile) ! Copy complete crystal
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
end subroutine save_internal_node
!
!*******************************************************************************
!
subroutine save_internal_backup(infile)
!-
!  Save the current structure as a backup with full settings
!
use crystal_mod
use discus_save_mod
use prop_para_func
use prop_para_mod
!
implicit none
!
character(len=*), intent(in) :: infile
!
character(len=PREC_STRING) :: line   ! Dummy line
integer                    :: length ! dummy length
!
call save_store_setting             ! Backup user "save" setting
call save_default_setting           ! Default to full saving
line       = 'ignore, all'          ! Ignore all properties
length     = 11
call property_select(line, length, sav_sel_prop)
line       = 'ignore, all'          ! Ignore all properties for global as well
length     = 11
call property_select(line, length,  cr_sel_prop)
!
!line = 'internal.' // infile(1:len_trim(infile))
line =                infile(1:len_trim(infile))
call save_internal(line)            !     thus this file name is unique
!
call save_restore_setting           ! Restore user save settings
!
end subroutine save_internal_backup
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
      INTEGER            :: n_nsite
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
         n_nsite = MAX(cr_ncatoms, SAV_MAXSCAT, MAXSCAT) 
         CALL alloc_save (  n_nscat, n_nsite )
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
      INTEGER            :: n_nsite
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT .OR. &
          SAV_T_MAXSCAT> SAV_MAXSCAT) THEN
         n_nscat = MAX(cr_nscat, SAV_MAXSCAT, SAV_T_MAXSCAT, MAXSCAT)
         n_nsite = MAX(cr_ncatoms, SAV_MAXSCAT, MAXSCAT) 
         CALL alloc_save (  n_nscat, n_nsite )
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
INTEGER            :: n_nsite
!                                                                       
maxw = MAX(MIN_PARA,MAXSCAT+1)
IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
   n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
   n_nsite = MAX(cr_ncatoms, SAV_MAXSCAT, MAXSCAT) 
   CALL alloc_save (  n_nscat, n_nsite )
   IF ( ier_num < 0 ) THEN
      RETURN
   ENDIF
ENDIF
!
SAV_MAXSCAT  = MAXSCAT
SAV_MAXSITE  = UBOUND(sav_lsite,1)
sav_latom(:) = .true.
sav_lsite(:) = .TRUE.
!
sav_sel_atom = .true.
!
sav_file    = 'crystal.stru'
!
sav_keyword = .true.
!
sav_w_scat  = .true.
sav_w_adp   = .true.
sav_w_occ   = .true.
sav_w_surf  = .false.
sav_w_magn  = .FALSE.
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
!
END SUBROUTINE save_default_setting
!
!*******************************************************************************
!
SUBROUTINE save_full_setting
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
INTEGER            :: n_nsite
!                                                                       
maxw = MAX(MIN_PARA,MAXSCAT+1)
IF( cr_nscat > SAV_MAXSCAT .or. MAXSCAT > SAV_MAXSCAT) THEN
   n_nscat = MAX(cr_nscat, SAV_MAXSCAT, MAXSCAT)
   n_nsite = MAX(cr_ncatoms, SAV_MAXSCAT, MAXSCAT) 
   CALL alloc_save (  n_nscat, n_nsite )
   IF ( ier_num < 0 ) THEN
      RETURN
   ENDIF
ENDIF
!
SAV_MAXSCAT  = MAXSCAT
SAV_MAXSITE  = UBOUND(sav_lsite,1)
sav_latom(:) = .true.
sav_lsite(:) = .TRUE.
!
sav_sel_atom = .true.
!
sav_file    = 'crystal.stru'
!
sav_keyword = .true.
!
sav_w_scat  = .true.
sav_w_adp   = .true.
sav_w_occ   = .true.
sav_w_surf  = .true.
sav_w_magn  = .true.
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
!
END SUBROUTINE save_full_setting
!
!*******************************************************************************
!
SUBROUTINE save_reset
!
USE discus_allocate_appl_mod
USE discus_save_mod

IMPLICIT NONE
!
CALL alloc_save(1, 1)
!
SAV_T_MAXSCAT = 1
SAV_MAXSCAT   = 1
SAV_MAXSITE   = 1
!
IF(ALLOCATED(sav_latom))   sav_latom(:)   = .TRUE. ! (0:MAXSCAT)
IF(ALLOCATED(sav_t_latom)) sav_t_latom(:) = .TRUE. ! (0:MAXSCAT)
IF(ALLOCATED(sav_lsite))   sav_lsite(:)   = .TRUE. ! (0:MAXSCAT)
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
sav_t_w_surf  = .FALSE.
sav_t_w_magn  = .FALSE.
sav_w_magn    = .FALSE.
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
   IF (sav_w_surf) THEN 
      WRITE (output_io, 3070) 'written' 
   ELSE 
      WRITE (output_io, 3070) 'omitted' 
   ENDIF 
   IF (sav_w_magn) THEN 
      WRITE (output_io, 3075) 'written' 
   ELSE 
      WRITE (output_io, 3075) 'omitted' 
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
 3070 FORMAT(' Surface  information: type, normal        : ',a7) 
 3075 FORMAT(' Magnetic moments                          : ',a7) 
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
