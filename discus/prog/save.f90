!+                                                                      
!     Menu for the structure saveing part of DISCUS. The user           
!     can define what structural information is saves to the structure  
!     file.                                                             
!     If the Menu is called with one or two parameters, the structure   
!     is written immedeately to the file defined by the first parameter.
!     In this case, the format is taken from the current settings.      
!                                                                       
!*****7*****************************************************************
      SUBROUTINE save_struc (zeile, lcomm) 
!-                                                                      
!     Main menu for generalized transformation operations               
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE modify_mod
      USE save_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'doact.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'macro.inc' 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA =  20 ! A command requires at least these no of parameters
      INTEGER maxw 
      LOGICAL lnew, lold 
      PARAMETER (lnew = .true., lold = .false.) 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
!
      CHARACTER ( LEN=* ) zeile 
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz, i 
      INTEGER lcomm, sav_flen 
      LOGICAL lend 
      REAL, DIMENSION(SAV_MAXSCAT) :: repl ! Dummy variable needed for atom_select
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      DATA sav_flen / 1 / 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      IF( cr_nscat > SAV_MAXSCAT) THEN
         CALL alloc_save ( cr_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!                                                                       
!     Interpret parameters used by 'save' command                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
      IF (ier_num.eq.0) then 
         IF (ianz.gt.0) then 
!                                                                       
!     ----Parameters were found, write file immedeatly, and return      
!         to calling routine                                            
!                                                                       
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.eq.0) then 
               sav_file = cpara (1) 
               sav_flen = lpara (1) 
!                                                                       
               IF (sav_keyword) then 
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
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/save' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.ne.' '.and.line (1:1) .ne.'#') then 
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') then 
                  IF (length.ge.2) then 
                     CALL file_kdo (line (2:length), length - 1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----list asymmetric unit 'asym'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
                  CALL show_asym 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----list atoms present in the crystal 'chem'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'chem', 2, lbef, 4) ) then 
                  CALL show_chem 
!                                                                       
!     ----Deselect which atoms are included in the wave 'dese'          
!                                                                       
               ELSEIF (str_comp (befehl, 'dese', 1, lbef, 4) ) then 
                   CALL atom_select (zeile, lp, 0, SAV_MAXSCAT, sav_latom, &
                   sav_sel_atom, .true., .false.)              
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num.eq.0) then 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num.eq.0) then 
!                       IF (werte (1) .eq. - 1) then 
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
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
                  CALL do_eval (zeile, lp) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     define format of output file 'format'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'format', 1, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (1) , 'keyword', 1, lpara (1) ,&
                     7) ) then                                          
                        sav_keyword = .true. 
                     ELSEIF (str_comp (cpara (1) , 'nokeywo', 1, lpara (&
                     1) , 7) ) then                                     
                        sav_keyword = .false. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  line = zeile 
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//line, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus save '//line, lp) 
                  ENDIF 
!                                                                       
!     ----Select range of atoms within crystal to be included 'incl'    
!                                                                       
               ELSEIF (str_comp (befehl, 'incl', 1, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           sav_start = nint (werte (1) ) 
                           sav_end = nint (werte (2) ) 
                        ENDIF 
                     ELSEIF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) then                                       
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
               ELSEIF (str_comp (befehl, 'omit', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) &
                     ) then                                             
                        sav_w_gene = .false. 
                        sav_w_ncell = .false. 
                        sav_w_symm = .false. 
                        sav_w_scat = .false. 
                        sav_w_adp = .false. 
                        sav_w_obje = .false. 
                        sav_w_doma = .false. 
                     ELSEIF (str_comp (cpara (1) , 'gene', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_gene = .false. 
                     ELSEIF (str_comp (cpara (1) , 'mole', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_mole = .false. 
                     ELSEIF (str_comp (cpara (1) , 'obje', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_obje = .false. 
                     ELSEIF (str_comp (cpara (1) , 'doma', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_doma = .false. 
                     ELSEIF (str_comp (cpara (1) , 'ncell', 1, lpara (1)&
                     , 5) ) then                                        
                        sav_w_ncell = .false. 
                     ELSEIF (str_comp (cpara (1) , 'symm', 2, lpara (1) &
                     , 4) ) then                                        
                        sav_w_symm = .false. 
                     ELSEIF (str_comp (cpara (1) , 'scat', 2, lpara (1) &
                     , 4) ) then                                        
                        sav_w_scat = .false. 
                     ELSEIF (str_comp (cpara (1) , 'adp', 1, lpara (1) ,&
                     3) ) then                                          
                        sav_w_adp = .false. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     define name of output file 'outfile'                              
!                                                                       
               ELSEIF (str_comp (befehl, 'outfile', 1, lbef, 7) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
                     IF (ier_num.eq.0) then 
                        sav_file = cpara (1) 
                        sav_flen = lpara (1) 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----run transformation 'run'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 1, lbef, 4) ) then 
                  IF (sav_keyword) then 
                     IF (str_comp (sav_file(1:8),'internal',8,8,8)) THEN
                        CALL save_internal (sav_file) 
                     ELSE 
                        CALL save_keyword (sav_file) 
                     ENDIF 
                  ELSE 
                     CALL save_nokeyword (sav_file) 
                  ENDIF 
!                                                                       
!     ----Select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'sele', 2, lbef, 4) ) then 
                   CALL atom_select (zeile, lp, 0, SAV_MAXSCAT, sav_latom, &
                   sav_sel_atom, .true., .true.)               
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num.eq.0) then 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num.eq.0) then 
!                       IF (werte (1) .eq. - 1) then 
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
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  IF (sav_flen.gt.0) then 
                     WRITE (output_io, 3000) sav_file (1:sav_flen) 
                  ELSE 
                     WRITE (output_io, 3000) ' ' 
                  ENDIF 
                  IF (sav_keyword) then 
                     WRITE (output_io, 3005) 
                     IF (sav_w_scat) then 
                        WRITE (output_io, 3008) 'written' 
                     ELSE 
                        WRITE (output_io, 3008) 'omitted' 
                     ENDIF 
                     IF (sav_w_adp) then 
                        WRITE (output_io, 3009) 'written' 
                     ELSE 
                        WRITE (output_io, 3009) 'omitted' 
                     ENDIF 
                     IF (sav_w_gene) then 
                        WRITE (output_io, 3010) 'written' 
                     ELSE 
                        WRITE (output_io, 3010) 'omitted' 
                     ENDIF 
                     IF (sav_w_symm) then 
                        WRITE (output_io, 3020) 'written' 
                     ELSE 
                        WRITE (output_io, 3020) 'omitted' 
                     ENDIF 
                     IF (sav_w_ncell) then 
                        WRITE (output_io, 3030) 'written' 
                     ELSE 
                        WRITE (output_io, 3030) 'omitted' 
                     ENDIF 
                     IF (sav_w_mole) then 
                        WRITE (output_io, 3040) 'written' 
                     ELSE 
                        WRITE (output_io, 3040) 'omitted' 
                     ENDIF 
                     IF (sav_w_obje) then 
                        WRITE (output_io, 3050) 'written' 
                     ELSE 
                        WRITE (output_io, 3050) 'omitted' 
                     ENDIF 
                     IF (sav_w_doma) then 
                        WRITE (output_io, 3060) 'written' 
                     ELSE 
                        WRITE (output_io, 3060) 'omitted' 
                     ENDIF 
!                                                                       
                     WRITE (output_io, 3090) 
                     WRITE (output_io, 3091) 
                     DO i = 0, cr_nscat 
                     IF (sav_latom (i) ) then 
                        WRITE (output_io, 3092) i, cr_at_lis (i) 
                     ENDIF 
                     ENDDO 
                     IF (sav_end.eq. - 1) then 
                        WRITE (output_io, 3080) 
                     ELSE 
                        WRITE (output_io, 3081) sav_start, sav_end 
                     ENDIF 
                  ELSE 
                     WRITE (output_io, 3006) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 2, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!                                                                       
!     ----write parameters from the keyword list  'write'               
!                                                                       
               ELSEIF (str_comp (befehl, 'write', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) &
                     ) then                                             
                        sav_w_gene = .true. 
                        sav_w_ncell = .true. 
                        sav_w_symm = .true. 
                        sav_w_scat = .true. 
                        sav_w_adp = .true. 
                        sav_w_obje = .true. 
                        sav_w_doma = .true. 
                     ELSEIF (str_comp (cpara (1) , 'gene', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_gene = .true. 
                     ELSEIF (str_comp (cpara (1) , 'mole', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_mole = .true. 
                     ELSEIF (str_comp (cpara (1) , 'obje', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_obje = .true. 
                     ELSEIF (str_comp (cpara (1) , 'doma', 1, lpara (1) &
                     , 4) ) then                                        
                        sav_w_doma = .true. 
                     ELSEIF (str_comp (cpara (1) , 'ncell', 1, lpara (1)&
                     , 5) ) then                                        
                        sav_w_ncell = .true. 
                     ELSEIF (str_comp (cpara (1) , 'symm', 2, lpara (1) &
                     , 4) ) then                                        
                        sav_w_symm = .true. 
                     ELSEIF (str_comp (cpara (1) , 'scat', 2, lpara (1) &
                     , 4) ) then                                        
                        sav_w_scat = .true. 
                     ELSEIF (str_comp (cpara (1) , 'adp', 1, lpara (1) ,&
                     3) ) then                                          
                        sav_w_adp = .true. 
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
      ENDDO 
!                                                                       
 9999 CONTINUE 
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
 3010 FORMAT(' Additional generator matrizes             : ',a7) 
 3020 FORMAT(' Additional symmetry operators             : ',a7) 
 3030 FORMAT(' Number of unit cells, atoms per unit cell : ',a7) 
 3040 FORMAT(' Molecule information: content etc.        : ',a7) 
 3050 FORMAT(' Object   information: content etc.        : ',a7) 
 3060 FORMAT(' Domain   information: content etc.        : ',a7) 
 3080 FORMAT(' Range of atoms from to    : All atoms included') 
 3081 FORMAT(' Range of atoms from to    : ',2(2x,i9)) 
 3090 FORMAT(' Selected atoms    :') 
 3091 FORMAT('                      type name') 
 3092 FORMAT(20x,2(4x,i2,1x,a4)) 
!                                                                       
      END SUBROUTINE save_struc                     
!*****7*****************************************************************
      SUBROUTINE save_nokeyword (strucfile) 
!+                                                                      
!     This subroutine saves the structure and/or the unit cell          
!     onto a file in the format used until DISCUS 3.0 i.e. without      
!     keywords.                                                         
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
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
      CALL oeffne (ist, strucfile, 'unknown', lread) 
      IF (ier_num.eq.0) then 
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
      USE config_mod 
      USE crystal_mod 
      USE gen_add_mod 
      USE molecule_mod 
      USE sym_add_mod 
      USE save_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER ( LEN=* ) strucfile 
      CHARACTER(31) fform 
      CHARACTER(15) C_MOLE ( - 4:4) 
      INTEGER ist, i, j, k, j_end, l1, l2, l 
      INTEGER i_start, i_end 
      INTEGER is, ie 
      LOGICAL lread 
      LOGICAL lwritten 
!                                                                       
      INTEGER len_str 
!                                                                       
!                                                                       
      DATA ist / 67 / 
      DATA C_MOLE / 'domain_fuzzy   ', 'domain_sphere  ', 'domain_cylind&
     &er', 'domain_cube    ', 'atoms          ', 'cube           ', 'cyl&
     &inder       ', 'sphere         ', 'cube           ' /             
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sav_start 
      i_end = sav_end 
      IF (sav_end.eq. - 1) i_end = cr_natoms 
!                                                                       
      lread = .false. 
      CALL oeffne (ist, strucfile, 'unknown', lread) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
!                                                                       
!-----      --Write new type of structur file                           
!                                                                       
      j = len_str (cr_name) 
      IF (j.eq.0) then 
         cr_name = ' ' 
         j = 1 
      ENDIF 
      WRITE (ist, 3000) cr_name (1:j) 
      IF (spcgr_ianz.eq.0) then 
         WRITE (ist, 3010) cr_spcgr 
      ELSEIF (spcgr_ianz.eq.1) then 
         WRITE (ist, 3011) cr_spcgr, spcgr_para 
      ENDIF 
      WRITE (ist, 3020) (cr_a0 (i), i = 1, 3), (cr_win (i), i = 1, 3) 
      IF (sav_w_scat) then 
         j = (cr_nscat - 1) / 7 
         DO i = 1, j 
         is = (i - 1) * 7 + 1 
         ie = is + 6 
         WRITE (ist, 3110) (cr_at_lis (k), k = is, ie) 
         ENDDO 
         IF (cr_nscat - j * 7.eq.1) then 
            WRITE (ist, 3111) cr_at_lis (cr_nscat) 
         ELSEIF (cr_nscat - j * 7.gt.1) then 
            WRITE (fform, 7010) cr_nscat - j * 7 - 1 
            WRITE (ist, fform) (cr_at_lis (i), i = j * 7 + 1, cr_nscat) 
         ENDIF 
      ENDIF 
      IF (sav_w_adp) then 
         j = (cr_nscat - 1) / 7 
         DO i = 1, j 
         is = (i - 1) * 7 + 1 
         ie = is + 6 
         WRITE (ist, 3120) (cr_dw (k), k = is, ie) 
         ENDDO 
         IF (cr_nscat - j * 7.eq.1) then 
            WRITE (ist, 3121) cr_dw (cr_nscat) 
         ELSEIF (cr_nscat - j * 7.gt.1) then 
            WRITE (fform, 7020) cr_nscat - j * 7 - 1 
            WRITE (ist, fform) (cr_dw (i), i = j * 7 + 1, cr_nscat) 
         ENDIF 
      ENDIF 
      IF (sav_w_gene) then 
         DO k = 1, gen_add_n 
         WRITE (ist, 3021) ( (gen_add (i, j, k), j = 1, 4), i = 1, 3),  &
         gen_add_power (k)                                              
         ENDDO 
      ENDIF 
      IF (sav_w_symm) then 
         DO k = 1, sym_add_n 
         WRITE (ist, 3022) ( (sym_add (i, j, k), j = 1, 4), i = 1, 3),  &
         sym_add_power (k)                                              
         ENDDO 
      ENDIF 
      IF (sav_w_ncell) then 
         WRITE (ist, 3030) cr_icc, cr_ncatoms 
      ENDIF 
      WRITE (ist, 3900) 
!                                                                       
!     Write content of objects                                          
!                                                                       
      IF (sav_w_obje) then 
         DO i = 1, mole_num_mole 
         IF (mole_char (i) .gt.MOLE_ATOM) then 
            WRITE (ist, 4000) 'object' 
            WRITE (ist, 4002) 'object', mole_type (i) 
            WRITE (ist, 4100) 'object', c_mole (mole_char (i) ) 
            WRITE (ist, 4200) 'object', mole_dens (i) 
            DO j = 1, mole_len (i) 
            k = mole_cont (mole_off (i) + j) 
            WRITE (ist, 4) cr_at_lis (cr_iscat (k) ), (cr_pos (l, k),   &
            l = 1, 3), cr_dw (cr_iscat (k) ), cr_prop (k)               
            ENDDO 
            WRITE (ist, 4900) 'object' 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Write content of micro domains                                    
!                                                                       
      IF (sav_w_doma) then 
         DO i = 1, mole_num_mole 
         IF (mole_char (i) .lt.MOLE_ATOM) then 
            WRITE (ist, 4000) 'domain' 
            WRITE (ist, 4002) 'domain', mole_type (i) 
            WRITE (ist, 4100) 'domain', c_mole (mole_char (i) ) 
            k = len_str (mole_file (i) ) 
            IF (k.gt.0) then 
               WRITE (ist, 4300) 'domain', mole_file (i) (1:k) 
            ELSE 
               WRITE (ist, 4300) 'domain' 
            ENDIF 
            WRITE (ist, 4400) 'domain', mole_fuzzy (i) 
            DO j = 1, mole_len (i) 
            k = mole_cont (mole_off (i) + j) 
            WRITE (ist, 4) cr_at_lis (cr_iscat (k) ), (cr_pos (l, k),   &
            l = 1, 3), cr_dw (cr_iscat (k) ), cr_prop (k)               
            ENDDO 
            WRITE (ist, 4900) 'domain' 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Write content of molecules                                        
!                                                                       
      IF (sav_w_mole) then 
         DO i = 1, mole_num_mole 
         IF (mole_char (i) .eq.MOLE_ATOM) then 
            WRITE (ist, 4000) 'molecule' 
            WRITE (ist, 4002) 'molecule', mole_type (i) 
            WRITE (ist, 4100) 'molecule', c_mole (mole_char (i) ) 
            DO j = 1, mole_len (i) 
            k = mole_cont (mole_off (i) + j) 
            WRITE (ist, 4) cr_at_lis (cr_iscat (k) ), (cr_pos (l, k),   &
            l = 1, 3), cr_dw (cr_iscat (k) ), cr_prop (k)               
            ENDDO 
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
      IF (sav_latom (cr_iscat (i) ) ) then 
         lwritten = .false. 
         IF (sav_w_obje) then 
            DO j = 1, mole_num_mole 
            IF (mole_char (j) .gt.MOLE_ATOM) then 
               DO k = 1, mole_len (j) 
               l = mole_cont (mole_off (j) + k) 
               lwritten = lwritten.or.l.eq.i 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
         IF (sav_w_doma) then 
            DO j = 1, mole_num_mole 
            IF (mole_char (j) .lt.MOLE_ATOM) then 
               DO k = 1, mole_len (j) 
               l = mole_cont (mole_off (j) + k) 
               lwritten = lwritten.or.l.eq.i 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
         IF (sav_w_mole) then 
            DO j = 1, mole_num_mole 
            IF (mole_char (j) .eq.MOLE_ATOM) then 
               DO k = 1, mole_len (j) 
               l = mole_cont (mole_off (j) + k) 
               lwritten = lwritten.or.l.eq.i 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
         IF (.not.lwritten) then 
            WRITE (ist, 4) cr_at_lis (cr_iscat (i) ), (cr_pos (j, i),   &
            j = 1, 3), cr_dw (cr_iscat (i) ), cr_prop (i)               
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     Write content of molecules ; OBSOLETE ! RBN                       
!           Molecules are now written just like domains                 
!                                                                       
!DBG      if(sav_w_mole) then                                           
!DBG        do i=1,mole_num_mole                                        
!DBG          if(mole_char(i).eq.MOLE_ATOM) then                        
!DBG            write(ist,4005) 'molecule',mole_type(i),i               
!DBG            write(ist,4100) 'molecule',c_mole(mole_char(i))         
!DBG            j_end = int(mole_len(i)/7.)                             
!DBG            do j=1,j_end                                            
!DBG              l1 = mole_off(i) + (j-1)*7 + 1                        
!DBG              l2 = l1 + 7 - 1                                       
!DBG              write(ist,4010) 'molecule',(',',mole_cont(k),k=l1,l2) 
!DBG            ENDDO                                                   
!DBG            if(j_end*7.lt.mole_len(i)) then                         
!DBG              l1= mole_off(i) + int((mole_len(i)-1)/7.)*7 + 1       
!DBG              l2= mole_off(i) + mole_len(i)                         
!DBG              write(ist,4010) 'molecule',(',',mole_cont(k),k=l1,l2) 
!DBG            endif                                                   
!DBG            write(ist,4900) 'molecule'                              
!DBG          endif                                                     
!DBG        ENDDO                                                       
!DBG      endif                                                         
      CLOSE (ist) 
!                                                                       
 3000 FORMAT    ('title ',a) 
 3010 FORMAT    ('spcgr ',a16) 
 3011 FORMAT    ('spcgr ',a16,',',i4) 
 3020 FORMAT    ('cell  ',5(f10.6,2x ),f10.6) 
 3021 FORMAT    ('gene  ',12(f9.6,','),i3) 
 3022 FORMAT    ('symm  ',12(f9.6,','),i3) 
 3030 FORMAT    ('ncell ',3(i8,','),i10) 
 3110 FORMAT    (('scat  ', a4,6(',',5x,a4))) 
 3111 FORMAT    (('scat  ', a4)) 
 3120 FORMAT    (('adp   ', f9.6,6(',',5x,f9.6))) 
 3121 FORMAT    (('adp   ', f9.6)) 
 3900 FORMAT    ('atoms ') 
                                                                        
 4000 FORMAT    (a) 
 4002 FORMAT    (a,' type,',i8) 
 4005 FORMAT    (a,' content,',i8,',',i8) 
 4010 FORMAT    (a,' atoms',7(a1,i8)) 
 4100 FORMAT    (a,' character,',a15) 
 4200 FORMAT    (a,' density  ,',f12.4) 
 4300 FORMAT    (a,' file     ,',a) 
 4400 FORMAT    (a,' fuzzy    ,',f12.4) 
 4900 FORMAT    (a,' end') 
    4 FORMAT (a4,3(1x,f14.6,','),4x,f10.6,',',i8) 
 7010 FORMAT    ('(''scat  '', a4  ,',i1,'('','',5x,a4  ))') 
 7020 FORMAT    ('(''adp   '', f9.6,',i1,'('','',5x,f9.6))') 
      END SUBROUTINE save_keyword                   
!********************************************************************** 
      SUBROUTINE save_internal (strucfile) 
!+                                                                      
!     This subroutine saves the structure and/or the unit cell          
!     onto a file. The format uses keyword description.                 
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE gen_add_mod 
      USE class_internal
      USE molecule_mod 
      USE sym_add_mod 
      USE save_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER ( LEN=* ), INTENT(IN) :: strucfile 
!
      INTEGER                         :: i_start, i_end 
      INTEGER                         :: istatus
!
      ALLOCATE(store_temp, STAT = istatus)        ! Allocate a temporary storage
      store_temp%strucfile = strucfile            ! Copy filename
!
      CALL store_add_node(store_root, store_temp) ! add this node to storage tree
      ALLOCATE(store_temp%crystal, STAT=istatus)  ! Allocate the crystal at this node
      CALL store_temp%crystal%alloc_arrays(cr_natoms,cr_nscat) ! Allocate the crystal arrays
      CALL store_temp%crystal%set_crystal_from_standard(strucfile) ! Copy complete crystal
!
!     An internal crystal has ALL headers saved, logical flags are used to indicate
!     whether they were supposed to be saved or not.
      CALL store_temp%crystal%set_crystal_save_flags (sav_w_scat, & 
           sav_w_adp, sav_w_gene, sav_w_symm,                     &
           sav_w_ncell, sav_w_obje, sav_w_doma, sav_w_mole)

write(*,*) ' SAVED the crystal'
!
      END SUBROUTINE save_internal
