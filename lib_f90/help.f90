!****7******************************************************************
!                                                                       
      SUBROUTINE do_hel (ein, length) 
!-                                                                      
!     This sublevel emulates the HELP function under VMS and            
!     controls the online help of DISCUS.                               
!                                                                       
!     Version : 1.3                                                     
!     Date    : 25 Oct 96                                               
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)      
!               Th. Proffen (proffen@pa.msu.edu)                        
!                                                                       
!-                                                                      
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
      include'envir.inc' 
      include'doact.inc' 
      include'macro.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER nviele, ihl, maxw 
      PARAMETER (nviele = 101, ihl = 77, maxw = 10) 
!                                                                       
      CHARACTER ( * ) ein 
      CHARACTER(80) line, zeile, dummy 
      CHARACTER(20) prom 
      CHARACTER(12) bef (maxw) 
      CHARACTER(4) befehl 
      CHARACTER(3) status 
      INTEGER i, ibef, il (maxw) 
      INTEGER length, ll, lbef, lp 
      LOGICAL lread, lende 
      LOGICAL found, fast, stay 
!                                                                       
      INTEGER len_str 
!                                                                       
!------ set up some variables                                           
!                                                                       
      line = ein 
      ier_num = 0 
      ier_typ = ER_NONE 
      status = 'old' 
      lread = .true. 
      lende = .false. 
      fast = .false. 
!                                                                       
!------ if length was given negative we do not want to enter help level 
!                                                                       
      IF (length.lt.0) then 
         stay = .false. 
         length = - length 
      ELSE 
         stay = .true. 
      ENDIF 
!                                                                       
!------ if we are in a loop or macro do nothing !!                      
!                                                                       
      IF (lmakro.or.lblock) return 
   10 CONTINUE 
!                                                                       
!------ - print help file entry                                         
!                                                                       
      IF (stay) call do_status (1, bef, ibef, il, maxw) 
      CALL oeffne (ihl, hlpfile, status, lread) 
      IF (ier_num.ne.0) return 
      CALL do_cap (line) 
      CALL do_trenn (line, bef, ibef, length, maxw) 
      CALL do_help_single (ihl, bef, ibef, maxw, found, fast) 
      CLOSE (ihl) 
!                                                                       
!------ - Entry found ?                                                 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
!------ - if there are no subcommands get up one level                  
!                                                                       
      IF (.not.found) then 
         IF (ibef.gt.1) then 
            ibef = ibef - 1 
            WRITE (line, 1000) (bef (i), i = 1, ibef) 
            length = len_str (line) 
         ENDIF 
      ENDIF 
!                                                                       
!------ - Want to enter help level ?                                    
!                                                                       
      IF (.not.stay) return 
!                                                                       
!     - Set reading mode back to normal                                 
!                                                                       
      fast = .false. 
!                                                                       
!------ - prompt user for input in online help                          
!                                                                       
   15 CONTINUE 
      CALL do_status (2, bef, ibef, il, maxw) 
      prom = prompt (1:len_str (prompt) ) //'/help' 
      CALL get_cmd (zeile, ll, befehl, lbef, dummy, lp, prom) 
!                                                                       
!------ - we will go up one level '..'                                  
!                                                                       
      IF (zeile (1:2) .eq.'..') then 
         IF (ibef.gt.1) then 
            ibef = ibef - 1 
            WRITE (line, 1000) (bef (i), i = 1, ibef) 
            length = len_str (line) 
         ENDIF 
         GOTO 15 
!                                                                       
!------ - leave help level ' '                                          
!                                                                       
      ELSEIF (ll.eq.0) then 
         lende = .true. 
!                                                                       
!     - repeat list of sublevels availble from current level            
!                                                                       
      ELSEIF (zeile (1:1) .eq.'?') then 
         fast = .true. 
!                                                                       
!------ - give help for this subcommand 'word'                          
!                                                                       
      ELSE 
         WRITE (line, 1000) (bef (i), i = 1, ibef), zeile 
         length = len_str (line) 
      ENDIF 
!                                                                       
      IF (.not.lende) goto 10 
!                                                                       
   20 CONTINUE 
!                                                                       
 1000 FORMAT    (8(a12,1x)) 
 2000 FORMAT    (a) 
      END SUBROUTINE do_hel                         
!****7******************************************************************
      SUBROUTINE do_status (iwhere, bef, ibef, il, maxw) 
!-                                                                      
!     Prints out status line in online help                             
!-                                                                      
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
!                                                                       
      INTEGER iwhere, ibef, maxw 
      INTEGER il (maxw) 
      CHARACTER ( * ) bef (maxw) 
!                                                                       
      CHARACTER(80) level 
      CHARACTER(6) prom 
      INTEGER i 
!                                                                       
      INTEGER len_str 
!                                                                       
!------ Build status lines                                              
!                                                                       
      IF (iwhere.eq.1) then 
         prom = prompt 
         CALL do_cap (prom) 
         WRITE (output_io, 2000) prom, version 
      ELSEIF (iwhere.eq.2) then 
         DO i = 1, ibef 
         il (i) = len_str (bef (i) ) 
         ENDDO 
         WRITE (level, 1000) (bef (i) (1:il (i) ), i = 1, ibef) 
         WRITE (output_io, 3000) level 
      ENDIF 
!                                                                       
 1000 FORMAT    (8(a,1x)) 
 2000 FORMAT    (' ----------- > ',A6,' online help ',25x,              &
     &                  'Version : ',a5,/)                              
 3000 FORMAT    (' ----------- > current help level : ',a40) 
      END SUBROUTINE do_status                      
!****7******************************************************************
      SUBROUTINE do_help_single (ihl, bef, ibef, maxw, found, fast) 
!-                                                                      
!     This routine emulates the VMS help function itself                
!-                                                                      
      IMPLICIT none 
!                                                                       
      include'envir.inc' 
      include'prompt.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER nviele 
      PARAMETER (nviele = 101) 
!                                                                       
      INTEGER ihl, ibef, maxw 
      CHARACTER ( * ) bef (maxw) 
!                                                                       
      CHARACTER(80) line 
      CHARACTER(80) zeile 
      CHARACTER(12) viele (nviele) 
      INTEGER lev, nl, nbef, ianz, laenge, ll 
      LOGICAL found, fast 
      REAL werte 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
!------ Loop over all helpfile levels until entry is found              
!                                                                       
      DO lev = 1, ibef 
    5 CONTINUE 
      found = .false. 
      READ (ihl, 1000, end = 10) zeile 
      ll = len_str (zeile) 
      line = zeile (1:ll) 
!                                                                       
      IF (line (1:1) .eq.'!') then 
         line (1:1) = ' ' 
         line (2:2) = ' ' 
      ENDIF 
!                                                                       
      CALL do_cap (line) 
      CALL hole_zahl (line (1:3), ianz, werte) 
      IF (ianz.eq.1) then 
         nl = int (werte) 
         IF (nl.eq.lev) then 
            laenge = len_str (bef (lev) ) 
            found = str_comp (line (4:4 + laenge-1), bef (lev), laenge, &
            laenge, laenge)                                             
            IF (.not.found) goto 5 
         ELSEIF (nl.gt.lev.or.nl.eq.0) then 
            GOTO 5 
         ELSEIF (nl.lt.lev) then 
            GOTO 10 
         ENDIF 
      ELSE 
         GOTO 5 
      ENDIF 
      ENDDO 
!                                                                       
!------ If entry is found, print text and get possible further          
!------ subentries                                                      
!                                                                       
   10 CONTINUE 
      nbef = 0 
      IF (found) then 
         IF (.not.fast) then 
            CALL lese_text (ihl) 
         ELSE 
            READ (ihl, 1000, end = 10) zeile 
            ll = len_str (zeile) 
            line = zeile (1:ll) 
!                                                                       
            IF (line (1:1) .eq.'!') then 
               line (1:1) = ' ' 
               line (2:2) = ' ' 
            ENDIF 
!                                                                       
         ENDIF 
         CALL weitere (ihl, ibef + 1, viele, nviele, nbef) 
         CALL schreib_viele (viele, nviele, nbef) 
         WRITE (output_io, * ) ' ' 
      ELSE 
         ier_typ = ER_IO 
         ier_num = - 5 
         CALL errlist 
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
!                                                                       
!     Mark whether subentries were found                                
!                                                                       
      found = nbef.gt.0 
!                                                                       
 1000 FORMAT  (a) 
      END SUBROUTINE do_help_single                 
!*****7*****************************************************************
      SUBROUTINE hole_zahl (zeile, ianz, werte) 
!+                                                                      
!     Gets numbers from command line and number of values of input.     
!-                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER ianz 
      REAL werte 
!                                                                       
      werte = 0.0 
      ianz = 0 
      READ (zeile, *, end = 98, err = 98) werte 
      ianz = 1 
      RETURN 
   98 CONTINUE 
      ianz = 0 
      END SUBROUTINE hole_zahl                      
!*****7*****************************************************************
      SUBROUTINE do_trenn (ein, bef, ibef, length, maxw) 
!+                                                                      
!     Separates the character string into single words                  
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER maxw, ibef 
      CHARACTER ( * ) ein, bef (maxw) 
      INTEGER ib, i, ianf, iend, length 
      LOGICAL lchar 
!                                                                       
      ib = 1 
      DO i = 1, 10 
      lchar = .true. 
      CALL locate (lchar, ein, ib, ianf, length) 
      IF (ianf.gt.length) goto 20 
      lchar = .false. 
      ib = ianf 
      CALL locate (lchar, ein, ib, iend, length) 
      bef (i) = ein (ianf:iend-1) 
      ib = iend 
      ENDDO 
   20 CONTINUE 
      ibef = i - 1 
      END SUBROUTINE do_trenn                       
!*****7*****************************************************************
      SUBROUTINE locate (lchar, ein, ib, ianf, length) 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) ein 
      INTEGER ianf, length, i, ib 
      LOGICAL lchar 
!                                                                       
      ianf = length + 1 
!                                                                       
      DO i = ib, length 
      IF (lchar.and.ein (i:i) .ne.' '.and.ein (i:i) .ne.',') then 
         ianf = i 
         GOTO 20 
      ELSEIF (.not.lchar.and. (ein (i:i) .eq.' '.OR.ein (i:i) .eq.',') )THEN                                                              
         ianf = i 
         GOTO 20 
      ENDIF 
      ENDDO 
   20 CONTINUE 
      END SUBROUTINE locate                         
!*****7*****************************************************************
      SUBROUTINE lese_text (ihl) 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'envir.inc' 
      include'prompt.inc' 
!                                                                       
      CHARACTER(80) line 
      CHARACTER(1) cdummy 
      INTEGER ihl, ilc, ll 
!                                                                       
      INTEGER len_str 
!                                                                       
   10 CONTINUE 
      DO ilc = 1, lines - 2 
      READ (ihl, 1000, end = 9010) line 
!                                                                       
      IF (line (1:1) .eq.'!') then 
         line (1:1) = ' ' 
         line (2:2) = ' ' 
      ENDIF 
!                                                                       
      IF (line (1:1) .ne.' '.and.line (1:1) .ne.achar (13) ) then 
         GOTO 20 
      ELSE 
         ll = len_str (line) 
         IF (ll.gt.0) then 
            WRITE (output_io, 1000) line (1:ll) 
         ELSE 
            WRITE (output_io, * ) 
         ENDIF 
      ENDIF 
      ENDDO 
      WRITE (output_io, 2000) 
      IF (output_io.eq.6) read ( *, 1000) cdummy 
      GOTO 10 
   20 CONTINUE 
!                                                                       
 9010 CONTINUE 
 1000 FORMAT  (a) 
 2000 FORMAT    (/,' ----------- > Press <RETURN> for next page ') 
      END SUBROUTINE lese_text                      
!*****7*****************************************************************
      SUBROUTINE weitere (ihl, nl, viele, nviele, nbef) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER nviele 
      CHARACTER ( * ) viele (nviele) 
      CHARACTER(80) line 
      CHARACTER(80) zeile 
      INTEGER ihl, nl, nbef, ianz 
      INTEGER ll 
      REAL werte 
!                                                                       
      CHARACTER(1) cdummy 
!                                                                       
      INTEGER len_str 
!                                                                       
      nbef = 0 
      BACKSPACE (ihl) 
   30 CONTINUE 
      READ (ihl, 1000, end = 9010) zeile 
      ll = len_str (zeile) 
      line = zeile (1:ll) 
      DO while (line (1:1) .eq.' '.or.line (1:1) .eq.achar (13) ) 
      READ (ihl, 1000, end = 9010) zeile 
      ll = len_str (zeile) 
      line = zeile (1:ll) 
      ENDDO 
      CALL hole_zahl (line, ianz, werte) 
      IF (ianz.eq.1) then 
         IF (int (werte) .eq.nl) then 
            nbef = nbef + 1 
            viele (nbef) = line (4:15) 
         ELSEIF (int (werte) .lt.nl) then 
            RETURN 
         ENDIF 
      ENDIF 
      GOTO 30 
 9010 CONTINUE 
 1000 FORMAT  (a) 
      END SUBROUTINE weitere                        
!*****7*****************************************************************
      SUBROUTINE schreib_viele (viele, nviele, nbef) 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
!                                                                       
      INTEGER nviele 
      CHARACTER ( * ) viele (nviele) 
      INTEGER nbef, i 
!                                                                       
      WRITE (output_io, 1000) (viele (i), i = 1, nbef) 
!                                                                       
 1000 FORMAT  (5(3x,a12)) 
      END SUBROUTINE schreib_viele                  
