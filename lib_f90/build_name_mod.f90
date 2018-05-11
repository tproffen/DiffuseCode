MODULE build_name_mod
!
CONTAINS
!
!*****7***************************************************************  
!
SUBROUTINE do_build_name (ianz, cpara, lpara, werte, MAXW, fpara) 
!-                                                                      
!     Creates a filename from several parameters. the first parameter   
!     is interpreted as format.                                         
!+                                                                      
      USE blanks_mod
      USE berechne_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, PARAMETER :: FORM_INTEGER    = 1
      INTEGER, PARAMETER :: FORM_INT_ZEROS  = 2
      INTEGER, PARAMETER :: FORM_REAL       = 3
      INTEGER, PARAMETER :: FORM_REAL_ZEROS = 4
      INTEGER, PARAMETER :: FORM_CHARACTER  = 5
      INTEGER, PARAMETER :: FORM_CHAR_ZEROS = 6
      INTEGER, PARAMETER :: NUM_FORM        = 6
!                                                                       
INTEGER                             , INTENT(INOUT) :: ianz 
INTEGER                             , INTENT(IN   ) :: MAXW 
CHARACTER(LEN=1024), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara
REAL               , DIMENSION(MAXW), INTENT(INOUT) :: werte
INTEGER                             , INTENT(IN)    :: fpara 
!                                                                       
      CHARACTER(1024) string 
!                                                                       
      CHARACTER(1) c_form (NUM_FORM) 
      CHARACTER(1024) fstring 
      CHARACTER(1024) line, number 
      INTEGER ftyp 
      INTEGER ind_ql, ind_qr 
      INTEGER ind_d, ind_i, ind_f 
      INTEGER ind_p, itot, idec 
      INTEGER ii, lp 
      INTEGER laenge 
      INTEGER lwert 
      INTEGER lsub 
      INTEGER npara 
      INTEGER pos 
      INTEGER ll, i 
      INTEGER fstring_l 
      INTEGER ising 
      REAL wert 
      LOGICAL lfloat 
      LOGICAL l_hyp 
      LOGICAL l_sin 
!                                                                       
      INTEGER len_str 
!                                                                       
      DATA c_form / 'd', 'D', 'f', 'F', 'c', 'C' / 
!                                                                       
      string = ' ' 
      lwert = 0
      wert = 0.
      fstring_l = 0
      IF (lpara (fpara) .le.1) then 
!                                                                       
!     --If first parameter is too short for a string, return            
!                                                                       
         RETURN 
      ENDIF 
!                                                                       
!     Remove leading blanks                                             
!                                                                       
      DO i = 1, ianz 
      ll = 1 
      line (1:1) = cpara (i) (ll:ll) 
      DO while (line (1:1) .eq.' ') 
      ll = ll + 1 
      line (1:1) = cpara (i) (ll:ll) 
      IF (ll.gt.lpara (i) ) goto 7575 
      ENDDO 
 7575 CONTINUE 
      string = cpara (i) (ll:lpara (i) ) 
      cpara (i) = string (1:lpara (i) - ll + 1) 
      lpara (i) = lpara (i) - ll + 1 
      ENDDO 
!                                                                       
!     Get position of left quotation mark into ind_ql                   
!                                                                       
      ind_ql = INDEX (cpara (fpara) (1:lpara (fpara) - 1) , '"') 
      IF (ind_ql.lt.1) then 
!                                                                       
!     --No quotation mark found, return                                 
!                                                                       
         RETURN 
      ENDIF 
!                                                                       
!     Get position of right quotation mark into ind_qr                  
!                                                                       
      ind_qr = ind_ql + INDEX(cpara(fpara)(ind_ql + 1:lpara(fpara)),'"')
      IF (ind_qr.lt.ind_ql + 1) then 
         ier_num = - 30 
         ier_typ = ER_FORT 
         cpara (fpara) = ' ' 
         lpara (fpara) = 1 
         RETURN 
      ENDIF 
!                                                                       
!     Find position of first format specifiers                          
!                                                                       
      ind_d = INDEX (cpara (fpara)(ind_ql + 1:ind_qr - 1) , '%') + ind_ql
!                                                                       
      laenge = 0 
      pos = ind_ql + 1 
      npara = fpara + 1 
!                                                                       
!     While there are any format specifiers left                        
!                                                                       
      DO while (ind_d.ge.pos) 
!                                                                       
!     --If there is space left for a format specifier                   
!                                                                       
      IF (ind_d+1.le.ind_qr - 1) then 
!                                                                       
!     ----Find first 'd' or 'f'                                         
!                                                                       
         ind_f = 9999 
         ftyp = 0 
         DO i = 1, NUM_FORM 
         ind_i = INDEX (cpara (fpara) (ind_d+1:ind_qr - 1), c_form (i) ) 
         IF (ind_i.gt.0.and.ind_i.lt.ind_f) then 
            ind_f = ind_i 
            ftyp = i 
         ENDIF 
         ENDDO 
         IF (ftyp.eq.0) then 
            ier_num = - 21 
            ier_typ = ER_FORT 
            cpara (fpara) = ' ' 
            lpara (fpara) = 1 
            RETURN 
         ENDIF 
         ind_i = ind_f 
!                                                                       
!     --We have an integer format                                       
!       The string between '%' and 'd' or 'D' is evaluated and          
!       written into the format variable fstring. Allows error checking 
!       and flexible formats                                            
!                                                                       
         IF (ftyp.eq.FORM_INTEGER.or.ftyp.eq.FORM_INT_ZEROS) then 
            ll = ind_i - 1 
            IF (ll.gt.0) then 
      line = '('//cpara (fpara)  (ind_d+1:ind_d+ind_i - 1) //')' 
               ll = ll + 2 
               CALL rem_bl (line (1:ll), ll) 
               ier_num = 0 
               wert = berechne (line, ll) 
               IF (ier_num.ne.0) then 
                  ier_msg (1) = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
               WRITE (fstring, 3000) int (wert) 
               ll = 15 
               CALL rem_bl (fstring (1:ll), ll) 
            ELSE 
               fstring = '*' 
               ll = 1 
            ENDIF 
!                                                                       
!     --We have a real format                                           
!       The string between '%' and 'f' or 'F' is evaluated and          
!       written into the format variable fstring. Allows error checking 
!       and flexible formats                                            
!                                                                       
         ELSEIF (ftyp.eq.FORM_REAL.or.ftyp.eq.FORM_REAL_ZEROS) then 
            ll = ind_f - 1 
            IF (ll.gt.0) then 
      line = '('//cpara (fpara)  (ind_d+1:ind_d+ind_f - 1) //')' 
               ll = ll + 2 
               CALL rem_bl (line (1:ll), ll) 
               ind_p = INDEX (line (1:ll) , '.') 
               IF (ind_p.gt.0) then 
                  number = line (1:ind_p - 1) //')' 
                  lp = ind_p 
                  itot = nint (berechne (number, lp) ) 
                  IF (ier_num.ne.0) then 
      ier_msg (1)  = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                     cpara (fpara) = ' ' 
                     lpara (fpara) = 1 
                     RETURN 
                  ENDIF 
                  IF (ind_p.lt.ll - 1) then 
                     number = '('//line (ind_p + 1:ll) 
                     lp = ll - ind_p + 1 
                     idec = nint (berechne (number, lp) ) 
                     IF (ier_num.ne.0) then 
      ier_msg (1)  = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                        cpara (fpara) = ' ' 
                        lpara (fpara) = 1 
                        RETURN 
                     ENDIF 
                  ELSE 
                     idec = 0 
                  ENDIF 
               ELSE 
                  itot = nint (berechne (line (1:ll), ll) ) 
                  IF (ier_num.ne.0) then 
      ier_msg (1)  = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                     cpara (fpara) = ' ' 
                     lpara (fpara) = 1 
                     RETURN 
                  ENDIF 
                  idec = 0 
               ENDIF 
               IF (ier_num.ne.0) return 
               WRITE (fstring, 3100) itot, idec 
               ll = 24 
               CALL rem_bl (fstring (1:ll), ll) 
            ELSE 
               fstring = '*' 
               ll = 1 
            ENDIF 
!                                                                       
!     --We have a character format                                      
!       The string between '%' and 'c' or 'C' is evaluated and          
!       written into the format variable fstring. Allows error checking 
!       and flexible formats                                            
!                                                                       
         ELSEIF (ftyp.eq.FORM_CHARACTER.or.ftyp.eq.FORM_CHAR_ZEROS) then
            ll = ind_i - 1 
            IF (ll.gt.0) then 
      line = '('//cpara (fpara)  (ind_d+1:ind_d+ind_i - 1) //')' 
               ll = ll + 2 
               CALL rem_bl (line (1:ll), ll) 
               ier_num = 0 
               wert = berechne (line, ll) 
               IF (ier_num.ne.0) then 
                  ier_msg (1) = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
               WRITE (fstring, 3200) int (wert) 
               fstring_l = int (wert) 
               ll = 15 
               CALL rem_bl (fstring (1:ll), ll) 
            ELSE 
               fstring = 'a' 
               ll = 1 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 21 
         ier_typ = ER_FORT 
         cpara (fpara) = ' ' 
         lpara (fpara) = 1 
         RETURN 
      ENDIF 
!                                                                       
!     --Copy strings left before the current format specifier.          
!                                                                       
      IF (ind_d.gt.pos) then 
         lsub = ind_d-pos 
         string (laenge+1:laenge+lsub) = cpara (fpara) (pos:ind_d-1) 
         laenge = laenge+lsub 
      ENDIF 
!                                                                       
!     --Write numerical value of parameter into string                  
!                                                                       
      IF (ianz.ge.npara) then 
         IF (.not. (ftyp.eq.FORM_CHARACTER.or.ftyp.eq.FORM_CHAR_ZEROS) )&
         then                                                           
!                                                                       
!     ------Calculate value of current parameter                        
!                                                                       
            ll = lpara (npara) 
            line = ' ' 
            line = '('//cpara (npara) (1:ll) //')' 
            ll = ll + 2 
            CALL rem_bl (line (1:ll), ll) 
            wert = berechne (line, ll) 
            IF (ier_num.ne.0) then 
               cpara (fpara) = ' ' 
               lpara (fpara) = 1 
               RETURN 
            ENDIF 
            werte (npara) = wert 
         ELSE 
!                                                                       
!     ------Calculate value of current character parameter              
!           if enclosed in " " take it as is, otherwise calculate       
!                                                                       
            ll = lpara (npara) 
            l_hyp = .false. 
            IF (cpara (npara) (1:1) .eq.'"') then 
               l_hyp = .true. 
               l_sin = .false. 
            ENDIF 
            IF (cpara (npara) (1:1) .eq.'''') then 
               l_hyp = .true. 
               l_sin = .true. 
            ENDIF 
            IF (l_hyp) then 
               IF (l_sin) then 
                  ising = INDEX (cpara (npara) (2:ll) , '''') + 1 
               ELSE 
                  ising = INDEX (cpara (npara) (2:ll) , '"') + 1 
               ENDIF 
               IF (ising.eq.ll) then 
                  cpara (npara) = cpara (npara) (2:ll - 1) 
                  ll = ll - 2 
                  lpara (npara) = ll 
               ELSEIF (ising.gt.1) then 
                  line (1:1) = '(' 
                  line (2:ising - 1) = cpara (npara) (2:ising - 1) 
                  line (ising:ll - 1) = cpara (npara) (ising + 1:ll) 
                  line (ll:ll) = ')' 
                  CALL berechne_char (line, ll) 
                  IF (line (1:1) .eq.''''.and.line (ll:ll) .eq.'''') then
                     cpara (npara) = line (2:ll) 
                     ll = ll - 2 
                  ELSE 
                     cpara (npara) = line (1:ll) 
                  ENDIF 
                  IF (ier_num.ne.0) then 
                     cpara (fpara) = ' ' 
                     lpara (fpara) = 1 
                     RETURN 
                  ENDIF 
               ELSE 
                  ier_num = - 30 
                  ier_typ = ER_FORT 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
            ELSE 
               line = ' ' 
               line = '('//cpara (npara) (1:ll) //')' 
               ll = ll + 2 
               CALL berechne_char (line, ll) 
               IF (line (1:1) .eq.''''.and.line (ll:ll) .eq.'''') then 
                  cpara (npara) = line (2:ll) 
                  ll = ll - 2 
               ELSE 
                  cpara (npara) = line (1:ll) 
               ENDIF 
               lpara (npara) = ll 
               IF (ier_num.ne.0) then 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
         IF (ftyp.eq.FORM_REAL) then 
            pos = ind_d+ind_f + 1 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) werte (npara) 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) werte (npara) 
            ENDIF 
            ii = len_str (number) 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_REAL_ZEROS) then 
            pos = ind_d+ind_f + 1 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) werte (npara) 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) werte (npara) 
            ENDIF 
            ii = len_str (number) 
            DO i = 1, ii 
            IF (number (i:i) .eq.' ') then 
               number (i:i) = '0' 
            ENDIF 
            ENDDO 
            DO i = 1, ii 
            IF (number (i:i) .eq.'-') then 
               number (i:i) = '0' 
               number (1:1) = '-' 
            ENDIF 
            ENDDO 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_INTEGER) then 
            pos = ind_d+ind_i + 1 
            ii = nint (werte (npara) ) 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) ii 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) ii 
            ENDIF 
            ii = len_str (number) 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_INT_ZEROS) then 
            pos = ind_d+ind_i + 1 
            ii = nint (werte (npara) ) 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) ii 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) ii 
            ENDIF 
            ii = len_str (number) 
            DO i = 1, ii 
            IF (number (i:i) .eq.' ') then 
               number (i:i) = '0' 
            ENDIF 
            ENDDO 
            DO i = 1, ii 
            IF (number (i:i) .eq.'-') then 
               number (i:i) = '0' 
               number (1:1) = '-' 
            ENDIF 
            ENDDO 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_CHARACTER) then 
            pos = ind_d+ind_i + 1 
            IF (fstring.eq.'a') then 
               number = cpara (npara) (1:lpara (npara) ) 
               lwert = lpara (npara) 
            ELSE 
               WRITE (number, fstring) cpara (npara) 
               lwert = fstring_l 
            ENDIF 
         ENDIF 
         IF (lwert.eq.0) then 
            ier_num = - 34 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
         string(laenge+1:LEN(STRING)) = ' '
         IF ( (number (1:1) .eq.''''.or.number (1:1) .eq.'"') .and.  &
              (number (lwert:lwert) .eq.''''.or.                     &
               number (lwert:lwert) .eq.'"') ) then
            string (laenge+1:laenge+lwert - 2) = number (2:lwert - 1) 
            laenge = laenge+lwert - 2 
         ELSE 
            string (laenge+1:laenge+lwert) = number (1:lwert) 
            laenge = laenge+lwert 
         ENDIF 
!                                                                       
!     ----Eliminate parameter from list                                 
!                                                                       
         DO i = npara, ianz - 1 
         cpara (i) = cpara (i + 1) 
         lpara (i) = lpara (i + 1) 
         ENDDO 
         cpara (ianz) = ' ' 
         lpara (ianz) = 0 
         ianz = ianz - 1 
      ELSE 
         ier_num = - 27 
         ier_typ = ER_IO 
         cpara (fpara) = ' ' 
         lpara (fpara) = 1 
         RETURN 
      ENDIF 
!                                                                       
!     --Any format specifiers left ?                                    
!                                                                       
      IF (pos.le.ind_qr - 1) then 
         ind_d = INDEX (cpara (fpara) (pos:ind_qr - 1) , '%') + pos - 1 
         lfloat = (cpara (fpara) (ind_d+1:ind_d+1) .eq.'f') 
      ELSE 
         ind_d = 0 
      ENDIF 
      ENDDO 
!                                                                       
!     Append rest of format specifying string to filename               
!                                                                       
      IF (ind_qr - 1.ge.pos) then 
         string (laenge+1:laenge+ind_qr - pos) = &
              cpara (fpara) (pos: ind_qr - 1)
         laenge = laenge+ind_qr - pos 
      ENDIF 
!                                                                       
      cpara (fpara) = string 
      lpara (fpara) = laenge 
!                                                                       
 3000 FORMAT    ('(i',i10,')') 
 3100 FORMAT    ('(f',i10,'.',i10,')') 
 3200 FORMAT    ('(a',i10,')') 
!
END SUBROUTINE do_build_name                  
!
!*****7***************************************************************  
!
END MODULE build_name_mod
