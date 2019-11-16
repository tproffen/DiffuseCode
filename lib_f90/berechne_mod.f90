MODULE berechne_mod
!
USE precision_mod
CONTAINS
!
!****7***************************************************************** 
!
REAL(KIND=PREC_DP) FUNCTION berechne (string, laenge) 
!-                                                                      
!     Calculates the value of the expression stored in string           
!+                                                                      
!     USE calc_intr_mod
USE charact_mod
USE do_read_number_mod
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE set_sub_generic_mod
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw =3
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: string 
INTEGER          , INTENT(INOUT) :: laenge
!
CHARACTER(LEN=1024) :: zeile, line, cpara (maxw) 
INTEGER          :: c
INTEGER :: lpara (maxw) 
!INTEGER :: max 
INTEGER :: ikla, iklz, ikla1, ikla2, ikl, ll, lll, ie 
INTEGER :: ikpa, ikpa1, ikpa2, ikp, ikpz, lp, ianz, i, ikom 
INTEGER :: omask ! , nmask might be needed later
REAL(KIND=PREC_DP)    :: werte (maxw) 
REAL(KIND=PREC_DP)    :: r 
LOGICAL  , DIMENSION(1024,0:1) :: lmask
!                                                                       
lmask = .TRUE.
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
berechne = 0D0
!
IF (laenge.eq.0.or.string.eq.' '.or.ier_num.ne.0) then 
   CONTINUE 
ELSE 
   CALL ersetz_variable (string, laenge, lmask, omask) 
   DO ie=2,laenge-1  !while (ie.ne.0)
      IF(string(ie:ie)=='E') THEN
         c = IACHAR(string(ie-1:ie-1))
         IF((zero<=c .and. c<=nine) .AND. (string(ie+1:ie+1)=='+' .OR. string(ie+1:ie+1)=='-')) THEN
            string (ie:ie) = 'e' 
         ENDIF
      ELSEIF(string(ie:ie)=='D') THEN
         c = IACHAR(string(ie-1:ie-1))
         IF((zero<=c .and. c<=nine) .AND. (string(ie+1:ie+1)=='+' .OR. string(ie+1:ie+1)=='-')) THEN
            string (ie:ie) = 'd' 
         ENDIF
      ENDIF
   ENDDO 
   ikla = INDEX (string, '(') 
   DO while (ikla.ne.0) 
      iklz = INDEX (string (ikla + 1:laenge) , ')') + ikla 
      IF (iklz.gt.ikla + 1) then 
         ikla2 = INDEX (string (ikla + 1:iklz) , '(') + ikla 
      ELSE 
         ikla2 = ikla 
      ENDIF 
      ikla1 = ikla 
      DO while (ikla2.lt.iklz.and.ikla2.gt.ikla1) 
         ikla1 = ikla2 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = INDEX (string (ikla1 + 1:iklz) , '(') + ikla1 
         ELSE 
            ikla2 = ikla1 
         ENDIF 
      ENDDO 
      ikl = max (ikla1, ikla2) 
      IF (ikl.ne.0) then 
!                                                                       
!     ----Print error messages if Parentheses are missing               
!                                                                       
         IF (iklz.eq.ikl) then 
            ier_num = -11 
            ier_typ = ER_FORT 
            RETURN 
         ELSEIF (iklz.eq.ikl + 1) then 
            ier_num = -12 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
         line = string (ikl + 1:iklz - 1) 
         ll = iklz - ikl - 1 
!                                                                       
!     ----Found a pair of parentheses, search for inner set             
!                                                                       
         ikpa = INDEX (line, '[') 
         DO while (ikpa.ne.0) 
            ikpz = INDEX (line (ikpa + 1:ll) , ']') + ikpa 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = INDEX (line (ikpa + 1:ikpz) , '[') + ikpa 
            ELSE 
               ikpa2 = ikpa 
            ENDIF 
            ikpa1 = ikpa 
            DO while (ikpa2.lt.ikpz.and.ikpa2.gt.ikpa1) 
               ikpa1 = ikpa2 
               IF (ikpz.gt.ikpa + 1) then 
                  ikpa2 = INDEX (line (ikpa1 + 1:ikpz) , '[') + ikpa1 
               ELSE 
                  ikpa2 = ikpa1 
               ENDIF 
            ENDDO 
            ikp = max (ikpa1, ikpa2) 
            IF (ikp.ne.0) then 
               IF (ikpz.gt.ikp + 1) then 
                  zeile = line (ikp + 1:ikpz - 1) 
                  lp = ikpz - ikp - 1 
                  cpara(:) = ' '
                  lpara(:) = 0
                  werte(:) = 0.0
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     DO i = 1, ianz 
                        CALL eval (cpara (i), lpara (i) ) 
                        IF (ier_num.ne.0) then 
                           RETURN 
                        ENDIF 
                        werte (i) = do_read_number (cpara (i), lpara (i) ) 
                        IF (ier_num.ne.0) then 
                           GOTO 999 
                        ENDIF 
                     ENDDO 
                     CALL p_ersetz_para (ikp, ikpz, line, ll, werte, maxw, ianz)
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                  ELSE 
                     RETURN 
                  ENDIF 
               ELSEIF (ikpz.eq.ikp + 1) then 
                  ier_num = - 10 
                  ier_typ = ER_FORT 
                  RETURN 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ENDIF 
            ikpa = INDEX (line, '[') 
         ENDDO 
!                                                                       
!                                                                       
         ikom = INDEX (line, ',') 
         IF (ikom.eq.0) then 
            CALL eval (line, ll) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF 
         IF (ikl.ge.3) then 
            CALL calc_intr (string, line, ikl, iklz, laenge, ll) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            zeile = ' ' 
            IF (ikl.gt.1) zeile (1:ikl - 1) = string (1:ikl - 1) 
            zeile (ikl:ikl + ll - 1) = line (1:ll) 
            lll = ikl + ll - 1 
            IF (iklz + 1.le.laenge) then 
               zeile (ikl + ll:ikl + ll + laenge-iklz - 1) =  &
                    string (iklz + 1:laenge)
               lll = ikl + ll + laenge-iklz - 1 
            ENDIF 
            string = zeile 
            laenge = lll 
         ENDIF 
      ENDIF 
      ikla = INDEX (string, '(') 
   ENDDO 
ENDIF 
r = do_read_number (string, laenge) 
berechne = do_read_number (string, laenge) 
  999 CONTINUE 
!                                                                       
END FUNCTION berechne                         
!
!****7***************************************************************** 
!
SUBROUTINE berechne_char (string, laenge) 
!-                                                                      
!     Calculates the value of the character expression stored in string 
!+                                                                      
!     USE calc_intr_mod
USE charact_mod
USE do_read_number_mod
USE do_variable_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE set_sub_generic_mod
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw =3
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: string 
INTEGER          , INTENT(INOUT) :: laenge
!                                                                       
CHARACTER(LEN=1024)                  :: zeile, line
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara !(maxw) 
CHARACTER(LEN=1024)                  :: substring 
INTEGER, DIMENSION(MAXW)             :: lpara !(maxw) 
INTEGER :: ikla, iklz, ikla1, ikla2, ikl, ll, lll
INTEGER :: ikpa, ikpa1, ikpa2, ikp, ikpz, lp, ianz, i, ikom 
INTEGER :: lsub 
INTEGER :: icol 
INTEGER :: iapo 
INTEGER, DIMENSION(2) ::  j! (2) 
INTEGER :: omask, nmask   ! Current location in mask
LOGICAL  , DIMENSION(1024,0:1) :: lmask
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!                                                                       
lmask = .TRUE.
!                                                                       
lll = 1
ier_num = 0 
ier_typ = ER_NONE
omask = 0
nmask = 1 
IF (laenge.eq.0.or.string.eq.' '.or.ier_num.ne.0) then 
   CONTINUE 
ELSE 
!
   CALL ersetz_variable (string, laenge, lmask, omask) 
   nmask = MOD(omask+1,2)   ! Shift nmask one index up
!
   ikla = INDEX_MASK (string, '(', lmask(1:LEN_TRIM(string),omask)) 
   parenth: DO while (ikla.ne.0) 
         iklz = INDEX_MASK (string (ikla + 1:laenge) , ')', lmask(ikla + 1:laenge,omask)) + ikla 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = INDEX_MASK (string (ikla + 1:iklz) , '(', lmask(ikla + 1:iklz, omask)) + ikla 
         ELSE 
            ikla2 = ikla 
         ENDIF 
         ikla1 = ikla 
         DO while (ikla2.lt.iklz.and.ikla2.gt.ikla1) 
         ikla1 = ikla2 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = INDEX_MASK (string (ikla1 + 1:iklz) , '(', lmask(ikla1 + 1:iklz, omask)) + ikla1 
         ELSE 
            ikla2 = ikla1 
         ENDIF 
         ENDDO 
         ikl = max (ikla1, ikla2) 
         IF (ikl.ne.0) then 
!                                                                       
!     ----Print error messages if Parentheses are missing               
!                                                                       
            IF (iklz.eq.ikl) then 
               ier_num = - 11 
               ier_typ = ER_FORT 
               RETURN 
            ELSEIF (iklz.eq.ikl + 1) then 
               ier_num = - 12 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            line = string (ikl + 1:iklz - 1) 
            ll = iklz - ikl - 1 
!                                                                       
!     ----Found a pair of parentheses, search for inner set             
!                                                                       
            ikpa = INDEX (line, '[') 
            DO while (ikpa.ne.0) 
            ikpz = INDEX (line (ikpa + 1:ll) , ']') + ikpa 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = INDEX (line (ikpa + 1:ikpz) , '[') + ikpa 
            ELSE 
               ikpa2 = ikpa 
            ENDIF 
            ikpa1 = ikpa 
            DO while (ikpa2.lt.ikpz.and.ikpa2.gt.ikpa1) 
            ikpa1 = ikpa2 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = INDEX (line (ikpa1 + 1:ikpz) , '[') + ikpa1 
            ELSE 
               ikpa2 = ikpa1 
            ENDIF 
            ENDDO 
            ikp = max (ikpa1, ikpa2) 
            IF (ikp.ne.0) then 
               IF (ikpz.gt.ikp + 1) then 
                  zeile = line (ikp + 1:ikpz - 1) 
                  lp = ikpz - ikp - 1 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     DO i = 1, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     CALL p_ersetz_para (ikp, ikpz, line, ll, werte, maxw, ianz)
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                  ELSE 
                     RETURN 
                  ENDIF 
               ELSEIF (ikpz.eq.ikp + 1) then 
                  ier_num = - 10 
                  ier_typ = ER_FORT 
                  RETURN 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ENDIF 
            ikpa = INDEX (line, '[') 
            ENDDO 
!                                                                       
!                                                                       
            ikom = INDEX (line, ',') 
            icol = INDEX (line, ':') 
            iapo = INDEX (line, '''') 
            IF (ikom.eq.0.and.icol.eq.0.and.iapo.eq.0 .AND. ikl<=1) then 
               ier_num = -43
               ier_typ = ER_FORT
               ier_msg(1) = 'Offending string: '//line(1:LEN_TRIM(line))
!              CALL eval (line, ll) 
!              IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
!           ELSEIF(ikom.eq.1.and.icol.eq.0.and.iapo.eq.0) then
!              CALL get_params (line, ianz, cpara, lpara, maxw, ll) 
!              CALL eval(cpara(1),lpara(1))
!              CALL eval(cpara(2),lpara(2))
!              line = cpara(1)(1:lpara(1))//','//cpara(2)(1:lpara(2))
!              ll   = lpara(1)+1+lpara(2)
!           ENDIF 
            IF (ikl.ge.3.and.icol.eq.0) then 
               CALL calc_intr (string, line, ikl, iklz, laenge, ll) 
               lmask(1:len_trim(string),nmask) = .FALSE.   ! NEEDS WORK 
               omask=MOD(omask+1,2)
               nmask=MOD(nmask+1,2)
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
            ELSE 
               IF (icol.ge.1.and.ikl.gt.1) THEN   ! We have a substring
                  IF(INDEX(string,'''')>1 .AND. string(ikl-1:ikl-1)=='''') THEN
                  cpara (1) = line (1:icol - 1) 
                  lpara (1) = icol - 1 
                  cpara (2) = line (icol + 1:ll) 
                  lpara (2) = ll - icol 
                  ianz = 2 
                  DO i = 1, ianz 
                  substring = '('//cpara (i) (1:lpara (i) ) //')' 
                  lsub = lpara (i) + 2 
                  j (i) = nint (berechne (substring, lsub) ) 
                  IF (ier_num.ne.0) then 
                     RETURN 
                  ENDIF 
                  ENDDO 
                  IF (ikl.gt.1) then 
                     IF (ikl.gt.2) then 
                        IF ( (string (2:2) .eq.''''.or.             &
                              string (2:2) .eq.'"'      ) .and.     &
                             (string (ikl - 1:ikl - 1) .eq.''''.or. &
                              string (ikl - 1:ikl - 1) .eq.'"') ) then
                           j (1) = j (1) + 1 
                           j (2) = j (2) + 1 
                        ENDIF 
                     ENDIF 
                     IF (j (1) .le.j(2) .and.1.le.j (1) .and.      &
                         j (1) .lt.ikl  .and.1.le.j (2) .and.      &
                         j (2) .lt.ikl                       ) then    
                        IF (string (1:1) .eq.'(') then 
                           zeile = string (j (1) + 1:j (2) + 1) 
                           lll = j (2) - j (1) + 1 
                        ELSE 
                           zeile = string (j (1) :j (2) ) 
                           lll = j (2) - j (1) + 1 
                        ENDIF 
                     ELSE 
                        ier_num = - 29 
                        ier_typ = ER_FORT 
                        string = ' ' 
                        laenge = 1 
                        RETURN 
                     ENDIF 
                  ENDIF 
                  string = zeile 
                  laenge = lll 
                  ELSE
                     ier_num = -43
                     ier_typ = ER_FORT
                     IF(ikl>2) THEN
                        ier_msg(1) = 'Offending string: '//string(2:ikl-1)
                     ENDIF
                     RETURN
                  ENDIF
               ELSE 
                  zeile = ' ' 
                  IF (ikl.gt.1) zeile (1:ikl - 1) = string (1:ikl - 1) 
                  zeile (ikl:ikl + ll - 1) = line (1:ll) 
                  lll = ikl + ll - 1 
                  IF (iklz + 1.le.laenge) then 
                     zeile (ikl + ll:ikl + ll + laenge-iklz - 1)        &
                     = string (iklz + 1:laenge)                         
                     lll = ikl + ll + laenge-iklz - 1 
                  ENDIF 
                  string = zeile 
                  laenge = lll 
               ENDIF 
            ENDIF 
         ENDIF 
         ikla = INDEX_MASK (string, '(', lmask(1:len_trim(string), omask))
   ENDDO  parenth
ENDIF 
999 CONTINUE 
!                                                                       
END SUBROUTINE berechne_char                  
!
!*****7**************************************************************** 
!
END MODULE berechne_mod
