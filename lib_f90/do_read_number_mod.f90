MODULE do_read_number_mod
!
CONTAINS
!
!*****7**************************************************************** 
!
REAL FUNCTION do_read_number (string, laenge) 
!-                                                                      
!     Reads a numerical value from the string. This is the VMS version. 
!     In order to recognize erroneous nonmunerical data in the last     
!     character of string, the format is one character longer than the  
!     actual; length of the parameter. Otherwise unix does not          
!     recognize an integer like e.g. "3a" as an error.                  
!                                                                       
!     The VMS version properly detects an error for the (more sensible) 
!     length equal to the actual length of the parameter                
!+                                                                      
      USE errlist_mod 
      USE charact_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: string 
      INTEGER          , INTENT(IN) :: laenge 
!                                                                       
      CHARACTER(LEN=1024) :: line 
      INTEGER :: i 
!                                                                       
      ier_num = - 1 
      ier_typ = ER_FORT 
      line = string 
!                                                                       
      IF (laenge.gt.0) then 
         i = iachar (string (laenge:laenge) ) 
         IF (zero.le.i.and.i.le.nine.or.i.eq.period.or.i.eq.blank1) then
            READ (line, *, end = 999, err = 999) do_read_number
!                                                                       
            ier_num = 0 
            ier_typ = ER_NONE 
         ELSE 
            do_read_number = 0. 
         ENDIF 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
END FUNCTION do_read_number                   
!
!****7***************************************************************** 
!
SUBROUTINE eval (line, ll) 
!-                                                                      
!       evaluates a line that has only the basic arithmetics            
!+                                                                      
USE blanks_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(LEN=1024), INTENT(INOUT) :: line
      INTEGER            , INTENT(INOUT) :: ll 
!
      INTEGER ipot, iz1, iz2, imal, idiv, iverk 
      INTEGER iexpo, iplus, iminus, ipp 
      LOGICAL lreal 
      REAL w1, w2, ww 
!                                                                       
CALL rem_bl(line,ll)
!...........Evaluate the exponentiation  '**'                           
ipot = INDEX (line, '**') 
DO while (ipot.ne.0) 
   CALL get_w1_w2 (w1, w2, line, ipot, iz1, iz2, ll, 2, lreal) 
   IF (ier_num.ne.0) then 
      RETURN 
   ENDIF 
   ww = w1**w2 
   CALL ersetz (line, iz1, iz2, ww, ipot, ll, 2, lreal) 
   ipot = INDEX (line, '**') 
ENDDO 
!...... Multiplication, division                                        
!
      imal = INDEX (line, '*') 
      idiv = INDEX (line, '/') 
      IF (imal.gt.0.and.idiv.gt.0) then 
         iverk = min (imal, idiv) 
      ELSEIF (imal.gt.0.and.idiv.eq.0) then 
         iverk = imal 
      ELSEIF (imal.eq.0.and.idiv.gt.0) then 
         iverk = idiv 
      ELSEIF (imal.eq.0.and.idiv.eq.0) then 
         iverk = 0 
      ENDIF 
!
loop_verk:DO while (iverk.ne.0) 
      CALL get_w1_w2 (w1, w2, line, iverk, iz1, iz2, ll, 1, lreal) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      IF (line (iverk:iverk) .eq.'*') then 
         ww = w1 * w2 
      ELSEIF (line (iverk:iverk) .eq.'/') then 
         IF (w2.ne.0) then 
            IF (lreal) then 
               ww = w1 / w2 
            ELSE 
               ww = float (int (w1) / int (w2) ) 
            ENDIF 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ENDIF 
      CALL ersetz (line, iz1, iz2, ww, iverk, ll, 1, lreal) 
      imal = INDEX (line, '*') 
      idiv = INDEX (line, '/') 
      IF (imal.gt.0.and.idiv.gt.0) then 
         iverk = min (imal, idiv) 
      ELSEIF (imal.gt.0.and.idiv.eq.0) then 
         iverk = imal 
      ELSEIF (imal.eq.0.and.idiv.gt.0) then 
         iverk = idiv 
      ELSEIF (imal.eq.0.and.idiv.eq.0) then 
         iverk = 0 
      ENDIF 
ENDDO  loop_verk
!...... addition,subtraction                                            
!...........search for all '+' that are not part of  xxxE+yyy           
      IF (line (1:1) .eq.'+') then 
         iplus = INDEX (line (2:ll) , '+') 
         IF (iplus.ne.0) iplus = iplus + 1 
      ELSE 
         iplus = INDEX (line, '+') 
      ENDIF 
      IF (iplus.gt.1) then 
         iexpo = INDEX (line, 'e') 
         DO while (iplus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iplus) 
         ipp = INDEX (line (iplus + 1:ll) , '+') 
         IF (ipp.eq.0) then 
            iplus = 0 
         ELSE 
            iplus = iplus + ipp 
            iexpo = iexpo + INDEX (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (line (1:1) .eq.'-') then 
         iminus = INDEX (line (2:ll) , '-') 
         IF (iminus.ne.0) iminus = iminus + 1 
      ELSE 
         iminus = INDEX (line, '-') 
      ENDIF 
      IF (iminus.gt.1) then 
         iexpo = INDEX (line, 'e') 
         DO while (iminus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iminus) 
         ipp = index (line (iminus + 1:ll) , '-') 
         IF (ipp.eq.0) then 
            iminus = 0 
         ELSE 
            iminus = iminus + ipp 
            iexpo = iexpo + INDEX (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (iplus.gt.0.and.iminus.gt.0) then 
         iverk = min (iplus, iminus) 
      ELSEIF (iplus.gt.0.and.iminus.eq.0) then 
         iverk = iplus 
      ELSEIF (iplus.eq.0.and.iminus.gt.0) then 
         iverk = iminus 
      ELSEIF (iplus.eq.0.and.iminus.eq.0) then 
         iverk = 0 
      ENDIF 
      DO while (iverk.gt.1) 
      CALL get_w1_w2 (w1, w2, line, iverk, iz1, iz2, ll, 1, lreal) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      IF (line (iverk:iverk) .eq.'+') then 
         ww = w1 + w2 
      ELSE 
         ww = w1 - w2 
      ENDIF 
      CALL ersetz (line, iz1, iz2, ww, iverk, ll, 1, lreal) 
      IF (line (1:1) .eq.'+') then 
         iplus = INDEX (line (2:ll) , '+') 
         IF (iplus.ne.0) iplus = iplus + 1 
      ELSE 
         iplus = INDEX (line, '+') 
      ENDIF 
      IF (iplus.gt.1) then 
         iexpo = INDEX (line, 'e') 
         DO while (iplus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iplus) 
         ipp = INDEX (line (iplus + 1:ll) , '+') 
         IF (ipp.eq.0) then 
            iplus = 0 
         ELSE 
            iplus = iplus + ipp 
            iexpo = iexpo + INDEX (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (line (1:1) .eq.'-') then 
         iminus = INDEX (line (2:ll) , '-') 
         IF (iminus.ne.0) iminus = iminus + 1 
      ELSE 
         iminus = INDEX (line, '-') 
      ENDIF 
      IF (iminus.gt.1) then 
         iexpo = INDEX (line, 'e') 
         DO while (iminus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iminus) 
         ipp = INDEX (line (iminus + 1:ll) , '-') 
         IF (ipp.eq.0) then 
            iminus = 0 
         ELSE 
            iminus = iminus + ipp 
            iexpo = iexpo + INDEX (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (iplus.gt.0.and.iminus.gt.0) then 
         iverk = min (iplus, iminus) 
      ELSEIF (iplus.gt.0.and.iminus.eq.0) then 
         iverk = iplus 
      ELSEIF (iplus.eq.0.and.iminus.gt.0) then 
         iverk = iminus 
      ELSEIF (iplus.eq.0.and.iminus.eq.0) then 
         iverk = 0 
      ENDIF 
      ENDDO 
END SUBROUTINE eval                           
!
!*****7**************************************************************** 
!
SUBROUTINE get_w1_w2 (w1, w2, line, iverk, iz1, iz2, ll, lverk, lreal)
!                                                                       
!     Get the two numbers in front and after an operator                
!                                                                       
      USE errlist_mod 
      USE search_string_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL             , INTENT(OUT)   :: w1, w2 
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      CHARACTER(1024) zeile 
      INTEGER          , INTENT(IN)    :: iverk
      INTEGER          , INTENT(OUT)   :: iz1
      INTEGER          , INTENT(OUT)   :: iz2
      INTEGER          , INTENT(IN)    :: ll
      INTEGER          , INTENT(IN)    :: lverk
!     INTEGER suche_vor2, suche_nach2 
      LOGICAL lreal 
!
INTEGER :: lll
!                                                                       
!                                                                       
! Catch a minus/plus that is intended as sign in front of a number
IF(iverk> 1) THEN
   IF(    line(iverk-1:iverk)=='--') THEN
      line(iverk-1:iverk) = '0+'
   ELSEIF(line(iverk-1:iverk)=='+-') THEN
      line(iverk-1:iverk) = '0-'
   ELSEIF(line(iverk-1:iverk)=='-+') THEN
      line(iverk-1:iverk) = '0-'
   ELSEIF(line(iverk-1:iverk)=='++') THEN
      line(iverk-1:iverk) = '0+'
   ENDIF
ENDIF
      lll = iverk - 1 
      iz1 = suche_vor2 (line (1:iverk - 1), lll) 
      zeile = line (iz1:iverk - 1) 
      w1 = do_read_number (zeile, iverk - iz1) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      lreal = INDEX (line (iz1:iverk - 1) , '.') .gt.0 
      lll = ll - (iverk + lverk) + 1 
! Catch a minus/plus that is intended as sign in front of a number
IF(lll>2) THEN 
   IF(    line(iverk+lverk:iverk+lverk+1)=='--') THEN
      line(iverk+lverk:iverk+lverk+1) = '  '
   ELSEIF(line(iverk+lverk:iverk+lverk+1)=='-+') THEN
      line(iverk+lverk:iverk+lverk+1) = ' -'
   ELSEIF(line(iverk+lverk:iverk+lverk+1)=='+-') THEN
      line(iverk+lverk:iverk+lverk+1) = ' -'
   ELSEIF(line(iverk+lverk:iverk+lverk+1)=='++') THEN
      line(iverk+lverk:iverk+lverk+1) = '  '
   ENDIF
ENDIF
      iz2 = suche_nach2 (line (iverk + lverk:ll), lll) 
      zeile = line (iverk + lverk:iverk + lverk + iz2 - 1) 
      w2 = do_read_number (zeile, iz2) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      lreal = lreal.or.INDEX (line (iverk + lverk:iverk + lverk + iz2 - &
      1) , '.') .gt.0                                                   
!                                                                       
      END SUBROUTINE get_w1_w2                      
!
!****7***************************************************************** 
!
SUBROUTINE ersetz (line, iz1, iz2, ww, iverk, ll, lverk, lreal) 
!+                                                                      
!     Replaces the corresponding part of line by the number ww          
!     Different format is used for real and integer numbers             
!                                                                       
      USE blanks_mod
USE precision_mod
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER          , INTENT(IN)    :: iz1
      INTEGER          , INTENT(IN)    :: iz2
      REAL             , INTENT(IN)    :: ww 
      INTEGER          , INTENT(IN)    :: iverk
      INTEGER          , INTENT(INOUT) :: ll
      INTEGER          , INTENT(IN)    :: lverk
      LOGICAL          , INTENT(IN)    :: lreal
      CHARACTER(5) form 
      CHARACTER(1024) zeile 
      INTEGER lw 
      INTEGER ltot 
!                                                                       
      zeile = ' ' 
      IF (iz1.gt.1) zeile (1:iz1 - 1) = line (1:iz1 - 1) 
      IF (lreal) then 
         lw = PREC_WIDTH
         WRITE (zeile (iz1:iz1 + lw - 1) , PREC_F_REAL) ww 
         zeile (iz1 + PREC_MANTIS:iz1 + PREC_MANTIS) = 'e' 
      ELSE 
         lw = int (alog (abs (ww) + 1.) ) + 2 
         WRITE (form, 1000) lw 
         WRITE (zeile (iz1:iz1 + lw - 1), form) int (ww) 
      ENDIF 
      IF (iverk + lverk + iz2.le.ll) then 
         ltot = (iz1 + lw) + (ll - iverk - lverk - iz2 + 1) - 1 
         IF (ltot.le.len (zeile) ) then 
            zeile (iz1 + lw:ltot) = line (iverk + lverk + iz2:ll) 
         ENDIF 
      ENDIF 
      line = zeile 
      ll = ll + lw - (lverk + iz2 + iverk - iz1 - 1) 
      CALL rem_bl (line, ll) 
!                                                                       
 1000 FORMAT    ('(i',i2,')') 
      END SUBROUTINE ersetz                         
!
!*****7**************************************************************** 
!
!
!****7***************************************************************** 
!
END MODULE do_read_number_mod
