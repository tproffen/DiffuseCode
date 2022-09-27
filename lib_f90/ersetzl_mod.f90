MODULE ersetzl_mod
!
private
!
public ersetzl
!
CONTAINS
!
!*****7**************************************************************** 
!
SUBROUTINE ersetzl (string, ikl, iklz, ww, lfunk, lll) 
!
!     Replaces the intrinsic logical function and its argument by the           
!     corresponding value ww
!
      USE blanks_mod
      IMPLICIT none 
!
      CHARACTER (LEN= * ) , INTENT(INOUT) ::string 
      INTEGER             , INTENT(IN)    :: ikl
      INTEGER             , INTENT(IN)    :: iklz
      LOGICAL             , INTENT(IN)    :: ww 
      INTEGER             , INTENT(IN)    :: lfunk
      INTEGER             , INTENT(INOUT) :: lll 
!
      CHARACTER(LEN(STRING)):: zeile 
      INTEGER ltot , laenge
!
      laenge = lll 
      zeile = ' ' 
      IF (ikl.gt.1) zeile (1:ikl - 1 - lfunk) = string (1:ikl - 1 - lfunk)
      WRITE (zeile (ikl - lfunk:ikl - lfunk     ) , '(L1     )') ww 
      lll = ikl - lfunk
      IF (iklz + 1.le.laenge) then 
         ltot = (ikl - lfunk + 1 ) + (laenge-iklz - 1 + 1) - 1 
         IF (ltot.le.len (zeile) ) then 
            zeile (ikl - lfunk + 1 :ltot) = string (iklz + 1:laenge) 
            lll = lll + laenge- (iklz + 1) + 1 
         ENDIF 
      ENDIF 
      string = zeile 
      CALL rem_bl (string, lll) 
END SUBROUTINE ersetzl                        
!
!*****7**************************************************************** 
!
END MODULE ersetzl_mod
