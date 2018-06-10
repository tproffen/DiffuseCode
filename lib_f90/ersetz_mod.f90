MODULE ersetz_mod
!
CONTAINS
!
!****7***************************************************************** 
!
      SUBROUTINE ersetz2 (string, ikl, iklz, ww, lfunk, lll) 
!                                                                       
!     Replaces the intrinsic function and its argument by the           
!     corresponding value ww                                            
!                                                                       
      USE blanks_mod
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1024) zeile 
      INTEGER ikl, iklz, laenge, lfunk, lll 
      INTEGER ltot 
      REAL ww 
!                                                                       
      laenge = lll 
      zeile = ' ' 
      IF (ikl.gt.1) zeile (1:ikl - 1 - lfunk) = string (1:ikl - 1 - lfunk)
      WRITE (zeile (ikl - lfunk:ikl - lfunk + 14) , '(e15.8e2)') ww 
      zeile (ikl - lfunk + 11:ikl - lfunk + 11) = 'e' 
      lll = ikl - lfunk + 14 
      IF (iklz + 1.le.laenge) then 
         ltot = (ikl - lfunk + 15) + (laenge-iklz - 1 + 1) - 1 
         IF (ltot.le.len (zeile) ) then 
            zeile (ikl - lfunk + 15:ltot) = string (iklz + 1:laenge) 
            lll = lll + laenge- (iklz + 1) + 1 
         ENDIF 
      ENDIF 
      string = zeile 
      CALL rem_bl (string, lll) 
      END SUBROUTINE ersetz2                        
!
!*****7**************************************************************** 
!
SUBROUTINE ersetzc (string, ikl, iklz, result, l_result, lfunk, lll)
!                                                                       
!     Replaces the intrinsic character function and its argument by the 
!     corresponding value line                                          
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*) , INTENT(INOUT) :: string 
INTEGER          , INTENT(IN)    :: ikl
INTEGER          , INTENT(IN)    :: iklz
CHARACTER(LEN=*) , INTENT(IN)    :: result 
INTEGER          , INTENT(IN)    :: l_result 
INTEGER          , INTENT(IN)    :: lfunk
INTEGER          , INTENT(INOUT) :: lll
!
CHARACTER(LEN=1024) :: zeile 
INTEGER             :: laenge
INTEGER             :: ltot 
!                                                                       
laenge = lll 
zeile = ' ' 
IF (ikl.gt.1) zeile (1:ikl - 1 - lfunk) = string (1:ikl - 1 - lfunk)
!
zeile (ikl - lfunk:ikl - lfunk + l_result + 1) = ''''//result(1:l_result)//''''
lll = ikl - lfunk + l_result + 1 
IF (iklz + 1.le.laenge) THEN 
   ltot = (ikl - lfunk + l_result) + (laenge-iklz - 1 + 1) + 1
   IF (ltot.le.len (zeile) ) THEN 
      zeile (ikl - lfunk + l_result+2:ltot) = string (iklz + 1: laenge)
   ENDIF 
   lll = lll + laenge- (iklz + 1) + 1  + 2
ENDIF 
lll = LEN_TRIM(zeile)
IF (zeile (1:1) .eq.'('.and.zeile (lll:lll) .eq.')') THEN 
   string = zeile (2:lll - 1) 
   lll = lll - 2 
ELSE 
   string = zeile 
ENDIF 
!
END SUBROUTINE ersetzc                        
!
!*****7**************************************************************** 
END MODULE ersetz_mod
