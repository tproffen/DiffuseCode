MODULE four_angles_mod
!
USE errlist_mod
!
CONTAINS
!
SUBROUTINE four_angles(ltop, length, angle_vh, ratio_vh, aver_vh, &
                                     angle_ht, ratio_ht, aver_ht, &
                                     angle_tv, ratio_tv, aver_tv)
!
USE diffuse_mod 
USE metric_mod
USE output_mod
use precision_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: ltop
REAL(KIND=PREC_DP), DIMENSION(3), INTENT(INOUT)::  length
REAL(KIND=PREC_DP)            , INTENT(OUT)::  angle_vh
REAL(KIND=PREC_DP)            , INTENT(OUT)::  ratio_vh
REAL(KIND=PREC_DP)            , INTENT(OUT)::   aver_vh
REAL(KIND=PREC_DP)            , INTENT(OUT)::  angle_ht
REAL(KIND=PREC_DP)            , INTENT(OUT)::  ratio_ht
REAL(KIND=PREC_DP)            , INTENT(OUT)::   aver_ht
REAL(KIND=PREC_DP)            , INTENT(OUT)::  angle_tv
REAL(KIND=PREC_DP)            , INTENT(OUT)::  ratio_tv
REAL(KIND=PREC_DP)            , INTENT(OUT)::   aver_tv
LOGICAL, PARAMETER :: lspace = .FALSE.
!
REAL(kind=PREC_DP), DIMENSION(3)             ::  hor
REAL(kind=PREC_DP), DIMENSION(3)             ::  ver
REAL(kind=PREC_DP), DIMENSION(3)             ::  top
REAL(kind=PREC_DP), DIMENSION(3)             ::  zero = (/0.0, 0.0, 0.0/)
!
!     Calculate lengths in Ang-1
!
hor(:) = vi(:,1)
ver(:) = vi(:,2)
top(:) = vi(:,3)
length(:) = 0.0
IF(inc(1)>1) length(1) = do_blen(lspace,hor,zero)
IF(inc(2)>1) length(2) = do_blen(lspace,ver,zero)
IF(inc(3)>1) length(3) = do_blen(lspace,top,zero)
CALL angle(angle_vh, inc(2), inc(1), ver, hor,        &
           length(2), length(1), extr_ord, extr_abs,  &
           ratio_vh, aver_vh)
IF(ltop .AND. inc(3)>1 .AND.length(3)>0.0) THEN
   CALL angle(angle_ht, inc(3), inc(1), hor, top,        &
              length(1), length(3), extr_abs, extr_top,  &
              ratio_ht, aver_ht)
   CALL angle(angle_tv, inc(3), inc(2), top, ver,        &
              length(3), length(2), extr_top, extr_ord,  &
              ratio_tv, aver_tv)
ELSE
   angle_ht = 90.0
   ratio_ht = 1.0
   aver_ht  = 1.0
   angle_tv = 90.0
   ratio_tv = 1.0
   aver_tv  = 1.0
ENDIF
!
END SUBROUTINE four_angles
!
!*******************************************************************************
!
SUBROUTINE angle(angle_vh, inc2, inc1, ver, hor , &
                 length2, length1, extr_ord, extr_abs,  &
                 ratio_vh, aver_vh)
!
USE metric_mod
use precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(OUT):: angle_vh
INTEGER           , INTENT(IN) :: inc1
INTEGER           , INTENT(IN) :: inc2
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN) :: hor
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN) :: ver
REAL(kind=PREC_DP), INTENT(IN) :: length1
REAL(kind=PREC_DP), INTENT(IN) :: length2
INTEGER           , INTENT(IN) :: extr_abs
INTEGER           , INTENT(IN) :: extr_ord
REAL(kind=PREC_DP), INTENT(OUT):: ratio_vh
REAL(kind=PREC_DP), INTENT(OUT):: aver_vh
!
REAL(kind=PREC_DP)   , DIMENSION(3) :: zero = (/0.0, 0.0, 0.0/)
LOGICAL, PARAMETER    :: lspace =.false.
!
angle_vh = 0.0
ratio_vh = 0.0 
aver_vh  = 0.0 
IF( inc1>1 .and. inc2>1 )THEN
   angle_vh = do_bang (lspace, hor, zero, ver) 
   ratio_vh = length2 / length1 
   IF (abs (hor(extr_abs)) > 0.0  .and. abs (ver(extr_ord)) >   0.0) then                                                  
      aver_vh = (length2 / ver(extr_ord) ) / (length1 / hor(extr_abs) ) 
   ELSE 
      ier_num = - 4 
      ier_typ = ER_FOUR 
   ENDIF 
ELSE 
ENDIF 
!
END SUBROUTINE angle
!
END MODULE four_angles_mod
