module lib_metric_mod
!-
!  Perform low level metric calculations equivalent to DISCUS : oro_blen; do_bang
!  Local routines here require input of metric/reciprocal metric tensors
!
!  lib_blen ( tensor, vector )
!  lib_bang ( tensor, vector_a, vector_b )
!  lib_tensor ( tensor, (a,b,c), (alpha, beta, gamma))
!  lib_angles() ! calculate aver from 3D data set content
!+
!
private
public lib_blen            ! vector length
public lib_bang            ! Angle between two vectors
public lib_vector_product  ! Vector product
public lib_d2r             ! Conversion direc <=> reciprocal
public lib_tensor          ! Build metric tensors
public lib_eps             ! Nuild epsilon tensor
public lib_angles          ! Calculate angles aver for KUPLOT/FOURIER
!
contains
!
!*******************************************************************************
!
pure function lib_blen(gten, u) result(blen)
!-
!  Calculate length of vactor u
!  direct     space: gten == metric tensor
!  reciprocal space: gten == reciprocal metric tensor
!+
use precision_mod
!
implicit none
!
real(kind=PREC_DP) :: blen
!
real(kind=PREC_DP), dimension(3,3), intent(in) :: gten
real(kind=PREC_DP), dimension(3)  , intent(in) :: u
!
blen = sqrt(dot_product(u, matmul(gten, u)))
!
end function lib_blen
!
!*******************************************************************************
!
pure function lib_bang(gten, u, v) result(bang)
!-
!  Calculate angle between vectors u, v
!  direct     space: gten == metric tensor
!  reciprocal space: gten == reciprocal metric tensor
!+
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: bang
!
real(kind=PREC_DP), dimension(3,3), intent(in) :: gten
real(kind=PREC_DP), dimension(3)  , intent(in) :: u
real(kind=PREC_DP), dimension(3)  , intent(in) :: v
!
real(kind=PREC_DP) :: xx
real(kind=PREC_DP) :: xy
real(kind=PREC_DP) :: yy
real(kind=PREC_DP) :: arg
!
bang = 0.0D0
!
xx = sqrt(dot_product(u, matmul(gten, u)))
xy =      dot_product(u, matmul(gten, v))
yy = sqrt(dot_product(v, matmul(gten, v)))
!
if(xx>0.0D0 .and. yy>0.0) then
   arg = xy / (xx * yy)
   if(abs(arg)<1.0D0) then
      bang = acos(arg)/rad
   elseif(arg>=1.0D0) then
      bang = 0.0D0
   elseif(arg<=-1.0D0) then
      bang = 180.0D0
   endif
else
   bang = -1000.0D0
endif
!
end function lib_bang
!
!*****7*****************************************************************
!
subroutine lib_vector_product(u, v, ww, eps, rten) 
!-                                                                      
!     calculates the VECTORPRODUCT in general triclinic space           
!     with  EPS and RTEN in direct space                                
!     with REPS and GTEN in reciprocal space                            
!+                                                                      
use precision_mod
!
implicit none 
!                                                                       
real(kind=PREC_DP), dimension(3),     intent(in)  :: u
real(kind=PREC_DP), dimension(3),     intent(in)  :: v
real(kind=PREC_DP), dimension(3),     intent(out) :: ww
real(kind=PREC_DP), dimension(3,3,3), intent(in)  :: eps
real(kind=PREC_DP), dimension(3,3),   intent(in)  :: rten
!
integer :: i, j, k, l 
!                                                                       
do i = 1, 3 
   ww (i) = 0.0D0 
   do j = 1, 3 
      do k = 1, 3 
         do l = 1, 3 
            ww (i) = ww (i) + eps (j, k, l) * u (k) * v (l) * rten (j, i) 
         enddo 
      enddo 
   enddo 
enddo 
!                                                                       
end subroutine lib_vector_product             
!
!*****7*****************************************************************
!
subroutine lib_d2r(gten, rten, u, v, w)
!-
! Convert a vector from direct to reciprocal space or vice versa
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: gten   ! Metric tensor
real(kind=PREC_DP), dimension(3,3), intent(in)  :: rten   ! Metric tensor in opposite space
real(kind=PREC_DP), dimension(3)  , intent(out) :: u      ! Input vector
real(kind=PREC_DP), dimension(3)  , intent(out) :: v      ! Output vector
real(kind=PREC_DP), dimension(3)  , intent(out) :: w      ! Output vector normalized
!
real(kind=PREC_DP) :: uu  ! length
real(kind=PREC_DP) :: vv  ! length
!
v = matmul(gten, u)             ! converted vector
uu = sqrt(dot_product(u, matmul(gten, u)))   ! Length input
vv = sqrt(dot_product(v, matmul(rten, v)))   ! Length output
if(uu>0.0D0 .and. vv>0.0D0) then
   w  = v/uu/vv                    ! Normalized result
endif
!
end subroutine lib_d2r             
!
!*****7*****************************************************************
!
subroutine lib_tensor(ten, vec, win) 
!+                                                                      
!     Calculates the metric tensor. Works both for direct and           
!     reciprocal metric tensor.                                         
!-                                                                      
use trig_degree_mod
use precision_mod
!
implicit none 
!                                                                       
integer, parameter :: IDIM = 3 
!                                                                       
real(kind=PREC_DP), dimension(IDIM, IDIM), intent(out) :: ten
real(kind=PREC_DP), dimension(IDIM)      , intent(in)  :: vec
real(kind=PREC_DP), dimension(IDIM)      , intent(in)  :: win
INTEGER :: i, j 
!                                                                       
DO i = 1, IDIM 
   DO j = 1, IDIM 
      IF(i /= j) then 
         ten(i, j) = vec(i) * vec(j) * cosd(win(6 - (i + j) ) ) 
      ELSE 
         ten(i, j) = vec(i) * vec(j) 
      ENDIF 
   ENDDO 
ENDDO 
!
end subroutine lib_tensor                         
!
!*****7*****************************************************************
!
subroutine lib_eps(ten, eps)
!-
!  calculate Epsilon tensor
!+
!
use matrix_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3)  , intent(in ) :: ten   ! Direct/reciprocal tensor
real(kind=PREC_DP), dimension(3,3,3), intent(out) :: eps   ! Direct/reciprocal epsilon tensor
!
real(kind=PREC_DP) :: vol
!
eps = 0.0D0
vol = sqrt(abs(determinant(ten)))
!
eps(1, 2, 3) = vol
eps(2, 3, 1) = vol
eps(3, 1, 2) = vol
eps(1, 3, 2) = -vol
eps(3, 2, 1) = -vol
eps(2, 1, 3) = -vol
!
end subroutine lib_eps
!
!*******************************************************************************
!
SUBROUTINE  lib_angles(ltop, length, &
   angle_vh, ratio_vh, aver_vh, &
                                     angle_ht, ratio_ht, aver_ht, &
                                     angle_tv, ratio_tv, aver_tv, &
   gten, inc, vi, use_coor)
!
!USE diffuse_mod 
!USE metric_mod
!USE output_mod
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
REAL(kind=PREC_DP), DIMENSION(3,3), intent(in) :: gten   ! Direct or reciprocal metric tensor
integer           , dimension(3)  , intent(in) :: inc
REAL(kind=PREC_DP), DIMENSION(3,3), intent(in) :: vi     ! Direct or reciprocal metric tensor
integer           , dimension(3)  , intent(in) :: use_coor  ! Use this index for coordinates
!
!
integer  :: extr_abs = 1
integer  :: extr_ord = 2
integer  :: extr_top = 3
REAL(kind=PREC_DP), DIMENSION(3)             ::  hor
REAL(kind=PREC_DP), DIMENSION(3)             ::  ver
REAL(kind=PREC_DP), DIMENSION(3)             ::  top
!
!     Calculate lengths in Ang-1
!
hor(:) = vi(:,1)
ver(:) = vi(:,2)
top(:) = vi(:,3)
length(:) = 0.0
extr_abs = use_coor(1)
extr_ord = use_coor(2)
extr_top = use_coor(3)
IF(inc(1)>1) length(1) = lib_blen(gten,hor)
IF(inc(2)>1) length(2) = lib_blen(gten,ver)
IF(inc(3)>1) length(3) = lib_blen(gten,top)
!write(*,'(a,3f12.8)') ' Tensor   1 ', gten(:,1)
!write(*,'(a,3f12.8)') ' Tensor   2 ', gten(:,2)
!write(*,'(a,3f12.8)') ' Tensor   3 ', gten(:,3)
!write(*,'(a,3f12.8)') ' vi       1 ', vi  (:,1)
!write(*,'(a,3f12.8)') ' vi       2 ', vi  (:,2)
!write(*,'(a,3f12.8)') ' vi       3 ', vi  (:,3)
!write(*,'(a,3f12.8)') ' Lengths vi ', length
!write(*,'(a,3i3)')    ' Use        ', extr_abs, extr_ord,extr_top
CALL lib_angle(angle_vh, inc(2), inc(1), ver, hor,        &
           length(2), length(1), extr_ord, extr_abs,  &
           ratio_vh, aver_vh, gten)
!write(*,'(a,3f12.5)') ' angle rat aver ', angle_vh,ratio_vh, aver_vh 
!write(*,'(a,3f12.5)') ' angle rat aver ', angle_vh,ratio_vh*(inc(2)-1)/(inc(1)-1)
IF(ltop .AND. inc(3)>1 .AND.length(3)>0.0) THEN
   CALL lib_angle(angle_ht, inc(3), inc(1), hor, top,        &
              length(1), length(3), extr_abs, extr_top,  &
              ratio_ht, aver_ht, gten)
   CALL lib_angle(angle_tv, inc(3), inc(2), top, ver,        &
              length(3), length(2), extr_top, extr_ord,  &
              ratio_tv, aver_tv, gten)
ELSE
   angle_ht = 90.0
   ratio_ht = 1.0
   aver_ht  = 1.0
   angle_tv = 90.0
   ratio_tv = 1.0
   aver_tv  = 1.0
ENDIF
!
END SUBROUTINE  lib_angles
!
!*******************************************************************************
!
SUBROUTINE lib_angle(angle_vh, inc2, inc1, ver, hor , &
                 length2, length1, extr_ord, extr_abs,  &
                 ratio_vh, aver_vh, &
                 gten)
!
!USE metric_mod
use errlist_mod
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
REAL(kind=PREC_DP), DIMENSION(3,3), intent(in) :: gten   ! Direct or reciprocal metric tensor
!
!
angle_vh = 0.0
ratio_vh = 0.0 
aver_vh  = 0.0 
IF( inc1>1 .and. inc2>1 )THEN
   angle_vh = lib_bang (gten, hor, ver) 
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
END SUBROUTINE lib_angle
!
!
!*******************************************************************************
!
end module lib_metric_mod
