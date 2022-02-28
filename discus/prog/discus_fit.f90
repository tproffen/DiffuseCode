MODULE fit_mod
!
! Module for least squares fits in DISCUS
!
PRIVATE
PUBLIC dis_fit_plane
!
CONTAINS
!
!###############################################################################
!
SUBROUTINE dis_fit_plane(NATOMS, list, hkl, dist)
!
USE crystal_mod
USE errlist_mod
USE matrix_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER                              , INTENT(IN)  :: NATOMS
INTEGER           , DIMENSION(NATOMS), INTENT(IN)  :: list
REAL(kind=PREC_DP), DIMENSION(3)     , INTENT(OUT) :: hkl
REAL(kind=PREC_DP)                   , INTENT(OUT) :: dist
!
INTEGER :: i
REAL(kind=PREC_DP), DIMENSION(3)                :: com
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: points
REAL(kind=PREC_DP)                              :: sum_xx, sum_yy, sum_zz
REAL(kind=PREC_DP)                              :: sum_xy, sum_xz, sum_yz
REAL(kind=PREC_DP)                              :: sum_yx, sum_zx, sum_zy
REAL(kind=PREC_DP)                              :: det_x, det_y, det_z
!
hkl(:) = 0
dist   = 0.0
!
ALLOCATE(points(3, NATOMS))
!
com(:) = 0.0
!
DO i=1,NATOMS
   com(1) = com(1) + cr_pos(1,list(i))
   com(2) = com(2) + cr_pos(2,list(i))
   com(3) = com(3) + cr_pos(3,list(i))
ENDDO
com(:) = com(:)/NATOMS
!
DO i=1,NATOMS
   points(1,i) = cr_pos(1,list(i)) - com(1)
   points(2,i) = cr_pos(2,list(i)) - com(2)
   points(3,i) = cr_pos(3,list(i)) - com(3)
ENDDO
!
! Calculate all cross products
!
sum_xx = 0.0
sum_yy = 0.0
sum_zz = 0.0
sum_xy = 0.0
sum_xz = 0.0
sum_yz = 0.0
DO i=1,NATOMS
  sum_xx = sum_xx + points(1,i) * points(1,i)
  sum_yy = sum_yy + points(2,i) * points(2,i)
  sum_zz = sum_zz + points(3,i) * points(3,i)
  sum_xy = sum_xy + points(1,i) * points(2,i)
  sum_xz = sum_xz + points(1,i) * points(3,i)
  sum_yz = sum_yz + points(2,i) * points(3,i)
ENDDO
sum_yx = sum_xy
sum_zx = sum_xz
sum_zy = sum_yz
det_z = sum_xx*sum_yy - sum_xy*sum_xy
det_y = sum_zz*sum_xx - sum_zx*sum_zx
det_x = sum_yy*sum_zz - sum_yz*sum_yz
IF(    ABS(det_z) >= ABS(det_x) .AND. ABS(det_z) >= ABS(det_y) ) THEN
   hkl(1) = (sum_yz*sum_xy - sum_xz*sum_yy)/det_z
   hkl(2) = (sum_xy*sum_xz - sum_xx*sum_yz)/det_z
   hkl(3) = 1.0
ELSEIF(ABS(det_y) >= ABS(det_z) .AND. ABS(det_y) >= ABS(det_x) ) THEN
   hkl(3) = (sum_xy*sum_zx - sum_zy*sum_xx)/det_y
   hkl(1) = (sum_zx*sum_zy - sum_zz*sum_xy)/det_y
   hkl(2) = 1.0
ELSEIF(ABS(det_x) >= ABS(det_y) .AND. ABS(det_x) >= ABS(det_z) ) THEN
   hkl(2) = (sum_zx*sum_xz - sum_yx*sum_zz)/det_x
   hkl(3) = (sum_yz*sum_yx - sum_yy*sum_zx)/det_x
   hkl(1) = 1.0
ENDIF
!
DEALLOCATE(points)
!
dist = hkl(1)*com(1) + hkl(2)*com(2) + hkl(3)*com(3)
if(dist<0.0) THEN
   hkl(:) = -hkl(:)
   dist   = -dist
ENDIF
!
!
END SUBROUTINE dis_fit_plane
!
!###############################################################################
!
!
!###############################################################################
!
END MODULE fit_mod
