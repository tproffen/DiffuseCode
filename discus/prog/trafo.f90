MODULE trafo_mod
!
CONTAINS
!*****7*****************************************************************
!                                                                       
SUBROUTINE trafo (hkl, u, xc, yc, zc, gmat, fmat, dist, eps, gten, rten)
!+                                                                      
!     Calculates the transformation matrix : crystal <> plotsection     
!     from lattice comstant and hkl of viewing vector. Also             
!     calculates the distance of the plotting section from the origin   
!     of the crystal.                                                   
!-                                                                      
use precision_mod
use matrix_mod
!
IMPLICIT none 
!                                                                       
real(kind=PREC_DP), dimension(3), intent(in) :: hkl
real(kind=PREC_DP), dimension(3), intent(in) :: u
real(kind=PREC_DP), dimension(3), intent(out) :: xc
real(kind=PREC_DP), dimension(3), intent(out) :: yc
real(kind=PREC_DP), dimension(3), intent(out) :: zc
real(kind=PREC_DP), dimension(3,3), intent(out) :: gmat
real(kind=PREC_DP), dimension(3,3), intent(out) :: fmat
real(kind=PREC_DP),                 intent(out) :: dist
real(kind=PREC_DP), dimension(3,3,3), intent(in) :: eps
real(kind=PREC_DP), dimension(3,3), intent(in) :: rten
real(kind=PREC_DP), dimension(3,3), intent(in) :: gten
!
INTEGER  :: inull2 (2), i, j, k, l, m, inull 
real(kind=PREC_DP) :: zcc, xcc
!                                                                       
!     Determine achses for transformation                               
!     Transform reciprocal vector HKL to real space vector ZC           
!                                                                       
!DO i = 1, 3 
!   zc (i) = 0.0 
!   DO j = 1, 3 
!      zc (i) = zc (i) + rten (i, j) * hkl (j) 
!   ENDDO 
!ENDDO 
zc = matmul(rten, hkl)
!                                                                       
!zcc = sqrt (quad (real(zc), real(zc), real(gten)) ) 
zcc = sqrt(dot_product(zc, matmul(gten, zc)))
IF (zcc.eq.0.0) return 
inull = 0 
inull2 (1) = 0 
inull2 (2) = 0 
DO i = 1, 3 
   zc (i) = zc (i) / zcc 
   IF (hkl (i) .eq.0.0) then 
      inull = inull + 1 
      inull2 (inull) = i 
   ENDIF 
   xc (i) = 0.0 
   yc (i) = 0.0 
ENDDO 
!                                                                       
!     Determine a direction as fundamental as possible for XC           
!                                                                       
IF (inull.eq.2.or.inull.eq.1) then 
   xc (inull2 (1) ) = 1.0 
ELSE 
   IF (hkl (1) .eq.hkl (2) ) then 
      xc (1) = hkl (1) 
      xc (2) = - hkl (2) 
   ELSEIF (hkl (1) .eq.hkl (3) ) then 
      xc (1) = hkl (1) 
      xc (3) = - hkl (3) 
   ELSEIF (hkl (2) .eq.hkl (3) ) then 
      xc (2) = hkl (2) 
      xc (3) = - hkl (3) 
   ELSE 
      xc (1) = hkl (2) 
      xc (2) = - hkl (1) 
   ENDIF 
ENDIF 
!xcc = sqrt (quad (real(xc), real(xc), real(gten)) ) 
xcc = sqrt(dot_product(xc, matmul(gten, xc)))
!                                                                       
!     normalize XC                                                      
!                                                                       
DO i = 1, 3 
   xc (i) = xc (i) / xcc 
ENDDO 
!                                                                       
!------ calculate yc as the vector product of zc and xc                 
!                                                                       
DO m = 1, 3 
   yc (m) = 0.0 
   DO j = 1, 3 
      DO k = 1, 3 
         DO l = 1, 3 
            yc (m) = yc (m) + eps (j, k, l) * zc (k) * xc (l) * rten (j, m) 
         ENDDO 
      ENDDO 
   ENDDO 
ENDDO 
!                                                                       
!     calculate indices of gmat                                         
!                                                                       
DO i = 1, 3 
   gmat (i, 1) = xc (i) 
   gmat (i, 2) = yc (i) 
   gmat (i, 3) = zc (i) 
ENDDO 
!                                                                       
!     calculate distance to origin                                      
!                                                                       
!CALL invmat (fmat, gmat) 
call matinv (gmat, fmat) 
dist = fmat(3, 1) * u(1) + fmat(3, 2) * u(2) + fmat(3, 3) * u(3)
!                                                                       
END SUBROUTINE trafo                          
!
!*****7*****************************************************************
!
!OLD SUBROUTINE trans_old (uc, gmat, up, idim) 
!+                                                                      
! IS OBSOLETE AND NO LONGER USED
!     Transforms a point in the crystal space into plot space           
!     and vice versa                                                    
!-                                                                      
!OLD USE precision_mod
!
!OLD IMPLICIT none 
!                                                                       
!OLD integer, intent(in) :: idim
!OLD REAL(KIND=PREC_DP), dimension(idim,idim), intent(in)  :: gmat
!OLD REAL(KIND=PREC_DP), dimension(idim)     , intent(in)  :: uc
!OLD REAL(KIND=PREC_DP), dimension(idim)     , intent(out) :: up
!
!OLD INTEGER i, j!, idim 
!                                                                       
!OLD DO i = 1, idim 
!OLD    up (i) = 0.0 
!OLD    DO j = 1, idim 
!OLD       up (i) = up (i) + gmat (i, j) * uc (j) 
!OLD    ENDDO 
!OLD ENDDO 
!                                                                       
!OLD END SUBROUTINE trans_old
!
END MODULE trafo_mod
