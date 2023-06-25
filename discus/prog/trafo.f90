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
call build_dmat
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
!
!*****7*****************************************************************
!
subroutine build_dmat
!-
!  Determine the transformation matrix D_ij as in eq. (48) in
!  Trueblood et all AC A, (1996), 52, 770-781
!  This is used to transform the anisotropic U_ic to cartesian space
!+
!
use crystal_mod
use metric_mod
!
use precision_mod
use matrix_mod, only:matinv
!
implicit none
! 
real(kind=PREC_DP), dimension(3) :: e1   ! Cartesian e_1
real(kind=PREC_DP), dimension(3) :: e2   ! Cartesian e_2
real(kind=PREC_DP), dimension(3) :: e3   ! Cartesian e_3
real(kind=PREC_DP), dimension(3) :: v    ! Dummy vector
!real(kind=PREC_DP), dimension(3,3) :: cr_dmat    ! Dummy vector
real(kind=PREC_DP) :: vv  ! dummy number
!integer :: i
!
!  Cartesian e2 is parallel b
e2(1) = 0.0D0
e2(2) = 1.0D0/cr_a0(2)
e2(3) = 0.0D0
!
!  Cartesian e3 is parallel c*
e3(1) = 0.0D0
e3(2) = 0.0D0
e3(3) = 1.0D0
v = matmul(cr_rten, e3)
vv = sqrt (skalpro (v, v, cr_gten) )
e3 = v / vv
!
!  Cartesian e1 is vector product e2 x e3
call vekprod (e2, e3, e1, cr_eps, cr_rten)
vv = sqrt (skalpro (e1, e1, cr_gten) )
e1 = e1 / vv
!write(*,'(a, 3(F12.6,2x) )') 'e (1,:): ', e1
!write(*,'(a, 3(F12.6,2x) )') 'e (2,:): ', e2
!write(*,'(a, 3(F12.6,2x) )') 'e (3,:): ', e3
!
v(1) = 1.00D0
v(2) = 0.00D0
v(3) = 0.00D0
cr_dmat(1,1) = skalpro(e1, v, cr_gten) * cr_ar(1)
cr_dmat(2,1) = skalpro(e2, v, cr_gten) * cr_ar(1)
cr_dmat(3,1) = skalpro(e3, v, cr_gten) * cr_ar(1)
!
v(1) = 0.00D0
v(2) = 1.00D0
v(3) = 0.00D0
cr_dmat(1,2) = skalpro(e1, v, cr_gten) * cr_ar(2)
cr_dmat(2,2) = skalpro(e2, v, cr_gten) * cr_ar(2)
cr_dmat(3,2) = skalpro(e3, v, cr_gten) * cr_ar(2)
!
v(1) = 0.00D0
v(2) = 0.00D0
v(3) = 1.00D0
cr_dmat(1,3) = skalpro(e1, v, cr_gten) * cr_ar(3)
cr_dmat(2,3) = skalpro(e2, v, cr_gten) * cr_ar(3)
cr_dmat(3,3) = skalpro(e3, v, cr_gten) * cr_ar(3)
!write(*,'(a, 3(F12.6,2x) )') 'DD(1,:): ', cr_dmat(1,:)
!write(*,'(a, 3(F12.6,2x) )') 'DD(2,:): ', cr_dmat(2,:)
!write(*,'(a, 3(F12.6,2x) )') 'DD(3,:): ', cr_dmat(3,:)
!read(*,*) i

call matinv(cr_dmat, cr_dimat)
!
end subroutine build_dmat
!
!*****7*****************************************************************
!
END MODULE trafo_mod
