MODULE spline_mod
!
CONTAINS
!
SUBROUTINE spline_prep(npkt, xpl, ypl, xmin, xmax, xstep, npkt_equi, xequi, yequi)
!
! Prepare and perform the spline operation on input arrays xpl, ypl
!
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                              , INTENT(IN)  :: npkt
REAL(KIND=PREC_SP), DIMENSION(1:npkt), INTENT(IN)  :: xpl
REAL(KIND=PREC_SP), DIMENSION(1:npkt), INTENT(IN)  :: ypl
REAL(KIND=PREC_SP)                   , INTENT(IN)  :: xmin
REAL(KIND=PREC_SP)                   , INTENT(IN)  :: xmax
REAL(KIND=PREC_SP)                   , INTENT(IN)  :: xstep
INTEGER                              , INTENT(IN)  :: npkt_equi
REAL(KIND=PREC_SP), DIMENSION(1:npkt_equi), INTENT(OUT) :: xequi
REAL(KIND=PREC_SP), DIMENSION(1:npkt_equi), INTENT(OUT) :: yequi
!
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: y2a
REAL(KIND=PREC_SP) :: xequ
REAL(KIND=PREC_SP) :: yequ
!
INTEGER :: ii
INTEGER :: all_status
!
!npkt_equi =     NINT((xmax-xmin)/xstep) + 1             
ALLOCATE(y2a (1:npkt),stat = all_status) ! Allocate array for calculated powder pattern
!ALLOCATE(xequi(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
!ALLOCATE(yequi(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
xequi = 0.0
yequi = 0.0
y2a  = 0.0
CALL spline (npkt, xpl, ypl, 1.e31, 1.e31, y2a)
DO ii = 1, npkt_equi
   xequ = xmin + (ii-1)*xstep
   CALL splint (npkt, xpl, ypl, y2a, xequ, yequ, ier_num)
   IF(ier_num/=0) THEN
!      DEALLOCATE( pow_tmp, stat = all_status)
!      DEALLOCATE( xpl, stat = all_status)
!      DEALLOCATE( ypl, stat = all_status)
      DEALLOCATE( y2a, stat = all_status)
!     DEALLOCATE( xequi, stat = all_status)
!     DEALLOCATE( yequi, stat = all_status)
      RETURN
   ENDIF
   xequi(ii) = xequ
   yequi(ii) = yequ
ENDDO
!npkt_equi = npkt_equi
DEALLOCATE(y2a, stat = all_status)
!
END SUBROUTINE spline_prep
!                                                                       
SUBROUTINE spline (n, x, y, yp1, ypn, y2) 
!
! Prepare the spline y2 values
!
!     USE kuplot_config 
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
!
IMPLICIT NONE
!
INTEGER           , INTENT(IN) :: n
REAL, DIMENSION(n), INTENT(IN) :: x
REAL, DIMENSION(n), INTENT(IN) :: y
REAL              , INTENT(IN) :: yp1
REAL              , INTENT(IN) :: ypn
REAL, DIMENSION(n), INTENT(OUT):: y2
!
!     PARAMETER (nmax = maxarray) 
!     DIMENSION x (n), y (n), y2 (n), u (nmax) 
REAL, DIMENSION(n) :: u
REAL    :: sig, p, qn, un
INTEGER :: i, k
!
IF(yp1 >  .99e30) THEN
   y2(1) = 0. 
   u (1) = 0. 
ELSE 
   y2(1) = - 0.5 
   u (1) = (3. / (x(2) - x(1)) ) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1)
ENDIF 
frst: DO i = 2, n - 1 
   sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1)) 
   p   = sig * y2(i-1) + 2. 
   y2(i) = (sig - 1.) / p 
   u (i) = (6. * ((y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1)) / (x(i)-x(i-1))) / (x(i+1)   &
         - x(i-1)) - sig * u(i - 1)) / p                          
ENDDO  frst
IF(ypn >  .99e30) THEN 
   qn = 0. 
   un = 0. 
ELSE 
   qn = 0.5 
   un = (3. / (x(n) - x(n-1))) * (ypn - (y(n)-y(n-1)) / (x(n)-x(n-1)))
ENDIF 
y2(n) = (un - qn * u(n-1)) / (qn*y2(n-1) + 1.) 
scnd: DO k = n - 1, 1, -1 
   y2(k) = y2(k) * y2(k+1) + u(k)
END DO scnd
!     RETURN 
END SUBROUTINE spline                         
!                                                                       
SUBROUTINE splint (n, xa, ya, y2a, x, y, ier) 
!
IMPLICIT NONE
!
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
!      DIMENSION xa (n), ya (n), y2a (n) 
INTEGER           , INTENT(IN)  :: n
REAL, DIMENSION(n), INTENT(IN)  :: xa
REAL, DIMENSION(n), INTENT(IN)  :: ya
REAL              , INTENT(IN)  :: x
REAL              , INTENT(OUT) :: y
REAL, DIMENSION(n), INTENT(OUT) :: y2a
INTEGER           , INTENT(OUT) :: ier
!
INTEGER :: klo, khi, k
REAL    :: h,a, b
!
klo = 1 
khi = n 
!1 IF(khi - klo > 1) THEN 
main: DO
   IF(khi - klo == 1) EXIT main 
   k = (khi + klo) / 2 
   IF(xa(k) > x) THEN 
      khi = k 
   ELSE 
      klo = k 
   ENDIF 
!   GOTO 1 
!ENDIF 
ENDDO main
h = xa (khi) - xa (klo) 
IF(h == 0.) THEN
   ier = -60
   RETURN
ENDIF
a = (xa (khi) - x) / h 
b = (x - xa (klo) ) / h 
y = a * ya(klo) + b * ya(khi) + ((a**3 - a) * y2a(klo)        &
      + (b**3 - b) * y2a (khi) ) * (h**2) / 6.                          
!     RETURN 
END SUBROUTINE splint                         
!
END MODULE spline_mod
