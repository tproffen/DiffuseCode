MODULE spline_mod
!
CONTAINS
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
