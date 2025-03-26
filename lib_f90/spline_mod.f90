MODULE spline_mod
!
CONTAINS
!
SUBROUTINE spline_prep(nlow,npkt, xpl_in, ypl_in, xmin, xmax, xstep, npkt_equi, xequi, yequi)
!
! Prepare and perform the spline operation on input arrays xpl, ypl
!
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                              , INTENT(IN)  :: nlow
INTEGER                              , INTENT(IN)  :: npkt
REAL(KIND=PREC_DP), DIMENSION(nlow:npkt), INTENT(IN)  :: xpl_in
REAL(KIND=PREC_DP), DIMENSION(nlow:npkt), INTENT(IN)  :: ypl_in
REAL(KIND=PREC_DP)                   , INTENT(IN)  :: xmin
REAL(KIND=PREC_DP)                   , INTENT(IN)  :: xmax
REAL(KIND=PREC_DP)                   , INTENT(IN)  :: xstep
INTEGER                              , INTENT(IN)  :: npkt_equi
REAL(KIND=PREC_DP), DIMENSION(nlow:npkt_equi), INTENT(OUT) :: xequi
REAL(KIND=PREC_DP), DIMENSION(nlow:npkt_equi), INTENT(OUT) :: yequi
!
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: xpl
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: ypl
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: y2a
REAL(KIND=PREC_DP) :: xequ
REAL(KIND=PREC_DP) :: yequ
!
integer :: ipkt
INTEGER :: ii
INTEGER :: all_status
!
ipkt = npkt-nlow+1
!write(*,*) ' SPLINE PREP ', ipkt, nlow, npkt
!npkt_equi =     NINT((xmax-xmin)/xstep) + 1             
if(nint((xmax-xmin)/xstep) +1 /= npkt_equi) then
  ! ERROR MESAGE
   return
endif
ALLOCATE(xpl (1:npkt-nlow+1),stat = all_status) ! Allocate array for calculated powder pattern
ALLOCATE(ypl (1:npkt-nlow+1),stat = all_status) ! Allocate array for calculated powder pattern
ALLOCATE(y2a (1:npkt-nlow+1),stat = all_status) ! Allocate array for calculated powder pattern
!ALLOCATE(xequi(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
!ALLOCATE(yequi(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
xpl(1:ipkt) = xpl_in(nlow:npkt)   !xpl_in(0:npkt-1)
ypl(1:ipkt) = ypl_in(nlow:npkt)   !ypl_in(0:npkt-1)
xequi = 0.0
yequi = 0.0
y2a  = 0.0
!write(*,*) ' PREPARE BD n', npkt+1-nlow, npkt_equi
!write(*,*) ' PREPARE BD x', lbound(xpl), ubound(xpl), xpl(1), xpl(2), xpl(3), maxval(xpl)
!write(*,*) ' PREPARE BD X', lbound(xpl), ubound(xpl), xpl(ipkt-2), xpl(ipkt-1), xpl(ipkt-0), maxval(xpl)
!write(*,*) ' PREPARE BD y', lbound(ypl), ubound(ypl), ypl(1), ypl(2), ypl(3), maxval(ypl)
!write(*,*) ' PREPARE BD 2', lbound(y2a), ubound(y2a)
!write(*,*) ' PREPARE xequ', lbound(xequi), ubound(xequi)
!write(*,*) ' PREPARE yequ', lbound(xequi), ubound(xequi)
!write(*,*) ' PREPARE BD n', npkt-1-nlow
!write(*,*) ' PREPARE xmin', xmin, xstep, xmin + (1 -1)*xstep
!CALL spline (npkt-2-nlow, xpl, ypl, 1.D31, 1.D31, y2a)
CALL spline (ipkt       , xpl, ypl, 1.D31, 1.D31, y2a)
!write(*,*) ' PREPARE BD 2', lbound(y2a), ubound(y2a), maxval(y2a)
DO ii = nlow, npkt_equi
   xequ = xmin + (ii-1)*xstep
!  CALL splint (npkt-2-nlow, xpl, ypl, y2a, xequ, yequ, ier_num)
   CALL splint (ipkt       , xpl, ypl, y2a, xequ, yequ, ier_num)
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
!write(*,*) ' POST xequi  ', xequi(0), xequi(1), xequi(2)
!write(*,*) ' POST yequi  ', yequi(0), yequi(1), yequi(2)
!
DEALLOCATE(y2a, stat = all_status)
deallocate(xpl)
deallocate(ypl)
!
END SUBROUTINE spline_prep
!                                                                       
SUBROUTINE spline (n, x, y, yp1, ypn, y2) 
!
! Prepare the spline y2 values
!
USE precision_mod
!     USE kuplot_config 
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
!
IMPLICIT NONE
!
INTEGER           , INTENT(IN) :: n
REAL(kind=PREC_DP), DIMENSION(n), INTENT(IN) :: x
REAL(kind=PREC_DP), DIMENSION(n), INTENT(IN) :: y
REAL(kind=PREC_DP)              , INTENT(IN) :: yp1
REAL(kind=PREC_DP)              , INTENT(IN) :: ypn
REAL(kind=PREC_DP), DIMENSION(n), INTENT(OUT):: y2
!
!     PARAMETER (nmax = maxarray) 
!     DIMENSION x (n), y (n), y2 (n), u (nmax) 
REAL(kind=PREC_DP), DIMENSION(n) :: u
REAL(kind=PREC_DP)    :: sig, p, qn, un
INTEGER :: i, k
!
!write(*,*) ' Boundaries n', n
!write(*,*) ' BOUNDARIES x', lbound(x), ubound(x), x(1), x(2), x(3)
!write(*,*) ' BOUNDARIES y', lbound(y), ubound(y), y(1), y(2), y(3)
!write(*,*) ' BOUNDARIES 2', lbound(y2), ubound(y2)
!write(*,*) ' yp1, ypn    ', yp1, ypn
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
!write(*,*) ' AFTER FRST ', p, x(3), x(1), y2(1), y2(2), y2(n-1) , maxval(y2)
!write(*,*) ' AFTER u    ',    u(1), u(2),  u(n-2), u(n-1), u(n) , maxval(u)
IF(ypn >  .99e30) THEN 
   qn = 0. 
   un = 0. 
ELSE 
   qn = 0.5 
   un = (3. / (x(n) - x(n-1))) * (ypn - (y(n)-y(n-1)) / (x(n)-x(n-1)))
ENDIF 
y2(n) = (un - qn * u(n-1)) / (qn*y2(n-1) + 1.) 
!write(*,*) ' PRIOR SCND ', qn, un,  x(1), y2(1), y2(2), y2(n-1), y(n), maxval(y2)
scnd: DO k = n - 1, 1, -1 
   y2(k) = y2(k) * y2(k+1) + u(k)
END DO scnd
!write(*,*) ' AFTER SCND ', p, x(3), x(1), y2(1), y2(2), y2(n-1), y2(n), maxval(y2)
!     RETURN 
END SUBROUTINE spline                         
!                                                                       
SUBROUTINE splint (n, xa, ya, y2a, x, y, ier) 
!
use errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
!IMPLICIT integer(i-n)
!IMPLICIT REAL (a-h, o-z)
!      DIMENSION xa (n), ya (n), y2a (n) 
INTEGER           , INTENT(IN)  :: n
REAL(kind=PREC_DP), DIMENSION(n), INTENT(IN)  :: xa
REAL(kind=PREC_DP), DIMENSION(n), INTENT(IN)  :: ya
REAL(kind=PREC_DP)              , INTENT(IN)  :: x
REAL(kind=PREC_DP)              , INTENT(OUT) :: y
REAL(kind=PREC_DP), DIMENSION(n), INTENT(IN ) :: y2a
INTEGER           , INTENT(OUT) :: ier
!
INTEGER :: klo, khi, k
REAL(kind=PREC_DP)    :: h,a, b
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
   ier = -55
   ier_typ = ER_FORT
   ier_msg(1) = 'Internal error in splint '
   ier_msg(2) = 'Please document and report to authors'
!write(*,*) ' ERROR IN SPLINE ', h, xa(khi), xa(klo), khi, klo, n
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
