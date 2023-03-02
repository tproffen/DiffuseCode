MODULE math_sup
!-
! some useful general math stuff
!
!+
! gcd ! Greatest common divisor
! cardano ! Solve a cubic equation: 0 = C*x^3 + Q*x^2 + L*x + const ; only real roots
!
INTERFACE gcd
   MODULE PROCEDURE gcd_two, gcd_three
END INTERFACE gcd
!
interface math_solve_poly
   module procedure math_cardano, math_quadratic, math_linear
end interface math_solve_poly
!
private
public gcd
public math_solve_poly
!
CONTAINS
!
!*******************************************************************************
!
INTEGER FUNCTION gcd_two(i1, i2)
!
! Determine the greatest common divisor
!
IMPLICIT NONE
INTEGER, INTENT(IN)  :: i1,i2
!
INTEGER :: j1,j2, k
!
k = 1
IF(IABS(i1) > IABS(i2)) THEN
   j1 = IABS(i1)
   j2 = IABS(i2)
   IF(j2==0) THEN
      gcd_two = j1
      RETURN
   ENDIF
ELSEIF(IABS(i1) < IABS(i2)) THEN
   j1 = IABS(i2)
   j2 = IABS(i1)
   IF(j2==0) THEN
      gcd_two = j1
      RETURN
   ENDIF
ELSE
   j1 = IABS(i1)
   j2 = IABS(i2)
   IF(j2==0) THEN
      gcd_two = 1
      RETURN
   ENDIF
ENDIF
!
search: DO
   k = MOD(j1,j2)
   IF(k == 0) EXIT search
   j1 = j2
   j2 = k
ENDDO search
!
gcd_two = j2
!
END FUNCTION gcd_two
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION gcd_three(i1,i2,i3)
!
! Determine the greatest common divisor. Version for three integers
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: i1
INTEGER, INTENT(IN) :: i2
INTEGER, INTENT(IN) :: i3
INTEGER             :: divisor
!
divisor   = gcd_two(i1, i2)
IF(divisor>1) THEN
   gcd_three = gcd_two(divisor, i3)
ELSE
   gcd_three = 1
   IF(i1==0 .AND. I2==0 .AND. IABS(I3)>0) THEN
      gcd_three = IABS(i3)
   ENDIF
ENDIF
!
END FUNCTION gcd_three
!
!*******************************************************************************
!
subroutine math_cardano(a, b, c, d, n, root)
!
! Determine the real roots of a cubic equation
! A x^3 + B x^2 + C x + D = 0
! with A, B, C, D real and 
! A /= 0
!
!+
!
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP)              , intent(in)  :: a     ! A x^3
real(kind=PREC_DP)              , intent(in)  :: b     ! B x^2
real(kind=PREC_DP)              , intent(in)  :: c     ! C x
real(kind=PREC_DP)              , intent(in)  :: d     ! D
integer                         , intent(out) :: n     ! number of distinct roots
real(kind=PREC_DP), dimension(3), intent(out) :: root
!
real(kind=prec_DP), parameter :: eps = 1.0D-5
real(kind=prec_DP) :: p
real(kind=prec_DP) :: q
real(kind=prec_DP) :: delta
real(kind=prec_DP) :: u
real(kind=prec_DP) :: v
real(kind=prec_DP) :: t1
real(kind=prec_DP) :: t2
real(kind=prec_DP) :: arg
!
!write(*,*) ' a,b,c,d', a,b,c,d
if(abs(a)<eps) then         ! Not a cubic equations, but a quadratic
   call math_quadratic(b, c, d, n, root)
   return
endif
p     = ( 9.0D0*a*c -3.0D0*b**2)/(9.0D0*a**2)
q     = (-9.0D0*a*b*c + 27.0D0*a**2*d + 2.0D0*b**3)/(27.0D0*a**3)
delta = (18.0D0*a*b*c*d - 4.0D0*a*c**3 - 27.0D0*a**2*d**2 + b**2*c**2 - 4.0D0*b**3*d)/(108.D0*a**4)
!write(*,*) ' p, q, d', p,q, delta
!
if(abs(delta)< EPS .and. abs(p)<EPS) then            ! Case A Triple solution
!write(*,*) ' CASE A '
   root(1) = -b/3.0D0
   root(2) = -b/3.0D0
   root(3) = -b/3.0D0
   n       = 1
elseif(abs(delta)< EPS .and. abs(p)>=EPS) then       ! Case B a single and a double solution
!write(*,*) ' CASE B '
   root(1) = (b**3 - 4.0D0*a*b*c + 9.0D0*a**2*d) / (3.0D0*a**2*c - a*b**2)
   root(2) = (b*c  - 9.0D0*a*d) / (6.0D0*a*c - 2.0D0*b**2)
   root(3) = (b*c  - 9.0D0*a*d) / (6.0D0*a*c - 2.0D0*b**2)
   n       = 2
elseif(delta<0.0D0) then                             ! Case C a single solution
!write(*,*) ' CASE C '
   t1   = -q/2. + sqrt(-delta)
   t2   = -q/2. - sqrt(-delta)
   u    =  ((abs(t1))**(1./3.))
   v    =  ((abs(t2))**(1./3.))
   if(t1<0) then
     u     = -u
   endif
   if(t2<0) then
     v     = -v
   endif
   root(1) = -b/a/3. + u + v
   root(2) = root(1)
   root(3) = root(1)
   n       = 1
elseif(delta>0.0D0) then                            ! Case D Three solutions
   arg =               (-0.5D0*q*sqrt(-27.D0/p**3))
!write(*,*) ' ARG 1 ', arg
   if(arg<-1.0D0) arg = -1.0D0
   if(arg> 1.0D0) arg =  1.0D0
   root(1) = -b/a/3.0D0 + sqrt(-4.D0/3.D0*p)*cos(1.D0/3.D0*acos(arg)          )
   root(2) = -b/a/3.0D0 - sqrt(-4.D0/3.D0*p)*cos(1.D0/3.D0*acos(arg)+ pi/3.0D0)
   root(3) = -b/a/3.0D0 - sqrt(-4.D0/3.D0*p)*cos(1.D0/3.D0*acos(arg)- pi/3.0D0)
!write(*,*) ' CASE D ', root
   n       = 3
endif
!
end subroutine math_cardano
!
!*******************************************************************************
!
subroutine math_quadratic(a, b, c, n, root)
!
! Determine the real roots of a quadratic equation
! A x^2 + B x   + C       = 0
! with A, B, C    real and 
! A /= 0
!
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP)              , intent(in)  :: a     ! A x^2
real(kind=PREC_DP)              , intent(in)  :: b     ! B x
real(kind=PREC_DP)              , intent(in)  :: c     ! C  
integer                         , intent(out) :: n     ! number of distinct roots
real(kind=PREC_DP), dimension(3), intent(out) :: root
!
real(kind=prec_DP), parameter :: eps = 1.0D-8
real(kind=prec_DP) :: p
real(kind=prec_DP) :: q
real(kind=prec_DP) :: arg
!
if(abs(a)<eps) then         ! Not a quadratic equation, but a linear
   call math_linear(b, c, n, root)
   return
endif
p = b/a
q = c/a
arg = 0.25*p*p - q
!
if(arg>=0.0D0) then      ! Two solutions
   root(1) = -0.5D0*p + sqrt(arg)
   root(2) = -0.5D0*p - sqrt(arg)
   root(3) = -0.5D0*p + sqrt(arg)
   n       = 2
elseif(arg < 0.0D0) then ! No real solutions
   n       = 0
else                     ! A single solution
   root(1) = -0.5D0*p + sqrt(arg)
   root(2) = -0.5D0*p + sqrt(arg)
   root(3) = -0.5D0*p + sqrt(arg)
   n       = 1
endif
!
end subroutine math_quadratic
!
!*******************************************************************************
!
subroutine math_linear(a, b, n, root)
!-
! Determine the real roots of a linear equation
! A x   + B               = 0
! with A, B       real and 
! A /= 0
!
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP)              , intent(in)  :: a     ! A x
real(kind=PREC_DP)              , intent(in)  :: b     ! B 
integer                         , intent(out) :: n     ! number of distinct roots
real(kind=PREC_DP), dimension(3), intent(out) :: root
!
real(kind=prec_DP), parameter :: eps = 1.0D-8
!
if(abs(a) < EPS) then    ! Constant, no solution
   n = 0
else
   root(1) = -b/a
   root(2) = -b/a
   root(3) = -b/a
   n       = 1
endif
!
end subroutine math_linear
!
!*******************************************************************************
!
END MODULE math_sup
