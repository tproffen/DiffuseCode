MODULE gamma_mod
!
!*******************************************************************************
!
CONTAINS
!
! ROUTINES for the GAMMA function 
! NumRec Ch 6
!
!
!*******************************************************************************
!
REAL FUNCTION gammaq(a, x)
!
! Incomplete gamma function
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: a
REAL, INTENT(IN) :: x
!
REAL :: gamser
REAL :: gamcf
REAL :: gln
!
gammaq = 0.0
!
IF(x < 0.0 .OR. a <= 0) THEN
   WRITE(*,*) ' COEFFICIENT <= 0! '
   RETURN
ENDIF
!
IF(x < a+1) THEN
   CALL GSER(a, x, gamser, gln)
   gammaq = 1. - gamser
ELSE
   CALL gcf(a, x, gamcf, gln)
   gammaq = gamcf
ENDIF
!
END FUNCTION gammaq
!
!*******************************************************************************
!
SUBROUTINE gser(a, x, gamser, gln)
!
! 
IMPLICIT NONE
!
REAL, INTENT(IN) :: a
REAL, INTENT(IN) :: x
REAL, INTENT(OUT) :: gamser
REAL, INTENT(OUT) :: gln
!
REAL   , PARAMETER :: EPS = 3.E-7
INTEGER, PARAMETER :: ITMAX = 1000
!
REAL :: ap
REAL :: sum
REAL :: del
INTEGER :: n
!
gln = gammln(a)
!
IF(x < 0.0) THEN
   WRITE(*,*) ' COEFFICIENT <= 0! '
   RETURN
ELSEIF(x == 0.0) THEN
   gamser = 0.0
   RETURN
ENDIF
ap  = a
sum = 1./a
del = sum
!
DO n=1,ITMAX
   ap  = ap + 1.
   del = del*x/ap
   sum = sum + del
   IF(ABS(del) < ABS(sum)*EPS) THEN
      gamser = sum*EXP(-x+a*ALOG(x)-gln)
      RETURN
   ENDIF
ENDDO
!
WRITE(*,*) ' a too large, too few iterations '
!
END SUBROUTINE gser
!
!*******************************************************************************
!
SUBROUTINE gcf(a, x, gamcf, gln)
!
! 
IMPLICIT NONE
!
REAL, INTENT(IN) :: a
REAL, INTENT(IN) :: x
REAL, INTENT(OUT) :: gamcf
REAL, INTENT(OUT) :: gln
!
REAL   , PARAMETER :: EPS = 3.E-7
INTEGER, PARAMETER :: ITMAX = 1000
!
REAL :: a0
REAL :: a1
REAL :: b0
REAL :: b1
REAL :: fac
REAL :: g
REAL :: gold
REAL :: an
REAL :: ana
REAL :: anf
INTEGER :: n
!
gln = gammln(a)
!
gold = 0.0
a0   = 1.0
a1   = x
b0   = 0.0
b1   = 1.0
fac  = 1.0
!
DO n=1,itmax
   an  = FLOAT(n)
   ana = an - a
   a0  = (a1+a0*ana)*fac
   b0  = (b1+b0*ana)*fac
   anf = an*fac
   a1  = x*a0 + anf*a1
   b1  = x*b0 + anf*b1
   IF(a1 /= 0.0) THEN
      fac = 1./a1
      g   = b1*fac
      IF(ABS((g-gold)/g) < EPS ) THEN
         gamcf = EXP(-x+a*ALOG(x)-gln)*g
         RETURN
      ENDIF
      gold = g
   ENDIF
ENDDO
!
END SUBROUTINE gcf
!
!*******************************************************************************
!
REAL FUNCTION gammln(xx)
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: xx
!
REAL(KIND=KIND(1.0d0)), DIMENSION(6) :: cof
REAL(KIND=KIND(1.0d0))               :: stp
REAL(KIND=KIND(1.0d0))               :: half
REAL(KIND=KIND(1.0d0))               :: one
REAL(KIND=KIND(1.0d0))               :: fpf
REAL(KIND=KIND(1.0d0))               :: x  
REAL(KIND=KIND(1.0d0))               :: tmp
REAL(KIND=KIND(1.0d0))               :: ser
!
INTEGER :: j
!
DATA cof / 76.180091730D0, -86.505320330D0,  24.01409822d0, &
           -1.231739516D0,   0.120858003D-2, -0.536382D-5  /
DATA stp / 2.50662827465D0 /
DATA half /0.5D0/
DATA one  /1.0D0/
DATA fpf  /5.5D0/
!
x   = xx - one
tmp = x + fpf
tmp = (x+half)*DLOG(tmp)-tmp
ser = one
!
DO j=1,6
   x = x + one
   ser = ser + cof(j)/x
ENDDO
!
gammln = REAL(tmp + DLOG(STP*SER))
!
!
END FUNCTION gammln
!
!*******************************************************************************
!
END MODULE gamma_mod
