MODULE gamma_mod
!
private
!
public gammaq
public func_beta
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
REAL(kind=PREC_DP) FUNCTION gammaq(a, x)
!
! Incomplete gamma function
!
use precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN) :: a
REAL(kind=PREC_DP), INTENT(IN) :: x
!
REAL(kind=PREC_DP) :: gamser
REAL(kind=PREC_DP) :: gamcf
REAL(kind=PREC_DP) :: gln
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
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN) :: a
REAL(kind=PREC_DP), INTENT(IN) :: x
REAL(kind=PREC_DP), INTENT(OUT) :: gamser
REAL(kind=PREC_DP), INTENT(OUT) :: gln
!
REAL(kind=PREC_DP)   , PARAMETER :: EPS = 3.E-7
INTEGER, PARAMETER :: ITMAX = 1000
!
REAL(kind=PREC_DP) :: ap
REAL(kind=PREC_DP) :: sum
REAL(kind=PREC_DP) :: del
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
      gamser = sum*EXP(-x+a*LOG(x)-gln)
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
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN) :: a
REAL(kind=PREC_DP), INTENT(IN) :: x
REAL(kind=PREC_DP), INTENT(OUT) :: gamcf
REAL(kind=PREC_DP), INTENT(OUT) :: gln
!
REAL(kind=PREC_DP)   , PARAMETER :: EPS = 3.E-7
INTEGER, PARAMETER :: ITMAX = 1000
!
REAL(kind=PREC_DP) :: a0
REAL(kind=PREC_DP) :: a1
REAL(kind=PREC_DP) :: b0
REAL(kind=PREC_DP) :: b1
REAL(kind=PREC_DP) :: fac
REAL(kind=PREC_DP) :: g
REAL(kind=PREC_DP) :: gold
REAL(kind=PREC_DP) :: an
REAL(kind=PREC_DP) :: ana
REAL(kind=PREC_DP) :: anf
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
   an  = real(n, kind=PREC_DP)
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
         gamcf = EXP(-x+a*LOG(x)-gln)*g
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
REAL(kind=PREC_DP) FUNCTION gammln(xx)
!
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN) :: xx
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
gammln = REAL(tmp + DLOG(STP*SER), kind=PREC_DP)
!
!
END FUNCTION gammln
!
!*******************************************************************************
!
REAL(kind=PREC_DP) function func_beta(z, w)
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: z
real(kind=PREC_DP), intent(in) :: w
!
func_beta = exp(gammln(z) + gammln(w) - gammln(z+w))
!
end function func_beta
!
!*******************************************************************************
!
END MODULE gamma_mod
