MODULE lib_random_func
!
! Routines for random number generation
!
CONTAINS
!                                                                       
!*****7*****************************************************************
!                                                                       
RECURSIVE REAL(KIND=KIND(1.0D0)) FUNCTION gasdev (sig) 
!-                                                                      
!     calculates a random number for gaussian distribution of           
!     sigma.                                                            
!+                                                                      
USE random_mod
USE precision_mod
!
IMPLICIT none 
SAVE
!                                                                       
REAL(KIND=PREC_DP), INTENT(IN) ::  sig
!                                                                       
REAL(KIND=PREC_DP) ::  v1, v2, r, fac, gset 
REAL(KIND=PREC_DP) :: r1
!                                                                       
!SAVE gset 
!                                                                       
IF (iset.eq.0) then 
1  CONTINUE
   CALL RANDOM_NUMBER(r1)
   v1 = 2.E0 * r1          - 1. 
   CALL RANDOM_NUMBER(r1)
   v2 = 2.E0 * r1          - 1. 
!  v1 = 2.E0 * ran1 (idum) - 1. 
!  v2 = 2.E0 * ran1 (idum) - 1. 
   r = v1**2 + v2**2 
   IF (r.ge.1.) goto 1 
   fac = sqrt ( - 2.0E0 * log (r) / r) 
   gset = v1 * fac 
   gasdev = v2 * fac 
   iset = 1 
ELSE 
   gasdev = gset 
   iset = 0 
ENDIF 
gasdev = gasdev * sig 
!
END FUNCTION gasdev                           
!*****7*****************************************************************
!                                                                       
REAL(KIND=KIND(1.0D0)) FUNCTION gasskew (sig, skew) 
!-                                                                      
!     calculates a random number for gaussian distribution of           
!     sigma and skewness 
!     skew =  0       : symmetric
!     skew =  0.99999 : skewed with right shoulder
!     skew = -0.99999 : skewed with left  shoulder
!+                                                                      
USE random_mod
USE precision_mod
IMPLICIT none 
!                                                                       
!                                                                       
REAL(KIND=PREC_DP), INTENT(in) :: sig
REAL(KIND=PREC_DP), INTENT(in) :: skew
!
REAL(KIND=PREC_DP)  :: v1, v2, v3
!REAL(KIND=PREC_DP)  :: gasdev
!                                                                       
v1      = gasdev(1.0D0)
v2      = gasdev(1.0D0)
v3      = skew*v1 + SQRT(1.0D0-skew**2)*v2
gasskew = v3 * SIGN(1.0D0 ,v1)*sig
!
END FUNCTION gasskew                           
!*****7*****************************************************************
!
REAL(KIND=KIND(1.0D0)) FUNCTION gaslim(sig,factor)
!
!     Calculates a gaussian distributed number, limited to +- factor*sig
!
USE precision_mod
REAL(KIND=PREC_DP), INTENT(IN) :: sig    ! Sigma of Gaussian distribution
REAL(KIND=PREC_DP), INTENT(IN) :: factor ! iLimit in multiples of sigma
!
REAL(KIND=PREC_DP)    :: x
INTEGER :: counter
!
!REAL(KIND=PREC_DP) :: gasdev
!
x = 0.0D0
counter = 0
main: DO 
   x = gasdev(sig)
   counter = counter + 1
   IF(ABS(x) <= factor*sig) EXIT main
   IF(counter>1000) THEN                ! Prevent infinite loop
      CALL RANDOM_NUMBER(x)
      x = (-2.0D0 + 4.D0*x) * sig   ! place into +-2*sigma interval
      EXIT main
   ENDIF
ENDDO main
gaslim = x
!
END FUNCTION gaslim
!*****7*****************************************************************
REAL(kind=PREC_DP) FUNCTION ran1 (idum)
!
!     kept for backwards compatibility, replaces old ran1 from Numerical recipes
!
use precision_mod
USE random_state_mod
!
INTEGER, INTENT(IN) :: idum   ! not needed, ...
!
REAL(kind=PREC_DP) :: r
!
CALL RANDOM_NUMBER(r)
ran1 = r
!
END FUNCTION ran1                             
!
!*****7*****************************************************************
!     REAL FUNCTION ran1 (idum) 
!                                                                       
!     USE random_state_mod
!     USE times_mod
!
!     IMPLICIT NONE 
!     INTEGER, INTENT(INOUT) :: idum
!                                                                       
!     INTEGER m1, m2, m3, ia1, ia2, ia3, ic1, ic2, ic3 
!     INTEGER                           j 
!     INTEGER iff, idum, ix1, ix2, ix3, j 
!     REAL rm1, rm2, r (97) 
!                                                                       
!     save         m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3,ix1,ix2,ix3,iff     
!     SAVE ix1, ix2, ix3, iff 
!     SAVE r 
!                                                                       
!       dimension r(97)                                                 
!                                                                       
!     PARAMETER (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 = 3.8580247e-6)
!     PARAMETER (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 = 7.4373773e-6)
!     PARAMETER (m3 = 243000, ia3 = 4561, ic3 = 51349)
!                                                                       
!     DATA iff / 0 / 
!                                                                       
!     IF (idum.lt.0.or.iff.eq.0) then 
!        iff = 1 
!        ix1 = mod (ic1 - idum, m1) 
!        ix1 = mod (ia1 * ix1 + ic1, m1) 
!        ix2 = mod (ix1, m2) 
!        ix1 = mod (ia1 * ix1 + ic1, m1) 
!        ix3 = mod (ix1, m3) 
!        DO 11 j = 1, 97 
!           ix1 = mod (ia1 * ix1 + ic1, m1) 
!           ix2 = mod (ia2 * ix2 + ic2, m2) 
!           r (j) = (float (ix1) + float (ix2) * rm2) * rm1 
!  11    END DO 
!        idum = 1 
!     ENDIF 
!     ix1 = mod (ia1 * ix1 + ic1, m1) 
!     ix2 = mod (ia2 * ix2 + ic2, m2) 
!     ix3 = mod (ia3 * ix3 + ic3, m3) 
!     j = 1 + (97 * ix3) / m3 
!     IF (j.gt.97.or.j.lt.1) then 
!
!        User probably provided erroneoous ix1, ix2, ix3, initialize
!
!        CALL  datum_intrinsic ()   !    by getting time since midnight
!        idum = - midnight
!        iff = 1 
!        ix1 = mod (ic1 - idum, m1) 
!        ix1 = mod (ia1 * ix1 + ic1, m1) 
!        ix2 = mod (ix1, m2) 
!        ix1 = mod (ia1 * ix1 + ic1, m1) 
!        ix3 = mod (ix1, m3) 
!        DO j = 1, 97 
!           ix1 = mod (ia1 * ix1 + ic1, m1) 
!           ix2 = mod (ia2 * ix2 + ic2, m2) 
!           r (j) = (float (ix1) + float (ix2) * rm2) * rm1 
!        ENDDO 
!        idum = 1 
!        ix1 = mod (ia1 * ix1 + ic1, m1) 
!        ix2 = mod (ia2 * ix2 + ic2, m2) 
!        ix3 = mod (ia3 * ix3 + ic3, m3) 
!        j = 1 + (97 * ix3) / m3 
!
!         WRITE ( *, * ) j, idum, iff 
!         WRITE ( *, * ) ix1, ix2, ix3 
!         STOP 
!     ENDIF 
!     ran1 = r (j) 
!     r (j) = (float (ix1) + float (ix2) * rm2) * rm1 
!     END FUNCTION ran1                             
!*****7*****************************************************************
!
REAL (kind=PREC_DP) FUNCTION bessj1 (x) 
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP) :: x 
REAL(kind=PREC_DP) :: ax, xx, z 
real(kind=PREC_DP) ::                                                   &
                      p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, r1, r2,   &
      r3, r4, r5, r6, s1, s2, s3, s4, s5, s6, y                         
SAVE p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, r1, r2, r3, r4, r5,  &
      r6, s1, s2, s3, s4, s5, s6                                        
DATA r1, r2, r3, r4, r5, r6 / 72362614232.d0, - 7895059235.d0,          &
      242396853.1d0, - 2972611.439d0, 15704.48260d0, - 30.16036606d0 /, &
      s1, s2, s3, s4, s5, s6 / 144725228442.d0, 2300535178.d0,          &
      18583304.74d0, 99447.43394d0, 376.9991397d0, 1.d0 /               
DATA p1, p2, p3, p4, p5 / 1.d0, .183105d-2, - .3516396496d-4,           &
      .2457520174d-5, - .240337019d-6 /, q1, q2, q3, q4, q5 /           &
      .04687499995d0, - .2002690873d-3, .8449199096d-5, - .88228987d-6, &
      .105787412d-6 /                                                   
!
IF (abs (x) .lt.8.) then 
         y = x**2 
         bessj1 = (x * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * &
         r6) ) ) ) ) / (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y *&
         s6) ) ) ) ))                                                    
ELSE 
         ax = abs (x) 
         z = 8. / ax 
         y = z**2 
         xx = ax - 2.356194491D0 
bessj1 = (sqrt(.636619772D0 / ax) * (cos(xx) * (p1 + y *                  &
         (p2 + y * (p3 + y * (p4 + y * p5) ) ) ) - z * sin (xx) *         &
         (q1 + y * (q2 + y * (q3 + y * (q4 + y * q5) ) ) ) ) * sign (1.D0,&
         x))                                                             
ENDIF 
!
END FUNCTION bessj1                           
!
!*****7*****************************************************************
REAL(KIND=KIND(1.0E0)) FUNCTION poidev (xm, idum) 
!                                                                       
USE wink_mod
USE precision_mod
!
IMPLICIT NONE 
!                                                                       
INTEGER, INTENT(IN) :: idum 
REAL(KIND=PREC_DP), INTENT(IN) :: xm
!
!U    USES gammln,ran1                                                  
REAL(KIND=PREC_DP) :: alxm, em, g, oldm, sq, t, y
!REAL(KIND=PREC_DP) :: gammln
!REAL               ran1
SAVE alxm, g, oldm, sq 
DATA oldm / -1.D0 / 
IF (xm.lt.12.) then 
   IF (xm.ne.oldm) then 
      oldm = xm 
      g = exp ( - xm) 
   ENDIF 
   em = -1.0D0 
   t = 1.0D0 
2  em = em + 1.D0
   t = t * ran1 (idum) 
   IF (t.gt.g) goto 2 
ELSE 
   IF (xm.ne.oldm) then 
      oldm = xm 
      sq = sqrt (2.D0 * xm) 
      alxm = log (xm) 
      g = xm * alxm - gammln (xm + 1.) 
   ENDIF 
1  y = tan (PI * ran1 (idum) ) 
   em = sq * y + xm 
   IF (em.lt.0.) goto 1 
   em = int (em) 
   t = 0.9E0 * (1. + y**2) * exp (em * alxm - gammln (em + 1.E0) - g)
   IF (ran1 (idum) .gt.t) goto 1 
ENDIF 
poidev = em 
!     RETURN 
END FUNCTION poidev                           
!*****7*****************************************************************
REAL(KIND=KIND(1.0E0)) FUNCTION gammln (xx) 
!
USE precision_mod
!
IMPLICIT none 
!                                                                       
REAL(KIND=PREC_DP), INTENT(IN) ::  xx 
!
INTEGER :: j 
REAL(KIND=PREC_DP) :: ser, stp, tmp, x, y, cof (6) 
SAVE cof, stp 
DATA cof, stp / 76.18009172947146d0,   -86.50532032941677d0,   &
                24.01409824083091d0,    -1.231739572450155d0,  &
                 0.1208650973866179d-2, -0.5395239384953d-5,   &
                 2.5066282746310005d0 /
x = xx 
y = x 
tmp = x + 5.5d0 
tmp = (x + 0.5d0) * log (tmp) - tmp 
ser = 1.000000000190015d0 
DO j = 1, 6 
   y = y + 1.d0 
   ser = ser + cof(j) / y 
ENDDO 
gammln = REAL(tmp + LOG(stp * ser / x) , kind=PREC_DP)
!
END FUNCTION gammln                           
!
!******************************************************************************
!
END MODULE lib_random_func
