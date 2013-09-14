!*****7*****************************************************************
!                                                                       
      REAL FUNCTION gasdev (sig) 
!-                                                                      
!     calculates a random number for gaussian distribution of           
!     sigma.                                                            
!+                                                                      
      USE random_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL sig, v1, v2, r, fac, gset 
      REAL ran1 
!                                                                       
      SAVE gset 
!                                                                       
      IF (iset.eq.0) then 
    1    v1 = 2. * ran1 (idum) - 1. 
         v2 = 2. * ran1 (idum) - 1. 
         r = v1**2 + v2**2 
         IF (r.ge.1.) goto 1 
         fac = sqrt ( - 2.0 * log (r) / r) 
         gset = v1 * fac 
         gasdev = v2 * fac 
         iset = 1 
      ELSE 
         gasdev = gset 
         iset = 0 
      ENDIF 
      gasdev = gasdev * sig 
      END FUNCTION gasdev                           
!*****7*****************************************************************
      REAL FUNCTION ran1 (idum) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER m1, m2, m3, ia1, ia2, ia3, ic1, ic2, ic3 
      INTEGER iff, idum, ix1, ix2, ix3, j 
      REAL rm1, rm2, r (97) 
!                                                                       
!     save         m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3,ix1,ix2,ix3,iff     
      SAVE ix1, ix2, ix3, iff 
      SAVE r 
!                                                                       
!       dimension r(97)                                                 
!                                                                       
      PARAMETER (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 =            &
      3.8580247e-6)                                                     
      PARAMETER (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 =            &
      7.4373773e-6)                                                     
      PARAMETER (m3 = 243000, ia3 = 4561, ic3 = 51349) 
!                                                                       
      DATA iff / 0 / 
!                                                                       
      IF (idum.lt.0.or.iff.eq.0) then 
         iff = 1 
         ix1 = mod (ic1 - idum, m1) 
         ix1 = mod (ia1 * ix1 + ic1, m1) 
         ix2 = mod (ix1, m2) 
         ix1 = mod (ia1 * ix1 + ic1, m1) 
         ix3 = mod (ix1, m3) 
         DO 11 j = 1, 97 
            ix1 = mod (ia1 * ix1 + ic1, m1) 
            ix2 = mod (ia2 * ix2 + ic2, m2) 
            r (j) = (float (ix1) + float (ix2) * rm2) * rm1 
   11    END DO 
         idum = 1 
      ENDIF 
      ix1 = mod (ia1 * ix1 + ic1, m1) 
      ix2 = mod (ia2 * ix2 + ic2, m2) 
      ix3 = mod (ia3 * ix3 + ic3, m3) 
      j = 1 + (97 * ix3) / m3 
      IF (j.gt.97.or.j.lt.1) then 
         WRITE ( *, * ) j, idum, iff 
         WRITE ( *, * ) ix1, ix2, ix3 
         STOP 
      ENDIF 
      ran1 = r (j) 
      r (j) = (float (ix1) + float (ix2) * rm2) * rm1 
      END FUNCTION ran1                             
!*****7*****************************************************************
      REAL FUNCTION bessj1 (x) 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL x 
      REAL ax, xx, z 
      DOUBLEPRECISION p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, r1, r2,   &
      r3, r4, r5, r6, s1, s2, s3, s4, s5, s6, y                         
      SAVE p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, r1, r2, r3, r4, r5,  &
      r6, s1, s2, s3, s4, s5, s6                                        
      DATA r1, r2, r3, r4, r5, r6 / 72362614232.d0, - 7895059235.d0,    &
      242396853.1d0, - 2972611.439d0, 15704.48260d0, - 30.16036606d0 /, &
      s1, s2, s3, s4, s5, s6 / 144725228442.d0, 2300535178.d0,          &
      18583304.74d0, 99447.43394d0, 376.9991397d0, 1.d0 /               
      DATA p1, p2, p3, p4, p5 / 1.d0, .183105d-2, - .3516396496d-4,     &
      .2457520174d-5, - .240337019d-6 /, q1, q2, q3, q4, q5 /           &
      .04687499995d0, - .2002690873d-3, .8449199096d-5, - .88228987d-6, &
      .105787412d-6 /                                                   
      IF (abs (x) .lt.8.) then 
         y = x**2 
         bessj1 = x * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * &
         r6) ) ) ) ) / (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y *&
         s6) ) ) ) )                                                    
      ELSE 
         ax = abs (x) 
         z = 8. / ax 
         y = z**2 
         xx = ax - 2.356194491 
         bessj1 = sqrt (.636619772 / ax) * (cos (xx) * (p1 + y *        &
         (p2 + y * (p3 + y * (p4 + y * p5) ) ) ) - z * sin (xx) *       &
         (q1 + y * (q2 + y * (q3 + y * (q4 + y * q5) ) ) ) ) * sign (1.,&
         x)                                                             
      ENDIF 
      END FUNCTION bessj1                           
!*****7*****************************************************************
      REAL FUNCTION poidev (xm, idum) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idum 
      REAL xm, PI 
      PARAMETER (PI = 3.141592654) 
!U    USES gammln,ran1                                                  
      REAL alxm, em, g, oldm, sq, t, y, gammln, ran1 
      SAVE alxm, g, oldm, sq 
      DATA oldm / - 1. / 
      IF (xm.lt.12.) then 
         IF (xm.ne.oldm) then 
            oldm = xm 
            g = exp ( - xm) 
         ENDIF 
         em = - 1 
         t = 1. 
    2    em = em + 1. 
         t = t * ran1 (idum) 
         IF (t.gt.g) goto 2 
      ELSE 
         IF (xm.ne.oldm) then 
            oldm = xm 
            sq = sqrt (2. * xm) 
            alxm = log (xm) 
            g = xm * alxm - gammln (xm + 1.) 
         ENDIF 
    1    y = tan (PI * ran1 (idum) ) 
         em = sq * y + xm 
         IF (em.lt.0.) goto 1 
         em = int (em) 
         t = 0.9 * (1. + y**2) * exp (em * alxm - gammln (em + 1.)      &
         - g)                                                           
         IF (ran1 (idum) .gt.t) goto 1 
      ENDIF 
      poidev = em 
      RETURN 
      END FUNCTION poidev                           
!*****7*****************************************************************
      REAL FUNCTION gammln (xx) 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx 
      INTEGER j 
      DOUBLEPRECISION ser, stp, tmp, x, y, cof (6) 
      SAVE cof, stp 
      DATA cof, stp / 76.18009172947146d0, - 86.50532032941677d0,       &
      24.01409824083091d0, - 1.231739572450155d0, .1208650973866179d-2, &
      - .5395239384953d-5, 2.5066282746310005d0 /                       
      x = xx 
      y = x 
      tmp = x + 5.5d0 
      tmp = (x + 0.5d0) * log (tmp) - tmp 
      ser = 1.000000000190015d0 
      DO 11 j = 1, 6 
         y = y + 1.d0 
         ser = ser + cof (j) / y 
   11 END DO 
      gammln = tmp + log (stp * ser / x) 
      RETURN 
      END FUNCTION gammln                           
