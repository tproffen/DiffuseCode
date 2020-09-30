MODULE four_strucf_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE four_strucf(iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!
!  Apparently a good OMP version. Very good serial version as well.
!  The outer loop is j=0,num(1)-1 and the threadprivate variable 
!  tcsfp is reduced to size (0:num(2)*num(3)-1)
!  To reduce the number of calculations within the Fourier loop, the
!  increments xinc*, oinc*, iinc*, jinc* are calculated once up front.
!
!  This algorithm is roughly twice as fast on an optimized serial compilation
!  compared to the original serial Fourier calculation. Thus its serial form
!  replaces the old serial algorithm as well.
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, jj
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
INTEGER (KIND=PREC_INT_LARGE)   :: num23,num123, num123_1
!
INTEGER :: IAND, ISHFT 
COMPLEX(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP
!
!------ zero fourier array                                              
!                                                                       
tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
num23    =        num(2)*num(3)
num123   = num(1)*num(2)*num(3)
num123_1 = num(1)*num(2)*num(3)-1
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
DO k = 1, nxat 
   xincu(k) = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
   xincv(k) = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
   xincw(k) = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
   oincu(k) = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
   oincv(k) = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
   oincw(k) = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
   iincu(k) = NINT (64 * I2PI * (xincu(k) - INT (xincu(k)) + 0.0d0) ) 
   iincv(k) = NINT (64 * I2PI * (xincv(k) - INT (xincv(k)) + 0.0d0) ) 
   iincw(k) = NINT (64 * I2PI * (xincw(k) - INT (xincw(k)) + 0.0d0) ) 
   jincu(k) = NINT (64 * I2PI * (oincu(k) - INT (oincu(k)) + 0.0d0) ) 
   jincv(k) = NINT (64 * I2PI * (oincv(k) - INT (oincv(k)) + 0.0d0) ) 
   jincw(k) = NINT (64 * I2PI * (oincw(k) - INT (oincw(k)) + 0.0d0) ) 
ENDDO
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
   IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
   ALLOCATE (tcsfp (0:       num(2)*num(3)-1)) !,0:nthreads-1))
   tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
      DO j = 0, num (1) - 1
         jj = j*num(2)*num(3)
         tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
         loop_k: DO k = 1, nxat 
            iarg0 =  lmn(1)*iincu(k) + lmn(2)*iincv(k) + lmn(3)*iincw(k) + &
                     lmn(4)*jincu(k) + lmn(5)*jincv(k) + lmn(6)*jincw(k) + &
                     j*iincu(k)
            ii = 0 
            DO i = 0, num (2) - 1
               iarg = iarg0 +              iincv(k)*i 
               DO h = 1, num (3) 
                  tcsfp(ii)= tcsfp(ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                  iarg     = iarg + iincw(k)
                  ii       = ii + 1 
               ENDDO 
            ENDDO 
         ENDDO  loop_k
         DO ii = 1,        num(2)*num(3)
            tcsf(jj+ii) = tcsf(jj+ii) + tcsfp(ii-1)
         ENDDO
      ENDDO 
   DEALLOCATE(tcsfp)
!
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
   DO  i = 1, num (1) * num (2) * num(3)
!  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
      tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
!  END FORALL
   END DO
ENDIF 
!                                                                       
END SUBROUTINE four_strucf
!
!*****7*****************************************************************
!
      SUBROUTINE four_strucf_original (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)        ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)        ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE)   ::                              jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
!
INTEGER :: IAND, ISHFT 
!
!------ zero fourier array                                              
!                                                                       
tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
DO k = 1, nxat 
!        xarg0 = xm (1)        * xat(k, 1) + xm (2)         * xat(k, 2) + xm (3)         * xat(k, 3)
   xincu = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
   xincv = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
   xincw = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
   oincu = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
   oincv = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
   oincw = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
!        iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 0.0d0) ) 
   iincu = nint (64 * I2PI * (xincu - int (xincu) + 0.0d0) ) 
   iincv = nint (64 * I2PI * (xincv - int (xincv) + 0.0d0) ) 
   iincw = nint (64 * I2PI * (xincw - int (xincw) + 0.0d0) ) 
   jincu = nint (64 * I2PI * (oincu - int (oincu) + 0.0d0) ) 
   jincv = nint (64 * I2PI * (oincv - int (oincv) + 0.0d0) ) 
   jincw = nint (64 * I2PI * (oincw - int (oincw) + 0.0d0) ) 
   iarg0 =  lmn(1)*iincu + lmn(2)*iincv + lmn(3)*iincw + &
            lmn(4)*jincu + lmn(5)*jincv + lmn(6)*jincw
   iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                 iadd      = ISHFT (iarg, - 6) 
!                 iadd      = IAND  (iadd, MASK) 
!                 tcsf (ii) = tcsf (ii) + cex (iadd, MASK) )
!                 iarg      = iarg + iincw
!                                                                       
   ii = 0 
!                                                                       
   DO j = 0, num (1) - 1
      DO i = 0, num (2) - 1
         iarg = iarg0 + iincu*j + iincv*i 
         DO h = 1, num (3) 
            ii       = ii + 1 
            tcsf(ii) = tcsf (ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )
            iarg     = iarg + iincw
         ENDDO 
      ENDDO 
   ENDDO 
ENDDO 
!
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
   DO  i = 1, num (1) * num (2) * num(3)
!  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
      tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
!  END FORALL
   END DO
ENDIF 
!                                                                       
END SUBROUTINE four_strucf_original                    
!
END MODULE four_strucf_mod
