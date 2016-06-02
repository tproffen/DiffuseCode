MODULE four_strucf_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE four_strucf (iscat, lform) 
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
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      REAL(PREC_DP)        :: xarg0, xincu, xincv , xincw
      REAL(PREC_DP)        ::        oincu, oincv , oincw
      INTEGER (KIND=16)   :: h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw
      INTEGER (KIND=16)   ::                              jincu, jincv, jincw
!
      INTEGER IAND, ISHFT 
!
!------ zero fourier array                                              
!                                                                       
      tcsf = cmplx (0.0D0, 0.0D0) 
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
                   tcsf(ii) = tcsf (ii) + cex (IAND  (ISHFT(iarg,-6), MASK) )
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
!        FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
            tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
!        END FORALL
         END DO
      ENDIF 
!                                                                       
      END SUBROUTINE four_strucf                    
END MODULE four_strucf_mod
