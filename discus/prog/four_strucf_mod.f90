MODULE four_strucf_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE four_strucf (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      REAL(8) xarg0, xincu, xincv , xincw
      INTEGER h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw, iadd 
!                                                                       
      INTEGER IAND, ISHFT 
!                                                                       
!------ zero fourier array                                              
!                                                                       
      DO i = 1, num (1) * num (2) *num (3)
         tcsf (i) = cmplx (0.0d0, 0.0d0) 
      ENDDO 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
      DO k = 1, nxat 
         xarg0 = xm (1) * xat(k, 1) + xm (2) * xat(k, 2) + xm  (3) * xat(k, 3)
         xincu = uin(1) * xat(k, 1) + uin(2) * xat(k, 2) + uin (3) * xat(k, 3)
         xincv = vin(1) * xat(k, 1) + vin(2) * xat(k, 2) + vin (3) * xat(k, 3)
         xincw = win(1) * xat(k, 1) + win(2) * xat(k, 2) + win (3) * xat(k, 3)
!                                                                       
         iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
         iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
         iincv = nint (64 * I2PI * (xincv - int (xincv) + 1.0d0) ) 
         iincw = nint (64 * I2PI * (xincw - int (xincw) + 1.0d0) ) 
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
                   ii        = ii + 1 
                   tcsf (ii) = tcsf (ii) + cex (IAND  (ISHFT(iarg,-6), MASK) )
                   iarg      = iarg + iincw
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
