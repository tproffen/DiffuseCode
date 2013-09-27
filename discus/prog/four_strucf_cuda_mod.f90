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
      USE crystal_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      REAL(8) xarg0, xincu, xincv , xincw
      INTEGER h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw, iadd 
!                                                                       
      INTEGER IAND, ISHFT 
      real cex_r(I2PI),cex_i(I2PI)
      real tcsf_r(MAXQXY),tcsf_i(MAXQXY)
      
      cex_r=real(cex)
      cex_i=aimag(cex)
      
      call cudastrucf(tcsf_r,tcsf_i,cex_r,cex_i,xat)
      
      tcsf=cmplx(tcsf_r,tcsf_i)



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
