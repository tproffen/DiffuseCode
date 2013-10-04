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
      
      call cudastrucf(tcsf,cex,xat,nxat,num,xm,win,vin,uin,cr_natoms)
      
      
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
