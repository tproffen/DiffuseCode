MODULE four_strucf_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE four_strucf (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE crystal_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      INTEGER i
!                                                                       
      
      call four_strucf_cuda(tcsf,cex,xat,nxat,num,xm,win,vin,uin,cr_natoms)
      
      
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) then 
         !DO  i = 1, num (1) * num (2) * num(3)                                           
!        FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
            !tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                            ! Neder's original code
            DO i = 1, num (1)
                  DO j = 1, num (2)
                        DO k = 1, num (3)
                              tcsf (i,j,k) = tcsf (i,j,k) * cfact (istl (i,j,k), iscat)          ! My declaration
                        ENDDO
                  ENDDO
            ENDDO
!        END FORALL
         !END DO
      ENDIF 
!                                                                       
      END SUBROUTINE four_strucf                    
END MODULE four_strucf_mod
