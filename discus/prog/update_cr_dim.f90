MODULE update_cr_dim_mod
!
CONTAINS
!********************************************************************** 
      SUBROUTINE update_cr_dim 
!-                                                                      
!     Updates the crystal dimensions to the current values              
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j 
!                                                                       
!     Set initial values                                                
!                                                                       
      IF (cr_natoms.gt.0) then 
         DO j = 1, 3 
         cr_dim (j, 1) = cr_pos (j, 1) 
         cr_dim (j, 2) = cr_pos (j, 1) 
         ENDDO 
      ELSE 
         DO j = 1, 3 
         cr_dim (j, 1) = 1.e10 
         cr_dim (j, 1) = - 1.e10 
         ENDDO 
      ENDIF 
!                                                                       
!     Update values from all atoms in crystal                           
!                                                                       
      DO i = 1, cr_natoms 
      DO j = 1, 3 
      cr_dim (j, 1) = amin1 (cr_dim (j, 1), cr_pos (j, i) ) 
      cr_dim (j, 2) = amax1 (cr_dim (j, 2), cr_pos (j, i) ) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE update_cr_dim                  
END MODULE update_cr_dim_mod
