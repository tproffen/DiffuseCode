MODULE update_cr_dim_mod
!
CONTAINS
!
!********************************************************************** 
!
SUBROUTINE update_cr_dim 
!-                                                                      
!     Updates the crystal dimensions to the current values              
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
!
IMPLICIT none 
!
INTEGER :: i, j     ! Counter variables
!
!     Set initial values
!
IF(cr_natoms > 0) THEN 
   DO j = 1, 3 
      cr_dim(j, 1) = cr_pos(j, 1) 
      cr_dim(j, 2) = cr_pos(j, 1) 
   ENDDO 
ELSE 
   DO j = 1, 3 
      cr_dim(j, 1) = 1.e10 
      cr_dim(j, 1) = - 1.e10 
   ENDDO 
ENDIF 
!
!     Update values from all atoms in crystal
!
DO i = 1, cr_natoms 
   DO j = 1, 3 
      cr_dim(j, 1) = MIN(cr_dim(j, 1), cr_pos(j, i)) 
      cr_dim(j, 2) = MAX(cr_dim(j, 2), cr_pos(j, i)) 
   ENDDO 
ENDDO 
!
END SUBROUTINE update_cr_dim                  
!
!********************************************************************** 
!
END MODULE update_cr_dim_mod
