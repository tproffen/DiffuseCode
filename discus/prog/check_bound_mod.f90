MODULE check_bound_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE check_bound (cell, offset, fp, ltype) 
!+                                                                      
!     This routine applies periodic boundaries if 'fp' is TRUE.         
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER cell (3) 
      REAL offset (3) 
      LOGICAL fp (3), ltype, lok (3) 
!                                                                       
      INTEGER i 
      INTEGER j 
!                                                                       
!------ Apply periodic boundaries                                       
!                                                                       
      DO i = 1, 3 
      IF (fp (i) ) then 
         offset (i) = 0.0 
!                                                                       
         IF (cell (i) .lt.1) then 
            IF (cr_icc (i) .eq.1) then 
               offset (i) = REAL(cell (i) - 1) 
               cell (i) = 1 
            ELSE 
               j = ( - cell (i) / cr_icc (i) + 1) * cr_icc (i) 
               cell (i) = cell (i) + j 
               offset (i) = - REAL(j) 
            ENDIF 
         ELSEIF (cell (i) .gt.cr_icc (i) ) then 
            IF (cr_icc (i) .eq.1) then 
               offset (i) = REAL(cell (i) - 1) 
               cell (i) = 1 
            ELSE 
               j = ( (cell (i) - 1) / cr_icc (i) ) * cr_icc (i) 
               cell (i) = cell (i) - j 
               offset (i) = REAL(j) 
            ENDIF 
         ENDIF 
!                                                                       
         lok (i) = .true. 
!                                                                       
!------ NO periodic boundaries                                          
!                                                                       
      ELSE 
         offset (i) = 0.0 
         lok (i) = (cell (i) .gt.0) .and. (cell (i) .le.cr_icc (i) ) 
      ENDIF 
      ENDDO 
!                                                                       
      ltype = lok (1) .and.lok (2) .and.lok (3) 
!                                                                       
      END SUBROUTINE check_bound                    
!*****7*****************************************************************
END MODULE check_bound_mod
