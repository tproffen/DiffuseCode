MODULE quad_mod
!
CONTAINS
!*****7*****************************************************************
      REAL FUNCTION quad (h, k, rten) 
!+                                                                      
!           Calculates the scalar product of h and k.                   
!           1/d**2 = h(i)*k(j)*rten(i,j)                                
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: idim = 3
!
      REAL, DIMENSION(3)  , INTENT(IN) :: h
      REAL, DIMENSION(3)  , INTENT(IN) :: k
      REAL, DIMENSION(3,3), INTENT(IN) :: rten
!                                                                       
      INTEGER i, j 
!                                                                       
      quad = 0.0 
      DO i = 1, idim 
         DO j = 1, idim 
            quad = quad+h (i) * k (j) * rten (i, j) 
         ENDDO 
      ENDDO 
      END FUNCTION quad                             
END MODULE quad_mod
