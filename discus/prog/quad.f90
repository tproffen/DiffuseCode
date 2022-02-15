MODULE quad_mod
!
! function quad has become obsolete
private
!
CONTAINS
!*****7*****************************************************************
!OLD      REAL kind=PREC_DP) FUNCTION quad (h, k, rten) 
!+                                                                      
!           Calculates the scalar product of h and k.                   
!           1/d**2 = h(i)*k(j)*rten(i,j)                                
!-                                                                      
!OLD      IMPLICIT none 
!                                                                       
!OLD      INTEGER, PARAMETER :: idim = 3
!
!OLD      REAL(KIND=PREC_DP), DIMENSION(3)  , INTENT(IN) :: h
!OLD      REAL(KIND=PREC_DP), DIMENSION(3)  , INTENT(IN) :: k
!OLD      REAL(KIND=PREC_DP), DIMENSION(3,3), INTENT(IN) :: rten
!                                                                       
!OLD      INTEGER i, j 
!                                                                       
!OLD      quad = 0.0 
!OLD      DO i = 1, idim 
!OLD         DO j = 1, idim 
!OLD            quad = quad+h (i) * k (j) * rten (i, j) 
!OLD         ENDDO 
!OLD      ENDDO 
!OLD      END FUNCTION quad                            
END MODULE quad_mod
