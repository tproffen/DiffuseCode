MODULE tensors_mod
!
!  All routines in tensors have become obsolete
!
private 
CONTAINS
!*****7*****************************************************************
!                                                                       
      SUBROUTINE matmulx (a, b, c) 
!+                                                                      
!     Matrixmultiplication a=b*c                                        
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k 
      REAL a (3, 3), b (3, 3), c (3, 3) 
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      a (i, j) = 0.0 
      DO k = 1, 3 
      a (i, j) = a (i, j) + b (i, k) * c (k, j) 
      ENDDO 
      ENDDO 
      ENDDO 
      END SUBROUTINE matmulx                        
!*****7*****************************************************************
      SUBROUTINE matmul4 (a, b, c) 
!+                                                                      
!     Matrixmultiplication a=b*c for (4x4) matrices.                    
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k 
      REAL a (4, 4), b (4, 4), c (4, 4) 
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      a (i, j) = 0.0 
      DO k = 1, 4 
      a (i, j) = a (i, j) + b (i, k) * c (k, j) 
      ENDDO 
      ENDDO 
      ENDDO 
      END SUBROUTINE matmul4                        
!*****7*****************************************************************
      SUBROUTINE transmat (mat, idim) 
!-                                                                      
!     Replaces a matrix by its transpose                                
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim, i, j 
      REAL mat (idim, idim), d 
!                                                                       
      DO i = 2, idim 
      DO j = 1, i - 1 
      d = mat (i, j) 
      mat (i, j) = mat (j, i) 
      mat (j, i) = d 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE transmat                       
!*****7*****************************************************************
      SUBROUTINE invmat (imat, a) 
!                                                                       
!     calculates the inverse matrix "imat" to input matrix "a"          
!                                                                       
      USE errlist_mod 
!
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL imat (3, 3), a (3, 3), det 
!                                                                       
      det = a (1, 1) * (a (2, 2) * a (3, 3) - a (2, 3) * a (3, 2) )     &
      + a (2, 1) * (a (3, 2) * a (1, 3) - a (1, 2) * a (3, 3) ) + a (3, &
      1) * (a (1, 2) * a (2, 3) - a (1, 3) * a (2, 2) )                 
                                                                        
      IF (abs (det) .gt.0.0) then 
         ier_num = 0 
         ier_typ = ER_NONE 
!                                                                       
         imat (1, 1) = (a (2, 2) * a (3, 3) - a (2, 3) * a (3, 2) )     &
         / det                                                          
         imat (1, 2) = - (a (1, 2) * a (3, 3) - a (1, 3) * a (3, 2) )   &
         / det                                                          
         imat (1, 3) = (a (1, 2) * a (2, 3) - a (1, 3) * a (2, 2) )     &
         / det                                                          
!                                                                       
         imat (2, 1) = - (a (2, 1) * a (3, 3) - a (3, 1) * a (2, 3) )   &
         / det                                                          
         imat (2, 2) = (a (1, 1) * a (3, 3) - a (1, 3) * a (3, 1) )     &
         / det                                                          
         imat (2, 3) = - (a (1, 1) * a (2, 3) - a (1, 3) * a (2, 1) )   &
         / det                                                          
!                                                                       
         imat (3, 1) = (a (2, 1) * a (3, 2) - a (3, 1) * a (2, 2) )     &
         / det                                                          
         imat (3, 2) = - (a (1, 1) * a (3, 2) - a (1, 2) * a (3, 1) )   &
         / det                                                          
         imat (3, 3) = (a (1, 1) * a (2, 2) - a (1, 2) * a (2, 1) )     &
         / det                                                          
!                                                                       
      ELSE 
         ier_num = - 1 
         ier_typ = ER_MATH 
      ENDIF 
!                                                                       
      END SUBROUTINE invmat                         
!*****7*****************************************************************
      SUBROUTINE invmat4 (matrix) 
!-                                                                      
!     inverts a 4*4 Symmetry operation                                  
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL matrix (4, 4) 
      INTEGER i, j 
      REAL a (3, 3), b (3, 3), t (3) 
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      a (i, j) = matrix (i, j) 
      ENDDO 
      t (i) = matrix (i, 4) 
      ENDDO 
!                                                                       
      CALL invmat (b, a) 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, 3 
         matrix (i, 4) = 0.0 
         DO j = 1, 3 
         matrix (i, j) = b (i, j) 
         matrix (i, 4) = matrix (i, 4) - b (i, j) * t (j) 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE invmat4                        
!*****7*****************************************************************
END MODULE tensors_mod
