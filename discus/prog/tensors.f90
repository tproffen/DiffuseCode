MODULE tensors_mod
!
!  All routines in tensors have become obsolete
!
private 
CONTAINS
!*****7*****************************************************************
!                                                                       
!OLD  SUBROUTINE matmulx (a, b, c) 
!+                                                                      
!     Matrixmultiplication a=b*c                                        
!-                                                                      
!OLD  IMPLICIT none 
!                                                                       
!OLD  INTEGER i, j, k 
!OLD  REAL a (3, 3), b (3, 3), c (3, 3) 
!                                                                       
!OLD  DO i = 1, 3 
!OLD  DO j = 1, 3 
!OLD  a (i, j) = 0.0 
!OLD  DO k = 1, 3 
!OLD  a (i, j) = a (i, j) + b (i, k) * c (k, j) 
!OLD  ENDDO 
!OLD  ENDDO 
!OLD  ENDDO 
!OLD  END SUBROUTINE matmulx                        
!*****7*****************************************************************
!OLD  SUBROUTINE matmul4 (a, b, c) 
!+                                                                      
!     Matrixmultiplication a=b*c for (4x4) matrices.                    
!-                                                                      
!OLD  IMPLICIT none 
!                                                                       
!OLD  INTEGER i, j, k 
!OLD  REAL a (4, 4), b (4, 4), c (4, 4) 
!                                                                       
!OLD  DO i = 1, 4 
!OLD  DO j = 1, 4 
!OLD  a (i, j) = 0.0 
!OLD  DO k = 1, 4 
!OLD  a (i, j) = a (i, j) + b (i, k) * c (k, j) 
!OLD  ENDDO 
!OLD  ENDDO 
!OLD  ENDDO 
!OLD  END SUBROUTINE matmul4                        
!*****7*****************************************************************
!OLD  SUBROUTINE transmat (mat, idim) 
!-                                                                      
!     Replaces a matrix by its transpose                                
!+                                                                      
!OLD  IMPLICIT none 
!                                                                       
!OLD  INTEGER idim, i, j 
!OLD  REAL mat (idim, idim), d 
!                                                                       
!OLD  DO i = 2, idim 
!OLD  DO j = 1, i - 1 
!OLD  d = mat (i, j) 
!OLD  mat (i, j) = mat (j, i) 
!OLD  mat (j, i) = d 
!OLD  ENDDO 
!OLD  ENDDO 
!                                                                       
!OLDD END SUBROUTINE transmat                       
!*****7*****************************************************************
!OLD  SUBROUTINE invmat (imat, a) 
!                                                                       
!     calculates the inverse matrix "imat" to input matrix "a"          
!                                                                       
!OLD  USE errlist_mod 
!
!OLD  IMPLICIT none 
!                                                                       
!                                                                       
!OLD  REAL imat (3, 3), a (3, 3), det 
!                                                                       
!OLD  det = a (1, 1) * (a (2, 2) * a (3, 3) - a (2, 3) * a (3, 2) )     &
!OLD  + a (2, 1) * (a (3, 2) * a (1, 3) - a (1, 2) * a (3, 3) ) + a (3, &
!OLD  1) * (a (1, 2) * a (2, 3) - a (1, 3) * a (2, 2) )                 
!OLD                                                                    
!OLD  IF (abs (det) .gt.0.0) then 
!OLD     ier_num = 0 
!OLD     ier_typ = ER_NONE 
!                                                                       
!OLD     imat (1, 1) = (a (2, 2) * a (3, 3) - a (2, 3) * a (3, 2) )     &
!OLD     / det                                                          
!OLD     imat (1, 2) = - (a (1, 2) * a (3, 3) - a (1, 3) * a (3, 2) )   &
!OLD     / det                                                          
!OLD     imat (1, 3) = (a (1, 2) * a (2, 3) - a (1, 3) * a (2, 2) )     &
!OLD     / det                                                          
!                                                                       
!OLD     imat (2, 1) = - (a (2, 1) * a (3, 3) - a (3, 1) * a (2, 3) )   &
!OLD     / det                                                          
!OLD     imat (2, 2) = (a (1, 1) * a (3, 3) - a (1, 3) * a (3, 1) )     &
!OLD     / det                                                          
!OLD     imat (2, 3) = - (a (1, 1) * a (2, 3) - a (1, 3) * a (2, 1) )   &
!OLD     / det                                                          
!                                                                       
!OLD     imat (3, 1) = (a (2, 1) * a (3, 2) - a (3, 1) * a (2, 2) )     &
!OLD     / det                                                          
!OLD     imat (3, 2) = - (a (1, 1) * a (3, 2) - a (1, 2) * a (3, 1) )   &
!OLD     / det                                                          
!OLD     imat (3, 3) = (a (1, 1) * a (2, 2) - a (1, 2) * a (2, 1) )     &
!OLD     / det                                                          
!                                                                       
!OLD  ELSE 
!OLD     ier_num = - 1 
!OLD     ier_typ = ER_MATH 
!OLD  ENDIF 
!                                                                       
!OLD  END SUBROUTINE invmat                         
!*****7*****************************************************************
!OLD  SUBROUTINE invmat4 (matrix) 
!-                                                                      
!     inverts a 4*4 Symmetry operation                                  
!+                                                                      
!OLD  USE errlist_mod 
!OLD  IMPLICIT none 
!                                                                       
!                                                                       
!OLD  REAL matrix (4, 4) 
!OLD  INTEGER i, j 
!OLD  REAL a (3, 3), b (3, 3), t (3) 
!                                                                       
!OLD  DO i = 1, 3 
!OLD  DO j = 1, 3 
!OLD  a (i, j) = matrix (i, j) 
!OLD  ENDDO 
!OLD  t (i) = matrix (i, 4) 
!OLD  ENDDO 
!                                                                       
!OLD  CALL invmat (b, a) 
!                                                                       
!OLD  IF (ier_num.eq.0) then 
!OLD     DO i = 1, 3 
!OLD     matrix (i, 4) = 0.0 
!OLD     DO j = 1, 3 
!OLD     matrix (i, j) = b (i, j) 
!OLD     matrix (i, 4) = matrix (i, 4) - b (i, j) * t (j) 
!OLD     ENDDO 
!OLD     ENDDO 
!OLD  ENDIF 
!                                                                       
!OLD  END SUBROUTINE invmat4                        
!*****7*****************************************************************
END MODULE tensors_mod
