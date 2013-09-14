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
      SUBROUTINE lattice (a0, ar, eps, gten, reps, rten, win, wrez, vol,&
      vr, lout)                                                         
!+                                                                      
!           Calculates lattice constants, metric and reciprocal metric  
!           tensor, permutation tensors and unit cell volume.           
!     It's done quite some old fashioned way, rather than calculating   
!     the direct metric tensor and its inverse.                         
!-                                                                      
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER i, i1, i2, j 
      LOGICAL lout 
      REAL a0 (3), ar (3), eps (3, 3, 3), gten (3, 3) 
      REAL reps (3, 3, 3), rten (3, 3) 
      REAL win (3), wrez (3), vol, vr 
      REAL cosa, cosb, cosg, cos1, cos2, sin1, sin2 
      REAL cosi, sind, cosd, acosd 
!                                                                       
      cosa = cosd (win (1) ) 
      cosb = cosd (win (2) ) 
      cosg = cosd (win (3) ) 
      vol = 1 - cosa * cosa - cosb * cosb - cosg * cosg 
      vol = vol + 2 * cosa * cosb * cosg 
!                                                                       
      IF (vol.gt.0.0) then 
         vol = sqrt (vol) * a0 (1) * a0 (2) * a0 (3) 
         vr = 1. / vol 
!                                                                       
!------ - calculate direct metric tensor                                
!                                                                       
         CALL tensor (gten, a0, win) 
!                                                                       
!------ - calculate reciprocal lattice constants                        
!                                                                       
         DO i = 1, 3 
         i1 = mod (i, 3) + 1 
         i2 = mod (i + 1, 3) + 1 
         ar (i) = a0 (i1) * a0 (i2) * sind (win (i) ) / vol 
         cos1 = cosd (win (i1) ) 
         cos2 = cosd (win (i2) ) 
         cosi = cosd (win (i) ) 
         sin1 = sind (win (i1) ) 
         sin2 = sind (win (i2) ) 
         wrez (i) = acosd ( (cos1 * cos2 - cosi) / (sin1 * sin2) ) 
         ENDDO 
!                                                                       
!------ - calculate reciprocal tensor                                   
!                                                                       
         CALL tensor (rten, ar, wrez) 
!                                                                       
!------ - calculate premutation tensors                                 
!                                                                       
         eps (1, 2, 3) = vol 
         eps (2, 3, 1) = vol 
         eps (3, 1, 2) = vol 
         eps (1, 3, 2) = - vol 
         eps (3, 2, 1) = - vol 
         eps (2, 1, 3) = - vol 
         reps (1, 2, 3) = vr 
         reps (2, 3, 1) = vr 
         reps (3, 1, 2) = vr 
         reps (1, 3, 2) = - vr 
         reps (3, 2, 1) = - vr 
         reps (2, 1, 3) = - vr 
!                                                                       
!------ - output ?                                                      
!                                                                       
         IF (lout) then 
            WRITE (output_io, 2001) (a0 (i), i = 1, 3), (win (i),       &
            i = 1, 3), vol                                              
            WRITE (output_io, 2002) ( (gten (i, j), j = 1, 3), i = 1, 3) 
            WRITE (output_io, 2003) (ar (i), i = 1, 3), (wrez (i),      &
            i = 1, 3), vr                                               
            WRITE (output_io, 2004) ( (rten (i, j), j = 1, 3), i = 1, 3) 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 35 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 2001 FORMAT     ( ' Lattice constants :'                               &
     &           ,/,4x,'a',10x,'b',10x,'c', 9x,                         &
     &           'alpha',6x,'beta',7x,'gamma', 6x,'volume',             &
     &           /,6(2X,F9.5),2X,G12.6)                                 
 2002 FORMAT     (/' Metric Tensor     :'/(3(' ',3(2X,F11.5)/))) 
 2003 FORMAT     ( ' Reciprocal Lattice constants :'                    &
     &           ,/,4x,'a*', 9x,'b*', 9x,'c*', 8x,                      &
     &           'alpha*',5x,'beta*',6x,'gamma*', 5x,'volume',          &
     &           /,6(2X,F9.5),2X,G12.6)                                 
 2004 FORMAT     (/' Reciprocal metric tensor     : '/                  &
     &            (3(' ',3(2X,F11.5)/)))                                
      END SUBROUTINE lattice                        
!*****7*****************************************************************
      SUBROUTINE tensor (ten, vec, win) 
!+                                                                      
!     Calculates the metric tensor. Works both for direct and           
!     reciprocal metric tensor.                                         
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
!                                                                       
      REAL ten (idim, idim), vec (idim), win (idim) 
      REAL cosd 
      INTEGER i, j 
!                                                                       
      DO i = 1, idim 
      DO j = 1, idim 
      IF (i.ne.j) then 
         ten (i, j) = vec (i) * vec (j) * cosd (win (6 - (i + j) ) ) 
      ELSE 
         ten (i, j) = vec (i) * vec (j) 
      ENDIF 
      ENDDO 
      ENDDO 
      END SUBROUTINE tensor                         
