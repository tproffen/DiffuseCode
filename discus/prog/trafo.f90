!*****7*****************************************************************
!                                                                       
      SUBROUTINE trafo (hkl, u, xc, yc, zc, gmat, fmat, dist, eps, gten,&
      reps, rten)                                                       
!+                                                                      
!     Calculates the transformation matrix : crystal <> plotsection     
!     from lattice comstant and hkl of viewing vector. Also             
!     calculates the distance of the plotting section from the origin   
!     of the crystal.                                                   
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER inull2 (2), i, j, k, l, m, inull 
      REAL eps (3, 3, 3), fmat (3, 3), gmat (3, 3), gten (3, 3) 
      REAL rten (3, 3), u (3), xc (3), yc (3), zc (3), dist 
      REAL reps (3, 3, 3), hkl (3), zcc, xcc, quad 
!                                                                       
!     Determine achses for transformation                               
!     Transform reciprocal vector HKL to real space vector ZC           
!                                                                       
      DO i = 1, 3 
      zc (i) = 0.0 
      DO j = 1, 3 
      zc (i) = zc (i) + rten (i, j) * hkl (j) 
      ENDDO 
      ENDDO 
!                                                                       
      zcc = sqrt (quad (zc, zc, gten) ) 
      IF (zcc.eq.0.0) return 
      inull = 0 
      inull2 (1) = 0 
      inull2 (2) = 0 
      DO i = 1, 3 
      zc (i) = zc (i) / zcc 
      IF (hkl (i) .eq.0.0) then 
         inull = inull + 1 
         inull2 (inull) = i 
      ENDIF 
      xc (i) = 0.0 
      yc (i) = 0.0 
      ENDDO 
!                                                                       
!     Determine a direction as fundamental as possible for XC           
!                                                                       
      IF (inull.eq.2.or.inull.eq.1) then 
         xc (inull2 (1) ) = 1.0 
      ELSE 
         IF (hkl (1) .eq.hkl (2) ) then 
            xc (1) = hkl (1) 
            xc (2) = - hkl (2) 
         ELSEIF (hkl (1) .eq.hkl (3) ) then 
            xc (1) = hkl (1) 
            xc (3) = - hkl (3) 
         ELSEIF (hkl (2) .eq.hkl (3) ) then 
            xc (2) = hkl (2) 
            xc (3) = - hkl (3) 
         ELSE 
            xc (1) = hkl (2) 
            xc (2) = - hkl (1) 
         ENDIF 
      ENDIF 
      xcc = sqrt (quad (xc, xc, gten) ) 
!                                                                       
!     normalize XC                                                      
!                                                                       
      DO i = 1, 3 
      xc (i) = xc (i) / xcc 
      ENDDO 
!                                                                       
!------ calculate yc as the vector product of zc and xc                 
!                                                                       
      DO m = 1, 3 
      yc (m) = 0.0 
      DO j = 1, 3 
      DO k = 1, 3 
      DO l = 1, 3 
      yc (m) = yc (m) + eps (j, k, l) * zc (k) * xc (l) * rten (j, m) 
      ENDDO 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     calculate indices of gmat                                         
!                                                                       
      DO i = 1, 3 
      gmat (i, 1) = xc (i) 
      gmat (i, 2) = yc (i) 
      gmat (i, 3) = zc (i) 
      ENDDO 
!                                                                       
!     calculate distance to origin                                      
!                                                                       
      CALL invmat (fmat, gmat) 
      dist = fmat (3, 1) * u (1) + fmat (3, 2) * u (2) + fmat (3, 3)    &
      * u (3)                                                           
!                                                                       
      END SUBROUTINE trafo                          
!*****7*****************************************************************
      SUBROUTINE trans (uc, gmat, up, idim) 
!+                                                                      
!     Transforms a point in the crystal space into plot space           
!     and vice versa                                                    
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER i, j, idim 
      REAL gmat (idim, idim), uc (idim), up (idim) 
!                                                                       
      DO i = 1, idim 
      up (i) = 0.0 
      DO j = 1, idim 
      up (i) = up (i) + gmat (i, j) * uc (j) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE trans                          
