MODULE qval_mod
!
CONTAINS
!*****7*****************************************************************
      REAL FUNCTION qval (i, value, ix, iy, laver) 
!-                                                                      
!     transforms the real and imaginary part of the Fourier transform   
!     into the desired output format                                    
!+                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      USE random_mod
      IMPLICIT none 
!                                                                       
      REAL DELTA 
      PARAMETER (DELTA = 0.000001) 
!                                                                       
!                                                                       
      INTEGER i, value, ix, iy 
      INTEGER k 
!                                                                       
      COMPLEX f 
      REAL h (3) 
!                                                                       
      REAL atan2d 
      REAL ran1 
      LOGICAL laver 
      qval = 0.0
!                                                                       
!------ Get values of F or <F>                                          
!                                                                       
      IF (laver) then 
         f = acsf (i) 
      ELSE 
         f = csf (i) 
      ENDIF 
!                                                                       
!     Calculate intensity 'intensity'                                   
!                                                                       
!------ We have to store dsi() here, because if lots are                
!------ used, csf() will only contain the values for the                
!------ last lot !!                                                     
!                                                                       
      IF (value.eq.1) then 
         IF (laver) then 
            qval = real (f * conjg (f) ) 
         ELSE 
            qval = dsi (i) 
         ENDIF 
!                                                                       
!     Calculate amplitude 'amplitude'                                   
!                                                                       
      ELSEIF (value.eq.2) then 
         qval = sqrt (real (f * conjg (f) ) ) 
!                                                                       
!     Calculate phase 'phase'                                           
!                                                                       
      ELSEIF (value.eq.3) then 
         IF (f.eq. (0, 0) ) then 
            qval = 0.0 
         ELSE 
            qval = atan2d (aimag (f), real (f) ) 
         ENDIF 
!                                                                       
!     Calculate real part 'real'                                        
!                                                                       
      ELSEIF (value.eq.4) then 
         qval = real (f) 
!                                                                       
!     Calculate imaginary part 'imaginary'                              
!                                                                       
      ELSEIF (value.eq.5) then 
         qval = aimag (f) 
!                                                                       
!     Calculate phase 'phase', random, except for integer hkl           
!                                                                       
      ELSEIF (value.eq.6) then 
         DO k = 1, 3 
         h (k) = out_eck (k, 1) + out_vi (k, 1) * float (ix - 1)        &
         + out_vi (k, 2) * float (iy - 1)                               
         ENDDO 
         IF (abs (h (1) - nint (h (1) ) ) .lt.DELTA.and.abs (h (2)      &
         - nint (h (2) ) ) .lt.DELTA.and.abs (h (3) - nint (h (3) ) )   &
         .lt.DELTA) then                                                
            IF (f.eq. (0, 0) ) then 
               qval = 0.0 
            ELSE 
               qval = atan2d (aimag (f), real (f) ) 
            ENDIF 
         ELSE 
            qval = (ran1 (idum) - 0.5) * 360. 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION qval                             
END MODULE qval_mod
