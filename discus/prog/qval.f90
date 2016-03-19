MODULE qval_mod
!
INTEGER, PARAMETER :: val_inten  =  1
INTEGER, PARAMETER :: val_ampli  =  2
INTEGER, PARAMETER :: val_phase  =  3
INTEGER, PARAMETER :: val_real   =  4
INTEGER, PARAMETER :: val_imag   =  5
INTEGER, PARAMETER :: val_ranph  =  6
INTEGER, PARAMETER :: val_sq     =  7
INTEGER, PARAMETER :: val_fq     =  8
INTEGER, PARAMETER :: val_f2aver =  9
INTEGER, PARAMETER :: val_faver2 = 10
INTEGER, PARAMETER :: val_faver  = 11
INTEGER, PARAMETER :: val_norm   = 12
CONTAINS
!*****7*****************************************************************
      REAL FUNCTION qval (i, value, ix, iy, laver) 
!-                                                                      
!     transforms the real and imaginary part of the Fourier transform   
!     into the desired output format                                    
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE output_mod 
      USE random_mod
      USE wink_mod
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
      REAL      :: q2     = 0.0
      REAL      :: faver2 = 0.0
      REAL      :: f2aver = 0.0
!                                                                       
      REAL atan2d 
      REAL ran1 
      LOGICAL laver 
      qval   = 0.0
      f2aver = 0.0
      faver2 = 0.0
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
      IF (value == val_inten) then 
         IF (laver) then 
            qval = real (f * conjg (f) ) 
         ELSE 
            qval = dsi (i) 
         ENDIF 
!                                                                       
!     Calculate amplitude 'amplitude'                                   
!                                                                       
      ELSEIF (value == val_ampli) then 
         qval = sqrt (real (f * conjg (f) ) ) 
!                                                                       
!     Calculate phase 'phase'                                           
!                                                                       
      ELSEIF (value == val_phase) then 
         IF (f.eq. (0, 0) ) then 
            qval = 0.0 
         ELSE 
            qval = atan2d (aimag (f), real (f) ) 
         ENDIF 
!                                                                       
!     Calculate real part 'real'                                        
!                                                                       
      ELSEIF (value == val_real) then 
         qval = real (f) 
!                                                                       
!     Calculate imaginary part 'imaginary'                              
!                                                                       
      ELSEIF (value == val_imag) then 
         qval = aimag (f) 
!                                                                       
!     Calculate phase 'phase', random, except for integer hkl           
!                                                                       
      ELSEIF (value == val_ranph) then 
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
!
!     Calculate S(Q) = I/<f>^2/N + 1-e^(q^2*u^2) normalized intensity plus inelastic part 
!
      ELSEIF (value == val_sq) then
         IF (laver) then 
            qval = real (f * conjg (f) ) 
         ELSE 
            qval = dsi (i) 
         ENDIF 
         DO k=1,cr_nscat
            faver2 = faver2 + (REAL(cfact_pure(istl(i),k)))*cr_amount(k)
         ENDDO
         faver2 = faver2**2
         q2   = zpi**2*(2*istl(i)*CFINC)**2
         qval = qval /faver2/ cr_n_real_atoms &
                +1.0 - exp(-q2*cr_u2aver)
!
!     Calculate E(Q) = I/<f^2>/N normalized intensity 
!
      ELSEIF (value == val_norm) then
         IF (laver) then 
            qval = real (f * conjg (f) ) 
         ELSE 
            qval = dsi (i) 
         ENDIF 
         DO k=1,cr_nscat
            f2aver = f2aver + (REAL(cfact_pure(istl(i),k)))**2*cr_amount(k)
         ENDDO
         q2   = zpi**2*(2*istl(i)*CFINC)**2
         qval = qval /f2aver/ cr_n_real_atoms 
!
!     Calculate average squared atomic form factor <f**2>
!
      ELSEIF (value == val_f2aver) then
         DO k=1,cr_nscat
            qval = qval + (REAL(cfact_pure(istl(i),k)))**2*cr_amount(k)
         ENDDO
         qval = qval / cr_n_real_atoms
!
!     Calculate average atomic form factor squared <f>**2
!
      ELSEIF (value == val_faver2) then
         DO k=1,cr_nscat
            qval = qval + REAL(cfact_pure(istl(i),k))*cr_amount(k)
         ENDDO
         qval = qval / cr_n_real_atoms
         qval = qval**2
!
!     Calculate average atomic form factor <f>
!
      ELSEIF (value == val_faver) then
         DO k=1,cr_nscat
            qval = qval + REAL(cfact_pure(istl(i),k))*cr_amount(k)
         ENDDO
         qval = qval / cr_n_real_atoms
      ENDIF 
!                                                                       
      END FUNCTION qval                             
END MODULE qval_mod
