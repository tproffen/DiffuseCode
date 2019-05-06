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
INTEGER, PARAMETER :: val_iq     = 13
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
      USE trig_degree_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
      REAL DELTA 
      PARAMETER (DELTA = 0.000001) 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN) :: i
      INTEGER, INTENT(IN) :: value
      INTEGER, INTENT(IN) :: ix
      INTEGER, INTENT(IN) :: iy 
      LOGICAL, INTENT(IN) :: laver
! 
      INTEGER k 
!                                                                       
      COMPLEX(KIND=KIND(0.0D0)) :: f 
      REAL h (3) 
      REAL      :: q2     = 0.0
      REAL      :: faver2 = 0.0
      REAL      :: f2aver = 0.0
      REAL      :: signum = 1.0
!                                                                       
!     REAL atan2d 
      REAL ran1 
      qval   = 0.0
      f2aver = 0.0
      faver2 = 0.0
!                                                                       
!------ Get values of F or <F>                                          
!                                                                       
      IF (laver) THEN 
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
      IF (value == val_inten) THEN 
         IF (laver) THEN 
            qval = REAL (f * CONJG (f) , KIND=KIND(1.0E0)) 
         ELSE 
            qval = REAL(dsi (i),KIND=KIND(0.0E0)) 
         ENDIF 
!                                                                       
!     Calculate amplitude 'amplitude'                                   
!                                                                       
      ELSEIF (value == val_ampli) THEN 
         qval = SQRT (REAL (f * CONJG (f) , KIND=KIND(1.0E0) )) 
!                                                                       
!     Calculate phase 'phase'                                           
!                                                                       
      ELSEIF (value == val_phase) THEN 
         IF (f.eq. (0, 0) ) THEN 
            qval = 0.0 
         ELSE 
            qval = REAL(atan2d (AIMAG (f), REAL (f) ) )
         ENDIF 
!                                                                       
!     Calculate real part 'real'                                        
!                                                                       
      ELSEIF (value == val_real) THEN 
         qval = REAL (f, KIND=KIND(1.0E0)) 
!                                                                       
!     Calculate imaginary part 'imaginary'                              
!                                                                       
      ELSEIF (value == val_imag) THEN 
         qval = REAL(AIMAG (f))
!                                                                       
!     Calculate phase 'phase', random, except for integer hkl           
!                                                                       
      ELSEIF (value == val_ranph) THEN 
         DO k = 1, 3 
            h (k) = out_eck (k, 1) + out_vi (k, 1) * REAL(ix - 1)     &
                                   + out_vi (k, 2) * REAL(iy - 1)
         ENDDO 
         IF (ABS (h (1) - NINT (h (1) ) ) .lt.DELTA.AND. &
             ABS (h (2) - NINT (h (2) ) ) .lt.DELTA.AND. &
             ABS (h (3) - NINT (h (3) ) ) .lt.DELTA) THEN                                                
            IF (f.eq. (0, 0) ) THEN 
               qval = 0.0 
            ELSE 
               qval = REAL(atan2d (AIMAG (f), REAL (f) ) )
            ENDIF 
         ELSE 
            qval = (ran1 (idum) - 0.5) * 360. 
         ENDIF 
!
!     Calculate S(Q) = I/<f>^2/N + 1-e^(q^2*u^2) normalized intensity plus inelastic part 
!
      ELSEIF (value == val_sq) THEN
         IF (laver) THEN 
            qval = REAL (f * CONJG (f), KIND=KIND(0.0E0) ) 
         ELSE 
            qval = REAL(dsi (i) , KIND=KIND(1.0E0))
         ENDIF 
         DO k=1,cr_nscat
            signum = 1.0
            IF(REAL(cfact_pure(1,k))<0.0) signum = -1.0
!           faver2 = faver2 + (REAL(cfact_pure(istl(i),k), KIND=KIND(0.0E0)))*cr_amount(k)
            faver2 = faver2 +                                    &
                     SQRT(DBLE (       cfact_pure(istl(i), k)  * &
                                conjg (cfact_pure(istl(i), k)))) &
                     *cr_amount(k)*signum
         ENDDO
         faver2 = faver2**2
         q2   = REAL(zpi**2*(REAL(2*istl(i),KIND=KIND(0.0D0))*CFINC)**2)
         qval = qval /faver2/ cr_n_real_atoms &
                +1.0 - EXP(-q2*cr_u2aver)
!
!     Calculate E(Q) = I/<f^2>/N normalized intensity 
!
      ELSEIF (value == val_norm) THEN
         IF (laver) THEN 
            qval = REAL (f * CONJG (f),KIND=KIND(0.0E0) ) 
         ELSE 
            qval = REAL(dsi (i), KIND=KIND(0.0E0))
         ENDIF 
         DO k=1,cr_nscat
!           f2aver = f2aver + REAL(cfact_pure(istl(i),k)**2,KIND=KIND(0.0E0))*cr_amount(k)
            f2aver = f2aver +                                 &
                      DBLE (       cfact_pure(istl(i), k)  *  &
                             conjg (cfact_pure(istl(i), k)))  &
                   *cr_amount(k)
         ENDDO
         q2   = REAL(zpi**2*(REAL(2*istl(i),KIND=KIND(0.0E0))*CFINC)**2)
         qval = qval /f2aver/ cr_n_real_atoms 
!
!     Calculate I(Q) = I/N normalized intensity 
!
      ELSEIF (value == val_iq) THEN
         IF (laver) THEN 
            qval = REAL (f * CONJG (f),KIND=KIND(0.0E0) ) 
         ELSE 
            qval = REAL(dsi (i), KIND=KIND(0.0E0))
         ENDIF 
         qval = qval /f2aver/ cr_n_real_atoms 
!
!     Calculate average squared atomic form factor <f**2>
!
      ELSEIF (value == val_f2aver) THEN
         DO k=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),k)**2,KIND=KIND(0.0E0))*cr_amount(k)
            qval = qval +                                     &
                      DBLE (       cfact_pure(istl(i), k)  *  &
                             conjg (cfact_pure(istl(i), k)))  &
                   *cr_amount(k)
         ENDDO
         qval = qval / cr_n_real_atoms
!
!     Calculate average atomic form factor squared <f>**2
!
      ELSEIF (value == val_faver2) THEN
         DO k=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),k),KIND=KIND(0.0E0))*cr_amount(k)
            signum = 1.0
            IF(REAL(cfact_pure(1,k))<0.0) signum = -1.0
            qval = qval +                                     &
                  SQRT(DBLE (       cfact_pure(istl(i), k)  * &
                             conjg (cfact_pure(istl(i), k)))) &
                  * cr_amount(k) * signum
         ENDDO
         qval = qval / cr_n_real_atoms
         qval = qval**2
!
!     Calculate average atomic form factor <f>
!
      ELSEIF (value == val_faver) THEN
         DO k=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),k),KIND=KIND(0.0E0))*cr_amount(k)
            qval = qval +                                      &
                   SQRT(DBLE (       cfact_pure(istl(i), k)  * &
                              conjg (cfact_pure(istl(i), k)))) &
                   * cr_amount(k) * signum
         ENDDO
         qval = qval / cr_n_real_atoms
      ENDIF 
!                                                                       
      END FUNCTION qval                             
END MODULE qval_mod
