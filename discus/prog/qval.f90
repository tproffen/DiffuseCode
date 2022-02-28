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
!INTEGER, PARAMETER :: val_norm   = 12
INTEGER, PARAMETER :: val_iq     = 13
INTEGER, PARAMETER :: val_pdf    = 14
INTEGER, PARAMETER :: val_3DPDF  = 15
INTEGER, PARAMETER :: val_3DBETA = 16
!
INTEGER, PARAMETER :: MAXVALS = 16
!
CHARACTER(LEN=14) :: cvalue (0:MAXVALS)
!
DATA cvalue / 'undefined     ', 'Intensity     ', 'Amplitude     ',&
              'Phase angle   ', 'Real Part     ', 'Imaginary Part',&
              'Random Phase  ', 'S(Q)          ', 'F(Q)          ',&
              'f2aver = <f^2>', 'faver2 = <f>^2', 'faver = <f>   ',&
              'Normal Inten  ', 'I(Q)          ', 'PDF           ',&
              '3DPDF         ', '3DBETA        '                   &
            /
CONTAINS
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION qval (i, value, ix, iy, laver) 
!-                                                                      
!     transforms the real and imaginary part of the Fourier transform   
!     into the desired output format                                    
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod 
USE output_mod 
USE lib_random_func
use precision_mod
USE random_mod
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), PARAMETER :: DELTA = 0.000001d0
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
REAL(kind=PREC_DP), dimension(3) :: h (3) 
REAL(kind=PREC_DP)      :: q2     = 0.0
REAL(kind=PREC_DP)      :: faver2 = 0.0
REAL(kind=PREC_DP)      :: f2aver = 0.0
REAL(kind=PREC_DP)      :: signum = 1.0
!                                                                       
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
            qval = REAL(f * CONJG (f) , KIND=KIND(1.0D0)) 
         ELSE 
            qval = REAL(dsi (i),KIND=KIND(0.0D0)) 
         ENDIF 
!
!     Calculate 3DPDF
!
      ELSEIF (value == val_3DPDF .or. value == val_3DBETA) THEN
         qval = REAL(rpdf(i), KIND=KIND(0.0D0))
!                                                                       
!     Calculate amplitude 'amplitude'                                   
!                                                                       
      ELSEIF (value == val_ampli) THEN 
         qval = SQRT (REAL(f * CONJG (f) , KIND=KIND(1.0D0) )) 
!                                                                       
!     Calculate phase 'phase'                                           
!                                                                       
      ELSEIF (value == val_phase) THEN 
         IF (f.eq. (0, 0) ) THEN 
            qval = 0.0D0 
         ELSE 
            qval = REAL(atan2d(DIMAG(f), REAL(f, kind=PREC_DP) ), kind=PREC_DP )
         ENDIF 
!                                                                       
!     Calculate real part 'real'                                        
!                                                                       
      ELSEIF (value == val_real) THEN 
         qval = REAL(f, KIND=KIND(1.0D0)) 
!                                                                       
!     Calculate imaginary part 'imaginary'                              
!                                                                       
      ELSEIF (value == val_imag) THEN 
         qval = REAL(DIMAG (f))
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
               qval = 0.0D0
            ELSE 
               qval = atan2d (DIMAG(f), REAL(f, kind=PREC_DP) )
            ENDIF 
         ELSE 
            qval = (ran1 (idum) - 0.5D0) * 360.D0 
         ENDIF 
!
!     Calculate S(Q) = I/<f>^2/N + 1-e^(q^2*u^2) normalized intensity plus inelastic part 
!
      ELSEIF (value == val_sq) THEN
         IF (laver) THEN 
            qval = f * CONJG(f)
         ELSE 
            qval = dsi(i)
         ENDIF 
         DO k=1,cr_nscat
            signum = 1.0
            IF(REAL(cfact_pure(1,k), kind=PREC_DP)<0.0) signum = -1.0D0
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
!     Calculate I(Q) = I/N normalized intensity 
!
      ELSEIF (value == val_iq) THEN
         IF (laver) THEN 
            qval = f * CONJG (f)
         ELSE 
            qval = dsi (i)
         ENDIF 
         qval = qval        / cr_n_real_atoms 
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
