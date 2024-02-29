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
INTEGER, PARAMETER :: val_f2averb= 17
INTEGER, PARAMETER :: val_faver2b= 18
INTEGER, PARAMETER :: val_faverb = 19
!
INTEGER, PARAMETER :: MAXVALS = 19
!
CHARACTER(LEN=17) :: cvalue (0:MAXVALS)
!
DATA cvalue / 'undefined        ', 'Intensity        ', 'Amplitude        ',&
              'Phase angle      ', 'Real Part        ', 'Imaginary Part   ',&
              'Random Phase     ', 'S(Q)             ', 'F(Q)             ',&
              'f2aver = <f^2>   ', 'faver2 = <f>^2   ', 'faver = <f>      ',&
              'Normal Inten     ', 'I(Q)             ', 'PDF              ',&
              '3DPDF            ', '3DBETA           ', 'f2averb=<f^2>*DBW',&
              'faver2b=<f>^2*DBW', 'faverb = <F>     '                      &
            /
CONTAINS
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION qval (i, j, k, value, ix, iy, laver)      ! We probably need to extend this function to 3D, so it recieves a vector of indices
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
INTEGER, INTENT(IN) :: i          ! Index along a*
INTEGER, INTENT(IN) :: j          ! Index along b*
INTEGER, INTENT(IN) :: k          ! Index along c*
INTEGER, INTENT(IN) :: value
INTEGER, INTENT(IN) :: ix
INTEGER, INTENT(IN) :: iy 
LOGICAL, INTENT(IN) :: laver
! 
INTEGER l 
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
         f = acsf (i,j,k)                                                                       ! Who is i here ? I think is an argument of the function, but how to extend to 3D ?
      ELSE 
         f = csf (i,j,k)                                                                        ! Who is i here ? I think is an argument of the function, but how to extend to 3D ?
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
            qval = REAL(dsi (i,j,k),KIND=KIND(0.0D0)) 
         ENDIF 
!
!     Calculate 3DPDF
!
      ELSEIF (value == val_3DPDF .or. value == val_3DBETA) THEN
         qval = REAL(rpdf(i,j,k), KIND=KIND(0.0D0))
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
         DO l = 1, 3 
            h (l) = out_eck (l, 1) + out_vi (l, 1) * REAL(ix - 1)     &
                                   + out_vi (l, 2) * REAL(iy - 1)
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
            qval = dsi(i,j,k)
         ENDIF 
         DO l=1,cr_nscat
            signum = 1.0
            IF(REAL(cfact_pure(1,l), kind=PREC_DP)<0.0) signum = -1.0D0
!           faver2 = faver2 + (REAL(cfact_pure(istl(i),l), KIND=KIND(0.0E0)))*cr_amount(l)
            faver2 = faver2 +                                    &
                     SQRT(DBLE (       cfact_pure(istl(i,j,1), l)  * &
                                conjg (cfact_pure(istl(i,j,1), l)))) &
                     *cr_amount(l)*signum
         ENDDO
         faver2 = faver2**2
         q2   = REAL(zpi**2*(REAL(2*istl(i,1,1),KIND=KIND(0.0D0))*CFINC)**2)
         qval = qval /faver2/ cr_n_real_atoms &
                +1.0 - EXP(-q2*cr_u2aver)
!
!     Calculate I(Q) = I/N normalized intensity 
!
      ELSEIF (value == val_iq) THEN
         IF (laver) THEN 
            qval = f * CONJG (f)
         ELSE 
            qval = dsi (i,j,k)
         ENDIF 
         qval = qval        / cr_n_real_atoms 
!
!     Calculate average squared atomic form factor <f**2>
!
      ELSEIF (value == val_f2aver) THEN
         DO l=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),l)**2,KIND=KIND(0.0E0))*cr_amount(l)
            qval = qval +                                     &
                      DBLE (       cfact_pure(istl(i,j,1), l)  *  &
                             conjg (cfact_pure(istl(i,j,1), l)))  &
                   *cr_amount(l)
         ENDDO
         qval = qval / cr_n_real_atoms
!
!     Calculate average atomic form factor squared <f>**2
!
      ELSEIF (value == val_faver2) THEN
         DO l=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),l),KIND=KIND(0.0E0))*cr_amount(l)
            signum = 1.0
            IF(REAL(cfact_pure(1,l))<0.0) signum = -1.0
            qval = qval +                                     &
                  SQRT(DBLE (       cfact_pure(istl(i,j,1), l)  * &
                             conjg (cfact_pure(istl(i,j,1), l)))) &
                  * cr_amount(l) * signum
         ENDDO
         qval = qval / cr_n_real_atoms
         qval = qval**2
!
!     Calculate average atomic form factor <f>
!
      ELSEIF (value == val_faver) THEN
         DO l=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),l),KIND=KIND(0.0E0))*cr_amount(l)
            qval = qval +                                      &
                   SQRT(DBLE (       cfact_pure(istl(i,j,1), l)  * &
                              conjg (cfact_pure(istl(i,j,1), l)))) &
                   * cr_amount(l) * signum
         ENDDO
         qval = qval / cr_n_real_atoms
!
!     Calculate average squared atomic form factor <f**2>
!
ELSEIF (value == val_f2averb) THEN
   DO l=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),l)**2,KIND=KIND(0.0E0))*cr_amount(l)
      qval = qval +                                     &
                DBLE (       cfact     (istl(i,j,1), l)  *  &
                       conjg (cfact     (istl(i,j,1), l)))  &
             *cr_amount(l)
   ENDDO
   qval = qval / cr_n_real_atoms
!
!     Calculate average atomic form factor squared <f>**2
!
ELSEIF (value == val_faver2b) THEN
   DO l=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),l),KIND=KIND(0.0E0))*cr_amount(l)
      signum = 1.0
      IF(REAL(cfact_pure(1,l))<0.0) signum = -1.0
      qval = qval +                                     &
            SQRT(DBLE (       cfact     (istl(i,j,1), l)  * &
                       conjg (cfact     (istl(i,j,1), l)))) &
            * cr_amount(l) * signum
   ENDDO
   qval = qval / cr_n_real_atoms
   qval = qval**2
!
!     Calculate average atomic form factor <f>
!
ELSEIF (value == val_faverb) THEN
   DO l=1,cr_nscat
!           qval = qval + REAL(cfact_pure(istl(i),l),KIND=KIND(0.0E0))*cr_amount(l)
      qval = qval +                                      &
             SQRT(DBLE (       cfact     (istl(i,j,1), l)  * &
                        conjg (cfact     (istl(i,j,1), l)))) &
             * cr_amount(l) * signum
   ENDDO
   qval = qval / cr_n_real_atoms
ENDIF 
!                                                                       
END FUNCTION qval                             
!
!*******************************************************************************
!
END MODULE qval_mod
