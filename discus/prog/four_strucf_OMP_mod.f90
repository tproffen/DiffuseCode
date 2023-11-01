MODULE four_strucf_mod
!
!
!  At least for a transition period, several versions of the Fourier 
!  algorithm are retained. The current serial and parallel algorithms
!  four_strucf_OMP and four-strucf_serial are both much faster than
!  the original algorithm is compiled with -O3 -ffast-math
!
!  four_strucf_omp     ! Best (?) OMP algorithm
!  four_strucf_serial  ! Best (?) serial algorithm
!  four_strucf_omp_a   ! Good     OMP algorithm, complete loop reversal
!  four_strucf_omp_b   ! Bad      OMP algorithm, Almost like original serial
!  four_strucf_serial_original !  Original algorithm, not as efficient as believed...
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE four_strucf (iscat, lform) 
!
!+
!  Interface to four_strucf_serial((iscat, lform)
!           and four_strucf_omp((iscat, lform)
!-
!
!$ USE OMP_LIB
USE parallel_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!
INTEGER :: tid
INTEGER :: nthreads
!
nthreads = 1
!$OMP PARALLEL PRIVATE(tid)
!$   tid = OMP_GET_THREAD_NUM()
!$   IF (tid == 0) THEN
!$      IF(par_omp_maxthreads == -1) THEN
!$         nthreads = MAX(1,MIN(par_omp_phys, OMP_GET_NUM_THREADS()))
!$      ELSE
!$         nthreads = MAX(1,MIN(par_omp_maxthreads, par_omp_phys, OMP_GET_NUM_THREADS()))
!$      ENDIF
!$   END IF
!$OMP END PARALLEL
IF(par_omp_use .AND. nthreads>1) THEN
   CALL four_strucf_omp  (iscat, lform)
ELSE
   CALL four_strucf_serial(iscat, lform)
ENDIF
END SUBROUTINE four_strucf
!
!*****7*****************************************************************
!
SUBROUTINE four_strucf_omp (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!
!  Apparently a good OMP version.
!  The outer loop is j=0,num(1)-1 and the threadprivate variable 
!  tcsfp is reduced to size (0:num(2)*num(3)-1)
!  To reduce the number of calculations within the Fourier loop, the
!  increments xinc*, oinc*, iinc*, jinc* are calculated once up front.
!
!  This algorithm is roughly twice as fast on an optimized serial compilation
!  compared to the original serial Fourier calculation. Thus its serial form
!  replaces the old serial algorithm as well.
!-                                                                      
!$ USE OMP_LIB
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, jj
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
INTEGER (KIND=PREC_INT_LARGE)   :: num23,num123, num123_1
!
INTEGER :: IAND, ISHFT 
! COMPLEX(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP (Neder's original code)
COMPLEX(KIND=PREC_DP), DIMENSION(  :,:), ALLOCATABLE, SAVE :: tcsfp   ! My declaration
!$OMP THREADPRIVATE(tcsfp)
!
!------ zero fourier array                                              
!
! Neder's original code
!                                                                       
!tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
tcsf = cmplx(0.0D0, 0.0D0, kind=kind(0.0D0))
!
num23    =        num(2)*num(3)
num123   = num(1)*num(2)*num(3)
num123_1 = num(1)*num(2)*num(3)-1
!
!****************************************************************************************************************
! Neder's original code
!****************************************************************************************************************
! !                                                                       
! !------ Loop over all atoms in 'xat'                                    
! !                                                                       
! !$OMP PARALLEL
! !$OMP DO SCHEDULE(STATIC)
! DO k = 1, nxat 
!    xincu(k) = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
!    xincv(k) = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
!    xincw(k) = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
!    oincu(k) = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
!    oincv(k) = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
!    oincw(k) = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
! !                                                                       
!    iincu(k) = NINT (64 * I2PI * (xincu(k) - INT (xincu(k)) + 0.0d0) ) 
!    iincv(k) = NINT (64 * I2PI * (xincv(k) - INT (xincv(k)) + 0.0d0) ) 
!    iincw(k) = NINT (64 * I2PI * (xincw(k) - INT (xincw(k)) + 0.0d0) ) 
!    jincu(k) = NINT (64 * I2PI * (oincu(k) - INT (oincu(k)) + 0.0d0) ) 
!    jincv(k) = NINT (64 * I2PI * (oincv(k) - INT (oincv(k)) + 0.0d0) ) 
!    jincw(k) = NINT (64 * I2PI * (oincw(k) - INT (oincw(k)) + 0.0d0) ) 
! ENDDO
! !$OMP END DO NOWAIT
! !$OMP END PARALLEL
! !                                                                       
! !------ Loop over all atoms in 'xat'                                    
! !                                                                       
! !$OMP PARALLEL PRIVATE(iarg, iarg0, ii, jj)
!    IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
!    ALLOCATE (tcsfp (0:       num(2)*num(3)-1)) !,0:nthreads-1))
!    tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
!    !$OMP DO SCHEDULE(STATIC)
!       DO j = 0, num (1) - 1
!          jj = j*num(2)*num(3)
!          tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
!          loop_k: DO k = 1, nxat 
!             iarg0 =  lmn(1)*iincu(k) + lmn(2)*iincv(k) + lmn(3)*iincw(k) + &
!                      lmn(4)*jincu(k) + lmn(5)*jincv(k) + lmn(6)*jincw(k) + &
!                      j*iincu(k)
!             ii = 0 
!             DO i = 0, num (2) - 1
!                iarg = iarg0 +              iincv(k)*i 
!                DO h = 1, num (3) 
!                   tcsfp(ii)= tcsfp(ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )
!                   iarg     = iarg + iincw(k)
!                   ii       = ii + 1 
!                ENDDO 
!             ENDDO 
!          ENDDO  loop_k
!          DO ii = 1,        num(2)*num(3)
!             tcsf(jj+ii) = tcsf(jj+ii) + tcsfp(ii-1)
!          ENDDO
!       ENDDO 
! !!       !$OMP END CRITICAL
!    !$OMP END DO
!    DEALLOCATE(tcsfp)
! !$OMP END PARALLEL
!****************************************************************************************************************
! My code ( Copilot generated )
!****************************************************************************************************************
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC)
DO k = 1, nxat 
   xincu(k) = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
   xincv(k) = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
   xincw(k) = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
   oincu(k) = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
   oincv(k) = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
   oincw(k) = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
   iincu(k) = NINT (64 * I2PI * (xincu(k) - INT (xincu(k)) + 0.0d0) ) 
   iincv(k) = NINT (64 * I2PI * (xincv(k) - INT (xincv(k)) + 0.0d0) ) 
   iincw(k) = NINT (64 * I2PI * (xincw(k) - INT (xincw(k)) + 0.0d0) ) 
   jincu(k) = NINT (64 * I2PI * (oincu(k) - INT (oincu(k)) + 0.0d0) ) 
   jincv(k) = NINT (64 * I2PI * (oincv(k) - INT (oincv(k)) + 0.0d0) ) 
   jincw(k) = NINT (64 * I2PI * (oincw(k) - INT (oincw(k)) + 0.0d0) ) 
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
!$OMP PARALLEL PRIVATE(iarg, iarg0, ii, jj)
   IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
   ALLOCATE (tcsfp (          1:num(2), 1:num(3))) !,0:nthreads-1))
   tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
   !$OMP DO SCHEDULE(STATIC)
      DO j = 1, num (1)! - 1
!        jj = j*num(2)*num(3)
         tcsfp         = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
         loop_k: DO k = 1, nxat 
            iarg0 =  lmn(1)*iincu(k) + lmn(2)*iincv(k) + lmn(3)*iincw(k) + &
                     lmn(4)*jincu(k) + lmn(5)*jincv(k) + lmn(6)*jincw(k) + &
                     (j-1)*iincu(k)
!           ii = 0 
            DO i = 1, num (2)! - 1
               iarg = iarg0 +              iincv(k)*(i-1) 
               DO h = 1, num (3) 
!                 tcsfp(iincu(k), iincv(k), iincw(k)) = tcsfp(iincu(k), iincv(k), iincw(k)) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                  tcsfp(  i,h) = tcsfp(  i,h) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                  iarg     = iarg + iincw(k)
!                 ii       = ii + 1 
               ENDDO 
            ENDDO 
         ENDDO  loop_k
         DO i  = 1,        num(2)
         DO h  = 1,        num(3)
!           tcsf(iincu(k), iincv(k), iincw(k)) = tcsf(iincu(k), iincv(k), iincw(k)) + tcsfp(iincu(k), iincv(k), iincw(k))
            tcsf(j,i,h) = tcsf(j,i,h) + tcsfp(  i,h)
         ENDDO
         ENDDO
      ENDDO 
!!       !$OMP END CRITICAL
   !$OMP END DO
   DEALLOCATE(tcsfp)
!$OMP END PARALLEL
!****************************************************************************************************************
!write(*,*) 'NUM ', num
!write(*,*) ' CEX  ', minval(real(cex )), maxval(real(cex ))
!write(*,*) ' TCSF ', minval(real(tcsf)), maxval(real(tcsf))
!****************************************************************************************************************
!
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
   !DO  i = 1, num (1) * num (2) * num(3)
!  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
      !tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                   ! Neder's original code
!  k = 0
   DO i = 1, num (1)
      DO j = 1, num (2)
         DO h = 1, num (3)
!           k = k + 1
            tcsf (i,j,h) = tcsf (i,j,h) * cfact(istl(i,j,h), iscat)          ! My declaration
         ENDDO
      ENDDO
   ENDDO   
!  END FORALL
   !END DO
ENDIF 
!                                                                       
END SUBROUTINE four_strucf_omp
!
! !**********************************************************************
! ! Neder's original code
! !**********************************************************************
! !
! SUBROUTINE four_strucf_serial (iscat, lform) 
! !+                                                                      
! !     Here the complex structure factor of 'nxat' identical atoms       
! !     from array 'xat' is computed.                                     
! !
! !     The phase "iarg0" is calculated via integer math as offset from 
! !     phase = 0 at hkl=0.
! !
! !  Apparently a good OMP version. Very good serial version as well.
! !  The outer loop is j=0,num(1)-1 and the threadprivate variable 
! !  tcsfp is reduced to size (0:num(2)*num(3)-1)
! !  To reduce the number of calculations within the Fourier loop, the
! !  increments xinc*, oinc*, iinc*, jinc* are calculated once up front.
! !
! !  This algorithm is roughly twice as fast on an optimized serial compilation
! !  compared to the original serial Fourier calculation. Thus its serial form
! !  replaces the old serial algorithm as well.
! !-                                                                      
! USE discus_config_mod 
! USE diffuse_mod 
! USE precision_mod
! !
! IMPLICIT none 
! !                                                                       
! INTEGER, INTENT(IN) :: iscat 
! LOGICAL, INTENT(IN) :: lform 
! !                                                                       
! REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        xincu, xincv , xincw
! REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        oincu, oincv , oincw
! INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               iincu, iincv, iincw
! INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               jincu, jincv, jincw
! INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, jj
! INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
! INTEGER (KIND=PREC_INT_LARGE)   :: num23,num123, num123_1
! !
! !INTEGER :: IAND, ISHFT 
! !  
! ! Neder's original code
! !
! !COMPLEX(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP
! !
! ! My declaration
! !
! COMPLEX(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP
! !
! !------ zero fourier array                                              
! !
! ! Neder's original code
! !                                                                       
! !tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
! !
! ! My declaration
! !
! DO i = 1, num (1)
!    DO j = 1, num (2)
!       DO k = 1, num (3)
!          tcsf(i,j,k) = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))
!       ENDDO
!    ENDDO
! ENDDO
! !
! num23    =        num(2)*num(3)
! num123   = num(1)*num(2)*num(3)
! num123_1 = num(1)*num(2)*num(3)-1
! !                                                                       
! !------ Loop over all atoms in 'xat'                                    
! !                                                                       
! DO k = 1, nxat 
!    xincu(k) = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
!    xincv(k) = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
!    xincw(k) = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
!    oincu(k) = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
!    oincv(k) = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
!    oincw(k) = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
! !                                                                       
!    iincu(k) = NINT (64 * I2PI * (xincu(k) - INT (xincu(k)) + 0.0d0) ) 
!    iincv(k) = NINT (64 * I2PI * (xincv(k) - INT (xincv(k)) + 0.0d0) ) 
!    iincw(k) = NINT (64 * I2PI * (xincw(k) - INT (xincw(k)) + 0.0d0) ) 
!    jincu(k) = NINT (64 * I2PI * (oincu(k) - INT (oincu(k)) + 0.0d0) ) 
!    jincv(k) = NINT (64 * I2PI * (oincv(k) - INT (oincv(k)) + 0.0d0) ) 
!    jincw(k) = NINT (64 * I2PI * (oincw(k) - INT (oincw(k)) + 0.0d0) ) 
! ENDDO
! !                                                           
! !------ Loop over all atoms in 'xat'                                    
! !                                                                       
!    IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
!    ALLOCATE (tcsfp (0:       num(2)*num(3)-1)) !,0:nthreads-1))
!    tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
!       DO j = 0, num (1) - 1
!          jj = j*num(2)*num(3)
!          tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
!          loop_k: DO k = 1, nxat 
!             iarg0 =  lmn(1)*iincu(k) + lmn(2)*iincv(k) + lmn(3)*iincw(k) + &
!                      lmn(4)*jincu(k) + lmn(5)*jincv(k) + lmn(6)*jincw(k) + &
!                      j*iincu(k)
!             ii = 0 
!             DO i = 0, num (2) - 1
!                iarg = iarg0 +              iincv(k)*i 
!                DO h = 1, num (3) 
!                   tcsfp(ii)= tcsfp(ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )
!                   iarg     = iarg + iincw(k)
!                   ii       = ii + 1 
!                ENDDO 
!             ENDDO 
!          ENDDO  loop_k
!          DO ii = 1,        num(2)*num(3)
!             tcsf(jj+ii) = tcsf(jj+ii) + tcsfp(ii-1)
!          ENDDO
!       ENDDO 
!    DEALLOCATE(tcsfp)
! !
! !------ Now we multiply with formfactor                                 
! !                                                                       
! IF (lform) then 
!    !DO  i = 1, num (1) * num (2) * num(3)
! !  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
!       !tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                            ! Neder's original code.
!    DO i = 1, num (1)
!       DO j = 1, num (2)
!          DO k = 1, num (3)
!             tcsf (i,j,k) = tcsf (i,j,k) * cfact (istl (i*j*k), iscat)          ! My declaration
!          ENDDO
!       ENDDO
!    ENDDO   
! !  END FORALL
! !   END DO
! ENDIF 
! !                                                                       
! END SUBROUTINE four_strucf_serial
!**********************************************************************
! My code ( Copilot generated )
!**********************************************************************
!
SUBROUTINE four_strucf_serial (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!
!  Apparently a good OMP version. Very good serial version as well.
!  The outer loop is j=0,num(1)-1 and the threadprivate variable 
!  tcsfp is reduced to size (0:num(2)*num(3)-1)
!  To reduce the number of calculations within the Fourier loop, the
!  increments xinc*, oinc*, iinc*, jinc* are calculated once up front.
!
!  This algorithm is roughly twice as fast on an optimized serial compilation
!  compared to the original serial Fourier calculation. Thus its serial form
!  replaces the old serial algreal_toniplorithm as well.
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
INTEGER (KIND=PREC_INT_LARGE)   :: num23,num123, num123_1
!
!INTEGER :: IAND, ISHFT 
!  
! Neder's original code
!
!COMPLEX(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP
!
! My declaration
!
COMPLEX(KIND=PREC_DP), DIMENSION(  :,:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP
!
!------ zero fourier array                                              
!
! Neder's original code
!                                                                       
!tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
!
! My declaration
!
!ALLOCATE(tcsf(num(1), num(2), num(3)))
tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))
!
num23    =        num(2)*num(3)
num123   = num(1)*num(2)*num(3)
num123_1 = num(1)*num(2)*num(3)-1
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
DO k = 1, nxat 
   xincu(k) = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
   xincv(k) = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
   xincw(k) = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
   oincu(k) = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
   oincv(k) = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
   oincw(k) = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
   iincu(k) = NINT (64 * I2PI * (xincu(k) - INT (xincu(k)) + 0.0d0) ) 
   iincv(k) = NINT (64 * I2PI * (xincv(k) - INT (xincv(k)) + 0.0d0) ) 
   iincw(k) = NINT (64 * I2PI * (xincw(k) - INT (xincw(k)) + 0.0d0) ) 
   jincu(k) = NINT (64 * I2PI * (oincu(k) - INT (oincu(k)) + 0.0d0) ) 
   jincv(k) = NINT (64 * I2PI * (oincv(k) - INT (oincv(k)) + 0.0d0) ) 
   jincw(k) = NINT (64 * I2PI * (oincw(k) - INT (oincw(k)) + 0.0d0) ) 
ENDDO
!                                                           
!------ Loop over all atoms in 'xat'                                    
!                                                                       
   IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
   ALLOCATE (tcsfp (             1: num(2)  , 1: num(3)  )) !,0:nthreads-1))
   tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
      DO j = 1, num (1)! - 1
         tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
         loop_k: DO k = 1, nxat 
            iarg0 =  lmn(1)*iincu(k) + lmn(2)*iincv(k) + lmn(3)*iincw(k) + &
                     lmn(4)*jincu(k) + lmn(5)*jincv(k) + lmn(6)*jincw(k) + &
                     (j-1)*iincu(k)
            ii = 0 
            DO i = 1, num (2)!- 1
               iarg = iarg0 +              iincv(k)*(i-1) 
               DO h = 1, num (3) 
!IAND  (ISHFT(iarg, shift), MASK)
!                 tcsfp(   i,h)= tcsfp(   i,h) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                  tcsf (j, i,h)= tcsf (j, i,h) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                  iarg     = iarg + iincw(k)
                  ii       = ii + 1 
               ENDDO 
            ENDDO 
         ENDDO  loop_k
!        DO ii = 1,        num(2)*num(3)
!           tcsf(j+1,ii/num(3)+1,MOD(ii,num(3))+1) = tcsf(j+1,ii/num(3)+1,MOD(ii,num(3))+1) + tcsfp(ii/num(3),MOD(ii,num(3)),j)
!        ENDDO
!        do i=1, num(2)
!        do h=1, num(3)
!           tcsf(j,i,h) = tcsf(j,i,h) + tcsfp(i,h)
!        enddo
!        enddo
      ENDDO 
   DEALLOCATE(tcsfp)
!
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
!ii = 0
   DO i = 1, num (1)
      DO j = 1, num (2)
         DO k = 1, num (3)
!            ii = ii + 1
            tcsf (i,j,k) = tcsf (i,j,k) * cfact (istl (i,j,k), iscat)          ! My declaration
         ENDDO
      ENDDO
   ENDDO   
ENDIF 
!                                                                       
END SUBROUTINE four_strucf_serial













































!*******************************************************************************
!*******************************************************************************
!
SUBROUTINE four_strucf_omp_a(iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!
!  OMP Version in which the Fourier loops has been replaced by a single linear loop
!  Indices j,i,h need to be calculated from loop index ii.
!  Scales nicely linear with thread number N by is not N times faster than serial,
!  As the total number of multiplications is higher.
!-                                                                      
!$ USE OMP_LIB
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)           , DIMENSION(nxat) ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE), DIMENSION(nxat) ::               jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE)   :: iarg
INTEGER (KIND=PREC_INT_LARGE)   :: num23,num123, num123_1
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
INTEGER :: ii, h,i,j, k
!INTEGER :: TID
!
!COMPLEX(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP     ! Neder's original code.
!
COMPLEX(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP      !    My declaration.
!
!$OMP THREADPRIVATE(tcsfp)
!
INTEGER :: IAND, ISHFT 
!
num23    =        num(2)*num(3)
num123   = num(1)*num(2)*num(3)
num123_1 = num(1)*num(2)*num(3)-1
!
!------ zero fourier array                                              
!                                                                       
!tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))                                                    ! Neder's original code.
!
DO i = 1, num (1)
   DO j = 1, num (2)
      DO k = 1, num (3)
         tcsf(i,j,k) = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))                                     ! My code.
      ENDDO
   ENDDO
ENDDO
!    
!****************************************************************************************************************
! PENDING !!!!!
!****************************************************************************************************************                                                                   
!------ Loop over all atoms in 'xat'                                    
!                                                                       
!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC)
DO k = 1, nxat 
   xincu(k) = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
   xincv(k) = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
   xincw(k) = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
   oincu(k) = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
   oincv(k) = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
   oincw(k) = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
   iincu(k) = NINT (64 * I2PI * (xincu(k) - INT (xincu(k)) + 0.0d0) ) 
   iincv(k) = NINT (64 * I2PI * (xincv(k) - INT (xincv(k)) + 0.0d0) ) 
   iincw(k) = NINT (64 * I2PI * (xincw(k) - INT (xincw(k)) + 0.0d0) ) 
   jincu(k) = NINT (64 * I2PI * (oincu(k) - INT (oincu(k)) + 0.0d0) ) 
   jincv(k) = NINT (64 * I2PI * (oincv(k) - INT (oincv(k)) + 0.0d0) ) 
   jincw(k) = NINT (64 * I2PI * (oincw(k) - INT (oincw(k)) + 0.0d0) ) 
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!$OMP PARALLEL PRIVATE(j, i, h, iarg)
      IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
      ALLOCATE (tcsfp (0:num(1)-1, 0:num(2)-1, 0:num(3)-1)) !,0:nthreads-1))
      tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
   !!!$   tid = OMP_GET_THREAD_NUM()
      !$OMP DO SCHEDULE(STATIC)
         DO ii = 0,num123_1
            j =          ii   / num23
            i = MOD(     ii   /num(3) , num(2))
            h = MOD(     ii           , num(3))
            DO k = 1, nxat 
              iarg  =  (j+lmn(1))*iincu(k) + (i+lmn(2))*iincv(k) + (h+lmn(3))*iincw(k) + &
                          lmn(4) *jincu(k) +    lmn(5) *jincv(k) +    lmn(6) *jincw(k)
              tcsfp(i,j,h) = tcsfp(i,j,h) + cex (IAND  (ISHFT(iarg, shift), MASK) )
            ENDDO 
         ENDDO 
      !$OMP END DO
      !$OMP CRITICAL
         DO i = 1, num(1)
            DO j = 1, num(2)
               DO h = 1, num(3)
                  tcsf(i,j,h) = tcsf(i,j,h) + tcsfp(i-1,j-1,h-1)
               ENDDO
            ENDDO
         ENDDO
      !$OMP END CRITICAL
      DEALLOCATE(tcsfp)
!$OMP END PARALLEL
!
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
   !DO  i = 1, num123
!  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
      !tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                         ! Neder's original code.
   DO i = 1, num (1)
      DO j = 1, num (2)
         DO k = 1, num (3)
            tcsf (i,j,k) = tcsf (i,j,k) * cfact (istl (i,j,k), iscat)          ! My declaration
         ENDDO
      ENDDO
   ENDDO
!  END FORALL
   !END DO
ENDIF 
!                                                                       
END SUBROUTINE four_strucf_omp_a
!
!****************************************************************************************************************
!****************************************************************************************************************
!
SUBROUTINE four_strucf_omp_b (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!
!  Simpleset OMP version.  Very close  to the serial algorithm.
!  The outer loop is still the one over all atoms k =1, nxat
!  Although tcsfp is made threadparallel, the calculation time increases
!  very fast if num(1)=num2)=num(3) > some 61 
!  Apparently the 'large' array tcsfp causes cache(?) issues.
!  tcsfp is of size (0:num(1)*num(2)*num(3))
!-                                                                      
!$ USE OMP_LIB
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)        ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)        ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE)   ::                              jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
INTEGER (KIND=PREC_INT_LARGE)   :: num23,num123, num123_1
!
INTEGER :: IAND, ISHFT 
!
!COMPLEX(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE, SAVE :: tcsfp     ! Partial structure factor from parallel OMP          ! Neder's original code.
!
COMPLEX(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: tcsfp  ! Partial structure factor from parallel OMP          ! My declaration.
!
!$OMP THREADPRIVATE(tcsfp)
!
!------ zero fourier array                                              
!                                                                       
!tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))                                    ! Neder's original code.                 
!
DO i = 1, num (1)
   DO j = 1, num (2)
      DO k = 1, num (3)
         tcsf(i,j,k) = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))                     ! My declaration.
      ENDDO
   ENDDO
ENDDO
!
num23    =        num(2)*num(3)
num123   = num(1)*num(2)*num(3)
num123_1 = num(1)*num(2)*num(3)-1
!                                     
!****************************************************************************************************************
! PENDING !!!!!
!****************************************************************************************************************                                  
!------ Loop over all atoms in 'xat'                                    
!                                                                       
!$OMP PARALLEL PRIVATE(xincu, xincv, xincw, iincu, iincv, iincw, &
!$OMP                  oincu, oincv, oincw, jincu, jincv, jincw, &
!$OMP                  iarg, iarg0, ii)
   IF(ALLOCATED(tcsfp)) DEALLOCATE(tcsfp)
   ALLOCATE (tcsfp (0:num(1)-1, 0:num(2)-1, 0:num(3)-1))
   tcsfp = CMPLX(0.0d0, 0.0d0, KIND=KIND(0.0D0))
   !$OMP DO SCHEDULE(STATIC)
      DO k = 1, nxat 
         xincu = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
         xincv = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
         xincw = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
         oincu = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
         oincv = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
         oincw = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
         iincu = NINT (64 * I2PI * (xincu - INT (xincu) + 0.0d0) ) 
         iincv = NINT (64 * I2PI * (xincv - INT (xincv) + 0.0d0) ) 
         iincw = NINT (64 * I2PI * (xincw - INT (xincw) + 0.0d0) ) 
         jincu = NINT (64 * I2PI * (oincu - INT (oincu) + 0.0d0) ) 
         jincv = NINT (64 * I2PI * (oincv - INT (oincv) + 0.0d0) ) 
         jincw = NINT (64 * I2PI * (oincw - INT (oincw) + 0.0d0) ) 
         iarg0 =  lmn(1)*iincu + lmn(2)*iincv + lmn(3)*iincw + &
                  lmn(4)*jincu + lmn(5)*jincv + lmn(6)*jincw
         iarg = iarg0 
!                                                                       
         ii = 0 
         DO j = 0, num (1) - 1
            DO i = 0, num (2) - 1
               iarg = iarg0 + iincu*j + iincv*i 
               DO h = 1, num (3) 
                  tcsfp(i,j,h)= tcsfp(i,j,h) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                  iarg     = iarg + iincw
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
   !$OMP END DO
   !$OMP CRITICAL
      DO k = 1, num(3)
         DO j = 1, num(2)
            DO i = 1, num(1)
               tcsf(i,j,k) = tcsf(i,j,k) + tcsfp(i-1,j-1,k-1)
            ENDDO
         ENDDO
      ENDDO

   !$OMP END CRITICAL
   DEALLOCATE(tcsfp)
!$OMP END PARALLEL
!
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
   !DO  i = 1, num (1) * num (2) * num(3)
!  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
      !tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                                     ! Neder's original code.
   DO i = 1, num (1)
      DO j = 1, num (2)
         DO k = 1, num (3)
            tcsf (i,j,k) = tcsf (i,j,k) * cfact (istl (i,j,k), iscat)          ! My declaration
         ENDDO
      ENDDO
   ENDDO
!  END FORALL
   !END DO
ENDIF 
!                                                                       
END SUBROUTINE four_strucf_omp_b
!
!****************************************************************************************************************
!****************************************************************************************************************
!
SUBROUTINE four_strucf_serial_original (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!-                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: iscat 
LOGICAL, INTENT(IN) :: lform 
!                                                                       
REAL(KIND=PREC_DP)        ::        xincu, xincv , xincw
REAL(KIND=PREC_DP)        ::        oincu, oincv , oincw
INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw
INTEGER (KIND=PREC_INT_LARGE)   ::                              jincu, jincv, jincw
INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
!
INTEGER :: IAND, ISHFT 
!
!------ zero fourier array                                              
!                                                                       
!tcsf = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))                                    ! Neder's original code.
!
DO i = 1, num (1)
   DO j = 1, num (2)
      DO k = 1, num (3)
         tcsf(i,j,k) = CMPLX(0.0D0, 0.0D0, KIND=KIND(0.0D0))                     ! My declaration.
      ENDDO
   ENDDO
ENDDO
!
!****************************************************************************************************************
! Neders original code
! !****************************************************************************************************************                 !                       
! !------ Loop over all atoms in 'xat'                                    
! !                                                                       
! DO k = 1, nxat 
! !        xarg0 = xm (1)        * xat(k, 1) + xm (2)         * xat(k, 2) + xm (3)         * xat(k, 3)
!    xincu = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
!    xincv = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
!    xincw = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
!    oincu = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
!    oincv = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
!    oincw = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
! !                                                                       
!    iincu = NINT (64 * I2PI * (xincu - INT (xincu) + 0.0d0) ) 
!    iincv = NINT (64 * I2PI * (xincv - INT (xincv) + 0.0d0) ) 
!    iincw = NINT (64 * I2PI * (xincw - INT (xincw) + 0.0d0) ) 
!    jincu = NINT (64 * I2PI * (oincu - INT (oincu) + 0.0d0) ) 
!    jincv = NINT (64 * I2PI * (oincv - INT (oincv) + 0.0d0) ) 
!    jincw = NINT (64 * I2PI * (oincw - INT (oincw) + 0.0d0) ) 
!    iarg0 =  lmn(1)*iincu + lmn(2)*iincv + lmn(3)*iincw + &
!             lmn(4)*jincu + lmn(5)*jincv + lmn(6)*jincw
!    iarg = iarg0 
! !                                                                       
! !------ - Loop over all points in Q. 'iadd' is the address of the       
! !------ - complex exponent table. 'IADD' divides out the 64 and         
! !------ - ISHFT acts as MOD so that the argument stays in the table     
! !------ - boundaries.                                                   
! !                 iadd      = ISHFT (iarg, - 6) 
! !                 iadd      = IAND  (iadd, MASK) 
! !                 tcsf (ii) = tcsf (ii) + cex (iadd, MASK) )
! !                 iarg      = iarg + iincw
! !                                                                       
!    ii = 0 
! !                                                                       
!    DO j = 0, num (1) - 1
!       DO i = 0, num (2) - 1
!          iarg = iarg0 + iincu*j + iincv*i 
!          DO h = 1, num (3) 
!             ii       = ii + 1 
!             tcsf(ii) = tcsf (ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )
!             iarg     = iarg + iincw
!          ENDDO 
!       ENDDO 
!    ENDDO 
! ENDDO 
!
!****************************************************************************************************************
! My code  ( Copilot generated )
!****************************************************************************************************************                 !------ Loop over all atoms in 'xat'                                    
!                                                                       
DO k = 1, nxat 
!        xarg0 = xm (1)        * xat(k, 1) + xm (2)         * xat(k, 2) + xm (3)         * xat(k, 3)
   xincu = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
   xincv = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
   xincw = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
   oincu = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
   oincv = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
   oincw = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
   iincu = NINT (64 * I2PI * (xincu - INT (xincu) + 0.0d0) ) 
   iincv = NINT (64 * I2PI * (xincv - INT (xincv) + 0.0d0) ) 
   iincw = NINT (64 * I2PI * (xincw - INT (xincw) + 0.0d0) ) 
   jincu = NINT (64 * I2PI * (oincu - INT (oincu) + 0.0d0) ) 
   jincv = NINT (64 * I2PI * (oincv - INT (oincv) + 0.0d0) ) 
   jincw = NINT (64 * I2PI * (oincw - INT (oincw) + 0.0d0) ) 
   iarg0 =  lmn(1)*iincu + lmn(2)*iincv + lmn(3)*iincw + &
            lmn(4)*jincu + lmn(5)*jincv + lmn(6)*jincw
   iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                 iadd      = ISHFT (iarg, - 6) 
!                 iadd      = IAND  (iadd, MASK) 
!                 tcsf (ii) = tcsf (ii) + cex (iadd, MASK) )
!                 iarg      = iarg + iincw
!                                                                       
   ii = 0 
!                                                                       
   DO j = 0, num (1) - 1
      DO i = 0, num (2) - 1
         iarg = iarg0 + iincu*j + iincv*i 
         DO h = 1, num (3) 
            ii       = ii + 1 
            tcsf(j+1, i+1, h) = tcsf(j+1, i+1, h) + cex (IAND  (ISHFT(iarg, shift), MASK) )
            iarg     = iarg + iincw
         ENDDO 
      ENDDO 
   ENDDO 
ENDDO 
!
!****************************************************************************************************************
!****************************************************************************************************************
!------ Now we multiply with formfactor                                 
!                                                                       
IF (lform) then 
   !DO  i = 1, num (1) * num (2) * num(3)
!  FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
   !   tcsf (i) = tcsf (i) * cfact (istl (i), iscat)                                         ! Neder's original code.
!   k = 0
   DO i = 1, num (1)
      DO j = 1, num (2)
         DO h = 1, num (3)
!            k = k + 1
!           tcsf (i,j,h) = tcsf (i,j,h) * cfact (istl (i*j*k), iscat)          ! My declaration
            tcsf (i,j,h) = tcsf (i,j,h) * cfact (istl (i,j,h), iscat)          ! My declaration
         ENDDO
      ENDDO
   ENDDO
!  END FORALL
   !END DO
ENDIF 
!                                                                       
END SUBROUTINE four_strucf_serial_original
!
END MODULE four_strucf_mod
