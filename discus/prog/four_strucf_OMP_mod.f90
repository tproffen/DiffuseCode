
#include "debug.h"

MODULE four_strucf_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE four_strucf (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
      USE omp_lib
      USE config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      REAL(KIND=8)        :: xarg0, xincu, xincv , xincw
      INTEGER             :: h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw, iadd 
!
      CHARACTER(len=16) :: omp_num_threads_env          ! Environment variable OMP_NUM_THREADS
      INTEGER :: omp_num_threads                        ! Number of threads
      INTEGER                              :: tid       ! Id of this thread
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: tcsfp     ! Partial structure factor from parallel OMP
!                                                                       
      INTEGER IAND, ISHFT 


      MSG('four_strucf OMP')
!
!------ zero fourier array                                              
!                                                                       
      tcsf = cmplx (0.0d0, 0.0d0) 

!     OMP_NUM_THREADS doesn't work for me. -Justin
!     Here, we handle this ourselves:
      CALL GET_ENVIRONMENT_VARIABLE("OMP_NUM_THREADS", &
                                     omp_num_threads_env)
      read(omp_num_threads_env, *), omp_num_threads
      CALL OMP_SET_NUM_THREADS(omp_num_threads)

!     Allocate, initialize tcsfp
      VAR(MAXQXY)
      ALLOCATE (tcsfp (1:MAXQXY,0:omp_num_threads-1))
      tcsfp = cmplx(0.0d0, 0.0d0)

      VAR(omp_num_threads)
      CALL OMP_SET_NUM_THREADS(omp_num_threads)
!$OMP PARALLEL PRIVATE(tid,k,xarg0,xincu,xincv,xincw,iincu,iincv,iincw,iarg,iarg0,ii,j,i,h)
!$OMP DO
!------ Loop over all atoms in 'xat'
      DO k = 1, nxat 
         tid = OMP_GET_THREAD_NUM()
         VAR(tid)
         xarg0 = xm (1) * xat(k, 1) + xm (2) * xat(k, 2) + xm  (3) * xat(k, 3)
         xincu = uin(1) * xat(k, 1) + uin(2) * xat(k, 2) + uin (3) * xat(k, 3)
         xincv = vin(1) * xat(k, 1) + vin(2) * xat(k, 2) + vin (3) * xat(k, 3)
         xincw = win(1) * xat(k, 1) + win(2) * xat(k, 2) + win (3) * xat(k, 3)
!                                                                       
         iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
         iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
         iincv = nint (64 * I2PI * (xincv - int (xincv) + 1.0d0) ) 
         iincw = nint (64 * I2PI * (xincw - int (xincw) + 1.0d0) ) 
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
                   ii            = ii + 1 
                   tcsfp(ii,tid) = tcsfp(ii,tid) + cex (IAND  (ISHFT(iarg,-6), MASK) )
                   iarg          = iarg + iincw
                ENDDO 
             ENDDO 
          ENDDO 
      ENDDO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
      tcsf = SUM(tcsfp, DIM=2)
      DEALLOCATE(tcsfp)
!                                                                       
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) then 
         DO  i = 1, num (1) * num (2) * num(3)
!        FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
            tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
!        END FORALL
         END DO
      ENDIF 
!

      MSG('four_strucf OMP DONE.')
      END SUBROUTINE four_strucf                    
END MODULE four_strucf_mod
