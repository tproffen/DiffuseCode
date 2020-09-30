MODULE parallel_mod
!-
!   Variables related to parallel processing using OMP
!+
LOGICAL   :: par_omp_use = .FALSE.     ! User does not want to use OMP
INTEGER   :: par_omp_maxthreads =  1   ! Maximum number of threads to use, 1
!
!*******************************************************************************
!
SUBROUTINE get_cores()
!
END SUBROUTINE get_cores()
!
!*******************************************************************************
!
END MODULE parallel_mod
