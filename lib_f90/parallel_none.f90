MODULE parallel_mod
!-
!   Variables related to parallel processing using CUDA
!+
LOGICAL   :: par_omp_use = .FALSE.     ! User does not want to use OMP
INTEGER   :: par_omp_maxthreads =  1   ! Maximum number of threads to use, 1
INTEGER   :: par_omp_phys       =  1   ! Maximum number of threads to use, 1
INTEGER   :: par_omp_logi       =  1   ! Maximum number of threads to use, 1
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE get_cores()
!
END SUBROUTINE get_cores
!
!*******************************************************************************
!
END MODULE parallel_mod
