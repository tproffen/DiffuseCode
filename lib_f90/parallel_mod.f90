MODULE parallel_mod
!-
!   Variables related to parallel processing using OMP
!+
LOGICAL   :: par_omp_use = .TRUE.      ! User wants to use OMP
INTEGER   :: par_omp_maxthreads = -1   ! Maximum number of threads to use, -1==use all available
!
END MODULE parallel_mod
