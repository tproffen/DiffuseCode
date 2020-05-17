MODULE refine_data_mod
!
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=PREC_STRING)           :: ref_load  = ' '  ! Load string data
CHARACTER(LEN=PREC_STRING)           :: ref_csigma = ' ' ! Load string Sigma's
INTEGER                              :: ref_kload = 0    ! Data set within KUPLOT
INTEGER                              :: ref_ksigma= 0    ! Sigma set within KUPLOT
INTEGER                              :: ref_kupl  = 0    ! Data set within KUPLOT that needs to be kept
INTEGER, DIMENSION(2)                :: ref_dim    ! Dimensions of data set
REAL   , DIMENSION(:,:), ALLOCATABLE :: ref_data   ! the actual data set
REAL   , DIMENSION(:,:), ALLOCATABLE :: ref_sigma  ! sigma at each data point
REAL   , DIMENSION(:  ), ALLOCATABLE :: ref_x      ! x-values of data set
REAL   , DIMENSION(:  ), ALLOCATABLE :: ref_y      ! y-values of data set
!
END MODULE refine_data_mod
