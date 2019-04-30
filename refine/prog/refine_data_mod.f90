MODULE refine_data_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1024)                  :: ref_load   ! Load string
INTEGER                              :: ref_kupl   ! Data set within KUPLOT
INTEGER, DIMENSION(2)                :: ref_dim    ! Dimensions of data set
REAL   , DIMENSION(:,:), ALLOCATABLE :: ref_data   ! the actual data set
REAL   , DIMENSION(:,:), ALLOCATABLE :: ref_weight ! weights at each data point
REAL   , DIMENSION(:  ), ALLOCATABLE :: ref_x      ! x-values of data set
REAL   , DIMENSION(:  ), ALLOCATABLE :: ref_y      ! y-values of data set
!
END MODULE refine_data_mod
