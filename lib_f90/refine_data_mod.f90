MODULE refine_data_mod
!
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=PREC_STRING)                                :: ref_load_u  = ' '  ! Load string data
CHARACTER(LEN=PREC_STRING)                                :: ref_csigma_u = ' ' ! Load string Sigma's
CHARACTER(LEN=PREC_STRING)                                :: ref_load  = ' '  ! Load string data
CHARACTER(LEN=PREC_STRING)                                :: ref_csigma = ' ' ! Load string Sigma's
integer                                                   :: ref_ndata = 0    ! Number of data set for REFINE
INTEGER                                                   :: ref_kload = 0    ! Data set within KUPLOT
INTEGER                                                   :: ref_ksigma= 0    ! Sigma set within KUPLOT
INTEGER                                                   :: ref_kupl  = 0    ! Data set within KUPLOT that needs to be kept
INTEGER                                                   :: ref_type  = 0    ! Data type => lib_data_types_mod
INTEGER                   , DIMENSION(3,1)                  :: ref_dim    ! Dimensions of data set
REAL(kind=PREC_DP)        , DIMENSION(:,:,:,:), ALLOCATABLE :: ref_data   ! the actual data set
REAL(kind=PREC_DP)        , DIMENSION(:,:,:,:), ALLOCATABLE :: ref_sigma  ! sigma at each data point
REAL(kind=PREC_DP)        , DIMENSION(:    ,:), ALLOCATABLE :: ref_x      ! x-values of data set
REAL(kind=PREC_DP)        , DIMENSION(:    ,:), ALLOCATABLE :: ref_y      ! y-values of data set
REAL(kind=PREC_DP)        , DIMENSION(:    ,:), ALLOCATABLE :: ref_z      ! z-values of data set
!
END MODULE refine_data_mod
