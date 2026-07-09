module refine_data_mod
!
use lib_data_struc_type_mod
use precision_mod
!
implicit none
!
character(LEN=PREC_STRING)                                :: ref_load_u  = ' '  ! Load string data
character(LEN=PREC_STRING)                                :: ref_csigma_u = ' ' ! Load string Sigma's
character(LEN=PREC_STRING)                                :: ref_load  = ' '  ! Load string data
character(LEN=PREC_STRING)                                :: ref_csigma = ' ' ! Load string Sigma's
character(len=PREC_STRING), dimension(:)      , allocatable :: ref_data_name    ! Data set names
integer                   , dimension(:)      , allocatable :: ref_data_number  ! Data set Rnumber in global
integer                                                   :: ref_ndata = 0    ! Number of data set for REFINE
integer                                                   :: ref_maxdata = 0  ! Max Number of data set for REFINE
integer                                                   :: ref_kload = 0    ! Data set within KUPLOT
integer                                                   :: ref_ksigma= 0    ! Sigma set within KUPLOT
integer                                                   :: ref_kupl  = 0    ! Data set within KUPLOT that needs to be kept
integer                                                   :: ref_type  = 0    ! Data type => lib_data_types_mod
real(kind=PREC_DP)        , dimension(:,:,:,:), allocatable :: ref_data   ! the actual data set
real(kind=PREC_DP)        , dimension(:,:,:,:), allocatable :: ref_sigma  ! sigma at each data point
real(kind=PREC_DP)        , dimension(      :), allocatable :: ref_weight ! Weight for each data set, defaults to 1 
real(kind=PREC_DP)        , dimension(:    ,:), allocatable :: ref_x      ! x-values of data set
real(kind=PREC_DP)        , dimension(:    ,:), allocatable :: ref_y      ! y-values of data set
real(kind=PREC_DP)        , dimension(:    ,:), allocatable :: ref_z      ! z-values of data set
!
type ref_ptr
   type(h5_data_struc), pointer :: data_ptr
end type ref_ptr 
!
type(ref_ptr), dimension(:), allocatable :: ref_data_ptr
type(ref_ptr), dimension(:), allocatable :: ref_calc_ptr
!
end module refine_data_mod
