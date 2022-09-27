module lib_data_struc_type_mod
!
use precision_mod
!
private
public h5_data_struc
!
type :: h5_data_struc
   integer                                               :: data_num       ! Current data set number
   character(len=PREC_STRING)                            :: infile         ! input file
   integer                                               :: layer=1        ! Current layer in data set
   LOGICAL                                               :: is_direct      ! Direct space == TRUE
   integer                                               :: ndims          ! Number of dimensions
   integer              , dimension(3)                   :: dims           ! Actual dimensions
   logical                                               :: is_grid        ! Data have grid like coordinates
   logical                                               :: has_dxyz       ! Data have uncertainties in xyz
   logical                                               :: has_dval       ! Data have uncertainties in value
   real(kind=PREC_DP)   , dimension(3,4)                 :: corners ! (3,4)
   real(kind=PREC_DP)   , dimension(3,3)                 :: vectors
   real(kind=PREC_DP)   , dimension(3)                   :: cr_a0
   real(kind=PREC_DP)   , dimension(3)                   :: cr_win
   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: x              ! Actual x-coordinates
   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: y              ! Actual y-coordinates
   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: z              ! Actual z-coordinates
   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: dx             ! Actual x-coordinates uncertainties
   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: dy             ! Actual y-coordinates uncertainties
   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: dz             ! Actual z-coordinates uncertainties
   real(kind=PREC_DP)   , dimension(:,:,:), allocatable  :: datamap        ! Actual diffraction data
   real(kind=PREC_DP)   , dimension(:,:,:), allocatable  :: sigma          ! Actual diffraction data uncertainties
   real(kind=PREC_DP)   , dimension(3)                   :: llims          ! Lower limits
   real(kind=PREC_DP)   , dimension(3)                   :: steps          ! steps in H, K, L
   real(kind=PREC_DP)   , dimension(3,3)                 :: steps_full     ! steps in H, K, L
   real(kind=PREC_DP)   , dimension(2)                   :: minmaxval      ! data extreme values
   real(kind=PREC_DP)   , dimension(3,2)                 :: minmaxcoor     ! coordinates extreme values
   type(h5_data_struc), pointer                          :: after
end type h5_data_struc
!
end module lib_data_struc_type_mod
