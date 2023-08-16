module lib_ik_mod
!-
!  Variables used by lib_math and refine for data exchange between global storage and local storage
!+
!
use precision_mod
!
implicit none
!
character(len=PREC_STRING)                         :: ik1_infile
character(len=PREC_STRING)                         :: ik2_infile
integer                                            :: ik1_node_number  ! Node in global data
integer                                            :: ik2_node_number  ! Node in global data
integer                                            :: ik1_data_type    ! Data type 
integer                                            :: ik2_data_type    ! Data type
integer                                            :: ik1_nlayer       ! Current layer in display
integer                                            :: ik2_nlayer       ! Current layer in display
logical                                            :: ik1_is_direct   ! Data are on direct / reciprocal scale
logical                                            :: ik2_is_direct   ! Data are on direct / reciprocal scale
logical                                            :: ik1_is_grid     ! Data are on direct / reciprocal scale
logical                                            :: ik2_is_grid     ! Data are on direct / reciprocal scale
logical                                            :: ik1_has_dxyz    ! Data are on direct / reciprocal scale
logical                                            :: ik2_has_dxyz    ! Data are on direct / reciprocal scale
logical                                            :: ik1_has_dval    ! Data are on direct / reciprocal scale
logical                                            :: ik2_has_dval    ! Data are on direct / reciprocal scale
logical                                            :: ik1_calc_coor   ! Need to calculate coordinates
logical                                            :: ik2_calc_coor   ! Need to calculate coordinates
integer, dimension(3)                              :: ik1_use_coor    ! Use these axes for coordinates
integer, dimension(3)                              :: ik2_use_coor    ! Use these axes for coordinates
integer                                            :: ik1_ndims        ! Number of dimensions
integer                                            :: ik2_ndims        ! Number of dimensions
integer, dimension(3)                              :: ik1_dims         ! Dimensions global array
integer, dimension(3)                              :: ik2_dims         ! Dimensions global array
real(kind=PREC_DP), dimension(3)                   :: ik1_a0           ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: ik2_a0           ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: ik1_win          ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: ik2_win          ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: ik1_llims        ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: ik2_llims        ! Lower limits global array
real(kind=PREC_DP), dimension(3,4)                 :: ik1_corners      ! Steps        global array
real(kind=PREC_DP), dimension(3,4)                 :: ik2_corners      ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: ik1_vectors      ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: ik2_vectors      ! Steps        global array
real(kind=PREC_DP), dimension(3)                   :: ik1_steps        ! Steps        global array
real(kind=PREC_DP), dimension(3)                   :: ik2_steps        ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: ik1_steps_full   ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: ik2_steps_full   ! Steps        global array
real(kind=PREC_DP), dimension(:)    , target, allocatable  :: ik1_x           ! Global data array for real
real(kind=PREC_DP), dimension(:)    , target, allocatable  :: ik1_y           ! Global data array for real
real(kind=PREC_DP), dimension(:)    , target, allocatable  :: ik1_z           ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: ik2_x           ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: ik2_y           ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: ik2_z           ! Global data array for real
real(kind=PREC_DP), dimension(:)    , target, allocatable  :: ik1_dx          ! Global data array for real
real(kind=PREC_DP), dimension(:)    , target, allocatable  :: ik1_dy          ! Global data array for real
real(kind=PREC_DP), dimension(:)    , target, allocatable  :: ik1_dz          ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: ik2_dx          ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: ik2_dy          ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: ik2_dz          ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), target, allocatable  :: ik1_data        ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: ik2_data        ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), target, allocatable  :: ik1_sigma       ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: ik2_sigma       ! Global data array for real
real(kind=PREC_DP), dimension(2)                   :: ik1_minmaxval    ! Steps        global array
real(kind=PREC_DP), dimension(2)                   :: ik2_minmaxval    ! Steps        global array
real(kind=PREC_DP), dimension(3,2)                 :: ik1_minmaxcoor   ! Steps        global array
real(kind=PREC_DP), dimension(3,2)                 :: ik2_minmaxcoor   ! Steps        global array
!
end module lib_ik_mod
