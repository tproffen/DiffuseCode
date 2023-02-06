module lib_data_struc_h5
!-
!  Contains routines to build the lib_f90 data structure, based on HDF5
!+
use lib_hdf5_params_mod
use lib_data_struc_type_mod
!
use hdf5_def_mod
use precision_mod
!
implicit none
!
private
public dgl5_new_node
public dgl5_set_node
public dgl5_get_node
public dgl5_set_pointer
public dgl5_find_node
public dgl5_copy_node
public dgl5_get_h5_is_ku
public dgl5_get_ku_is_h5
public dgl5_get_layer
public dgl5_get_height
public dgl5_get_direct
public dgl5_get_is_grid
public dgl5_get_ndims
public dgl5_get_infile
public dgl5_get_dims
public dgl5_get_lattice
public dgl5_get_llims
public dgl5_get_steps
public dgl5_get_minmax
public dgl5_get_minmaxcoor
public dgl5_get_map
public dgl5_get_tmap
public dgl5_get_x
public dgl5_get_y
public dgl5_get_z
public dgl5_get_data
public dgl5_get_number
public dgl5_set_is_grid
public dgl5_set_llims
public dgl5_set_steps
public dgl5_set_h5_is_ku
public dgl5_set_ku_is_h5
public dgl5_set_layer
public dgl5_set_map
public dgl5_set_x
public dgl5_set_y
public dgl5_set_z
public data2local
public local2data
public fft2data
public dgl5_reset
!
integer, parameter ::  maxkurvtot = 200      ! Relic from KUPLOT
!
integer, dimension(MAXKURVTOT)                        :: h5_ku_is_h5 = 0   ! Pointer from kuplot number to h5 number
integer, dimension(MAXKURVTOT)                        :: h5_h5_is_ku = 0   ! Pointer from h5 number to kuplot number
integer                                               :: h5_number   = 0   ! Currently loaded h5 data sets
!
!type :: h5_data_struc
!   integer                                               :: data_num       ! Current data set number
!   character(len=PREC_STRING)                            :: infile         ! input file
!   integer                                               :: layer=1        ! Current layer in data set
!   LOGICAL                                               :: direct         ! Direct space == TRUE
!   integer                                               :: ndims          ! Number of dimensions
!   integer              , dimension(3)                   :: dims           ! Actual dimensions
!   logical                                               :: is_grid        ! Data have grid like coordinates
!!   logical                                               :: has_dxyz       ! Data have uncertainties in xyz
!   logical                                               :: has_dval       ! Data have uncertainties in value
!   real(kind=PREC_DP)   , dimension(3,4)                 :: corners ! (3,4)
!   real(kind=PREC_DP)   , dimension(3,3)                 :: vectors
!   real(kind=PREC_DP)   , dimension(3)                   :: cr_a0
!   real(kind=PREC_DP)   , dimension(3)                   :: cr_win
!   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: x              ! Actual x-coordinates
!   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: y              ! Actual y-coordinates
!   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: z              ! Actual z-coordinates
!   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: dx             ! Actual x-coordinates uncertainties
!   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: dy             ! Actual y-coordinates uncertainties
!   real(kind=PREC_DP)   , dimension(:)    , allocatable  :: dz             ! Actual z-coordinates uncertainties
!   real(kind=PREC_DP)   , dimension(:,:,:), allocatable  :: datamap        ! Actual diffraction data
!   real(kind=PREC_DP)   , dimension(:,:,:), allocatable  :: sigma          ! Actual diffraction data uncertainties
!   real(kind=PREC_DP)   , dimension(3)                   :: llims          ! Lower limits
!   real(kind=PREC_DP)   , dimension(3)                   :: steps          ! steps in H, K, L
!   real(kind=PREC_DP)   , dimension(3,3)                 :: steps_full     ! steps in H, K, L
!   real(kind=PREC_DP)   , dimension(2)                   :: minmaxval      ! data extreme values
!   real(kind=PREC_DP)   , dimension(3,2)                 :: minmaxcoor     ! coordinates extreme values
!   type(h5_data_struc), pointer                          :: after
!end type h5_data_struc
!
type(h5_data_struc), pointer                          :: h5_root => NULL()
type(h5_data_struc), pointer                          :: h5_temp => NULL()
type(h5_data_struc), pointer                          :: h5_find => NULL()
!
!
contains
!
!*******************************************************************************
!
subroutine dgl5_new_node
!-
! Create a new node
!+
!
if(associated(h5_root)) then                                ! A root node exists
   h5_temp => h5_root                                       ! Point to current node
   find_node: do while(associated(h5_temp%after))           ! Does next node exist?
      h5_temp => h5_temp%after                              ! Next node exists, point to this next node
   enddo find_node
   allocate(h5_temp%after)                                  ! Create next node
   h5_temp => h5_temp%after                                 ! Point to Current working node
else
   allocate(h5_root)                                        ! Create first node
   nullify(h5_root%after)
   h5_temp => h5_root                                       ! Point to Current working node
endif
! Work on current node
nullify(h5_temp%after)
h5_temp%data_num = h5_number + 1                         ! Increment the data number
h5_number = h5_number + 1                                   ! Increment the global data number
!
!write(*,*) ' ROOT NW', associated(h5_root), ASSOCIATED(h5_temp), h5_temp%data_num, h5_root%data_num, h5_number
!
end subroutine dgl5_new_node
!
!*******************************************************************************
!
subroutine dgl5_copy_node(old, new)
!-
!  Copies old node to new node
!+
integer, intent(in)  :: old
integer, intent(out) :: new
integer :: ier_num
integer :: ier_typ
!
call dgl5_find_node(old, ier_num, ier_typ)
call dgl5_new_node
h5_temp%infile     = h5_find%infile         ! input file
h5_temp%layer      = h5_find%layer          ! Current layer in data set
h5_temp%is_direct     = h5_find%is_direct         ! Direct space == TRUE
h5_temp%ndims      = h5_find%ndims          ! Number of dimensions
h5_temp%dims       = h5_find%dims           ! Actual dimensions
h5_temp%is_grid    = h5_find%is_grid        ! Data are on periodic grid
h5_temp%has_dxyz   = h5_find%has_dxyz       ! Data have uncertainties on coordinates
h5_temp%has_dval   = h5_find%has_dval       ! Data have uncertainties in valu
h5_temp%corners    = h5_find%corners        ! 
h5_temp%vectors    = h5_find%vectors        ! 
h5_temp%cr_a0      = h5_find%cr_a0
h5_temp%cr_win     = h5_find%cr_win 
h5_temp%x          = h5_find%x              ! Coordinates x
h5_temp%y          = h5_find%y              ! Coordinates x
h5_temp%z          = h5_find%z              ! Coordinates x
h5_temp%dx         = h5_find%dx             ! Coordinates x
h5_temp%dy         = h5_find%dy             ! Coordinates x
h5_temp%dz         = h5_find%dz             ! Coordinates x
h5_temp%datamap    = h5_find%datamap        ! Actual diffraction data
h5_temp%sigma      = h5_find%sigma          ! Actual diffraction data
h5_temp%llims      = h5_find%llims          ! Lower limits
h5_temp%steps      = h5_find%steps          ! steps in H, K, L
h5_temp%steps_full = h5_find%steps_full     ! steps in H, K, L
h5_temp%minmaxval  = h5_find%minmaxval      ! Data extrema
h5_temp%minmaxcoor = h5_find%minmaxcoor     ! Data extrema
!
new = h5_temp%data_num
!
end subroutine dgl5_copy_node
!
!*******************************************************************************
!
subroutine dgl5_set_node(l_infile, l_layer, l_direct, l_ndims, l_dims, &
                   l_is_grid,l_has_dxyz, l_has_dval, l_corners, l_vectors, l_cr_a0, l_cr_win, &
                   l_x, l_y, l_z, l_dx, l_dy, l_dz,                 &
                   l_data, l_sigma, l_llims, l_steps,  &
                   l_steps_full)
!-
!  Place the temporary values into the current dgl5 node
!+
character(len=*)                        , intent(in) :: l_infile         ! Input file
integer                                 , intent(in) :: l_layer          ! Current layer in data set
LOGICAL                                 , intent(in) :: l_direct         ! Direct space == TRUE
integer                                 , intent(in) :: l_ndims          ! Number of dimensions
integer           , dimension(3)        , intent(in) :: l_dims           ! Actual dimensions
logical                                 , intent(in) :: l_is_grid        ! Data on periodic grid
logical                                 , intent(in) :: l_has_dxyz       ! Data on periodic grid
logical                                 , intent(in) :: l_has_dval       ! Data on periodic grid
real(kind=PREC_DP), dimension(3,4)      , intent(in) :: l_corners        ! 
real(kind=PREC_DP), dimension(3,3)      , intent(in) :: l_vectors        ! 
real(kind=PREC_DP), dimension(3)        , intent(in) :: l_cr_a0          ! 
real(kind=PREC_DP), dimension(3)        , intent(in) :: l_cr_win         ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(in) :: l_x              ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(in) :: l_y              ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(in) :: l_z              ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(in) :: l_dx             ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(in) :: l_dy             ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(in) :: l_dz             ! 
real(kind=PREC_DP), dimension(l_dims(1), l_dims(2), l_dims(3)), intent(in):: l_data           ! Actual diffraction data
real(kind=PREC_DP), dimension(:, :, :), allocatable           , intent(in):: l_sigma          ! Actual diffraction data
real(kind=PREC_DP), dimension(3)        , intent(in) :: l_llims          ! Lower limits
real(kind=PREC_DP), dimension(3)        , intent(in) :: l_steps          ! steps in H, K, L
real(kind=PREC_DP), dimension(3,3)      , intent(in) :: l_steps_full     ! steps in H, K, L
!
h5_temp%infile   = l_infile         ! input file
h5_temp%layer    = l_layer          ! Current layer in data set
h5_temp%is_direct   = l_direct         ! Direct space == TRUE
h5_temp%ndims    = l_ndims          ! Number of dimensions
h5_temp%dims     = l_dims           ! Actual dimensions
h5_temp%is_grid  = l_is_grid        ! Actual dimensions
h5_temp%has_dxyz = l_has_dxyz       ! Actual dimensions
h5_temp%has_dval = l_has_dval       ! Actual dimensions
h5_temp%corners  = l_corners        ! Actual dimensions
h5_temp%vectors  = l_vectors        ! Actual dimensions
h5_temp%cr_a0    = l_cr_a0          ! Lattice parameters
h5_temp%cr_win   = l_cr_win         ! lattice angles
!allocate(h5_temp%data(h5_temp%dims(1), h5_temp%dims(2), h5_temp%dims(3)))
!if(.not.l_is_grid) then
h5_temp%x        = l_x              ! Actual x-coordinates
h5_temp%y        = l_y              ! Actual y-coordinates
h5_temp%z        = l_z              ! Actual z-coordinates
!endif
!if(l_has_dxyz) then
   h5_temp%dx        = l_dx              ! Actual x-coordinates uncertainties
   h5_temp%dy        = l_dy              ! Actual y-coordinates uncertainties
   h5_temp%dz        = l_dz              ! Actual z-coordinates uncertainties
!endif
h5_temp%datamap  = l_data           ! Actual diffraction data
if(l_has_dval) then
   h5_temp%sigma    = l_sigma          ! Actual diffraction data
endif
h5_temp%llims    = l_llims          ! Lower limits
h5_temp%steps    = l_steps          ! steps in H, K, L
h5_temp%steps_full    = l_steps_full          ! steps in H, K, L
h5_temp%minmaxval(1) = minval(l_data)
h5_temp%minmaxval(2) = maxval(l_data)
h5_temp%minmaxcoor(1,1) = minval(l_x)
h5_temp%minmaxcoor(1,2) = maxval(l_x)
h5_temp%minmaxcoor(2,1) = minval(l_y)
h5_temp%minmaxcoor(2,2) = maxval(l_y)
h5_temp%minmaxcoor(3,1) = minval(l_z)
h5_temp%minmaxcoor(3,2) = maxval(l_z)
!write(*,*) ' SET NODE NUM, NDIMS, DIMS ', h5_temp%data_num, h5_temp%ndims, h5_temp%dims
!write(*,*) ' X in NODE ', minval(h5_temp%x), maxval(h5_temp%x)

!
end subroutine dgl5_set_node
!
!*******************************************************************************
!
subroutine dgl5_get_node(l_infile, l_layer, l_direct, l_ndims, l_dims, &
                   l_is_grid,l_has_dxyz, l_has_dval, l_corners, l_vectors, l_cr_a0, l_cr_win, &
                   l_x, l_y, l_z, l_dx, l_dy, l_dz,                 &
                   l_data, l_sigma, l_llims, l_steps,  &
                   l_steps_full,l_minmaxval, l_minmaxcoor)
!-
!  Place the temporary values into the current dgl5 node
!+
character(len=*)                             , intent(out) :: l_infile         ! Input file
integer                                      , intent(out) :: l_layer          ! Current layer in data set
LOGICAL                                      , intent(out) :: l_direct         ! Direct space == TRUE
integer                                      , intent(out) :: l_ndims          ! Number of dimensions
integer           , dimension(3)             , intent(out) :: l_dims           ! Actual dimensions
logical                                      , intent(out) :: l_is_grid        ! Data on periodic grid
logical                                      , intent(out) :: l_has_dxyz       ! Data on periodic grid
logical                                      , intent(out) :: l_has_dval       ! Data on periodic grid
real(kind=PREC_DP), dimension(3,4)           , intent(out) :: l_corners        ! 
real(kind=PREC_DP), dimension(3,3)           , intent(out) :: l_vectors        ! 
real(kind=PREC_DP), dimension(3)             , intent(out) :: l_cr_a0          ! 
real(kind=PREC_DP), dimension(3)             , intent(out) :: l_cr_win         ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: l_x              ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: l_y              ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: l_z              ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: l_dx             ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: l_dy             ! 
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: l_dz             ! 
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(out):: l_data           ! Actual diffraction data
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(out):: l_sigma          ! Actual diffraction data
real(kind=PREC_DP), dimension(3)             , intent(out) :: l_llims          ! Lower limits
real(kind=PREC_DP), dimension(3)             , intent(out) :: l_steps          ! steps in H, K, L
real(kind=PREC_DP), dimension(3,3)           , intent(out) :: l_steps_full     ! steps in H, K, L
real(kind=PREC_DP), dimension(2)             , intent(out) :: l_minmaxval        ! steps in H, K, L
real(kind=PREC_DP), dimension(3,2)           , intent(out) :: l_minmaxcoor       ! steps in H, K, L
!
l_infile   = h5_temp%infile         ! input file
l_layer    = h5_temp%layer          ! Current layer in data set
l_direct   = h5_temp%is_direct         ! Direct space == TRUE
l_ndims    = h5_temp%ndims          ! Number of dimensions
l_dims     = h5_temp%dims           ! Actual dimensions
l_is_grid  = h5_temp%is_grid        ! Actual dimensions
l_has_dxyz = h5_temp%has_dxyz       ! Actual dimensions
l_has_dval = h5_temp%has_dval       ! Actual dimensions
l_corners  = h5_temp%corners        ! Actual dimensions
l_vectors  = h5_temp%vectors        ! Actual dimensions
l_cr_a0    = h5_temp%cr_a0          ! Lattice parameters
l_cr_win   = h5_temp%cr_win         ! lattice angles
!allocate(l_data(h5_temp%dims(1), h5_temp%dims(2), h5_temp%dims(3)))
!if(.not. h5_temp%is_grid) then
allocate(l_x(1:h5_temp%dims(1)))
allocate(l_y(1:h5_temp%dims(2)))
allocate(l_z(1:h5_temp%dims(3)))
   l_x        = h5_temp%x              ! Actual x-coordinates
   l_y        = h5_temp%y              ! Actual y-coordinates
   l_z        = h5_temp%z              ! Actual z-coordinates
!endif
!if(h5_temp%has_dxyz) then
   l_dx        = h5_temp%dx              ! Actual x-coordinates uncertainties
   l_dy        = h5_temp%dy              ! Actual y-coordinates uncertainties
   l_dz        = h5_temp%dz              ! Actual z-coordinates uncertainties
!endif
allocate(l_data(h5_temp%dims(1), h5_temp%dims(2), h5_temp%dims(3)))
l_data     = h5_temp%datamap        ! Actual diffraction data
if(h5_temp%has_dval) then
   l_sigma    = h5_temp%sigma          ! Actual diffraction data
endif
l_llims        = h5_temp%llims          ! Lower limits
l_steps        = h5_temp%steps          ! steps in H, K, L
l_steps_full   = h5_temp%steps_full          ! steps in H, K, L
l_minmaxval(1) = minval(h5_temp%datamap)
l_minmaxval(2) = maxval(h5_temp%datamap)
l_minmaxcoor(1,1) = minval(h5_temp%x)
l_minmaxcoor(1,2) = maxval(h5_temp%x)
l_minmaxcoor(2,1) = minval(h5_temp%y)
l_minmaxcoor(2,2) = maxval(h5_temp%y)
l_minmaxcoor(3,1) = minval(h5_temp%z)
l_minmaxcoor(3,2) = maxval(h5_temp%z)
!write(*,*) ' get_node   ' , allocated(l_x), allocated(h5_temp%dx), h5_temp%data_num
!write(*,*) ' X in NODE ', minval(h5_temp%x), maxval(h5_temp%x)
!write(*,*) ' X in NODE ', minval(l_x), maxval(l_x)
!
end subroutine dgl5_get_node
!
!*******************************************************************************
!
subroutine dgl5_set_pointer(izz, ier_num, ier_typ, node_number)
!-
!  Find the node associated to kuplot data set number izz
!+
integer, intent(in ) :: izz
integer, intent(out) :: ier_num
integer, intent(out) :: ier_typ
integer, intent(out) :: node_number
!
!write(*,*) ' ROOT ST', associated(h5_root), ASSOCIATED(h5_temp), node_number, h5_number
if(.NOT. associated(h5_root)) then
   ier_num = -74         ! Root node does not exist !
   ier_typ =   6         ! ER_APPL
   return
endif
!
h5_temp => h5_root
!write(*,*) ' DATA NUMB ', h5_root%data_num, h5_temp%data_num 
find_node: do            ! Search for node
!write(*,*) ' SEARCHING ', izz, h5_ku_is_h5(izz), h5_temp%data_num, h5_ku_is_h5(izz) == h5_temp%data_num
   if(h5_ku_is_h5(izz) == h5_temp%data_num) then
      node_number = h5_temp%data_num
      exit find_node
   else
      if(associated(h5_temp%after)) then    ! A next node exists 
         h5_temp => h5_temp%after
      else
         ier_num = -74         ! Root node does not exist !
         ier_typ =   6         ! ER_APPL
         return
      endif
   endif
enddo find_node
!write(*,*) ' NODE number ', izz,h5_temp%data_num, h5_temp%infile(1:40)
!write(*,*) ' NODE minmax ', h5_temp%minmaxval
!
end subroutine dgl5_set_pointer
!
!*******************************************************************************
!
subroutine dgl5_find_node(node_number, ier_num, ier_typ)
!-
!  Find the node with node_number
!+
integer, intent(in)  :: node_number
integer, intent(out) :: ier_num
integer, intent(out) :: ier_typ
!
if(.NOT. associated(h5_root)) then
   ier_num = -74         ! Root node does not exist !
   ier_typ =   6         ! ER_APPL
   return
endif
h5_find => h5_root
find_node: do            ! Search for node
!write(*,*) 'dgl5_find_node ', node_number, h5_find%data_num
   if(node_number == h5_find%data_num) then
      exit find_node
   else
      if(associated(h5_find%after)) then    ! A next node exists 
         h5_find => h5_find%after
      else
         ier_num = -74         ! Root node does not exist !
         ier_typ =   6         ! ER_APPL
         return
      endif
   endif
enddo find_node
!write(*,*) 'dgl5_FIND_node ', node_number, h5_find%data_num
!
end subroutine dgl5_find_node
!
!*******************************************************************************
!
integer function dgl5_get_layer()
!
implicit none
!
dgl5_get_layer = h5_temp%layer
!
end function dgl5_get_layer
!
!*******************************************************************************
!
integer function dgl5_get_number()
!
implicit none
!
dgl5_get_number = h5_number
!
end function dgl5_get_number
!
!*******************************************************************************
!
integer function dgl5_get_ndims()
!
implicit none
!
dgl5_get_ndims = h5_temp%ndims
!
end function dgl5_get_ndims
!
!*******************************************************************************
!
integer function dgl5_get_ku_is_h5(izz)
!
implicit none
!
integer,intent(in) :: izz
!
dgl5_get_ku_is_h5 = h5_ku_is_h5(izz)
!
end function dgl5_get_ku_is_h5
!
!*******************************************************************************
!
integer function dgl5_get_h5_is_ku(inumber)
!
implicit none
!
integer,intent(in) :: inumber
!
dgl5_get_h5_is_ku = h5_h5_is_ku(inumber)
!
end function dgl5_get_h5_is_ku
!
!*******************************************************************************
!
subroutine dgl5_set_layer(h5_layer)
!
implicit none
!
integer, intent(in) :: h5_layer
!
h5_temp%layer = h5_layer
!
end subroutine dgl5_set_layer
!
!*******************************************************************************
!
LOGICAL function dgl5_get_direct()
!
implicit none
!
dgl5_get_direct = h5_temp%is_direct
!
end function dgl5_get_direct
!
!*******************************************************************************
!
LOGICAL function dgl5_get_is_grid()
!
implicit none
!
dgl5_get_is_grid = h5_temp%is_grid
!
end function dgl5_get_is_grid
!
!*******************************************************************************
!
real function dgl5_get_height()
!
implicit none
!
dgl5_get_height = h5_temp%llims(3) + (h5_temp%layer-1)*h5_temp%steps(3)
!
end function dgl5_get_height
!
!*******************************************************************************
!
subroutine dgl5_get_dims(idata, dims)
!
implicit none
!
integer,               intent(in)  :: idata
integer, dimension(3), intent(out) :: dims
!
dims = h5_temp%dims
!
end subroutine dgl5_get_dims
!
!*******************************************************************************
!
subroutine dgl5_get_lattice(idata, lattice)
!
implicit none
!
integer,                          intent(in)  :: idata
real(kind=PREC_DP), dimension(6), intent(out) :: lattice
!
lattice(1:3) = h5_temp%cr_a0
lattice(4:6) = h5_temp%cr_win
!
end subroutine dgl5_get_lattice
!
!*******************************************************************************
!
subroutine dgl5_get_llims(idata, llims)
!
use precision_mod
!
implicit none
!
integer,               intent(in)  :: idata
real(kind=PREC_DP), dimension(3), intent(out) :: llims
!
llims = h5_temp%llims
!
end subroutine dgl5_get_llims
!
!*******************************************************************************
!
subroutine dgl5_get_minmax(idata, minmax)
!
use precision_mod
!
implicit none
!
integer,               intent(in)  :: idata
real(kind=PREC_DP), dimension(2), intent(out) :: minmax
!
minmax = h5_temp%minmaxval
!
end subroutine dgl5_get_minmax
!
!*******************************************************************************
!
subroutine dgl5_get_minmaxcoor(idata, minmax)
!
use precision_mod
!
implicit none
!
integer,               intent(in)  :: idata
real(kind=PREC_DP), dimension(3,2), intent(out) :: minmax
!
minmax = h5_temp%minmaxcoor
!
end subroutine dgl5_get_minmaxcoor
!
!*******************************************************************************
!
subroutine dgl5_get_steps(idata, steps)
!
!use hdf5_def_mod
use precision_mod
!
implicit none
!
integer,               intent(in)  :: idata
real(kind=PREC_DP), dimension(3,3), intent(out) :: steps
!
if(ALL(yd_present(YD_step_sizes_abs:YD_step_sizes_TOP))) then
   steps = h5_temp%steps_full
else
   steps = 0.0D0
   steps(1,1) = h5_temp%steps(1)
   steps(2,2) = h5_temp%steps(2)
   steps(3,3) = h5_temp%steps(3)
endif
!
end subroutine dgl5_get_steps
!
!*******************************************************************************
!
subroutine dgl5_get_x(dims, x)
!
implicit none
!
integer,            dimension(3),       intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1)), intent(out) :: x
!
x = h5_temp%x
!
end subroutine dgl5_get_x
!
!*******************************************************************************
!
subroutine dgl5_get_y(dims, y)
!
implicit none
!
integer,            dimension(3),       intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1)), intent(out) :: y
!
y = h5_temp%y
!
end subroutine dgl5_get_y
!
!*******************************************************************************
!
subroutine dgl5_get_z(dims, z)
!
implicit none
!
integer,            dimension(3),       intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1)), intent(out) :: z
!
z = h5_temp%z
!
end subroutine dgl5_get_z
!
!*******************************************************************************
!
subroutine dgl5_get_map(dims, odata)
!
implicit none
!
integer,            dimension(3),                         intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1), dims(2), dims(3)), intent(out) :: odata
!
integer :: i,j,k
!
do i=1, dims(1)
   do j=1, dims(2)
      do k=1, dims(3)
         odata(i,j,k) = h5_temp%datamap(i,j,k)
      enddo
   enddo
enddo
!
end subroutine dgl5_get_map
!
!*******************************************************************************
!
subroutine dgl5_get_tmap(dims, odata)
!-
! Get the matrix in transposed form, which is the regular Fortran style
!+
!
implicit none
!
integer,            dimension(3),                         intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1), dims(2), dims(3)), intent(out) :: odata
!
integer :: i,j,k
!
do i=1, dims(1)
   do j=1, dims(2)
      do k=1, dims(3)
         odata(i,j,k) = h5_temp%datamap(k,j,i)
      enddo
   enddo
enddo
!
end subroutine dgl5_get_tmap
!
!*******************************************************************************
!
real(kind=PREC_DP) function dgl5_get_data(i,j,k)
!
!  Get a single data point from the data structure
!
integer, intent(in) :: i
integer, intent(in) :: j
integer, intent(in) :: k
!
dgl5_get_data = h5_temp%datamap(i,j,k)
!
end function dgl5_get_data
!
!*******************************************************************************
!
subroutine dgl5_get_infile(idata, infile)
!
! Get the file name for the current data set
!
integer         , intent(in)  :: idata
character(len=*), intent(out) :: infile
!
infile = h5_temp%infile
!
end subroutine dgl5_get_infile
!
!*******************************************************************************
!
subroutine dgl5_set_is_grid(is_grid)
!
implicit none
!
logical, intent(in) :: is_grid
!
h5_temp%is_grid = is_grid
!
end subroutine dgl5_set_is_grid
!
!*******************************************************************************
!
subroutine dgl5_set_ku_is_h5(izz,ku_is_h5)
!
implicit none
!
integer, intent(in) :: izz
integer, intent(in) :: ku_is_h5
!
h5_ku_is_h5(izz) = ku_is_h5
!
end subroutine dgl5_set_ku_is_h5
!
!*******************************************************************************
!
subroutine dgl5_set_h5_is_ku(inumber, h5_is_ku)
!
implicit none
!
integer, intent(in) :: inumber
integer, intent(in) :: h5_is_ku
!
h5_h5_is_ku(inumber) = h5_is_ku
!
end subroutine dgl5_set_h5_is_ku
!
!*******************************************************************************
!
subroutine dgl5_set_llims(idata, llims)
!
use precision_mod
!
implicit none
!
integer,                          intent(in) :: idata
real(kind=PREC_DP), dimension(3), intent(in) :: llims
!
h5_temp%llims = llims
!
end subroutine dgl5_set_llims
!
!*******************************************************************************
!
subroutine dgl5_set_steps(idata, steps_full)
!
!use hdf5_def_mod
use precision_mod
!
implicit none
!
integer,               intent(in)  :: idata
real(kind=PREC_DP), dimension(3,3), intent(in) :: steps_full
!
h5_temp%steps_full = steps_full
h5_temp%steps(1)   = steps_full(1,1)
h5_temp%steps(3)   = steps_full(2,2)
h5_temp%steps(2)   = steps_full(3,3)
!
end subroutine dgl5_set_steps
!
!*******************************************************************************
!
subroutine dgl5_set_x(dims, x)
!
implicit none
!
integer,            dimension(3),       intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1)), intent(in) :: x
!
h5_temp%x = x
h5_temp%minmaxcoor(1,1) = minval(x)
h5_temp%minmaxcoor(1,2) = maxval(x)
!
end subroutine dgl5_set_x
!
!*******************************************************************************
!
subroutine dgl5_set_y(dims, y)
!
implicit none
!
integer,            dimension(3),       intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1)), intent(in) :: y
!
h5_temp%y = y
h5_temp%minmaxcoor(2,1) = minval(y)
h5_temp%minmaxcoor(2,2) = maxval(y)
!
end subroutine dgl5_set_y
!
!*******************************************************************************
!
subroutine dgl5_set_z(dims, z)
!
implicit none
!
integer,            dimension(3),       intent(in)  :: dims
real(kind=PREC_DP), dimension(dims(1)), intent(in) :: z
!
h5_temp%z = z
h5_temp%minmaxcoor(3,1) = minval(z)
h5_temp%minmaxcoor(3,2) = maxval(z)
!
end subroutine dgl5_set_z
!
!*******************************************************************************
!
subroutine dgl5_set_map(dims, odata)
!
implicit none
!
integer,            dimension(3),                         intent(in) :: dims
real(kind=PREC_DP), dimension(dims(1), dims(2), dims(3)), intent(in) :: odata
!
integer :: i,j,k
!
do i=1, dims(1)
   do j=1, dims(2)
      do k=1, dims(3)
         h5_temp%datamap(i,j,k)= odata(i,j,k) 
      enddo
   enddo
enddo
h5_temp%minmaxval(1) = minval(odata)
h5_temp%minmaxval(2) = maxval(odata)
!
end subroutine dgl5_set_map
!
!*******************************************************************************
!
subroutine data2local(ik, ier_num, ier_typ, node_number,                        &
                      infile, layer, is_direct, ndims, dims,                    &
                      is_grid, has_dxyz, has_dval, corners, vectors, a0, win,   &
                      x, y, z, dx, dy, dz,                                      &
                      odata, sigma, llims, steps,                               &
                      steps_full, minmaxval, minmaxcoor)
!-
! Transfer a data set "ik" into the N-Dim local array odata
!+
!
!use lib_data_struc_h5
use precision_mod
!
implicit none
!
integer                                                 , intent(in )   :: ik          ! Kuplot data set number
integer                                                 , intent(out)   :: ier_num     ! Error number
integer                                                 , intent(out)   :: ier_typ     ! error type
integer                                                 , intent(inout) :: node_number ! node in global array
character(len=*)                                        , intent(out)    :: infile      ! File name
integer                                                 , intent(inout) :: layer       ! Current layer to display
logical                                                 , intent(out)    :: is_direct   ! Direct or reciprocal pspace
!integer                                                 , intent(inout) :: nlayer      ! Current layer to display
integer                                                 , intent(out)    :: ndims       ! Number of dimensions
integer                  , dimension(3)                 , intent(out) :: dims        ! Actual x coordinates
logical                                                 , intent(out)    :: is_grid     ! Data on periodic grid 
logical                                                 , intent(out)    :: has_dxyz    ! Data on periodic grid 
logical                                                 , intent(out)    :: has_dval    ! Data on periodic grid 
real(kind=PREC_DP)       , dimension(3,4)               , intent(out)    :: corners    ! Lattice parameter
real(kind=PREC_DP)       , dimension(3,3)               , intent(out)    :: vectors    ! Lattice parameter
real(kind=PREC_DP)       , dimension(3)                 , intent(out)    :: a0          ! Lattice parameter
real(kind=PREC_DP)       , dimension(3)                 , intent(out)    :: win         ! Lattice angles
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: x           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(out) :: y           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(out) :: z           ! Actual y coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(out) :: dx           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(out) :: dy           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(out) :: dz           ! Actual y coordinates
real(kind=PREC_DP)       , dimension(:,:,:), allocatable, intent(out) :: odata       ! Actual data
real(kind=PREC_DP)       , dimension(:,:,:), allocatable, intent(out) :: sigma       ! Actual data
real(kind=PREC_DP)       , dimension(3)                 , intent(out)    :: llims       ! Lower left coordinates
real(kind=PREC_DP)       , dimension(3  )               , intent(out)    :: steps       ! Increments along axes
real(kind=PREC_DP)       , dimension(3,3)               , intent(out)    :: steps_full  ! Increments along axes
real(kind=PREC_DP)       , dimension(2  )               , intent(out)    :: minmaxval   ! Increments along axes
real(kind=PREC_DP)       , dimension(3,2)               , intent(out)    :: minmaxcoor  ! Increments along axes

call dgl5_set_pointer(ik, ier_num, ier_typ, node_number)
!write(*,*) ' DGL  ik, node : ', ik, node_number, ier_num, ier_typ
layer = dgl5_get_layer()
if(ier_num/=0) return
!
call dgl5_get_node(infile, layer, is_direct, ndims, dims,                    &
                      is_grid, has_dxyz, has_dval, corners, vectors, a0, win,   &
                      x, y, z, dx, dy, dz,                                      &
                      odata, sigma, llims, steps,                               &
                      steps_full, minmaxval, minmaxcoor)
!write(*,*) ' DATA2local ' , allocated(x)
!
end subroutine data2local
!
!*******************************************************************************
!
subroutine local2data(ik, ier_num, ier_typ, node_number,                        &
                      infile, layer, is_direct, ndims, dims,                    &
                      is_grid, has_dxyz, has_dval, corners, vectors, a0, win,   &
                      x, y, z, dx, dy, dz,                                      &
                      odata, sigma, llims, steps,                               &
                      steps_full)
!-
! Transfer a local array into the global storage
! If a node at node_number exists, put there, else create a new node
!+
!
!use lib_data_struc_h5
!use lib_hdf5_params_mod
use precision_mod
!
implicit none
!
integer                                                 , intent(in )   :: ik          ! Kuplot data set number
integer                                                 , intent(out)   :: ier_num     ! Error number
integer                                                 , intent(out)   :: ier_typ     ! error type
integer                                                 , intent(inout) :: node_number ! node in global array
character(len=*)                                        , intent(in)    :: infile      ! File name
integer                                                 , intent(inout) :: layer       ! Current layer to display
logical                                                 , intent(in)    :: is_direct   ! Direct or reciprocal pspace
!integer                                                 , intent(inout) :: nlayer      ! Current layer to display
integer                                                 , intent(in)    :: ndims       ! Number of dimensions
integer                  , dimension(3)                 , intent(in) :: dims        ! Actual x coordinates
logical                                                 , intent(in)    :: is_grid     ! Data on periodic grid 
logical                                                 , intent(in)    :: has_dxyz    ! Data on periodic grid 
logical                                                 , intent(in)    :: has_dval    ! Data on periodic grid 
real(kind=PREC_DP)       , dimension(3,4)               , intent(in)    :: corners    ! Lattice parameter
real(kind=PREC_DP)       , dimension(3,3)               , intent(in)    :: vectors    ! Lattice parameter
real(kind=PREC_DP)       , dimension(3)                 , intent(in)    :: a0          ! Lattice parameter
real(kind=PREC_DP)       , dimension(3)                 , intent(in)    :: win         ! Lattice angles
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: x           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: y           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: z           ! Actual y coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: dx           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: dy           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: dz           ! Actual y coordinates
real(kind=PREC_DP)       , dimension(dims(1),dims(1),dims(1)), intent(in) :: odata       ! Actual data
real(kind=PREC_DP)       , dimension(:,:,:), allocatable, intent(in) :: sigma       ! Actual data
real(kind=PREC_DP)       , dimension(3)                 , intent(in)    :: llims       ! Lower left coordinates
real(kind=PREC_DP)       , dimension(3  )               , intent(in)    :: steps       ! Increments along axes
real(kind=PREC_DP)       , dimension(3,3)               , intent(in)    :: steps_full  ! Increments along axes
!
integer :: i
!write(*,*) ' OLD NODE_NUMBER ', node_number, ik
!
if(node_number>0) then
!write(*,*) ' FIND NODE ', node_number
   call dgl5_find_node(node_number, ier_num, ier_typ)
!write(*,*) ' GOT  NODE ', node_number
   if(ier_num/=0) return
   h5_temp => h5_find
else
   CALL dgl5_new_node
endif
!write(*,*) ' POINTING TO     ',h5_temp%data_num 
node_number = dgl5_get_number()
node_number = h5_temp%data_num
!write(*,*) ' POINTING TO     ',h5_temp%data_num 
!write(*,*) ' NEW NODE_NUMBER ', node_number!, minval(odata), maxval(odata)
if(allocated( x)) deallocate( x)
if(allocated( y)) deallocate( y)
if(allocated( z)) deallocate( z)
if(allocated(dx)) deallocate(dx)
if(allocated(dy)) deallocate(dy)
if(allocated(dz)) deallocate(dz)
allocate( x(1:dims(1)))
allocate( y(1:dims(2)))
allocate( z(1:dims(3)))
allocate(dx(1:dims(1)))
allocate(dy(1:dims(2)))
allocate(dz(1:dims(3)))
do i=1, dims(1)
  x(i) = corners(1,1) + (i-1)*steps_full(1,1)
enddo
do i=1, dims(2)
  y(i) = corners(2,2) + (i-1)*steps_full(2,2)
enddo
do i=1, dims(3)
  z(i) = corners(3,3) + (i-1)*steps_full(3,3)
enddo
dx = 0.0D0
dy = 0.0D0
dz = 0.0D0
!
call dgl5_set_ku_is_h5(ik, node_number)
call dgl5_set_h5_is_ku(node_number, ik)
call dgl5_set_node(infile, layer, is_direct, ndims, dims, &
                   is_grid, has_dxyz, has_dval, corners, vectors,  &
                   a0, win, x, y, z, dx, dy, dz, &
                   odata, sigma, llims, steps,  &
                   steps_full)
!write(*,*) ' KUPLOT DATA SET STORED AS NODE ', dgl5_get_ku_is_h5(ik)
!write(*,*) ' NODE has KUPLOT data set       ', dgl5_get_h5_is_ku(node_number)
!write(*,*) ' ik = 1 : ', dgl5_get_ku_is_h5(1)
!write(*,*) ' ik = 1 : ', dgl5_get_ku_is_h5(2)
!write(*,*) ' ik = 1 : ', dgl5_get_ku_is_h5(3)

!write(*,*) ' node 1 : ', dgl5_get_h5_is_ku(1)
!write(*,*) ' node 2 : ', dgl5_get_h5_is_ku(2)
!write(*,*) ' node 3 : ', dgl5_get_h5_is_ku(3)
deallocate( x)
deallocate( y)
deallocate( z)
deallocate(dx)
deallocate(dy)
deallocate(dz)
!
end subroutine local2data
!
!*******************************************************************************
!
subroutine   fft2data(ik, ier_num, ier_typ, node_number,                        &
                              layer, is_direct, ndims, dims,                    &
                      is_grid, has_dxyz, has_dval, corners, vectors, a0, win,   &
                      x, y, z, dx, dy, dz,                                      &
                      rl_data, rl_sigma, llims, steps,                          &
                      steps_full, im_data, im_sigma)
!-
! Transfer a FFT local array pair into the global storage
!+
!
!use lib_data_struc_h5
!use lib_hdf5_params_mod
use precision_mod
!
implicit none
!
integer                                                 , intent(in )   :: ik          ! Kuplot data set number
integer                                                 , intent(out)   :: ier_num     ! Error number
integer                                                 , intent(out)   :: ier_typ     ! error type
integer                                                 , intent(inout) :: node_number ! node in global array
integer                                                 , intent(inout) :: layer       ! Current layer to display
logical                                                 , intent(inout) :: is_direct   ! Direct or reciprocal pspace
!integer                                                 , intent(inout) :: nlayer      ! Current layer to display
integer                                                 , intent(in)    :: ndims       ! Number of dimensions
integer                  , dimension(3)                 , intent(in)    :: dims        ! Actual x coordinates
logical                                                 , intent(inout) :: is_grid     ! Data on periodic grid 
logical                                                 , intent(inout) :: has_dxyz    ! Data on periodic grid 
logical                                                 , intent(inout) :: has_dval    ! Data on periodic grid 
real(kind=PREC_DP)       , dimension(3,4)               , intent(inout) :: corners     ! Corners ll, lr, ul, tl
real(kind=PREC_DP)       , dimension(3,3)               , intent(inout) :: vectors     ! Increment vectors
real(kind=PREC_DP)       , dimension(3)                 , intent(in)    :: a0          ! Lattice parameter
real(kind=PREC_DP)       , dimension(3)                 , intent(in)    :: win         ! Lattice angles
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: x           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: y           ! Actual x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: z           ! Actual y coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: dx          ! Sigma  x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: dy          ! Sigma  x coordinates
real(kind=PREC_DP)       , dimension(:), allocatable    , intent(inout) :: dz          ! Sigma  y coordinates
real(kind=PREC_DP)       , dimension(dims(1),dims(1),dims(1)), intent(in   ) :: rl_data   ! Actual data real part
real(kind=PREC_DP)       , dimension(:,:,:), allocatable     , intent(in   ) :: rl_sigma  ! Actual sigma real part
real(kind=PREC_DP)       , dimension(3)                 , intent(in)    :: llims       ! Lower left coordinates
real(kind=PREC_DP)       , dimension(3  )               , intent(in)    :: steps       ! Increments along axes
real(kind=PREC_DP)       , dimension(3,3)               , intent(inout) :: steps_full  ! Increments along axes
real(kind=PREC_DP)       , dimension(dims(1),dims(1),dims(1)), intent(in   ) :: im_data   ! Actual data imaginary part
real(kind=PREC_DP)       , dimension(:,:,:), allocatable     , intent(in   ) :: im_sigma  ! Actual sigma imaginary part
!
character(len=PREC_STRING) :: rl_file
character(len=PREC_STRING) :: im_file
integer :: im
!
im = ik + 1
rl_file = 'Fourier_real'
im_file = 'Fourier_imag'
is_direct = .not. is_direct  ! swap direct versus reciprocal space
is_grid   = .true.
has_dxyz  = .false.
has_dval  = .false.
!
steps_full      = 0.0D0
steps_full(1,1) = steps(1)
steps_full(2,2) = steps(2)
steps_full(3,3) = steps(3)
corners(:,1)    = llims
corners(:,2)    = corners(:,1) + (dims(1)-1)*steps_full(:,1)
corners(:,3)    = corners(:,1) + (dims(2)-1)*steps_full(:,2)
corners(:,4)    = corners(:,1) + (dims(3)-1)*steps_full(:,3)
vectors(:,1)    = steps_full(:,1)
vectors(:,2)    = steps_full(:,2)
vectors(:,3)    = steps_full(:,3)
!
node_number = 0
call local2data(ik, ier_num, ier_typ, node_number,                        &
                     rl_file, layer, is_direct, ndims, dims,                    &
                      is_grid, has_dxyz, has_dval, corners, vectors, a0, win,   &
                      x, y, z, dx, dy, dz,                                      &
                      rl_data, rl_sigma, llims, steps,                               &
                      steps_full)
!
node_number = 0
call local2data(im, ier_num, ier_typ, node_number,                        &
                      im_file, layer, is_direct, ndims, dims,                    &
                      is_grid, has_dxyz, has_dval, corners, vectors, a0, win,   &
                      x, y, z, dx, dy, dz,                                      &
                      im_data, im_sigma, llims, steps,                               &
                      steps_full)
!
end subroutine   fft2data
!
!*******************************************************************************
!
subroutine dgl5_reset
!
type(h5_data_struc), pointer :: h5_current => NULL()
!
if(associated(h5_root)) then       ! A storage does exist
   h5_temp => h5_root
   if(allocated(h5_temp%datamap)) deallocate(h5_temp%datamap)
   find_node: do 
      if(associated(h5_temp%after)) then   ! A next node exists
         h5_current => h5_temp             ! Point to current
         h5_temp    => h5_temp%after       ! Point to next node
         deallocate(h5_current)            ! Clean up current node
      else
         h5_current => h5_temp             ! Point to current
         deallocate(h5_current)            ! Clean up current node
         exit find_node                    ! We are done
      endif
   enddo find_node
endif
nullify(h5_temp)
nullify(h5_root)
h5_number   = 0
h5_h5_is_ku = 0
h5_ku_is_h5 = 0
!
end subroutine dgl5_reset
!
!*******************************************************************************
!
end module lib_data_struc_h5
