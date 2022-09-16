module lib_math_mod
!
!  General data manipluation routines that work on the data in the global storage
!
!  calc_h5_coord_global   Calculations on the coordinates and their uncertainties
!  calc_h5_val_global     Calculations on the values/sigmas 
!  edge2zero              Calculates average values at data set edge and subtracts from data
!  filter_minabs          Applies a filter: minimum of (abs(i), abs(j)... )
!  fft_1D_global          Performs a FFT on 2D data
!  fft_2D_global          Performs a FFT on 2D data
!  fft_3D_global          Performs a FFT on 2D data
!  kmath_h5_global        Perform math on two data sets 
!  rvalue_h5_global       R-value calculation
!  do_merge_global        Merge data sets
!  
!*******************************************************************************
!
use precision_mod
!
character(len=PREC_STRING)                         :: ik1_infile
character(len=PREC_STRING)                         :: ik2_infile
integer                                            :: ik1_node_number  ! Node in global data
integer                                            :: ik2_node_number  ! Node in global data
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
real(kind=PREC_DP), dimension(3)                   :: ft_steps        ! Steps after Fourier
real(kind=PREC_DP), dimension(3)                   :: ft_start        ! Start after Fourier
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: res_rl_data        ! result data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: res_im_data        ! result data array for real
real(kind=PREC_DP), dimension(3) :: steps
!
!
private
public calc_h5_coord_global
public calc_h5_val_global
public edge2zero
public filter_minabs
public fft_1D_global
public fft_2D_global
public fft_3D_global
public kmath_h5_global
public rvalue_h5_global
public do_merge_global
!
contains
!
!*******************************************************************************
!
subroutine calc_h5_coord_global(ik, unt, oper, MAXW, werte, ianz) !, node_number, ndims, dims)
!-
!  Manipulate coordinates
!+
!
use errlist_mod
use lib_data_struc_h5
!
implicit none
!
integer                            , intent(in) :: ik     ! Kuplot data set number 
character(len=*)                   , intent(in) :: unt    ! Qualifier 'wx', 'wy', 'wz', 'wi', 'dx', dy'
character(len=*)                   , intent(in) :: oper   ! Calculation command like 'add', 'mul'...
integer                            , intent(in) :: MAXW   ! Dimension of werte
real(kind=PREC_DP), dimension(MAXW), intent(in) :: werte  ! Value do be added, multiplied etc.
integer                            , intent(in) :: ianz   ! Number parameters in werte
!
integer                                   :: jdim     ! Coordinate number to work on
logical                                   :: lcoord   ! Working on coordinates, not on sigmas
real(kind=PREC_DP) :: factor, summand !, thresh
real(kind=PREC_DP), dimension(:), pointer :: p2coord  ! Pointer to Coordinates to work on
!
call data2local(ik      , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
!  Get the appropriate coordinates
!
jdim = 0
   lcoord = .true.
if(unt=='WX') then
   p2coord => ik1_x
   jdim = 1
   lcoord = .true.
elseif(unt=='WY') then
   p2coord => ik1_y
   jdim = 2
   lcoord = .true.
elseif(unt=='WZ') then
   p2coord => ik1_z
   jdim = 3
   lcoord = .true.
elseif(unt=='DX') then
   p2coord => ik1_dx
   jdim = 1
   lcoord = .false.
elseif(unt=='DY') then
   p2coord => ik1_dy
   jdim = 2
   lcoord = .false.
elseif(unt=='DZ') then
   p2coord => ik1_dz
   jdim = 3
   lcoord = .false.
endif
!
! Perform actual modification
!
if(oper=='ADD') then
   summand = 0.0D0
   if(ianz==2) summand = werte(2)
   p2coord     = p2coord       + summand
elseif(oper.eq.'EXP') then
   p2coord       = exp(p2coord)
   if(lcoord) then
      ik1_is_grid = .false.                 ! No longer an equidistant grid
   endif
elseif(oper.eq.'INV') then
   where(p2coord/=0.0D0)
      p2coord = 1.0D0/p2coord
   end where
   if(lcoord) then
      ik1_is_grid = .false.                 ! No longer an equidistant grid
   endif
elseif(oper.eq.'LOG') then
   where(p2coord>=0.0D0)
      p2coord = log(p2coord)
   end where
   if(lcoord) then
      ik1_is_grid = .false.                 ! No longer an equidistant grid
   endif
elseif(oper=='MUL') then
   factor = 1.0D0
   if(ianz==2) factor = werte(2)
   p2coord = p2coord * factor
   if(lcoord) then
      ik1_llims(jdim)           = ik1_llims(jdim)           * factor
      ik1_steps(jdim)           = ik1_steps(jdim)           * factor
      ik1_steps_full(jdim,jdim) = ik1_steps_full(jdim,jdim) * factor
   endif
elseif(oper.eq.'SQR') then
   where(p2coord>=0.0D0)
      p2coord = sqrt(p2coord)
   end where
   if(lcoord) then
      ik1_is_grid = .false.                 ! No longer an equidistant grid
   endif
elseif(oper.eq.'SQU') then
   p2coord = p2coord **2
   if(lcoord) then
      ik1_is_grid = .false.                 ! No longer an equidistant grid
   endif
elseif(oper.eq.'ABS') then
   p2coord       = abs(p2coord)
   if(lcoord) then
      ik1_is_grid = .false.                 ! No longer an equidistant grid
   endif
else
   ier_num = -6
   ier_typ = ER_COMM
endif
!
if(lcoord) then
   ik1_llims(jdim) = minval(p2coord)
endif
!
call local2data(ik, ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_nlayer,  &
     ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
     ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
     ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma, ik1_llims, ik1_steps,  &
     ik1_steps_full)
!
nullify(p2coord)
call math_clean
!
end subroutine calc_h5_coord_global
!
!*****7*****************************************************************
!
subroutine calc_h5_val_global(ik, unt, oper, MAXW, werte, ianz)
!+                                                                      
!     Calculations for a (3D) data set stored in the (HDF5)-type general storage
!-
!
use errlist_mod
use lib_data_struc_h5
use lib_random_func
use precision_mod
!
implicit none
!
integer                            , intent(in) :: ik     ! Kuplot data set number 
character(len=*)                   , intent(in) :: unt    ! Qualifier 'val', 'sig'
character(len=*)                   , intent(in) :: oper   ! Calculation command like 'add', 'mul'...
integer                            , intent(in) :: MAXW   ! Dimension of werte
real(kind=PREC_DP), dimension(MAXW), intent(in) :: werte  ! Value do be added, multiplied etc.
integer                            , intent(in) :: ianz   ! Number parameters in werte
!
integer :: i,j,k
real(kind=PREC_DP) :: factor, summand !, thresh
real(kind=PREC_DP), dimension(:,:,:), pointer :: p2data  ! Pointer to data to work on
!
call data2local(ik      , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
if(unt=='VAL') then
   p2data => ik1_data
elseif(unt=='SIG') then
   p2data => ik1_sigma
endif
!                                                                       
if(oper=='ADD') then
   summand = 0.0D0
   if(ianz==2) summand = werte(2)
   p2data = p2data + summand
elseif(oper.eq.'EXP') then
   p2data = exp(p2data)
elseif(oper.eq.'INV') then
   where(p2data/=0.0D0)
      p2data = 1.0D0/p2data
   end where
elseif(oper.eq.'LOG') then
   where(p2data>=0.0D0)
      p2data = log(p2data)
   end where
elseif(oper=='MUL') then
   factor = 1.0D0
   if(ianz==2) factor = werte(2)
   p2data = p2data * factor
elseif(oper.eq.'SQR') then
   where(p2data>=0.0D0)
      p2data = sqrt(p2data)
   end where
elseif(oper.eq.'SQU') then
   p2data = p2data **2
elseif(oper.eq.'SIG') then
   if(allocated(ik1_sigma)) then
     do k=1, ik1_dims(3)
     do j=1, ik1_dims(2)
     do i=1, ik1_dims(1)
     ik1_data(i,j,k) = ik1_data(i,j,k) + gasdev(ik1_sigma(i,j,k))
     enddo
     enddo
     enddo
   else
      ier_num = -6
      ier_typ = ER_COMM
   endif
else
   ier_num = -6
   ier_typ = ER_COMM
endif
!
if(ier_num==0) then
call local2data(ik, ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_nlayer,  &
     ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
     ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
     ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma, ik1_llims, ik1_steps,  &
     ik1_steps_full)
endif
!
nullify(p2data)
call math_clean
!
end subroutine calc_h5_val_global
!
!*******************************************************************************
!
subroutine edge2zero(line, length)
!-
!  Get the average value around the edge and subtract this from all values
!+
!
use ber_params_mod
use errlist_mod
use get_params_mod
use lib_data_struc_h5
use precision_mod
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: line     ! Input command line
integer         , intent(inout) :: length   ! Input command line length
!
integer, parameter :: EDGE   = -1
integer, parameter :: SPHERE = -2
!
integer :: ik     ! Data set number
!
integer, parameter :: MAXW = 3
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
integer :: ianz           ! Number of parameters
real(kind=PREC_DP)        , dimension(MAXW) :: werte    ! Calculated values
real(kind=PREC_DP)        , dimension(3   ) :: vector   ! Calculated values
integer                   , dimension(3   ) :: inull    ! Calculated values
integer                                     :: radius   ! Calculated values
integer :: i,j,k
integer :: ii,jj,kk, ll
!
logical           , dimension(:,:,:), allocatable  :: ed_mask      ! Mask field edges
logical           , dimension(:,:,:), allocatable  :: ze_mask      ! Mask field zeros
!
integer             :: nval
integer             :: mask
real(kind=PREC_DP)  :: av_val
!
integer, parameter :: NOPTIONAL = 3
integer, parameter :: O_DATA    = 1
integer, parameter :: O_MASK    = 2
integer, parameter :: O_RADIUS  = 3
character(len=          6), dimension(NOPTIONAL) :: oname     ! Optional parameter names
character(len=PREC_STRING), dimension(NOPTIONAL) :: opara     ! Optional parameter strings returned
integer                   , dimension(NOPTIONAL) :: loname    ! Lenght opt. para name
integer                   , dimension(NOPTIONAL) :: lopara    ! Lenght opt. para name returned
logical                   , dimension(NOPTIONAL) :: lpresent  ! opt. para present
real(kind=PREC_DP)        , dimension(NOPTIONAL) :: owerte    ! Calculated values
integer                   , parameter            :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'data ', 'mask', 'radius'/
data loname /  4     ,  4    ,  6/
opara  =  (/ '1.000000000000', 'edge          ', '[-1.0,0.0,0.0]' /)
lopara =  (/  14             ,  4              ,  14              /)
owerte =  (/  1.0            ,  -1.0           ,  0.0             /)
!
call get_params (line, ianz, cpara, lpara, MAXW, length)
if(ier_num /= 0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num /= 0) return
if(.not.lpresent(O_DATA)) then
   if(ianz==1) then
      call ber_params(ianz, cpara, lpara, werte, MAXW)
      if(ier_num /= 0) return
      owerte(O_DATA) = werte(1)
   else
      ier_num = -6
      ier_typ = ER_COMM
      return
   endif
endif
!
ik = nint(owerte(O_DATA))
!
mask = EDGE
if(lpresent(O_MASK)) then
   if(opara(O_MASK)=='edge') then
      mask = EDGE
   elseif(opara(O_MASK)=='sphere') then
      mask = SPHERE
      if(lpresent(O_RADIUS)) then
         i = index(opara(O_RADIUS),':')
         cpara(1) = opara(O_RADIUS)(i+1:)
         lpara(1) = len_trim(cpara(1))
         call get_optional_multi(MAXW, cpara(1), lpara(1), vector, ianz)
      else
         ier_num = -58
         ier_num = ER_FORT
         return
      endif
   else
      ianz = 1
      cpara(1) = opara(O_MASK)
      lpara(1) = lopara(O_MASK)
      call ber_params(ianz, cpara, lpara, werte, MAXW)
      if(ier_num /= 0) return
      mask = nint(werte(1))
   endif
!  if(opara(O_MASK)=='edge') then
!     mask = EDGE
!  elseif(opara(O_MASK)=='sphere') then
!     mask = SPHERE
!  else
!  endif
else
   mask = EDGE
endif
!
! Get data field
!
call data2local(ik      , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
! Get/ Prepare mask field
!
allocate(ed_mask(ik1_dims(1), ik1_dims(2), ik1_dims(3)))
allocate(ze_mask(ik1_dims(1), ik1_dims(2), ik1_dims(3)))
ed_mask = .false.
ze_mask = .false.
nval    = 1
!
if(mask>0) then
!
!  Get mask field
!
   call data2local(mask    , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
        ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
        ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
        ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
        ik2_llims, ik2_steps, ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
   if(ik1_ndims/=ik2_ndims .or. any(ik1_dims/=ik2_dims)) then
      ier_num = -56
      ier_typ = ER_FORT
   endif
   where(nint(ik2_data)==-1)
     ed_mask = .true.
   end where
   nval = count(ed_mask)
   where(nint(ik2_data)==0 )
     ze_mask = .true.
   end where
!
elseif(mask==EDGE) then
   if(ik1_ndims==1) then
      ed_mask(1          ,1,1) = .true.
      ed_mask(ik1_dims(1),1,1) = .true.
      nval = 2
   elseif(ik1_ndims==2) then
      ed_mask(:          ,1          ,1) = .true.
      ed_mask(:          ,ik1_dims(2),1) = .true.
      ed_mask(1          ,:          ,1) = .true.
      ed_mask(ik1_dims(1),:          ,1) = .true.
      nval = 2*(ik1_dims(1)+ik1_dims(2)) - 4
   elseif(ik1_ndims==3) then
      ed_mask(:          ,:          ,1          ) = .true.
      ed_mask(:          ,:          ,ik1_dims(3)) = .true.
      ed_mask(:          ,1          ,:          ) = .true.
      ed_mask(:          ,ik1_dims(2),:          ) = .true.
      ed_mask(1          ,:          ,:          ) = .true.
      ed_mask(ik1_dims(1),:          ,:          ) = .true.
      nval = 2*(ik1_dims(1)*(ik1_dims(2)-1) + (ik1_dims(1)-1)*ik1_dims(3) +             &
                ik1_dims(2)*(ik1_dims(3)-1)                             ) - 8
   endif
elseif(mask==SPHERE) then
   radius   = (int(vector(1)/ik1_steps_full(1,1)))
   inull(1) = nint(-ik1_llims(1)/ik1_steps_full(1,1))+1
   inull(2) = nint(-ik1_llims(2)/ik1_steps_full(2,2))+1
   inull(3) = nint(-ik1_llims(3)/ik1_steps_full(3,3))+1
   do k=1, ik1_dims(3)
      kk = (k-inull(3))**2
      do j=1, ik1_dims(2)
         jj = (j-inull(2))**2
         do i=1, ik1_dims(1)
            ii = (i-inull(1))**2
            ll = radius-nint(sqrt(1.0D0*(ii + jj + kk)))
            if(ll>=0) then 
               if(ll< 3) then
                  ed_mask(i,j,k) = .true.
               endif
            else
               ze_mask(i,j,k) = .true.
            endif
         enddo
      enddo
   enddo
   nval = count(ed_mask)
endif
!
cond_error: if(ier_num==0) then
   av_val = 0.0D0
   av_val = sum(ik1_data, ed_mask) / real(nval,kind=PREC_DP)
   ik1_data = ik1_data - av_val
! Zero outside range
   where(ze_mask)
      ik1_data = 0.0D0
   end where
!
   call local2data(ik , ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_nlayer,  &
        ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
        ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
        ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma, ik1_llims, ik1_steps,  &
        ik1_steps_full)
endif cond_error
!
deallocate(ze_mask)
deallocate(ed_mask)
!
call math_clean
!
end subroutine edge2zero
!
!*******************************************************************************
!
subroutine filter_minabs(line, length)
!-
!  Get the average value around the edge and subtract this from all values
!+
!
use errlist_mod
use get_params_mod
use lib_data_struc_h5
use precision_mod
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: line     ! Input command line
integer         , intent(inout) :: length   ! Input command line length
!
!
integer, parameter :: MAXW = 4
integer, parameter :: MAXU = 20
!
character(len=PREC_STRING)                  :: string   ! 
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
integer :: ianz           ! Number of parameters
integer :: janz           ! Number of parameters
integer :: l
integer :: ik             ! Data set number
!
integer                                            :: ndims        ! Number of dimensions
integer                                            :: f_ndims        ! Number of dimensions
integer                                            :: node_number  ! Number of dimensions
integer, dimension(3)                              :: dims         ! Dimensions global array
integer, dimension(3)                              :: f_dims         ! Dimensions global array
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: fil_data        ! Global data array for result
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: res_data        ! Global data array for result
!
real(kind=PREC_DP), dimension(MAXU) :: uwerte   ! Calculated values
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_DATA    = 1
integer, parameter :: O_FILT    = 2
character(len=          6), dimension(NOPTIONAL) :: oname     ! Optional parameter names
character(len=PREC_STRING), dimension(NOPTIONAL) :: opara     ! Optional parameter strings returned
integer                   , dimension(NOPTIONAL) :: loname    ! Lenght opt. para name
integer                   , dimension(NOPTIONAL) :: lopara    ! Lenght opt. para name returned
logical                   , dimension(NOPTIONAL) :: lpresent  ! opt. para present
real(kind=PREC_DP)        , dimension(NOPTIONAL) :: owerte    ! Calculated values
integer                   , parameter            :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'data ' , 'filter'/
data loname /  4      ,  6      /
opara  =  (/ '1.0000' ,'[1]   '/)
lopara =  (/  6       , 3      /)
owerte =  (/  1.0     , 1.0    /)
!
call get_params (line, ianz, cpara, lpara, MAXW, length)
if(ier_num /= 0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num /= 0) return
!
string = opara(O_FILT)
length = lopara(O_FILT)
call get_optional_multi(MAXU, string, length, uwerte, janz)
if(ier_num/=0) return
!
ik = nint(owerte(O_DATA))
!
!  Get data to check dimensions and to initialize the result array
!
call dgl5_set_pointer(ik, ier_num, ier_typ, node_number)
ndims = dgl5_get_ndims()
call dgl5_get_dims(node_number, dims)
allocate(res_data(dims(1), dims(2), dims(3)))
allocate(fil_data(dims(1), dims(2), dims(3)))
res_data = HUGE(1.0)
!
loop_filter:do l=1, janz              ! Loop over all filters
   ik = nint(uwerte(l))
   call dgl5_set_pointer(ik, ier_num, ier_typ, node_number)
   if(ier_num/=0) exit loop_filter
!
   f_ndims =  dgl5_get_ndims()
   call dgl5_get_dims(node_number, f_dims)
   if(ndims/=f_ndims .or. ALL(dims/=f_dims)) then
      ier_num = -56
      ier_typ = ER_FORT
      exit loop_filter
   endif
   call dgl5_get_map(dims, fil_data)
   res_data = min(res_data, abs(fil_data))    ! Apply min(abs()) Filter
enddo loop_filter
!
if(ier_num==0) then
   ik = nint(owerte(O_DATA))
   call dgl5_set_pointer(ik, ier_num, ier_typ, node_number)
   call dgl5_get_map(dims, fil_data)   ! Get original data into fil_data
   where(fil_data<0.0d0)
      fil_data = -1.0D0*res_data              ! Negative data
   else where
      fil_data =        res_data              ! Positive data
   end where
   call dgl5_set_map(dims, fil_data)
endif
!
deallocate(res_data)
deallocate(fil_data)
call math_clean
!
end subroutine filter_minabs
!
!*******************************************************************************
!
subroutine fft_1D_global(idata, odata, xscale, yscale, ipad)
!-
!  Calculate 1D FFT of data set idata
!+
!
use lib_data_struc_h5
use errlist_mod
use map_1dtofield
use singleton
use wink_mod
!
implicit none
!
integer, dimension(2), INTENT(IN)  :: idata     ! Data set number to transform
integer, dimension(2), INTENT(OUT) :: odata     ! Data set number to transform
real(kind=PREC_DP)   , intent(in)  :: xscale
real(kind=PREC_DP)   , intent(in)  :: yscale
integer, dimension(3), intent(in)  :: ipad      ! Pad to these dimensions
!
integer :: i
integer :: kdat
integer :: length             ! Data set length
integer :: istart             ! first index in case of padding
integer, dimension(3) :: num  ! DATA set dimensions
integer, dimension(3) :: dsort
complex(kind=kind(0.0D0)) , dimension(:), allocatable  :: k_data   ! The Kuplot in data set)
complex(kind=kind(0.0D0)) , dimension(:), allocatable  :: pattern  ! The Curve to be FFT'd
logical, dimension(3) :: leven            ! Make original pattern symmetric around 0,0,0
integer, dimension(3) :: peven            ! Make original pattern symmetric around 0,0,0
!
!
if(idata(1)>0) then                             !  User provided real par
   call data2local(idata(1), ier_num, ier_typ, ik1_node_number, ik1_infile,     &
        ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
        ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
        ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
        ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
   length = ik1_dims(1)
   kdat   = 1
   if(idata(2)==0) then
      ik2_ndims = ik1_ndims
      ik2_dims  = ik1_dims
      allocate(ik2_data(ik2_dims(1), ik2_dims(2),ik2_dims(3)))
      ik2_data = 0.0D0
   endif
endif
if(idata(2)>0) then                             !  User provided imag part
   call data2local(idata(2), ier_num, ier_typ, ik2_node_number, ik2_infile,     &
        ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
        ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
        ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
        ik2_llims, ik2_steps, ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
   length = ik2_dims(1)
   kdat   = 2
   if(idata(1)==0) then
      ik1_ndims = ik2_ndims
      ik1_dims  = ik2_dims
      allocate(ik1_data(ik1_dims(1), ik1_dims(2),ik1_dims(3)))
      ik1_data = 0.0D0
   else
      if(ik1_ndims/=ik2_ndims .or. any(ik1_dims/=ik2_dims)) then
         ier_num = -56
         ier_typ = ER_FORT
         call math_clean
         return
      endif
   endif
endif
!
leven    = .false.
peven    = 0
if(mod(ik1_dims(1),2)==0) then      ! Even number of pixels
  leven(1) = .true.
  peven(1) = 1
endif
!
num    = 1
num(1) = max(ipad(1), ik1_dims(1)         )       ! Dimensions are increased to padding if needed
istart = nint((num(1)-ik1_dims(1)+0.1D0)/2.0D0)   ! Calculate start index prior to symmetry padding
num(1) = max(ipad(1), ik1_dims(1)+peven(1))       ! Dimensions are increased for symmetry padding if needed
dsort(1) = 1                       ! 1D does not need sort
dsort(2) = 2
dsort(3) = 3
allocate(k_data (num(1)))          ! Complex array for FFT
allocate(pattern(num(1)))          ! Complex array for FFT in shifted sequence           
!
k_data = cmplx(0.0D0, 0.0D0)
do i=1, ik1_dims(1)
   k_data(istart+i) = CMPLX(ik1_data(i,1,1), ik2_data(i,1,1))
enddo
if(leven(1)) then 
   k_data(istart+ik1_dims(1)+1) = k_data(istart+1)
endif
!
call maptofftfd(num, dsort, k_data, pattern)           ! Shifts origin to (1,1,1)
!
pattern = fft(pattern) / SQRT(REAL(num(1))) * yscale   ! Actual fft routine
!
call mapfftfdtoline(num, dsort, k_data, pattern)       ! Shift backinto discus sequence
!
!  Convert limits (currently the assumption is that the origininal data are -x...0...+x
!  
ft_start(1)   = -0.5D0/ik1_steps(1)  
ft_start(2:3) =  0.00D0
ft_steps(1)   =  1.0D0/(ik1_steps(1)  *real(num(1)  -1,kind=PREC_DP))
ft_steps(2:3) =  0.0D0
!
allocate(res_rl_data(num(1),      1, 1))
allocate(res_im_data(num(1),      1, 1))
res_rl_data(:,1,1) =  real(k_data(:),kind=PREC_DP)
res_im_data(:,1,1) = dimag(k_data(:))
!
ik1_dims  = num          ! Copy current dimensions 
steps(1) = ft_steps(1)       ! New steps 
steps(2) = 0.0D0
steps(3) = 0.0D0
ik1_llims(1) = ft_start(1)    ! New lower limits
ik1_llims(2) = 0.0D0
ik1_llims(3) = 0.0D0
ik1_nlayer      = 1   ! Dummy layer number
!
! Store real/imag pair into global storage
!
call fft2data(odata(1), ier_num, ier_typ, ik1_node_number, ik1_nlayer,          &
     ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
     ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
     ik1_z, ik1_dx, ik1_dy, ik1_dz, res_rl_data, ik1_sigma, ik1_llims, steps,   &
     ik1_steps_full, res_im_data, ik2_sigma)
!
deallocate(k_data )
deallocate(pattern)
deallocate(res_rl_data)
deallocate(res_im_data)
call math_clean
!
end subroutine fft_1D_global
!
!*****7**************************************************************** 
!
subroutine fft_2D_global(idata, odata, xscale, yscale, ipad)
!-
!  Calculate 2D FFT of data sets idata; 
!  If ipad > dimension of data pad with zeros
!+
!
use lib_data_struc_h5
use errlist_mod
use map_1dtofield
use singleton
use wink_mod
!
implicit none
!
integer, dimension(2), INTENT(IN)  :: idata     ! Input  Data set number to transform
integer, dimension(2), INTENT(in)  :: odata     ! Output Data set number to transform
real(kind=PREC_DP)   , intent(in)  :: xscale
real(kind=PREC_DP)   , intent(in)  :: yscale
integer, dimension(3), intent(in)  :: ipad      ! Pad to these dimensions
!
integer :: i,j
integer :: kdat
integer :: length             ! Data set length == nx * ny
integer, dimension(3) :: num  ! DATA set dimensions
integer, dimension(3) :: dsort
complex(kind=kind(0.0D0)) , dimension(:,:), allocatable  :: k_data   ! The Kuplot in data set)
complex(kind=kind(0.0D0)) , dimension(:,:), allocatable  :: pattern  ! The Curve to be FFT'd
integer :: istart  ! Entry at which data start in k_data due to padding
integer :: jstart  ! Entry at which data start in k_data due to padding
logical, dimension(3) :: leven            ! Make original pattern symmetric around 0,0,0
integer, dimension(3) :: peven            ! Make original pattern symmetric around 0,0,0
!
!
if(idata(1)>0) then                             !  User provided real part
   call data2local(idata(1), ier_num, ier_typ, ik1_node_number, ik1_infile,     &
        ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
        ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
        ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
        ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
   kdat   = 1
   if(idata(2)==0) then
      ik2_ndims = ik1_ndims
      ik2_dims  = ik1_dims
      allocate(ik2_data(ik2_dims(1), ik2_dims(2),ik2_dims(3)))
      ik2_data = 0.0D0
   endif
endif
if(idata(2)>0) then                             !  User provided imag part
   call data2local(idata(2), ier_num, ier_typ, ik2_node_number, ik2_infile,     &
        ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
        ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
        ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
        ik2_llims, ik2_steps, ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
   kdat   = 2
   if(idata(1)==0) then
      ik1_ndims = ik2_ndims
      ik1_dims  = ik2_dims
      ik1_steps = ik2_steps
      allocate(ik1_data(ik1_dims(1), ik1_dims(2),ik1_dims(3)))
      ik1_data = 0.0D0
   else
      if(ik1_ndims/=ik2_ndims .or. any(ik1_dims/=ik2_dims)) then
         ier_num = -56
         ier_typ = ER_FORT
         call math_clean
         return
      endif
   endif
endif
!
!  Determine sequence of array dimensions 
!
if(idata(1)>0) then              ! Real is provided
   kdat   = 1
else                             ! Imag only
   kdat   = 2
endif
!
leven    = .false.
peven    = 0
if(mod(ik1_dims(1),2)==0) then      ! Even number of pixels
  leven(1) = .true.
  peven(1) = 1
endif
if(mod(ik1_dims(2),2)==0) then      ! Even number of pixels
  leven(2) = .true.
  peven(2) = 1
endif
!
num(1) = max(ipad(1), ik1_dims(1)         )  ! Pad with zeros if requested
num(2) = max(ipad(2), ik1_dims(2)         )  ! Pad with zeros if requested
num(3) = 1
istart = nint((num(1)-ik1_dims(1)+0.1D0)/2.0D0)
jstart = nint((num(2)-ik1_dims(2)+0.1D0)/2.0D0)
num(1) = max(ipad(1), ik1_dims(1)+peven(1))  ! Pad with symmetry if even numbered
num(2) = max(ipad(2), ik1_dims(2)+peven(2))  ! Pad with symmetry if even numbered
num(3) = 1
length = num(1)*num(2)
!
dsort(1)      = MAXLOC(num, 1)
num(dsort(1)) = -num(dsort(1))
dsort(2)      = MAXLOC(num, 1)
num(dsort(2)) = -num(dsort(2))
dsort(3)      = MAXLOC(num, 1)
num(dsort(3)) = -num(dsort(3))
num = -num
!
allocate(k_data (num(1), num(2)))
k_data = cmplx(0.0D0,0.0D0)
allocate(pattern(num(dsort(1)), num(dsort(2)) ))
do j=1, ik1_dims(2)
   do i=1, ik1_dims(1)
      k_data(istart+i,jstart+j) = CMPLX(ik1_data(i,j,1), ik2_data(i,j,1))
   enddo
enddo
if(leven(1)) then 
   do j=1, ik1_dims(2)
      k_data(istart+ik1_dims(1)+1,jstart+j) = k_data(istart+1, jstart+ik1_dims(2)+1+peven(2)-j)
   enddo
endif
if(leven(2)) then 
   do i=1, ik1_dims(1)
      k_data(istart+i, jstart+ik1_dims(2)+1) = k_data(istart+ik1_dims(1)+1+peven(1)-i, jstart+1)
   enddo
endif
!
call maptofftfd(num, dsort, k_data, pattern)
!
pattern = fft(pattern) / SQRT(REAL(num(1)*num(2))) * yscale
!
call mapfftfdtoline(num, dsort, k_data, pattern)
!
ft_start(1:2) = -0.5D0/ik1_steps(1:2)
ft_start(3)   =  0.00D0
ft_steps(1:2) =  1.0D0/(ik1_steps(1:2)*real(num(1:2)-1,kind=PREC_DP))
ft_steps(3)   =  0.0D0
!
allocate(res_rl_data(num(1), num(2), 1))
allocate(res_im_data(num(1), num(2), 1))
res_rl_data(:,:,1) =  real(k_data(:,:),kind=PREC_DP)
res_im_data(:,:,1) = dimag(k_data(:,:))
!
ik1_dims  = num          ! Copy current dimensions 
steps     = ft_steps
ik1_llims = ft_start
ik1_nlayer      = 1   ! Dummy layer number
!
! Store real/imag pair into global storage
!
ik1_node_number = 0
call fft2data(odata(1), ier_num, ier_typ, ik1_node_number, ik1_nlayer,          &
     ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
     ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
     ik1_z, ik1_dx, ik1_dy, ik1_dz, res_rl_data, ik1_sigma, ik1_llims, steps,   &
     ik1_steps_full, res_im_data, ik2_sigma)
!
!
!
deallocate(k_data)
deallocate(pattern)
deallocate(res_rl_data)
deallocate(res_im_data)
call math_clean
!
end subroutine fft_2D_global
!
!*****7**************************************************************** 
!
SUBROUTINE fft_3D_global(idata, odata, xscale, yscale, ipad)
!-
!  Calculate 3D FFT of data set idata
!+
!
use lib_data_struc_h5
!
USE errlist_mod
USE map_1dtofield
USE prompt_mod
USE singleton
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(2), INTENT(IN)  :: idata     ! Data set number to transform
INTEGER, DIMENSION(2), INTENT(OUT) :: odata     ! Data set number to transform
real(kind=PREC_DP)   , intent(in)  :: xscale
real(kind=PREC_DP)   , intent(in)  :: yscale
integer, dimension(3), intent(in)  :: ipad      ! Pad to these dimensions
!
INTEGER :: i, j, k
INTEGER :: kdat
INTEGER :: length             ! Data set length == nx * ny
INTEGER, DIMENSION(3) :: num  ! DATA set dimensions
INTEGER, DIMENSION(3) :: dsort
complex(KIND=KIND(0.0D0)) , DIMENSION(:,:,:), ALLOCATABLE  :: k_data   ! The Kuplot in data set)
COMPLEX(KIND=KIND(0.0D0)) , DIMENSION(:,:,:), ALLOCATABLE  :: pattern  ! The Curve to be FFT'd
integer :: istart  ! Entry at which data start in k_data due to padding
integer :: jstart  ! Entry at which data start in k_data due to padding
integer :: kstart  ! Entry at which data start in k_data due to padding
logical, dimension(3) :: leven            ! Make original pattern symmetric around 0,0,0
integer, dimension(3) :: peven            ! Make original pattern symmetric around 0,0,0
!
!
if(idata(1)>0) then                             !  User provided real part
   call data2local(idata(1), ier_num, ier_typ, ik1_node_number, ik1_infile,     &
        ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
        ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
        ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
        ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
   kdat   = 1
   if(idata(2)==0) then
      ik2_ndims = ik1_ndims
      ik2_dims  = ik1_dims
      ik2_llims = ik1_llims
      ik2_steps = ik1_steps
      allocate(ik2_data(ik2_dims(1), ik2_dims(2),ik2_dims(3)))
      ik2_data = 0.0D0
   endif
endif
if(idata(2)>0) then                             !  User provided imag part
   call data2local(idata(2), ier_num, ier_typ, ik2_node_number, ik2_infile,     &
        ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
        ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
        ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
        ik2_llims, ik2_steps, ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
   kdat   = 2
   if(idata(1)==0) then
      ik1_ndims = ik2_ndims
      ik1_dims  = ik2_dims
      ik1_llims = ik2_llims
      ik1_steps = ik2_steps
      allocate(ik1_data(ik1_dims(1), ik1_dims(2),ik1_dims(3)))
      ik1_data = 0.0D0
   else
      if(ik1_ndims/=ik2_ndims .or. any(ik1_dims/=ik2_dims)) then
         ier_num = -56
         ier_typ = ER_FORT
         call math_clean
         return
      endif
   endif
endif
!
!  Determine sequence of array dimensions 
!
IF(idata(1)>0) THEN              ! Real is provided
   kdat   = 1
ELSE                             ! Imag only
   kdat   = 2
ENDIF
!
peven = 0
leven = .false.
!
if(mod(ik1_dims(1),2)==0) then      ! Even number of pixels
  leven(1) = .true.
  peven(1) = 1
endif
if(mod(ik1_dims(2),2)==0) then      ! Even number of pixels
  leven(2) = .true.
  peven(2) = 1
endif
if(mod(ik1_dims(3),2)==0) then      ! Even number of pixels
  leven(3) = .true.
  peven(3) = 1
endif
!
num(1) = max(ipad(1), ik1_dims(1)         )  ! Pad with zeros if requested
num(2) = max(ipad(2), ik1_dims(2)         )  ! Pad with zeros if requested
num(3) = max(ipad(3), ik1_dims(3)         )  ! Pad with zeros if requested
!
istart = nint((num(1)-ik1_dims(1)+0.1D0)/2.0D0)
jstart = nint((num(2)-ik1_dims(2)+0.1D0)/2.0D0)
kstart = nint((num(3)-ik1_dims(3)+0.1D0)/2.0D0)
!
num(1) = max(ipad(1), ik1_dims(1)+peven(1))  ! Pad with zeros if requested
num(2) = max(ipad(2), ik1_dims(2)+peven(2))  ! Pad with zeros if requested
num(3) = max(ipad(3), ik1_dims(3)+peven(3))  ! Pad with zeros if requested
length = num(1)*num(2)*num(3)
!
dsort(1)      = MAXLOC(num, 1)
num(dsort(1)) = -num(dsort(1))
dsort(2)      = MAXLOC(num, 1)
num(dsort(2)) = -num(dsort(2))
dsort(3)      = MAXLOC(num, 1)
num(dsort(3)) = -num(dsort(3))
num = -num
!
!
ALLOCATE(k_data (num(1)       , num(2)       ,num(3)         ))
ALLOCATE(pattern(num(dsort(1)), num(dsort(2)), num(dsort(3)) ))
!
k_data = cmplx(0.0D0,0.0D0)        ! Initialize for padding
do k=1, ik1_dims(3)
   do j=1, ik1_dims(2)
      do i=1, ik1_dims(1)
         k_data(istart+i,jstart+j,kstart+k) = CMPLX(ik1_data(i,j,k), ik2_data(i,j,k))
      enddo
   enddo
enddo
!
!Make centrosymmetric
if(leven(1)) then
! (100)
do k=1, ik1_dims(3)
   do j=1, ik1_dims(2)
      k_data(istart+ik1_dims(1)+1, jstart+j                      , kstart+k                      ) =  &
      k_data(istart+1           , jstart+ik1_dims(2)+1+peven(2)-j, kstart+ik1_dims(3)+1+peven(3)-k)
   enddo
enddo
endif
!
if(leven(1)) then
! (010)
do k=1, ik1_dims(3)
   do i=1, ik1_dims(1)
      k_data(istart+i                      , jstart+ik1_dims(2)+1, kstart+k                      ) =  &
      k_data(istart+ik1_dims(1)+1+peven(1)-i, jstart+1           , kstart+ik1_dims(3)+1+peven(3)-k)
   enddo
enddo
endif
!
if(leven(1)) then
! (001)
do j=1, ik1_dims(2)
   do i=1, ik1_dims(1)
      k_data(istart+i                      , jstart+j                      , kstart+ik1_dims(3)+1) =  &
      k_data(istart+ik1_dims(1)+1+peven(1)-j, jstart+ik1_dims(2)+1+peven(2)-j, kstart+1)
   enddo
enddo
endif
!
!
CALL maptofftfd(num, dsort, k_data, pattern)
!
pattern = fft(pattern) / SQRT(REAL(num(1)*num(2)*num(3))) * yscale
!
CALL mapfftfdtoline(num, dsort, k_data, pattern)
!
!  Convert limits, copy into result array
!
ft_start = -0.5D0/ik1_steps
ft_steps =  1.0D0/(ik1_steps*real(num-1,kind=PREC_DP))
!
allocate(res_rl_data(num(1), num(2), 1))
allocate(res_im_data(num(1), num(2), 1))
res_rl_data =  real(k_data,kind=PREC_DP)
res_im_data = dimag(k_data)
!
ik1_dims  = num    ! Copy current dimensions
steps    = ft_steps
ik1_llims = ft_start
ik1_node_number = 0
!
! Replace the current image with the central layer
!
IF(MOD(num(3),2)==0) THEN
  ik1_nlayer = num(3)/2
ELSE
  ik1_nlayer = (num(3)+1)/2
ENDIF
!                 ! Copy result into global storage
call fft2data(odata(1), ier_num, ier_typ, ik1_node_number, ik1_nlayer,          &
     ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
     ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
     ik1_z, ik1_dx, ik1_dy, ik1_dz, res_rl_data, ik1_sigma, ik1_llims, steps,   &
     ik1_steps_full, res_im_data, ik2_sigma)
!
DEALLOCATE(k_data)
DEALLOCATE(pattern)
deallocate(res_rl_data)
deallocate(res_im_data)
call math_clean
!
end subroutine fft_3D_global
!
!*****7*****************************************************************
!
subroutine kmath_h5_global(ik1, ik2, oper, ik3)
!-
!  Perform data set products / additions for data in global storage
!+
!
use errlist_mod
use lib_data_struc_h5
use precision_mod
use prompt_mod
!
implicit none
!
integer         , intent(in) :: ik1
integer         , intent(in) :: ik2
character(Len=*), intent(in) :: oper
integer         , intent(in) :: ik3
!
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: res_data        ! Global data array for real
!
call data2local(ik1     , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
call data2local(ik2     , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
     ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
     ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
     ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
     ik2_llims, ik2_steps,  ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
!
if(ik1_ndims/=ik2_ndims .or. any(ik1_dims/=ik2_dims)) then
   ier_num = -56
   ier_typ = ER_FORT
   call math_clean
   return
endif
!
allocate(res_data(ik1_dims(1), ik1_dims(2), ik1_dims(3)))
!
if(oper=='ADD') then
   res_data = ik1_data + ik2_data
!
elseif(oper=='SUB') then
   res_data = ik1_data - ik2_data
!
elseif(oper=='MUL') then
   res_data = ik1_data * ik2_data
!
elseif(oper=='DIV') then
   res_data = 0.0D0
   where(ik2_data/=0.0D0)
      res_data = ik1_data /ik2_data
   end where
!
endif
!
ik1_node_number = 0
call local2data(ik3, ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_nlayer,  &
     ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
     ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
     ik1_z, ik1_dx, ik1_dy, ik1_dz, res_data, ik1_sigma, ik1_llims, ik1_steps,  &
     ik1_steps_full)
!
deallocate(res_data)
call math_clean
!
end subroutine kmath_h5_global
!
!*****7*****************************************************************
!
subroutine rvalue_h5_global(ik1, ik2, iweight, bck_k, rval, wrval)
!-
!  Perform R-value calculations
!+
!
use errlist_mod
use lib_data_struc_h5
use lib_weights_mod
use precision_mod
!use prompt_mod
!
implicit none
!
integer           , intent(in) :: ik1
integer           , intent(in) :: ik2
integer           , intent(in)  :: iweight
real(kind=PREC_DP), intent(out) :: rval 
real(kind=PREC_DP), intent(out) :: wrval 
real(kind=PREC_DP), intent(in)  :: bck_k
!                                                                       
integer :: ii,jj,kk
real(kind=PREC_DP) :: sumrz, sumrn 
real(kind=PREC_DP) :: sumwrz, sumwrn 
real(kind=PREC_DP) :: wght 
real(kind=PREC_DP) :: a_ik1, a_ik2 

!
call data2local(ik1     , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
call data2local(ik2     , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
     ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
     ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
     ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
     ik2_llims, ik2_steps,  ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
!
if(ik1_ndims/=ik2_ndims .or. any(ik1_dims/=ik2_dims)) then
   ier_num = -56
   ier_typ = ER_FORT
   call math_clean
   return
endif
!
sumrz  = 0.0D0
sumrn  = 0.0D0
sumwrz = 0.0D0
sumwrn = 0.0D0
!
do kk=1,ik1_dims(3)
   do jj=1,ik1_dims(2)
      do ii=1,ik1_dims(1)
         a_ik1 = ik1_data(ii,jj,kk)
         a_ik2 = ik2_data(ii,jj,kk)
         sumrz = sumrz + abs (a_ik1 -a_ik2)
         sumrn = sumrn + abs (a_ik1) 
         wght = r_wichtung (a_ik1  , 1.0, iweight, a_ik2, bck_k) 
         sumwrz = sumwrz + wght * (a_ik1   - a_ik2   ) **2 
         sumwrn = sumwrn + wght * (a_ik1   ) **2 
      enddo
   enddo
enddo
!
rval = sumrz / sumrn 
wrval = sqrt (sumwrz / sumwrn) 
!                                                                       
call math_clean
!
end subroutine rvalue_h5_global
!
!*****7*****************************************************************
!
subroutine do_merge_global (zeile, lp, ik3) 
!+                                                                      
!     Merge different data sets                                         
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
use lib_data_struc_h5
!OUSE kuplot_config 
!OUSE kuplot_mod 
!Ouse kuplot_extrema_mod
!Ouse kuplot_show_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
INTEGER          , INTENT(INOUT) :: ik3
!
INTEGER, parameter :: maxw = 30
!PARAMETER (maxw = MAXKURVTOT) 
!                                                                       
CHARACTER(LEN=PREC_STRING), dimension(:), allocatable :: cpara !(maxw) 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_DP) :: mdelta, mmin, mmax 
REAL(kind=PREC_DP) ::  mdeltay, mminy, mmaxy 
INTEGER :: lpara (maxw)
INTEGER :: ianz, ik, ibin, i, ip, ntot
integer :: ibiny
INTEGER :: ntoty
INTEGER :: ix, iy       ! Dummy loop indices
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: npt ! number of data points
LOGICAL :: lvalid, ladd, lall 
!                                                                       
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
!IF (iz.gt.maxkurvtot) then 
!   ier_num = - 1 
!   ier_typ = ER_APPL 
!   RETURN 
!ENDIF 
!
allocate(cpara(maxw))
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) then
   deallocate(cpara)
   return 
endif
!                                                                       
ladd = str_comp (cpara (ianz) , 'add', 1, lpara (ianz) , 3) 
IF (ladd) ianz = ianz - 1 
!                                                                       
lall = str_comp (cpara (ianz) , 'all', 1, lpara (ianz) , 3) 
IF (lall) then 
   werte (1) = - 1 
   deallocate(cpara)
ELSE 
!                                                                       
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   deallocate(cpara)
   IF (ier_num.ne.0) then
      return 
   endif
ENDIF 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF(nint(werte(1)) /= -1 .and. ianz < 2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
IF(ianz == 1) then 
   IF(nint(werte(1)) == -1) then 
      DO i = 1, dgl5_get_number()
         werte(i) = i 
      ENDDO 
      ianz = dgl5_get_number()
   ENDIF 
else
   if(nint(maxval(werte))>dgl5_get_number()) then
      ier_num = -4
      ier_typ = ER_APPL
      return
   endif
ENDIF 
!                                                                       
!------ First we get extend of new data set (first one defines DELTA)   
!                                                                       
ik = nint(werte(1)) 
!
call data2local(ik      , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
if(ik1_ndims==1) then                    ! 1-D data sets
!
   mdelta = (ik1_x(ik1_dims(1)) - ik1_x(1) )/real(ik1_dims(1)-1,kind=PREC_DP)
   mmin   = ik1_minmaxcoor(1,1)
   mmax   = ik1_minmaxcoor(1,2)
!
   lvalid = .true. 
   loop_valid_1d: do i = 1, ianz 
      ik = nint(werte(i))
      call data2local(ik      , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
           ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
           ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
           ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
           ik2_llims, ik2_steps,  ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
      lvalid = lvalid .and. ik2_ndims==1
      if(.not.lvalid) exit loop_valid_1d
      mmin   = min(mmin, ik2_minmaxcoor(1,1))
      mmax   = max(mmax, ik2_minmaxcoor(1,2))
   enddo loop_valid_1d
!
   IF (.not.lvalid) then 
      ier_num = - 4 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!
   ntot = nint((mmax - mmin) / mdelta) + 1 
!   
   deallocate(ik1_x)
   deallocate(ik1_y)
   deallocate(ik1_z)
   if(allocated(ik1_x)) deallocate(ik1_x)
   if(allocated(ik1_y)) deallocate(ik1_y)
   if(allocated(ik1_z)) deallocate(ik1_z)
   if(allocated(ik1_data))  deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   ik1_dims(1) = ntot
   ik1_dims(2) = 1
   ik1_dims(3) = 1
   allocate(ik1_x(1:ik1_dims(1)))
   allocate(ik1_y(1))
   allocate(ik1_z(1))
   allocate(ik1_data(ik1_dims(1), ik1_dims(2), ik1_dims(3)))
   allocate(ik1_sigma(ik1_dims(1), ik1_dims(2), ik1_dims(3)))
   allocate(npt      (ik1_dims(1), ik1_dims(2), ik1_dims(3)))
   ik1_x = 0.0D0
   ik1_y = 0.0D0
   ik1_z = 0.0D0
   npt   = 0
!
   loop_merge_1d: do i=1, ianz    ! Merge values of all data
      ik = nint(werte(i))
      call data2local(ik      , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
           ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
           ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
           ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
           ik2_llims, ik2_steps,  ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
      do ip = 1, ik2_dims(1)
         ibin = nint((ik1_x(ip) - mmin) / mdelta) + 1
         ik1_data(ibin,1,1) = ik1_data(ibin,1,1) + ik2_data(ip,1,1)
         if(ik2_has_dval) ik1_sigma(ibin,1,1) = ik1_sigma(ibin,1,1) + ik2_sigma(ibin,1,1)**2
         npt   (ibin,1,1) = npt   (ibin,1,1) + 1
      enddo
   enddo loop_merge_1d
   if(.not.ladd) then                      !Merge data, not just add
      do ibin=1, ntot
         if(npt(ibin,1, 1)>0) then 
            ik1_data(ibin,1,1)  = ik1_data(ibin,1,1)/npt(ibin,1,1)
            ik1_sigma(ibin,1,1) = sqrt(ik1_data(ibin,1,1))/npt(ibin,1,1)
         endif
      enddo
   endif
   do ibin=1, ntot
      ik1_x(ibin) = mmin + (ibin - 1) * mdelta 
   ENDDO 
   deallocate(npt)
   ik1_is_grid = .true.
   ik1_infile  = 'merged.dat'
!
   ik1_node_number = 0
   call local2data(ik3, ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_nlayer,  &
        ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
        ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
        ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma, ik1_llims, ik1_steps,  &
        ik1_steps_full)
!
   call math_clean
!
elseif(ik1_ndims==2) then                ! 2-D data sets
!
   mdelta  = (ik1_x(ik1_dims(1)) - ik1_x(1) )/real(ik1_dims(1)-1,kind=PREC_DP)
   mdeltay = (ik1_y(ik1_dims(2)) - ik1_y(1) )/real(ik1_dims(2)-1,kind=PREC_DP)
   mmin    = ik1_minmaxcoor(1,1)
   mminy   = ik1_minmaxcoor(2,1)
   mmax    = ik1_minmaxcoor(1,2)
   mmaxy   = ik1_minmaxcoor(2,2)
!
   lvalid = .true. 
   loop_valid_2d: do i = 1, ianz 
      ik = nint(werte(i))
      call data2local(ik      , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
           ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
           ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
           ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
           ik2_llims, ik2_steps,  ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
      lvalid = lvalid .and. ik2_ndims==2
      if(.not.lvalid) exit loop_valid_2d
      mmin   = min(mmin , ik2_minmaxcoor(1,1))
      mminy  = min(mminy, ik2_minmaxcoor(2,1))
      mmax   = max(mmax , ik2_minmaxcoor(1,2))
      mmaxy  = max(mmaxy, ik2_minmaxcoor(2,2))
   enddo loop_valid_2d
!
   IF (.not.lvalid) then 
      ier_num = - 4 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!
   ntot  = nint((mmax  - mmin)  / mdelta ) + 1 
   ntoty = nint((mmaxy - mminy) / mdeltay) + 1 
!   
   deallocate(ik1_x)
   deallocate(ik1_y)
   deallocate(ik1_z)
   if(allocated(ik1_x)) deallocate(ik1_x)
   if(allocated(ik1_y)) deallocate(ik1_y)
   if(allocated(ik1_z)) deallocate(ik1_z)
   if(allocated(ik1_data))  deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   ik1_dims(1) = ntot
   ik1_dims(2) = 1
   ik1_dims(3) = 1
   allocate(ik1_x(1:ik1_dims(1)))
   allocate(ik1_y(1:ik1_dims(2)))
   allocate(ik1_z(1))
   allocate(ik1_data( ik1_dims(1), ik1_dims(2), ik1_dims(3)))
   allocate(ik1_sigma(ik1_dims(1), ik1_dims(2), ik1_dims(3)))
   allocate(npt      (ik1_dims(1), ik1_dims(2), ik1_dims(3)))
   ik1_x = 0.0D0
   ik1_y = 0.0D0
   ik1_z = 0.0D0
   npt   = 0
   ik1_data  = 0.0D0
   ik1_sigma = 0.0D0
!
   loop_merge_2d: do i=1, ianz    ! Merge values of all data
      ik = nint(werte(i))
      call data2local(ik      , ier_num, ier_typ, ik2_node_number, ik2_infile,     &
           ik2_nlayer, ik2_is_direct, ik2_ndims, ik2_dims, ik2_is_grid,            &
           ik2_has_dxyz, ik2_has_dval, ik2_corners, ik2_vectors, ik2_a0, ik2_win,  &
           ik2_x, ik2_y, ik2_z, ik2_dx, ik2_dy, ik2_dz, ik2_data, ik2_sigma,       &
           ik2_llims, ik2_steps,  ik2_steps_full, ik2_minmaxval, ik2_minmaxcoor)
      do iy = 1, ik2_dims(2)
      do ix = 1, ik2_dims(1)
         ibin  = nint((ik1_x(ix) - mmin ) / mdelta ) + 1
         ibiny = nint((ik1_x(iy) - mminy) / mdeltay) + 1
         ik1_data(ibin,ibiny,1) = ik1_data(ibin,ibiny,1) + ik2_data(ix,iy,1)
         if(ik2_has_dval) ik1_sigma(ibin,ibiny,1) = ik1_sigma(ibin,ibiny,1) + ik2_sigma(ibin,ibiny,1)**2
         npt(ibin, ibiny, 1) = npt(ibin, ibiny, 1) + 1
      enddo
      enddo
   enddo loop_merge_2d
   if(.not.ladd) then                      !Merge data, not just add
      do ibin=1, ntot
         do ibiny=1, ntoty
         if(npt(ibin, ibiny, 1)>0) then 
            ik1_data(ibin,ibiny,1)  = ik1_data(ibin,ibiny,1)/npt(ibin, ibiny, 1)
            ik1_sigma(ibin,ibiny,1) = sqrt(ik1_data(ibin,ibiny,1))/npt(ibin, ibiny, 1)
         endif
         enddo
      enddo
   endif
   do ibin=1, ntot
      ik1_x(ibin) = mmin + (ibin - 1) * mdelta 
   ENDDO 
   do ibiny=1, ntoty
      ik1_y(ibiny) = mminy + (ibiny - 1) * mdeltay 
   ENDDO 
   deallocate(npt)
   ik1_is_grid = .true.
   ik1_infile = 'merged.dat'
!
   ik1_node_number = 0
   call local2data(ik3, ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_nlayer,  &
        ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid, ik1_has_dxyz,             &
        ik1_has_dval, ik1_corners, ik1_vectors, ik1_a0, ik1_win, ik1_x, ik1_y,     &
        ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma, ik1_llims, ik1_steps,  &
        ik1_steps_full)
!
   call math_clean
   deallocate(npt)
ENDIF
!
!                                                                       
END SUBROUTINE do_merge_global
!
!*******************************************************************************
!
subroutine math_clean
!-
!  deallocate all local arrays
!+
if(allocated(ik1_x    )) deallocate(ik1_x    )
if(allocated(ik1_y    )) deallocate(ik1_y    )
if(allocated(ik1_z    )) deallocate(ik1_z    )
if(allocated(ik1_dx   )) deallocate(ik1_dx   )
if(allocated(ik1_dy   )) deallocate(ik1_dy   )
if(allocated(ik1_dz   )) deallocate(ik1_dz   )
if(allocated(ik1_data )) deallocate(ik1_data )
if(allocated(ik1_sigma)) deallocate(ik1_sigma)
!
if(allocated(ik2_x    )) deallocate(ik2_x    )
if(allocated(ik2_y    )) deallocate(ik2_y    )
if(allocated(ik2_z    )) deallocate(ik2_z    )
if(allocated(ik2_dx   )) deallocate(ik2_dx   )
if(allocated(ik2_dy   )) deallocate(ik2_dy   )
if(allocated(ik2_dz   )) deallocate(ik2_dz   )
if(allocated(ik2_data )) deallocate(ik2_data )
if(allocated(ik2_sigma)) deallocate(ik2_sigma)
!
end subroutine math_clean
!
!*******************************************************************************
!
end module lib_math_mod
