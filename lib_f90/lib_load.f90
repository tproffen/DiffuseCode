module lib_load_mod
!
!  Load data sets into the global 3D data structure
!
private
public gen_load_hklf4
!
contains
!
!*******************************************************************************
!
subroutine gen_load_hklf4(infile, node_number, lout)
!-
! Load a SHELX type HKLF4 file into memory
!+
use errlist_mod
use lib_data_struc_h5
use lib_data_types_mod
use precision_mod
use prompt_mod
!
implicit none
!
character(len=*), intent(in) :: infile
integer, intent(out) :: node_number
logical, intent(in ) :: lout
!
integer, parameter :: IRD = 33
!
character(len=PREC_STRING) :: line
integer               :: i          ! Dummy loop index
integer               :: ios        ! I/O status
integer, dimension(3) :: minhkl     ! Max values for hkl
integer, dimension(3) :: maxhkl     ! Max values for hkl
integer, dimension(3) ::  ihkl      ! An individual hkl
integer, dimension(3) ::  jhkl      ! An individual hkl shifted by minhkl
integer               :: nhkl       ! number reflections
integer               :: nobs       ! number reflections > 3sigma
logical               :: lda        ! File is present / absent
real(kind=PREC_DP), dimension(2)    :: inte_range       ! intensity for reflection j
real(kind=PREC_DP), dimension(2)    ::  sig_range       ! intensity for reflection j
real(kind=PREC_DP)    :: inte       ! intensity for reflection j
real(kind=PREC_DP)    :: sig        ! sigma     for reflection j
real(kind=PREC_DP), dimension(:,:,:), allocatable :: observed
real(kind=PREC_DP), dimension(:,:,:), allocatable :: sigma
integer           , dimension(:,:,:), allocatable :: weight
!
integer :: d5_layer
integer :: d5_nndims
integer :: d5_data_type
integer,dimension(3) :: d5_dims
integer,dimension(3) :: d5_use_coor
logical :: d5_direct
logical :: d5_is_grid
logical :: d5_has_dxyz
logical :: d5_has_dval
logical :: d5_calc_coor
real(kind=PREC_DP)   , dimension(3)                   :: d5_llims          ! left lower bottom corner
real(kind=PREC_DP)   , dimension(3)                   :: d5_steps          ! simple steps along axes 
real(kind=PREC_DP)   , dimension(3,3)                 :: d5_steps_full     ! full steps in H, K, L
real(kind=PREC_DP)   , dimension(3,3)                 :: d5_vectors        ! full steps in H, K, L
real(kind=PREC_DP)   , dimension(3,4)                 :: d5_corners        ! steps in H, K, L
real(kind=PREC_DP)   , dimension(6)                   :: d5_unit           ! steps in H, K, L
real(kind=PREC_DP)   , dimension(:), allocatable      :: d5_x              ! h coordinates
real(kind=PREC_DP)   , dimension(:), allocatable      :: d5_y              ! k coordinates
real(kind=PREC_DP)   , dimension(:), allocatable      :: d5_z              ! l coordinates
real(kind=PREC_DP)   , dimension(:), allocatable      :: d5_dx             ! h coordinates
real(kind=PREC_DP)   , dimension(:), allocatable      :: d5_dy             ! k coordinates
real(kind=PREC_DP)   , dimension(:), allocatable      :: d5_dz             ! l coordinates
!
ios = 0
inquire(file=infile,exist=lda)
if(.not.lda) then
   ier_num = -1
   ier_typ = ER_IO
   return
endif
!
open(unit=IRD, file=infile, iostat=ios)
if(ios/=0) then
   ier_num = -2
   ier_typ = ER_IO
   return
endif
!
!  Read once to get maximum abs(hkl)
!
inte_range(1) = +huge(1.0D0)
inte_range(2) = -huge(1.0D0)
 sig_range(1) = +huge(1.0D0)
 sig_range(2) = -huge(1.0D0)
minhkl =  9999
maxhkl = -9999
nhkl = 0
nobs = 0
read(IRD, '(a)', iostat=ios) line
init_loop: do
   if(is_iostat_end(ios)) exit init_loop
   read(line, '(3i4,2f8.2)', iostat=ios) ihkl
   if(is_iostat_end(ios)) exit init_loop
   if(all(ihkl==0)) exit init_loop
   nhkl = nhkl + 1
   minhkl(1) = min(minhkl(1),     ihkl(1))
   minhkl(2) = min(minhkl(2),     ihkl(2))
   minhkl(3) = min(minhkl(3),     ihkl(3))
   maxhkl(1) = max(maxhkl(1),     ihkl(1))
   maxhkl(2) = max(maxhkl(2),     ihkl(2))
   maxhkl(3) = max(maxhkl(3),     ihkl(3))
   read(IRD, '(a)', iostat=ios) line
enddo init_loop
close(IRD)
!
allocate(observed( maxhkl(1)-minhkl(1)+1, maxhkl(2)-minhkl(2)+1, maxhkl(3)-minhkl(3)+1))
allocate(sigma   ( maxhkl(1)-minhkl(1)+1, maxhkl(2)-minhkl(2)+1, maxhkl(3)-minhkl(3)+1))
allocate(weight  ( maxhkl(1)-minhkl(1)+1, maxhkl(2)-minhkl(2)+1, maxhkl(3)-minhkl(3)+1))
!
observed = 0.0D0
sigma    = -10000.0D0    ! Flag all points as missing
weight   = 0
!
open(unit=IRD, file=infile, iostat=ios)
read_loop: do i=1, nhkl
   read(IRD, '(a)', iostat=ios) line
   if(is_iostat_end(ios)) exit read_loop
   read(line, '(3i4,2f8.2)', iostat=ios) ihkl, inte, sig
   if(is_iostat_end(ios)) exit read_loop
   if(all(ihkl==0)) exit read_loop
   sig = max(0.001D0,sig)
!
!if(weight(ihkl(1), ihkl(2), ihkl(3))/=0) then  ! new reflection
!write(output_io,'(a, i6,2x,3i4)') 'DOUBLE AT ', i, ihkl
!endif
    jhkl = ihkl - minhkl + 1
    inte_range(1) = min(inte_range(1), inte)
     sig_range(1) = min( sig_range(1),  sig)
    inte_range(2) = max(inte_range(2), inte)
     sig_range(2) = max( sig_range(2),  sig)
    if(inte>3.0D0*abs(sig)) nobs = nobs + 1
    observed(jhkl(1), jhkl(2), jhkl(3)) = observed(jhkl(1), jhkl(2), jhkl(3)) + inte
    if(weight  (jhkl(1), jhkl(2), jhkl(3)) ==0) then 
       sigma   (jhkl(1), jhkl(2), jhkl(3)) = abs(sig)
    else
       sigma   (jhkl(1), jhkl(2), jhkl(3)) = sigma   (jhkl(1), jhkl(2), jhkl(3)) + abs(sig)
    endif
    weight  (jhkl(1), jhkl(2), jhkl(3)) = weight  (jhkl(1), jhkl(2), jhkl(3)) + 1
!
enddo read_loop
!
!write(output_io, '(a, 4i6)') 'Maximum weight ', maxval(weight)
where(weight>0)
   observed = observed/weight
   sigma    = sigma   /weight
   weight   = 1
else where
   sigma = -1.0D0
end where
!write(output_io, '(a, 4i6)') 'Maximum weight ', maxval(weight)
if(lout) then
   write(output_io,*)
   write(output_io,'(a)')         ' Intensity statistics ' 
   write(output_io,'(a, 3i6)')    '  # reflections', nhkl
   write(output_io,'(a, 3i6)')    '  # I>3sigma   ', nobs
   write(output_io,'(a, 3i6)')    '  Minimum hkl  ', minhkl
   write(output_io,'(a, 3i6)')    '  Maximum hkl  ', maxhkl
   write(output_io,'(a, 2f12.4)') '  Intensity range ', inte_range
   write(output_io,'(a, 2f12.4)') '  Sigma     range ',  sig_range
   write(output_io,*)
endif
!
close(IRD)
!
d5_layer  = max(abs(lbound(observed,3)),-minhkl(3)+1)
d5_direct = .false.
d5_nndims = 3
d5_dims   = (maxhkl-minhkl) + 1
d5_llims   = minhkl
d5_steps   = 1.0D0
d5_steps_full      = 0.0D0
d5_steps_full(1,1) = 1.0D0
d5_steps_full(2,2) = 1.0D0
d5_steps_full(3,3) = 1.0D0
d5_is_grid= .true.
d5_has_dxyz=.false.
d5_has_dval=.true.
d5_calc_coor=.false.
d5_use_coor = (/1,2,3/)
d5_corners(:,1) = minhkl(:)
d5_corners(:,2) = minhkl(:)
d5_corners(:,3) = minhkl(:)
d5_corners(:,4) = minhkl(:)
d5_corners(1,2) = maxhkl(1)
d5_corners(2,3) = maxhkl(2)
d5_corners(3,4) = maxhkl(3)
d5_vectors      = d5_steps_full
d5_unit(1:3)    = 1.0D0
d5_unit(4:6)    =90.0D0
allocate(d5_x(d5_dims(1)))
allocate(d5_y(d5_dims(2)))
allocate(d5_z(d5_dims(3)))
allocate(d5_dx(d5_dims(1)))
allocate(d5_dy(d5_dims(2)))
allocate(d5_dz(d5_dims(3)))
do i=1,d5_dims(1)
   d5_x(i) = minhkl(1) + real((i-1),kind=PREC_DP)
enddo
do i=1,d5_dims(2)
   d5_y(i) = minhkl(2) + real((i-1),kind=PREC_DP)
enddo
do i=1,d5_dims(3)
   d5_z(i) = minhkl(3) + real((i-1),kind=PREC_DP)
enddo
d5_dx = 0.0D0
d5_dy = 0.0D0
d5_dz = 0.0D0
call dgl5_new_node
node_number = dgl5_get_number()
d5_data_type = H5_BRAGG_I
call dgl5_set_node(   infile, d5_data_type, d5_layer, d5_direct, d5_nndims, d5_dims ,         &
                   d5_is_grid, d5_has_dxyz, d5_has_dval, d5_calc_coor, d5_use_coor, &
                   d5_corners, d5_vectors,&
                   d5_unit(1:3), d5_unit(4:6), d5_x, d5_y, d5_z, d5_dx, d5_dy,  &
                   d5_dz,      observed              ,    sigma, d5_llims,      &
                   d5_steps, d5_steps_full)
!
deallocate(observed)
deallocate(sigma)
deallocate(weight)
deallocate(d5_x)
deallocate(d5_y)
deallocate(d5_z)
deallocate(d5_dx)
deallocate(d5_dy)
deallocate(d5_dz)
!
end subroutine gen_load_hklf4
!
!*******************************************************************************
!
end module lib_load_mod
