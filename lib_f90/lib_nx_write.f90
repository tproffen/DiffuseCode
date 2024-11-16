module nx_write_mod
!-
!  Subroutines to write into a NEXUS file
!  DiffuseDevelopers Software convention
!+
private
!
public nx_write_scattering  ! Subroutine; write the provided scattering data into Nexus file
!
contains
!
!*******************************************************************************
!
subroutine nx_write_scattering(outfile, out_inc, out_eck, out_vi,               &
           abs_is_hkl, ord_is_hkl, top_is_hkl,                                  &
           cr_a0, cr_win, qvalues, content, is_space, radiation,                &
           space_group_name, symmetry_applied, symmetry_n_mat, symmetry_mat,    &
           ier_num)
!-
!  Main writing routine for Nexus Format DiffuseDevelopers Scattering format
!+
!
use precision_mod
!
use lib_forpython_mod
use forpy_mod
!
implicit none
!
character(len=*)                  , intent(in)  :: outfile     ! Output file name
integer           , dimension(3)  , intent(in)  :: out_inc     ! Dimensions of qvalue
real(kind=PREC_DP), dimension(3,4), intent(in)  :: out_eck     ! (:,1) LowerLeft
real(kind=PREC_DP), dimension(3,3), intent(in)  :: out_vi      ! (:,1) Abscissa
!                                                              ! (:,2) Ordinate
!                                                              ! (:,3) Top axis
integer                           , intent(in)  :: abs_is_hkl  ! Abscissa is 1=h 2=k 3=l
integer                           , intent(in)  :: ord_is_hkl  ! Ordinate is "
integer                           , intent(in)  :: top_is_hkl  ! Top axis is "
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_a0       ! Lattice params (a, b, c)
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_win      ! Lattice params (alpha, beta, gamma)
real(kind=PREC_DP), dimension(out_inc(1), out_inc(2), out_inc(3)), intent(in) :: qvalues  ! Actual data array
character(len=*)                  , intent(in)  :: content     ! The data contain this information 'intensity' .....
character(len=*)                  , intent(in)  :: is_space    ! The data are in 'reciprocal', 'patterson', 'direct' space
character(len=*)                  , intent(in)  :: radiation   ! 'xray', 'neutron', 'electron'
character(len=*)                  , intent(in)  :: space_group_name ! Hermann-Mauguin Symbol 
integer                           , intent(in)  :: symmetry_applied ! Data conform to symmetry Yes=1 / no = 0
integer                           , intent(in)  :: symmetry_n_mat   ! Number of symmetry operations
real(kind=PREC_DP), dimension(3,4,symmetry_n_mat), intent(in) :: symmetry_mat ! Actual symmetry matrices (real space)
!                                                                (x)'   (11, 12, 13, 14)   (x)
!                                                                (y)' = (21, 22, 23, 24) * (y)
!                                                                (z)'   (31, 32, 33, 34)   (z)
!                                                                (1)'   ( 0,  0,  0,  1)   (1)
!                                                                                         (11, 12, 13)
!                                                                (h, k, l)' = (h, k, l) * (21, 22, 23)
!                                                                                         (31, 32, 33)
integer                           , intent(out) :: ier_num     ! Error number
!
type(tuple)     :: p_args              ! Tuple of arguments for write_diffuse_scattering
type(object )   :: p_outfile           ! Output filename in python interface
type(object )   :: p_program           ! Output program  in python interface
type(object )   :: p_author            ! Output author   in python interface
type(object )   :: p_qvalues_dim       ! Data dimensionality in python interface
type(object )   :: p_radiation         ! radiation type 
type(object )   :: p_space             ! reciprocal, direct space
type(object )   :: p_content           ! intensity, 3d-delta-PDF, etc
type(object )   :: p_axes              ! axes labels ["h", "k","l"] etc
type(ndarray)   :: p_qvalues           ! Intensities     in python interface
type(ndarray)   :: p_lower_limits      ! Left_lower_bottom corner in python  interface
type(ndarray)   :: p_step_vectors      ! Increment vi along abscissa python  interface
type(ndarray)   :: p_h_indices         ! H indices       in python interface
type(ndarray)   :: p_k_indices         ! K indices       in python interface
type(ndarray)   :: p_l_indices         ! L indices       in python interface
type(ndarray)   :: p_unit_cell_length  ! Unit cell length in python interface
type(ndarray)   :: p_unit_cell_angle   ! Unit cell angles in python interface
type(object )   :: p_space_group       ! Space group name in python interface
type(object )   :: p_symmetry_applied  ! Space group symmetry has been applied to data 
type(object )   :: p_symmetry_n_mat    ! Number of symmetry matrices
type(ndarray)   :: p_symmetry_mat      ! Space group symmetry matrices 
!
type(module_py) :: interface_module    ! python script name
type(list)      :: paths_to_module     ! python script path
type(object)    :: return_value        ! forpy return value
!
character(len=1), dimension(3), parameter :: chkl = (/'h', 'k', 'l'/)
character(len=1), dimension(3), parameter :: cxyz = (/'x', 'y', 'z'/)
character(len=1), dimension(3), parameter :: cuvw = (/'u', 'v', 'w'/)
character(len=PREC_STRING)    :: axes ! ==  '["h", "k", "l"]' or similar permutation
integer :: i, j, k   ! Dummy indices
integer :: qvalues_dim  ! Data dimensionality
!real(KIND=PREC_DP), dimension(:,:,:), allocatable :: values        ! Values in C sequence
real(kind=PREC_DP), dimension(out_inc(abs_is_hkl)) :: h_indices      ! Actual H indices along a*
real(kind=PREC_DP), dimension(out_inc(ord_is_hkl)) :: k_indices      ! Actual K indices along b*
real(kind=PREC_DP), dimension(out_inc(top_is_hkl)) :: l_indices      ! Actual L indices along c*
!
write(*,*) ' WRITE NX '
!
!  Determine data dimensionality^
qvalues_dim = 0     ! Assume zero dimensional data 
do i=1, 3
   if(out_inc(i)>1) qvalues_dim = qvalues_dim + 1  ! Increment for each dimension
enddo
!
!  Set proper axes names at abscissa, ordinate, top_axis
!
if(is_space=='reciprocal') then
   axes      = '["h", "k", "l"]'
   axes( 3: 3) = chkl(abs_is_hkl)    ! Place the correct 'h', 'k', or 'l'
   axes( 8: 8) = chkl(ord_is_hkl)
   axes(13:13) = chkl(top_is_hkl)
elseif(is_space=='direct') then
   axes      = '["x", "y", "z"]'
   axes( 3: 3) = cxyz(abs_is_hkl)    ! Place the correct 'x', 'y', or 'z'
   axes( 8: 8) = cxyz(ord_is_hkl)
   axes(13:13) = cxyz(top_is_hkl)
elseif(is_space=='patterson') then
   axes      = '["u", "v", "w"]'
   axes( 3: 3) = cxyz(abs_is_hkl)    ! Place the correct 'u', 'v', or 'w'
   axes( 8: 8) = cxyz(ord_is_hkl)
   axes(13:13) = cxyz(top_is_hkl)
endif
!
write(*,*) ' Wrote original   array 321 ', qvalues(3,2,1)
write(*,*)
write(*,*) ' Wrote original   array llb ', qvalues(1         ,1         ,1         )
write(*,*) ' Wrote original   array rlb ', qvalues(out_inc(1),1         ,1         )
write(*,*) ' Wrote original   array lub ', qvalues(1         ,out_inc(2),1         )
write(*,*) ' Wrote original   array rub ', qvalues(out_inc(1),out_inc(2),1         )
!
write(*,*) ' Wrote original   array llt ', qvalues(1         ,1         ,out_inc(3))
write(*,*) ' Wrote original   array rlt ', qvalues(out_inc(1),1         ,out_inc(3))
write(*,*) ' Wrote original   array lut ', qvalues(1         ,out_inc(2),out_inc(3))
write(*,*) ' Wrote original   array rut ', qvalues(out_inc(1),out_inc(2),out_inc(3))
!
! Write full indices H, K, L along axes
!
write(*,*) ' abs_is_hkl ', abs_is_hkl, ord_is_hkl, top_is_hkl
do i=1, out_inc(abs_is_hkl)
   h_indices(i) = out_eck(1,1) + (i-1)*out_vi(abs_is_hkl,1)
enddo
do i=1, out_inc(ord_is_hkl)
   k_indices(i) = out_eck(2,1) + (i-1)*out_vi(ord_is_hkl,2)
enddo
do i=1, out_inc(top_is_hkl)
   l_indices(i) = out_eck(3,1) + (i-1)*out_vi(top_is_hkl,3)
enddo
!
ier_num = 0
call forpy_start()
!
!
ier_num = cast(          p_outfile,     outfile(1:len_trim(outfile)) )
ier_num = cast(          p_program,     'DISCUS 6.17.02')
ier_num = cast(          p_author,      'rbn')
ier_num = cast(          p_qvalues_dim, qvalues_dim)
ier_num = ndarray_create(p_qvalues,     qvalues)
ier_num = cast(          p_radiation,   radiation(1:len_trim(radiation)) )
ier_num = cast(          p_space,       is_space(1:len_trim(is_space)) )
ier_num = cast(          p_content,     content(1:len_trim(content)) )
ier_num = ndarray_create(p_lower_limits,out_eck(:,1))
ier_num = ndarray_create(p_step_vectors,(out_vi))
ier_num = cast          (p_axes        ,axes(1:len_trim(axes)))
ier_num = ndarray_create(p_h_indices,   h_indices)
ier_num = ndarray_create(p_k_indices,   k_indices)
ier_num = ndarray_create(p_l_indices,   l_indices)
ier_num = ndarray_create(p_unit_cell_length, cr_a0)
ier_num = ndarray_create(p_unit_cell_angle,  cr_win)
ier_num = cast(          p_space_group, space_group_name)
ier_num = cast(          p_symmetry_applied, symmetry_applied)
ier_num = cast(          p_symmetry_n_mat  , symmetry_n_mat  )
ier_num = ndarray_create(p_symmetry_mat    , symmetry_mat    )
!
!  Collect the 20 arguments into a tuple: p_args
!
ier_num = tuple_create(p_args, 20)
ier_num = p_args%setitem( 0, p_outfile)
ier_num = p_args%setitem( 1, p_program)
ier_num = p_args%setitem( 2, p_author)
ier_num = p_args%setitem( 3, p_qvalues_dim)
ier_num = p_args%setitem( 4, p_qvalues)
ier_num = p_args%setitem( 5, p_radiation)
ier_num = p_args%setitem( 6, p_space)
ier_num = p_args%setitem( 7, p_content)
ier_num = p_args%setitem( 8, p_lower_limits)
ier_num = p_args%setitem( 9, p_step_vectors)
ier_num = p_args%setitem(10, p_axes)
ier_num = p_args%setitem(11, p_h_indices)
ier_num = p_args%setitem(12, p_k_indices)
ier_num = p_args%setitem(13, p_l_indices)
ier_num = p_args%setitem(14, p_unit_cell_length)
ier_num = p_args%setitem(15, p_unit_cell_angle)
ier_num = p_args%setitem(16, p_space_group)
ier_num = p_args%setitem(17, p_symmetry_applied)
ier_num = p_args%setitem(18, p_symmetry_n_mat)
ier_num = p_args%setitem(19, p_symmetry_mat)
!ier_num = print_py(p_args)
!
! Append current directory to paths
!
ier_num = get_sys_path(paths_to_module)
ier_num = paths_to_module%append('.')
ier_num = import_py(interface_module, 'write_diffuse_scattering')
write(*,*) ' ier_num C ', ier_num
ier_num = call_py(return_value, interface_module, 'write_diffuse_scattering', p_args)
write(*,*) ' ier_num D ', ier_num
!
call p_outfile%destroy
call p_program%destroy
call p_author%destroy
call p_qvalues_dim%destroy
call p_qvalues%destroy
call p_radiation%destroy
call p_space%destroy
call p_content%destroy
call p_lower_limits%destroy
call p_step_vectors%destroy
call p_axes%destroy
call p_h_indices%destroy
call p_k_indices%destroy
call p_l_indices%destroy
call p_unit_cell_length%destroy
call p_unit_cell_angle%destroy
call p_space_group%destroy
call p_args%destroy
call p_symmetry_applied%destroy
call p_symmetry_n_mat%destroy
call p_symmetry_mat%destroy
call interface_module%destroy
!
end subroutine nx_write_scattering
!
!*******************************************************************************
end module nx_write_mod
