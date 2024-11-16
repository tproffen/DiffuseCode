module lib_nx_read_mod
!-
!  Load an Nexus File DiffuseDevelopers Common Format into the lib_f90 data structure
!+
!
private
public nx_read_scattering
!
contains
!
!*******************************************************************************
!
subroutine nx_read_scattering(infile, &
file_type, file_version, file_date, file_program, file_author, &
data_dim, data_shape, data, &
signal, radiation, is_space, content, &
lower_limits, step_vectors, axes, &
indices_abs, indices_ord, indices_top, &
unit_cell_lengths, unit_cell_angles, &
space_group, symmetry_applied, symmetry_n_mat, symmetry_mat,    &
ier_num  &
) 
!-
! Read a NEXUS HDF5 file in the Diffuse Developers common file format via FORPY
!+
!
use lib_forpython_mod
use forpy_mod
use precision_mod
!
use iso_fortran_env, only: real64
use iso_c_binding, only:C_CHAR
!
implicit none
!
character(LEN=1024)                       , intent(in)   :: infile
character(len=PREC_STRING)                , intent(out)  :: file_type   ! File type should be "Disorder scattering"
character(len=PREC_STRING)                , intent(out)  :: file_version! File version 
character(len=PREC_STRING)                , intent(out)  :: file_date   ! File date 
character(len=PREC_STRING)                , intent(out)  :: file_program! File program
character(len=PREC_STRING)                , intent(out)  :: file_author ! File author
integer                                   , intent(out)  :: data_dim          ! Data is 1, 2, 3 -dimensional
integer                   , dimension(3)  , intent(out)  :: data_shape        ! Actual dimensions along three axes
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(out) :: data ! Actual data array
character(len=PREC_STRING)                , intent(out)  :: signal            ! Data contain 'measurement', 'background' ,,, ,,,
character(len=PREC_STRING)                , intent(out)  :: radiation         ! Data were measured/calculated for this radiation
character(len=PREC_STRING)                , intent(out)  :: is_space          ! Data are 'reciprocal or 'direct'
character(len=PREC_STRING)                , intent(out)  :: content           ! Data contain this value type = 'Intensity', 3DPDF' ...
real(kind=PREC_DP)        , dimension(3)  , intent(out)  :: lower_limits      ! Left lower bottom corner
real(kind=PREC_DP)        , dimension(3,3), intent(out)  :: step_vectors      ! (:,1) along abscissa, (:,2) ordinate, (:,3) top axis
character(len=PREC_STRING)                , intent(out)  :: axes              ! string with ["h", "k","l"] or so
real(kind=PREC_DP)        , dimension(:    ), allocatable, intent(out) :: indices_abs ! Actual coordinates along abscissa
real(kind=PREC_DP)        , dimension(:    ), allocatable, intent(out) :: indices_ord ! Actual coordinates along ordinate
real(kind=PREC_DP)        , dimension(:    ), allocatable, intent(out) :: indices_top ! Actual coordinates along top axis
real(kind=PREC_DP)        , dimension(3)  , intent(out)  :: unit_cell_lengths ! Unit cell parameters (a, b, c)
real(kind=PREC_DP)        , dimension(3)  , intent(out)  :: unit_cell_angles  ! Unit cell parameters (alpha, beta, gamma)
character(len=PREC_STRING)                , intent(out)  :: space_group       ! Hermann-Mauguin Symbol
integer                                   , intent(out)  :: symmetry_applied  ! Data conform to symmetry
integer                                   , intent(out)  :: symmetry_n_mat    ! Data conform to symmetry
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(out) :: symmetry_mat ! Actual Symmetry matrices
integer                                   , intent(out)  :: ier_num           ! Error number
!
character(len=:), allocatable :: return_string
!
type(object )   :: p_infile           ! Output filename in python interface
!
type(tuple)     :: p_args             ! Tuple of arguments for  read_diffuse_scattering
type(module_py) :: interface_module   ! python script name
type(list)      :: paths_to_module    ! python script path
type(object)    :: return_value       ! forpy return value
type(tuple)     :: returned_tuple     ! All results as tuple
type(object)    :: temp, temp2        ! temporary python object 
type(tuple)     :: temp_tuple         ! temporary python tuple
type(ndarray)   :: temp_arr           ! temporary numpy array
!
real(kind=real64     ), dimension(:,:,:), pointer :: matrix_3d  ! Pointer to result of get_data
real(kind=real64     ), dimension(:,:)  , pointer :: matrix_2d  ! Pointer to result of get_data
real(kind=real64     ), dimension(:)    , pointer :: matrix_1d  ! Pointer to result of get_data
!
integer :: iarg
integer :: ih,ik,il   ! Dummy loop arrays
integer :: i ,j ,k    ! Dummy loop arrays
!
! TEMPORARY DEBUG
!character(kind=C_CHAR, len=:), allocatable :: dname
!
!
call forpy_start()                   ! If needed, start forpy
!
ier_num = cast(p_infile, infile(1:len_trim(infile)) )
ier_num = tuple_create(p_args, 1)
ier_num = p_args%setitem( 0, p_infile)
!
! Append current directory to paths
!
ier_num = get_sys_path(paths_to_module)
ier_num = paths_to_module%append('.')!
ier_num = import_py(interface_module, 'read_diffuse_scattering')
ier_num = call_py(return_value, interface_module, 'read_diffuse_scattering', p_args)
ier_num = cast(returned_tuple, return_value)   ! Get multiple returned objects
!
!  Get file type 
!
iarg = 0
ier_num = returned_tuple%getitem(temp, iarg)      ! getitem does not allow getitem(p_arr, o)
ier_num = cast(return_string, temp)
file_type = return_string(1:len_trim(return_string))
call temp%destroy
!
!  Get file version
!
iarg = 1
ier_num = returned_tuple%getitem(temp, iarg)      ! getitem does not allow getitem(p_arr, o)
ier_num = cast(return_string, temp)
file_version = return_string(1:len_trim(return_string))
call temp%destroy
!
!  Get file date
!
iarg = 2
ier_num = returned_tuple%getitem(temp, iarg)      ! getitem does not allow getitem(p_arr, o)
ier_num = cast(return_string, temp)
file_date = return_string(1:len_trim(return_string))
call temp%destroy
!
!  Get file program
!
iarg = 3
ier_num = returned_tuple%getitem(temp, iarg)      ! getitem does not allow getitem(p_arr, o)
ier_num = cast(return_string, temp)
file_program = return_string(1:len_trim(return_string))
call temp%destroy
!
!  Get file author
!
iarg = 4
ier_num = returned_tuple%getitem(temp, iarg)      ! getitem does not allow getitem(p_arr, o)
ier_num = cast(return_string, temp)
file_author = return_string(1:len_trim(return_string))
call temp%destroy
!
!  Get data dimension
!
iarg = 5
ier_num = returned_tuple%getitem(temp, iarg)      !
ier_num = cast(data_dim, temp)
call temp%destroy
!
!  Get actual dimensions
!
iarg = 6
ier_num = returned_tuple%getitem(temp, iarg)      !
ier_num = cast(temp_tuple, temp)
call temp%destroy
ier_num = temp_tuple%getitem(temp2, 0)      !
ier_num = cast(data_shape(1), temp2)
call temp2%destroy
ier_num = temp_tuple%getitem(temp2, 1)      !
ier_num = cast(data_shape(2), temp2)
call temp2%destroy
ier_num = temp_tuple%getitem(temp2, 2)      !
ier_num = cast(data_shape(3), temp2)
call temp2%destroy
!
! Get actual data array
!
iarg = 7
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_3d,'C')
allocate(data(data_shape(1), data_shape(2), data_shape(3)))
do ih = 1, data_shape(1)
do ik = 1, data_shape(2)
do il = 1, data_shape(3)
   data(ih, ik, il) = matrix_3d(il,ik,ih)
enddo
enddo
enddo
nullify(matrix_3d)
call temp%destroy
call temp_arr%destroy
!
!write(*,*) ' Read  original   array 321 ', d5_data(3,2,1)
!write(*,*)
!write(*,*) ' Read  original   array llb ', d5_data(1         ,1         ,1         )
!write(*,*) ' Read  original   array rlb ', d5_data(data_shape(1),1         ,1         )
!write(*,*) ' Read  original   array lub ', d5_data(1         ,data_shape(2),1         )
!write(*,*) ' Read  original   array rub ', d5_data(data_shape(1),data_shape(2),1         )
!
!write(*,*) ' Read  original   array llt ', d5_data(1         ,1         ,data_shape(3))
!write(*,*) ' Read  original   array rlt ', d5_data(data_shape(1),1         ,data_shape(3))
!write(*,*) ' Read  original   array lut ', d5_data(1         ,data_shape(2),data_shape(3))
!write(*,*) ' Read  original   array rut ', d5_data(data_shape(1),data_shape(2),data_shape(3))
!call get_data(iarg, returned_tuple, ergebnis, 'C')
!
! get signal
!
iarg = 8
ier_num = returned_tuple%getitem(temp, iarg)     ! Get item into temporary object
ier_num = cast(return_string, temp)
signal = return_string(1:len_trim(return_string))
call temp%destroy
!
! get radiation
!
iarg = 9
ier_num = returned_tuple%getitem(temp, iarg)     ! Get item into temporary object
ier_num = cast(return_string, temp)
radiation = return_string(1:len_trim(return_string))
call temp%destroy
!
! get space
!
iarg = 10
ier_num = returned_tuple%getitem(temp, iarg)     ! Get item into temporary object
ier_num = print_py(temp)
ier_num = cast(return_string, temp)
is_space = return_string(1:len_trim(return_string))
call temp%destroy
!
! get content == value type 
!
iarg = 11
ier_num = returned_tuple%getitem(temp, iarg)     ! Get item into temporary object
ier_num = print_py(temp)
ier_num = cast(return_string, temp)
content = return_string(1:len_trim(return_string))
call temp%destroy
!
! Get lower limits
!
iarg = 12
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
lower_limits = matrix_1d
deallocate(matrix_1d)
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
! Get step vectors
!
iarg = 13
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_2d,'C')
step_vectors = transpose(matrix_2d)
call temp%destroy
call temp_arr%destroy
!
!  Get axes string
!
iarg = 14
ier_num = returned_tuple%getitem(temp, iarg)      ! getitem does not allow getitem(p_arr, o)
ier_num = cast(return_string, temp)
axes = return_string(1:len_trim(return_string))
call temp%destroy
!
! Get indices along abscissa
!
iarg = 15
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
allocate(indices_abs(ubound(matrix_1d,1)))
indices_abs = matrix_1d
!deallocate(matrix_1d)
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
! Get indices along ordinate
!
iarg = 16
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
allocate(indices_ord(ubound(matrix_1d,1)))
indices_ord = matrix_1d
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
! Get indices along top_axis
!
iarg = 17
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
allocate(indices_top(ubound(matrix_1d,1)))
indices_top = matrix_1d
call temp%destroy
call temp_arr%destroy
nullify(matrix_1d)
!
! Get unit cell parameters
!
iarg = 18
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
unit_cell_lengths = matrix_1d
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
iarg = 19
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
unit_cell_angles = matrix_1d
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
! Get space group symbol
!
iarg = 20
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(return_string, temp)
space_group = return_string(1:len_trim(return_string))
call temp%destroy
!
! Get information if symmetry has been applied
!
iarg = 21
ier_num = returned_tuple%getitem(temp, iarg)      !
ier_num = cast(symmetry_applied, temp)
call temp%destroy
!
! Get number of symmetry matrices
!
iarg = 22
ier_num = returned_tuple%getitem(temp, iarg)      !
ier_num = cast(symmetry_n_mat, temp)
call temp%destroy
!
! Get symmetry matrices
!
iarg = 23
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_3d,'C')
allocate(symmetry_mat(3, 3, symmetry_n_mat))
symmetry_mat = 0.0_PREC_DP
do ih = 1, 3
   do ik = 1, 3
      do il = 1, symmetry_n_mat
         symmetry_mat(ih, ik, il) = matrix_3d(il,ik,ih)
      enddo
   enddo
enddo
nullify(matrix_3d)
!
call temp%destroy
call temp_arr%destroy
!
end subroutine nx_read_scattering
!
!*******************************************************************************
!
end module lib_nx_read_mod
