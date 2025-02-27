module lib_nx_read_mod
!-
!  Load an Nexus File DiffuseDevelopers Common Format into the lib_f90 data structure
!+
!
private
public nx_read_scattering
public nx_read_structure
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
!integer :: i ,j ,k    ! Dummy loop arrays
!
! TEMPORARY DEBUG
!character(kind=C_CHAR, len=:), allocatable :: dname
!
!
call forpy_start(ier_num)                   ! If needed, start forpy
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
subroutine nx_read_structure(python_script_dir, infile, unit_cell_lengths, unit_cell_angles,    &
                             metric_tensor,  &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat, unit_cells, number_of_types, types_names,      &
                             types_ordinal, types_charge, types_isotope, number_of_atoms, &
                             atom_ID, atom_type, atom_pos, atom_unit_cell, atom_site,     &
                             atom_property, crystal_flags, crystal_meta,                  &
                             anisotropic_adp, molecules, average_struc, types_occupancy,  &
                             ier_num)
!-
! Read a NEXUS HDF5 structure file in the Diffuse Developers common file format via FORPY
!+
!
use lib_forpython_mod
use lib_nx_transfer_mod
use forpy_mod
use precision_mod
!
use iso_fortran_env, only: real64
use iso_c_binding, only:C_CHAR
!
implicit none
!
character(len=*)                                         , intent(in)  :: python_script_dir  ! Full path to directory with python script
                                                                     ! '/home/reinhard/.local/share/diffuse_scattering'
character(len=*)                                         , intent(in)  :: infile
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
real(kind=PREC_DP)        , dimension(3,3)               , intent(out) :: metric_tensor      ! Direct space metric tensor
character(len=32)                                        , intent(out) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(out) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(out) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(out) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(out) :: symmetry_mat       ! Actual Symmetry matrices
integer                   , dimension(3,3  )             , intent(out) :: unit_cells         ! Number of unit cells
!
integer                                                  , intent(out) :: number_of_types    ! Crystal has this many atom types
character(len=4)          , dimension(:),     allocatable, intent(out) :: types_names        ! Name as "O", "O2-" etc
integer                   , dimension(:),     allocatable, intent(out) :: types_ordinal      ! Ordinal number of chemical element
integer                   , dimension(:),     allocatable, intent(out) :: types_charge       ! Atom types have this charge
integer                   , dimension(:),     allocatable, intent(out) :: types_isotope      ! Atom type is this isotope or zero
integer                                                  , intent(out) :: number_of_atoms    ! Crystal contains this many actual atoms
integer                   , dimension(:    ), allocatable, intent(out) :: atom_ID            ! Atom ID 
integer                   , dimension(:    ), allocatable, intent(out) :: atom_type          ! Atom is of this type
real(kind=PREC_DP)        , dimension(:,:  ), allocatable, intent(out) :: atom_pos           ! Atom is at these fractional coordinates
integer                   , dimension(:,:  ), allocatable, intent(out) :: atom_unit_cell     ! Atom is in this unit cell
integer                   , dimension(:    ), allocatable, intent(out) :: atom_site          ! Atom is on this site in its unit cell
integer                   , dimension(:    ), allocatable, intent(out) :: atom_property      ! Atom has this property flag
logical                   , dimension(2,6)               , intent(out) :: crystal_flags      ! Flags , see "c_flags"
character(len=PREC_STRING), dimension(  5)               , intent(out) :: crystal_meta       ! Metadata, see "c_meta"
type(anis_adp_type)                                      , intent(out) :: anisotropic_adp    ! Info on anisotropic ADP ==> lib_nx_transfer.f90
type(molecule_data)                                      , intent(out) :: molecules          ! Info on molecules       ==> lib_nx_transfer.f90
type(average_structure)                                  , intent(out) :: average_struc      ! Info on average struct  ==> lib_nx_transfer.f90
real(kind=PREC_DP)        , dimension(:),     allocatable, intent(out) :: types_occupancy    ! This atoms type has an occupancy of value
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
!
character(len=:), allocatable :: return_string
!
character(len=24), dimension(6) :: c_flags
character(len=26), dimension(5) :: c_meta
integer :: i, j
integer :: ih, ik, il
integer :: iarg
logical :: lflag
!
character(len=:)                    , allocatable :: str_value
integer               , dimension(:  )  , pointer :: matrix_1d_i  ! Pointer to integer result of get_data
integer               , dimension(:,:)  , pointer :: matrix_2d_i  ! Pointer to integer result of get_data
integer               , dimension(:,:)  , pointer :: matrix_2d_it ! Pointer to integer result of get_data
!integer               , dimension(:,:,:), pointer :: matrix_3d_i  ! Pointer to integer result of get_data
real(kind=real64     ), dimension(:,:,:), pointer :: matrix_3d  ! Pointer to result of get_data
real(kind=real64     ), dimension(:,:)  , pointer :: matrix_2d  ! Pointer to result of get_data
real(kind=real64     ), dimension(:,:)  , pointer :: matrix_2d_t! Pointer to result of get_data
real(kind=real64     ), dimension(:)    , pointer :: matrix_1d  ! Pointer to result of get_data
!
type(object )   :: p_infile           ! Output filename in python interface
type(tuple)     :: p_args             ! Tuple of arguments for  read_diffuse_scattering
type(module_py) :: interface_module   ! python script name
type(list)      :: paths_to_module    ! python script path
type(object)    :: return_value       ! forpy return value
type(tuple)     :: returned_tuple     ! All results as tuple
type(object)    :: temp, temp2, temp3 ! temporary python object 
type(tuple )    :: temp_tuple         ! temporary python tuple object 
type(list  )    :: temp_list          ! temporary python list object 
type(dict  )    :: temp_dict          ! temporary python dircionary object 
type(ndarray)   :: temp_arr           ! temporary numpy array
!
data c_flags /'is_super_structure      ', &
              'is_asymmetric_unit      ', &
              'is_periodic_x           ', &
              'is_periodic_y           ', &
              'is_periodic_z           ', &
              'is_homogeneous          '  &
             /
!
data c_meta  /'audit_author_name         ', &
              'audit_conform_dict_name   ', &
              'audit_conform_dict_version', &
              'audit_creation_date       ', &
              'audit_creation_method     '  &
             /
!
!
call forpy_start(ier_num)                   ! If needed, start forpy
!
!write(*,*) ' INFILE >',infile(1:len_trim(infile)),'<' 
ier_num = cast(p_infile, infile(1:len_trim(infile)) )
call err_print()
ier_num = tuple_create(p_args, 1)
call err_print()
ier_num = p_args%setitem( 0, p_infile)
call err_print()
!ier_num = print_py(p_args)
!call err_print()
!
! Append current directory to paths
!
ier_num = get_sys_path(paths_to_module)
call err_print()
ier_num = paths_to_module%append(python_script_dir(1:len_trim(python_script_dir))) !'.')!
call err_print()
ier_num = import_py(interface_module, 'nexus_structure')
call err_print()
!write(*,*) ' EXECUTE MACRO '
ier_num = call_py(return_value, interface_module, 'read_diffuse_structure', p_args)
call err_print()
!write(*,*) ' CAST RETURN VALUES'
ier_num = cast(returned_tuple, return_value)   ! Get multiple returned objects
call err_print()
!read(*,*) iarg
!ier_num= -1
!
! Get unit_cell_lengths
!
iarg =  0
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
unit_cell_lengths = matrix_1d
deallocate(matrix_1d)
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
! Get unit_cell_angles
!
iarg =  1
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_1d,'C')
unit_cell_angles = matrix_1d
deallocate(matrix_1d)
nullify(matrix_1d)
call temp%destroy
call temp_arr%destroy
!
! Get metric tensor
!
iarg = 2
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_2d,'C')        ! Copy into Fortran 3x3 matrix
metric_tensor = matrix_2d
deallocate(matrix_2d)
nullify(matrix_2d)
call temp%destroy
call temp_arr%destroy
!
! Get Hermann Mauguin Symbol
!
iarg =  3
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(return_string, temp)
symmetry_H_M = return_string(1:len_trim(return_string))
call temp%destroy
!
! Get Space group origin
!
iarg =  4
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(symmetry_origin, temp)
call temp%destroy
!
! Get Space group permutation 'abc'
!
iarg =  5
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(return_string, temp)
symmetry_abc = return_string(1:len_trim(return_string))
call temp%destroy
!
! Get number of symmetry matrices
!
iarg =  6
ier_num = returned_tuple%getitem(temp, iarg)      !
ier_num = cast(symmetry_n_mat, temp)
call temp%destroy
!
! Get symmetry matrices
!
iarg =  7
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_3d,'C')
allocate(symmetry_mat(4, 4, symmetry_n_mat))
symmetry_mat = 0.0_PREC_DP
symmetry_mat(1,1,:) = 1.0_PREC_DP
do ih = 1, 3
   do ik = 1, 4
      do il = 1, symmetry_n_mat
         symmetry_mat(ih, ik, il) = matrix_3d(il,ik,ih)
      enddo
   enddo
enddo
deallocate(matrix_3d)
nullify(matrix_3d)
call temp%destroy
call temp_arr%destroy
!
! Get Number of unit cells
!
iarg =  8
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
ier_num = temp_arr%get_data(matrix_2d_i,'C')
unit_cells = matrix_2d_i
call temp%destroy
call temp_arr%destroy
deallocate(matrix_2d_i)
nullify(matrix_2d_i)
!
! Get number of atom types
!
iarg =  9
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(number_of_types, temp)
call temp%destroy
!
! Get and interpret atom_type structure
!
iarg = 10
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_list, temp)                   ! Recast the temporary object into a list
!
allocate(types_names(number_of_types))
allocate(types_ordinal(number_of_types))
allocate(types_charge (number_of_types))
allocate(types_isotope(number_of_types))
types_names   = ' '
types_ordinal = 0
do i=0,number_of_types-1
   ier_num = temp_list%getitem(temp2, i)          ! Get item for current atom type
   ier_num = cast(temp_tuple, temp2)              ! Recast the temporary object into a list
!
   ier_num = temp_tuple%getitem(temp3, 0)         ! Get atom name
   ier_num = cast(str_value, temp3)               ! Cast atom naem into fortran string
   types_names(i+1) = str_value(1:len_trim(str_value))  ! Copy atom_name into Fortran array
   deallocate(str_value)
   call temp3%destroy
!
   ier_num = temp_tuple%getitem(temp3, 1)         ! Get atom ordinal number
   ier_num = cast(j, temp3)                       ! Cast into a number
   types_ordinal(i+1) = j                         !Copy ordinal number into Fortran array
   call temp3%destroy
!
   ier_num = temp_tuple%getitem(temp3, 2)         ! Get atom ordinal number
   ier_num = cast(j, temp3)                       ! Cast into a number
   types_charge (i+1) = j                         ! Copy charge 
   call temp3%destroy
!
   ier_num = temp_tuple%getitem(temp3, 3)         ! Get atom ordinal number
   ier_num = cast(j, temp3)                       ! Cast into a number
   types_isotope(i+1) = j                         !Copy ordinal number into Fortran array
   call temp3%destroy
!
   call temp2%destroy
   call temp_tuple%destroy
enddo
call temp_list%destroy
call temp%destroy
!
! Get number of atoms
!
iarg = 11
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(number_of_atoms, temp)
call temp%destroy
!
! Get and interpret atom_data structure; integer part
!
iarg = 12
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
if(number_of_atoms==1) then                       ! Single atom, need to perform get_data without 'C'
   ier_num = temp_arr%get_data(matrix_2d_it   )
   allocate(matrix_2d_i(number_of_atoms,6))
   !write(*,*) ' GOT ATOM DATA F ', ier_num, shape(matrix_2d_it), shape(matrix_2d_i)
   matrix_2d_i = transpose(matrix_2d_it)
   deallocate(matrix_2d_it)
   nullify(matrix_2d_it)
else
   ier_num = temp_arr%get_data(matrix_2d_i,'C')
   !write(*,*) ' GOT ATOM DATA C ', ier_num, shape(matrix_2d_i)
endif
allocate(atom_ID       (  number_of_atoms))
allocate(atom_type     (  number_of_atoms))
allocate(atom_unit_cell(3,number_of_atoms))
allocate(atom_site     (  number_of_atoms))
! 
do i=1,number_of_atoms
   atom_ID       (  i) = matrix_2d_i(i,1)
   atom_type     (  i) = matrix_2d_i(i,2)
   atom_unit_cell(:,i) = matrix_2d_i(i,3:5)
   atom_site     (  i) = matrix_2d_i(i,6)
enddo
call temp%destroy
call temp_arr%destroy
deallocate(matrix_2d_i)
nullify(matrix_2d_i)
!
! Get and interpret atom_data structure; float part == positions
!
iarg = 13
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
ier_num = cast(temp_arr, temp)                    ! Now copy temp into numpy array
if(number_of_atoms==1) then                       ! Single atom, need to perform get_data without 'C'
   ier_num = temp_arr%get_data(matrix_2d_t  )
   allocate(matrix_2d(number_of_atoms,3))
   !write(*,*) ' GOT ATOM POS  F ', ier_num, shape(matrix_2d_t), shape(matrix_2d)
   matrix_2d = transpose(matrix_2d_t)
   deallocate(matrix_2d_t)
   nullify(matrix_2d_t)
else
   ier_num = temp_arr%get_data(matrix_2d,'C')
   !write(*,*) ' GOT ATOM POS  C ', ier_num, shape(matrix_2d)
endif
allocate(atom_pos      (3,number_of_atoms))
!  unit_cells = matrix_2d_i
do i=1,number_of_atoms
   atom_pos      (:,i) = matrix_2d(i,:)
enddo
call temp%destroy
call temp_arr%destroy
deallocate(matrix_2d)
nullify(matrix_2d)
!
! Get Flags
!
iarg = 14
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
crystal_flags = .false.
if(.not.is_none(temp)) then
   ier_num = cast(temp_dict, temp)
   do i=0,5
      ier_num = temp_dict%getitem(lflag, c_flags(i+1)(1:len_trim(c_flags(i+1))))
      if(ier_num==0) then
         crystal_flags(1, i+1) = lflag
         crystal_flags(2, i+1) = .true.
      else
         crystal_flags(1, i+1) = .false.
         crystal_flags(2, i+1) = .false.
      endif
   enddo
   call temp_dict%destroy
!else
!   write(*,*) ' Status flags are None'
endif
call temp%destroy
!
! Get Metadata
!
iarg = 15
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
crystal_meta  = ' '
if(.not.is_none(temp)) then
   ier_num = cast(temp_dict, temp)
   do i=0,4
      ier_num = temp_dict%getitem(str_value, c_meta(i+1)(1:len_trim(c_meta(i+1))))
      if(ier_num==0) then
        crystal_meta(i+1)  = str_value(1:len_trim(str_value))
      else
        crystal_meta(i+1)  = ' '
      endif
      deallocate(str_value)
   enddo
   call temp_dict%destroy
endif
call temp%destroy
!
! Get and interpret average structure
!
iarg = 16
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
if(is_none(temp)) then
   !write(*,'(a)') ' Average Structure is NONE'
   continue
else
   continue
   !write(*,'(a)') ' Average Structure exists'
   ier_num = cast(temp_tuple, temp)               ! Now copy temp into temporary tuple
!
   ier_num = temp_tuple%getitem(temp2, 0)         ! Get item into temporary objec
   ier_num = cast(average_struc%aver_n_atoms, temp2)
   call temp2%destroy
   !write(*,*) ' AVERAGE contains ', average_struc%aver_n_atoms
!
   if(average_struc%aver_n_atoms>0) then                ! We have atoms in the average structure
      allocate(average_struc%atom_type  (   average_struc%aver_n_atoms))
      allocate(average_struc%position   (3, average_struc%aver_n_atoms))
      allocate(average_struc%occupancy  (   average_struc%aver_n_atoms))
      allocate(average_struc%anis_adp   (7, average_struc%aver_n_atoms))
      allocate(average_struc%site_number(   average_struc%aver_n_atoms))
!
      ier_num = temp_tuple%getitem(temp2, 1)         ! Get item into temporary objec
      ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
      ier_num = temp_arr%get_data(matrix_1d_i)
      average_struc%atom_type = matrix_1d_i          ! Get atom types
      deallocate(matrix_1d_i)
      nullify(matrix_1d_i)
      call temp2%destroy
!
      ier_num = temp_tuple%getitem(temp2, 2)         ! Get item into temporary object
      ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
      if(average_struc%aver_n_atoms==1) then         ! Just one atom, need to do get_item without 'C'
         ier_num = temp_arr%get_data(matrix_2d_t)
         allocate(matrix_2d(average_struc%aver_n_atoms,3))
      !write(*,*) ' GOT ADPs      F ', ier_num, shape(matrix_2d_t), shape(matrix_2d)
         matrix_2d = transpose(matrix_2d_t)
         deallocate(matrix_2d_t)
         nullify(matrix_2d_t)
      else
         ier_num = temp_arr%get_data(matrix_2d, 'C')
      !write(*,*) ' GOT ADPs      C ', ier_num, shape(matrix_2d)
      endif
      average_struc%position = transpose(matrix_2d)
      deallocate(matrix_2d)
      nullify(matrix_2d)
      call temp2%destroy
      call temp_arr%destroy
!
      ier_num = temp_tuple%getitem(temp2, 3)         ! Get item into temporary objec
      ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
      ier_num = temp_arr%get_data(matrix_1d)
      average_struc%occupancy = matrix_1d            ! Get atom types
      deallocate(matrix_1d)
      nullify(matrix_1d)
      call temp2%destroy
!
      ier_num = temp_tuple%getitem(temp2, 4)         ! Get item into temporary object
      ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
      if(average_struc%aver_n_atoms==1) then         ! Just one atom, need to do get_item without 'C'
         ier_num = temp_arr%get_data(matrix_2d_t)
         allocate(matrix_2d(average_struc%aver_n_atoms,7))
      !write(*,*) ' GOT ADPs      F ', ier_num, shape(matrix_2d_t), shape(matrix_2d)
         matrix_2d = transpose(matrix_2d_t)
         deallocate(matrix_2d_t)
         nullify(matrix_2d_t)
      else
         ier_num = temp_arr%get_data(matrix_2d, 'C')
      !write(*,*) ' GOT ADPs      C ', ier_num, shape(matrix_2d)
      endif
      average_struc%anis_adp = transpose(matrix_2d)
      deallocate(matrix_2d)
      nullify(matrix_2d)
      call temp2%destroy
      call temp_arr%destroy
!
      ier_num = temp_tuple%getitem(temp2, 5)         ! Get item into temporary objec
      ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
      ier_num = temp_arr%get_data(matrix_1d_i)
      average_struc%site_number = matrix_1d_i        ! Get atom types
      deallocate(matrix_1d_i)
      nullify(matrix_1d_i)
      call temp2%destroy
   endif
endif
call temp%destroy
!
! Get and interpret property flags
!
iarg = 18
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
!write(*,*) ' NUMBER OF ATOMS ', number_of_atoms
allocate(atom_property (  number_of_atoms))
!write(*,*) ' NUMBER OF ATOMS ', number_of_atoms
atom_property = 1                                 ! Default to simplest property
if(is_none(temp)) then
   !write(*,'(a)') ' Property flags are not provided'
!   write(*,*)     ' PROPERTY ', ubound(atom_property)   
else
   ier_num = cast(temp_arr, temp)                 ! Now copy temp into numpy array
!   write(*,*) ' CAST PROPERTY to TEMP ', ier_num
   ier_num = temp_arr%get_data(matrix_1d_i) !,'C')
!   write(*,*) ' GOT DATA   FROM  TEMP ', ier_num
!   write(*,*) ' PROPERTY ', ubound(matrix_1d_i), ubound(atom_property), number_of_atoms
   atom_property = matrix_1d_i
   deallocate(matrix_1d_i)
   nullify(matrix_1d_i)
   call temp_arr%destroy
endif
call temp%destroy
!
! Get and interpret ADP parameters
!
iarg = 19
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
!write(*,*) ' GOT INTO TEMP ', ier_num
if(is_none(temp)) then
!   write(*,'(a)') ' ADP parameters are not provided '
   anisotropic_adp%anis_n_type = 0
else
   ier_num = cast(temp_tuple, temp)               ! Now copy temp into temporary tuple
!
   ier_num = temp_tuple%getitem(temp2, 0)         ! Get item into temporary objec
   ier_num = cast(i, temp2)
   call temp2%destroy
!
   ier_num = temp_tuple%getitem(temp2, 1)         ! Get item into temporary objec
   ier_num = cast(anisotropic_adp%anis_n_type, temp2)
   call temp2%destroy
!write(*,*) ' Number of UIJ ', i, anisotropic_adp%anis_n_type
!
   ier_num = temp_tuple%getitem(temp2, 2)         ! Get item into temporary objec
   ier_num = cast(anisotropic_adp%anis_n_atom, temp2)
   call temp2%destroy
!write(*,*) ' Number of ADPs', anisotropic_adp%anis_n_atom
!
   ier_num = temp_tuple%getitem(temp2, 3)         ! Get item into temporary object
   ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
   if(anisotropic_adp%anis_n_type==1) then        ! Just one type, need to do get_item without 'C'
      ier_num = temp_arr%get_data(matrix_2d_t)
      allocate(matrix_2d(anisotropic_adp%anis_n_type,7))
   !write(*,*) ' GOT ADPs      F ', ier_num, shape(matrix_2d_t), shape(matrix_2d)
      matrix_2d = transpose(matrix_2d_t)
      deallocate(matrix_2d_t)
      nullify(matrix_2d_t)
   else
      ier_num = temp_arr%get_data(matrix_2d, 'C')
   !write(*,*) ' GOT ADPs      C ', ier_num, shape(matrix_2d)
   endif
   anisotropic_adp%anis_n_type = ubound(matrix_2d,1)
   allocate(anisotropic_adp%anis_adp(7, anisotropic_adp%anis_n_type ))
   anisotropic_adp%anis_adp = transpose(matrix_2d)
   deallocate(matrix_2d)
   nullify(matrix_2d)
   call temp2%destroy
   call temp_arr%destroy
!
   ier_num = temp_tuple%getitem(temp2, 4)         ! Get item into temporary object
   ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
   ier_num = temp_arr%get_data(matrix_1d_i)
!  call      err_print()
   allocate(anisotropic_adp%atom_index(number_of_atoms))
   anisotropic_adp%atom_index = matrix_1d_i
   deallocate(matrix_1d_i)
   nullify(matrix_1d_i)
   call temp2%destroy
   call temp_arr%destroy
   call temp_tuple%destroy
endif
call temp%destroy
!
! Get and interpret ADP index list
!
!iarg = 20
!ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
!if(is_none(temp)) then
!   write(*,'(a)') ' ADP index list is not provided'
!else
!   ier_num = cast(temp_arr, temp)                 ! Now copy temp into numpy array
!   ier_num = temp_arr%get_data(matrix_1d_i)
!   allocate(atom_adp_index(number_of_atoms))
!   write(*,*) ' ADP_INDEX dimension ', ubound(matrix_1d_i,1), number_of_atoms
!   atom_adp_index = matrix_1d_i
!   deallocate(matrix_1d_i)
!   nullify(matrix_1d_i)
!endif
!
! Get and interpret Molecular info
!
iarg = 20
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
if(is_none(temp)) then
   molecules%number_moles = 0
   molecules%number_types = 0
else
   ier_num = cast(temp_tuple, temp)               ! Molecular info is contained inside a tuple 
!
!  Determine molecule number
!
   ier_num = temp_tuple%getitem(temp2, 0)         ! Get item into temporary object
   ier_num = cast(i, temp2)
   molecules%number_moles = i
   call temp2%destroy
!write(*,*) ' READ MOLE_NUMBER ', molecules%number_moles
!
!  Determine molecule type number
!
   ier_num = temp_tuple%getitem(temp2, 1)         ! Get item into temporary object
   ier_num = cast(i, temp2)
   molecules%number_types = i
   call temp2%destroy
!write(*,*) ' READ MOLE_TYPES  ', molecules%number_types
!
!  Get integer info for each molecule
!
   ier_num = temp_tuple%getitem(temp2, 2)         ! Get item into temporary object
call err_print()
   ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
call err_print()
   if(molecules%number_moles==1) then             ! Just one molecule
      ier_num = temp_arr%get_data(matrix_2d_i)    ! A (3,1) Matrix needs to be read as default
call err_print()
!write(*,*) ' MATRIX f  2 ', ier_num, ' SHAPE ', shape(matrix_2d_i)
      allocate(molecules%mole_int(3, molecules%number_moles))
      molecules%mole_int =          (matrix_2d_i)
   elseif(molecules%number_moles>1) then          ! More than one molecule
      ier_num = temp_arr%get_data(matrix_2d_i ,'C') ! A (3, N) matrix can be read as C-style
call err_print()
!write(*,*) ' MATRIX c  2 ', ier_num, ' SHAPE ', shape(matrix_2d_i)
      allocate(molecules%mole_int(3, molecules%number_moles))
      molecules%mole_int = transpose(matrix_2d_i)
   endif
   deallocate(matrix_2d_i)
   nullify(matrix_2d_i)
   call temp2%destroy
   call temp_arr%destroy
!write(*,*) ' INTEGER MATRIX ', shape(molecules%mole_int)
!
!  Get real valued info for each molecule type
!
   ier_num = temp_tuple%getitem(temp2, 3)         ! Get item into temporary object
   ier_num = cast(temp_arr, temp2)                ! Now copy temp into numpy array
   if(molecules%number_types==1) then              ! Just one molecule type
      ier_num = temp_arr%get_data(matrix_2d )      ! Cast as 2d matrix in default order
!write(*,*) ' MATRIX f  3 ', ier_num, ' SHAPE ', shape(matrix_2d)
      allocate(molecules%mole_real(3, molecules%number_types))
      molecules%mole_real =          (matrix_2d )
   elseif(molecules%number_types>1) then           ! More than one molecule type
      ier_num = temp_arr%get_data(matrix_2d ,'C')  ! Cast as 2d matrix in C order
!write(*,*) ' MATRIX c  3 ', ier_num, ' SHAPE ', shape(matrix_2d)
      allocate(molecules%mole_real(ubound(matrix_2d,1), ubound(matrix_2d,2)))
      molecules%mole_real = transpose(matrix_2d)
   endif
   deallocate(matrix_2d)
   nullify(matrix_2d)
   call temp2%destroy
   call temp_arr%destroy
!write(*,*) ' REAL    MATRIX ', shape(molecules%mole_real)
!
!  Get molecule content for each molecule
!
   ier_num = temp_tuple%getitem(temp2, 4)      ! Get item into temporary object
   ier_num = cast(temp_arr, temp2)             ! Now copy temp into numpy array
   j = maxval(molecules%mole_int(3,:))         ! Maximum molecule length
   if(molecules%number_moles==1 .or. j==1) then          ! Just one molecule or at most one atom per molecule
      ier_num = temp_arr%get_data(matrix_2d_i) ! ,'C')
!write(*,*) ' MATRIX f  4 ', ier_num, ' SHAPE ', shape(matrix_2d_i), shape(temp_arr)
      allocate(molecules%atom_index(ubound(matrix_2d_i,1), molecules%number_moles))
      molecules%atom_index =          (matrix_2d_i)
!  elseif(molecules%number_moles>1) then        ! More than one molecule
   else
      ier_num = temp_arr%get_data(matrix_2d_i ,'C')
!write(*,*) ' MATRIX c  4 ', ier_num, ' SHAPE ', shape(matrix_2d_i)
!write(*,*) ' MATRIX c  4 ', ier_num, ' PSHAP ', shape(temp_arr)
      allocate(molecules%atom_index(ubound(matrix_2d_i,2), ubound(matrix_2d_i,1 )))
      molecules%atom_index = transpose(matrix_2d_i)
   endif
   deallocate(matrix_2d_i)
   nullify(matrix_2d_i)
!write(*,*) ' INDEX   MATRIX ', shape(molecules%atom_index)
!
   call temp2%destroy
   call temp_arr%destroy
   call temp_tuple%destroy
endif
call temp%destroy
!
! Get and interpret Occupancy info
!
iarg = 21
ier_num = returned_tuple%getitem(temp, iarg)      ! Get item into temporary object
if(is_none(temp)) then
   allocate(types_occupancy(number_of_types))
   types_occupancy = 1.0_PREC_DP
else
   ier_num = cast(temp_arr, temp)             ! Now copy temp into numpy array
   ier_num = temp_arr%get_data(matrix_1d,'C')
   allocate(types_occupancy(number_of_types))
   types_occupancy = matrix_1d
!
   deallocate(matrix_1d)
   nullify(matrix_1d)
   call temp_arr%destroy
   call temp%destroy
endif
!
!
!write(*,'(a,3f12.6)') ' UNIT_CELL_LENGTHS   ', unit_cell_lengths
!write(*,'(a,3f12.6)') ' UNIT_CELL_ANGLES    ', unit_cell_angles
!write(*,'(a,3f12.6)') ' Metric Tensor (1,:) ', metric_tensor(1,:)
!write(*,'(a,3f12.6)') ' Metric Tensor (1,:) ', metric_tensor(2,:)
!write(*,'(a,3f12.6)') ' Metric Tensor (1,:) ', metric_tensor(3,:)
!write(*,'(a,a     )') ' Symmetry H_M        ', symmetry_H_M(1:len_trim(symmetry_H_M))
!write(*,'(a,i8)'    ) ' Space group origin  ', symmetry_origin
!write(*,'(a,a     )') ' Symmetry abc        ', symmetry_abc(1:len_trim(symmetry_abc))
!write(*,'(a,i8)'    ) ' Number of Symmetry  ', symmetry_n_mat 
!do i=1,symmetry_n_mat
!write(*,'(a,i8    )') ' Symmetry matrix No. ', i
!write(*,'(a,4f12.6)') ' Symmetry matrix     ', symmetry_mat(1,:,i)
!write(*,'(a,4f12.6)') ' Symmetry matrix     ', symmetry_mat(2,:,i)
!write(*,'(a,4f12.6)') ' Symmetry matrix     ', symmetry_mat(3,:,i)
!enddo
!write(*,'(a,3i8   )') ' Unit cells          ', unit_cells(1,:)
!write(*,'(a,3i8   )') ' Unit cells          ', unit_cells(2,:)
!write(*,'(a,3i8   )') ' Unit cells          ', unit_cells(3,:)
!write(*,'(a,i8)'    ) ' Number of types     ', number_of_types
!write(*,'(a,20(a4,a1))') ' Atom names        ', (types_names(i), ' ',i=1, number_of_types)
!write(*,'(a,20(i4,a1))') ' Atom ordinal      ', (types_ordinal(i), ' ',i=1, number_of_types)
!write(*,'(a,20(i4,a1))') ' Atom charge       ', (types_charge (i), ' ',i=1, number_of_types)
!write(*,'(a,20(i4,a1))') ' Atom isotope      ', (types_isotope(i), ' ',i=1, number_of_types)
!write(*,'(a,20(f8.4,a1 ))') ' Atom occupancy    ', (types_occupancy(i), ' ',i=1, number_of_types)
!write(*,'(a,i8)'    ) ' Number of atoms     ', number_of_atoms
!do i=1, number_of_atoms
!write(*,'(a,3i8,3f12.6,5i8)') 'Atom:              ', i, atom_id(i), atom_type(i), atom_pos(:,i), atom_property(i), atom_unit_cell(:,i), atom_site(i)
!enddo
!do i=1, 6
!write(*,'(2a,2l2    )') ' Crystal Flags ', c_flags(i), crystal_flags(:,i)
!enddo
!do i=1, 5
!write(*,'(a, a26, 2a )') ' Meta data     ', c_meta(i)(1:len_trim(c_meta(i))), ': ',crystal_meta(i)(1:len_trim(crystal_meta(i)))
!enddo
!write(*,'(a,i8)'    ) ' Number of ADP type  ', anisotropic_adp%anis_n_type
!do i=1, anisotropic_adp%anis_n_type
!write(*,'(a,6f10.6,2x, f9.6)') ' ADP Uij       ', anisotropic_adp%anis_adp(:,i)
!enddo
!write(*,'(a, 10i5)') ' Atoms have ADP type : ', anisotropic_adp%atom_index
!write(*,'(a,2i5)')   ' Mole: number, types   ', molecules%number_moles, molecules%number_types
!do i=1, molecules%number_moles
!write(*,'(a,4i5)'  )    ' Mole: Nr; Ty, CH, LEN ', i, molecules%mole_int(:,i)
!write(*,'(a,5x,3f9.6)') '           Ueqv, CL, CQ',    molecules%mole_real(:,molecules%mole_int(1,i))
!write(*,'(a,5x,20i5)' ) '           content     ',    molecules%atom_index(:,i)
!enddo

call p_infile%destroy             ! Output filename in python interface
call      err_print()
call p_args%destroy             ! Output filename in python interface
call      err_print()
call return_value%destroy
call      err_print()
!call returned_tuple%destroy
call interface_module%destroy
call      err_print()
!
end subroutine nx_read_structure
!
!*******************************************************************************
!
end module lib_nx_read_mod
