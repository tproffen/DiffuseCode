module nx_write_mod
!-
!  Subroutines to write into a NEXUS file
!  DiffuseDevelopers Software convention
!+
use precision_mod
!
private
!
public nx_write_scattering  ! Subroutine; write the provided scattering data into Nexus file
public nx_write_structure   ! Subroutine; write the provided crystal structure to Nexus file
!
contains
!
!*******************************************************************************
!
subroutine nx_write_scattering(outfile, program_version, author,                &
           out_inc, out_eck, out_vi,               &
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
character(len=*)                  , intent(in)  :: program_version ! ==  'DISCUS 6.18.01' or similar
character(len=*)                  , intent(in)  :: author          ! ==  'R.B. Neder    ' or similar
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
integer :: i         ! Dummy indices
integer :: qvalues_dim  ! Data dimensionality
real(kind=PREC_DP), dimension(out_inc(abs_is_hkl)) :: h_indices      ! Actual H indices along a*
real(kind=PREC_DP), dimension(out_inc(ord_is_hkl)) :: k_indices      ! Actual K indices along b*
real(kind=PREC_DP), dimension(out_inc(top_is_hkl)) :: l_indices      ! Actual L indices along c*
!
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
!write(*,*) ' Wrote original   array 321 ', qvalues(3,2,1)
!write(*,*)
!write(*,*) ' Wrote original   array llb ', qvalues(1         ,1         ,1         )
!write(*,*) ' Wrote original   array rlb ', qvalues(out_inc(1),1         ,1         )
!write(*,*) ' Wrote original   array lub ', qvalues(1         ,out_inc(2),1         )
!write(*,*) ' Wrote original   array rub ', qvalues(out_inc(1),out_inc(2),1         )
!
!write(*,*) ' Wrote original   array llt ', qvalues(1         ,1         ,out_inc(3))
!write(*,*) ' Wrote original   array rlt ', qvalues(out_inc(1),1         ,out_inc(3))
!write(*,*) ' Wrote original   array lut ', qvalues(1         ,out_inc(2),out_inc(3))
!write(*,*) ' Wrote original   array rut ', qvalues(out_inc(1),out_inc(2),out_inc(3))
!
! Write full indices H, K, L along axes
!
!write(*,*) ' abs_is_hkl ', abs_is_hkl, ord_is_hkl, top_is_hkl
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
!
ier_num = 0
call forpy_start()
!
!
ier_num = cast(          p_outfile,     outfile(1:len_trim(outfile)) )
ier_num = cast(          p_program,     program_version(1:len_trim(program_version)) )
ier_num = cast(          p_author,      author(1:len_trim(author)) )
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
!write(*,*) ' ier_num C ', ier_num
ier_num = call_py(return_value, interface_module, 'write_diffuse_scattering', p_args)
!write(*,*) ' ier_num D ', ier_num
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
!
subroutine nx_write_structure(python_script_dir, outfile, program_version, author,    &
           unit_cell_lengths, unit_cell_angles, metric_tensor,                  &
           symmetry_H_M, symmetry_origin, symmetry_abc,                         &
           symmetry_n_mat, symmetry_mat, unit_cells,                            &
           MAX_TYPES, types_names, types_ordinal, types_charge, types_isotope,  &
           MAX_ATOMS, atom_id, atom_type, atom_pos, atom_unit_cell, atom_site,  &
           status_flags, anis_n,                                                &
           ier_num,                                                             &
           property_flags,                                                      & ! Optional
           magnetic_spin,                                                       & ! Optional
           anis_adp,                                                            & ! Optional
           molecules,                                                           & ! Optional
           types_occupancy,                                                     & ! Optional
           average_struc                                                        & ! Optional
           )
!-
!  Main writing routine for Nexus Format DiffuseDevelopers Structure  format
!+
!
use precision_mod
!
use lib_forpython_mod
use lib_nx_transfer_mod
use forpy_mod
!
implicit none
!
character(len=*)                  , intent(in)  :: python_script_dir ! Full path to directory with python script
                                                                     ! '/home/reinhard/.local/share/diffuse_scattering'
character(len=*)                  , intent(in)  :: outfile           ! Output file name
character(len=*)                  , intent(in)  :: program_version   ! ==  'DISCUS 6.18.01' or similar
character(len=*)                  , intent(in)  :: author            ! ==  'R.B. Neder    ' or similar
real(kind=PREC_DP), dimension(3)  , intent(in)  :: unit_cell_lengths ! Lattice params (a, b, c)
real(kind=PREC_DP), dimension(3)  , intent(in)  :: unit_cell_angles  ! Lattice params (alpha, beta, gamma)
real(kind=PREC_DP), dimension(3,3), intent(in)  :: metric_tensor     ! Direct space metric tensor
character(len=*)                  , intent(in)  :: symmetry_H_M      ! Hermann-Mauguin Symbol 
integer                           , intent(in)  :: symmetry_origin   ! Space group origin 1 or 2
character(len=*)                  , intent(in)  :: symmetry_abc      ! Permutation 'abc', 'cba' ...
integer                           , intent(in)  :: symmetry_n_mat    ! Number of symmetry operations
real(kind=PREC_DP), dimension(3,4,symmetry_n_mat), intent(in) :: symmetry_mat ! Actual symmetry matrices (direct space)
!                                                                (x)'   (11, 12, 13, 14)   (x)
!                                                                (y)' = (21, 22, 23, 24) * (y)
!                                                                (z)'   (31, 32, 33, 34)   (z)
!                                                                (1)'   ( 0,  0,  0,  1)   (1) <== not stored!
!                                                                                         (11, 12, 13)
!                                                                (h, k, l)' = (h, k, l) * (21, 22, 23)
!                                                                                         (31, 32, 33)
integer           , dimension(3,3)        , intent(in)  :: unit_cells     ! Crystal consists of this many unit cells
integer                                   , intent(in)  :: MAX_TYPES      ! Number of atom types
character(len=*)  , dimension(MAX_TYPES)  , intent(in)  :: types_names    ! Atom types have these names
integer           , dimension(MAX_TYPES)  , intent(in)  :: types_ordinal  ! Atom types have these ordinal numberes
integer           , dimension(MAX_TYPES)  , intent(in)  :: types_charge   ! Atom types have this charge
integer           , dimension(MAX_TYPES)  , intent(in)  :: types_isotope  ! Atom types are these isotope numbers
integer                                   , intent(in)  :: MAX_ATOMS      ! Number of atoms
integer           , dimension(MAX_ATOMS)  , intent(in)  :: atom_id        ! Atoms are references by this ID
integer           , dimension(MAX_ATOMS)  , intent(in)  :: atom_type      ! Atoms have this atom_type 
real(kind=PREC_DP), dimension(3,MAX_ATOMS), intent(in)  :: atom_pos       ! Atoms are on these fractional coordinates
integer           , dimension(3,MAX_ATOMS), intent(in)  :: atom_unit_cell ! Atoms are in these unit cells
integer           , dimension(MAX_ATOMS)  , intent(in)  :: atom_site      ! Atoms are on this site in a unit cell
logical           , dimension(6)          , intent(in)  :: status_flags   ! Flags for:
!              'is_super_structure   : Crystal consists of many unit cells
!              'is_asymmetric_unit   : No symmetry operations have been applied, this is an asymmetric unit
!              'is_periodic_x        : Periodic boundary conditions can be applied along x, y, z
!              'is_periodic_y        :
!              'is_periodic_z        :
!              'is_homogeneous       : The crystal should be homogeneous, an average structure can be calculated
integer                                   , intent(in)           :: ANIS_N
!
integer           , dimension(  MAX_ATOMS), intent(in), optional :: property_flags ! Property flags for each atom
real(kind=PREC_DP), dimension(3,MAX_ATOMS), intent(in), optional :: magnetic_spin  ! Magnetic spin vector for  atom
type(anis_adp_type)                       , intent(in), optional :: anis_adp       ! Anisotropic ADPs for for type
type(molecule_data)                       , intent(in), optional :: molecules      ! Molecule info
type(average_structure)                   , intent(in), optional :: average_struc  ! Average structure 
real(kind=PREC_DP), dimension(MAX_TYPES)  , intent(in), optional :: types_occupancy! Atom types have these occupancies
             
!
integer                                   , intent(out) :: ier_num     ! Error number
!
!
type(tuple)     :: p_args              ! Tuple of arguments for write_diffuse_scattering
type(object )   :: p_outfile           ! Output filename in python interface
type(object )   :: p_program           ! Output program  in python interface
type(object )   :: p_author            ! Output author   in python interface
type(ndarray)   :: p_unit_cell_lengths ! Unit cell length in python interface
type(ndarray)   :: p_unit_cell_angles  ! Unit cell angles in python interface
type(ndarray)   :: p_metric_tensor     ! Metric tensor    in python interface
type(object )   :: p_symmetry_H_M      ! Space group name in python interface
type(object )   :: p_symmetry_origin   ! Space group symmetry has been applied to data 
type(object )   :: p_symmetry_abc      ! Space group symmetry has been applied to data 
type(object )   :: p_symmetry_n_mat    ! Number of symmetry matrices
type(ndarray)   :: p_symmetry_mat      ! Space group symmetry matrices 
type(ndarray)   :: p_unit_cells        ! Crystal has this many unit cells in python interface
type(object )   :: p_types_number      ! Crystal has these Number atom types in python interface
type(tuple  )   :: p_types_names       ! Crystal has these atom type names   in python interface
type(ndarray)   :: p_types_ordinal     ! Crystal has these atom type ordinal numbers in python interface
type(ndarray)   :: p_types_charge      ! Crystal has these atom type ID      numbers in python interface
type(ndarray)   :: p_types_isotope     ! Crystal has these atom type isotope numbers in python interface
type(object )   :: p_number_of_atoms   ! Crystal has these Number atoms       in python interface
type(ndarray)   :: p_atom_type         ! Atom is referenced by this ID number in python interface
type(ndarray)   :: p_atom_ID           ! Atom is of this type                 in python interface
type(ndarray)   :: p_atom_pos          ! Atom is at this fractional coordin.  in python interface
type(ndarray)   :: p_atom_unit_cell    ! Atom is in this unit_cell            in python interface
type(ndarray)   :: p_atom_site         ! Atom is on this site                 in python interface
type(dict   )   :: p_flags             ! Dictionary of Mandatory flags        in python interface
type(ndarray)   :: p_property_flags    ! Atom has this property_flag          in python interface
type(ndarray)   :: p_magnetic_spin     ! Atom has this magnetic spin vector   in python interface
type(ndarray)   :: p_anis_adp          ! Type has this anisotropic ADP        in python interface
type(ndarray)   :: p_anis_index        ! Type has this anisotropic ADP        in python interface
type(tuple)     :: p_anisotropic_adp   ! Anisotropic    combined
type(ndarray)   :: p_molecules_int     ! Molecular structure
type(ndarray)   :: p_molecules_real    ! Molecular structure
type(ndarray)   :: p_molecules_index   ! Molecular structure
type(tuple)     :: p_molecules         ! Molecular info combined
type(ndarray)   :: p_average_type      ! Molecular structure
type(ndarray)   :: p_average_pos       ! Molecular structure
type(ndarray)   :: p_average_occ       ! Molecular structure
type(ndarray)   :: p_average_adp       ! Molecular structure
type(ndarray)   :: p_average_site      ! Molecular structure
type(tuple)     :: p_average           ! Molecular info combined
type(ndarray)   :: p_occupancy         ! Occupancies per atom type
!
type(dict)      :: kwargs              ! Keyword specified optional arguments
!
type(module_py) :: interface_module    ! python script name
type(list)      :: paths_to_module     ! python script path
type(object)    :: return_value        ! forpy return value
!
integer :: n_arg   ! Number of arguments to python script "write_diffuse_structure"
integer :: i
!real(kind=PREC_DP), dimension(:,:), allocatable :: anis_adp       ! Anisotropic ADPs for for type
!
!write(*,*) 'WRITE_NEXUS_STRUCTURE '
!write(*,'(3a, i3)') ' OUTFILE >', outfile(1:len_trim(outfile)),'<', len_trim(outfile)
!write(*,'(3a, i3)') ' PROGRAM >', program_version(1:len_trim(program_version)),'<', len_trim(program_version)
!write(*,'(3a, i3)') ' AUTHOR  >', author(1:len_trim(author)),'<', len_trim(author)
!write(*,'(a,i3)')   ' Symm_n ', symmetry_n_mat
!write(*,'(a,3i3)')   ' SYMM MAT ', ubound(symmetry_mat)
!write(*,'(4f5.1)') symmetry_mat(1,:,1)
!write(*,'(4f5.1)') symmetry_mat(2,:,1)
!write(*,'(4f5.1)') symmetry_mat(3,:,1)
!write(*,'(a,3i4  )')   ' Unit CELL1', unit_cells(1,:)
!write(*,'(a,3i4  )')   ' Unit CELL2', unit_cells(2,:)
!write(*,'(a,3i4  )')   ' Unit CELL3', unit_cells(3,:)
!write(*,*) ' TYPES_NAMES   ', ubound(types_names), len(types_names)
!write(*,*) ' TYPES_ORDIN   ', ubound(types_ordinal)
!write(*,*) ' TYPES_TYPE    ', ubound(types_charge)
!write(*,*) ' TYPES_ISOTOPE ', ubound(types_isotope)
!
!types_names_string = 
!
ier_num = 0
call forpy_start()
!
ier_num = cast(          p_outfile,     outfile(1:len_trim(outfile)) )
ier_num = cast(          p_program,     program_version(1:len_trim(program_version)) )
ier_num = cast(          p_author,      author(1:len_trim(author)) )
ier_num = ndarray_create(p_unit_cell_lengths,  unit_cell_lengths)
ier_num = ndarray_create(p_unit_cell_angles,   unit_cell_angles )
ier_num = ndarray_create(p_metric_tensor,      metric_tensor)
ier_num = cast(          p_symmetry_H_M    , symmetry_H_M   )
ier_num = cast(          p_symmetry_origin , symmetry_origin)
ier_num = cast(          p_symmetry_abc    , symmetry_abc   )
ier_num = cast(          p_symmetry_n_mat  , symmetry_n_mat )
ier_num = ndarray_create(p_symmetry_mat    , symmetry_mat   )
ier_num = ndarray_create(p_unit_cells      , unit_cells     )
ier_num = cast(          p_types_number    , MAX_TYPES      )
!
ier_num = tuple_create(p_types_names, MAX_TYPES)
do i=1,MAX_TYPES
   ier_num = p_types_names%setitem((i-1), types_names(i))
enddo
!
ier_num = ndarray_create(p_types_ordinal   , types_ordinal  )
ier_num = ndarray_create(p_types_charge    , types_charge   )
ier_num = ndarray_create(p_types_isotope   , types_isotope  )
ier_num = cast(          p_number_of_atoms , MAX_ATOMS      )
ier_num = ndarray_create(p_atom_ID         , atom_id        )
ier_num = ndarray_create(p_atom_type       , atom_type      )
ier_num = ndarray_create(p_atom_pos        , atom_pos       )
ier_num = ndarray_create(p_atom_unit_cell  , atom_unit_cell )
ier_num = ndarray_create(p_atom_site       , atom_site      )
!
ier_num = dict_create(p_flags)
ier_num = p_flags%setitem('is_super_structure', status_flags(1))
ier_num = p_flags%setitem('is_asymmetric_unit', status_flags(2))
ier_num = p_flags%setitem('is_periodic_x'     , status_flags(3))
ier_num = p_flags%setitem('is_periodic_y'     , status_flags(4))
ier_num = p_flags%setitem('is_periodic_z'     , status_flags(5))
ier_num = p_flags%setitem('is_homogeneous'    , status_flags(6))
!
n_arg = 23   ! We have 23 mandatory arguments
!
!  Collect the arguments into a tuple: p_args
!
ier_num = tuple_create(p_args, n_arg)
ier_num = p_args%setitem( 0, p_outfile)
ier_num = p_args%setitem( 1, p_program)
ier_num = p_args%setitem( 2, p_author)
ier_num = p_args%setitem( 3, p_unit_cell_lengths)
ier_num = p_args%setitem( 4, p_unit_cell_angles)
ier_num = p_args%setitem( 5, p_metric_tensor)
ier_num = p_args%setitem( 6, p_symmetry_H_M)
ier_num = p_args%setitem( 7, p_symmetry_origin)
ier_num = p_args%setitem( 8, p_symmetry_abc)
ier_num = p_args%setitem( 9, p_symmetry_n_mat)
ier_num = p_args%setitem(10, p_symmetry_mat)
ier_num = p_args%setitem(11, p_unit_cells)
ier_num = p_args%setitem(12, p_types_number    )
ier_num = p_args%setitem(13, p_types_names     )
ier_num = p_args%setitem(14, p_types_ordinal   )
ier_num = p_args%setitem(15, p_types_charge    )
ier_num = p_args%setitem(16, p_types_isotope   )
ier_num = p_args%setitem(17, p_number_of_atoms )
ier_num = p_args%setitem(18, p_atom_id         )
ier_num = p_args%setitem(19, p_atom_type       )
ier_num = p_args%setitem(20, p_atom_pos        )
ier_num = p_args%setitem(21, p_atom_unit_cell  )
ier_num = p_args%setitem(22, p_atom_site       )
!
! Optional keyword specified arguments
!
ier_num = dict_create(kwargs)                 ! Create dictionary for keyword arguments

ier_num = kwargs%setitem("status_flags", p_flags           )
!
! Property Flags
!
if(present(property_flags)) then              ! property_flags are present
   ier_num = ndarray_create(p_property_flags   , property_flags  )
   ier_num = kwargs%setitem("property_flags", p_property_flags  )
endif
!
! Magnetic spin vectors
!
if(present(magnetic_spin)) then               ! Magnetic spins are present
   ier_num = ndarray_create(p_magnetic_spin   , magnetic_spin  )
   ier_num = kwargs%setitem("magnetic_spin", p_magnetic_spin   )
endif
!
! Anisotropic ADP
!
if(present(anis_adp     )) then               ! anis_adp     s are present
!write(*,*) ' ANIS_ADP       present', ubound(anis_adp%anis_adp), anis_adp%anis_n_type
!write(*,*) ' ANIS_ADP ', anis_adp%anis_adp(:,1)
!write(*,*) ' ANIS_INDX', anis_adp%atom_index(:)
   ier_num = ndarray_create(p_anis_adp  , (anis_adp%anis_adp))
!write(*,*) ' NDARRAY anis_adp ', ier_num
   ier_num = ndarray_create(p_anis_index, anis_adp%atom_index      )
!write(*,*) ' NDARRAY anis_ind ', ier_num
   !ier_num = print_py(p_anis_adp)
   n_arg = 2
   ier_num = tuple_create(p_anisotropic_adp, n_arg)
   ier_num = p_anisotropic_adp%setitem( 0, p_anis_adp  )
   ier_num = p_anisotropic_adp%setitem( 1, p_anis_index)
   ier_num = kwargs%setitem("anisotropic_adp"     , p_anisotropic_adp)
   !ier_num = kwargs%setitem("anisotropic_adp"     , p_anis_adp        )
   !ier_num = kwargs%setitem("anisotropic_index"   , p_anis_index      )
endif
!
! Molecule Information
!
if(present(molecules)) then     ! Crystal contains molecules
   if(molecules%number_moles>0) then         ! Crystal does actually contain molecules
   ier_num = ndarray_create(p_molecules_int   , (molecules%mole_int ))
   ier_num = ndarray_create(p_molecules_real  , (molecules%mole_real))
   ier_num = ndarray_create(p_molecules_index , (molecules%atom_index))
   n_arg = 5
   ier_num = tuple_create(p_molecules, n_arg)
   ier_num = p_molecules%setitem( 0, molecules%number_moles)
   ier_num = p_molecules%setitem( 1, molecules%number_types)
   ier_num = p_molecules%setitem( 2, p_molecules_int)
   ier_num = p_molecules%setitem( 3, p_molecules_real)
   ier_num = p_molecules%setitem( 4, p_molecules_index)
   ier_num = kwargs%setitem("molecules"   , p_molecules      )
   endif
endif
!
! Occupancies
!
if(present(types_occupancy)) then  ! Crystal contains occupancy information
   ier_num = ndarray_create(p_occupancy, types_occupancy)
   ier_num = kwargs%setitem("occupancy"   , p_occupancy      )
endif
!
! Average structure
!
if(present(average_struc)) then   ! Average structure info is known
   if(average_struc%aver_n_atoms > 0) then  ! Average structure contains atoms
      ier_num = ndarray_create(p_average_type, average_struc%atom_type)
      ier_num = ndarray_create(p_average_pos,  average_struc%position)
      ier_num = ndarray_create(p_average_occ , average_struc%occupancy)
      ier_num = ndarray_create(p_average_adp,  average_struc%anis_adp)
      ier_num = ndarray_create(p_average_site, average_struc%site_number)
      n_arg = 5
      ier_num = tuple_create(p_average, n_arg)
      ier_num = p_average%setitem( 0, p_average_type)
      ier_num = p_average%setitem( 1, p_average_pos )
      ier_num = p_average%setitem( 2, p_average_occ )
      ier_num = p_average%setitem( 3, p_average_adp )
      ier_num = p_average%setitem( 4, p_average_site)
      ier_num = kwargs%setitem("average_structure", p_average      )
   endif 
endif
!
! Append current directory to paths
!
ier_num = get_sys_path(paths_to_module)
!write(*,*) ' ier_num A ', ier_num
call err_print()
!write(*,*) ' PATH ', python_script_dir(1:len_trim(python_script_dir))
ier_num = paths_to_module%append(python_script_dir(1:len_trim(python_script_dir))) !'.')
!write(*,*) ' ier_num B ', ier_num
call err_print()
!ier_num = print_py(paths_to_module)
!call err_print()
ier_num = import_py(interface_module, 'nexus_structure')
!write(*,*) ' ier_num C ', ier_num
call err_print()
ier_num = call_py(return_value, interface_module, 'write_diffuse_structure', p_args, kwargs)
!write(*,*) ' ier_num D ', ier_num
call err_print()
!
call p_outfile%destroy
call p_program%destroy
call p_author%destroy
call p_unit_cell_lengths%destroy
call p_unit_cell_angles%destroy
call p_metric_tensor%destroy
call p_symmetry_H_M%destroy
call p_symmetry_origin%destroy
call p_symmetry_abc%destroy
call p_symmetry_n_mat%destroy
call p_symmetry_mat%destroy
call p_unit_cells%destroy
call p_types_number%destroy
call p_types_names%destroy
call p_types_ordinal%destroy
call p_types_charge%destroy
call p_types_isotope%destroy
call p_number_of_atoms%destroy
call p_atom_id%destroy
call p_atom_type%destroy
call p_atom_pos%destroy
call p_atom_unit_cell%destroy
call p_atom_site%destroy
!
end subroutine nx_write_structure
!
!*******************************************************************************
!
end module nx_write_mod
