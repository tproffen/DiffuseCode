module unified_write_mod
!
!   Module to read the unified structure / scattering via the h5fortran library
!
private
public unified_write_structure
public unified_write_data
!
contains
!
!*******************************************************************************
!
subroutine unified_write_structure( infile, unit_cell_lengths, unit_cell_angles,          &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat, unit_cells, number_of_types, types_names,      &
                             types_ordinal, types_charge, types_isotope, number_of_atoms, &
                             atom_type, atom_pos, atom_unit_cell, atom_site,              &
                             NMSG, ier_num, ier_msg,                                      &
                             optional_intended,                                              &
                             atom_property, crystal_flags, crystal_meta,                  &
                             anisotropic_adp, molecules, average_struc, magnetic_spins,   &
                             types_occupancy)
!
! low level read
!
use lib_h5fortran_mod
use lib_nx_transfer_mod
use lib_unified_chars_mod
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
character(len=*)                                         , intent(in) :: infile
real(kind=PREC_DP)        , dimension(3)                 , intent(in) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(in) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
character(len=32)                                        , intent(in) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(in) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(in) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(in) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(3,4,symmetry_n_mat), intent(in) :: symmetry_mat       ! Actual Symmetry matrices
integer                   , dimension(3    )             , intent(in) :: unit_cells         ! Number of unit cells
!
integer                                                  , intent(in) :: number_of_types    ! Crystal has this many atom types
character(len=4)          , dimension(number_of_types)   , intent(in) :: types_names        ! Name as "O", "O2-" etc
integer                   , dimension(number_of_types)   , intent(in) :: types_ordinal      ! Ordinal number of chemical element
integer                   , dimension(number_of_types)   , intent(in) :: types_charge       ! Atom types have this charge
integer                   , dimension(number_of_types)   , intent(in) :: types_isotope      ! Atom type is this isotope or zero
integer                                                  , intent(in) :: number_of_atoms    ! Crystal contains this many actual atoms
integer                   , dimension(  number_of_atoms) , intent(in) :: atom_type          ! Atom is of this type
real(kind=PREC_DP)        , dimension(3,number_of_atoms) , intent(in) :: atom_pos           ! Atom is at these fractional coordinates
integer                   , dimension(3,number_of_atoms) , intent(in) :: atom_unit_cell     ! Atom is in this unit cell
integer                   , dimension(  number_of_atoms) , intent(in) :: atom_site          ! Atom is on this site in its unit cell
!
integer                                                  , intent(in ) :: NMSG               ! Dimension of ier_msg
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*)          , dimension(NMSG)              , intent(out) :: ier_msg            ! Error messages
!
logical                   , dimension(8)                 , intent(in)  :: optional_intended     ! This optional flags should be present
!
integer                   , dimension(  number_of_atoms) , intent(in) , optional :: atom_property      ! Atom has this property flag
logical                   , dimension(2,6)               , intent(in ), optional :: crystal_flags      ! Flags , see "c_flags"
character(len=PREC_STRING), dimension(  5)               , intent(in ), optional :: crystal_meta       ! Metadata, see "c_meta"
type(anis_adp_type)                                      , intent(in ), optional :: anisotropic_adp    ! Info on anisotropic ADP ==> lib_nx_transfer.f90
type(molecule_data)                                      , intent(in ), optional :: molecules          ! Info on molecules       ==> lib_nx_transfer.f90
type(average_structure)                                  , intent(in ), optional :: average_struc      ! Info on average struct  ==> lib_nx_transfer.f90
real(kind=PREC_DP)        , dimension(3, number_of_atoms), intent(in ), optional :: magnetic_spins     ! Atom has these magnetic spins
real(kind=PREC_DP)        , dimension(   number_of_types), intent(in ), optional :: types_occupancy    ! This atoms type has an occupancy of value
!
!
type(hdf5_file) :: h5f           ! The h5fortran library object
!
integer :: i
!
character(len=256) :: string
integer, dimension(3) :: cdim
!
!
call h5f_reset
call h5f_open(h5f, infile, 'w', NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f%create_group('/entry')
call h5f%create_group('/entry/data')
!
! Unit cell parameters
!
cdim(1) = ubound(unit_cell_lengths,1)
call h5f_write(h5f, 'unit_cell_lengths', cdim(1), unit_cell_lengths, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'unit_cell_angles' , cdim(1), unit_cell_angles , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Symmetry information
call h5f_write(h5f, 'symmetry_space_group_name_H-M' , symmetry_H_M , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'space_group_origin'            , symmetry_origin, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'symmetry_space_group_abc'      , symmetry_abc   , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'space_group_symop_number'      , symmetry_n_mat , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!Symmetry matrices / transposed write
cdim(1) = 3
cdim(2) = 4
cdim(3) = symmetry_n_mat
call h5f_write(h5f, 'space_group_symop_operation_mat' , cdim, symmetry_mat, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Number of unit cells along x, y, c
call h5f_write(h5f, 'unit_cells' , ubound(unit_cells, 1), unit_cells, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Number of atom types
call h5f_write(h5f, 'number_of_types' , number_of_types, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Number of atoms
call h5f_write(h5f, 'number_of_atoms' , number_of_atoms, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Atom names, Ordinal, Charge , Isotope
string = types_names(1)(1:len_trim(types_names(1)))     ! Concatenate all atom names into a string
do i=2, number_of_types
   string = string(1:len_trim(string)) // ';' // types_names(i)(1:len_trim(types_names(i)))
enddo
call h5f_write(h5f, 'types_names' , string, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_write(h5f, 'types_ordinal' , number_of_types, types_ordinal, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'types_charge'  , number_of_types, types_charge , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'types_isotope' , number_of_types, types_isotope, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Atom: type, coordinate, unit cell, site_number
!
call h5f_write(h5f, 'atom_type' , number_of_atoms, atom_type, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
cdim(1) = 3
cdim(2) = number_of_atoms
call h5f_write(h5f, 'atom_position' , cdim(1:2), atom_pos, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_write(h5f, 'atom_unit_cell' , cdim(1:2), atom_unit_cell, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_write(h5f, 'atom_site_number' , number_of_atoms, atom_site, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!  OPTIONAL PARTS
!
! OPTIONAL STATUS_FLAGS
!
if(present(crystal_flags)) then     ! Crystal flags are present
   if(optional_intended(o_crystal_flags)) then
      cdim(1) = 2
      cdim(2) = ubound(crystal_flags,2)
      call unified_write_c_flags(h5f, cdim(1:2), crystal_flags)
   endif
endif
!
if(present(crystal_meta)) then
   if(optional_intended(o_crystal_meta)) then
      call unified_write_meta(h5f, N_META, crystal_meta)
   endif
endif
!
! OPTIONAL_AVERAGE_STRUCTURE
!
if(present(average_struc)) then
   if(optional_intended(o_average_struc)) then
   i = average_struc%aver_n_atoms
   call h5f_write(h5f, 'average_number',            average_struc%aver_n_atoms, NMSG, ier_num, ier_msg)
   call h5f_write(h5f, 'average_type'  , i,         average_struc%atom_type   , NMSG, ier_num, ier_msg)
   cdim(1) = 3
   cdim(2) = i
   call h5f_write(h5f, 'average_pos'   , cdim(1:2), average_struc%position    , NMSG, ier_num, ier_msg)
   call h5f_write(h5f, 'average_occ'   , i,         average_struc%occupancy   , NMSG, ier_num, ier_msg)
   cdim(1) = 7
   cdim(2) = i
   call h5f_write(h5f, 'average_adp'   , cdim(1:2), average_struc%anis_adp    , NMSG, ier_num, ier_msg)
   call h5f_write(h5f, 'average_site'  , i,         average_struc%site_number , NMSG, ier_num, ier_msg)
   endif
endif
!
! OPTIONAL OCCUPANCY
!
if(present(types_occupancy)) then
   if(optional_intended(o_types_occupancy)) then
      call h5f_write(h5f, 'occupancy', number_of_types, types_occupancy, NMSG, ier_num, ier_msg)
   endif
endif
!
! OPTIONAL PROPERTY
!
if(present(atom_property)) then
   if(optional_intended(o_atom_property)) then
      call h5f_write(h5f, 'property_flags', number_of_atoms, atom_property, NMSG, ier_num, ier_msg)
      if(ier_num/=0) return
   endif
endif
!
! OPTIONAL MOLECULES
!
if(present(molecules)) then
   if(optional_intended(o_molecules)) then
   call h5f_write(h5f, 'molecules_number', molecules%number_moles, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   call h5f_write(h5f, 'molecules_types', molecules%number_types, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   cdim(1) = 3
   cdim(2) = molecules%number_moles
   call h5f_write(h5f, 'molecules_int', cdim(1:2), molecules%mole_int, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   cdim(1) = 3
   cdim(2) = molecules%number_types
   call h5f_write(h5f, 'molecules_real', cdim(1:2), molecules%mole_real, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   i = maxval(molecules%mole_int(3,:))    ! Maximum molecule length
   cdim(1) = i
   cdim(2) = molecules%number_moles
   call h5f_write(h5f, 'molecules_index', cdim(1:2), molecules%atom_index, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   endif
endif
!
! OPTIONAL MAGNETIC SPINS
!
if(present(magnetic_spins)) then
   if(optional_intended(o_magnetic_spins)) then
   cdim(1) = 3
   cdim(2) = number_of_atoms
   call h5f_write(h5f, 'magnetic_spins', cdim(1:2), magnetic_spins, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   endif
endif
!
! OPTIONAL ANISOTROPIC ADP
!
if(present(anisotropic_adp)) then
   if(optional_intended(o_anisotropic_adp)) then
!write(*,*) ' LOW LEVEL '
!write(*,*) ' ADP BOUND ', ubound(anisotropic_adp%anisotropic_adp)
!write(*,'(a,12i4)')           ' ADP index ', anisotropic_adp%atom_index(1:12)
!write(*,'(a,12i4)')           ' ADP number', anisotropic_adp%anisotropic_n_type
!do i=1, anisotropic_adp%anisotropic_n_type
!write(*,'(a,6f10.6,a,f10.6)') ' ADP uij   ', anisotropic_adp%anisotropic_adp(1:6,i), ' || ',&
!anisotropic_adp%anisotropic_adp(7,i)
!enddo
   call h5f_write(h5f, 'anisotropic_number', anisotropic_adp%anisotropic_n_type, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
!  anisotropic_adp%anisotropic_n_atom = number_of_atoms
!
   call h5f_write(h5f, 'anisotropic_is_iso', anisotropic_adp%anisotropic_n_type, anisotropic_adp%anisotropic_is_iso, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   cdim(1) = 7
   cdim(2) = anisotropic_adp%anisotropic_n_type
   call h5f_write(h5f, 'anisotropic_adp', cdim(1:2), anisotropic_adp%anisotropic_adp, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   call h5f_write(h5f, 'anisotropic_index', number_of_atoms, anisotropic_adp%atom_index, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   endif
endif
!
call h5f%close()
call h5f_reset
!
end subroutine unified_write_structure
!
!*******************************************************************************
!
subroutine unified_write_data(datafile, unit_cell_lengths, unit_cell_angles,               &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat,                                                &
                             data_type_experiment, &
                             data_type_style     , &
                             data_type_content   , &
                             data_type_reciprocal, &
                             data_type_with_bragg, &
                             data_type_symmetrized, &
                             data_rad_radiation , data_rad_symbol, data_rad_length, &
                             data_dimension  , &
                             data_abs_is_hkl     , &
                             data_ord_is_hkl     , &
                             data_top_is_hkl     , &
                             data_corner     , &
                             data_vector     , &
                             data_values     , &
                             NMSG, ier_num, ier_msg,                                      &
                                            crystal_meta                   &
                             )
!
! low level read
!
use lib_h5fortran_mod
use lib_nx_transfer_mod
use lib_unified_chars_mod
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
character(len=*)                                         , intent(in) :: datafile
real(kind=PREC_DP)        , dimension(3)                 , intent(in) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(in) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
character(len=32)                                        , intent(in) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(in) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(in) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(in) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(3,4,symmetry_n_mat), intent(in) :: symmetry_mat       ! Actual Symmetry matrices
!
character(len=*)                                         , intent(in) :: data_type_experiment ! 'experimental';  'calculated'
character(len=*)                                         , intent(in) :: data_type_style      ! 'powder_diffraction', 'powder_pdf', 'single_diffraction', 'single_pdf' ...
character(len=*)                                         , intent(in) :: data_type_content    ! 'intensity', '3d-delta-pdf', 'amplitide', 'density' ...
character(len=*)                                         , intent(in) :: data_type_reciprocal ! 'reciprocal', 'patterson', 'direct'
character(len=*)                                         , intent(in) :: data_type_with_bragg ! 'bragg_present'; 'bragg_subtracted'
character(len=*)                                         , intent(in) :: data_type_symmetrized! 'none'; 'point'; 'laue'; 'space', 'sym_elem'
character(len=*)                                         , intent(in) :: data_rad_radiation   ! 'xray'; 'neutron'; 'electron'
character(len=*)                                         , intent(in) :: data_rad_symbol      ! CU, CUA1, CU12
real(kind=PREC_DP)        , dimension(3)                 , intent(in) :: data_rad_length      ! Numerical value
integer                   , dimension(3)                 , intent(in) :: data_dimension   ! Data have dimensions along (HKL) / (UVW)
integer                                                  , intent(in) :: data_abs_is_hkl      ! Abscissa is 1=h 2=k 3=l
integer                                                  , intent(in) :: data_ord_is_hkl      ! Ordinate is 1=h 2=k 3=l
integer                                                  , intent(in) :: data_top_is_hkl      ! top-axis is 1=h 2=k 3=l
real(kind=PREC_DP)        , dimension(3   )              , intent(in) :: data_corner          ! Lower left bottom corner in fractional coordinates
real(kind=PREC_DP)        , dimension(3, 3)              , intent(in) :: data_vector          ! Increment vectors abs: (:,1); ord: (:,2); top: (:,3)
real(kind=PREC_DP)        , dimension(data_dimension(1), &
                                      data_dimension(2), &
                                      data_dimension(3)), intent(in) :: data_values  ! Actual data array

!
integer                                                  , intent(in ) :: NMSG               ! Dimension of ier_msg
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*)          , dimension(NMSG)              , intent(out) :: ier_msg            ! Error messages
!
!logical                   , dimension(8)                 , intent(in)  :: optional_intended     ! These optional flags should be present
!
character(len=PREC_STRING), dimension(  5)               , intent(in ), optional :: crystal_meta       ! Metadata, see "c_meta"
!
!
type(hdf5_file) :: h5f           ! The h5fortran library object
!
integer :: i
!
character(len=1), dimension(3), parameter :: chkl = (/'h', 'k', 'l'/)
character(len=1), dimension(3), parameter :: cxyz = (/'x', 'y', 'z'/)
character(len=1), dimension(3), parameter :: cuvw = (/'u', 'v', 'w'/)
character(len=PREC_STRING)    :: axes ! ==  '["h", "k", "l"]' or similar permutation
character(len=256) :: string
integer, dimension(3) :: cdim
!integer                :: data_dim_number      ! Data have 0,1,2,3 dimensions
integer, dimension(11) :: data_info
!
!data_dim_number = 0
!do i=1, 3
!   if(data_dimension(i)>1) data_dim_number = data_dim_number +1  ! Increment for each dimension
!enddo
!
!  Set proper axes names at abscissa, ordinate, top_axis
!
if(data_type_reciprocal=='reciprocal') then
   axes      = '["h", "k", "l"]'
   axes( 3: 3) = chkl(data_abs_is_hkl)    ! Place the correct 'h', 'k', or 'l'
   axes( 8: 8) = chkl(data_ord_is_hkl)
   axes(13:13) = chkl(data_top_is_hkl)
elseif(data_type_reciprocal=='direct') then
   axes      = '["x", "y", "z"]'
   axes( 3: 3) = cxyz(data_abs_is_hkl)    ! Place the correct 'x', 'y', or 'z'
   axes( 8: 8) = cxyz(data_ord_is_hkl)
   axes(13:13) = cxyz(data_top_is_hkl)
elseif(data_type_reciprocal=='patterson') then
   axes      = '["u", "v", "w"]'
   axes( 3: 3) = cxyz(data_abs_is_hkl)    ! Place the correct 'u', 'v', or 'w'
   axes( 8: 8) = cxyz(data_ord_is_hkl)
   axes(13:13) = cxyz(data_top_is_hkl)
endif
!
data_info = 0           ! Generic integer array
!
call h5f_reset
call h5f_open(h5f, datafile, 'w', NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f%create_group('/entry')
call h5f%create_group('/entry/data')
!
! Unit cell parameters
!
cdim(1) = ubound(unit_cell_lengths,1)
call h5f_write(h5f, 'unit_cell_lengths', cdim(1), unit_cell_lengths, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'unit_cell_angles' , cdim(1), unit_cell_angles , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Symmetry information
call h5f_write(h5f, 'symmetry_space_group_name_H-M' , symmetry_H_M , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'space_group_origin'            , symmetry_origin, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'symmetry_space_group_abc'      , symmetry_abc   , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_write(h5f, 'space_group_symop_number'      , symmetry_n_mat , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!Symmetry matrices / transposed write
cdim(1) = 3
cdim(2) = 4
cdim(3) = symmetry_n_mat
call h5f_write(h5f, 'space_group_symop_operation_mat' , cdim, symmetry_mat, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!cdim(1) = ubound(data_info,1)
!call h5f_write(h5f, 'data_info', cdim(1), data_info, NMSG, ier_num, ier_msg)
!if(ier_num/=0) return
!
! Data type experiment/calculated  information
call h5f_write(h5f, 'data_type_experiment' , data_type_experiment , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data type style  information
call h5f_write(h5f, 'data_type_style' , data_type_style , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data type content information
call h5f_write(h5f, 'data_type_content' , data_type_content , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data type space information
call h5f_write(h5f, 'data_type_reciprocal' , data_type_reciprocal , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data type with_bragg information
call h5f_write(h5f, 'data_type_with_bragg' , data_type_with_bragg , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data type symmetrized
call h5f_write(h5f, 'data_type_symmetrized' , data_type_symmetrized , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data radiation
call h5f_write(h5f, 'data_radiation', data_rad_radiation, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data radiation symbol
call h5f_write(h5f, 'data_rad_symbol', data_rad_symbol, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data radiation wave length
cdim(1) = 3
call h5f_write(h5f, 'data_rad_length', cdim(1), data_rad_length, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data dimension
cdim(1) = 3
call h5f_write(h5f, 'data_dimension' , cdim(1), data_dimension, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Data axes types
data_info( 1) = data_abs_is_hkl       ! Abscissa is 0=unknown 1=h 2=k 3=l
data_info( 2) = data_ord_is_hkl       ! Ordinate is 0=unknown 1=h 2=k 3=l
data_info( 3) = data_top_is_hkl       ! top-axis is 0=unknown 1=h 2=k 3=l
cdim(1) = 3
call h5f_write(h5f, 'data_axes' , cdim(1), data_info(1:3), NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Lower_left_bottom corner / transposed write
cdim(1) = 3
call h5f_write(h5f, 'data_corner' , cdim(1), data_corner, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Increment vectors / transposed write
cdim(1) = 3
cdim(2) = 3
cdim(3) = 3
call h5f_write(h5f, 'data_increment_vector' , cdim, data_vector, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Actual data       / transposed write
cdim(1) = data_dimension(1)
cdim(2) = data_dimension(2)
cdim(3) = data_dimension(3)
call h5f_write(h5f, 'data_values' , cdim, data_values, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!
if(present(crystal_meta)) then
!  if(optional_intended(o_crystal_meta)) then
      call unified_write_meta(h5f, N_META, crystal_meta)
!  endif
endif
!
call h5f%close()
call h5f_reset
!
!
end subroutine unified_write_data
!
!*******************************************************************************
!
subroutine unified_write_c_flags(h5f, cdim, crystal_flags)
!-
!  Write the optional crystal_flags
!+
!
use lib_unified_chars_mod
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                     , intent(in)  :: h5f
integer, dimension(2)               , intent(in)  :: cdim
logical, dimension(cdim(1), cdim(2)), intent(in)  :: crystal_flags
!
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: attribute
integer :: i
!
do i=1, cdim(2)
   if(crystal_flags(1,i)) then
      if(crystal_flags(2,i)) then
         string = 'is_true'
      else
         string = 'is_false'
   attribute = 'status_flag_' // c_flags(i)(1:len_trim(c_flags(i)))
      endif
      attribute = 'status_flag_' // c_flags(i)(1:len_trim(c_flags(i)))
      call h5f%writeattr('/entry/data/',attribute(1:len_trim(attribute)), string)
   endif
enddo
!
end subroutine unified_write_c_flags
!
!*******************************************************************************
!
subroutine unified_write_meta(h5f, NN_META, crystal_meta)
!-
!  Write the optional crystal_meta information
!+
!
use lib_unified_chars_mod
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                     , intent(in)  :: h5f
integer                             , intent(in)  :: NN_META
character(len=*), dimension(NN_META), intent(in) :: crystal_meta
!
character(len=PREC_STRING) :: meta
integer :: i
!logical :: lexist
!
do i=1, N_META
   meta   = '/entry/data/' // c_meta(i)
   call h5f%write(meta,  crystal_meta(i))
enddo
!
end subroutine unified_write_meta
!
!*******************************************************************************
!
end module unified_write_mod
