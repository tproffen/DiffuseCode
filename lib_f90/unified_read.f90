module unified_read_mod
!
!   Module to read the unified structure / scattering via the h5fortran library
!
private
public unified_read_structure
public unified_read_data
!
contains
!
!*******************************************************************************
!
subroutine unified_read_structure( infile, unit_cell_lengths, unit_cell_angles,           &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat, unit_cells, number_of_types, types_names,      &
                             types_ordinal, types_charge, types_isotope, number_of_atoms, &
                             atom_type, atom_pos, atom_unit_cell, atom_site,              &
                             atom_property, crystal_flags, crystal_meta,                  &
                             anisotropic_adp, molecules, average_struc, magnetic_spins,   &
                             types_occupancy,                                             &
                             l_property, l_anisotropic_adp, l_molecules, l_average_struc, &
                             l_magnetic_spins, l_types_occupancy,                         &
                             NMSG, ier_num, ier_msg)
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
character(len=*)                                         , intent(in)  :: infile
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
character(len=32)                                        , intent(out) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(out) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(out) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(out) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(out) :: symmetry_mat       ! Actual Symmetry matrices
integer                   , dimension(3    )             , intent(out) :: unit_cells         ! Number of unit cells
!
integer                                                  , intent(out) :: number_of_types    ! Crystal has this many atom types
character(len=4)          , dimension(:),     allocatable, intent(out) :: types_names        ! Name as "O", "O2-" etc
integer                   , dimension(:),     allocatable, intent(out) :: types_ordinal      ! Ordinal number of chemical element
integer                   , dimension(:),     allocatable, intent(out) :: types_charge       ! Atom types have this charge
integer                   , dimension(:),     allocatable, intent(out) :: types_isotope      ! Atom type is this isotope or zero
integer                                                  , intent(out) :: number_of_atoms    ! Crystal contains this many actual atoms
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
real(kind=PREC_DP)        , dimension(:,:  ), allocatable, intent(out) :: magnetic_spins     ! Atom has these magnetic spins
real(kind=PREC_DP)        , dimension(:),     allocatable, intent(out) :: types_occupancy    ! This atoms type has an occupancy of value
logical                                                  , intent(out) :: l_property
logical                                                  , intent(out) :: l_anisotropic_adp
logical                                                  , intent(out) :: l_molecules
logical                                                  , intent(out) :: l_average_struc
logical                                                  , intent(out) :: l_magnetic_spins
logical                                                  , intent(out) :: l_types_occupancy
integer                                                  , intent(in ) :: NMSG               ! Dimension of ier_msg
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*)          , dimension(NMSG)              , intent(out) :: ier_msg            ! Error messages
!
type(hdf5_file) :: h5f
!
integer :: i,j, k
!
character(len=256) :: string
integer, dimension(3) :: cdim
logical :: lexist
!
call h5f_reset
call h5f_open(h5f, infile, 'r', NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Unit cell parameters
!
cdim(1) = ubound(unit_cell_lengths,1)
call h5f_read(h5f, 'unit_cell_lengths', cdim(1), unit_cell_lengths, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'unit_cell_angles' , cdim(1), unit_cell_angles , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Symmetry information
call h5f_read(h5f, 'symmetry_space_group_name_H-M' , symmetry_H_M , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'space_group_origin'            , symmetry_origin, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'symmetry_space_group_abc'      , symmetry_abc   , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'space_group_symop_number'      , symmetry_n_mat , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!Symmetry matrices / transposed read
allocate(symmetry_mat(3,4,symmetry_n_mat))
cdim(1) = 3
cdim(2) = 4
cdim(3) = symmetry_n_mat
call h5f_read(h5f, 'space_group_symop_operation_mat' , cdim, symmetry_mat, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Number of unit cells
call h5f_read(h5f, 'unit_cells' , ubound(unit_cells, 1), unit_cells, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Number of atom types
call h5f_read(h5f, 'number_of_types' , number_of_types, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Number of atoms
call h5f_read(h5f, 'number_of_atoms' , number_of_atoms, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Atom names, Ordinal, Charge , Isotope
allocate(types_names(1:number_of_types))
allocate(types_ordinal(1:number_of_types))
allocate(types_charge(1:number_of_types))
allocate(types_isotope(1:number_of_types))
call h5f_read(h5f, 'types_names' , string, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
! Extract individual names
j = 1
do i=1,number_of_types-1
  k = index(string(j:len_trim(string)),';')
  types_names(i) = string(j:j+k-2)
  j = j + k 
enddo
types_names(number_of_types) = string(j:len_trim(string))
!
call h5f_read(h5f, 'types_ordinal' , number_of_types, types_ordinal, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'types_charge'  , number_of_types, types_charge , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'types_isotope' , number_of_types, types_isotope, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Atom: type, coordinate, unit cell, site_number
!
allocate(atom_type     (   number_of_atoms))
allocate(atom_pos      (3, number_of_atoms))
allocate(atom_unit_cell(3, number_of_atoms))
allocate(atom_site     (   number_of_atoms))
!
!
call h5f_read(h5f, 'atom_type' , number_of_atoms, atom_type, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
cdim(1) = 3
cdim(2) = number_of_atoms
call h5f_read(h5f, 'atom_position' , cdim(1:2), atom_pos, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'atom_unit_cell' , cdim(1:2), atom_unit_cell, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'atom_site_number' , number_of_atoms, atom_site, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!  OPTIONAL PARTS
!
!  Default to all missing
l_property = .false.
l_anisotropic_adp = .false.
l_molecules = .false.
l_average_struc = .false.
l_magnetic_spins = .false.
l_types_occupancy = .false.
!
! OPTIONAL STATUS_FLAGS
!
cdim(1) = 2
cdim(2) = ubound(crystal_flags,2)
crystal_flags  = .false.
do i=1,N_C_FLAGS
  call unified_read_c_flags(h5f, i, c_flags(i), cdim(1:2), crystal_flags)
enddo
!
crystal_meta = ' '
do i=1,N_META
  call unified_read_meta(h5f, i, c_meta(i), N_META, crystal_meta)
enddo
!
! OPTIONAL_AVERAGE_STRUCTURE
!
lexist =  h5f%exist('/entry/data/average_number')
if(lexist) then                             ! Average structure exists
   l_average_struc = .true.
   call h5f_read(h5f, 'average_number',            average_struc%aver_n_atoms, NMSG, ier_num, ier_msg)
   i = average_struc%aver_n_atoms
   if(ier_num/=0) return
   allocate(average_struc%atom_type  (   average_struc%aver_n_atoms))
   allocate(average_struc%position   (3, average_struc%aver_n_atoms))
   allocate(average_struc%occupancy  (   average_struc%aver_n_atoms))
   allocate(average_struc%anis_adp   (7, average_struc%aver_n_atoms))
   allocate(average_struc%site_number(   average_struc%aver_n_atoms))
!
   call h5f_read(h5f, 'average_type'  , i,         average_struc%atom_type   , NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
   cdim(1) = 3
   cdim(2) = i
   call h5f_read(h5f, 'average_pos'   , cdim(1:2), average_struc%position    , NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
   call h5f_read(h5f, 'average_occ'   , i,         average_struc%occupancy   , NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
   cdim(1) = 7
   cdim(2) = i
   call h5f_read(h5f, 'average_adp'   , cdim(1:2), average_struc%anis_adp    , NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
   call h5f_read(h5f, 'average_site'  , i,         average_struc%site_number , NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
endif
!
! OPTIONAL OCCUPANCY
!
lexist =  h5f%exist('/entry/data/occupancy')
if(lexist) then                               ! Occupancy exists
   l_types_occupancy = .true.
   allocate(types_occupancy(number_of_types))
   call h5f_read(h5f, 'occupancy', number_of_types, types_occupancy, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
endif
!
! OPTIONAL PROPERTY
!
lexist =  h5f%exist('/entry/data/property_flags')
if(lexist) then                               ! Occupancy exists
   l_property = .true.
   allocate(atom_property(number_of_atoms))
   call h5f_read(h5f, 'property_flags', number_of_atoms, atom_property, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
endif
!
! OPTIONAL MOLECULES
!
lexist =  h5f%exist('/entry/data/molecules_number')
if(lexist) then
   l_molecules = .true.
   call h5f_read(h5f, 'molecules_number', molecules%number_moles, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   call h5f_read(h5f, 'molecules_types', molecules%number_types, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   allocate(molecules%mole_int(3,molecules%number_moles))
   cdim(1) = 3
   cdim(2) = molecules%number_moles
   call h5f_read(h5f, 'molecules_int', cdim(1:2), molecules%mole_int, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   cdim(1) = 3
   cdim(2) = molecules%number_types
   allocate(molecules%mole_real(3,molecules%number_moles))
   call h5f_read(h5f, 'molecules_real', cdim(1:2), molecules%mole_real, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   i = maxval(molecules%mole_int(3,:))    ! Maximum molecule length
   cdim(1) = i
   cdim(2) = molecules%number_moles
   allocate(molecules%atom_index(i,molecules%number_moles))
   call h5f_read(h5f, 'molecules_index', cdim(1:2), molecules%atom_index, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
endif
!
! OPTIONAL MAGNETIC SPINS
!
lexist =  h5f%exist('/entry/data/magnetic_spins')
if(lexist) then                               ! Occupancy exists
   l_magnetic_spins = .true.
   allocate(magnetic_spins(3, number_of_atoms))
   cdim(1) = 3
   cdim(2) = number_of_atoms
   call h5f_read(h5f, 'magnetic_spins', cdim(1:2), magnetic_spins, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
endif
!
! OPTIONAL ANISOTROPIC ADP
!
lexist =  h5f%exist('/entry/data/anisotropic_number')
if(lexist) then
   l_anisotropic_adp = .true.
   call h5f_read(h5f, 'anisotropic_number', anisotropic_adp%anisotropic_n_type, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   anisotropic_adp%anisotropic_n_atom = number_of_atoms
   allocate(anisotropic_adp%anisotropic_is_iso(anisotropic_adp%anisotropic_n_type))
   allocate(anisotropic_adp%anisotropic_adp(7, number_of_atoms))
   allocate(anisotropic_adp%atom_index(number_of_atoms))
!
   call h5f_read(h5f, 'anisotropic_is_iso', anisotropic_adp%anisotropic_n_type, anisotropic_adp%anisotropic_is_iso, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   cdim(1) = 7
   cdim(2) = anisotropic_adp%anisotropic_n_type
   call h5f_read(h5f, 'anisotropic_adp', cdim(1:2), anisotropic_adp%anisotropic_adp, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
   call h5f_read(h5f, 'anisotropic_index', number_of_atoms, anisotropic_adp%atom_index, NMSG, ier_num, ier_msg)
   if(ier_num/=0) return
!
endif
!
call h5f%close()
call h5f_reset
!
end subroutine unified_read_structure
!
!*******************************************************************************
!
subroutine unified_read_data( infile, unit_cell_lengths, unit_cell_angles,           &
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
!-
!  Read the Diffuse Developers common file format for diffraction data
!+
!
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
character(len=*)                                         , intent(in)  :: infile
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
character(len=32)                                        , intent(out) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(out) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(out) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(out) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(out) :: symmetry_mat       ! Actual Symmetry matrices
!
character(len=*)                                         , intent(out) :: data_type_experiment ! Experimental = 0, Calculated = 1
character(len=*)                                         , intent(out) :: data_type_style      ! 'powder', 'powder_pdf', 'single_diffraction', 'single_pdf' ...
character(len=*)                                         , intent(out) :: data_type_content    ! 'intensity', '3d-delta-pdf', 'amplitide', 'density' ...
character(len=*)                                         , intent(out) :: data_type_reciprocal ! 'reciprocal', 'patterson', 'direct'
character(len=*)                                         , intent(out) :: data_type_with_bragg ! With Bragg   = 0, No Bragg   = 1
character(len=*)                                         , intent(out) :: data_type_symmetrized! No = 0; Point=1; Laue=2; Space=3;SymElem=4
character(len=*)                                         , intent(out) :: data_rad_radiation   ! Xray=1; Neutron=2; Electron=3
character(len=*)                                         , intent(out) :: data_rad_symbol      ! CU, CUA1, CU12
real(kind=PREC_DP)        , dimension(3)                 , intent(out) :: data_rad_length      ! Numerical value
integer                   , dimension(3)                 , intent(out) :: data_dimension   ! Data have dimensions along (HKL) / (UVW)
integer                                                  , intent(out) :: data_abs_is_hkl      ! Abscissa is 1=h 2=k 3=l
integer                                                  , intent(out) :: data_ord_is_hkl      ! Ordinate is 1=h 2=k 3=l
integer                                                  , intent(out) :: data_top_is_hkl      ! top-axis is 1=h 2=k 3=l
real(kind=PREC_DP)        , dimension(3   )              , intent(out) :: data_corner          ! Lower left bottom corner in fractional coordinates
real(kind=PREC_DP)        , dimension(3, 3)              , intent(out) :: data_vector          ! Increment vectors abs: (:,1); ord: (:,2); top: (:,3)
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(out) :: data_values  ! Actual data array

!
integer                                                  , intent(in ) :: NMSG               ! Dimension of ier_msg
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*)          , dimension(NMSG)              , intent(out) :: ier_msg            ! Error messages
!
!logical                   , dimension(8)                 , intent(in)  :: optional_intended     ! These optional flags should be present
!
character(len=PREC_STRING), dimension(  5)               , intent(in ), optional :: crystal_meta       ! Metadata, see "c_meta"
!
type(hdf5_file) :: h5f
integer, dimension(3) :: cdim
integer, dimension(11):: data_info
!
call h5f_reset
call h5f_open(h5f, infile, 'r', NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Unit cell parameters
!
cdim(1) = ubound(unit_cell_lengths,1)
call h5f_read(h5f, 'unit_cell_lengths', cdim(1), unit_cell_lengths, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'unit_cell_angles' , cdim(1), unit_cell_angles , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
! Symmetry information
call h5f_read(h5f, 'symmetry_space_group_name_H-M' , symmetry_H_M , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'space_group_origin'            , symmetry_origin, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'symmetry_space_group_abc'      , symmetry_abc   , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
call h5f_read(h5f, 'space_group_symop_number'      , symmetry_n_mat , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!Symmetry matrices / transposed read
allocate(symmetry_mat(3,4,symmetry_n_mat))
cdim(1) = 3
cdim(2) = 4
cdim(3) = symmetry_n_mat
call h5f_read(h5f, 'space_group_symop_operation_mat' , cdim, symmetry_mat, NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
!cdim(1) = 11
!call h5f_read(h5f, 'data_info'                     , cdim(1), data_info , NMSG, ier_num, ier_msg)
!if(ier_num/=0) return
!!data_dim_number       = data_info( 1) ! Data dimension 0,1,2,3
!data_type_experiment  = data_info( 2) ! Experimental = 0, Calculated = 1
!data_type_with_bragg  = data_info( 3) ! With Bragg   = 0, No Bragg   = 1
!data_type_symmetrized = data_info( 4) ! No = 0; Point=1; Laue=2; Space=3;SymElem=4
!data_rad_radiation    = data_info( 5) ! Xray=0; Neutron=1; Electron=2
!data_dimension(1)     = data_info( 6) ! Data have these dimensions along (HKL) / (UVW) abscissa
!data_dimension(2)     = data_info( 7) ! Data have these dimensions along (HKL) / (UVW) ordinate
!data_dimension(3)     = data_info( 8) ! Data have these dimensions along (HKL) / (UVW) top_axis
!
call h5f_read(h5f, 'data_type_experiment' , data_type_experiment , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_type_style' , data_type_style , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_type_content' , data_type_content , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_type_reciprocal' , data_type_reciprocal , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_type_with_bragg' , data_type_with_bragg , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_type_symmetrized' , data_type_symmetrized , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_radiation' , data_rad_radiation , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
call h5f_read(h5f, 'data_rad_symbol' , data_rad_symbol , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
cdim(1) = 3
call h5f_read(h5f, 'data_rad_length' , cdim(1), data_rad_length , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
cdim(1) = 3
call h5f_read(h5f, 'data_dimension' , cdim(1), data_dimension , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
cdim(1) =  3
call h5f_read(h5f, 'data_axes' , cdim(1), data_info , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
data_abs_is_hkl       = data_info( 1) ! Abscissa is 1=h 2=k 3=l
data_ord_is_hkl       = data_info( 2) ! Ordinate is 1=h 2=k 3=l
data_top_is_hkl       = data_info( 3) ! top-axis is 1=h 2=k 3=l
!
cdim(1) = 3
call h5f_read(h5f, 'data_corner' , cdim(1), data_corner , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
cdim(:) = 3
call h5f_read(h5f, 'data_increment_vector' , cdim   , data_vector , NMSG, ier_num, ier_msg)
if(ier_num/=0) return
!
allocate(data_values(data_dimension(1), data_dimension(2), data_dimension(3)))
cdim    = data_dimension
call h5f_read(h5f, 'data_values' , cdim   , data_values , NMSG, ier_num, ier_msg)
if(ier_num/=0) return

!
call h5f%close()
call h5f_reset
!
end subroutine unified_read_data
!
!*******************************************************************************
!
subroutine unified_read_c_flags(h5f, i, c_flag, cdim, crystal_flags)
!-
!  Check for existence and read the optional crystal_flags
!+
!
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                     , intent(in)  :: h5f
integer                             , intent(in)  :: i
character(len=*)                    , intent(in)  :: c_flag
integer, dimension(2)               , intent(in)  :: cdim
logical, dimension(cdim(1), cdim(2)), intent(out) :: crystal_flags
!
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: attribute
logical :: lexist
!
attribute = 'status_flag_' // c_flag(1:len_trim(c_flag))
lexist = h5f%exist_attr('/entry/data',attribute)
crystal_flags(:,i) = .false.
if(lexist) then
   crystal_flags(1,i) = .true.
   call h5f%readattr('/entry/data/',attribute(1:len_trim(attribute)), string)
   crystal_flags(2,i) = string == 'is_true'
endif
!
end subroutine unified_read_c_flags
!
!*******************************************************************************
!
subroutine unified_read_meta(h5f, i, c_meta, N_META, crystal_meta)
!-
!  Check for existence and read the optional crystal_meta information
!+
!
use precision_mod
!
use h5fortran, only: hdf5_file
!
implicit none
!
type(hdf5_file)                     , intent(in)  :: h5f
integer                             , intent(in)  :: i
character(len=*)                    , intent(in)  :: c_meta
integer                             , intent(in)  :: N_META
character(len=*), dimension(N_META), intent(out) :: crystal_meta
!
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: meta
logical :: lexist
!
crystal_meta(i) = ' '
meta   = '/entry/data/' // c_meta
lexist = h5f%exist(meta  )
if(lexist) then     ! Metadata exist
   call h5f%read(meta,  string)
   crystal_meta(i) = string 
endif
!
end subroutine unified_read_meta
!
!*******************************************************************************
!
end module unified_read_mod
