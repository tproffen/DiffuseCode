module unified_dump_mod
!-
!  For debug purposes write a summary of the unified structure / data
!+
!
private
public unified_dump_structure
public unified_dump_data
!
!*******************************************************************************
!
contains
!
!*******************************************************************************
!
subroutine unified_dump_structure( unit_cell_lengths, unit_cell_angles,           &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat, unit_cells, number_of_types, types_names,      &
                             types_ordinal, types_charge, types_isotope, coordinate_unit, &
                             number_of_atoms, &
                             atom_type, atom_pos, atom_unit_cell, atom_site,              &
                             atom_property, crystal_flags, crystal_meta,                  &
                             anisotropic_adp, molecules, average_struc, magnetic_spins,   &
                             types_occupancy,                                             &
                             l_property, l_anisotropic_adp, l_molecules, l_average_struc, &
                             l_magnetic_spins, l_types_occupancy,                         &
                             NMSG, ier_num, ier_msg)
!
! Mostly for debug purposes dump the structure to screen; might be  long
!
!use lib_h5fortran_mod
use lib_nx_transfer_mod
use lib_unified_chars_mod
use precision_mod
!
!
implicit none
!
real(kind=PREC_DP)        , dimension(3)                 , intent(in ) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(in ) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
character(len=32)                                        , intent(in ) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(in ) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(in ) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(in ) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(3,4,symmetry_n_mat), intent(in ) :: symmetry_mat       ! Actual Symmetry matrices
integer                   , dimension(3    )             , intent(in ) :: unit_cells         ! Number of unit cells
!
integer                                                  , intent(in ) :: number_of_types    ! Crystal has this many atom types
character(len=4)          , dimension(number_of_types)   , intent(in ) :: types_names        ! Name as "O", "O2-" etc
integer                   , dimension(number_of_types)   , intent(in ) :: types_ordinal      ! Ordinal number of chemical element
integer                   , dimension(number_of_types)   , intent(in ) :: types_charge       ! Atom types have this charge
integer                   , dimension(number_of_types)   , intent(in ) :: types_isotope      ! Atom type is this isotope or zero
integer                                                  , intent(in ) :: number_of_atoms    ! Crystal contains this many actual atoms
character(len=32)                                        , intent(in ) :: coordinate_unit    ! Units for coordinates 'basecell_fractional', isupercell_fractional'
integer                   , dimension(  number_of_atoms) , intent(in ) :: atom_type          ! Atom is of this type
real(kind=PREC_DP)        , dimension(3,number_of_atoms) , intent(in ) :: atom_pos           ! Atom is at these fractional coordinates
integer                   , dimension(3,number_of_atoms) , intent(in ) :: atom_unit_cell     ! Atom is in this unit cell
integer                   , dimension(  number_of_atoms) , intent(in ) :: atom_site          ! Atom is on this site in its unit cell
integer                   , dimension(  number_of_atoms) , intent(in ) :: atom_property      ! Atom has this property flag
logical                   , dimension(2,6)               , intent(in ) :: crystal_flags      ! Flags , see "c_flags"
character(len=PREC_STRING), dimension(  5)               , intent(in ) :: crystal_meta       ! Metadata, see "c_meta"
type(anis_adp_type)                         , optional   , intent(in ) :: anisotropic_adp    ! Info on anisotropic ADP ==> lib_nx_transfer.f90
type(molecule_data)                         , optional   , intent(in ) :: molecules          ! Info on molecules       ==> lib_nx_transfer.f90
type(average_structure)                     , optional   , intent(in ) :: average_struc      ! Info on average struct  ==> lib_nx_transfer.f90
real(kind=PREC_DP)        , dimension(:,:  ),              optional, intent(in ) :: magnetic_spins     ! Atom has these magnetic spins
real(kind=PREC_DP)        , dimension(number_of_types)   , optional, intent(in ) :: types_occupancy    ! This atoms type has an occupancy of value
logical                                                  , intent(in ) :: l_property
logical                                                  , intent(in ) :: l_anisotropic_adp
logical                                                  , intent(in ) :: l_molecules
logical                                                  , intent(in ) :: l_average_struc
logical                                                  , intent(in ) :: l_magnetic_spins
logical                                                  , intent(in ) :: l_types_occupancy
integer                                                  , intent(in ) :: NMSG               ! Dimension of ier_msg
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*)          , dimension(NMSG)              , intent(out) :: ier_msg            ! Error messages
!
integer :: i
!
write(*,*)
write(*,'(a)') '      Partial dump of Crystal structure '
write(*,*)
write(*,'(a)') '   Crystal Meta information'
do i= 1, N_META
   write(*,'(a26, 1x, a)') c_meta(i), crystal_meta(i)(1:len_trim(crystal_meta(i)))
enddo
write(*,*)
write(*,'(a)') 'Crystal Flag              Status; is_intended'
do i= 1, N_C_FLAGS
   write(*,'(a26, 1x, l3, 10x,l3)') c_flags(i), crystal_flags(:,i)
enddo
write(*,*)
!
write(*,'(a,3f12.6)')  ' Unit cell params ', unit_cell_lengths
write(*,'(a,3f12.6)')  ' Unit cell angles ', unit_cell_angles
!
write(*,'(a, a     )') ' Symmetry H_M     ', symmetry_H_M
write(*,'(a, i5    )') ' Symmetry Origin  ', symmetry_origin
write(*,'(a, a     )') ' Symmetry abc     ', symmetry_abc
write(*,'(a, i5    )') ' Symmetry n MAT   ', symmetry_n_mat
write(*,'(a, 4f11.6)') ' Symmetry 1 fst   ', symmetry_mat(1,:,1)
write(*,'(a, 4f11.6)') ' Symmetry 1 scnd  ', symmetry_mat(2,:,1)
write(*,'(a, 4f11.6)') ' Symmetry 1 lst   ', symmetry_mat(3,:,1)
if(symmetry_n_mat>10) then
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'fst  ',symmetry_mat(1,:,symmetry_n_mat/2)
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'scnd ',symmetry_mat(2,:,symmetry_n_mat/2)
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'lst  ',symmetry_mat(3,:,symmetry_n_mat/2)
endif
write(*,'(a, 4f11.6)') ' Symmetry N fst   ', symmetry_mat(1,:,symmetry_n_mat)
write(*,'(a, 4f11.6)') ' Symmetry N scnd  ', symmetry_mat(2,:,symmetry_n_mat)
write(*,'(a, 4f11.6)') ' Symmetry N lst   ', symmetry_mat(3,:,symmetry_n_mat)
write(*,*)
!
write(*,'(a, i8)') ' Number of atom types ', number_of_types
write(*,'(a)') '  Number  Name  Ordinal  Charge  Isotope'
do i=1, number_of_types
write(*,'(i8, 1x, a5, 3i8)') i, types_names(i), types_ordinal(i), types_charge(1),  types_isotope(i)
enddo
!
write(*,*)
write(*,'(a,a)') 'Atom coordinates: ', coordinate_unit
!
write(*,*)
write(*,'(a, i8)') ' Number of atoms      ', number_of_atoms
write(*,'(a)') '  Number  Name    Type               Position                       Unit cell         Site  Property'
i=1
write(*,'(i8, a5, i8, 3f12.6, 5i8)') i, types_names(atom_type(i)), atom_type(i), atom_pos(:,i), atom_unit_cell(:,i), atom_site(i), atom_property(i)
i=number_of_atoms
write(*,'(i8, a5, i8, 3f12.6, 5i8)') i, types_names(atom_type(i)), atom_type(i), atom_pos(:,i), atom_unit_cell(:,i), atom_site(i), atom_property(i)
write(*,*)
!
if(l_anisotropic_adp) then
write(*,'(a,i8)') ' Number of ADPs ', anisotropic_adp%anisotropic_n_type
write(*,'(a   )') '   Number ISO    U11         U22         U33         U23         U13         U11         Ueqv'
do i= 1, anisotropic_adp%anisotropic_n_type
   write(*,'(i8, l3, 7f12.6)') i, anisotropic_adp%anisotropic_is_iso(i), anisotropic_adp%anisotropic_adp(:,i)
enddo
else
   write(*,'(a)') ' No anisotropic ADP are present '
endif
write(*,*)
!
if(l_average_struc) then
write(*,'(a,i8)') ' Number of Atoms in average structure ', average_struc%aver_n_atoms
do i=1, average_struc%aver_n_atoms
   write(*,'(i8, i8, 11f12.6,i8)') i, average_struc%atom_type(i), average_struc%position(:,i), &
   average_struc%occupancy(i), average_struc%anis_adp(:,i), average_struc%site_number(i)
enddo
else
   write(*,'(a)') ' No average structure info '
endif
write(*,*)
!
if(l_molecules) then
write(*,'(a,i8)') ' Number of  molecules', molecules%number_moles
write(*,'(a,i8)') '   Number   Type  Character N_atoms   Biso       Corrlin    Corrquad'
write(*,'(a,i8)') '   Atoms in the molecule '
  i=1
  write(*,'(4i8, 3f12.6)') i, molecules%mole_int(:,i), molecules%mole_real(:,i)
  write(*,*) molecules%atom_index(:,i)
  i=molecules%number_moles
else
   write(*,'(a)') ' No molecules present'
endif
write(*,*)
!
write(*,*) ' Magnetic spins '
if(l_magnetic_spins) then
   write(*,'(i8, 3f12.6)') i, magnetic_spins(:,1)
else
   write(*,'(a)') ' No magnetic spins present'
endif
write(*,*)
write(*,*)
!
end subroutine unified_dump_structure
!
!*******************************************************************************
!
subroutine unified_dump_data(unit_cell_lengths, unit_cell_angles,           &
                             symmetry_H_M, symmetry_origin, symmetry_abc, symmetry_n_mat, &
                             symmetry_mat,                                                &
                             data_type_experiment, &
                             data_type_style     , &
                             data_type_axes      , &
                             data_type_content   , &
                             data_type_reciprocal, &
                             data_type_with_bragg, &
                             data_type_symmetrized, &
                             data_type_number     , &
                             data_rad_radiation , data_rad_symbol, data_rad_length, &
                             data_dimension  , &
                             data_abs_is_hkl     , &
                             data_ord_is_hkl     , &
                             data_top_is_hkl     , &
                             coordinate_unit , &
                             data_corner     , &
                             data_vector     , &
                             data_values     , &
                             NMSG, ier_num, ier_msg,                                      &
                                            crystal_meta,                  &
                             data_imag              &
                            ) 
!-
!   Dump partion data info
!+
!
!use lib_nx_transfer_mod
use lib_unified_chars_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP)        , dimension(3)                 , intent(in ) :: unit_cell_lengths  ! Unit cell parameters a, b, c
real(kind=PREC_DP)        , dimension(3)                 , intent(in ) :: unit_cell_angles   ! Unit cell angles alpha, beta, gammy
character(len=32)                                        , intent(in ) :: symmetry_H_M       ! Hermann-Mauguin Symbol
integer                                                  , intent(in ) :: symmetry_origin    ! Space group origin 1 / 2
character(len=3)                                         , intent(in ) :: symmetry_abc       ! Permutation for orthorhombic space groups
integer                                                  , intent(in ) :: symmetry_n_mat     ! Number of symmetry elements
real(kind=PREC_DP)        , dimension(3,4,symmetry_n_mat), intent(in ) :: symmetry_mat       ! Actual Symmetry matrices
!
character(len=*)                                         , intent(in) :: data_type_experiment ! 'experimental';  'calculated'
character(len=*)                                         , intent(in) :: data_type_style      ! 'powder_diffraction', 'powder_pdf', 'single_diffraction', 'single_pdf' ...
character(len=*)                                         , intent(in) :: data_type_axes       ! 'hkl', 'Q', 'theta', 'theta', 'uvw', 'xyz', 'r' ...
character(len=*)                                         , intent(in) :: data_type_content    ! 'intensity', '3d-delta-pdf', 'amplitide', 'density' ...
character(len=*)                                         , intent(in) :: data_type_reciprocal ! 'reciprocal', 'patterson', 'direct'
character(len=*)                                         , intent(in) :: data_type_with_bragg ! 'bragg_present'; 'bragg_subtracted'
character(len=*)                                         , intent(in) :: data_type_symmetrized! 'none'; 'point'; 'laue'; 'space', 'sym_elem'
character(len=*)                                         , intent(in) :: data_type_number     ! 'none'; 'point'; 'laue'; 'space', 'sym_elem'
character(len=*)                                         , intent(in) :: data_rad_radiation   ! 'xray'; 'neutron'; 'electron'
character(len=*)                                         , intent(in) :: data_rad_symbol      ! CU, CUA1, CU12
real(kind=PREC_DP)        , dimension(3)                 , intent(in) :: data_rad_length      ! Numerical value
integer                   , dimension(3)                 , intent(in) :: data_dimension   ! Data have dimensions along (HKL) / (UVW)
integer                                                  , intent(in) :: data_abs_is_hkl      ! Abscissa is 1=h 2=k 3=l
integer                                                  , intent(in) :: data_ord_is_hkl      ! Ordinate is 1=h 2=k 3=l
integer                                                  , intent(in) :: data_top_is_hkl      ! top-axis is 1=h 2=k 3=l
character(len=32)                                        , intent(in) :: coordinate_unit    ! Units for coordinates 'basecell_fractional', isupercell_fractional'
real(kind=PREC_DP)        , dimension(3   )              , intent(in) :: data_corner          ! Lower left bottom corner in fractional coordinates
real(kind=PREC_DP)        , dimension(3, 3)              , intent(in) :: data_vector          ! Increment vectors abs: (:,1); ord: (:,2); top: (:,3)
real(kind=PREC_DP)        , dimension(data_dimension(1), &
                                      data_dimension(2), &
                                      data_dimension(3)), intent(in) :: data_values  ! Actual data array

!
integer                                                  , intent(in ) :: NMSG               ! Dimension of ier_msg
integer                                                  , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*)          , dimension(NMSG)              , intent(out) :: ier_msg            ! Error messages
character(len=PREC_STRING), dimension(  5)               , intent(in ) :: crystal_meta       ! Metadata, see "c_meta"
real(kind=PREC_DP)        , dimension(:,:,:), allocatable, intent(in ) :: data_imag
!
integer :: i   ! Dummy loop index
!
write(*,*)
write(*,'(a)') '      Partial dump of unified data '
write(*,*)
write(*,'(a)') '   Crystal Meta information'
do i= 1, N_META
   write(*,'(a26, 1x, a)') c_meta(i), crystal_meta(i)(1:len_trim(crystal_meta(i)))
enddo
write(*,*)
!
write(*,'(a,3f12.6)')  ' Unit cell params ', unit_cell_lengths
write(*,'(a,3f12.6)')  ' Unit cell angles ', unit_cell_angles
!
write(*,'(a, a     )') ' Symmetry H_M     ', symmetry_H_M
write(*,'(a, i5    )') ' Symmetry Origin  ', symmetry_origin
write(*,'(a, a     )') ' Symmetry abc     ', symmetry_abc
write(*,'(a, i5    )') ' Symmetry n MAT   ', symmetry_n_mat
write(*,'(a, 4f11.6)') ' Symmetry 1 fst   ', symmetry_mat(1,:,1)
write(*,'(a, 4f11.6)') ' Symmetry 1 scnd  ', symmetry_mat(2,:,1)
write(*,'(a, 4f11.6)') ' Symmetry 1 lst   ', symmetry_mat(3,:,1)
if(symmetry_n_mat>10) then
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'fst  ',symmetry_mat(1,:,symmetry_n_mat/2)
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'scnd ',symmetry_mat(2,:,symmetry_n_mat/2)
write(*,'(a, i3, a5, 4f11.6)') ' Symmetry ', symmetry_n_mat/2, 'lst  ',symmetry_mat(3,:,symmetry_n_mat/2)
endif
write(*,'(a, 4f11.6)') ' Symmetry N fst   ', symmetry_mat(1,:,symmetry_n_mat)
write(*,'(a, 4f11.6)') ' Symmetry N scnd  ', symmetry_mat(2,:,symmetry_n_mat)
write(*,'(a, 4f11.6)') ' Symmetry N lst   ', symmetry_mat(3,:,symmetry_n_mat)
write(*,*)
write(*,'(a, a     )') ' Data experiment  ', data_type_experiment
write(*,'(a, a     )') ' Data with  bragg ', data_type_with_bragg
write(*,'(a, a     )') ' Data symmetrized ', data_type_symmetrized
write(*,'(a, a     )') ' Data number      ', data_type_number
write(*,'(a, a     )') ' Data axes        ', data_type_axes
write(*,'(a, a     )') ' Data radiation   ', data_rad_radiation
write(*,'(a, 3i5   )') ' Data dimensions  ', data_dimension
write(*,'(a, 3i5   )') ' Data abs,ord,top ', data_abs_is_hkl, data_ord_is_hkl, data_top_is_hkl
write(*,'(a, a     )') ' Data style       ', data_type_style
write(*,'(a, a     )') ' Data content     ', data_type_content
write(*,'(a, a     )') ' Data symbol      ', data_rad_symbol
write(*,'(a, 3f11.6)') ' Data length      ', data_rad_length
!
write(*,*)
write(*,'(a,a)')       'Axes coordinates: ', coordinate_unit
write(*,*)
write(*,'(a, 3f11.6)') ' Data corner      ', data_corner(:)
write(*,*)
write(*,'(a, 3f11.6)') ' Data abscissa    ', data_vector(:,1)
write(*,'(a, 3f11.6)') ' Data ordinate    ', data_vector(:,2)
write(*,'(a, 3f11.6)') ' Data top_axis    ', data_vector(:,3)
write(*,*)
write(*,*) ' DATA_VALUES', data_values(1:2,1,1)
write(*,*) ' DATA_MINMAX', minval(data_values), maxval(data_values)
if(data_type_number == 'complex' .and. allocated(data_imag)) then
write(*,*)
write(*,*) ' DATA_IMAG  ', data_imag  (1:2,1,1)
write(*,*) ' DATA_MINMAX', minval(data_imag  ), maxval(data_imag  )
endif
write(*,*)
!
end subroutine unified_dump_data
!
!*******************************************************************************
!
end module unified_dump_mod
