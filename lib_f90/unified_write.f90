module unified_write_mod
!
!   Module to read the unified structure / scattering via the h5fortran library
!
private
public unified_write_structure
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
!
end subroutine unified_write_structure
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
