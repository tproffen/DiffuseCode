module lib_nx_transfer_mod
!
!  Contains types for the transfer to/from nx_write / read
!
use precision_mod
!
private 
!
type :: anis_adp_type
   integer :: anisotropic_n_type                                        ! Number of ADP type
   integer :: anisotropic_n_atom                                        ! Number of atoms in crystal
   integer           , dimension(:)  , allocatable :: anisotropic_is_iso! ADP type is isotropic(==0) or anisotropic(==1)
   real(kind=PREC_DP), dimension(:,:), allocatable :: anisotropic_adp   ! (:,i) = (U11, U22, U33, U23, U13, U12, Ueqv) for type j
   integer           , dimension(:)  , allocatable :: atom_index ! Atom i has this ADP type
end type anis_adp_type
!
type :: molecule_data
   integer                                         :: number_moles  ! Number of molecules in crystal = K
   integer                                         :: number_types  ! Number of molecule types       = J
   integer           , dimension(:,:), allocatable :: mole_int      ! (3,K) 
                                                                    ! (1, k) Molecule type of this molecule
                                                                    ! (2, k) Molecule character of this molecule 0=atoms
                                                                    ! (3, k) Number of atoms in this molecule
   real(kind=PREC_DP), dimension(:,:), allocatable :: mole_real     ! (3,J)
                                                                    ! (1, j) Uiso for this molecule type
                                                                    ! (2, j) Corrlin  for this molecule type
                                                                    ! (2, j) Corrquad for this molecule type
   integer           , dimension(:,:), allocatable :: atom_index    ! (O,K)  Molecule k contains atom 1 to o
                                                                    ! O is the maximum number of atoms per molecule
end type molecule_data
!
type :: average_structure
   integer :: aver_n_atoms
   integer           , dimension(:),   allocatable :: atom_type   ! Atom is of this type
   real(kind=PREC_DP), dimension(:,:), allocatable :: position    ! Fractional coordinates (x,y,z) for each atom
   real(kind=PREC_DP), dimension(:),   allocatable :: occupancy   ! Occupancy of this atom [0,1]
   real(kind=PREC_DP), dimension(:,:), allocatable :: anis_adp    ! (:,i) = (U11, U22, U33, U23, U13, U12, Ueqv) for type j
   integer           , dimension(:),   allocatable :: site_number ! Atom is on this site in the unit cell 
end type average_structure

!
public anis_adp_type
public molecule_data
public average_structure
!
end module lib_nx_transfer_mod
