module lib_unified_chars_mod
!
integer, parameter :: N_C_FLAGS = 6
integer, parameter :: N_META    = 5
character(len=24), dimension(N_C_FLAGS) :: c_flags
character(len=26), dimension(N_META   ) :: c_meta
!
integer,parameter :: o_atom_property   = 1! Atom has this property flag
integer,parameter :: o_crystal_flags   = 2! Flags , see "c_flags"
integer,parameter :: o_crystal_meta    = 3! Metadata, see "c_meta"
integer,parameter :: o_anisotropic_adp = 4! Info on anisotropic ADP ==> lib_nx_transfer.f90
integer,parameter :: o_molecules       = 5! Info on molecules       ==> lib_nx_transfer.f90
integer,parameter :: o_average_struc   = 6! Info on average struct  ==> lib_nx_transfer.f90
integer,parameter :: o_magnetic_spins  = 7! Atom has these magnetic spins
integer,parameter :: o_types_occupancy = 8! This atoms type has an occupancy of value
!
data c_flags /'is_super_structure      ', &
              'is_asymmetric_unit      ', &
              'is_periodic_x           ', &
              'is_periodic_y           ', &
              'is_periodic_z           ', &
              'is_homogeneous          '  &
             /
!
data c_meta  /'audit_conform_dict_name'    , &
              'audit_conform_dict_version' , &
              'audit_creation_date'        , &
              'audit_creation_method'      , &
              'audit_author_name'            &
             /
end module lib_unified_chars_mod
