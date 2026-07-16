module lib_unified_chars_mod
!
implicit none
!
logical :: l_dump = .false.
!
integer, parameter :: N_C_FLAGS        = 6        ! Number of Crystal flags
integer, parameter :: N_META           = 5        ! Number of meta entries
integer, parameter :: N_FIELDS_DATA    = 23       ! Number of entries in HDF5 file  scattering data
integer, parameter :: N_FIELDS_CRYSTAL = 36       ! Number of entries in HDF5 file  crystal structure
!
character(len=24), dimension(N_C_FLAGS       ), parameter :: c_flags  &
       =    (/'is_super_structure      ', &
              'is_asymmetric_unit      ', &
              'is_periodic_x           ', &
              'is_periodic_y           ', &
              'is_periodic_z           ', &
              'is_homogeneous          '  &
             /)
!
character(len=26), dimension(N_META          ), parameter :: c_meta &
       =    (/'audit_conform_dict_name   ' , &
              'audit_conform_dict_version' , &
              'audit_creation_date       ' , &
              'audit_creation_method     ' , &
              'audit_author_name         '   &
             /)
!
character(len=32), dimension(N_FIELDS_CRYSTAL), parameter :: c_fields_crystal &
       =    (/ &
              'unit_cell_lengths              ', &
              'unit_cell_angles               ', &
              'symmetry_space_group_name_H-M  ', &
              'space_group_origin             ', &
              'symmetry_space_group_abc       ', &
              'space_group_symop_number       ', &
              'space_group_symop_operation_mat', &
              'unit_cells                     ', &
              'number_of_types                ', &
              'number_of_atoms                ', &
              'types_names                    ', &
              'types_ordinal                  ', &
              'types_charge                   ', &
              'types_isotope                  ', &
              'atom_type                      ', &
              'atom_position                  ', &
              'atom_unit_cell                 ', &
              'atom_site_number               ', &
              'average_number                 ', &
              'average_type                   ', &
              'average_pos                    ', &
              'average_occ                    ', &
              'average_adp                    ', &
              'average_site                   ', &
              'occupancy                      ', &
              'property_flags                 ', &
              'molecules_number               ', &
              'molecules_types                ', &
              'molecules_int                  ', &
              'molecules_real                 ', &
              'molecules_index                ', &
              'magnetic_spins                 ', &
              'anisotropic_number             ', &
              'anisotropic_is_iso             ', &
              'anisotropic_adp                ', &
              'anisotropic_index              '  &
             /)
!
character(len=32), dimension(N_FIELDS_DATA   ), parameter :: c_fields_data   &
       =    (/ &
              'unit_cell_lengths              ', &
              'unit_cell_angles               ', &
              'symmetry_space_group_name_H-M  ', &
              'space_group_origin             ', &
              'symmetry_space_group_abc       ', &
              'space_group_symop_number       ', &
              'space_group_symop_operation_mat', &
              'data_type_experiment           ', &
              'data_type_style                ', &
              'data_type_axes                 ', &
              'data_type_content              ', &
              'data_type_reciprocal           ', &
              'data_type_with_bragg           ', &
              'data_type_symmetrized          ', &
              'data_type_number               ', &
              'data_radiation                 ', &
              'data_rad_symbol                ', &
              'data_rad_length                ', &
              'data_dimension                 ', &
              'data_axes                      ', &
              'data_corner                    ', &
              'data_increment_vector          ', &
              'data_values                    '  &
             /)
!
character(len=32), dimension(2               ), parameter :: c_dictionary_names    &  !(1=structure, 2 = data)
                    = (/'Disorder structure   '   , &
                        'Disorder unified data'  &
                       /)
!
logical          , dimension(N_C_FLAGS       ) :: l_flags
logical          , dimension(N_META          ) :: l_meta
logical          , dimension(N_FIELDS_DATA   ) :: l_fields_data
logical          , dimension(N_FIELDS_CRYSTAL) :: l_fields_crystal
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
integer, parameter :: INDX_unit_cell_lengths               =  1
integer, parameter :: INDX_unit_cell_angles                =  2
integer, parameter :: INDX_symmetry_space_group_name_H_M   =  3
integer, parameter :: INDX_space_group_origin              =  4
integer, parameter :: INDX_symmetry_space_group_abc        =  5
integer, parameter :: INDX_space_group_symop_number        =  6
integer, parameter :: INDX_space_group_symop_operation_mat =  7
integer, parameter :: INDX_data_type_experiment            =  8
integer, parameter :: INDX_data_type_style                 =  9
integer, parameter :: INDX_data_type_axes                  = 10
integer, parameter :: INDX_data_type_content               = 11
integer, parameter :: INDX_data_type_reciprocal            = 12
integer, parameter :: INDX_data_type_with_bragg            = 13
integer, parameter :: INDX_data_type_symmetrized           = 14
integer, parameter :: INDX_data_type_number                = 15
integer, parameter :: INDX_data_radiation                  = 16
integer, parameter :: INDX_data_rad_symbol                 = 17
integer, parameter :: INDX_data_rad_length                 = 18
integer, parameter :: INDX_data_dimension                  = 19
integer, parameter :: INDX_data_axes                       = 20
integer, parameter :: INDX_data_corner                     = 21
integer, parameter :: INDX_data_increment_vector           = 22
integer, parameter :: INDX_data_values                     = 23
!
!integer, parameter :: INDX_unit_cell_lengths              =  1
!integer, parameter :: INDX_unit_cell_angles               =  2
!integer, parameter :: INDX_symmetry_space_group_name_H-M  =  3
!integer, parameter :: INDX_space_group_origin             =  4
!integer, parameter :: INDX_symmetry_space_group_abc       =  5
!integer, parameter :: INDX_space_group_symop_number       =  6
!integer, parameter :: INDX_space_group_symop_operation_mat=  7 
integer, parameter :: INDX_unit_cells                     =  8
integer, parameter :: INDX_number_of_types                =  9
integer, parameter :: INDX_number_of_atoms                = 10
integer, parameter :: INDX_types_names                    = 11
integer, parameter :: INDX_types_ordinal                  = 12
integer, parameter :: INDX_types_charge                   = 13
integer, parameter :: INDX_types_isotope                  = 14
integer, parameter :: INDX_atom_type                      = 15
integer, parameter :: INDX_atom_position                  = 16
integer, parameter :: INDX_atom_unit_cell                 = 17
integer, parameter :: INDX_atom_site_number               = 18
integer, parameter :: INDX_average_number                 = 19
integer, parameter :: INDX_average_type                   = 20
integer, parameter :: INDX_average_pos                    = 21
integer, parameter :: INDX_average_occ                    = 22
integer, parameter :: INDX_average_adp                    = 23
integer, parameter :: INDX_average_site                   = 24
integer, parameter :: INDX_occupancy                      = 25
integer, parameter :: INDX_property_flags                 = 26
integer, parameter :: INDX_molecules_number               = 27
integer, parameter :: INDX_molecules_types                = 28
integer, parameter :: INDX_molecules_int                  = 29
integer, parameter :: INDX_molecules_real                 = 30
integer, parameter :: INDX_molecules_index                = 31
integer, parameter :: INDX_magnetic_spins                 = 32
integer, parameter :: INDX_anisotropic_number             = 33
integer, parameter :: INDX_anisotropic_is_iso             = 34
integer, parameter :: INDX_anisotropic_adp                = 35
integer, parameter :: INDX_anisotropic_index              = 36
!
end module lib_unified_chars_mod
