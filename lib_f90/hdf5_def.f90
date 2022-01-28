module hdf5_def_mod
!-
!  Basic definitions for HDF5 file format YELL, DISCUS
!+
!
implicit none
!
integer, parameter :: YD_ND =10        ! Number of data sets
integer, parameter :: YD_PROGRAM        =  1
integer, parameter :: YD_data           =  2
integer, parameter :: YD_format         =  3
integer, parameter :: YD_is_direct      =  4
integer, parameter :: YD_lower_limits   =  5
integer, parameter :: YD_step_sizes     =  6
integer, parameter :: YD_unit_cell      =  7
integer, parameter :: YD_step_sizes_abs =  8
integer, parameter :: YD_step_sizes_ord =  9
integer, parameter :: YD_step_sizes_top = 10
!                                                                   ! Data set names
character(len=128), dimension(YD_ND), parameter :: yd_datasets = &
(/ 'PROGRAM       ', 'data          ', 'format        ', 'is_direct     ',      &
   'lower_limits  ', 'step_sizes    ', 'unit_cell     ', 'step_sizes_abs',      &
   'step_sizes_ord', 'step_sizes_top'                                           &
/)
!                                                                   ! Used in YELL / DISCUS
integer            , dimension(YD_ND,2), parameter :: yd_req =                  &
reshape( (/ 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, &
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), (/YD_ND, 2/))
!                                                                   ! Data set was read
logical            , dimension(YD_ND)              :: yd_present = .false.
!
end module hdf5_def_mod
