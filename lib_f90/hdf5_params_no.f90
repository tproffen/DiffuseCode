module lib_hdf5_params_mod
!
! copy hdf5 parameters into DISCUS parameters
! needed  just in case a user does not have HDF5
! Version without HDF5
!
!
integer, parameter :: LIB_HSIZE_T = 8   ! == HSIZE_T
integer, parameter :: LIB_HID_T   = 8   ! == HID_T
!
end module lib_hdf5_params_mod
