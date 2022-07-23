module lib_hdf5_params_mod
!
! copy hdf5 parameters into DISCUS parameters
! needed  just in case a user does not have HDF5
!
use hdf5
!
integer, parameter :: LIB_HSIZE_T = HSIZE_T
integer, parameter :: LIB_HID_T   = HID_T
!
end module lib_hdf5_params_mod
