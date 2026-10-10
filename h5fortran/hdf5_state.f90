module hdf5_state
!-
!  Global variable beyond h5fortran to indicate that HDF5 had been opendd once
!  
!+
private
public hdf5_set_state
public hdf5_get_state
public hdf5_final_close
!
logical :: hdf5_was_opened = .false.
!
contains
!
!*******************************************************************************
!
subroutine hdf5_set_state
!-
!  Used once at inition call to hdf5_open
!
  hdf5_was_opened = .true.
!
end subroutine hdf5_set_state
!
!*******************************************************************************
!
logical function hdf5_get_state()
!-
!  returns current value of hdf5_was_opened
!+
hdf5_get_state = hdf5_was_opened
!
end function hdf5_get_state
!
!*******************************************************************************
!
subroutine hdf5_final_close
!
! performs the final HDF5 closure, will only be performed at final program termination
!
use hdf5, only: H5close_f
!
integer :: ierr
!
call H5close_f(ierr)
!
end subroutine hdf5_final_close
!
!*******************************************************************************
!
end module hdf5_state
