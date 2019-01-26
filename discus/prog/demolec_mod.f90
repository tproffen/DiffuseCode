MODULE demolec_mod
!
!
INTEGER                            :: DEM_MAX_MOLETYPE   ! MAX Number of molecule types
INTEGER, DIMENSION(2)              :: dem_incl = 0       ! molecule range is included
LOGICAL, DIMENSION(:), ALLOCATABLE :: dem_lmoletype      ! molecule types are selected
!
!*******************************************************************************
!
END MODULE demolec_mod
