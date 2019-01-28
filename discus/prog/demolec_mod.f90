MODULE demolec_mod
!
!
INTEGER                            :: DEM_MAX_MOLETYPE   ! MAX Number of molecule types
INTEGER                            :: DEM_MAX_ATOMTYPE   ! MAX Number of molecule types
INTEGER, DIMENSION(2)              :: dem_molerange = 0  ! molecule range is included
INTEGER, DIMENSION(2)              :: dem_atomrange = (/1,-1/)  ! atom range is included
INTEGER, DIMENSION(0:1)            :: dem_sel_prop  = (/0,0/)   ! Property selection mask
LOGICAL, DIMENSION(:), ALLOCATABLE :: dem_lmoletype      ! molecule types are selected
LOGICAL, DIMENSION(:), ALLOCATABLE :: dem_latomtype      ! atom     types are selected
!
!*******************************************************************************
!
END MODULE demolec_mod
