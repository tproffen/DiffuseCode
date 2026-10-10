module OUTPUT_FORMAT_mod
!-
!  Output format specifiers 
!+
!
implicit none
!
integer, parameter :: OUTPUT_FORMAT_stan   =  0
integer, parameter :: OUTPUT_FORMAT_post   =  1
integer, parameter :: OUTPUT_FORMAT_PGM    =  2
integer, parameter :: OUTPUT_FORMAT_GNU    =  3
integer, parameter :: OUTPUT_FORMAT_PPM    =  4
integer, parameter :: OUTPUT_FORMAT_POWDER =  5
integer, parameter :: OUTPUT_FORMAT_SHELX  =  6
integer, parameter :: OUTPUT_FORMAT_LIST5  =  7
integer, parameter :: OUTPUT_FORMAT_LIST9  =  8
integer, parameter :: OUTPUT_FORMAT_3D     =  9
integer, parameter :: OUTPUT_FORMAT_NEXUS  = 10
integer, parameter :: OUTPUT_FORMAT_VTK    = 11
integer, parameter :: OUTPUT_FORMAT_MRC    = 12
integer, parameter :: OUTPUT_FORMAT_HDF5   = 13
integer, parameter :: OUTPUT_FORMAT_XPLOR  = 14
integer, parameter :: OUTPUT_FORMAT_VESTA  = 15
integer, parameter :: OUTPUT_FORMAT_COMMON = 16
!
end module OUTPUT_FORMAT_mod
