module lib_data_types_mod
!
use precision_mod
!
public 
!
integer, parameter :: H5_UNKNOWN   =  0 
integer, parameter :: H5_1D_GEN    =  1  ! Generic 1-D line
integer, parameter :: H5_2D_GEN    =  2  ! Generic 2-D plane
integer, parameter :: H5_3D_GEN    =  3  ! Generic 3-D plane
integer, parameter :: H5_1D_RECI   =  4  ! Reciprocal 1-Line
integer, parameter :: H5_2D_RECI   =  5  ! Reciprocal 2-D plane
integer, parameter :: H5_3D_RECI   =  6  ! Reciprocal 3-D plane
integer, parameter :: H5_1D_DIRECT =  7  ! Direct     1-Line
integer, parameter :: H5_2D_DIRECT =  8  ! Direct     2-D plane
integer, parameter :: H5_3D_DIRECT =  9  ! Direct     3-D plane
integer, parameter :: H5_BRAGG_I   = 10  ! Bragg data hkl, inte, sigma
integer, parameter :: H5_BRAGG_SYM = 11  ! Bragg data hkl, inte, sigma, Symmetry averaged
integer, parameter :: H5_POWDER_I  = 12  ! Powder diffraction Intensity
integer, parameter :: H5_POWDER_SQ = 13  ! Powder diffraction Intensity
integer, parameter :: H5_POWDER_FQ = 14  ! Powder diffraction Intensity
!
end module lib_data_types_mod
