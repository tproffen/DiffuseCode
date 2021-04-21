# File:   FindHDF5.cmake
# Author: Reinhard Neder
# Date:   2021
#
####################################################################################################
#
#  Find a HDF5 fortran library
#
#  Defines:
#   Variable: HDF5_FORTRAN_FOUND          - TRUE, if the package has been completely found
#   Variable: HDF5_FORTRAN_SHARED_LIBRARY - List of full paths to libraries required for using HDF5_FORTRAN
#   Variable: HDF5_FORTRAN_INCLUDE_DIR    - List of full paths to include directories required for using HDF5_FORTRAN
#
####################################################################################################
#
INCLUDE( DefaultSearchPaths )
INCLUDE( FindPackageHandleStandardArgs )
#
FILE(GLOB HDF5_LIB_SEARCH_PATHS "/usr/lib/x86_64-linux-gnu/hdf5/serial*"
                                "/lib/x86_64-linux-gnu/hdf5/serial*/"
                                "/usr/lib*"
                                "/usr/lib64*"
                                "/usr/local/Cellar/hdf5*")

FILE(GLOB HDF5_INC_SEARCH_PATHS "/usr/include/hdf5/serial*"
                                "/usr/include*"
                                "/usr/lib64/gfortran/modules*"
                                "/usr/local/Cellar/hdf5*")

FIND_LIBRARY(HDF5_FORTRAN_LIBRARY NAME hdf5_fortran PATHS ${HDF5_LIB_SEARCH_PATHS} )

FIND_PATH(HDF5_FORTRAN_INCLUDE hdf5.mod PATHS ${HDF5_INC_SEARCH_PATHS} )

####################################################################################################
#### EVALUATE SEARCH RESULTS
####################################################################################################

FIND_PACKAGE_HANDLE_STANDARD_ARGS( HDF5
	HDF5_FORTRAN_LIBRARY
	HDF5_FORTRAN_INCLUDE
)

IF(HDF5_FOUND)
	SET (HDF5_FORTRAN_INCLUDE_PATH ${HDF5_FORTRAN_INCLUDE} )
	SET (HDF5_FORTRAN_SHARED_LIBRARY ${HDF5_FORTRAN_LIBRARY})
ENDIF()
#message ( STATUS  HDF5_FORTRAN_SHARED_LIBRARY  ${HDF5_FORTRAN_LIBRARY})
#message ( STATUS  HDF5_FORTRAN_INCLUDE_PATH ${HDF5_FORTRAN_INCLUDE} )


