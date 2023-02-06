# File:   FindFINUFFT.cmake
# Author: Reinhard Neder
# Date:   2022--2023
#
# Copyright (C) 2022--2023 by the authors. All rights reserved.
#
####################################################################################################
#
# Find a FINUFFT installation
#
#
# Defines:
#   Variable: FINUFFT_LIB          - List of full paths to libraries required for using FINUFFT
#
###################################################################################################

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

if( ${FFTW_USE_STATIC_LIBS} )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )
endif()

INCLUDE( DefaultSearchPaths )
INCLUDE( FindPackageHandleStandardArgs )

FILE(GLOB FINUFFT_SEARCH_PATHS "/usr/local/lib" "/usr/lib" "/usr/lib/finufft" )

find_library(
   FINUFFT_LIB
   NAMES "libfinufft.so" "libfinufft.dylib"
   PATHS ${FINUFFT_SEARCH_PATHS}
   PATH_SUFFIXES "lib" "lib64" "finufft"
   NO_DEFAULT_PATH
)
#
#MESSAGE( STATUS "Looking for FINUFFT: at ${FINUFFT_SEARCH_PATHS} Names ${NAMES}")
#
#--------------------------------------- components

if (FINUFFT_LIB)
#  MESSAGE( STATUS " YES ? ${FINUFFT_LIB} ??")
   set(FINUFFT_FOUND TRUE)
   set(FINUFFT_LIB ${FINUFFT_LIB})
else ()
   set(FINUFFT_FOUND FALSE)
   MESSAGE( STATUS " NO? ${FINUFFT_LIB} ??")
endif ()

#--------------------------------------- end components

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FINUFFT
    REQUIRED_VARS FINUFFT_LIB
    HANDLE_COMPONENTS
    )

mark_as_advanced(
   FINUFFT_LIB
)

