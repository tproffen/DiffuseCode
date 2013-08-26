# File:   FindNexus.cmake
# Author: Reinhard Neder
# Date:   2013--2013
# 
# based on:
### File:   FindPGPLOT.cmake
### Author: Christian Hoffmann
### Date:   2007--2009
###
### This file is part of proprietary software of the
###   Simulation and Optimization Workgroup
###   Interdisciplinary Center for Scientific Computing (IWR)
###   University of Heidelberg, Germany.
### Copyright (C) 2007--2009 by the authors. All rights reserved.
#
####################################################################################################
#
# Find a NEXUS installation
#
# This file fulfils the CMAKE modules guidelines:
#   http://www.cmake.org/cgi-bin/viewcvs.cgi/Modules/readme.txt?root=CMake&view=markup
# PLEASE READ THERE BEFORE CHANGING SOMETHING!
#
# Defines:
#   Variable: NEXUS_FOUND        - TRUE, if the package has been completely found
#   Variable: NEXUS_INCLUDE_DIRS - List of full paths to include directories required for using NEXUS
#   Variable: NEXUS_MODULE_DIRS  - List of full paths to fortran90 module directories required for using NEXUS
#   Variable: NEXUS_LIBRARIES    - List of full paths to libraries required for using NEXUS
#   Function: USE_NEXUS()      - Convenience function for preparing CMake for usage of NEXUS.
#
###################################################################################################

INCLUDE( DefaultSearchPaths )
INCLUDE( FindPackageHandleStandardArgs )

#FIND_PACKAGE( HDF5 REQUIRED )  # I might add a search for HDF5 later on
MESSAGE( STATUS "Looking for NEXUS: " )

####################################################################################################
#### SEARCH PACKAGE COMPONENTS
####################################################################################################
SET( NEXUS_DIR "" CACHE PATH "Path to an NEXUS installation" )

SET( NEXUS_LIBRARY_SEARCH_DIRS
	${NEXUS_DIR}
	${DEFAULT_PACKAGE_DIRS}
	${DEFAULT_LIBRARY_DIRS}
)

#
# Find platform-dependend libs Just a template, if needed....
#
#IF( WIN32 OR CYGWIN )
#	MESSAGE( STATUS "Looking for NEXUS: WIN32 speacial features library" )
#	FIND_LIBRARY( GRWIN_LIBRARY
#		NAMES GrWin
#		PATHS	${NEXUS_LIBRARY_SEARCH_DIRS}
#		PATH_SUFFIXES lib  lib64 LIB
#	)
#
#	IF( GRWIN_LIBRARY )
#		SET( GRAPHIC_LIBRARIES ${GRWIN_LIBRARY} )
#		SET( GRAPHICS_LIBS_FOUND TRUE )
#	ENDIF()
#
#ELSE( WIN32 OR CYGWIN )
#	FIND_PACKAGE( X11 REQUIRED )
#	SET( GRAPHIC_LIBRARIES ${X11_LIBRARIES} )
#	SET( GRAPHICS_LIBS_FOUND TRUE )
#ENDIF()


#
# Find the components of NEXUS itself
#
MESSAGE( STATUS "Looking for NEXUS: include directory" )
FIND_PATH( NEXUS_INCLUDE_DIR napi.h napif.inc nxmodule.mod nxumodule.mod
	PATHS
		${NEXUS_DIR}
		${DEFAULT_PACKAGE_DIRS}
		${DEFAULT_INCLUDE_DIRS}
	PATH_SUFFIXES include INCLUDE
)
	MESSAGE( STATUS "Looking for NEXUS: include directory (${NEXUS_INCLUDE_DIR})" )

MESSAGE( STATUS "Looking for NEXUS: library directory " )
FIND_PATH( NEXUS_LIBRARY_DIR libNeXus.so
	PATHS
		${NEXUS_DIR}
		${DEFAULT_PACKAGE_DIRS}
		${DEFAULT_INCLUDE_DIRS}
	PATH_SUFFIXES lib lib64 LIB LIB64
)
	MESSAGE( STATUS "Looking for NEXUS: library directory (${NEXUS_LIBRARY_DIR})" )

MESSAGE( STATUS "Looking for NEXUS: nexus library" )
FIND_LIBRARY( NEXUS_NEXUS_LIBRARY
	NAMES NeXus
	PATHS	${NEXUS_LIBRARY_SEARCH_DIRS}
	PATH_SUFFIXES lib  lib64 LIB
)
MESSAGE( STATUS "Looking for NEXUS:  nexus library (${NEXUS_NEXUS_LIBRARY})" )

MESSAGE( STATUS "Looking for NEXUS: F90 library" )
FIND_LIBRARY( NEXUS_F90_LIBRARY
	NAMES NeXus90
	PATHS ${NEXUS_LIBRARY_SEARCH_DIRS}
	PATH_SUFFIXES lib  lib64 LIB
)
MESSAGE( STATUS "Looking for NEXUS: F90 library (${NEXUS_F90_LIBRARY})" )


####################################################################################################
#### EVALUATE SEARCH RESULTS
####################################################################################################
FIND_PACKAGE_HANDLE_STANDARD_ARGS( NEXUS DEFAULT_MSG
	NEXUS_NEXUS_LIBRARY
	NEXUS_F90_LIBRARY
	NEXUS_LIBRARY_DIR
	NEXUS_INCLUDE_DIR
)


IF( NEXUS_FOUND )
	SET( NEXUS_INCLUDE_DIRS ${NEXUS_INCLUDE_DIR} )
	SET( NEXUS_LIBRARIES
		"${NEXUS_NEXUS_LIBRARY}"
		"${NEXUS_F90_LIBRARY}"
	)
	SET( NEXUS_LIBRARY_PATH 
		"${NEXUS_LIBRARY_DIR}"
	)
	SET( NEXUS_INCLUDE_PATH 
		"${NEXUS_INCLUDE_DIR}"
	)

	# Function making NEXUS ready to be used.
	FUNCTION( USE_NEXUS )
		IF( NOT NEXUS_USED )
			INCLUDE_DIRECTORIES( ${NEXUS_INCLUDE_DIRS} )
			SET( NEXUS_USED TRUE )
		ENDIF()
		MESSAGE( STATUS "Using NEXUS in ${CMAKE_CURRENT_SOURCE_DIR}" )
	ENDFUNCTION( USE_NEXUS )

ENDIF()

IF( NEXUS_FOUND OR Nexus_FIND_QUIETLY )
	MARK_AS_ADVANCED(
		NEXUS_DIR
		NEXUS_INCLUDE_DIR
		NEXUS_NEXUS_LIBRARY
		NEXUS_F90_LIBRARY
	)
ENDIF()

