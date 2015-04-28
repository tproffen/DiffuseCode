# File:   FindPGPLOT.cmake
# Author: Christian Hoffmann
# Date:   2007--2009
#
# This file is part of proprietary software of the
#   Simulation and Optimization Workgroup
#   Interdisciplinary Center for Scientific Computing (IWR)
#   University of Heidelberg, Germany.
# Copyright (C) 2007--2009 by the authors. All rights reserved.
#
####################################################################################################
#
# Find a PGPLOT installation
#
# This file fulfils the CMAKE modules guidelines:
#   http://www.cmake.org/cgi-bin/viewcvs.cgi/Modules/readme.txt?root=CMake&view=markup
# PLEASE READ THERE BEFORE CHANGING SOMETHING!
#
# Defines:
#   Variable: PGPLOT_FOUND        - TRUE, if the package has been completely found
#   Variable: PGPLOT_INCLUDE_DIRS - List of full paths to include directories required for using PGPLOT
#   Variable: PGPLOT_LIBRARIES    - List of full paths to libraries required for using PGPLOT
#   Function: USE_PGPLOT()      - Convenience function for preparing CMake for usage of PGPLOT.
#
###################################################################################################

INCLUDE( DefaultSearchPaths )
INCLUDE( FindPackageHandleStandardArgs )

FIND_PACKAGE( PNG REQUIRED )
MESSAGE( STATUS "Looking for PGPLOT: " )

####################################################################################################
#### SEARCH PACKAGE COMPONENTS
####################################################################################################
SET( PGPLOT_DIR "" CACHE PATH "Path to an PGPLOT installation" )

SET( PGPLOT_LIBRARY_SEARCH_DIRS
	${PGPLOT_DIR}
	${DEFAULT_PACKAGE_DIRS}
	${DEFAULT_LIBRARY_DIRS}
        /usr/local/pgplot
)

#MESSAGE( STATUS "Looking for PGPLOT: tcl library" )
#FIND_LIBRARY( TCL_LIBRARY
#	NAMES tcl tcl85 tcl8.5 tcl84 tcl8.4 tcl83 tcl8.3 tcl82 tcl80
#	PATHS	${PGPLOT_LIBRARY_SEARCH_DIRS}
#	PATH_SUFFIXES lib  lib64 LIB
#)
#MESSAGE( STATUS "Looking for PGPLOT: tcl library" TCL_LIBRARY )

#
# Find platform-dependend graphic libs
#
##64## IF( WIN32 OR CYGWIN )
##64## 	MESSAGE( STATUS "Looking for PGPLOT: GrWin library" )
##64## 	FIND_LIBRARY( GRWIN_LIBRARY
##64## 		NAMES GrWin
##64## 		PATHS	${PGPLOT_LIBRARY_SEARCH_DIRS}
##64## 		PATH_SUFFIXES lib  lib64 LIB
##64## 	)
##64## 	MESSAGE( STATUS "Looking for PGPLOT: GrWin library" GRWIN_LIBRARY )
##64## 
##64## #	MESSAGE( STATUS "Looking for PGPLOT: gdi32 library" )
##64## #	FIND_LIBRARY( GDI32_LIBRARY
##64## #		NAMES gdi32
##64## #		PATHS	${PGPLOT_LIBRARY_SEARCH_DIRS}
##64## #		PATH_SUFFIXES lib  lib64 LIB
##64## #	)
##64## #	MESSAGE( STATUS "Looking for PGPLOT: gdi32 library" GDI32_LIBRARY )
##64## #
##64## #	MESSAGE( STATUS "Looking for PGPLOT: th32 library" )
##64## #	FIND_LIBRARY( TH32_LIBRARY
##64## #		NAMES th32
##64## #		PATHS	${PGPLOT_LIBRARY_SEARCH_DIRS}
##64## #		PATH_SUFFIXES lib  lib64 LIB
##64## #	)
##64## #	MESSAGE( STATUS "Looking for PGPLOT: th32 library" TH32_LIBRARY )
##64## #
##64## #	MESSAGE( STATUS "Looking for PGPLOT: imagehlp library" )
##64## #	FIND_LIBRARY( IMAGEHLP_LIBRARY
##64## #		NAMES imagehlp
##64## #		PATHS	${PGPLOT_LIBRARY_SEARCH_DIRS}
##64## #		PATH_SUFFIXES lib  lib64 LIB
##64## #	)
##64## 	MESSAGE( STATUS "Looking for PGPLOT: imagehlp library" IMAGEHLP_LIBRARY )
##64## 
##64## 	IF( GRWIN_LIBRARY )
##64## 		SET( GRAPHIC_LIBRARIES ${GRWIN_LIBRARY} )
##64## 		SET( GRAPHICS_LIBS_FOUND TRUE )
##64## 	ENDIF()
##64## 
##64## ELSE( WIN32 OR CYGWIN )
	FIND_PACKAGE( X11 REQUIRED )
	SET( GRAPHIC_LIBRARIES ${X11_LIBRARIES} )
	SET( GRAPHICS_LIBS_FOUND TRUE )
##64## ENDIF()


#
# Find the components of PGPLOT itself
#
MESSAGE( STATUS "Looking for PGPLOT: include directory" )
FIND_PATH( PGPLOT_INCLUDE_DIR cpgplot.h
	PATHS
		${PGPLOT_DIR}
		${DEFAULT_PACKAGE_DIRS}
		${DEFAULT_INCLUDE_DIRS}
                /usr/local/pgplot
	PATH_SUFFIXES include INCLUDE
)
	MESSAGE( STATUS "Looking for PGPLOT: include directory (${PGPLOT_INCLUDE_DIR})" )

MESSAGE( STATUS "Looking for PGPLOT: pgplot library" )
FIND_LIBRARY( PGPLOT_PGPLOT_LIBRARY
	NAMES pgplot
	PATHS	${PGPLOT_LIBRARY_SEARCH_DIRS}
	PATH_SUFFIXES lib  lib64 LIB
)
MESSAGE( STATUS "Looking for PGPLOT: pgplot library (${PGPLOT_PGPLOT_LIBRARY})" )

MESSAGE( STATUS "Looking for PGPLOT: cpgplot library" )
FIND_LIBRARY( PGPLOT_CPGPLOT_LIBRARY
	NAMES cpgplot
	PATHS ${PGPLOT_LIBRARY_SEARCH_DIRS}
	PATH_SUFFIXES lib  lib64 LIB
)
MESSAGE( STATUS "Looking for PGPLOT: cpgplot library (${PGPLOT_CPGPLOT_LIBRARY})" )

####################################################################################################
#### EVALUATE SEARCH RESULTS
####################################################################################################
FIND_PACKAGE_HANDLE_STANDARD_ARGS( PGPLOT DEFAULT_MSG
	PGPLOT_PGPLOT_LIBRARY
	PGPLOT_INCLUDE_DIR
	PGPLOT_CPGPLOT_LIBRARY
	GRAPHICS_LIBS_FOUND
)

IF( PGPLOT_FOUND )
	SET( PGPLOT_INCLUDE_DIRS ${PGPLOT_INCLUDE_DIR} )
	SET( PGPLOT_LIBRARIES
		"${PGPLOT_PGPLOT_LIBRARY}"
		"${PGPLOT_CPGPLOT_LIBRARY}"
		"${GRAPHIC_LIBRARIES}"
		"${TCL_LIBRARY}"
	)

	# Function making PGPLOT ready to be used.
	FUNCTION( USE_PGPLOT )
		IF( NOT PGPLOT_USED )
			INCLUDE_DIRECTORIES( ${PGPLOT_INCLUDE_DIRS} )
			SET( PGPLOT_USED TRUE )
		ENDIF()
		MESSAGE( STATUS "Using PGPLOT in ${CMAKE_CURRENT_SOURCE_DIR}" )
	ENDFUNCTION( USE_PGPLOT )

ENDIF()

IF( PGPLOT_FOUND OR PGPLOT_FIND_QUIETLY )
	MARK_AS_ADVANCED(
		PGPLOT_DIR
		PGPLOT_INCLUDE_DIR
		PGPLOT_PGPLOT_LIBRARY
		PGPLOT_CPGPLOT_LIBRARY
		TCL_LIBRARY
		GRWIN_LIBRARY
		GDI32_LIBRARY
		TH32_LIBRARY
		IMAGEHLP_LIBRARY
	)
ENDIF()

