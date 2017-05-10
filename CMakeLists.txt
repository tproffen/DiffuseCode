# Toplevel for Diffuse Project

cmake_minimum_required(VERSION 2.8.4)
project(DIFFUSE)
enable_language(Fortran)

# Detect if Intel or Portland Fortran compiler is used.
if (${CMAKE_Fortran_COMPILER} MATCHES "ifort")
  set(IFORT TRUE)
elseif (${CMAKE_Fortran_COMPILER} MATCHES "pgf")
  set(PGF TRUE)
endif (${CMAKE_Fortran_COMPILER} MATCHES "ifort")


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # gfortran 2008 doesn't exist until gcc 4.6
  exec_program(${CMAKE_C_COMPILER} ARGS "-dumpversion" OUTPUT_VARIABLE _gcc_version_info)
  if ( _gcc_version_info VERSION_LESS 4.6)
    message (SEND_ERROR "Too old of a version GNU  ${_gcc_version_info} (Need >=4.6).")
  endif( _gcc_version_info VERSION_LESS 4.6)
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (DIFFUSE_VERSION_MAJOR 5)
set (DIFFUSE_VERSION_MINOR 11)
set (DIFFUSE_VERSION_PATCH 0)

# Readline required for all programs

find_package(Readline REQUIRED)

# Option to build a PYTHON interface
option(DIFFUSE_PYTHON "DIFFUSE_PYTHON" OFF)
if (DIFFUSE_PYTHON)
  find_package(NumPy)
endif (DIFFUSE_PYTHON)

# Set common libraries

set (DIFFUSE_LIBS ${READLINE_LIBRARY})

##64## if (WIN32 OR CYGWIN)
##64##   set (DIFFUSE_LIBS ${DIFFUSE_LIBS} ws2_32)
##64## endif (WIN32 OR CYGWIN)

# Set Compiler Flags

set (CMAKE_CXX_FLAGS "-fPIC")
set (CMAKE_C_FLAGS "-fPIC")
set (CMAKE_Fortran_FLAGS "-fPIC -cpp -O3")
#set (CMAKE_Fortran_FLAGS "-fPIC -cpp -fcheck=all -g -fbacktrace -Warray-bounds -pedantic")

# Set different compiler flags for Portland compiler
if (PGF)
  set (CMAKE_Fortran_FLAGS "-fPIC -fast -Mbackslash -Mpreprocess")
  set (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") # To fix link flags
endif (PGF)

# Debug mode - turn on warnings and other options
option(DEBUG "DEBUG" OFF)
if (DEBUG)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall")
endif (DEBUG)

# Go into directories

add_subdirectory(lib_f90) # this must be first for date.inc
add_subdirectory(mixscat/prog)
add_subdirectory(discus/prog)
add_subdirectory(diffev/prog)
add_subdirectory(kuplot/prog)
add_subdirectory(suite/prog)
add_subdirectory(python)
