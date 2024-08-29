# Toplevel for Diffuse Project

cmake_minimum_required(VERSION 3.12.0)
#cmake_minimum_required(VERSION 2.8.4)
project(DIFFUSE)
enable_language(Fortran)

if ( NOT DEFINED ENV{FC} )
    message(FATAL_ERROR "The environment variable FC is not set
Please set FC to the full path to fortran compiler before running CMake/ccmake
For example use: export FC=$(which gfortran) or export FC=$(which ifort) ")
else ()
   message (status " FC set to $ENV{FC} ")
endif()
set (CMAKE_Fortran_COMPILER $ENV{FC})
message ( status " TOPLEVEL COMPILER is ${CMAKE_Fortran_COMPILER}" )

# Detect if Intel or Portland Fortran compiler is used.
if (${CMAKE_Fortran_COMPILER} MATCHES "ifort")
  set(IFORT TRUE)
elseif (${CMAKE_Fortran_COMPILER} MATCHES "pgf")
  set(PGF TRUE)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(IFORT TRUE)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(ISGNU TRUE)
endif (${CMAKE_Fortran_COMPILER} MATCHES "ifort")


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # gfortran 2008 doesn't exist until gcc 4.6
  exec_program(${CMAKE_Fortran_COMPILER} ARGS "-dumpversion" OUTPUT_VARIABLE _gf90_version_info)
  if ( _gf90_version_info VERSION_LESS 4.6)
    message (SEND_ERROR "Too old of a version GNU  ${_gf90_version_info} (Need >=4.6).")
  endif( _gf90_version_info VERSION_LESS 4.6)
endif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (DIFFUSE_VERSION_MAJOR 6)
set (DIFFUSE_VERSION_MINOR 17)
set (DIFFUSE_VERSION_PATCH 02)

# Readline required for all programs

find_package(Readline REQUIRED)
find_package(X11 REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(FFTW REQUIRED)

option(DIFFUSE_SHARED_LIB "If on build shared libraries" OFF)
#
option(DISCUS_OMP    "Switch for compiling DISCUS_SUITE with OpenMP" ON)
if (DISCUS_OMP)
  if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
      set(OpenMP_C "${CMAKE_C_COMPILER}" CACHE STRING "" FORCE)
      set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
      set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
      set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
      set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(OpenMP_CXX "${CMAKE_CXX_COMPILER}" CACHE STRING "" FORCE)
    set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
    set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
    set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
    set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
  endif()
  endif()
  find_package(OpenMP REQUIRED)
endif (DISCUS_OMP)

option(DISCUS_CUDA   "Switch for compiling DISCUS_SUITE with CUDA"   OFF)
if (DISCUS_CUDA)
  find_package(CUDA REQUIRED)
endif (DISCUS_CUDA)

# Option to build a PYTHON interface
option(DIFFUSE_PYTHON "DIFFUSE_PYTHON" OFF)

# Option to build with HDF5 support 
option(DIFFUSE_HDF5 "DIFFUSE_HDF5" ON )

if ( DIFFUSE_HDF5 )
  find_package(HDF5 COMPONENTS Fortran REQUIRED)
  set (LINK_LIBS ${LINK_LIBS} "${HDF5_LIBRARIES}" )
  INCLUDE_DIRECTORIES ("${HDF5_INCLUDE_DIRS}")
  message (status " HDF5 is ON  ${DIFFUSE_HDF5} ")
# message (status " LINK_LIBS ${LINK_LIBS}")
# message (status " HDF5_LIBRARY ${HDF5_FORTRAN_SHARED_LIBRARY}" )
# message (status " HDF5_INCLUDE ${HDF5_FORTRAN_INCLUDE_PATH}" )
else ()
  message (status " HDF5 is OFF ${DIFFUSE_HDF5} ")
endif ()

# Set common libraries

set (DIFFUSE_LIBS ${READLINE_LIBRARY})
set (DIFFUSE_LIBS ${DIFFUSE_LIBS} ${FFTW_DOUBLE_LIB} ${FFTW_DOUBLE_OPENMP_LIB} )

# Set Compiler Flags

set (CMAKE_CXX_FLAGS "-fPIC -funroll-loops -march=native -fcx-limited-range -std=c++14 -fopenmp")
set (CMAKE_C_FLAGS "-fPIC")
#
if ( ISGNU )
# set (CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -std=legacy -fno-automatic")
  set (CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -ffast-math -std=legacy -ffree-line-length-400")
elseif ( IFORT )
  set (CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 ")
endif ()

# Set different compiler flags for Portland compiler
if (PGF)
  set (CMAKE_Fortran_FLAGS "-fPIC -fast -Mbackslash -Mpreprocess")
  set (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") # To fix link flags
endif (PGF)

# Debug mode - turn on warnings and other options
option(DEBUG "DEBUG" OFF)
if (DEBUG)
  if (ISGNU)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
#   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall")
#   set (CMAKE_Fortran_FLAGS "-fPIC -cpp -fcheck=all -g -Og -fbacktrace -Warray-bounds -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -std=f2008 -fimplicit-none  -pedantic -ffree-line-length-400")
#   set (CMAKE_Fortran_FLAGS "-fPIC -cpp -fcheck=all -g -Og -fbacktrace -Warray-bounds -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -Wunused-variable -Wunused-function -Wunused-value -Wunused-dummy-argument -Wmaybe-uninitialized -Wconversion -fwhole-file -fimplicit-none  -pedantic -ffree-line-length-400")
    set (CMAKE_Fortran_FLAGS "-fPIC -cpp -fcheck=all -g -Og -fbacktrace -Warray-bounds -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -Wunused-variable -Wunused-function -Wunused-value -fwhole-file -fimplicit-none -pedantic -Wuninitialized -ffree-line-length-400")
  elseif (IFORT)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
    set (CMAKE_Fortan_FLAGS "-fPIC -cpp -check=all -g -O0 -backtrace")
  endif (ISGNU)
endif (DEBUG)

# Go into directories

add_subdirectory(finufft)
add_subdirectory(lib_f90) # this must be first for date.inc
add_subdirectory(kuplot/prog)
add_subdirectory(diffev/prog)
add_subdirectory(refine/prog)
add_subdirectory(discus/prog)
add_subdirectory(experi/prog)
add_subdirectory(suite/prog)

# Add Python directory if Python option is selected
if (DIFFUSE_PYTHON)
  add_subdirectory(python)
endif (DIFFUSE_PYTHON)

FILE (GLOB manfiles "${DIFFUSE_SOURCE_DIR}/Manual/*.pdf")
install (FILES ${manfiles} DESTINATION share)
FILE (WRITE ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt    "Installation: ${CMAKE_INSTALL_PREFIX}/bin/")
if(CYGWIN)
   FILE (APPEND ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt "\nPlatform    : Cygwin")
elseif(APPLE)
   FILE (APPEND ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt "\nPlatform    : Apple")
elseif(UNIX)
   FILE (APPEND ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt "\nPlatform    : Unix")
endif(CYGWIN)
FILE (APPEND ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt "\nManual      : ${CMAKE_INSTALL_PREFIX}/share/")
FILE (APPEND ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt "\nHelpdir     : ${CMAKE_INSTALL_PREFIX}/share/")
FILE (APPEND ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt "\nVersion     : ${DIFFUSE_VERSION_MAJOR}.${DIFFUSE_VERSION_MINOR}.${DIFFUSE_VERSION_PATCH}\n")
install (FILES ${DIFFUSE_BINARY_DIR}/DiscusSuite.txt DESTINATION share)
#FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/suite/prog/sysmac/*.mac")
#install (FILES ${files} DESTINATION share/suite)

