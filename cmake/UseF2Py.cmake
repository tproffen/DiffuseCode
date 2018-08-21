# +-----------------------------------------------------------------------------+
# |   Copyright (C) 2011-2015                                                   |
# |   Original by Marcel Loose (loose <at> astron.nl) 2011-2013                 |
# |   Modified by Chris Kerr (chris.kerr <at> mykolab.ch) 2013-2015             |
# |                                                                             |
# |   This program is free software; you can redistribute it and/or modify      |
# |   it under the terms of the GNU General Public License as published by      |
# |   the Free Software Foundation; either version 2 of the License, or         |
# |   (at your option) any later version.                                       |
# |                                                                             |
# |   This program is distributed in the hope that it will be useful,           |
# |   but WITHOUT ANY WARRANTY; without even the implied warranty of            |
# |   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             |
# |   GNU General Public License for more details.                              |
# |                                                                             |
# |   You should have received a copy of the GNU General Public License         |
# |   along with this program; if not, write to the                             |
# |   Free Software Foundation, Inc.,                                           |
# |   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 |
# +-----------------------------------------------------------------------------+

find_package(PythonInterp REQUIRED)

if (NOT F2PY_SUFFIX)
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or sysconfig.get_config_var('SO'))"
                  OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
                  RESULT_VARIABLE FOUND_PYTHON_EXT_SUFFIX)
  if (NOT ${FOUND_PYTHON_EXT_SUFFIX} EQUAL 0)
    set (F2PY_SUFFIX "" CACHE STRING "Suffix added by F2Py to the module name to get the output file name." )
    message(FATAL_ERROR "Unable to determine file extension of compiled Python modules - specify it with F2PY_SUFFIX")
  endif (NOT ${FOUND_PYTHON_EXT_SUFFIX} EQUAL 0)
  set (F2PY_SUFFIX ${PYTHON_EXT_SUFFIX} CACHE INTERNAL "the F2PY extension")
endif (NOT F2PY_SUFFIX)

## Path to the f2py executable
find_program(F2PY_EXECUTABLE NAMES "f2py${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
                                   "f2py-${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}"
                                   "f2py${PYTHON_VERSION_MAJOR}"
                                   "f2py"
                             REQUIRED)


## -----------------------------------------------------------------------------
## Macro to generate a Python interface module from one or more Fortran sources
##
## Usage: add_f2py_module(<module-name> <src1>..<srcN> DESTINATION <install-dir>
##
macro (add_f2py_module _name)
  # Parse arguments.
  string(REGEX REPLACE ";?DESTINATION.*" "" _srcs "${ARGN}")
  string(REGEX MATCH "DESTINATION;.*" _dest_dir "${ARGN}")
  string(REGEX REPLACE "^DESTINATION;" "" _dest_dir "${_dest_dir}")

  # Sanity checks.
  if(_srcs MATCHES "^$")
    message(FATAL_ERROR "add_f2py_module: no source files specified")
  endif(_srcs MATCHES "^$")

  # Get the compiler-id and map it to compiler vendor as used by f2py.
  # Currently, we only check for GNU, but this can easily be extended. 
  # Cache the result, so that we only need to check once.
  if(NOT F2PY_FCOMPILER)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
        set(_fcompiler "gnu95")
      else(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
        set(_fcompiler "gnu")
      endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    else(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set(_fcompiler "F2PY_FCOMPILER-NOTFOUND")
    endif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(F2PY_FCOMPILER ${_fcompiler} CACHE STRING
      "F2PY: Fortran compiler type by vendor" FORCE)
    if(NOT F2PY_FCOMPILER)
      message(STATUS "[F2PY]: Could not determine Fortran compiler type. "
                     "Troubles ahead!")
    endif(NOT F2PY_FCOMPILER)
  endif(NOT F2PY_FCOMPILER)

  # Set f2py compiler options: compiler vendor and path to Fortran77/90 compiler.
  if(F2PY_FCOMPILER)
    set(_fcompiler_opts "--fcompiler=${F2PY_FCOMPILER}")
    list(APPEND _fcompiler_opts "--f77exec=${CMAKE_Fortran_COMPILER}" "--f77flags='${CMAKE_Fortran_FLAGS} -DF2PY'")
    if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
      list(APPEND _fcompiler_opts "--f90exec=${CMAKE_Fortran_COMPILER}" "--f90flags='${CMAKE_Fortran_FLAGS} -DF2PY'")
    endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  endif(F2PY_FCOMPILER)
  # Make the source filenames absolute.
  set(_abs_srcs)
  foreach(_src ${_srcs})
    get_filename_component(_abs_src ${_src} ABSOLUTE)
    list(APPEND _abs_srcs ${_abs_src})
  endforeach(_src ${_srcs})

  # Get a list of the include directories.
  # The f2py --include_paths option, used when generating a signature file,
  # needs a colon-separated list. The f2py -I option, used when compiling
  # the sources, must be repeated for every include directory.
  get_directory_property(_inc_dirs INCLUDE_DIRECTORIES)
  string(REPLACE ";" ":" _inc_paths "${_inc_dirs}")
  set(_inc_opts)
  foreach(_dir ${_inc_dirs})
    list(APPEND _inc_opts "-I${_dir}")
  endforeach(_dir)
  # Define the command to generate the Fortran to Python interface module. The
  # output will be a shared library that can be imported by python.
  if ( "${_srcs}" MATCHES "^[^;]*\\.pyf;" )
    add_custom_command(OUTPUT "${_name}${F2PY_SUFFIX}"
      COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name}
              --build-dir "${CMAKE_CURRENT_BINARY_DIR}/f2py-${_name}"
              ${_fcompiler_opts} ${_inc_opts} -c ${_abs_srcs}
      DEPENDS ${_srcs}
      COMMENT "[F2PY] Building Fortran to Python interface module ${_name}")

  else ( "${_srcs}" MATCHES "^[^;]*\\.pyf;" )
    add_custom_command(OUTPUT "${_name}${F2PY_SUFFIX}"
      COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name} -h ${_name}.pyf
              --build-dir "${CMAKE_CURRENT_BINARY_DIR}/f2py-${_name}"
              --include-paths ${_inc_paths} --overwrite-signature ${_abs_srcs}
      COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name} -c "${CMAKE_CURRENT_BINARY_DIR}/f2py-${_name}/${_name}.pyf"
              --build-dir "${CMAKE_CURRENT_BINARY_DIR}/f2py-${_name}"
              ${_fcompiler_opts} ${_inc_opts} ${_abs_srcs} ${F2PY_LIBS}
      DEPENDS ${_srcs}
      COMMENT "[F2PY] Building Fortran to Python interface module ${_name}")

  endif ( "${_srcs}" MATCHES "^[^;]*\\.pyf;" )
  


  # Add a custom target <name> to trigger the generation of the python module.
  add_custom_target(${_name} ALL DEPENDS "${_name}${F2PY_SUFFIX}")

  if(NOT (_dest_dir MATCHES "^$" OR _dest_dir MATCHES ";"))
    # Install the python module
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${_name}${F2PY_SUFFIX}"
            DESTINATION ${_dest_dir})
  endif(NOT (_dest_dir MATCHES "^$" OR _dest_dir MATCHES ";"))


endmacro (add_f2py_module)