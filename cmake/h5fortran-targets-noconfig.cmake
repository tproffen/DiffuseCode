#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "h5fortran::h5fortran" for configuration ""
set_property(TARGET h5fortran::h5fortran APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(h5fortran::h5fortran PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "Fortran"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libh5fortran.a"
  )

list(APPEND _cmake_import_check_targets h5fortran::h5fortran )
list(APPEND _cmake_import_check_files_for_h5fortran::h5fortran "${_IMPORT_PREFIX}/lib/libh5fortran.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
