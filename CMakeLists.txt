# Toplevel for Diffuse Project

cmake_minimum_required(VERSION 2.8.4)
project(DIFFUSE)
enable_language(Fortran)

set (DIFFUSE_VERSION_MAJOR 5)
set (DIFFUSE_VERSION_MINOR 0)
set (DIFFUSE_VERSION_PATCH 0)

# Go into directories

add_subdirectory(lib_f90)
add_subdirectory(mixscat/prog)
add_subdirectory(discus/prog)
add_subdirectory(diffev/prog)

add_custom_command (
  OUTPUT date.inc
  COMMAND ${DIFFUSE_SOURCE_DIR}/lib_f90/makedate "${DIFFUSE_VERSION_MAJOR}.${DIFFUSE_VERSION_MINOR}.${DIFFUSE_VERSION_PATCH}" > ${DIFFUSE_BINARY_DIR}/lib_f90/date.inc
  )
add_custom_target(dateinc DEPENDS date.inc)
add_dependencies(mixscat dateinc)
add_dependencies(discus dateinc)
add_dependencies(diffev dateinc)

# Packaging
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "DISCUS Program Package")
set (CPACK_PACKAGE_CONTACT "Thomas Proffen, Reinhard Neder")
set (CPACK_PACKAGE_VERSION_MAJOR "${DIFFUSE_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${DIFFUSE_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${DIFFUSE_VERSION_PATCH}")

include (CPack)

