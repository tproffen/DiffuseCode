# Toplevel for Diffuse Project

cmake_minimum_required(VERSION 2.8.4)
project(DIFFUSE)
enable_language(Fortran)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (DIFFUSE_VERSION_MAJOR 5)
set (DIFFUSE_VERSION_MINOR 0)
set (DIFFUSE_VERSION_PATCH 0)

# Readline required for all programs

find_package(Readline REQUIRED)

# Set common libraries

set (DIFFUSE_LIBS ${READLINE_LIBRARY})

if (WIN32 OR CYGWIN)
  set (DIFFUSE_LIBS ${DIFFUSE_LIBS} ws2_32)
endif (WIN32 OR CYGWIN)


# Go into directories

add_subdirectory(lib_f90)
add_subdirectory(mixscat/prog)
add_subdirectory(discus/prog)
add_subdirectory(diffev/prog)
add_subdirectory(kuplot/prog)

add_custom_command (
  OUTPUT date.inc
  COMMAND ${DIFFUSE_SOURCE_DIR}/lib_f90/makedate "${DIFFUSE_VERSION_MAJOR}.${DIFFUSE_VERSION_MINOR}.${DIFFUSE_VERSION_PATCH}" > ${DIFFUSE_BINARY_DIR}/lib_f90/date.inc
  )
add_custom_target(dateinc DEPENDS date.inc)
add_dependencies(mixscat dateinc)
add_dependencies(discus dateinc)
add_dependencies(diffev dateinc)
add_dependencies(kuplot dateinc)

