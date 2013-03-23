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

add_subdirectory(lib_f90) # this must be first for date.inc
add_subdirectory(mixscat/prog)
add_subdirectory(discus/prog)
add_subdirectory(diffev/prog)
add_subdirectory(kuplot/prog)
