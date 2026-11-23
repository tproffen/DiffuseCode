# Use the interpreter found by the top-level find_package(Python3) for the
# older Find modules below, so that NumPy headers and f2py come from the same
# Python (under pip/scikit-build-core: the isolated build environment).
if (Python3_EXECUTABLE AND NOT PYTHON_EXECUTABLE)
  set(PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
endif ()
# The suffix of Python extension modules (".cpython-312-x86_64-linux-gnu.so")
# and the include folder. These came from FindPythonExtensions, which needs
# distutils (removed in Python 3.12); sysconfig and FindPython3 give the same.
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
  OUTPUT_VARIABLE PYTHON_EXTENSION_MODULE_SUFFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT PYTHON_EXTENSION_MODULE_SUFFIX)
  message(FATAL_ERROR "Could not get the extension suffix from ${PYTHON_EXECUTABLE}")
endif ()
set(PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
find_package(NumPy)
find_package(MPI)
find_package(F2PY REQUIRED)

message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")

set(PYLINK_LIBS discus_suite_all diffev_all refine_all experi_all discus_all kuplot_all)
set(PYLINK_LIBS ${PYLINK_LIBS} lib_f90 lib_f90c finufft_lib)
set(PYLINK_LIBS ${PYLINK_LIBS} ${DIFFUSE_LIBS} ${MPI_Fortran_LIBRARIES})
set(PYLINK_LIBS ${PYLINK_LIBS} ${HDF5_LIBRARIES})
if (DIFFUSE_GRAPHICS)
  set(PYLINK_LIBS ${PYLINK_LIBS} ${PGPLOT_PGPLOT_LIBRARY} ${X11_X11_LIB} ${PNG_LIBRARY})
endif (DIFFUSE_GRAPHICS)

include_directories(${DIFFUSE_BINARY_DIR}/suite/prog/)
include_directories(${DIFFUSE_BINARY_DIR}/diffev/prog/)
include_directories(${DIFFUSE_BINARY_DIR}/discus/prog/)
include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog/)
include_directories(${DIFFUSE_BINARY_DIR}/refine/prog/)
include_directories(${DIFFUSE_BINARY_DIR}/experi/prog/)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90/)
include_directories(${DIFFUSE_BINARY_DIR}/finufft/)
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(f2py_module_name "suite_python")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/python/suite_python.f90")

set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  # --quiet: no progress report and no notes about the DISCUS modules that
  # suite_python.f90 uses internally (f2py does not need them)
  COMMAND ${F2PY_EXECUTABLE} --quiet
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

     message(">>>>>>>>>>>>>>> ${PYLINK_LIBS}")

target_link_libraries(${generated_module_file} ${PYLINK_LIBS})

target_include_directories(${generated_module_file} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS})
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

# pip install (scikit-build-core sets SKBUILD): put the module inside the
# diffuse package. A plain CMake install keeps the old place, lib/.
if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION diffuse COMPONENT python)
else ()
  install(TARGETS ${generated_module_file} DESTINATION lib COMPONENT python)
endif ()

# --- C interface: libdiffuse_core (WP2 option B, migration/WP2_Binding_Decision.md)
# A plain shared library with bind(C) functions (python/diffuse_capi.f90,
# declared in python/diffuse_core.h), for ctypes and other languages. It
# needs no Python: the Python library is taken out of the link list, so the
# library does not depend on the Python version. Built next to the f2py
# module while both exist (Fortran modules in their own folder, as
# suite_python.f90 is compiled into both).
set(CORE_LINK_LIBS ${PYLINK_LIBS})
if (Python3_LIBRARIES)
  list(REMOVE_ITEM CORE_LINK_LIBS ${Python3_LIBRARIES})
endif ()
add_library(diffuse_core SHARED
            ${fortran_src_file}
            "${CMAKE_SOURCE_DIR}/python/diffuse_capi.f90")
set_target_properties(diffuse_core PROPERTIES
                      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/diffuse_core_mod"
                      PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/python/diffuse_core.h")
target_link_libraries(diffuse_core ${CORE_LINK_LIBS})
if (APPLE)
  set_target_properties(diffuse_core PROPERTIES INSTALL_NAME_DIR "@rpath")
endif ()
# COMPONENT goes with each artifact kind (after PUBLIC_HEADER it would cover
# only the header, and pip installs only the python component).
if (SKBUILD)
  install(TARGETS diffuse_core LIBRARY DESTINATION diffuse COMPONENT python
          PUBLIC_HEADER DESTINATION diffuse COMPONENT python)
else ()
  install(TARGETS diffuse_core LIBRARY DESTINATION lib COMPONENT python
          PUBLIC_HEADER DESTINATION include COMPONENT python)
endif ()

