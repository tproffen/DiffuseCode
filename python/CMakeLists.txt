find_package(PythonExtensions REQUIRED)
find_package(NumPy)
find_package(MPI)
find_package(F2PY REQUIRED)

message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")

set(PYLINK_LIBS discus_suite_all diffev_all refine_all discus_all kuplot_all experi_all)
set(PYLINK_LIBS ${PYLINK_LIBS} lib_f90 lib_f90c finufft_lib)
set(PYLINK_LIBS ${PYLINK_LIBS} ${DIFFUSE_LIBS} ${MPI_Fortran_LIBRARIES})
set(PYLINK_LIBS ${PYLINK_LIBS} ${PGPLOT_PGPLOT_LIBRARY} ${X11_X11_LIB} ${PNG_LIBRARY} ${HDF5_LIBRARIES})

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
set(fortran_src_file "${CMAKE_SOURCE_DIR}/python/suite_python.f90"
	             "${CMAKE_SOURCE_DIR}/python/suite_python_mod.f90")

set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
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

install(TARGETS ${generated_module_file} DESTINATION lib)

