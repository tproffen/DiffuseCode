# The Python package's compiled part: libdiffuse_core, a plain shared library
# with bind(C) functions (python/diffuse_capi.f90 on top of the suite module in
# python/diffuse_suite.f90, declared in python/diffuse_core.h). Python loads
# it with ctypes (python/diffuse/_core.py); C, C++, Julia and others can call
# it as well. It needs no Python, NumPy or f2py at build time and does not
# depend on the Python version (migration/WP2_Binding_Decision.md, option B).

find_package(MPI)

set(CORE_LINK_LIBS discus_suite_all diffev_all refine_all experi_all discus_all kuplot_all)
set(CORE_LINK_LIBS ${CORE_LINK_LIBS} lib_f90 lib_f90c finufft_lib)
set(CORE_LINK_LIBS ${CORE_LINK_LIBS} ${DIFFUSE_LIBS} ${MPI_Fortran_LIBRARIES})
set(CORE_LINK_LIBS ${CORE_LINK_LIBS} ${HDF5_LIBRARIES})
if (DIFFUSE_GRAPHICS)
  set(CORE_LINK_LIBS ${CORE_LINK_LIBS} ${PGPLOT_PGPLOT_LIBRARY} ${X11_X11_LIB} ${PNG_LIBRARY})
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

add_library(diffuse_core SHARED
            "${CMAKE_SOURCE_DIR}/python/diffuse_suite.f90"
            "${CMAKE_SOURCE_DIR}/python/diffuse_capi.f90")
set_target_properties(diffuse_core PROPERTIES
                      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/diffuse_core_mod"
                      PUBLIC_HEADER "${CMAKE_SOURCE_DIR}/python/diffuse_core.h")
target_link_libraries(diffuse_core ${CORE_LINK_LIBS})
if (APPLE)
  set_target_properties(diffuse_core PROPERTIES INSTALL_NAME_DIR "@rpath")
endif ()

# pip install (scikit-build-core sets SKBUILD): into the diffuse package.
# A plain CMake install: lib/ and include/.
# COMPONENT goes with each artifact kind (after PUBLIC_HEADER it would cover
# only the header, and pip installs only the python component).
if (SKBUILD)
  install(TARGETS diffuse_core LIBRARY DESTINATION diffuse COMPONENT python
          PUBLIC_HEADER DESTINATION diffuse COMPONENT python)
else ()
  install(TARGETS diffuse_core LIBRARY DESTINATION lib COMPONENT python
          PUBLIC_HEADER DESTINATION include COMPONENT python)
endif ()
