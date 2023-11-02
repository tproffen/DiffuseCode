# REFINE Build

option(DIFFEV_MPI    "DIFFEV_MPI"    ON)

find_package(MPI)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
include_directories(${DIFFUSE_SOURCE_DIR}/diffev/prog)
include_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
include_directories(${DIFFUSE_SOURCE_DIR}/kuplot/prog)
include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
link_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

set (SOURCES_ALL refine_add_param.f90 refine_allocate_appl.f90 refine_blk_appl.f90 
             refine_branch.f90 refine_constraint.f90 refine_control.f90
             refine_current.f90 refine_data_mod.f90 refine_do_exit.f90 
             refine_err_appl.f90 refine_fit_erg.f90 refine_fit_mache_kdo.f90 
             refine_fit_set_sub_mod.f90 refine_fix.f90 refine_kdo.f90
             refine_load.f90 refine_loop.f90 refine_loop_mpi.f90
             refine_mac_mod.f90 
             refine_random_mod.f90
             refine_reset.f90 refine_run.f90 refine_show.f90 
             refine_set.f90 refine_setup.f90 refine_setup_sub.f90 refine_set_param.f90
             refine_upd_par.f90)

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

if (MPI_Fortran_FOUND AND DIFFEV_MPI)
  message (status " building REFINE with MPI")
  set (SOURCES_ALL ${SOURCES_ALL} ${DIFFUSE_SOURCE_DIR}/diffev/prog/with_mpi.f90)
  set (LIBS_OPT ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_LIBARY_PATH})

else (MPI_Fortran_FOUND AND DIFFEV_MPI)
  message (status " building REFINE without MPI")
  set (SOURCES_ALL ${SOURCES_ALL} ${DIFFUSE_SOURCE_DIR}/diffev/prog/no_mpi.f90)

endif (MPI_Fortran_FOUND AND DIFFEV_MPI)

if (DIFFUSE_SHARED_LIB)
  add_library (refine_all SHARED ${SOURCES_ALL})
else (DIFFUSE_SHARED_LIB)
  add_library (refine_all ${SOURCES_ALL})
endif (DIFFUSE_SHARED_LIB)

add_dependencies ( refine_all lib_f90 diffev_all kuplot_all)

set (LIBS refine_all ${LIBS} refine_all ${LIBS_OPT})

add_custom_target(refinehlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/refine/prog/appl_ref.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/refine/prog/appl_ref.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/refine/prog/refine.hlp)

add_dependencies(refine_all refinehlp)

#install (TARGETS refine DESTINATION bin)
#install (FILES ${DIFFUSE_BINARY_DIR}/refine/prog/refine.hlp DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/refine/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/refine)

