# DIFFEV Build

option(DIFFEV_MPI    "DIFFEV_MPI"    ON)

find_package(MPI)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

set (SOURCE_MAIN diffev.f90 )
set (SOURCES_ALL diffev_allocate_appl.f90 diffev_blk_appl.f90 compare.f90 diffev_config.f90
             constraint.f90 create_trial.f90 diff_evol.f90 diffev_add_param.f90
             diffev_setup.f90 diffev_loop.f90 diffev_branch.f90
             diffev_err_appl.f90 diffev_do_exit.f90 diffev_execute_cost.f90 initialise.f90 diffev_kdo.f90
             diffev_refine.f90 diffev_reset.f90 diffev_random.f90 diffev_loop_mpi.f90
             diffev_distrib_mod.f90
             population.f90 diffev_show.f90 support_diffev_mod.f90
             run_mpi_mod.f90 triple_perm.f90 diffev_upd_par.f90)

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

if (MPI_FOUND AND DIFFEV_MPI)
  message (status " building DIFFEV with MPI")
  set (SOURCES_ALL ${SOURCES_ALL} with_mpi.f90)
  set (LIBS_OPT ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_LIBARY_PATH})

else (MPI_FOUND AND DIFFEV_MPI)
  message (status " building DIFFEV without MPI")
  set (SOURCES_ALL ${SOURCES_ALL} no_mpi.f90)

endif (MPI_FOUND AND DIFFEV_MPI)

add_library (diffev_all ${SOURCES_ALL})
add_dependencies ( diffev_all lib_f90)

set (LIBS diffev_all ${LIBS} diffev_all ${LIBS_OPT})

add_executable(diffev ${SOURCE_MAIN} )
target_link_libraries (diffev ${LIBS})

add_custom_target(diffevhlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp)

add_dependencies(diffev diffevhlp)

install (TARGETS diffev DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/diffev/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/diffev)

