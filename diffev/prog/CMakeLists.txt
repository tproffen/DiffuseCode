# DIFFEV Build

option(DIFFEV_MPI "DIFFEV_MPI" ON)

find_package(MPI)
include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

set (SOURCES allocate_appl.f90 blk_appl.f90 compare.f90 config.f90
             constraint.f90 create_trial.f90 diff_evol.f90 diffev.f90
             err_appl.f90 exit.f90 initialise.f90 kdo.f90
             population.f90 show.f90 support_diffev_mod.f90
             run_mpi_mod.f90 triple_perm.f90 upd_par.f90)

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

if (MPI_FOUND AND DIFFEV_MPI)
  set (SOURCES ${SOURCES} with_mpi.f90)
  set (LIBS ${LIBS} ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_LIBARY_PATH})

else (MPI_FOUND AND DIFFEV_MPI)
  set (SOURCES ${SOURCES} no_mpi.f90)

endif (MPI_FOUND AND DIFFEV_MPI)

add_executable(diffev ${SOURCES})
target_link_libraries (diffev ${LIBS})

add_custom_command (
  OUTPUT diffev.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp
  )

add_custom_target(diffevhlp DEPENDS diffev.hlp)
add_dependencies(diffev diffevhlp)

install (TARGETS diffev DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/diffev/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/diffev)

