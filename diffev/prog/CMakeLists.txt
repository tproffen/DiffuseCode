# DIFFEV Build

option(DIFFEV_MPI    "DIFFEV_MPI"    ON)

find_package(MPI)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

set (SOURCE_MAIN diffev.f90 )
set (SOURCES_ALL diffev_allocate_appl.f90 diffev_blk_appl.f90 compare.f90 diffev_config.f90
             constraint.f90 create_trial.f90 diff_evol.f90
             diffev_setup.f90 diffev_loop.f90 diffev_branch.f90
             diffev_err_appl.f90 diffev_do_exit.f90 diffev_execute_cost.f90 initialise.f90 diffev_kdo.f90
             population.f90 diffev_show.f90 support_diffev_mod.f90
             run_mpi_mod.f90 triple_perm.f90 diffev_upd_par.f90)

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

if (MPI_FOUND AND DIFFEV_MPI)
  set (SOURCES_ALL ${SOURCES_ALL} with_mpi.f90)
  set (LIBS_OPT ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_LIBARY_PATH})

else (MPI_FOUND AND DIFFEV_MPI)
  set (SOURCES_ALL ${SOURCES_ALL} no_mpi.f90)

endif (MPI_FOUND AND DIFFEV_MPI)

add_library (diffev_all ${SOURCES_ALL})
add_dependencies ( diffev_all lib_f90)

set (LIBS diffev_all ${LIBS} diffev_all ${LIBS_OPT})

add_executable(diffev ${SOURCE_MAIN} )
target_link_libraries (diffev ${LIBS})

add_custom_target(diffevhlp 
  DEPENDS      ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp
               ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp)

if (NUMPY_FOUND AND DIFFUSE_PYTHON)
  add_custom_command (
    OUTPUT lib_diffev.so
    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/prompt_mod.mod     ${DIFFUSE_BINARY_DIR}/diffev/prog/
    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/errlist_mod.mod    ${DIFFUSE_BINARY_DIR}/diffev/prog/
    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/class_macro_internal.mod      ${DIFFUSE_BINARY_DIR}/diffev/prog/
    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/param_mod.mod      ${DIFFUSE_BINARY_DIR}/diffev/prog/
    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/set_sub_generic_mod.mod      ${DIFFUSE_BINARY_DIR}/diffev/prog/
    COMMAND cp ${DIFFUSE_SOURCE_DIR}/lib_f90/send_get.f90       ${DIFFUSE_BINARY_DIR}/diffev/prog/
    COMMAND cp ${DIFFUSE_SOURCE_DIR}/diffev/prog/diffev_py.f90  ${DIFFUSE_BINARY_DIR}/diffev/prog
    COMMAND f2py -c -m lib_diffev 
            ${DIFFUSE_LIBS}
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90c 
            -L${DIFFUSE_BINARY_DIR}/diffev/prog -ldiffev_all 
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90 
            -L${DIFFUSE_BINARY_DIR}/diffev/prog -ldiffev_all 
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90c 
            ${LIBS_OPT}
            ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev_py.f90
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/errlist_mod.mod
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/class_macro_internal.mod
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/param_mod.mod
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/prompt_mod.mod
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/set_sub_generic_mod.mod
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev_py.f90
    COMMAND rm ${DIFFUSE_BINARY_DIR}/diffev/prog/send_get.f90
    )

  add_custom_target (diffev_python DEPENDS lib_diffev.so )
  add_dependencies  (diffev_python diffev_all lib_f90c)

  add_dependencies(diffev diffevhlp diffev_python)

else  (NUMPY_FOUND AND DIFFUSE_PYTHON)

  add_dependencies(diffev diffevhlp)

endif (NUMPY_FOUND AND DIFFUSE_PYTHON)

install (TARGETS diffev DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/diffev/prog/diffev.hlp DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/diffev/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/diffev)

