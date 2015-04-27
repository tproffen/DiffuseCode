# DISCUS_SUITE Build

option(SUITE_MPI    "SUITE_MPI"    ON)

find_package(MPI)
find_package(NEXUS)
find_package(PGPLOT REQUIRED)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/suite/prog)
include_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)
include_directories(${DIFFUSE_BINARY_DIR}/discus/prog)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${DIFFUSE_BINARY_DIR}/discus/prog)
link_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
link_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

set (SOURCE_MAIN suite.f90)
set (SOURCES_ALL suite_setup_mod.f90 suite_branch.f90 suite_loop.f90
	suite_err_appl.f90 suite_exit.f90 suite_execute_cost.f90 
	suite_init_mod.f90 suite_kdo.f90 
	suite_upd_par.f90 )

set (LIBS diffev_all discus_all kuplot_all lib_f90 lib_f90c ${DIFFUSE_LIBS} ${PNG_LIBRARY} )

if (MPI_FOUND AND SUITE_MPI)
#  set (SOURCES_ALL ${SOURCES_ALL} with_mpi.f90)
  set (LIBS_OPT ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_LIBARY_PATH})

else (MPI_FOUND AND SUITE_MPI)
  set (SOURCES_ALL ${SOURCES_ALL} ../../diffev/prog/no_mpi.f90)

endif (MPI_FOUND AND SUITE_MPI)

# NEXUS support ?

if (NEXUS_FOUND)
#  set (SOURCES_ALL ${SOURCES_ALL} nexus.f90)
  set (LIBS_OPT     ${LIBS_OPT} ${NEXUS_LIBRARIES})
  include_directories(${NEXUS_INCLUDE_PATH})
  link_directories(${NEXUS_LIBARY_PATH})

else (NEXUS_FOUND)
#  set (SOURCES_ALL ${SOURCES_ALL} nexus_no.f90)
  set (LIBS_OPT )

endif (NEXUS_FOUND)

# Get the correct PGPLOT libraries

if (WIN32 OR CYGWIN)
  set (LIBS_GRAPHIC ${PGPLOT_PGPLOT_LIBRARY} ${GRWIN_LIBRARY}
                     -mwindows -Wl,--subsystem,console)
#  set (SOURCES_ALL ${SOURCES_ALL} init_win.f90)
  set (PGSERVER ${PGPLOT_DIR}/GrWnd.exe)

else (WIN32 OR CYGWIN)
  set (LIBS_GRAPHIC ${X11_X11_LIB} ${PGPLOT_PGPLOT_LIBRARY})
#  set (SOURCES_ALL ${SOURCES_ALL} init_x11.f90)
  set (PGSERVER ${PGPLOT_DIR}/pgxwin_server)
endif (WIN32 OR CYGWIN)

set (LIBS_GRAPHIC ${LIBS_GRAPHIC} ${PNG_LIBRARY})


add_library (discus_suite_all ${SOURCES_ALL})
add_dependencies ( discus_suite_all diffev_all discus_all kuplot_all lib_f90)

set (LIBS discus_suite_all diffev_all discus_all kuplot_all ${LIBS_GRAPHIC} ${LIBS} ${LIBS_OPT} discus_suite_all )


add_executable(discus_suite ${SOURCE_MAIN} )
target_link_libraries (discus_suite ${LIBS})

add_custom_command (
  OUTPUT suite.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/suite/prog/appl_suite.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/suite/prog/suite.hlp
  )

add_custom_target(suitehlp DEPENDS suite.hlp)

install (TARGETS discus_suite DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/suite/prog/suite.hlp DESTINATION share)
#FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/suite/prog/sysmac/*.mac")
#install (FILES ${files} DESTINATION share/suite)


