# DISCUS_SUITE Build

find_package(MPI)
find_package(NEXUS)
find_package(PGPLOT REQUIRED)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_SOURCE_DIR}/diffev/prog)
include_directories(${DIFFUSE_SOURCE_DIR}/refine/prog)
include_directories(${DIFFUSE_SOURCE_DIR}/kuplot/prog)
include_directories(${DIFFUSE_SOURCE_DIR}/discus/prog)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/suite/prog)
include_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
include_directories(${DIFFUSE_BINARY_DIR}/refine/prog)
include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)
include_directories(${DIFFUSE_BINARY_DIR}/discus/prog)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${DIFFUSE_BINARY_DIR}/discus/prog)
link_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
link_directories(${DIFFUSE_BINARY_DIR}/refine/prog)
link_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)
link_directories(${PGPLOT_LIBRARIES})

set (SOURCE_MAIN suite.f90)
set (SOURCES_ALL suite_setup_mod.f90 suite_branch.f90 suite_loop.f90
	suite_err_appl.f90 suite_exit.f90 suite_execute_cost.f90 
	suite_init_mod.f90 suite_kdo.f90  suite_set_sub_mod.f90
        suite_parallel.f90 suite_upd_par.f90 suite_top.f90 )
set( SOURCES_LOCAL suite_out_file.f90 suite_pdf_file.f90 suite_powder_file.f90)

set (LIBS diffev_all refine_all discus_all kuplot_all lib_f90 lib_f90c ${DIFFUSE_LIBS} ${PNG_LIBRARY} )

if( DISCUS_OMP)
  message (status " building SUITE with openmp")

#  set (SOURCES_ALL ${SOURCES_ALL} four_strucf_OMP_mod.f90)
  # Add OpenMP flag depending on compiler
  if (PGF)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
  elseif (IFORT)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
  else (PGF) # Else assume gfortran
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  endif (PGF)
else(DISCUS_OMP)
  message (status " building SUITE without openmp")
endif(DISCUS_OMP)

if (MPI_Fortran_FOUND AND DIFFEV_MPI)
  message (status " building SUITE with MPI")
  set (LIBS_OPT ${MPI_Fortran_LIBRARIES})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  link_directories(${MPI_LIBARY_PATH})

else (MPI_Fortran_FOUND AND DIFFEV_MPI)
  message (status " building SUITE without MPI")

endif (MPI_Fortran_FOUND AND DIFFEV_MPI)

# NEXUS support ?

if (NEXUS_FOUND)
#  set (SOURCES_ALL ${SOURCES_ALL} nexus.f90)
  set (LIBS_OPT     ${LIBS_OPT} ${NEXUS_LIBRARIES})
  include_directories(${NEXUS_INCLUDE_PATH})
  link_directories(${NEXUS_LIBARY_PATH})

else (NEXUS_FOUND)
#  set (SOURCES_ALL ${SOURCES_ALL} nexus_no.f90)
#  set (LIBS_OPT )

endif (NEXUS_FOUND)

set (LIBS_GRAPHIC ${LIBS_GRAPHIC} ${PNG_LIBRARY})

# Get the correct PGPLOT libraries

##64## if (WIN32 OR CYGWIN)
##64##   set (LIBS_GRAPHIC ${PGPLOT_PGPLOT_LIBRARY} ${GRWIN_LIBRARY}
##64##                      -mwindows -Wl,--subsystem,console)
##64## #  set (SOURCES_ALL ${SOURCES_ALL} init_win.f90)
##64##   set (PGSERVER ${PGPLOT_DIR}/GrWnd.exe)

##64## else (WIN32 OR CYGWIN)
  set (LIBS_GRAPHIC ${X11_X11_LIB} ${PGPLOT_PGPLOT_LIBRARY})
#  set (SOURCES_ALL ${SOURCES_ALL} init_x11.f90)
  set (PGSERVER ${PGPLOT_DIR}/pgxwin_server)
##64## endif (WIN32 OR CYGWIN)



add_library (discus_suite_all ${SOURCES_ALL})
add_dependencies ( discus_suite_all diffev_all refine_all discus_all kuplot_all lib_f90 suitehlp)

set (LIBS ${LIBS} discus_suite_all diffev_all refine_all discus_all kuplot_all ${LIBS_GRAPHIC} ${LIBS} ${LIBS_OPT} discus_suite_all )

add_library (discus_suite_full ${SOURCES_ALL} ${SOURCES_LOCAL})
add_dependencies ( discus_suite_full  discus_suite_all diffev_all refine_all discus_all kuplot_all lib_f90  suitehlp)

add_executable(discus_suite ${SOURCE_MAIN} ${SOURCES_LOCAL} )
target_link_libraries (discus_suite ${LIBS} ${LINK_LIBS})

add_custom_target(suitehlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/suite/prog/appl_suite.hlp
               ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp
               ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp
               ${DIFFUSE_SOURCE_DIR}/kuplot/prog/appl_kup.hlp
               ${DIFFUSE_SOURCE_DIR}/refine/prog/appl_ref.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/suite/prog/appl_suite.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp   
               ${DIFFUSE_SOURCE_DIR}/diffev/prog/appl_dif.hlp
               ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp
               ${DIFFUSE_SOURCE_DIR}/kuplot/prog/appl_kup.hlp
               ${DIFFUSE_SOURCE_DIR}/refine/prog/appl_ref.hlp
               >
               ${DIFFUSE_BINARY_DIR}/suite/prog/suite.hlp)

install (TARGETS discus_suite DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/suite/prog/suite.hlp DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/suite/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/suite)


