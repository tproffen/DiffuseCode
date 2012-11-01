# KUPLOT Build

option(KUPLOT_NEXUS "KUPLOT_NEXUS" ON)

find_package(PGPLOT REQUIRED)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${PGPLOT_LIBRARIES})

set (SOURCES blk_appl.f90 calc.f90 color.f90 draw.f90 
        err_appl.f90 exit.f90 fit.f90 frame.f90 gsas.f90 
	kdo.f90 kuplot.f90 load.f90 math.f90 para.f90 
        plot.f90 save.f90 show.f90 upd_par.f90)

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS} ${PNG_LIBRARY})

# Get the correct PGPLOT libraries

if (WIN32 OR CYGWIN) 
  set (LIBS ${LIBS} ${PGPLOT_PGPLOT_LIBRARY} ${GRWIN_LIBRARY}
                     -mwindows -Wl,--subsystem,console)
  set (SOURCES ${SOURCES} init_win.f90)
  set (PGSERVER ${PGPLOT_DIR}/GrWnd.exe)

else (WIN32 OR CYGWIN) 
  set (LIBS ${LIBS} ${X11_X11_LIB} ${PGPLOT_PGPLOT_LIBRARY})
  set (SOURCES ${SOURCES} init_x11.f90)
  set (PGSERVER ${PGPLOT_DIR}/pgxwin_server)
endif (WIN32 OR CYGWIN) 

# NEXUS support ?

if (NEXUS_FOUND AND KUPLOT_NEXUS)
  set (SOURCES ${SOURCES} nexus.f90)
  set (LIBS ${LIBS} ${NEXUS_LIBRARIES})
  include_directories(${NEXUS_INCLUDE_PATH})
  link_directories(${NEXUS_LIBARY_PATH})

else (NEXUS_FOUND AND KUPLOT_NEXUS)
  set (SOURCES ${SOURCES} nonexus.f90)

endif (NEXUS_FOUND AND KUPLOT_NEXUS)

add_executable(kuplot ${SOURCES})
target_link_libraries (kuplot ${LIBS})

add_custom_command (
  OUTPUT kuplot.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/kuplot/prog/appl_kup.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/kuplot/prog/kuplot.hlp
  )

add_custom_target(kuplothlp DEPENDS kuplot.hlp)
add_dependencies(kuplot kuplothlp)

install (TARGETS kuplot DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/kuplot/prog/kuplot.hlp DESTINATION bin)
install (FILES ${PGPLOT_DIR}/grfont.dat DESTINATION bin)
install (PROGRAMS ${PGSERVER} DESTINATION bin)

FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/kuplot/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION bin/mac/kuplot)

# Install libraries needed

get_target_property(MY_BINARY_LOCATION kuplot LOCATION)
configure_file(
   "${CMAKE_SOURCE_DIR}/cmake/dependencies.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake"
   @ONLY
   )
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake")

