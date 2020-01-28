# KUPLOT Build

find_package(PGPLOT REQUIRED)
find_package(NEXUS)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${PGPLOT_LIBRARIES})

set (SOURCE_MAIN kuplot.f90 )
set (SOURCES_ALL kuplot_blk_appl.f90 calc.f90 color.f90 draw.f90 
	kuplot_err_appl.f90 kuplot_exit.f90 fit.f90 frame.f90 gsas.f90 
	fit_macro.f90 fit_mache_kdo.f90 fit_params_mod.f90 fit_set_sub_mod.f90
        kuplot_diffev_plot.f90
	kuplot_kdo.f90 koordinate_mod.f90 kuplot_loop.f90 
	kuplot_setup.f90 kuplot_branch.f90 load.f90 math.f90 para.f90 
	kuplot_plot.f90 kuplot_save.f90 kuplot_show.f90 kuplot_upd_par.f90
	kuplot_config.f90 kuplot_mod.f90 kuplot_branch.f90
	kuplot_words_mod.f90 kuplot_top.f90 kuplot_2dmap_mod.f90 kuplot_2dmap.f90
	kuplot_3dmap_mod.f90  kuplot_3dmap_draw.f90)

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS} ${PNG_LIBRARY})

# Get the correct PGPLOT libraries

##64## if (WIN32 OR CYGWIN) 
##64##   set (LIBS_GRAPHIC ${PGPLOT_PGPLOT_LIBRARY} ${GRWIN_LIBRARY}
##64##                      -mwindows -Wl,--subsystem,console)
##64##   set (SOURCES_ALL ${SOURCES_ALL} init_win.f90)
##64##   set (PGSERVER ${PGPLOT_DIR}/GrWnd.exe)
##64## 
##64## else (WIN32 OR CYGWIN) 
  set (LIBS_GRAPHIC ${X11_X11_LIB} ${PGPLOT_PGPLOT_LIBRARY})
  set (SOURCES_ALL ${SOURCES_ALL} init_x11.f90)
  set (PGSERVER ${PGPLOT_DIR}/pgxwin_server)
##64## endif (WIN32 OR CYGWIN) 

set (LIBS_GRAPHIC ${LIBS_GRAPHIC} ${PNG_LIBRARY})

# NEXUS support ?

if (NEXUS_FOUND)
  set (SOURCES_ALL ${SOURCES_ALL} kuplot_nexus.f90 ../../nexus/NXmodule.f90 ../../nexus/NXUmodule.f90)
  set (LIBS_OPT         ${NEXUS_LIBRARIES})
  include_directories(${NEXUS_INCLUDE_PATH})
  link_directories(${NEXUS_LIBARY_PATH})

else (NEXUS_FOUND)
  set (SOURCES_ALL ${SOURCES_ALL} kuplot_nonexus.f90)
  set (LIBS_OPT )

endif (NEXUS_FOUND)

add_library (kuplot_all ${SOURCES_ALL})
add_dependencies (kuplot_all lib_f90)

set (LIBS kuplot_all ${LIBS} kuplot_all ${LIBS_GRAPHIC} ${LIBS_OPT} )

#add_executable(kuplot ${SOURCE_MAIN})
#target_link_libraries (kuplot ${LIBS})

add_custom_target(kuplothlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/kuplot/prog/appl_kup.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/kuplot/prog/appl_kup.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/kuplot/prog/kuplot.hlp)

#add_dependencies(kuplot kuplothlp)

#install (TARGETS kuplot DESTINATION bin)
#install (FILES ${DIFFUSE_BINARY_DIR}/kuplot/prog/kuplot.hlp DESTINATION share)
install (FILES ${PGPLOT_DIR}/grfont.dat DESTINATION share)
install (PROGRAMS ${PGSERVER} DESTINATION bin)

FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/kuplot/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/kuplot)

