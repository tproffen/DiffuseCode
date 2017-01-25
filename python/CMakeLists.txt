find_package(NumPy  REQUIRED)

if (NUMPY_FOUND AND DIFFUSE_PYTHON)

  find_package(PNG REQUIRED)
  find_package(PGPLOT REQUIRED)
  find_package(NEXUS)

# NEXUS support ?

  if (NEXUS_FOUND)
    set (SOURCES_ALL ${SOURCES_ALL} kuplot_nexus.f90)
    set (LIBS_OPT         ${NEXUS_LIBRARIES})
    include_directories(${NEXUS_INCLUDE_PATH})
    link_directories(${NEXUS_LIBARY_PATH})

  else (NEXUS_FOUND)
    set (SOURCES_ALL ${SOURCES_ALL} kuplot_nonexus.f90)
    set (LIBS_OPT )

  endif (NEXUS_FOUND)

  include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
  include_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
  include_directories(${DIFFUSE_BINARY_DIR}/discus/prog)
  include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)
  include_directories(${DIFFUSE_BINARY_DIR}/suite/prog)

  link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
  link_directories(${DIFFUSE_BINARY_DIR}/diffev/prog)
  link_directories(${DIFFUSE_BINARY_DIR}/discus/prog)
  link_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)
  link_directories(${DIFFUSE_BINARY_DIR}/suite/prog)
  link_directories(${PGPLOT_LIBRARIES})

# Operating system dependent Library paths 
  if ( CYGWIN )
    set ( LIB_PATH /usr/lib/ )
  elseif ( UNIX )
    set ( LIB_PATH /usr/lib64/ )
  endif ( CYGWIN )

  add_custom_command (
    OUTPUT lib_discus_suite.so
    COMMAND cp ${DIFFUSE_SOURCE_DIR}/python/suite_python_support.f90        ${DIFFUSE_BINARY_DIR}/python/
    COMMAND cp ${DIFFUSE_SOURCE_DIR}/python/suite_python.f90    ${DIFFUSE_BINARY_DIR}/python/
    COMMAND f2py3 -c -I${DIFFUSE_BINARY_DIR}/lib_f90/
                     -I${DIFFUSE_BINARY_DIR}/diffev/prog/
                     -I${DIFFUSE_BINARY_DIR}/discus/prog/
                     -I${DIFFUSE_BINARY_DIR}/kuplot/prog/
                     -I${DIFFUSE_BINARY_DIR}/suite/prog/
                     -I${DIFFUSE_BINARY_DIR}/python/
            -m lib_discus_suite
            -L${LIB_PATH} -lreadline
            -L${LIB_PATH} -lX11
            -L${LIB_PATH} -lpng
            ${LIBS_OPT} 
            -L${PGPLOT_DIR}     -lcpgplot
            -L${PGPLOT_DIR}     -lpgplot
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90c
            -L${DIFFUSE_BINARY_DIR}/discus/prog -ldiscus_all
            -L${DIFFUSE_BINARY_DIR}/diffev/prog -ldiffev_all
            -L${DIFFUSE_BINARY_DIR}/kuplot/prog -lkuplot_all
            -L${DIFFUSE_BINARY_DIR}/suite/prog -ldiscus_suite_all
            -L${DIFFUSE_BINARY_DIR}/suite/prog -ldiscus_suite_full
            -L${LIB_PATH} -lreadline
            -L${LIB_PATH} -lX11
            -L${LIB_PATH} -lpng
            ${LIBS_OPT} 
            -L${PGPLOT_DIR}     -lcpgplot
            -L${PGPLOT_DIR}     -lpgplot
            -L${DIFFUSE_BINARY_DIR}lib_f90 -llib_f90
            -L${DIFFUSE_BINARY_DIR}lib_f90 -llib_f90c
            -L${DIFFUSE_BINARY_DIR}discus/prog -ldiscus_all
            -L${DIFFUSE_BINARY_DIR}diffev/prog -ldiffev_all
            -L${DIFFUSE_BINARY_DIR}kuplot/prog -lkuplot_all
            -L${DIFFUSE_BINARY_DIR}suite/prog -ldiscus_suite_all
            -L${DIFFUSE_BINARY_DIR}suite/prog -ldiscus_suite_full
            ${DIFFUSE_BINARY_DIR}/python/suite_python_support.f90
            ${DIFFUSE_BINARY_DIR}/python/suite_python.f90
            -L${LIB_PATH} -lreadline
            -L${LIB_PATH} -lX11
            -L${LIB_PATH} -lpng
            ${LIBS_OPT} 
            -L${PGPLOT_DIR}     -lcpgplot
            -L${PGPLOT_DIR}     -lpgplot
            -L${DIFFUSE_BINARY_DIR}lib_f90 -llib_f90
            -L${DIFFUSE_BINARY_DIR}lib_f90 -llib_f90c
            -L${DIFFUSE_BINARY_DIR}discus/prog -ldiscus_all
            -L${DIFFUSE_BINARY_DIR}diffev/prog -ldiffev_all
            -L${DIFFUSE_BINARY_DIR}kuplot/prog -lkuplot_all
            -L${DIFFUSE_BINARY_DIR}suite/prog -ldiscus_suite_all
            -L${DIFFUSE_BINARY_DIR}suite/prog -ldiscus_suite_full
            -L${LIB_PATH} -lreadline
            -L${LIB_PATH} -lX11
            -L${LIB_PATH} -lpng
            ${LIBS_OPT} 
            -L${PGPLOT_DIR}     -lcpgplot
            -L${PGPLOT_DIR}     -lpgplot
            -L${DIFFUSE_BINARY_DIR}lib_f90 -llib_f90
            -L${DIFFUSE_BINARY_DIR}lib_f90 -llib_f90c
            -L${DIFFUSE_BINARY_DIR}discus/prog -ldiscus_all
            -L${DIFFUSE_BINARY_DIR}diffev/prog -ldiffev_all
            -L${DIFFUSE_BINARY_DIR}kuplot/prog -lkuplot_all
            -L${DIFFUSE_BINARY_DIR}suite/prog -ldiscus_suite_all
            -L${DIFFUSE_BINARY_DIR}suite/prog -ldiscus_suite_full
    COMMAND rm ${DIFFUSE_BINARY_DIR}/python/suite_python*.f90
    )
  add_custom_target (suite_python DEPENDS lib_discus_suite.so)
endif (NUMPY_FOUND AND DIFFUSE_PYTHON)
