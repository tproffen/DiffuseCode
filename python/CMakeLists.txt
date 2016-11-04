if (NUMPY_FOUND AND DIFFUSE_PYTHON)
  add_custom_command (
    OUTPUT lib_discus_suite.so
    COMMAND cp ${DIFFUSE_SOURCE_DIR}/python/suite_python_support.f90        ${DIFFUSE_BINARY_DIR}/python/
    COMMAND cp ${DIFFUSE_SOURCE_DIR}/python/suite_python.f90    ${DIFFUSE_BINARY_DIR}/python/
    COMMAND f2py3 -c -I${DIFFUSE_BINARY_DIR}/lib_f90/
                     -I${DIFFUSE_BINARY_DIR}/diffev/prog/
                     -I${DIFFUSE_BINARY_DIR}/discus/prog/
                     -I${DIFFUSE_BINARY_DIR}/kuplot/prog/
                     -I${DIFFUSE_BINARY_DIR}/suite/prog/
            -m lib_discus_suite
            ${DIFFUSE_LIBS}
            ${LIBS_OPT}
            -L${DIFFUSE_BINARY_DIR}/suite/prog -ldiscus_suite_full
            -L${DIFFUSE_BINARY_DIR}/suite/prog -ldiscus_suite_all
            -L${DIFFUSE_BINARY_DIR}/kuplot/prog -lkuplot_all
            -L${DIFFUSE_BINARY_DIR}/discus/prog -ldiscus_all
            -L${DIFFUSE_BINARY_DIR}/diffev/prog -ldiffev_all
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90
            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90c
            -L${DIFFUSE_BINARY_DIR}/suite/prog -ldiscus_suite_full
            ${LIBS_OPT}
            ${DIFFUSE_LIBS}
            -L/usr/lib64/ -lNeXus90
            -L/usr/lib64/ -lX11
            -L/usr/lib64/ -lpng
            -L/usr/local/pgplot/ -lcpgplot
            -L/usr/local/pgplot/ -lpgplot
            ${DIFFUSE_BINARY_DIR}/python/suite_python_support.f90
            ${DIFFUSE_BINARY_DIR}/python/suite_python.f90
    COMMAND rm ${DIFFUSE_BINARY_DIR}/python/suite_python.f90
    )
  add_custom_target (suite_python DEPENDS lib_discus_suite.so)
  add_dependencies  (suite_python suite_full                 )
endif (NUMPY_FOUND AND DIFFUSE_PYTHON)
