include(UseF2Py)

find_package(NumPy  REQUIRED)
find_package(PGPLOT REQUIRED)

# Sources for Python interface 

set (SOURCES suite_python_support.f90 suite_python.f90)

# Libraries we need
		
set (F2PY_LIBS ${F2PY_LIBS} -L${DIFFUSE_BINARY_DIR}/kuplot/prog/ -lkuplot_all -L/usr/local/pgplot -lpgplot -lpng -lX11)
set (F2PY_LIBS ${F2PY_LIBS} -L${DIFFUSE_BINARY_DIR}/diffev/prog/ -ldiffev_all)
set (F2PY_LIBS ${F2PY_LIBS} -L${DIFFUSE_BINARY_DIR}/discus/prog/ -ldiscus_all)
set (F2PY_LIBS ${F2PY_LIBS} -L${DIFFUSE_BINARY_DIR}/suite/prog/ -ldiscus_suite_full)
set (F2PY_LIBS ${F2PY_LIBS} -L${DIFFUSE_BINARY_DIR}/lib_f90/ -llib_f90 -llib_f90c -lreadline)

# Build the module

add_f2py_module(suite_python ${SOURCES})

