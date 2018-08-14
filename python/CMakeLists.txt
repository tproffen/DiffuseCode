# Need to clean up and define all of those at root level

find_package(NumPy  REQUIRED)
find_package(PGPLOT REQUIRED)

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

# Sources for Python interface 

set (SOURCES suite_python_support.f90 suite_python.f90)
		
add_custom_target(objs
	DEPENDS ${SOURCES}
	COMMAND ${CMAKE_Fortran_COMPILER} -fPIC -c
        -I${DIFFUSE_BINARY_DIR}/lib_f90/
        -I${DIFFUSE_BINARY_DIR}/diffev/prog/
        -I${DIFFUSE_BINARY_DIR}/discus/prog/
        -I${DIFFUSE_BINARY_DIR}/kuplot/prog/
        -I${DIFFUSE_BINARY_DIR}/suite/prog/
        -I${DIFFUSE_BINARY_DIR}/python/
        ${SOURCES}
)

add_custom_target(kindmap
	DEPENDS ${DIFFUSE_SOURCE_DIR}/python/kind_map
	COMMAND cp ${DIFFUSE_SOURCE_DIR}/python/kind_map .
)

add_custom_target(suite_python
	DEPENDS kindmap objs
	COMMAND f90wrap -k kind_map -m suite_python suite_python.f90
    COMMAND f2py-f90wrap -c -m _suite_python
		f90wrap_suite_python.f90
		suite_python_support.o suite_python.o
		-I${DIFFUSE_BINARY_DIR}/lib_f90/
		-I${DIFFUSE_BINARY_DIR}/diffev/prog/
		-I${DIFFUSE_BINARY_DIR}/discus/prog/
		-I${DIFFUSE_BINARY_DIR}/kuplot/prog/
		-I${DIFFUSE_BINARY_DIR}/suite/prog/
		-I${DIFFUSE_BINARY_DIR}/python/
		-L${DIFFUSE_BINARY_DIR}/lib_f90
		-L${DIFFUSE_BINARY_DIR}/kuplot/prog
        -L${DIFFUSE_BINARY_DIR}/diffev/prog
        -L${DIFFUSE_BINARY_DIR}/discus/prog
        -L${DIFFUSE_BINARY_DIR}/suite/prog
		-L/usr/local/pgplot/
		-L${PGPLOT_LIBRARIES}
        -lkuplot_all -ldiffev_all -ldiscus_all -ldiscus_suite_full -llib_f90 -llib_f90c
		-lreadline -lpgplot -lcpgplot -lpng -lX11 ${X11_X11_LIB} ${PGPLOT_PGPLOT_LIBRARY}
    )
