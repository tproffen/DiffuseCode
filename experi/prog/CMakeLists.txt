# EXPERI build

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
include_directories(${DIFFUSE_SOURCE_DIR}/discus/prog)
include_directories(${DIFFUSE_BINARY_DIR}/discus/prog)
include_directories(${DIFFUSE_SOURCE_DIR}/kuplot/prog)
include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${DIFFUSE_BINARY_DIR}/discus/prog)
link_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

set (SOURCES_ALL  experi_kdo.f90 experi_loop.f90 experi_reset.f90
                  experi_setup.f90 experi_setup_sub.f90
    )

set (LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

if (DIFFUSE_SHARED_LIB)
  add_library (experi_all SHARED ${SOURCES_ALL})
else (DIFFUSE_SHARED_LIB)
  add_library (experi_all ${SOURCES_ALL})
endif (DIFFUSE_SHARED_LIB)

add_dependencies ( experi_all lib_f90 discus_all kuplot_all)

add_custom_target(experihlp
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/experi/prog/appl_exp.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/experi/prog/appl_exp.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp >
               ${DIFFUSE_BINARY_DIR}/experi/prog/experi.hlp)

add_dependencies(experi_all experihlp)
if (DIFFUSE_SHARED_LIB)
  install (TARGETS experi_all DESTINATION lib)
endif (DIFFUSE_SHARED_LIB)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/experi/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/experi)
