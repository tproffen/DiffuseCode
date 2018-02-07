# MIXSCAT Build

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

set (SOURCES blk_appl.f90 err_appl.f90 exit.f90 
             fit.f90 kdo.f90 lazy.f load.f90 
             mixscat.f90 save.f90 set.f90 show.f90 
             upd_par.f90 weights.f90
             config_mod.f90 mixscat_mod.f90)

add_executable(mixscat ${SOURCES})
target_link_libraries (mixscat lib_f90 lib_f90c ${DIFFUSE_LIBS}) 

add_custom_target(mixscathlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/mixscat/prog/appl_mix.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/mixscat/prog/appl_mix.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/mixscat/prog/mixscat.hlp)

add_dependencies(mixscat mixscathlp)

install (TARGETS mixscat DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/mixscat/prog/mixscat.hlp DESTINATION share)
file(GLOB files "${DIFFUSE_SOURCE_DIR}/mixscat/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/mixscat)

