# DISCUS Build

option(DISCUS_NEXUS "DISCUS_NEXUS" OFF)

find_package(NEXUS)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

###	
set( SOURCES  addfile.f90 allocate_appl_mod.f90 atom_env_mod.f90
	blk_appl.f90 bv_data_mod.f90 
	chem.f90 chem_mod.f90 config_mod.f90
	conn_def_mod.f90 conn_mod.f90 crystal_mod.f90 
	class_atom.f90 class_crystal.f90 class_internal.f90
	class_dc_def.f90 
	deco_mod.f90 
	debye_mod.f90 diffuse_mod.f90
	discus.f90 domain.f90 domain_mod.f90 domaindis_mod.f90 element_data_mod.f90
	err_appl.f90 exit.f90 external.f90 external_mod.f90 extrmc.f90 fourier.f90 fourier_sup.f90
	gen_add_mod.f90 generate_mod.f90 graphic.f90 insert.f90 insert_mod.f90 intens_mod.f90
	inter_readstru.f90
	interpret.f90 inverse_mod.f90 kdo.f90 mc.f90 mc_mod.f90 metric.f90 micro_mod.f90
	mmc.f90 mmc_mod.f90 modify_func_mod.f90 modify_mod.f90 mole_env_mod.f90
	molecule_mod.f90 output_mod.f90 patters.f90 patters_mod.f90 
	place_molecule.f90  read_internal_mod.f90
	pdf.f90 pdf_mod.f90
	plot.f90 plot_mod.f90 powder.f90 powder_mod.f90 powder_scat_mod.f90 powder_write_mod.f90
	prop_para_mod.f90
	recipro_mod.f90 rmc.f90 rmc_mod.f90 save.f90 save_mod.f90 shear.f90 shear_mod.f90
	show.f90 spcgr_apply.f90 spcgr_mod.f90 spcgr_setup.f90 stack.f90 stack_cr_mod.f90 stack_mod.f90
	structur.f90 surface_mod.f90 sym_add_mod.f90 symm.f90 symm_mod.f90 tensors.f90
	thermal.f90 trafo.f90 transfrm.f90 transfrm_mod.f90 unitcell_mod.f90 upd_par.f90
	waves.f90 waves_mod.f90 wyckoff_mod.f90 
  )

set( LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

# NEXUS support ?

if (NEXUS_FOUND AND DISCUS_NEXUS)
  set (SOURCES ${SOURCES} nexus.f90)
  set (LIBS ${LIBS} ${NEXUS_LIBRARIES})
  include_directories(${NEXUS_INCLUDE_PATH})
  link_directories(${NEXUS_LIBARY_PATH})

else (NEXUS_FOUND AND DISCUS_NEXUS)
  set (SOURCES ${SOURCES} nexus_no.f90)

endif (NEXUS_FOUND AND DISCUS_NEXUS)

add_executable (discus ${SOURCES})
#target_link_libraries (discus lib_f90 lib_f90c ${DIFFUSE_LIBS}) 
target_link_libraries (discus ${LIBS}) 

add_custom_command (
  OUTPUT discus.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp
  )

add_custom_target(discushlp DEPENDS discus.hlp)
add_dependencies(discus discushlp)

install (TARGETS discus DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp DESTINATION share)
install (FILES ${DIFFUSE_SOURCE_DIR}/discus/prog/color.map DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/discus/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/discus)
