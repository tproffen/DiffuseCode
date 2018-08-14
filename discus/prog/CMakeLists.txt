# DISCUS Build

option(DISCUS_OMP    "Switch for compiling DISCUS with OpenMP" OFF)
if (DISCUS_OMP)
  find_package(OpenMP REQUIRED)
endif (DISCUS_OMP)

option(DISCUS_CUDA   "Switch for compiling DISCUS with CUDA"   OFF)
if (DISCUS_CUDA)
  find_package(CUDA REQUIRED)
endif (DISCUS_CUDA)

find_package(NEXUS)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)

###	
set( SOURCE_MAIN discus.f90 )
set( SOURCES_ALL  addfile.f90 discus_allocate_appl_mod.f90 atom_env_mod.f90
	atom_name.f90 discus_blk_appl.f90 bv_data_mod.f90 
	chem.f90 chem_mod.f90 chem_aver_mod.f90 chem_symm_mod.f90 discus_config_mod.f90
	discus_branch.f90 discus_export.f90 discus_nipl_header.f90
        check_blen_mod.f90 check_bound_mod.f90 check_user_prop_mod.f90
	conn_def_mod.f90 conn_mod.f90 conn_sup_mod.f90 conn_type_mod.f90 crystal_mod.f90 
	class_atom.f90 class_crystal.f90 class_internal.f90
	celltoindex.f90
	deco_mod.f90 discus_mrc.f90
	debye_mod.f90 diffuse_mod.f90
	discus_loop.f90 discus_setup.f90 domain.f90 domain_mod.f90 domaindis_mod.f90 element_data_mod.f90
	discus_err_appl.f90 discus_exit.f90 external.f90 external_mod.f90 extrmc.f90 fourier.f90 fourier_sup.f90
        discus_fit.f90 discus_pointgrp.f90 discus_reset_all.f90
        do_find_mod.f90
        four_angles_mod.f90 fourier_lmn_mod.f90
	gen_add_mod.f90 generate_mod.f90 graphic.f90 insert.f90 insert_mod.f90 intens_mod.f90
        get_iscat_mod.f90
	inter_readstru.f90
	interpret.f90 inverse_mod.f90 discus_kdo.f90 mc_mod.f90 metric.f90 micro_mod.f90
	mmc.f90 mmc_mod.f90 modify_func_mod.f90 modify_mod.f90 mole_env_mod.f90
	molecule_alloc.f90 molecule_func.f90 molecule_mod.f90 output_mod.f90 patters.f90 patters_mod.f90 
	place_molecule.f90  read_internal_mod.f90
	pdf.f90 pdf_mod.f90
	discus_plot.f90 discus_plot_mod.f90 discus_plot_init.f90 powder.f90 powder_mod.f90 powder_scat_mod.f90 powder_write_mod.f90
	powder_tables.f90 prop_char_mod.f90 prop_para_func.f90 prop_para_mod.f90 refine_mod.f90
	recipro_mod.f90 rmc.f90 rmc_mod.f90 rmc_symm.f90 rmc_sup_mod.f90 discus_save.f90 discus_save_mod.f90 shear.f90 shear_mod.f90
	discus_show.f90 spcgr_apply.f90 spcgr_mod.f90 spcgr_setup.f90 stack.f90 stack_cr_mod.f90 stack_mod.f90
        stack_rese.f90 storage.f90
	structur.f90 surface_func_mod.f90 surface_mod.f90 sym_add_mod.f90 symm.f90 symm_mod.f90 symm_sup.f90 tensors.f90
	thermal.f90 trafo.f90 transfrm.f90 trans_sup.f90 transfrm_mod.f90 trans_cart_mod.f90
	unitcell_mod.f90 discus_upd_par.f90
	update_cr_dim.f90 vtk_mod.f90
	waves.f90 waves_mod.f90 wyckoff_mod.f90 quad.f90 qval.f90 utilities.f90 zone.f90
  )
set( SOURCES_LOCAL discus_pdf_file.f90 discus_powder_file.f90 discus_out_file.f90)
#set( SOURCES ${SOURCE_MAIN} ${SOURCES_ALL} )

set( LIBS lib_f90 lib_f90c ${DIFFUSE_LIBS})

# NEXUS support ?

if (NEXUS_FOUND)
  set (SOURCES_ALL ${SOURCES_ALL} discus_nexus.f90)
  set (LIBS_OPT     ${NEXUS_LIBRARIES})
  include_directories(${NEXUS_INCLUDE_PATH})
  link_directories(${NEXUS_LIBARY_PATH})

else (NEXUS_FOUND)
  set (SOURCES_ALL ${SOURCES_ALL} discus_nexus_no.f90)
  set (LIBS_OPT )

endif (NEXUS_FOUND)

# OpenMP ?

if (CUDA_FOUND)
  message (status " building DISCUS with cuda")

  set (SOURCES_ALL ${SOURCES_ALL} four_strucf_cuda_mod.f90)
  set (CUDA_NVCC_FLAGS -O3)
  cuda_compile(FOUR_STRUCF_CUDA_O four_strucf_cuda.cu)
  set (SOURCE_MAIN ${SOURCE_MAIN} ${FOUR_STRUCF_CUDA_O})
  set (LIBS ${LIBS} cudart)

elseif( DISCUS_OMP)
  message (status " building DISCUS with openmp")

  set (SOURCES_ALL ${SOURCES_ALL} four_strucf_OMP_mod.f90)
  # Add OpenMP flag depending on compiler
  if (PGF)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
  elseif (IFORT)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
  else (PGF) # Else assume gfortran
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  endif (PGF)

else (CUDA_FOUND)

  message (status " building DISCUS serial ")
  set (SOURCES_ALL ${SOURCES_ALL} four_strucf_mod.f90)

endif (CUDA_FOUND)

add_library (discus_all ${SOURCES_ALL} )
add_dependencies (discus_all lib_f90 )

set ( LIBS discus_all ${LIBS} discus_all ${LIBS_OPT} ) 

add_executable (discus ${SOURCE_MAIN} ${SOURCES_LOCAL}) 

target_link_libraries (discus ${LIBS}) 


add_custom_target(discushlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp)

# Numpy ?

#if (NUMPY_FOUND AND  DISCUS_PYTHON)
#  add_custom_command (
#    OUTPUT lib_discus.so
#    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/prompt_mod.mod     ${DIFFUSE_BINARY_DIR}/discus/prog/
#    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/errlist_mod.mod    ${DIFFUSE_BINARY_DIR}/discus/prog/
#    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/class_macro_internal.mod      ${DIFFUSE_BINARY_DIR}/discus/prog/
#    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/param_mod.mod      ${DIFFUSE_BINARY_DIR}/discus/prog/
#    COMMAND cp ${DIFFUSE_BINARY_DIR}/lib_f90/set_sub_generic_mod.mod       ${DIFFUSE_BINARY_DIR}/discus/prog/
#    COMMAND cp ${DIFFUSE_SOURCE_DIR}/lib_f90/send_get.f90       ${DIFFUSE_BINARY_DIR}/discus/prog/
#    COMMAND cp ${DIFFUSE_SOURCE_DIR}/discus/prog/discus_py.f90  ${DIFFUSE_BINARY_DIR}/discus/prog
#    COMMAND f2py -c -m lib_discus 
#            ${DIFFUSE_LIBS}
#            ${LIBS_OPT}
#            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90c 
#            -L${DIFFUSE_BINARY_DIR}/discus/prog -ldiscus_all 
#            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90 
#            -L${DIFFUSE_BINARY_DIR}/discus/prog -ldiscus_all 
#            -L${DIFFUSE_BINARY_DIR}/lib_f90 -llib_f90c 
#            ${DIFFUSE_BINARY_DIR}/discus/prog/discus_py.f90
#            ${SOURCES_LOCAL}
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/errlist_mod.mod
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/prompt_mod.mod
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/param_mod.mod
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/set_sub_generic_mod.mod
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/class_macro_internal.mod
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/discus_py.f90
#    COMMAND rm ${DIFFUSE_BINARY_DIR}/discus/prog/send_get.f90
#    )
#  add_custom_target (discus_python DEPENDS lib_discus.so)
#  add_dependencies  (discus_python discus_all lib_f90c  )
#
#  add_dependencies(discus discushlp discus_python)
#
#else  (NUMPY_FOUND AND  DISCUS_PYTHON)

  add_dependencies(discus discushlp)

#endif (NUMPY_FOUND AND  DISCUS_PYTHON)



install (TARGETS discus DESTINATION bin)
install (FILES ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp DESTINATION share)
install (FILES ${DIFFUSE_SOURCE_DIR}/discus/prog/color.map DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/discus/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/discus)
