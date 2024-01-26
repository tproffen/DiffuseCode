# DISCUS Build

#option(DISCUS_OMP    "Switch for compiling DISCUS with OpenMP" OFF)
#if (DISCUS_OMP)
#  find_package(OpenMP REQUIRED)
#endif (DISCUS_OMP)
#
#option(DISCUS_CUDA   "Switch for compiling DISCUS with CUDA"   OFF)
#if (DISCUS_CUDA)
#  find_package(CUDA REQUIRED)
#endif (DISCUS_CUDA)

find_package(NEXUS)

include_directories(${DIFFUSE_SOURCE_DIR}/lib_f90)
include_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
include_directories(${DIFFUSE_SOURCE_DIR}/kuplot/prog)
include_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

link_directories(${DIFFUSE_BINARY_DIR}/finufft)
link_directories(${DIFFUSE_BINARY_DIR}/lib_f90)
link_directories(${DIFFUSE_BINARY_DIR}/kuplot/prog)

###	
set( SOURCE_MAIN discus.f90 )
set( SOURCES_ALL  addfile.f90 discus_allocate_appl_mod.f90 atom_env_mod.f90
	atom_line.f90 atom_name.f90 discus_blk_appl.f90 discus_bragg.f90 bv_data_mod.f90 
	chem.f90 chem_mod.f90 chem_aver_mod.f90 chem_symm_mod.f90 chem_multi.f90 chem_neig_multi.f90
        discus_config_mod.f90
	discus_branch.f90 discus_crystal.f90 discus_export.f90 discus_nipl_header.f90
        check_blen_mod.f90 check_bound_mod.f90 check_user_prop_mod.f90
	conn_def_mod.f90 conn_mod.f90 conn_sup_mod.f90 conn_type_mod.f90 crystal_mod.f90 
	class_atom.f90 class_crystal.f90 class_internal.f90
	celltoindex.f90
	deco_mod.f90 demolec.f90 demolec_mod.f90 discus_mrc.f90
	debye_mod.f90 diffuse_mod.f90 discus_discamb.f90 discus_estimate.f90 discus_exp2pdf.f90
        discus_exp2pdf_data.f90 discus_exp2pdf_load.f90 discus_exp2pdf_run.f90 discus_exp2pdf_supp.f90
	discus_loop.f90 discus_reduce.f90 discus_setup.f90 discus_set_sub.f90 
        discus_s_out_file.f90 discus_s_pdf_file.f90 discus_s_powder_file.f90
        domain.f90 domain_mod.f90 domaindis_mod.f90  domain_irreg_mod.f90
	discus_err_appl.f90 discus_exit.f90 external.f90 external_mod.f90 extrmc.f90 fourier.f90 fourier_sup.f90
        discus_fft.f90 discus_fit.f90 discus_pointgrp.f90 discus_prep_anis.f90 discus_reset_all.f90
        discus_powder_fft.f90 discus_symmetrize.f90 discus_3dpdf.f90 discus_xplor.f90
        do_find_mod.f90 find_top.f90
        four_angles_mod.f90 fourier_conv.f90 fourier_finufft.f90 fourier_lmn_mod.f90 fourier_reset.f90
	gen_add_mod.f90 generate_mod.f90 graphic.f90 insert.f90 insert_mod.f90 intens_mod.f90
        get_iscat_mod.f90 
	inter_readstru.f90
	interpret.f90 inverse_mod.f90 discus_kdo.f90 discus_kdo_common.f90 mc_mod.f90 metric.f90 micro_mod.f90
	mmc.f90 mmc_mod.f90 mmc_basic.f90 mmc_mole.f90 modify_func_mod.f90 modify_mod.f90 mole_env_mod.f90
	molecule_alloc.f90 molecule_func.f90 molecule_mod.f90 output_mod.f90 patters.f90 patters_mod.f90 
        phases_mod.f90 phases_set_mod.f90 phases_stack.f90
	place_molecule.f90  read_internal_mod.f90
	perioditize.f90 pdf.f90 pdf_mod.f90
	discus_plot.f90 discus_plot_mod.f90 discus_plot_init.f90  discus_plot_export.f90
	powder.f90 powder_mod.f90 powder_scat_mod.f90 powder_write_mod.f90
	powder_tables.f90 powder_pdf_hist.f90 
        prop_char_mod.f90 prop_para_func.f90 prop_para_mod.f90 refine_mod.f90
	recipro_mod.f90 rmc.f90 rmc_mod.f90 rmc_symm.f90 rmc_sup_mod.f90 discus_save.f90 
        discus_save_mod.f90 shear.f90 shear_mod.f90
	discus_show.f90 spcgr_apply.f90 spcgr_mod.f90 spcgr_setup.f90 stack.f90 stack_cr_mod.f90 stack_mod.f90
        stack_rese.f90 storage.f90
	structur.f90 surface_func_mod.f90 surface_mod.f90 sym_add_mod.f90 symm.f90 symm_mod.f90 symm_sup.f90 tensors.f90
	thermal.f90 trafo.f90 transfrm.f90 trans_sup.f90 transfrm_mod.f90 trans_cart_mod.f90
	unitcell_mod.f90 discus_upd_par.f90
	update_cr_dim.f90 vtk_mod.f90
	waves.f90 waves_mod.f90 wyckoff_mod.f90 quad.f90 qval.f90 utilities.f90 zone.f90
        private.f90
  )
#set( SOURCES_LOCAL discus_pdf_file.f90 discus_powder_file.f90 discus_out_file.f90)
#set( SOURCES ${SOURCE_MAIN} ${SOURCES_ALL} )

set( LIBS lib_f90 lib_f90c kuplot_all ${DIFFUSE_LIBS})

if (DIFFUSE_HDF5)
  if (DEFINED ENV{HDF5_INC_DIR})
    include_directories($ENV{HDF5_INC_DIR})
  endif ()
  set (SOURCES_ALL ${SOURCES_ALL} hdf_write.f90)
  message (DISCUS " HDF5 is ON  ${DIFFUSE_HDF5} ")
else (DIFFUSE_HDF5)
  set (SOURCES_ALL ${SOURCES_ALL} hdf_write_no.f90)
  message (DISCUS " HDF5 is OFF ${DIFFUSE_HDF5} ")
endif ()

# NEXUS support ?

if (NEXUS_FOUND)
  set (SOURCES_ALL ${SOURCES_ALL} discus_nexus.f90 ../../nexus/NXmodule.f90 ../../nexus/NXUmodule.f90)
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
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
  else (PGF) # Else assume gfortran
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  endif (PGF)

else (CUDA_FOUND)

  message (status " building DISCUS serial ")
  set (SOURCES_ALL ${SOURCES_ALL} four_strucf_mod.f90)

endif (CUDA_FOUND)
set (SOURCES_ALL ${SOURCES_ALL} "${DIFFUSE_SOURCE_DIR}/finufft/include/finufft_mod.f90")

if (DIFFUSE_SHARED_LIB)
  add_library (discus_all SHARED ${SOURCES_ALL})
else (DIFFUSE_SHARED_LIB)
  add_library (discus_all ${SOURCES_ALL})
endif (DIFFUSE_SHARED_LIB)

add_dependencies (discus_all kuplot_all lib_f90 )

set ( LIBS discus_all ${LIBS} discus_all ${LIBS_OPT} ) 

#add_executable (discus ${SOURCE_MAIN} ${SOURCES_LOCAL}) 

#target_link_libraries (discus ${LIBS}) 


add_custom_target(discushlp 
  DEPENDS      ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp
  COMMAND cat  ${DIFFUSE_SOURCE_DIR}/discus/prog/appl_dis.hlp 
               ${DIFFUSE_SOURCE_DIR}/lib_f90/lib_f90.hlp > 
               ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp)

#add_dependencies(discus discushlp)

#install (TARGETS discus DESTINATION bin)
#install (FILES ${DIFFUSE_BINARY_DIR}/discus/prog/discus.hlp DESTINATION share)
install (FILES ${DIFFUSE_SOURCE_DIR}/discus/prog/color.map DESTINATION share)
FILE(GLOB files "${DIFFUSE_SOURCE_DIR}/discus/prog/sysmac/*.mac")
install (FILES ${files} DESTINATION share/discus)
