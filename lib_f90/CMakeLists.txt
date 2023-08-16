# Create the lib_90 library for the command language

set (FSOURCES allocate_generic.f90 lib_f90_alloc.f90 lib_f90_default.f90
              appl_unix.f90 arrays_mod.f90 blanks_mod.f90 berechne_mod.f90 
              ber_params_mod.f90 blockdat.f90 build_name_mod.f90
              calc_expr_mod.f90  calc_intr_mod.f90 cmdline_args_mod.f90
              constants_mod.f90
              get_params_mod.f90 count_col_mod.f90 
              data_struc.f90 data_struc_type.f90 data_types.f90
              define_variable_mod.f90 dkdo.f90 do_echo_mod.f90 do_eval_mod.f90
              do_execute_mod.f90 do_read_number_mod.f90 do_replace_expr.f90
              do_show_mod.f90 do_set_mod.f90
              do_string_alloc_mod.f90 do_variable_mod.f90 do_wait_mod.f90
              element_data_mod.f90
              errlist.f90 ersetz_mod.f90 ersetzl_mod.f90 exit_mod.f90 exit_para.f90
              fast_fourier_mod.f90
              fput_mod.f90 help.f90 jsu_readline.f90 kdo_all_mod.f90 learn.f90 
              gauss_lorentz_pseudo.f90
              gamma_mod.f90 gaussj_mod.f90 gen_mpi_mod.f90 global_data.f90
              hdf5_def.f90 
              length_string_mod.f90 lib_f90_config.f90 lib_f90_fftw3.f90 
              lib_conv.f90 lib_functions.f90
              lib_get_var_type.f90  lib_ik.f90 lib_lanczos.f90 lib_load.f90
              lib_math.f90 lib_metric.f90 lib_f90_profile.f90 lib_trans.f90 
              lib_use_coor.f90 lib_weights.f90
              lib_write_debug.f90 lib_upd_mod.f90 
              macro.f90  matrix_mod.f90 map_1dtofield.f90 
              mrc_data.f90
              operating_mod.f90 reserved_mod.f90
              random.f90 random_state_mod.f90 str_comp_mod.f90 
              string_convert_mod.f90
              sup.f90 set_sub_generic_mod.f90 variable_test.f90 
              charact_mod.f90 debug_mod.f90 doact_mod.f90 doexec_mod.f90 
              doloop_mod.f90 envir_mod.f90 errlist_mod.f90 learn_mod.f90
              terminal_mod.f90 mpi_slave_mod.f90 dummy_loop_mpi.f90
              class_macro.f90 macro_internal.f90 macro_mod.f90 
              math_sup.f90 trig_degree_mod.f90
              param_mod.f90 profile_tof.f90 prompt_mod.f90 random_mod.f90 
              search_string_mod.f90 sorting.f90 take_param_mod.f90
              sine_table.f90 spline_mod.f90 string_extract_mod.f90
              times_mod.f90 variable_mod.f90 wink_mod.f90 precision_mod.f90 
              memory.f90 support_mod.f90 
              refine_params.f90)
#             sockets_mod.f90

#set (CSOURCES FCreadline.c csocket.c csup.c)
set (CSOURCES FCreadline.c csup.c)

# Will change the sys_*.f90 file depending on which
# compiler is used. Default to gfortran.
if (PGF)
  set( FSOURCES ${FSOURCES} sys_pgi.f90 )
elseif (IFORT)
  set( FSOURCES ${FSOURCES} sys_intel.f90 )
else (PGF)
  set( FSOURCES ${FSOURCES} sys_gfortran.f90 )
endif (PGF)

# Operating system dependent fortran code
if ( CYGWIN ) 
  set ( FSOURCES ${FSOURCES} op_windows.f90 )
  message ( status " OP SYSTEM is Cygwin" )
elseif ( APPLE )
  set ( FSOURCES ${FSOURCES} op_macos.f90 )
  message ( status " OP SYSTEM is MACOS" )
elseiF( UNIX ) 
  set ( FSOURCES ${FSOURCES} op_linux.f90 )
  message ( status " OP SYSTEM is UNIX" )
else ( CYGWIN )
  set ( FSOURCES ${FSOURCES} op_linux.f90 )
  message ( status " OP SYSTEM is unknown" )
endif ( CYGWIN ) 

# OpenMP ?

if (CUDA_FOUND)
  message (status " building LIBF90 with cuda")
  set ( FSOURCES ${FSOURCES} parallel_CUDA.f90)
elseif( DISCUS_OMP)
  message (status " building LIBF90 with openmp")
  set ( FSOURCES ${FSOURCES} parallel_mod.f90)
  # Add OpenMP flag depending on compiler
  if (PGF)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mp")
  elseif (IFORT)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
  else (PGF) # Else assume gfortran
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  endif (PGF)
else (CUDA_FOUND)
  message (status " building LIBF90 serial ")
  set (FSOURCES ${FSOURCES} parallel_none.f90)
endif (CUDA_FOUND)
#
if (DIFFUSE_HDF5)
  if (DEFINED ENV{HDF5_INC_DIR})
    include_directories($ENV{HDF5_INC_DIR})
  endif ()
  set (FSOURCES ${FSOURCES} hdf_write.f90 hdf5_read.f90 hdf5_params.f90)
  message (LIB_F90 " HDF5 is ON  ${DIFFUSE_HDF5} ")
else (DIFFUSE_HDF5)
  set (FSOURCES ${FSOURCES} hdf_write_no.f90 hdf5_params_no.f90)
  message (LIB_F90 " HDF5 is OFF ${DIFFUSE_HDF5} ")
endif ()

include_directories(${FFTW_INCLUDE_DIRS})


# This macro modified from http://www.cmake.org/pipermail/cmake/2010-November/041034.html
# The date is returned in ISO8601 format
MACRO (TODAY RESULT)
    IF(UNIX OR CYGWIN OR LINUX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/(....).*" "\\3-\\2-\\1" ${RESULT} ${${RESULT}})
    ELSEIF (WIN32)
        EXECUTE_PROCESS(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/(....).*" "\\3-\\2-\\1" ${RESULT} ${${RESULT}})
    ELSE (UNIX OR CYGWIN OR LINUX)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} "0000-00-00")
    ENDIF (UNIX OR CYGWIN OR LINUX)
ENDMACRO (TODAY)

TODAY(BUILD_DATE)
message (status " BUILD_DATE=${BUILD_DATE} ")
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/date.inc.template
                 ${CMAKE_CURRENT_BINARY_DIR}/date.inc )
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/version.inc.template
                 ${CMAKE_CURRENT_BINARY_DIR}/version.inc )

IF (WIN32 OR CYGWIN )
  set (CSOURCES ${CSOURCES} win32-glob.c)
ENDIF (WIN32 OR CYGWIN)

add_library (lib_f90  ${FSOURCES})
add_library (lib_f90c ${CSOURCES})
#
