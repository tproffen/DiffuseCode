MODULE refine_setup_mod
!
CONTAINS
!
!*****7*****************************************************************
SUBROUTINE refine_setup(standalone)
!                                                                       
!     This routine makes inital setup of REFINE                         
!                                                                       
!
USE refine_allocate_appl
USE refine_blk_appl
!
USE appl_env_mod
USE cmdline_args_mod
USE gen_mpi_mod
USE prompt_mod 
USE lib_f90_default_mod
USE random_state_mod
!
IMPLICIT none 
!
LOGICAL, INTENT(IN) :: standalone
!
INTEGER, PARAMETER  :: np = 1
!REAL, DIMENSION(np) :: werte = 0.0
INTEGER, DIMENSION(np) :: iwerte = 0
!                                                                       
      include'date.inc' 
CHARACTER(LEN=13)  :: is_debug
LOGICAL            :: lend 
!                                                                       
lend              = .false. 
blank             = ' ' 
pname             = 'refine' 
pname_cap         = 'REFINE' 
prompt            = pname 
prompt_status     = PROMPT_ON 
prompt_status_old = PROMPT_ON 
!                                                                       
!CALL ini_ran (np, werte) 
IF(random_linit) CALL ini_ran_ix (np, iwerte, 0) 
!
!     Call initial default allocation
!
IF(standalone) CALL lib_alloc_default
!MAXPOP     = 0
!MAXDIMX    = 0
!MAX_CONSTR = 0
CALL refine_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL refine_initarrays 
IF(standalone) CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
IF(standalone) CALL appl_env (lstandalone) !, gen_mpi_myid)
!                                                                       
!------ Write starting screen                                           
!                                                                       
version   = aktuell 
!
IF(standalone) THEN
   IF(cdebug=='ON') THEN
      is_debug = 'DEBUG VERSION'
   ELSE
      is_debug = '             '
   ENDIF
   WRITE ( *, 1000) version, is_debug, cdate
   CALL write_appl_env (lstandalone, gen_mpi_myid)
ENDIF
!                                                                       
!     try to read default file                                          
!                                                                       
CALL refine_autodef 
!
!     Define Slave/stand alone status
!
!IF(standalone) THEN
!   pop_result_file_rd = .true.
!   pop_trial_file_wrt = .true.
!ELSE
!   pop_result_file_rd = .false.
!   pop_trial_file_wrt = .false.
!ENDIF
!                                                                       
!     Check for command line parameters                                 
!                                                                       
IF(standalone) CALL cmdline_args(gen_mpi_myid)
!
lsetup_done = .true.
!                                                                       
1000 FORMAT (/,                                                              &
     10x,59('*'),/,                                                          &
     10x,'*',15x,'R E F I N E   Version ',a10,10x,'*',/,                     &
     10x,'*',22(' '),a13,22(' '),'*',/                                       &
     10x,'*         Created : ',a35,3x,'*',/,                                &
     10x,'*',57('-'),'*',/,                                                  &
     10x,'* (c) R.B. Neder  ','(reinhard.neder@fau.de)                 *',/, &
     10x,59('*'),/)                                        
END SUBROUTINE refine_setup                          
!
SUBROUTINE refine_set_sub
!
! Sets the specific REFINE interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE refine_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE refine_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_errlist_appl
   END SUBROUTINE refine_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   USE precision_mod
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE refine_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_upd_para (ctype, ww, maxw, wert, ianz, cstring)
!
   USE precision_mod
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL(KIND=PREC_DP), INTENT(IN   )    :: wert
   CHARACTER (LEN=* ), INTENT(IN   )    :: cstring
!
   END SUBROUTINE refine_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_get_var_type(line,length, var_is_type)
!
   CHARACTER(LEN=*)     , INTENT(IN)  :: line
   INTEGER              , INTENT(IN)  :: length
   INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
   END SUBROUTINE refine_get_var_type
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   USE precision_mod
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL(KIND=PREC_DP) , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE refine_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_calc_intr_log_spec(string, length)
!
   IMPLICIT NONE
   CHARACTER(LEN=*) , INTENT(INOUT) :: string
   INTEGER          , INTENT(INOUT) :: length
!
   END SUBROUTINE refine_calc_intr_log_spec
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE refine_validate_var_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_branch(zeile, length, lreset)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
   END SUBROUTINE refine_branch
END INTERFACE
!
!INTERFACE
!   SUBROUTINE refine_loop_mpi(prog_n, prog_l, mac_n, mac_l, out_n, out_l, repeat, nindiv)
!!
!CHARACTER (LEN=*), INTENT(IN) :: prog_n
!CHARACTER (LEN=*), INTENT(IN) :: mac_n
!CHARACTER (LEN=*), INTENT(IN) :: out_n
!INTEGER          , INTENT(IN) :: prog_l
!INTEGER          , INTENT(IN) :: mac_l
!INTEGER          , INTENT(IN) :: out_l
!LOGICAL          , INTENT(IN) :: repeat
!INTEGER          , INTENT(IN) :: nindiv
!!
!   END SUBROUTINE refine_loop_mpi
!END INTERFACE
!
p_mache_kdo         => refine_mache_kdo
p_errlist_appl      => refine_errlist_appl
p_ersetz_para       => refine_ersetz_para
p_upd_para          => refine_upd_para
p_calc_intr_spec    => refine_calc_intr_spec
p_calc_intr_log_spec=> refine_calc_intr_log_spec
p_validate_var_spec => refine_validate_var_spec
p_branch            => refine_branch
!p_loop_mpi          => refine_loop_mpi
p_get_var_type      => refine_get_var_type
!
END SUBROUTINE refine_set_sub
!
!SUBROUTINE refine_set_sub_cost
!!
!! Sets the specific REFINE interfaces for the cost calculation function
!!
!USE set_sub_generic_mod
!!
!INTERFACE
!   SUBROUTINE refine_execute_cost( repeat,    &
!                            prog_len,         &
!                            prog  ,  prog_l , &
!                            mac_len,          &
!                            mac   ,  mac_l  , &
!                            direc_len,        &
!                            direc ,  direc_l, &
!!                            kid   ,  indiv  , &
!!                            n_rvalue_i, n_rvalue_o, NRVAL,  &
!                            rvalue, l_rvalue, &
!                            output_len,       &
!                            output, output_l, &
!                            generation, member, &
!                            children, parameters, &
!                                    nindiv  , &
!                            trial_n,         &
!                            trial_v, NTRIAL, &
!                            l_get_random_state,     &
!                            rd_nseeds,rd_seeds,     &
!                            l_first_job,            &
!                            ierr )
!!
!   IMPLICIT NONE
!   LOGICAL                , INTENT(IN) :: repeat
!   INTEGER                , INTENT(IN) :: prog_len
!   INTEGER                , INTENT(IN) :: prog_l
!   INTEGER                , INTENT(IN) :: mac_len
!   INTEGER                , INTENT(IN) :: mac_l
!   INTEGER                , INTENT(IN) :: direc_len
!   INTEGER                , INTENT(IN) :: direc_l
!   INTEGER                , INTENT(IN) :: output_len
!   INTEGER                , INTENT(IN) :: output_l
!   CHARACTER(LEN=prog_len  ), INTENT(IN) :: prog
!   CHARACTER(LEN=mac_len   ), INTENT(IN) :: mac
!   CHARACTER(LEN=direc_len ), INTENT(IN) :: direc
!   INTEGER                , INTENT(IN) :: kid
!   INTEGER                , INTENT(IN) :: indiv
!   CHARACTER(LEN=output_len), INTENT(IN) :: output
!   INTEGER                , INTENT(IN ):: n_rvalue_i
!   INTEGER                , INTENT(OUT):: n_rvalue_o
!   INTEGER                , INTENT(IN ):: NRVAL
!   REAL, DIMENSION(0:NRVAL)    , INTENT(OUT):: rvalue
!   LOGICAL                , INTENT(OUT):: l_rvalue
!   INTEGER                , INTENT(IN) :: generation
!   INTEGER                , INTENT(IN) :: member
!   INTEGER                , INTENT(IN) :: children
!   INTEGER                , INTENT(IN) :: parameters
!   INTEGER                , INTENT(IN) :: nindiv
!   INTEGER                , INTENT(IN) :: NTRIAL
!   CHARACTER(LEN=16),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_n
!   REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
!   LOGICAL                , INTENT(IN)  :: l_get_random_state
!   INTEGER                , INTENT(OUT) :: rd_nseeds
!   INTEGER, DIMENSION(64) , INTENT(OUT) :: rd_seeds
!   LOGICAL                , INTENT(IN ) :: l_first_job
!   INTEGER                , INTENT(OUT):: ierr 
!!
!   END SUBROUTINE refine_execute_cost
!END INTERFACE
!!
!p_execute_cost      => refine_execute_cost
!!
!END SUBROUTINE refine_set_sub_cost
!
END MODULE refine_setup_mod
