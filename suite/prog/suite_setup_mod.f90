MODULE suite_setup_mod
!
CONTAINS
!
!*****7*****************************************************************
SUBROUTINE setup_suite
!                                                                       
!     This routine makes inital setup of DISCUS_SUITE                         
!                                                                       
!USE allocate_appl
!USE blk_appl
!
USE run_mpi_mod
USE prompt_mod
USE lib_f90_default_mod
!
IMPLICIT none
!                                                                       
      include'date.inc'
LOGICAL                        :: lend
INTEGER, PARAMETER  :: np = 1
REAL, DIMENSION(np) :: werte = 0.0
!                                                                       
lend              = .false.
blank             = ' '
pname             = 'suite'
pname_cap         = 'SUITE'
prompt            = pname
prompt_status     = PROMPT_ON
prompt_status_old = PROMPT_ON
!                                                                       
CALL ini_ran (np, werte)
!
!     Call initial default allocation
!
CALL lib_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
!CALL initarrays
CALL init_sysarrays
!                                                                       
!     get envirmonment information                                      
!                                                                       
CALL appl_env (.TRUE., run_mpi_myid)
!
IF(run_mpi_myid==0) THEN
!                                                                       
!------ Write starting screen                                           
!                                                                       
version   = aktuell
WRITE ( *, 1000) version, cdate
CALL write_appl_env (.TRUE., run_mpi_myid)
ENDIF
!                                                                       
!     try to read default file                                          
!                                                                       
!CALL autodef
!                                                                       
!     Check for command line parameters                                 
!                                                                       
CALL cmdline_args(run_mpi_myid)
!                                                                       
lsetup_done = .true.
!
!write(fname,2000) run_mpi_myid
!2000 format('MEMORY.',I4.4) 
!open(101,file=fname, status='unknown')
!write(101,'(a)') '#INFO            vmpeak     vmsize     vmhwm    vmrss    vmpte'
!write(101,'(a)') '#'
!                                                                       
1000 FORMAT (/,                                                              &
     10x,59('*'),/,                                                          &
     10x,'*', 9x,'D I S C U S - S U I T E  Version ',a6, 9x,'*',/,           &
     10x,'*',57(' '),'*',/                                                   &
     10x,'*         Created : ',a35,3x,'*',/,                                &
     10x,'*',57('-'),'*',/,                                                  &
     10x,'* (c) R.B. Neder  ','(reinhard.neder@fau.de)                 *',/, &
     10x,59('*'),/)
END SUBROUTINE setup_suite
!
SUBROUTINE suite_set_sub
!
! Sets the specific SUITE  interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE suite_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE suite_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_errlist_appl
   END SUBROUTINE suite_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE suite_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE suite_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE suite_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE suite_validate_var_spec 
END INTERFACE

!
p_mache_kdo         => suite_mache_kdo
p_errlist_appl      => suite_errlist_appl
p_ersetz_para       => suite_ersetz_para
p_upd_para          => suite_upd_para
p_calc_intr_spec    => suite_calc_intr_spec
p_validate_var_spec => suite_validate_var_spec
!
END SUBROUTINE suite_set_sub
!
SUBROUTINE  suite_set_sub_cost
!
! Sets the specific SUITE  interfaces for the cost calculation function
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE suite_execute_cost( repeat,    &
                            prog_len,         &
                            prog  ,  prog_l , &
                            mac_len,          &
                            mac   ,  mac_l  , &
                            direc_len,        &
                            direc ,  direc_l, &
                            kid   ,  indiv  , &
                            rvalue, l_rvalue, &
                            output_len,       &
                            output, output_l, &
                            generation, member, &
                            children, parameters, &
                                    nindiv  , &
                            trial_n,              &
                            trial_v,  NTRIAL,     &
                            l_get_random_state,     &
                            rd_nseeds,rd_seeds,     &
                            ierr )
!
   IMPLICIT NONE
   LOGICAL                , INTENT(IN) :: repeat
   INTEGER                , INTENT(IN) :: prog_len
   INTEGER                , INTENT(IN) :: prog_l
   INTEGER                , INTENT(IN) :: mac_len
   INTEGER                , INTENT(IN) :: mac_l
   INTEGER                , INTENT(IN) :: direc_len
   INTEGER                , INTENT(IN) :: direc_l
   INTEGER                , INTENT(IN) :: output_len
   INTEGER                , INTENT(IN) :: output_l
   CHARACTER(LEN=prog_len), INTENT(IN) :: prog
   CHARACTER(LEN=mac_len ), INTENT(IN) :: mac
   CHARACTER(LEN=direc_len ), INTENT(IN) :: direc
   INTEGER                , INTENT(IN) :: kid
   INTEGER                , INTENT(IN) :: indiv
   CHARACTER(LEN=output_len), INTENT(IN) :: output
   REAL                   , INTENT(OUT):: rvalue
   LOGICAL                , INTENT(OUT):: l_rvalue
   INTEGER                , INTENT(IN) :: generation
   INTEGER                , INTENT(IN) :: member
   INTEGER                , INTENT(IN) :: children
   INTEGER                , INTENT(IN) :: parameters
   INTEGER                , INTENT(IN) :: nindiv
   INTEGER                , INTENT(IN) :: NTRIAL
   CHARACTER(LEN=16),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_n
   REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
   LOGICAL                , INTENT(IN)  :: l_get_random_state
   INTEGER                , INTENT(OUT) :: rd_nseeds
   INTEGER, DIMENSION(64) , INTENT(OUT) :: rd_seeds
   INTEGER                , INTENT(OUT):: ierr 
!
   END SUBROUTINE suite_execute_cost
END INTERFACE
!
p_execute_cost      =>  suite_execute_cost
!
END SUBROUTINE suite_set_sub_cost
!
SUBROUTINE suite_set_sub_branch
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE suite_branch(zeile, length)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
!
   END SUBROUTINE suite_branch
END INTERFACE
!
p_branch            =>  suite_branch 
!
END SUBROUTINE suite_set_sub_branch
!
END MODULE suite_setup_mod
