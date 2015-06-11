MODULE diffev_setup_mod
!
CONTAINS
!
!*****7*****************************************************************
SUBROUTINE diffev_setup(standalone)
!                                                                       
!     This routine makes inital setup of DIFFEV                         
!                                                                       
USE diffev_allocate_appl
USE diffev_blk_appl
USE constraint
USE population
!
USE prompt_mod 
!
IMPLICIT none 
!
LOGICAL, INTENT(IN) :: standalone
!                                                                       
      include'date.inc' 
LOGICAL                        :: lend 
!                                                                       
lend              = .false. 
blank             = ' ' 
pname             = 'diffev' 
pname_cap         = 'DIFFEV' 
prompt            = pname 
prompt_status     = PROMPT_ON 
prompt_status_old = PROMPT_ON 
!                                                                       
CALL ini_ran (0) 
!                                                                       
!------ Write starting screen                                           
!                                                                       
version   = aktuell 
IF(standalone) WRITE ( *, 1000) version, cdate
!
!     Call initial default allocation
!
MAXPOP     = 0
MAXDIMX    = 0
MAX_CONSTR = 0
CALL diffev_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL diffev_initarrays 
IF(standalone) CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
CALL appl_env (lstandalone)
!                                                                       
!     try to read default file                                          
!                                                                       
CALL diffev_autodef 
!
!     Define Slave/stand alone status
!
IF(standalone) THEN
   pop_result_file_rd = .true.
   pop_trial_file_wrt = .true.
ELSE
   pop_result_file_rd = .false.
   pop_trial_file_wrt = .false.
ENDIF
!                                                                       
!     Check for command line parameters                                 
!                                                                       
IF(standalone) CALL cmdline_args
!
lsetup_done = .true.
!                                                                       
1000 FORMAT (/,                                                              &
     10x,59('*'),/,                                                          &
     10x,'*',15x,'D I F F E V   Version ',a6,14x,'*',/,                      &
     10x,'*',57(' '),'*',/                                                   &
     10x,'*         Created : ',a35,3x,'*',/,                                &
     10x,'*',57('-'),'*',/,                                                  &
     10x,'* (c) R.B. Neder  ','(reinhard.neder@fau.de)                 *',/, &
     10x,59('*'),/)                                        
END SUBROUTINE diffev_setup                          
!
SUBROUTINE diffev_set_sub
!
! Sets the specific DIFFEV interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE diffev_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE diffev_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_errlist_appl
   END SUBROUTINE diffev_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE diffev_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE diffev_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE diffev_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE diffev_validate_var_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_branch(zeile, length)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
!
   END SUBROUTINE diffev_branch
END INTERFACE
!
p_mache_kdo         => diffev_mache_kdo
p_errlist_appl      => diffev_errlist_appl
p_ersetz_para       => diffev_ersetz_para
p_upd_para          => diffev_upd_para
p_calc_intr_spec    => diffev_calc_intr_spec
p_validate_var_spec => diffev_validate_var_spec
p_branch            => diffev_branch
!
END SUBROUTINE diffev_set_sub
!
SUBROUTINE diffev_set_sub_cost
!
! Sets the specific DIFFEV interfaces for the cost calculation function
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE diffev_execute_cost( repeat,    &
                            prog  ,  prog_l , &
                            mac   ,  mac_l  , &
                            direc ,  direc_l, &
                            kid   ,  indiv  , &
                            rvalue, l_rvalue, &
                            output, output_l, &
                            generation, member, &
                            children, parameters, &
                            trial_v, NTRIAL, &
                            ierr )
!
   IMPLICIT NONE
   LOGICAL                , INTENT(IN) :: repeat
   INTEGER                , INTENT(IN) :: prog_l
   INTEGER                , INTENT(IN) :: mac_l
   INTEGER                , INTENT(IN) :: direc_l
   INTEGER                , INTENT(IN) :: output_l
   CHARACTER(LEN=prog_l  ), INTENT(IN) :: prog
   CHARACTER(LEN=mac_l   ), INTENT(IN) :: mac
   CHARACTER(LEN=direc_l ), INTENT(IN) :: direc
   INTEGER                , INTENT(IN) :: kid
   INTEGER                , INTENT(IN) :: indiv
   CHARACTER(LEN=output_l), INTENT(IN) :: output
   REAL                   , INTENT(OUT):: rvalue
   LOGICAL                , INTENT(OUT):: l_rvalue
   INTEGER                , INTENT(IN) :: generation
   INTEGER                , INTENT(IN) :: member
   INTEGER                , INTENT(IN) :: children
   INTEGER                , INTENT(IN) :: parameters
   INTEGER                , INTENT(IN) :: NTRIAL
   REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
   INTEGER                , INTENT(OUT):: ierr 
!
   END SUBROUTINE diffev_execute_cost
END INTERFACE
!
p_execute_cost      => diffev_execute_cost
!
END SUBROUTINE diffev_set_sub_cost
!
END MODULE diffev_setup_mod
