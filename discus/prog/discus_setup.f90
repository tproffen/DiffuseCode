MODULE discus_setup_mod
!
CONTAINS
!
!
SUBROUTINE discus_setup (standalone)
!                                                                       
!     This routine makes inital setup of DISCUS                         
!                                                                       
      USE discus_allocate_appl_mod
      USE discus_init_mod
!
      USE appl_env_mod
      USE cmdline_args_mod
      USE errlist_mod
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
!
      pname      = 'discus'
      pname_cap  = 'DISCUS'
!                                                                       
      blank   = ' '
      prompt  = pname
      prompt_status = PROMPT_ON
      prompt_status_old = PROMPT_ON
!                                                                       
!CALL ini_ran (np,werte) 
CALL ini_ran_ix (np,iwerte) 
!
!     Call initial default allocation
IF(standalone) CALL lib_alloc_default
!
      CALL discus_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL discus_initarrays 
IF(standalone) CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
CALL appl_env (lstandalone, 0)
!                                                                       
!     try to read default file                                          
!                                                                       
CALL discus_autodef 
!                                                                       
!------ Write starting screen                                           
!                                                                       
version = aktuell 
IF(standalone) THEN
   WRITE ( *, 1000) version, cdate 
   CALL write_appl_env (lstandalone, 0)
ENDIF
!                                                                       
!     Check for command line parameters                                 
!                                                                       
IF(standalone) CALL cmdline_args(0) 
!
CALL no_error
lsetup_done = .true.
!
 1000 FORMAT (/,10x,59('*'),/,                                      &
              10x,'*',15x,'D I S C U S   Version ',a10,10x,'*',/,   &
              10x,'*',57(' '),'*',/                                 &
     &        10x,'*         Created : ',a35,3x,'*',/,              &
              10x,'*',57('-'),'*',/,                                &
     &        10x,'* (c) R.B. Neder  ',                             &
     &        '(reinhard.neder@fau.de)                 *',/,        &
     &        10x,'*     Th. Proffen ',                             &
     &        '(tproffen@ornl.gov)                     *',/,        &
     &        10x,59('*'),/,                                        &
              10x,'*',57(' '),'*',/,                                &
     &        10x,'* For information on current changes',           &
     &            ' type: help News',6x,'*',/,                      &
     &        10x,'*',57(' '),'*',/,10x,59('*'),/                   &
     &                     )                                            
END SUBROUTINE discus_setup                          
!
SUBROUTINE discus_set_sub
!
! Sets the specific DIFFEV interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE discus_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE discus_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_errlist_appl
   END SUBROUTINE discus_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE discus_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE discus_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_get_var_type(line,length, var_is_type)
!
   CHARACTER(LEN=*)     , INTENT(IN)  :: line
   INTEGER              , INTENT(IN)  :: length
   INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
   END SUBROUTINE discus_get_var_type
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE discus_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_calc_intr_log_spec(string, length)
!
   IMPLICIT NONE
   CHARACTER(LEN=*) , INTENT(INOUT) :: string
   INTEGER          , INTENT(INOUT) :: length
!
   END SUBROUTINE discus_calc_intr_log_spec
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE discus_validate_var_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_branch(zeile, length)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
!
   END SUBROUTINE discus_branch
END INTERFACE
!
p_mache_kdo         => discus_mache_kdo
p_errlist_appl      => discus_errlist_appl
p_ersetz_para       => discus_ersetz_para
p_upd_para          => discus_upd_para
p_calc_intr_spec    => discus_calc_intr_spec
p_calc_intr_log_spec=> discus_calc_intr_log_spec
p_validate_var_spec => discus_validate_var_spec
p_branch            => discus_branch
p_loop_mpi          => dummy_loop_mpi
p_get_var_type      => discus_get_var_type
!
END SUBROUTINE discus_set_sub
!
END MODULE discus_setup_mod
