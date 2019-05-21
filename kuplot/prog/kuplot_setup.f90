MODULE kuplot_setup_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE kuplot_setup (standalone)
!                                                                       
!     This routine makes inital setup of KUPLOT                         
!                                                                       
      USE prompt_mod 
      USE lib_f90_default_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE cmdline_args_mod
      USE appl_env_mod
      USE random_state_mod
      IMPLICIT none 
!
      LOGICAL, INTENT(IN) :: standalone
!
      INTEGER, PARAMETER  :: np = 1
!     REAL, DIMENSION(np) :: werte = 0.0
      INTEGER, DIMENSION(np) :: iwerte = 0
!                                                                       
      include'date.inc' 
      CHARACTER(LEN=13)  :: is_debug
!                                                                       
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
!                                                                       
      prompt = pname 
      blank = ' ' 
      prompt_status = PROMPT_ON 
      prompt_status_old = PROMPT_ON 
!                                                                       
!     CALL ini_ran (np, werte) 
      CALL ini_ran_ix (np, iwerte) 
!                                                                       
!     Call initialization routines                                      
!                                                                       
      IF(standalone) CALL lib_alloc_default
      CALL kuplot_initarrays 
      IF(standalone) CALL init_sysarrays 
      CALL appl_env (lstandalone,0)
!                                                                       
!------ Write starting screen                                           
!                                                                       
      version = aktuell 
      IF(standalone) THEN
         IF(cdebug=='ON') THEN
            is_debug = 'DEBUG VERSION'
         ELSE
            is_debug = '             '
         ENDIF
         WRITE ( *, 1000) version, is_debug, cdate 
         CALL write_appl_env (lstandalone,0)
      ENDIF
      CALL kuplot_auto_def 
      IF(standalone) CALL cmdline_args (0)
      CALL no_error
!
      lsetup_done = .true.
!                                                                       
 1000 FORMAT (/,10x,59('*'),/,10x,'*',15x,                              &
     &        'K U P L O T   Version ',                                 &
     &        a10,10x,'*',/,10x,'*',22(' '),a13,22(' '),'*',/           &
     &        10x,'*         Created : ',a35,3x,'*',/,10x,'*',          &
     &        57('-'),'*',/,10x,'* (c) Th. Proffen ',                   &
     &        '(tproffen@ornl.gov)                     *',/,            &
     &        10x,'*     R.B. Neder  ',                                 &
     &        '(reinhard.neder@fau.de)                 *',/,            &
     &        10x,59('*'),/,                                            &
     &        10x,'* GSAS code: (c) Allen C. Larson and',               &
     &        ' Robert B. Von Dreele *',/,                              &
     &        10x,59('*'),/)                                            
      END SUBROUTINE kuplot_setup                          
!
SUBROUTINE kuplot_set_sub
!
! Sets the specific DIFFEV interfaces four routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
USE prompt_mod
!
INTERFACE
   SUBROUTINE kuplot_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE kuplot_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_errlist_appl
   END SUBROUTINE kuplot_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
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
   END SUBROUTINE kuplot_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_upd_para (ctype, ww, maxw, wert, ianz)
!
   USE precision_mod
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL(KIND=PREC_DP), INTENT(IN   )    :: wert
!
   END SUBROUTINE kuplot_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_get_var_type(line,length, var_is_type)
!
   CHARACTER(LEN=*)     , INTENT(IN)  :: line
   INTEGER              , INTENT(IN)  :: length
   INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
   END SUBROUTINE kuplot_get_var_type
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   USe precision_mod
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL(KIND=PREC_DP) , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE kuplot_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_calc_intr_log_spec(string, length)
!
   IMPLICIT NONE
   CHARACTER(LEN=*) , INTENT(INOUT) :: string
   INTEGER          , INTENT(INOUT) :: length
!
   END SUBROUTINE kuplot_calc_intr_log_spec
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE kuplot_validate_var_spec 
END INTERFACE

INTERFACE
   SUBROUTINE kuplot_branch(zeile, length, lreset)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
   END SUBROUTINE kuplot_branch
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_top(zeile)
!
   CHARACTER (LEN=*), INTENT(IN) :: zeile
   END SUBROUTINE kuplot_top
END INTERFACE
!
p_mache_kdo         => kuplot_mache_kdo
p_errlist_appl      => kuplot_errlist_appl
p_ersetz_para       => kuplot_ersetz_para
p_upd_para          => kuplot_upd_para
p_calc_intr_spec    => kuplot_calc_intr_spec
p_calc_intr_log_spec=> kuplot_calc_intr_log_spec
p_validate_var_spec => kuplot_validate_var_spec
p_branch            => kuplot_branch
p_loop_mpi          => dummy_loop_mpi
p_get_var_type      => kuplot_get_var_type
IF(lstandalone) THEN
p_top               => kuplot_top
ENDIF
!
END SUBROUTINE kuplot_set_sub
!
END MODULE kuplot_setup_mod
