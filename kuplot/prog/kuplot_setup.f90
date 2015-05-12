MODULE kuplot_setup_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE kuplot_setup (standalone)
!                                                                       
!     This routine makes inital setup of KUPLOT                         
!                                                                       
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
      IMPLICIT none 
!
      LOGICAL, INTENT(IN) :: standalone
!                                                                       
      include'date.inc' 
!                                                                       
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
!                                                                       
      prompt = pname 
      blank = ' ' 
      prompt_status = PROMPT_ON 
      prompt_status_old = PROMPT_ON 
!                                                                       
      CALL ini_ran (0) 
!                                                                       
!------ Write starting screen                                           
!                                                                       
      version = aktuell 
      IF(standalone) WRITE ( *, 1000) version, cdate 
!                                                                       
!     Call initialization routines                                      
!                                                                       
      CALL kuplot_initarrays 
      IF(standalone) CALL init_sysarrays 
      CALL appl_env 
      CALL kuplot_auto_def 
      IF(standalone) CALL cmdline_args 
      CALL no_error
!
      lsetup_done = .true.
!                                                                       
 1000 FORMAT (/,10x,59('*'),/,10x,'*',15x,                              &
     &        'K U P L O T   Version ',                                 &
     &        a10,10x,'*',/,10x,'*',57x,'*',/                           &
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
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE kuplot_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE kuplot_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE kuplot_calc_intr_spec 
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
   SUBROUTINE kuplot_branch(zeile, length)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
!
   END SUBROUTINE kuplot_branch
END INTERFACE
!
p_mache_kdo         => kuplot_mache_kdo
p_errlist_appl      => kuplot_errlist_appl
p_ersetz_para       => kuplot_ersetz_para
p_upd_para          => kuplot_upd_para
p_calc_intr_spec    => kuplot_calc_intr_spec
p_validate_var_spec => kuplot_validate_var_spec
p_branch            => kuplot_branch
!
END SUBROUTINE kuplot_set_sub
!
END MODULE kuplot_setup_mod
