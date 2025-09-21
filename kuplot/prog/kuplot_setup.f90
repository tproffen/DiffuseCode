MODULE kuplot_setup_mod
!
CONTAINS
!*****7*****************************************************************
!
SUBROUTINE kuplot_setup
!                                                                       
!     This routine makes inital setup of KUPLOT                         
!                                                                       
USE prompt_mod 
USE lib_f90_default_mod
USE lib_init_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_blk_mod
USE cmdline_args_mod
USE appl_env_mod
USE lib_errlist_func
USE random_state_mod
!
IMPLICIT none 
!
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
      IF(random_linit) CALL ini_ran_ix (np, iwerte, 0) 
!                                                                       
!     Call initialization routines                                      
!                                                                       
      CALL kuplot_initarrays 
!                                                                       
!------ Write starting screen                                           
!                                                                       
      version = aktuell 
      CALL kuplot_auto_def 
      CALL no_error
!
      lsetup_done = .true.
!                                                                       
!1000 FORMAT (/,10x,59('*'),/,10x,'*',15x,                              &
!    &        'K U P L O T   Version ',                                 &
!    &        a10,10x,'*',/,10x,'*',22(' '),a13,22(' '),'*',/           &
!    &        10x,'*         Created : ',a35,3x,'*',/,10x,'*',          &
!    &        57('-'),'*',/,10x,'* (c) Th. Proffen ',                   &
!    &        '(tproffen@ornl.gov)                     *',/,            &
!    &        10x,'*     R.B. Neder  ',                                 &
!    &        '(reinhard.neder@fau.de)                 *',/,            &
!    &        10x,59('*'),/,                                            &
!    &        10x,'* GSAS code: (c) Allen C. Larson and',               &
!    &        ' Robert B. Von Dreele *',/,                              &
!    &        10x,59('*'),/)                                            
END SUBROUTINE kuplot_setup                          
!
!*****7*****************************************************************
!
END MODULE kuplot_setup_mod
