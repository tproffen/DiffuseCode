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
USE discus_reset_all_mod
!
USE appl_env_mod
USE cmdline_args_mod
USE errlist_mod
USE prompt_mod 
USE lib_errlist_func
USE lib_f90_default_mod
USE lib_init_mod
USE random_state_mod
!
IMPLICIT none 
!
LOGICAL, INTENT(IN) :: standalone
!
INTEGER, PARAMETER  :: np = 1
INTEGER, DIMENSION(np) :: iwerte = 0
!                                                                       
      include'date.inc' 
      CHARACTER(LEN=13)  :: is_debug
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
IF(random_linit) CALL ini_ran_ix (np,iwerte, 0) 
!
!     Call initial default allocation
IF(standalone) CALL lib_alloc_default
!
      CALL discus_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL discus_initarrays 
CALL discus_reset_all
IF(standalone) CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
IF(standalone) CALL appl_env (lstandalone) !, 0)
!                                                                       
!     try to read default file                                          
!                                                                       
CALL discus_autodef 
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
              10x,'*',22(' '),a13,22(' '),'*',/                     &
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
!
END MODULE discus_setup_mod
