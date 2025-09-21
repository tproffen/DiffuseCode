MODULE discus_setup_mod
!
CONTAINS
!
!
SUBROUTINE discus_setup
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
IF(random_linit) CALL ini_ran_ix (np,iwerte, 0) 
!
!     Call initial default allocation
!
CALL discus_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL discus_initarrays 
CALL discus_reset_all
!                                                                       
!     try to read default file                                          
!                                                                       
CALL discus_autodef 
!                                                                       
version = aktuell 
!                                                                       
!
CALL no_error
lsetup_done = .true.
!
END SUBROUTINE discus_setup                          
!
!
END MODULE discus_setup_mod
