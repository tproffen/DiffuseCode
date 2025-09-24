MODULE diffev_setup_mod
!
CONTAINS
!
!*****7*****************************************************************
SUBROUTINE diffev_setup
!                                                                       
!     This routine makes inital setup of DIFFEV                         
!                                                                       
USE diffev_allocate_appl
USE diffev_blk_appl
USE constraint
USE population
USE gen_mpi_mod
!
USE appl_env_mod
USE cmdline_args_mod
USE prompt_mod 
USE lib_f90_default_mod
USE lib_init_mod
USE random_state_mod
!
IMPLICIT none 
!
INTEGER, PARAMETER  :: np = 1
INTEGER, DIMENSION(np) :: iwerte = 0
!                                                                       
      include'date.inc' 
CHARACTER(LEN=13)  :: is_debug
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
!CALL ini_ran (np, werte) 
IF(random_linit) CALL ini_ran_ix (np, iwerte, 0) 
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
!                                                                       
!------ Write starting screen                                           
!                                                                       
version   = aktuell 
!
!     try to read default file                                          
!                                                                       
CALL diffev_autodef 
!
!
lsetup_done = .true.
!                                                                       
END SUBROUTINE diffev_setup                          
!
!
END MODULE diffev_setup_mod
