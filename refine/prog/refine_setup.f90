MODULE refine_setup_mod
!
CONTAINS
!
!*****7*****************************************************************
SUBROUTINE refine_setup !(standalone)
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
USE lib_init_mod
USE random_state_mod
!
IMPLICIT none 
!
!LOGICAL, INTENT(IN) :: standalone
!
INTEGER, PARAMETER  :: np = 1
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
IF(random_linit) CALL ini_ran_ix (np, iwerte, 0) 
!
!MAXPOP     = 0
!MAXDIMX    = 0
!MAX_CONSTR = 0
CALL refine_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL refine_initarrays 
!                                                                       
version   = aktuell 
!
!     try to read default file                                          
!                                                                       
CALL refine_autodef 
!
lsetup_done = .true.
!                                                                       
END SUBROUTINE refine_setup                          
!
!*******************************************************************************
!
END MODULE refine_setup_mod
