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
LOGICAL, INTENT(IN) :: standalone
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
IF(standalone) CALL lib_alloc_default
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
IF(standalone) CALL appl_env (lstandalone) !, gen_mpi_myid)
!                                                                       
!------ Write starting screen                                           
!                                                                       
version   = aktuell 
!
IF(standalone) THEN
   IF(cdebug=='ON') THEN
      is_debug = 'DEBUG VERSION'
   ELSE
      is_debug = '             '
   ENDIF
   WRITE ( *, 1000) version, is_debug, cdate
   CALL write_appl_env (lstandalone, gen_mpi_myid)
ENDIF
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
IF(standalone) CALL cmdline_args(gen_mpi_myid)
!
lsetup_done = .true.
!                                                                       
1000 FORMAT (/,                                                              &
     10x,59('*'),/,                                                          &
     10x,'*',15x,'D I F F E V   Version ',a10,10x,'*',/,                     &
     10x,'*',22(' '),a13,22(' '),'*',/                                       &
     10x,'*         Created : ',a35,3x,'*',/,                                &
     10x,'*',57('-'),'*',/,                                                  &
     10x,'* (c) R.B. Neder  ','(reinhard.neder@fau.de)                 *',/, &
     10x,59('*'),/)                                        
END SUBROUTINE diffev_setup                          
!
!
END MODULE diffev_setup_mod
