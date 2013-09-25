MODULE diffev_setup_mod
!
CONTAINS
!
SUBROUTINE diffev_setup
!                                                                       
USE diffev_mpi_mod
USE run_mpi_mod
USE doact_mod
USE errlist_mod 
USE learn_mod 
USE macro_mod 
USE prompt_mod 
!                                                                       
IMPLICIT none 
!
!*****7*****************************************************************
!                                                                       
!                                                                       
CHARACTER (LEN=1024)           :: line, zeile 
CHARACTER (LEN=4)              :: befehl 
LOGICAL                        :: lend 
INTEGER                        :: laenge, lp, lbef 
!
INTEGER, PARAMETER             :: master = 0 ! MPI ID of MASTER process
!                                                                       
run_mpi_myid      = 0
lend      = .false. 
blank     = ' ' 
pname     = 'diffev' 
pname_cap = 'DIFFEV' 
prompt            = pname 
prompt_status     = PROMPT_ON 
prompt_status_old = PROMPT_ON 
!
!------ Setting up variables and print start screen                     
!                                                                       
CALL setup
CALL run_mpi_init
lsetup_done = .true.
!                                                                       
END SUBROUTINE diffev_setup
!*****7*****************************************************************
SUBROUTINE setup 
!                                                                       
!     This routine makes inital setup of DIFFEV                         
!                                                                       
USE allocate_appl
USE blk_appl
USE constraint
USE population
!
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
      include'date.inc' 
!                                                                       
CALL ini_ran (0) 
!                                                                       
!------ Write starting screen                                           
!                                                                       
version   = aktuell 
WRITE ( *, 1000) version, cdate
!
!     Call initial default allocation
!
MAXPOP     = 0
MAXDIMX    = 0
MAX_CONSTR = 0
CALL alloc_appl
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL initarrays 
CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
CALL appl_env 
!                                                                       
!     try to read default file                                          
!                                                                       
CALL autodef 
!                                                                       
!     Check for command line parameters                                 
!                                                                       
CALL cmdline_args
!                                                                       
1000 FORMAT (/,                                                              &
     10x,59('*'),/,                                                          &
     10x,'*',15x,'D I F F E V   Version ',a6,14x,'*',/,                      &
     10x,'*',57(' '),'*',/                                                   &
     10x,'*         Created : ',a35,3x,'*',/,                                &
     10x,'*',57('-'),'*',/,                                                  &
     10x,'* (c) R.B. Neder  ','(reinhard.neder@krist.uni-erlangen.de)  *',/, &
     10x,59('*'),/)                                        
END SUBROUTINE setup                          
END MODULE diffev_setup_mod
