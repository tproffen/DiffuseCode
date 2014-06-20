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
USE diffev_mpi_mod
USE population
USE run_mpi_mod
!
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
      include'date.inc' 
LOGICAL                        :: lend 
INTEGER, PARAMETER             :: master = 0 ! MPI ID of MASTER process
!                                                                       
run_mpi_myid      = 0
lend              = .false. 
blank             = ' ' 
pname             = 'diffev' 
pname_cap         = 'DIFFEV' 
prompt            = pname 
prompt_status     = PROMPT_ON 
prompt_status_old = PROMPT_ON 
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
CALL diffev_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL diffev_initarrays 
CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
CALL appl_env 
!                                                                       
!     try to read default file                                          
!                                                                       
CALL diffev_autodef 
!                                                                       
!     Check for command line parameters                                 
!                                                                       
CALL cmdline_args
!
CALL run_mpi_init
lsetup_done = .true.
!                                                                       
1000 FORMAT (/,                                                              &
     10x,59('*'),/,                                                          &
     10x,'*',15x,'D I F F E V   Version ',a6,14x,'*',/,                      &
     10x,'*',57(' '),'*',/                                                   &
     10x,'*         Created : ',a35,3x,'*',/,                                &
     10x,'*',57('-'),'*',/,                                                  &
     10x,'* (c) R.B. Neder  ','(reinhard.neder@krist.uni-erlangen.de)  *',/, &
     10x,59('*'),/)                                        
END SUBROUTINE diffev_setup                          
!
SUBROUTINE diffev_set_sub
!
! Sets the specific DIFFEV interfaces four routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE diffev_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE diffev_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_errlist_appl
   END SUBROUTINE diffev_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE diffev_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE diffev_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE diffev_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE diffev_validate_var_spec 
END INTERFACE

!
p_mache_kdo         => diffev_mache_kdo
p_errlist_appl      => diffev_errlist_appl
p_ersetz_para       => diffev_ersetz_para
p_upd_para          => diffev_upd_para
p_calc_intr_spec    => diffev_calc_intr_spec
p_validate_var_spec => diffev_validate_var_spec
!
END SUBROUTINE diffev_set_sub
!
END MODULE diffev_setup_mod
