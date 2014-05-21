MODULE setup_mod
!
CONTAINS
!
!
SUBROUTINE setup 
!                                                                       
!     This routine makes inital setup of DISCUS                         
!                                                                       
      USE discus_allocate_appl_mod
      USE discus_init_mod
!
      USE errlist_mod
      USE prompt_mod 
!
IMPLICIT none 
!                                                                       
      include'date.inc' 
!
      pname      = 'discus'
      pname_cap  = 'DISCUS'
!                                                                       
      blank   = ' '
      prompt  = pname
      prompt_status = PROMPT_ON
      prompt_status_old = PROMPT_ON
!                                                                       
CALL ini_ran (0) 
!                                                                       
!------ Write starting screen                                           
!                                                                       
version = aktuell 
WRITE ( *, 1000) version, cdate 
!
!     Call initial default allocation
!
      CALL discus_alloc_default
!                                                                       
!     Call initialization routine.                                      
!                                                                       
CALL discus_initarrays 
CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
CALL appl_env 
!                                                                       
!     try to read default file                                          
!                                                                       
CALL discus_autodef 
!                                                                       
!     Check for command line parameters                                 
!                                                                       
CALL cmdline_args 
!
CALL no_error
lsetup_done = .true.
!
 1000 FORMAT (/,10x,59('*'),/,                                      &
              10x,'*',15x,'D I S C U S   Version ',a10,10x,'*',/,   &
              10x,'*',57(' '),'*',/                                 &
     &        10x,'*         Created : ',a35,3x,'*',/,              &
              10x,'*',57('-'),'*',/,                                &
     &        10x,'* (c) R.B. Neder  ',                             &
     &        '(reinhard.neder@fau.der)                *',/,        &
     &        10x,'*     Th. Proffen ',                             &
     &        '(tproffen@ornl.gov)                     *',/,        &
     &        10x,59('*'),/,                                        &
              10x,'*',57(' '),'*',/,                                &
     &        10x,'* For information on current changes',           &
     &            ' type: help News',6x,'*',/,                      &
     &        10x,'*',57(' '),'*',/,10x,59('*'),/                   &
     &                     )                                            
END SUBROUTINE setup                          
!
SUBROUTINE discus_set_sub
!
! Sets the specific DIFFEV interfaces four routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE discus_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE discus_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_errlist_appl
   END SUBROUTINE discus_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE discus_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE discus_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE discus_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE discus_validate_var_spec 
END INTERFACE

!
p_mache_kdo         => discus_mache_kdo
p_errlist_appl      => discus_errlist_appl
p_ersetz_para       => discus_ersetz_para
p_upd_para          => discus_upd_para
p_calc_intr_spec    => discus_calc_intr_spec
p_validate_var_spec => discus_validate_var_spec
!
END SUBROUTINE discus_set_sub
!
END MODULE setup_mod
