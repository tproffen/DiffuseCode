
#include "debug.h"

MODULE setup_mod
!
CONTAINS
!
!
SUBROUTINE setup 
!                                                                       
!     This routine makes inital setup of DISCUS                         
!                                                                       
      USE allocate_appl_mod
      USE init_mod
!
      USE errlist_mod
      USE prompt_mod 
!
IMPLICIT none 
!                                                                       
      include'date.inc' 

MSG('setup_mod.setup()...')

      pname      = 'discus'
      pname_cap  = 'DISCUS'
!                                                                       
      blank   = ' '
      prompt  = pname
      prompt_status = PROMPT_ON
      prompt_status_old = PROMPT_ON
write (*,*) 'lsetup_done: ', lsetup_done
CALL ini_ran (0)
write (*,*) 'lsetup_done: ', lsetup_done
!                                                                       
!------ Write starting screen                                           
!                                                                       
version = aktuell 
WRITE ( *, 1000) version, cdate 
!
!     Call initial default allocation
!
    MSG('calling alloc_default...')

      CALL alloc_default
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
CALL no_error
VAR(lsetup_done)
lsetup_done = .true.
VAR(lsetup_done)
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
END MODULE setup_mod
