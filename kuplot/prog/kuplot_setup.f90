MODULE setup_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE setup 
!                                                                       
!     This routine makes inital setup of KUPLOT                         
!                                                                       
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
      IMPLICIT none 
!                                                                       
      include'date.inc' 
!                                                                       
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
!                                                                       
      prompt = pname 
      blank = ' ' 
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
!     Call initialization routines                                      
!                                                                       
      CALL initarrays 
      CALL init_sysarrays 
      CALL appl_env 
      CALL auto_def 
      CALL cmdline_args 
      CALL no_error
!
      lsetup_done = .true.
!                                                                       
 1000 FORMAT (/,10x,59('*'),/,10x,'*',15x,                              &
     &        'K U P L O T   Version ',                                 &
     &        a10,10x,'*',/,10x,'*',57x,'*',/                           &
     &        10x,'*         Created : ',a35,3x,'*',/,10x,'*',          &
     &        57('-'),'*',/,10x,'* (c) Th. Proffen ',                   &
     &        '(tproffen@ornl.gov)                     *',/,            &
     &        10x,'*     R.B. Neder  ',                                 &
     &        '(reinhard.neder@fau.de)                 *',/,            &
     &        10x,59('*'),/,                                            &
     &        10x,'* GSAS code: (c) Allen C. Larson and',               &
     &        ' Robert B. Von Dreele *',/,                              &
     &        10x,59('*'),/)                                            
      END SUBROUTINE setup                          
END MODULE setup_mod
