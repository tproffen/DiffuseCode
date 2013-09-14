PROGRAM discus 
!                                                                       
   USE doact_mod 
   USE errlist_mod 
   USE learn_mod 
   USE macro_mod 
   USE prompt_mod 

   IMPLICIT none 
!
!*****7*****************************************************************
!                                                                       
!     Main program for DISCUS                                           
!                                                                       
!     This is the main program for DISCUS. It sets up most              
!     variables and calls the loop interpreting the commands.           
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)      
!               Th. Proffen (tproffen@ornl.gov)                         
!                                                                       
!*****7*****************************************************************
!                                                                       
CHARACTER (LEN=1024) :: line, zeile 
CHARACTER (LEN=4)    :: befehl 
LOGICAL              :: lend 
INTEGER              :: laenge, lp, lbef 
!                                                                       
pname      = 'discus'
pname_cap  = 'DISCUS'
!                                                                       
lend    = .false. 
blank   = ' ' 
prompt  = pname 
prompt_status = PROMPT_ON 
prompt_status_old = PROMPT_ON 
!                                                                       
!------ Setting up variables and print start screen                     
!                                                                       
CALL setup 
CALL no_error 
!                                                                       
!------ This is the main loop: reading commands ..                      
!                                                                       
main: DO while (.not.lend) 
   CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt) 
ok: IF (ier_num.eq.0.and.laenge.gt.0) then 
!                                                                       
!     - If not a comment continue                                       
!                                                                       
      IF (.not. (line (1:1) .eq.'#'.or.line (1:1) .eq.'!') ) then 
!                                                                       
!     - execute command                                                 
!                                                                       
         IF (line (1:3) .eq.'do '.OR.line (1:2) .eq.'if') then 
            CALL do_loop (line, lend, laenge) 
         ELSE 
            CALL mache_kdo (line, lend, laenge) 
         ENDIF 
      ENDIF 
   ENDIF ok
!                                                                       
!     - Handle error message                                            
!                                                                       
fehler: IF (ier_num.ne.0) then 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) then 
         IF (lmakro) then 
            CALL macro_close 
            prompt_status = PROMPT_ON 
         ENDIF 
         lblock = .false. 
         CALL no_error 
      ENDIF 
   ENDIF fehler
ENDDO main
!                                                                       
 9999 CONTINUE 
CALL do_exit 
!                                                                       
END PROGRAM discus                            
!*****7*****************************************************************
SUBROUTINE setup 
!                                                                       
!     This routine makes inital setup of DISCUS                         
!                                                                       
      USE allocate_appl_mod
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
version = aktuell 
WRITE ( *, 1000) version, cdate 
!
!     Call initial default allocation
!
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
 1000 FORMAT (/,10x,59('*'),/,                                      &
              10x,'*',15x,'D I S C U S   Version ',a10,10x,'*',/,   &
              10x,'*',57(' '),'*',/                                 &
     &        10x,'*         Created : ',a35,3x,'*',/,              &
              10x,'*',57('-'),'*',/,                                &
     &        10x,'* (c) R.B. Neder  ',                             &
     &        '(reinhard.neder@krist.uni-erlangen.de)  *',/,        &
     &        10x,'*     Th. Proffen ',                             &
     &        '(tproffen@ornl.gov)                     *',/,        &
     &        10x,59('*'),/,                                        &
              10x,'*',57(' '),'*',/,                                &
     &        10x,'* For information on current changes',           &
     &            ' type: help News',6x,'*',/,                      &
     &        10x,'*',57(' '),'*',/,10x,59('*'),/                   &
     &                     )                                            
END SUBROUTINE setup                          
