      PROGRAM kuplot 
!                                                                       
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE macro_mod 
      USE prompt_mod 
      IMPLICIT none 
!*****7*****************************************************************
!       This is the universal plot program KUPLOT. It sets up most      
!     variables and calls the loop interpreting the commands.           
!*****7*****************************************************************
!                                                                       
!                                                                       
      CHARACTER(1024) line, zeile 
      CHARACTER(4) befehl 
      LOGICAL lend 
      INTEGER lbef, lp, ll 
!                                                                       
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
!                                                                       
      prompt = pname 
      lend = .false. 
      blank = ' ' 
      prompt_status = PROMPT_ON 
      prompt_status_old = PROMPT_ON 
!                                                                       
      CALL setup 
      CALL no_error 
!                                                                       
      DO while (.not.lend) 
   10 CONTINUE 
      CALL get_cmd (line, ll, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0.and.ll.gt.0) then 
         IF (line.eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!------ --- Execute command                                             
!                                                                       
         IF (befehl (1:3) .eq.'do '.or.befehl (1:2) .eq.'if') then 
            CALL do_loop (line, lend, ll) 
         ELSE 
            CALL mache_kdo (line, lend, ll) 
         ENDIF 
      ENDIF 
!                                                                       
!     - Handle error message                                            
!                                                                       
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            lblock = .false. 
            CALL no_error 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!------ END of PROGRAM                                                  
!                                                                       
 9999 CONTINUE 
      CALL do_exit 
!                                                                       
 2000 FORMAT    (a) 
!                                                                       
      END PROGRAM kuplot                            
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
!                                                                       
 1000 FORMAT (/,10x,59('*'),/,10x,'*',15x,                              &
     &        'K U P L O T   Version ',                                 &
     &        a10,10x,'*',/,10x,'*',57x,'*',/                           &
     &        10x,'*         Created : ',a35,3x,'*',/,10x,'*',          &
     &        57('-'),'*',/,10x,'* (c) Th. Proffen ',                   &
     &        '(tproffen@ornl.gov)                     *',/,            &
     &        10x,'*     R.B. Neder  ',                                 &
     &        '(reinhard.neder@krist.uni-erlangen.de)  *',/,            &
     &        10x,59('*'),/,                                            &
     &        10x,'* GSAS code: (c) Allen C. Larson and',               &
     &        ' Robert B. Von Dreele *',/,                              &
     &        10x,59('*'),/)                                            
      END SUBROUTINE setup                          
