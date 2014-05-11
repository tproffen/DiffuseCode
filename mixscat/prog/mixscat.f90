      PROGRAM mixsca 
!                                                                       
      USE doact_mod
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal 
      USE prompt_mod 
      IMPLICIT none 
!*****7*****************************************************************
!                                                                       
!     Main program for MIXSCAT                                          
!                                                                       
!*****7*****************************************************************
!                                                                       
      CHARACTER(1024) line, zeile 
      CHARACTER(4) befehl 
      LOGICAL lend 
      INTEGER laenge, lp, lbef 
!
      EXTERNAL :: mixscat_mache_kdo
!                                                                       
      pname             = 'mixscat'
      pname_cap         = 'MIXSCAT'
!                                                                       
      lend = .false. 
      blank = ' ' 
      prompt = pname 
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
      DO while (.not.lend) 
      CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0.and.laenge.gt.0) then 
!                                                                       
!     - If not a comment continue                                       
!                                                                       
         IF (.not. (line (1:1) .eq.'#'.or.line (1:1) .eq.'!') ) then 
!                                                                       
!     - execute command                                                 
!                                                                       
            IF (line (1:3) .eq.'do '.OR.line (1:2) .eq.'if') then 
               CALL do_loop (line, lend, laenge, mixscat_mache_kdo) 
            ELSE 
               CALL mixscat_mache_kdo (line, lend, laenge) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!     - Handle error message                                            
!                                                                       
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (lmakro) then 
            CALL macro_close 
            prompt_status = PROMPT_ON 
         ENDIF 
         lblock = .false. 
         CALL no_error 
      ENDIF 
      ENDDO 
!                                                                       
 9999 CONTINUE 
!                                                                       
      IF (output_io.ne.OUTPUT_SCREEN) then 
         CLOSE (output_io) 
      ENDIF 
!                                                                       
      END PROGRAM mixsca                            
!*****7*****************************************************************
      SUBROUTINE setup 
!                                                                       
!     This routine makes inital setup                                   
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      include'date.inc' 
!                                                                       
      CALL ini_ran (0) 
!                                                                       
!------ Write starting screen                                           
!                                                                       
      version=aktuell
      WRITE ( *, 1000) version, cdate 
!                                                                       
!     Call initialization routine                                       
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
!     try to read command line arguments                                
!                                                                       
      CALL cmdline_args 
!                                                                       
 1000 FORMAT   (/,10x,59('*'),/,10x,'*',15x,                            &
     &            'M I X S C A T   Version ',                           &
     &            a10,8x,'*',/,10x,'*',57(' '),'*',/                    &
     &            10x,'*         Created : ',a35,3x,'*',/,10x,'*',      &
     &            57('-'),'*',/,                                        &
     &                 10x,'* by Caroline Wurden, Kathar',              &
     &            'ine Page, Anna Llobet and     *',/,                  &
     &                 10x,'*    Thomas Proffen - Lujan ',              &
     &            ' Neutron Scattering Center    *',/,                  &
     &            10x,59('*'),/)                                        
      END SUBROUTINE setup                          
