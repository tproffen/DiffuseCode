PROGRAM mixsca 
!                                                                       
      USE doact_mod
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal 
      USE prompt_mod 
      USE lib_f90_default_mod
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
      CALL lib_alloc_default
      CALL setup 
      CALL mixscat_set_sub
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
               CALL do_loop (line, lend, laenge) 
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
!     Call initialization routine                                       
!                                                                       
      CALL mixscat_initarrays 
      CALL init_sysarrays 
!                                                                       
!     get envirmonment information                                      
!                                                                       
      CALL appl_env (.TRUE.,0)
!                                                                       
!------ Write starting screen                                           
!                                                                       
      version=aktuell
      WRITE ( *, 1000) version, cdate 
      CALL write_appl_env (.TRUE.,0)
!                                                                       
!     try to read default file                                          
!                                                                       
      CALL mixscat_autodef 
!                                                                       
!     try to read command line arguments                                
!                                                                       
      CALL cmdline_args (0)
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
!
SUBROUTINE mixscat_set_sub
!
! Sets the specific DIFFEV interfaces four routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE mixscat_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE mixscat_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE mixscat_errlist_appl
   END SUBROUTINE mixscat_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE mixscat_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE mixscat_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE mixscat_upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE mixscat_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE mixscat_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE mixscat_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE mixscat_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE mixscat_validate_var_spec 
END INTERFACE

!
p_mache_kdo         => mixscat_mache_kdo
p_errlist_appl      => mixscat_errlist_appl
p_ersetz_para       => mixscat_ersetz_para
p_upd_para          => mixscat_upd_para
p_calc_intr_spec    => mixscat_calc_intr_spec
p_validate_var_spec => mixscat_validate_var_spec
!
END SUBROUTINE mixscat_set_sub
!
