MODULE discus_kdo_common_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE discus_kdo_common(befehl, lbef, line, length, zeile, lp, hlp_sec, &
                             lend, success)
!
! Runs the common set of commands that occur in any menu
!
USE discus_show_menu
USE class_macro_internal
USE errlist_mod
USE do_eval_mod
USE do_show_mod
USE do_wait_mod
USE get_params_mod
USE sup_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: befehl    ! The command to execute
INTEGER         , INTENT(IN)    :: lbef      ! length of command word
CHARACTER(LEN=*), INTENT(INOUT) :: line      ! Command line for macro processing
INTEGER         , INTENT(INOUT) :: length    ! length of command line
CHARACTER(LEN=*), INTENT(INOUT) :: zeile     ! Arguments after "befehl"
INTEGER         , INTENT(INOUT) :: lp        ! Length of argument string
CHARACTER(LEN=*), INTENT(IN)    :: hlp_sec   ! help entry to be looked up
LOGICAL         , INTENT(OUT)   :: lend      ! The "befehl" was "exit"
LOGICAL         , INTENT(OUT)   :: success   ! The "befehl" corresponded to a general command
!
INTEGER, PARAMETER :: MAXP = 5
!
CHARACTER(LEN=1024) :: string
CHARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
REAL               , DIMENSION(MAXP) :: werte
INTEGER            , DIMENSION(MAXP) :: lpara
INTEGER                              :: ianz
!
LOGICAL :: str_comp
!
success = .FALSE.
lend    = .FALSE.
!
!----execute a macro file
!
IF(befehl(1:1) == '@') THEN
   IF(length >= 2) THEN
      CALL file_kdo(line(2:length), length - 1)
      success = .TRUE.
   ELSE
      ier_num = -13
      ier_typ = ER_MAC
   ENDIF
!
!----list asymmetric unit 'asym'                                   
!
ELSEIF(str_comp(befehl, 'asym', 2, lbef, 4) ) THEN
   CALL show_asym
   success = .TRUE.
!
!----continues a macro 'continue'
!
ELSEIF(str_comp (befehl, 'continue', 2, lbef, 8) ) THEN
   CALL macro_continue (zeile, lp)
   success = .TRUE.
!                                                                       
!     ----list atoms present in the crystal 'chem'                      
!                                                                       
ELSEIF(str_comp (befehl, 'chem', 2, lbef, 4) ) THEN
   CALL show_chem
   success = .TRUE.
!
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!
ELSEIF(str_comp (befehl, 'echo', 2, lbef, 4) ) THEN
   CALL echo (zeile, lp)
   success = .TRUE.
!
!---Evaluate an expression, just for interactive check 'eval'     
!
ELSEIF(str_comp (befehl, 'eval', 2, lbef, 4) ) THEN
   CALL do_eval (zeile, lp, .TRUE.)
   success = .TRUE.
!
!     ----exit 'exit'                                                   
!
ELSEIF(str_comp (befehl, 'exit', 2, lbef, 4) ) THEN
   lend = .true.
   success = .TRUE.
!
!     ----help 'help','?'                                               
!
ELSEIF(str_comp (befehl, 'help', 2, lbef, 4) .OR.         &
       str_comp (befehl, '?   ', 1, lbef, 4)      ) THEN
   IF(str_comp (zeile, 'errors', 2, lp, 6) ) THEN
      lp = lp + 7
      CALL do_hel ('discus '//zeile, lp)
      success = .TRUE.
   ELSE
      string = 'discus ' // hlp_sec // ' ' // zeile
      lp = LEN_TRIM(string)
      CALL do_hel (string, lp)
      success = .TRUE.
   ENDIF
!
!------ - Do a generic show
!
ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN
   CALL get_params (zeile, ianz, cpara, lpara, MAXP, lp)
   IF(ier_num==0) THEN
      IF(ianz > 0) THEN
         CALL do_show_generic (cpara, lpara, MAXP)
         success = .TRUE.
      ENDIF
   ENDIF
!
!------- -Operating System Kommandos 'syst'                             
!
ELSEIF(str_comp (befehl, 'syst', 2, lbef, 4) ) THEN
   IF (zeile /= ' ') THEN
      CALL do_operating (zeile (1:lp), lp)
      success = .TRUE.
   ELSE
      ier_num = - 6
      ier_typ = ER_COMM
   ENDIF
!
!-----waiting for user input                                    
!
ELSEIF(str_comp (befehl, 'wait', 3, lbef, 4) ) THEN
   CALL do_input (zeile, lp)
   success = .TRUE.
ENDIF
!
END SUBROUTINE discus_kdo_common
!
!*******************************************************************************
!
END MODULE discus_kdo_common_mod
