module exp2pdf_menu
!
CONTAINS
!
SUBROUTINE exp2pdf 
!                                                                       
!+                                                                      
!  Convert an experimental powder diffraction pattern to a PDF
!-                                                                      
use discus_config_mod 
use exp2pdf_load_mod
use exp2pdf_supp_mod
use exp2pdf_run_mod
!
use calc_expr_mod
use doact_mod 
use errlist_mod 
use kdo_all_mod
use lib_errlist_func
use lib_length
use lib_macro_func
use class_macro_internal
use prompt_mod 
use sup_mod
use precision_mod
use str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MIN_PARA = 99  ! A command requires at least these no of parameters
INTEGER :: maxw 
!                                                                       
CHARACTER(len=5)           :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=PREC_STRING) :: line, zeile
INTEGER :: lp, length, lbef 
INTEGER :: indxg
LOGICAL :: lend
!
!                                                                       
maxw = MIN_PARA
lend = .false. 
CALL no_error 
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/exp2pdf' 
!                                                                       
loop_main: DO while (.not.lend)                 ! Main exp2pdf loop
!
   if_error: IF (ier_num.ne.0) THEN 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in exp2pdf menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close 
                  prompt_status = PROMPT_ON 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (lblock) THEN 
            ier_num = - 11 
            ier_typ = ER_COMM 
            prompt_status = PROMPT_ON 
            prompt = orig_prompt
            RETURN 
         ENDIF 
         CALL no_error 
         lmakro_error = .FALSE.
         sprompt = ' '
      ENDIF 
   ENDIF if_error
!
   CALL get_cmd(line, length, befehl, lbef, zeile, lp, prompt) 
   IF (ier_num /= 0) cycle loop_main
!
   IF(line == ' '      .or. line(1:1) == '#' .or. &
      line == char(13) .or. line(1:1) == '!'        ) cycle loop_main
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=') 
   IF (indxg.ne.0                                              &
        .AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
        .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
        .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                    str_comp (befehl, '?   ', 2, lbef, 4) )    &
        .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluate an expression and assign the value to a variable   
!                                                                       
            CALL do_math(line, indxg, length) 
      cycle loop_main
   endif
!
!--EXP2PDF specific commands
!
!  --- Experimental data file
!
   if(str_comp (befehl, 'back', 2, lbef, 4)) then
      call exp2pdf_load(2, zeile, lp)
!
!  --- Experimental composition 'comp'
!
   elseif(str_comp (befehl, 'comp', 2, lbef, 4)) then
      call exp2pdf_composition(zeile, lp)
!
!  --- Experimental data file
!
   elseif(str_comp (befehl, 'data', 2, lbef, 4)) then
      call exp2pdf_load(1, zeile, lp)
!
!     ----exit 'exit'
!
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then
      lend = .true.
!
!     ----Set output file and PDF range
!
   ELSEIF (str_comp (befehl, 'outf', 2, lbef, 4) .or.              &
           str_comp (befehl, 'outputfile', 2, lbef, 10) ) then
      call exp2pdf_outfile(zeile, lp)
!
!     ----Set radiation 'radiation'
!
   ELSEIF (str_comp (befehl, 'poly', 2, lbef, 4) ) then
      call exp2pdf_poly(zeile, lp)
!
!     ----Set radiation 'radiation'
!
   ELSEIF (str_comp (befehl, 'radiation', 2, lbef, 9) ) then
      call exp2pdf_radiation(zeile, lp)
!
!     ----Reset         'reset'     
!
   ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) then
      call exp2pdf_reset
!
!     ----perform task 'run'
!
   ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then
      call exp2pdf_run(zeile, lp)
!
!     ----Set qmax limits 'qmax'
!
   ELSEIF (str_comp (befehl, 'limits', 2, lbef, 6) ) then
      call exp2pdf_qmax(zeile, lp)
   ELSE 
!                                                                       
!     ----Unknown exp2pdf command, try general commands
!                                                                       
      call kdo_all(befehl, lbef, zeile, lp)
   ENDIF 
!
ENDDO loop_main
!
prompt = orig_prompt
!                                                                       
END SUBROUTINE exp2pdf                          
!
!*******************************************************************************
!
end module exp2pdf_menu
