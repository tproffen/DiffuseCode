MODULE sup_mod
!
CONTAINS
!****7***************************************************************** 
!                                                                       
!     This file contains several subroutines for command language       
!     handling and other common support routines.                       
!                                                                       
!*****7*****************************************************************
!*****7*****************************************************************
!
SUBROUTINE get_cmd (line, ll, befehl, lbef, zeile, lp, prom) 
!+                                                                      
!     This subroutine gets a command for processing. If it              
!     is keyboard input and the program was compiled with               
!     READLINE defined, you will have basic line editing                
!     functions.                                                        
!-                                                                      
USE charact_mod
USE debug_mod 
USE do_execute_mod
USE doact_mod 
USE errlist_mod 
USE jsu_readline
USE learn_mod 
USE lib_length
USE lib_macro_func
use macro_mod
USE class_macro_internal 
USE precision_mod
USE prompt_mod 
USE support_mod
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER  :: ILRN = 33
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line
INTEGER          , INTENT(OUT)   :: ll
CHARACTER (LEN=*), INTENT(OUT)   :: befehl
CHARACTER (LEN=*), INTENT(OUT)   :: zeile
INTEGER          , INTENT(OUT)   :: lp
CHARACTER (LEN=*), INTENT(IN )   :: prom 
!
CHARACTER(LEN=PREC_STRING)       :: input
CHARACTER(LEN=PREC_STRING)       :: as_typed
CHARACTER(60) bprom 
INTEGER :: lbef, indxb , indxt
INTEGER il, jl
INTEGER :: lt
LOGICAL lreg 
LOGICAL :: ldone
!                                                                       
INTEGER socket_accept 
INTEGER socket_get 
INTEGER socket_send 
!                                                                       
ldone = .FALSE.      ! Learn has not yet been written
input  = ' ' 
line   = ' ' 
zeile  = ' ' 
befehl = ' ' 
!                                                                       
ll   = 0 
lp   = 0 
lbef = 0
ier_num = 0
ier_typ = ER_NONE 
!                                                                       
IF (lblock) THEN 
   CALL do_execute (lreg, input, ll) 
   IF (ier_num.ne.0.or..not.lreg) RETURN 
ELSEIF (lmakro .and. .not.lblock_dbg) THEN 
   CALL macro_read (input, ll) 
   IF (ier_num.ne.0) RETURN 
elseif(.not.lblock .and. lmakro .and. lblock_dbg .and. lmakro_dbg) then
   CALL macro_read (input, ll) 
   IF (ier_num.ne.0) RETURN 
ELSE 
   ier_ctrlc = .FALSE.
   ier_rep   = .FALSE.
!                                                                       
!     --Normal mode, if status is PROMPT_OFF or PROMPT_REDIRECT         
!---- --we assume non interactive input and use 'normal' FORTRAN      
!---- --READ to avoid EOF problems.                                   
!                                                                       
!        On old Red Hat systems uses of readline for the first input
!        line caused problems. Seems not to be an issue any longer
   IF(linteractive) THEN
      first_input = .false.
      input = ' '
      ll    = 0
      IF (prompt_status.eq.PROMPT_ON.and..not.first_input) THEN 
         bprom = ' '//prom (1:len_str(prom)) //' > ' 
!                                                                       
!     ----call the c-routine that enables command history & line editing
!                                                                       
         CALL iso_readline (input,bprom) 
         ll=len_str(input)
         CALL learn_save(input, ll, llearn, lmakro, ILRN)
         CALL remove_comment(input, ll) 
         IF(ll>0) THEN
            is_cont1: DO WHILE(input(ll:ll)=='&')
               ll = ll - 1
               bprom = ' '//prom (1:len_str(prom)) //'/continuation > '
               CALL iso_readline(as_typed, bprom)
               lt = len_str(as_typed)
               CALL learn_save(as_typed, lt, llearn, lmakro, ILRN)
               CALL remove_comment(as_typed, lt) 
               input = input(1:ll) // ' ' // as_typed(1:lt)
               ll = ll + 1 + lt
               IF(ll>LEN(input)-100) THEN
                  ier_num = -16
                  ier_typ = ER_COMM
                  RETURN
               ENDIF
            ENDDO is_cont1
         ENDIF
!                                                                       
!------ --otherwise use normal READ                                     
!                                                                       
      ELSE 
         CALL do_prompt (prom) 
         READ ( *, 2000, end = 990, err = 995) input 
         first_input = .FALSE. 
!
         ll=len_str(input)
         CALL learn_save(input, ll, llearn, lmakro, ILRN)
         CALL remove_comment (input, ll) 
         IF(ll>0) THEN
            is_cont2: DO WHILE(input(ll:ll)=='&')
               ll = ll - 1
               bprom = ' '//prom (1:len_str(prom)) //'/continuation > '
               CALL do_prompt (bprom) 
               READ ( *, 2000, end = 990, err = 995) as_typed 
               lt = len_str(as_typed)
               CALL learn_save(as_typed, lt, llearn, lmakro, ILRN)
               CALL remove_comment(as_typed, lt) 
               input = input(1:ll) // ' ' // as_typed(1:lt)
               ll = ll + 1 + lt
               IF(ll>LEN(input)-100) THEN
                  ier_num = -16
                  ier_typ = ER_COMM
                  RETURN
               ENDIF
            ENDDO is_cont2
         ENDIF
      ENDIF 
   ELSE 
      input = input_gui
   ENDIF 
!                                                                       
   ll = len_str (input) 
   IF (prompt_status.eq.PROMPT_REDIRECT) THEN 
      IF (ll.gt.0) THEN 
         WRITE (output_io, 2000) input (1:ll) 
      ELSE 
         WRITE (output_io, 2000) 
      ENDIF 
   ENDIF 
ENDIF 
!                                                                       
!     For commands of significant length, remove leading 'white'        
!     blanks, get the command.                                          
!     Return comment lines "# ..." unchanged.                           
!                                                                       
!
!     - Learn command for learn sequence                                
!
IF (ll.ne.0) THEN 
   IF (input (1:1) .ne.'#'.and.input (1:1) .ne.'!') THEN 
      ll = len_str (input) 
      CALL remove_comment (input, ll) 
!                                                                       
!------ --- Remove leading blanks from 'line'                           
!                                                                       
      il = 1 
      jl = len_str (input) 
      DO while ( (input (il:il) .eq.' '.or.input (il:il) .eq.TAB) &
                 .and.il.le.jl)                                              
         il = il + 1 
      ENDDO 
      line = input (il:jl) 
      ll = len_str (line) 
!                                                                       
!     - The maximum number of significant characters depends on the     
!     - length of the character constant befehl.                        
!                                                                       
      lbef = len (befehl) 
      indxb = index (line, ' ') 
      indxt = index (line, TAB) 
      if(indxb==0) indxb = len(line) + 2      ! If no blank, set beyond length of befehl
      if(indxt==0) indxt = len(line) + 2      ! If no tab,   set beyond length of befehl
      lbef = min (indxb - 1, indxt-1, lbef) 
      befehl = line (1:lbef) 
!                                                                       
!     - command parameters start at the first character following       
!------ - the blank                                                     
!                                                                       
      indxb = min(indxb, indxt)   ! Both have been set to non zero values
      IF (indxb + 1.le.ll) THEN 
         zeile = line (indxb + 1:ll) 
         lp = ll - indxb 
      ENDIF 
   ELSE 
      line = input 
   ENDIF 
ENDIF 
!                                                                       
!     Normal return                                                     
!                                                                       
RETURN 
!                                                                       
!     EOF in input                                                      
!                                                                       
  990 CONTINUE 
      line = 'exit' 
      ll = 4 
      befehl = 'exit' 
      lbef = 4 
      zeile = ' ' 
      lp = 0 
      RETURN 
!                                                                       
!     Error in input                                                    
!                                                                       
  995 CONTINUE 
      ier_num = - 9 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 2000 FORMAT     (a) 
!
END SUBROUTINE get_cmd                        
!
!*****7*****************************************************************
!
SUBROUTINE learn_save(input, ll, llearn, lmakro, ILRN)
!
USE lib_length
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: input
INTEGER         , INTENT(IN) :: ll
LOGICAL         , INTENT(IN) :: llearn
LOGICAL         , INTENT(IN) :: lmakro
INTEGER         , INTENT(IN) :: ILRN
!
INTEGER :: lbef
!
lbef = MIN(LEN(input),LEN_TRIM(input)) 
IF (llearn.and..not.str_comp (input, 'lend', 3, lbef, 4)   &
          .and..not.str_comp (input, 'mouse', 3, lbef, 5)  &
          .and..not.lmakro) THEN
   IF(ll.gt.0) THEN 
      WRITE(ILRN, '(a)') input(1:len_str (input) ) 
   ELSE 
      WRITE(ILRN, '(a)') 
   ENDIF 
ENDIF 
!
END SUBROUTINE learn_save
!
!*****7*****************************************************************
!
SUBROUTINE remove_comment (line, ll) 
!                                                                       
!     removes trailing in line comments                                 
!                                                                       
USE charact_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(INOUT) ::ll 
!                                                                       
INTEGER :: i 
LOGICAL :: quote 
LOGICAL :: search 
!                                                                       
search = .true. 
DO i = 1, ll 
   quote = line (i:i) .eq.'"'.or.line (i:i) .eq.'''' 
   IF (quote) THEN 
      search = .not.search 
   ENDIF 
   IF (search) THEN 
      IF (line (i:i) .eq.'#'.or.line (i:i) .eq.'!') THEN 
         line (i:ll) = ' ' 
         ll = i - 1 
      ENDIF 
   ENDIF 
ENDDO 
ll = LEN_TRIM(line)
!                                                                       
END SUBROUTINE remove_comment                 
!                                                                       
!*****7***********************************************************      
!
SUBROUTINE do_prompt (prom) 
!*                                                                      
!     This routine prints the prompt on the screen                      
!-                                                                      
USE lib_length
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: prom 
!
!                                                                       
IF(prompt_status == PROMPT_ON.OR.prompt_status == PROMPT_REDIRECT) THEN
   WRITE (output_io, '(1X,A,'' > '')',advance='no') prom (1:len_str (prom) ) 
ENDIF 
!                                                                       
END SUBROUTINE do_prompt
!
!*****7***********************************************************      
!
END MODULE sup_mod
