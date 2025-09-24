MODULE do_if_mod
!
private
!
public do_loop
!
CONTAINS
!+                                                                      
!     Handling of blockstructures and all math functions                
!     and parameter updating and referencing is done by                 
!     routines in this file.                                            
!                                                                       
!*****7*****************************************************************
SUBROUTINE do_loop (line, lend, length) 
!+                                                                      
!     All commands included in a block structure are read and stored in 
!     character array. Executable commands are parsed to mach_kdo.      
!     If an illegal nesting of block structures is detected, an         
!     error flag is returned and the block structure is not executed    
!     at all.                                                           
!-                                                                      
USE doact_mod 
USE doexec_mod 
USE do_execute_mod
USE errlist_mod 
USE lib_macro_func
USE class_macro_internal 
!                                                                       
IMPLICIT none 
!
CHARACTER(LEN=*), INTENT(INOUT) :: line 
LOGICAL         , INTENT(INOUT) :: lend
INTEGER         , INTENT(INOUT) :: length
!                                                                       
!-----      read first block structure command                          
!                                                                       
CALL do_do_init(line, length)
!
cond_error: IF(ier_num ==0 )then
!                                                                       
!.....read all commands                                                 
!                                                                       
   lblock_read = .true. 
   DO WHILE (level.gt. - 1) 
      CALL do_insert_line!(jlevel) ! Moved to separate subroutine
      IF(ier_num /= 0) exit cond_error
   ENDDO 
      lblock_read = .false. 
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
!-----      execute the block structure                                 
!                                                                       
   CALL do_execute_block(lend)
!
endif cond_error
!                                                                       
lblock = .false. 
IF (ier_num.ne.0) THEN 
   IF(lmakro .AND. lmacro_close) THEN
      CALL macro_close(-1)
   ENDIF 
ENDIF 
do_macro = ' '
!                                                                       
!                                                                       
!2000 FORMAT    (a1) 
!3000 FORMAT    ('Erroneous line in block structure') 
!3100 FORMAT    (a41) 
END SUBROUTINE do_loop                        
!
!*****7*****************************************************************
!
SUBROUTINE do_do_init(line, length)
!                                                                       
!-----      read first block structure command                          
!                                                                       
USE doact_mod 
USE doexec_mod 
!
USE blanks_mod
USE errlist_mod 
USE class_macro_internal 
USE lib_length
USE lib_macro_func
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line 
!OGICAL         , INTENT(INOUT) :: lend 
INTEGER         , INTENT(INOUT) :: length
! 
CHARACTER(LEN=20) :: prom 
CHARACTER(LEN=3)  :: cprom (0:3) 
INTEGER :: i
!
!
DATA cprom / '/do', '/if', '/do', '/do' / 
!                                                                       
ier_num = 0
!
      IF (line (1:2) .eq.'do'.and.INDEX (line, 'while') .ne.0) then 
         jlevel (0) = 2 
         i = length - 3 
!        CALL rem_bl (line (4:length), i) 
         CALL rem_insig_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'do'.and.INDEX (line, '=') .ne.0) then 
         jlevel (0) = 0 
         i = length - 3 
         CALL rem_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'do'.and.length.eq.2) then 
         jlevel (0) = 3 
      ELSEIF (line (1:2) .eq.'if'.and.INDEX (line, 'then') .ne.0) then 
         jlevel (0) = 1 
!        CALL rem_bl (line, length) 
         CALL rem_insig_bl (line, length) 
      ELSE 
         ier_num = -31 
         ier_typ = ER_FORT 
         WRITE (ier_msg (1), 3000) 
         WRITE (ier_msg (2), 3100) line (1:41) 
      ENDIF 
!
cond_error:if(ier_num==0) then
   prom = prompt (1:len_str (prompt) ) //cprom (jlevel (0) ) 
   DO i = 0, maxlev 
      nlevel (i) = - 1 
   ENDDO 
   level = 0 
   nlevel (level) = 0 
   do_comm (0, 0) = line 
   do_leng (0, 0) = length 
   IF(lmakro) do_macro(0,0) = &
      mac_tree_active%active%macrofile(1:LEN_TRIM(mac_tree_active%active%macrofile))
   lblock_read = .TRUE. 
endif cond_error
!
IF (ier_num.ne.0) then 
   lblock_read = .FALSE. 
   lblock      = .FALSE. 
   IF(lmakro .AND. lmacro_close) THEN
      CALL macro_close(-1)
   ENDIF 
   do_macro = ' '
ENDIF 
!
 3000 FORMAT    ('Erroneous line in block structure') 
 3100 FORMAT    (a41) 
!
END SUBROUTINE do_do_init
!
!*****7*****************************************************************
!
SUBROUTINE do_insert_line!(jlevel)
!
! Inserts an individual line into the array
!
USE ber_params_mod
USE blanks_mod
USE doact_mod 
USE doexec_mod 
USE errlist_mod 
USE get_params_mod
USE learn_mod 
USE lib_length
USE lib_macro_func
USE class_macro_internal 
USE precision_mod
USE prompt_mod 
USE set_sub_generic_mod
USE sup_mod
USE mpi_slave_mod
USE lib_macro_func
USE str_comp_mod
USE take_param_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: MAXW = 5
!
CHARACTER(LEN=PREC_STRING)           :: zeile 
CHARACTER(LEN=PREC_STRING)           :: line 
CHARACTER(LEN=PREC_STRING)           :: outfile
CHARACTER(LEN=20)                    :: prom 
CHARACTER(LEN=4)                     :: befehl 
CHARACTER(LEN=3)                     :: cprom (0:3) 
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
CHARACTER(LEN=PREC_STRING)                  :: prog_n
CHARACTER(LEN=PREC_STRING)                  :: mac_n
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
INTEGER                              :: i, length, lbef, lp
INTEGER                              :: indxm
INTEGER                              :: length_m
INTEGER                              :: out_l
INTEGER                              :: prog_l
INTEGER                              :: mac_l
INTEGER                              :: nindiv
!INTEGER, DIMENSION(0:maxlev)         :: jlevel ! (0:maxlev) 
LOGICAL                              :: repeat
!
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 4
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
DATA oname  / 'partial', 'repeat' , 'logfile', 'compute'  /
DATA loname /  7       ,  6       ,  7       ,  7  /
!
DATA cprom / '/do', '/if', '/do', '/do' / 
!
opara  =  (/ '0.000000', '1.000000', 'none    ', 'parallel' /)   ! Always provide fresh default values
lopara =  (/  8        ,  8        ,  8        ,  8         /)
owerte =  (/  0.0      ,  1.0      ,  0.0      ,  1.0       /)
!
!                                                                       
lblock_read = .true. 
nlevel (level) = nlevel (level) + 1 
IF (nlevel (level) .gt.maxcom) then 
   ier_num = - 15 
   ier_typ = ER_FORT 
ENDIF 
!                                                                       
cond_error: if(ier_num==0) then
   loop_macro: do
!                                                                       
      prom = prompt (1:len_str (prompt) ) //cprom (jlevel (level) ) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF(length==0) THEN
          nlevel (level) = nlevel (level) - 1 
          exit cond_error
      ENDIF
      IF (line (1:2) .eq.'do'.and.INDEX (line, 'while') .ne.0) then 
         i = length - 3 
!        CALL rem_bl (line (4:length), i) 
         CALL rem_insig_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'do'.and.INDEX (line, '=') .ne.0) then 
         i = length - 3 
         CALL rem_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'if') then 
!        CALL rem_bl (line, length) 
         CALL rem_insig_bl (line, length) 
      ENDIF 
!                                                                       
!------ execute a macro file                                            
!                                                                       
      IF (line (1:1) .eq.'@') then 
         ier_num = 0 
         ier_typ = ER_NONE 
         IF (length.ge.2) then 
            line(1:length-1) = line(2:length)
            length = length - 1
            CALL file_kdo (line, length) 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_MAC 
            exit cond_error
         ENDIF 
      else
         exit loop_macro
      endif
   enddo loop_macro
!
      IF(line(1:6)=='branch') THEN
         indxm = INDEX(line,'-macro')
         IF(indxm>0) THEN   ! there is a macro name specified, 
            zeile = line(indxm+7:length)
            length_m = length - indxm - 6
            CALL file_kdo(zeile, length_m)
            line = line(1:indxm-2)
            length = indxm-2
         ENDIF
      ELSEIF(line(1:7)=='run_mpi' .AND. .NOT.mpi_active ) THEN 
!
!        run  a slave from diffev within discus_suite at NO MPI
!
        IF(pname=='diffev' .OR. pname=='suite') THEN
!
            opara  =  (/ '0.000000', '1.000000', 'none    ', 'parallel' /)   ! Always provide fresh default values
            lopara =  (/  8        ,  8        ,  8        ,  8         /)
            owerte =  (/  0.0      ,  1.0      ,  0.0      ,  1.0       /)
            CALL get_params (zeile, ianz, cpara, lpara, maxw, length)
            CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
            IF(ier_num == 0) THEN
               IF ( ianz >= 2 .and. ianz <= 4 ) THEN
                  IF ( ianz == 4 ) THEN
                     outfile = cpara(4)(1:100)! Target for program output 
                     out_l   = lpara(4)
                     cpara(4) = '0'
                     lpara(4) = 1
                  ELSE
                     IF(opara(3)=='none') THEN
                        outfile = '/dev/null'   ! Default output
                        out_l   = 9
                     ELSE
                        outfile = opara(3)(1:lopara(3))
                        out_l   = lopara(3)
                     ENDIF
                  ENDIF
                  prog_n = cpara(1)
                  prog_l = lpara(1)
                  mac_n  = cpara(2)
                  mac_l  = lpara(2)
                  cpara(1:2) = '0'
                  lpara(1:2) = 1
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  repeat = .false.
                  IF (ier_num.eq.0) then
                     IF(ianz>2) THEN
                     IF ( nint(werte(3)) > 0 ) THEN
                        repeat = .true.             ! repeat = false means no repetition
                     ENDIF
                     nindiv = max(1,nint(werte(3))) ! nindiv is at least 1

                     IF ( nint(werte(3)) > 0 ) THEN
                        repeat = .true.             ! repeat = false means no repetition
                        IF (str_comp (opara(4), 'parallel', 4, lopara(4), 8) ) THEN
                           repeat = .TRUE.
                        ELSE
                           repeat = .FALSE.
                        ENDIF
                     ELSE !
                        repeat = .FALSE.
                     ENDIF
                     nindiv = max(1,nint(werte(3))) ! nindiv is at least 1
                     ELSE
                        IF (str_comp (opara(4), 'parallel', 4, lopara(4), 8) ) THEN
                           repeat = .TRUE.
                        ELSE
                           repeat = .FALSE.
                        ENDIF
                        nindiv = max(1,nint(owerte(2))) ! nindiv is at least 1
                     ENDIF
                  ENDIF
                  CALL p_loop_mpi(prog_n, prog_l, mac_n, mac_l, &
                               outfile, out_l, repeat, nindiv)
                  IF(repeat) THEN
!                      WRITE(zeile,'(a,a)') mac_n(1:mac_l),' ., i[205], i[206]'
!                     length_m = mac_l + 18
                      WRITE(zeile,'(a,a)') mac_n(1:mac_l),' ., REF_KID, REF_INDIV'
                     length_m = mac_l + 22
                  ELSE
!                     WRITE(zeile,'(a,a)') mac_n(1:mac_l),' ., i[205]'
!                     length_m = mac_l + 10
                     WRITE(zeile,'(a,a)') mac_n(1:mac_l),' ., REF_KID'
                     length_m = mac_l + 11
                  ENDIF
   
                  CALL file_kdo(zeile, length_m)
                   level_mpi = level
                  nlevel_mpi = nlevel(level)
                  line = 'run_mpi '//prog_n(1:prog_l)//','//mac_n(1:mac_l)//', DOLOOP'
                  length = 17 + prog_l + mac_l
               ENDIF
            ENDIF
         ELSE
            ier_num = -9
            ier_typ = ER_COMM
            ier_msg(4) = 'Found run_mpi command within a loop, NO MPI'
            ier_msg(5) = 'Loop not started within DIFFEV or SUITE'
            ier_msg(6) = 'Start loop with run_mpi in DIFFEV or SUITE'
         ENDIF
      ENDIF 
!                                                                       
      do_comm (nlevel (level), level) = line 
      do_leng (nlevel (level), level) = length 
      ldostart (level) = .true. 
      IF(lmakro) THEN
         do_macro(nlevel (level), level) =  &
           mac_tree_active%active%macrofile(1:LEN_TRIM(mac_tree_active%active%macrofile))
      ENDIF
!
      IF(line(1:3) ==  'end') THEN
         IF (INDEX (line, 'if') >   0) THEN
            CALL rem_insig_bl(line, length)
         ELSEIF (INDEX (line, 'do') >   0) THEN
            CALL rem_insig_bl(line, length)
         ENDIF
      ENDIF
!                                                                       
!     sort commands into proper levels                                  
!                                                                       
      IF (line (1:3) .eq.'do '.or.line (1:2) .eq.'if') then 
         IF (level.lt.maxlev) then 
            do_comm (nlevel (level) , level) = '&&' 
            do_leng (nlevel (level), level) = 6 
            WRITE (do_comm (nlevel (level) , level) (3:6) , '(i4)')     &
            nlevel (level + 1) + 1                                      
            level = level + 1 
            nlevel (level) = nlevel (level) + 1 
            do_comm (nlevel (level), level) = line 
            do_leng (nlevel (level), level) = length 
            IF (line(1:2).eq.'do'.and.INDEX(line,'while').ne.0) then
               jlevel (level) = 2 
            ELSEIF (line (1:2) .eq.'do'.and.INDEX (line, '=') .ne.0) then 
               jlevel (level) = 0 
            ELSEIF (line (1:2) .eq.'do'.and.length.eq.2) then 
               jlevel (level) = 3 
            ELSEIF (line (1:2).eq.'if'.and.INDEX(line,'then').ne.0) then
               jlevel (level) = 1 
            ELSE 
               ier_num = - 31 
               ier_typ = ER_FORT 
               exit cond_error
            ENDIF 
         ELSE 
            ier_num = - 16 
            ier_typ = ER_FORT 
            exit cond_error
         ENDIF 
      ELSEIF (line (1:5) .eq.'enddo') then 
         IF (jlevel (level) .eq.0.and.length.eq.5) then 
            level = level - 1 
         ELSEIF (jlevel (level) .eq.2.and.length.eq.5) then 
            level = level - 1 
         ELSEIF (jlevel (level) .eq.3.and.INDEX (line, 'until') .ne.0)  then
            level = level - 1 
         ELSE 
            ier_num = - 19 
            ier_typ = ER_FORT 
            exit cond_error
         ENDIF 
      ELSEIF (line (1:5) .eq.'endif') then 
         IF (jlevel (level) .eq.1) then 
            level = level - 1 
         ELSE 
            ier_num = - 19 
            ier_typ = ER_FORT 
            exit cond_error
         ENDIF 
      ENDIF 
endif cond_error

!
IF (ier_num.ne.0) then 
   WRITE (ier_msg (1), 3000) 
   WRITE (ier_msg (2), 3100) line (1:41) 
   IF(lmakro .AND. lmacro_close) THEN
      CALL macro_close(-1)
   ENDIF 
ENDIF 
 3000 FORMAT    ('Erroneous line in block structure') 
 3100 FORMAT    (a41) 
!
END SUBROUTINE do_insert_line
!
!*****7*****************************************************************
!
SUBROUTINE do_execute_block(lend)
!
!  Runs through all commands stored in the block array and executes the loop/if
!
USE doact_mod 
USE doexec_mod 
      USE do_execute_mod
USE doloop_mod 
USE errlist_mod 
USE learn_mod 
USE class_macro_internal 
USE prompt_mod 
USE lib_errlist_func
USE lib_macro_func
use macro_mod
USE set_sub_generic_mod
USE str_comp_mod
USE sup_mod
USE mpi_slave_mod
USE terminal_mod
!                                                                       
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(OUT) :: lend
!                                                                       
CHARACTER(LEN=PREC_STRING) :: line 
CHARACTER(LEN=PREC_STRING)                  :: zeile 
CHARACTER(LEN=4)                     :: befehl 
INTEGER :: length, lp, lbef
LOGICAL :: lreg 
!
!
!      DO i = 0, maxlev 
ilevel   (0:maxlev) = - 1 
ltest    (0:maxlev) = .false. 
ldostart (0:maxlev) = .true. 
!     ENDDO 
level = 0 
lblock = .true. 
!                                                                       
!-----      as long as there are commands and no error ...              
!                                                                       
main: DO WHILE (level.gt. - 1.and. (                                    &
      ier_num.eq.0.or.ier_num.ne.0.and.ier_sta.eq.ER_S_LIVE) )          
!                                                                       
!     increment the command array, follow up with block commands        
!                                                                       
   CALL do_execute (lreg, line, length) 
   IF (ier_num.ne.0.AND.ier_sta.ne.ER_S_LIVE) THEN 
      EXIT main               ! Error finish loop
   ELSEIF (ier_num.ne.0.AND.ier_sta.eq.ER_S_LIVE) THEN 
      CALL errlist 
   ENDIF 
!                                                                       
!  Regular command, call mach_kdo or file_kdo                        
!                                                                       
   IF (lreg) then 
      IF (.NOT. (level.eq.0.AND.ilevel (level) .eq.0) ) THEN 
         IF (str_comp (line (1:4) , 'stop', 4, length, 4) ) THEN 
            WRITE(output_io, '(a,a1)') '------ > do/if halted, continue with cont ...', achar(7) 
            WRITE(output_io, '(a,a1)') '------ > return to the current section/menu  ', achar(7) 
            WRITE(output_io, '(a,a1)') '------ > prior to the continue command       ', achar(7) 
            prompt_stop = prompt
            lblock_dbg = .true. 
            if(lmakro) lmakro_dbg = .true.
            lblock = .false. 
            line = '#' 
            length = 1 
!                                                                       
!     ----Continuous loop until debug mode is switched off              
!                                                                       
            DO WHILE (lblock_dbg) 
               CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt)
               IF (line (1:1) .eq.'@') THEN 
                  line(1:length-1) = line(2:length)
                  length = length - 1
                  CALL file_kdo (line, length) 
               ELSE 
                  CALL p_mache_kdo (line, lend, length) 
               ENDIF 
               IF (ier_num.ne.0.AND.ier_sta.ne.ER_S_LIVE) THEN 
                  EXIT main               ! Error finish loop
               ELSEIF (ier_num.ne.0.AND.ier_sta.eq.ER_S_LIVE) THEN 
                  CALL errlist 
               ENDIF 
            ENDDO 
            IF (.NOT.lblock) THEN 
               EXIT main               ! Error finish loop
            ENDIF 
            line = '#' 
            length = 1 
         ENDIF 
         IF (line (1:1) .eq.'@') THEN 
            line(1:length-1) = line(2:length)
            length = length - 1
            CALL file_kdo (line, length) 
         ELSE 
            CALL p_mache_kdo (line, lend, length) 
         ENDIF 
         IF (ier_num.ne.0.AND.ier_sta.ne.ER_S_LIVE) THEN 
            EXIT main               ! Error finish loop
         ELSEIF (ier_num.ne.0.AND.ier_sta.eq.ER_S_LIVE) THEN 
            CALL errlist 
         ENDIF 
      ENDIF 
   ENDIF 
   if(.not.lblock) then
      if(ier_num==0) then
         write(output_io,'(a)') ' Block terminated by user'
      endif
      exit main
   endif
ENDDO main
!
!
IF (ier_num.ne.0) then 
   WRITE (ier_msg (4), 3000) 
   WRITE (ier_msg (5), 3100) line (1:41) 
   WRITE (ier_msg (6), 3200) 'Error in block line: ',ilevel(level),' Nested', level+1
   CALL errlist
   IF(do_macro(ilevel(level),level)/=' ') THEN
      WRITE (output_io,'(a,'' ***MAC *** '',a,a)') TRIM(color_err), &
      do_macro(ilevel(level),level)(1:LEN_TRIM(do_macro(ilevel(level),level))), &
      TRIM(color_fg)
   ENDIF
!  ier_num = -5
!  ier_typ = ER_FORT
   IF(lmakro .AND. lmacro_close) THEN
      CALL macro_close(-1)
   ENDIF 
ENDIF 
!                                                                       
!                                                                       
!2000 FORMAT    (a1) 
 3000 FORMAT    ('Erroneous line in block structure') 
 3100 FORMAT    (a41) 
 3200 FORMAT    (a, i5, a, i3)
!
END SUBROUTINE do_execute_block
!
!*****7**************************************************************** 
!O!      INTEGER FUNCTION bindex (string, char) 
!O!!                                                                       
!O!!     Searches for 'char' in 'string' from the back                     
!O!!                                                                       
!O!      IMPLICIT none 
!O!!                                                                       
!O!      CHARACTER ( * ) string 
!O!      CHARACTER(1) char 
!O!!                                                                       
!O!      INTEGER il, ik 
!O!      INTEGER len_str 
!O!!                                                                       
!O!      ik = 0 
!O!      il = len_str (string) 
!O!!                                                                       
!O!      DO while (il.ge.1.and.ik.eq.0) 
!O!      IF (string (il:il) .eq.char) then 
!O!         ik = il 
!O!      ELSE 
!O!         il = il - 1 
!O!      ENDIF 
!O!      ENDDO 
!O!!                                                                       
!O!      bindex = ik 
!O!!                                                                       
!O!      END FUNCTION bindex                           
!*****7**************************************************************** 
END MODULE do_if_mod
