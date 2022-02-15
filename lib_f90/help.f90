MODULE lib_help
!
CONTAINS
!
!****7******************************************************************
!                                                                       
SUBROUTINE do_hel (ein, inlength) 
!-                                                                      
!     This sublevel emulates the HELP function under VMS and            
!     controls the online help of DISCUS.                               
!                                                                       
!     Version : 1.3                                                     
!     Date    : 25 Oct 96                                               
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)      
!               Th. Proffen (proffen@pa.msu.edu)                        
!                                                                       
!-                                                                      
USE doact_mod 
USE envir_mod 
USE errlist_mod 
USE class_macro_internal 
USE lib_errlist_func
USE lib_length
USE prompt_mod 
USE string_convert_mod
USE sup_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(IN) :: ein         ! input string
INTEGER         , INTENT(IN) :: inlength    ! input string length
!                                                                       
INTEGER, PARAMETER :: IHL = 77
INTEGER, PARAMETER :: MAXW = 10
!                                                                       
CHARACTER(LEN=80) :: line, zeile, dummy 
CHARACTER(LEN=40) :: prom 
CHARACTER(LEN=12), DIMENSION(MAXW) :: bef
CHARACTER(LEN=4)  :: befehl 
CHARACTER(LEN=3)  :: status 
INTEGER  :: i, ibef
INTEGER, DIMENSION(MAXW) :: il
INTEGER  :: ll, lbef, lp 
INTEGER  :: length
LOGICAL  :: lread, lende 
LOGICAL  :: found, fast, stay 
!                                                                       
!                                                                       
!------ set up some variables                                           
!                                                                       
line = ein 
length = inlength
ier_num = 0 
ier_typ = ER_NONE 
status = 'old' 
lread = .true. 
lende = .false. 
fast = .false. 
!                                                                       
!------ if length was given negative we do not want to enter help level 
!                                                                       
IF (length <  0) then 
   stay = .false. 
   length = - length 
ELSE 
   stay = .true. 
ENDIF 
!                                                                       
!------ if we are in a loop or macro do nothing !!                      
!                                                                       
IF (lmakro.OR.lblock) RETURN 
10 CONTINUE 
main_loop: DO
!                                                                       
!------ - print help file entry                                         
!                                                                       
   IF (stay) CALL do_status (1, bef, ibef, il, maxw) 
   CALL oeffne (ihl, hlpfile, status) 
   IF (ier_num /= 0) RETURN 
   CALL do_cap (line) 
   CALL do_trenn (line, bef, ibef, length, maxw) 
   CALL do_help_single (ihl, bef, ibef, maxw, found, fast) 
   CLOSE (ihl) 
!                                                                       
!------ - Entry found ?                                                 
!                                                                       
   IF (ier_num /= 0) RETURN 
!                                                                       
!------ - if there are no subcommands get up one level                  
!                                                                       
   IF (.NOT.found) THEN 
      IF (ibef >  1) THEN 
         ibef = ibef - 1 
         WRITE (line, 1000) (bef (i), i = 1, ibef) 
         length = len_str (line) 
      ENDIF 
   ENDIF 
!                                                                       
!------ - Want to enter help level ?                                    
!                                                                       
   IF (.NOT.stay) RETURN 
!                                                                       
!     - Set reading mode back to normal                                 
!                                                                       
   fast = .false. 
!                                                                       
!------ - prompt user for input in online help                          
!                                                                       
   15 CONTINUE 
   CALL do_status (2, bef, ibef, il, maxw) 
   prom = prompt (1:len_str (prompt) ) //'/help' 
   CALL get_cmd (zeile, ll, befehl, lbef, dummy, lp, prom) 
!                                                                       
!------ - we will go up one level '..'                                  
!                                                                       
   IF(zeile(1:2)  == '..') THEN 
      IF (ibef >  1) THEN 
         ibef = ibef - 1 
         WRITE(line, 1000) (bef (i), i = 1, ibef) 
         length = len_str (line) 
      ENDIF 
      GOTO 15 
!                                                                       
!------ - leave help level ' '                                          
!                                                                       
   ELSEIF(ll == 0) THEN 
      lende = .true. 
!                                                                       
!     - repeat list of sublevels availble from current level            
!                                                                       
   ELSEIF(zeile(1:1)  == '?') THEN 
      fast = .true. 
!                                                                       
!------ - give help for this subcommand 'word'                          
!                                                                       
   ELSE 
      WRITE (line, 1000) (bef (i), i = 1, ibef), zeile 
      length = len_str (line) 
   ENDIF 
!                                                                       
   IF (lende) EXIT main_loop
ENDDO main_loop
!                                                                       
 1000 FORMAT    (8(a12,1x)) 
END SUBROUTINE do_hel                         
!
!****7******************************************************************
!
SUBROUTINE do_status (iwhere, bef, ibef, il, MAXW) 
!-                                                                      
!     Prints out status line in online help                             
!-                                                                      
USE prompt_mod 
USE string_convert_mod
USE lib_length
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER                          , INTENT(IN) :: iwhere
INTEGER                          , INTENT(IN) :: MAXW 
CHARACTER(LEN=*), DIMENSION(MAXW), INTENT(IN) :: bef
INTEGER                          , INTENT(IN) :: ibef
INTEGER, DIMENSION(MAXW)         , INTENT(OUT) :: il
!                                                                       
CHARACTER(LEN=80) :: level 
CHARACTER(LEN=6) ::prom 
INTEGER :: i 
!                                                                       
!                                                                       
!------ Build status lines                                              
!                                                                       
IF(iwhere == 1) THEN 
   prom = prompt(1:LEN(prom))
   CALL do_cap (prom) 
   WRITE(output_io, 2000) prom, version 
ELSEIF(iwhere == 2) THEN 
   DO i = 1, ibef 
      il(i) = len_str(bef(i)) 
   ENDDO 
   WRITE(level, 1000) (bef(i)(1:il(i)), i = 1, ibef) 
   WRITE(output_io, 3000) level 
ENDIF 
!                                                                       
 1000 FORMAT    (8(a,1x)) 
 2000 FORMAT    (' ----------- > ',A6,' online help ',15x,              &
     &                  'Version : ',a10,/)                              
 3000 FORMAT    (' ----------- > current help level : ',a40) 
!
END SUBROUTINE do_status                      
!
!****7******************************************************************
!
SUBROUTINE do_help_single (ihl, bef, ibef, maxw, found, fast) 
!-                                                                      
!     This routine emulates the VMS help function itself                
!-                                                                      
USE envir_mod 
USE errlist_mod 
USE lib_errlist_func
USE lib_length
use precision_mod
USE prompt_mod 
USE string_convert_mod
USE str_comp_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: ihl
INTEGER, INTENT(IN) :: ibef
INTEGER, INTENT(IN) :: MAXW
CHARACTER(LEN=*), DIMENSION(MAXW), INTENT(IN) ::  bef
LOGICAL, INTENT(OUT) :: found
LOGICAL, INTENT(IN)  :: fast
!                                                                       
INTEGER, PARAMETER :: nviele = 201
!                                                                       
CHARACTER(LEN=80) :: line 
CHARACTER(LEN=80) :: zeile 
CHARACTER(LEN=12), DIMENSION(nviele) :: viele
INTEGER :: lev, nl, nbef, ianz, laenge, ll 
REAL(kind=PREC_DP) :: werte 
!                                                                       
!                                                                       
!------ Loop over all helpfile levels until entry is found              
!                                                                       
main: DO lev = 1, ibef 
    5 CONTINUE 
      found = .false. 
      READ (ihl, 1000, end = 10) zeile 
      ll = len_str (zeile) 
      line = zeile (1:ll) 
!                                                                       
      IF (line (1:1) .eq.'!') then 
         line (1:1) = ' ' 
         line (2:2) = ' ' 
      ENDIF 
!                                                                       
      CALL do_cap (line) 
      CALL hole_zahl (line (1:3), ianz, werte) 
      IF (ianz.eq.1) then 
         nl = int (werte) 
         IF (nl.eq.lev) then 
            laenge = len_str (bef (lev) ) 
            found = str_comp (line (4:4 + laenge-1), bef (lev), laenge, &
            laenge, laenge)                                             
            IF (.not.found) goto 5 
         ELSEIF (nl.gt.lev.or.nl.eq.0) then 
            GOTO 5 
         ELSEIF (nl.lt.lev) then 
            GOTO 10 
         ENDIF 
      ELSE 
         GOTO 5 
      ENDIF 
ENDDO main
!                                                                       
!------ If entry is found, print text and get possible further          
!------ subentries                                                      
!                                                                       
   10 CONTINUE 
      nbef = 0 
      IF (found) then 
         IF (.not.fast) then 
            CALL lese_text (ihl) 
         ELSE 
            READ (ihl, 1000, end = 10) zeile 
            ll = len_str (zeile) 
            line = zeile (1:ll) 
!                                                                       
            IF (line (1:1) .eq.'!') then 
               line (1:1) = ' ' 
               line (2:2) = ' ' 
            ENDIF 
!                                                                       
         ENDIF 
         CALL weitere (ihl, ibef + 1, viele, nviele, nbef) 
         CALL schreib_viele (viele, nviele, nbef) 
         WRITE (output_io, * ) ' ' 
      ELSE 
         ier_typ = ER_IO 
         ier_num = - 5 
         CALL errlist 
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
!                                                                       
!     Mark whether subentries were found                                
!                                                                       
      found = nbef.gt.0 
!                                                                       
 1000 FORMAT  (a) 
END SUBROUTINE do_help_single                 
!
!*****7*****************************************************************
!
SUBROUTINE hole_zahl (zeile, ianz, werte) 
!+                                                                      
!     Gets numbers from command line and number of values of input.     
!-                                                                      
use precision_mod
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(IN)  :: zeile 
INTEGER         , INTENT(OUT) :: ianz 
REAL(kind=PREC_DP), INTENT(OUT) :: werte 
!                                                                       
werte = 0.0 
ianz = 0 
READ(zeile, *, end = 98, err = 98) werte 
ianz = 1 
RETURN 
   98 CONTINUE 
ianz = 0 
!
END SUBROUTINE hole_zahl                      
!
!*****7*****************************************************************
!
SUBROUTINE do_trenn (ein, bef, ibef, length, MAXW) 
!+                                                                      
!     Separates the character string into single words                  
!-                                                                      
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN)  :: MAXW
INTEGER, INTENT(OUT) :: ibef 
CHARACTER(LEN=*)                 , INTENT(IN)  :: ein
CHARACTER(LEN=*), DIMENSION(MAXW), INTENT(OUT) :: bef
INTEGER, INTENT(IN) :: length 
!
INTEGER :: ib, i, ianf, iend
LOGICAL :: lchar 
!                                                                       
ib = 1 
main_loop: DO i = 1, 10 
   lchar = .true. 
   CALL locate(lchar, ein, ib, ianf, length) 
   IF (ianf >  length) EXIT main_loop
   lchar = .false. 
   ib = ianf 
   CALL locate(lchar, ein, ib, iend, length) 
   bef(i) = ein(ianf:iend-1) 
   ib = iend 
ENDDO  main_loop
!
ibef = i - 1 
!
END SUBROUTINE do_trenn                       
!
!*****7*****************************************************************
!
SUBROUTINE locate (lchar, ein, ib, ianf, length) 
!                                                                       
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: lchar 
CHARACTER(LEN=*), INTENT(IN) :: ein 
INTEGER, INTENT(IN)  :: ib 
INTEGER, INTENT(OUT) :: ianf
INTEGER, INTENT(IN)  :: length
!
INTEGER :: i
!                                                                       
ianf = length + 1 
!                                                                       
main_loop: DO i = ib, length 
   IF (lchar.AND.ein (i:i)  /= ' '.AND.ein (i:i)  /= ',') THEN
      ianf = i 
      EXIT main_loop
   ELSEIF(.NOT.lchar.AND. (ein(i:i)  == ' '.OR.ein(i:i)  == ',') ) THEN
      ianf = i 
      EXIT main_loop
   ENDIF 
ENDDO main_loop
!
END SUBROUTINE locate                         
!
!*****7*****************************************************************
!
SUBROUTINE lese_text (ihl) 
!                                                                       
USE envir_mod 
USE lib_length
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ihl
!
CHARACTER(LEN=PREC_STRING) :: clines
CHARACTER(LEN=80)  ::line 
CHARACTER(LEN=1)  ::cdummy 
INTEGER :: ilc, ll 
!                                                                       
!                                                                       
main_loop: DO
   CALL GET_ENVIRONMENT_VARIABLE('LINES', clines)
   IF(clines==' ') THEN
      lines = 20
   ELSE
      READ(clines,*) lines
   ENDIF
   reading: DO ilc = 1, lines - 6 
      READ(ihl, 1000, end = 9010) line 
      IF(line(1:2) == '!%' ) CYCLE reading
!                                                                       
      IF(line (1:1) .eq.'!') THEN 
         line (1:1) = ' ' 
         line (2:2) = ' ' 
      ENDIF 
!                                                                       
      IF(line (1:1)  /= ' '.AND.line (1:1)  /= achar (13) ) THEN 
         EXIT main_loop
      ELSE 
         ll = len_str (line) 
         IF(ll >  0) THEN 
            WRITE(output_io, 1000) line (1:ll) 
         ELSE 
            WRITE(output_io, * ) 
         ENDIF 
      ENDIF 
   ENDDO  reading
   WRITE(output_io, 2000) 
   IF(output_io == 6) read ( *, 1000) cdummy 
ENDDO main_loop
!                                                                       
 9010 CONTINUE 
 1000 FORMAT  (a) 
 2000 FORMAT    (/,' ----------- > Press <RETURN> for next page ') 
!
END SUBROUTINE lese_text                      
!
!*****7*****************************************************************
!
SUBROUTINE weitere (ihl, nl, viele, nviele, nbef) 
!                                                                       
USE lib_length
use precision_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: ihl
INTEGER, INTENT(IN) :: nviele
CHARACTER(LEN=*), DIMENSION(NVIELE), INTENT(OUT) :: viele
INTEGER, INTENT(IN) :: nl
INTEGER, INTENT(OUT) :: nbef
!
CHARACTER(LEN=80) :: line 
CHARACTER(LEN=80) :: zeile 
INTEGER :: ianz 
INTEGER :: ll 
REAL(kind=PREC_DP) :: werte 
!                                                                       
!                                                                       
nbef = 0 
BACKSPACE (ihl) 
!
main_LOOP: DO
   READ(ihl, 1000, END = 9010) zeile 
   IF(zeile(1:2) == '!%' ) CYCLE main_LOOP
   ll = len_str(zeile) 
   line = zeile(1:ll) 
   inner_loop: DO WHILE(line(1:1)  == ' '.or.line(1:1)  == ACHAR(13)) 
      READ(ihl, 1000, END = 9010) zeile 
      IF(zeile(1:2) == '!%' ) CYCLE inner_LOOP
      ll = len_str(zeile) 
      line = zeile(1:ll) 
   ENDDO inner_loop
   CALL hole_zahl(line, ianz, werte) 
   IF(ianz == 1) then 
      IF(INT(werte)  == nl) then 
         nbef = nbef + 1 
         viele (nbef) = line(4:15) 
      ELSEIF(INT(werte)  <  nl) then 
         RETURN 
      ENDIF 
   ENDIF 
ENDDO main_loop
!
 9010 CONTINUE 
 1000 FORMAT  (a) 
!
END SUBROUTINE weitere                        
!
!*****7*****************************************************************
!
SUBROUTINE schreib_viele (viele, nviele, nbef) 
!                                                                       
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: nviele 
INTEGER, INTENT(IN) :: nbef 
CHARACTER(LEN=*), DIMENSION(nviele), INTENT(IN) :: viele
!
INTEGER :: i 
!                                                                       
WRITE(output_io, 1000) (viele (i), i = 1, nbef) 
!                                                                       
 1000 FORMAT  (5(3x,a12)) 
!
END SUBROUTINE schreib_viele                  
!
!*****7*****************************************************************
!
SUBROUTINE do_manual(line, length)
!
USE envir_mod 
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE prompt_mod 
USE take_param_mod
USE string_convert_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER  :: N_VIEWER = 5
INTEGER, PARAMETER  :: W_VIEWER = 5
INTEGER, PARAMETER  :: MAXW     = 2
INTEGER, PARAMETER  :: MAX_MAN  = 6
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: string, zeile
CHARACTER(LEN=    PREC_STRING           ) :: command, message
CHARACTER(LEN=128 ), DIMENSION(N_VIEWER)  :: pdf_viewer  ! List of possible viewers
CHARACTER(LEN=128 ), DIMENSION(W_VIEWER)  :: win_test_v  ! List of possible viewers
CHARACTER(LEN=128 ), DIMENSION(W_VIEWER)  :: win_viewer  ! List of possible viewers
CHARACTER(LEN= 7  ), DIMENSION(MAX_MAN )  :: c_manual    ! list of possibel sections
CHARACTER(LEN= 7  )                       :: manual      ! actual section
CHARACTER(LEN=PREC_STRING)                :: viewer      ! actual viewer
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara   !parameters
INTEGER            , DIMENSION(MAXW) :: lpara   !parameters
!
INTEGER, PARAMETER :: NOPTIONAL = 2
CHARACTER(LEN=    PREC_STRING           ), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Length opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Length opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER             :: ianz
INTEGER             :: laenge, i
INTEGER             :: ierror
LOGICAL             :: lexist
!
!
DATA oname  / 'section','viewer' /
DATA loname /  7       , 6 /
!
DATA pdf_viewer / 'qpdfview  ', 'evince    ', 'xpdf      ',   &
                  'okular    ', 'acroread  '/ 
DATA win_test_v /                                             &
'/cygdrive/c/Program Files (x86)/Foxit Software/Foxit Reader/FoxitReader.exe   ', &
'/cygdrive/c/Program Files (x86)/STDU Viewer/STDUViewerApp.exe   ', &
'/cygdrive/c/Program Files/SumatraPDF/SumatraPDF.exe      ', &
'/cygdrive/c/Program Files/Internet Explorer/iexplore.exe',  &
'/cygdrive/c/Program Files/Mozilla Firefox/firefox.exe   '   &
/
DATA win_viewer /                                             &
'/cygdrive/c/Program\ Files\ \(x86\)/Foxit\ Software/Foxit\ Reader/FoxitReader.exe', &
'/cygdrive/c/Program\ Files\ \(x86\)/STDU\ Viewer/STDUViewerApp.exe', &
'/cygdrive/c/Program\ Files/SumatraPDF/SumatraPDF.exe      ', &
'/cygdrive/c/Program\ Files/Internet\ Explorer/iexplore.exe', &
'/cygdrive/c/Program\ Files/Mozilla\ Firefox/firefox.exe   '  &
/
DATA c_manual / 'suite'  ,'discus' , 'diffev' ,'kuplot',      &
                'package','refine'  /
!
opara (1) = pname         ! Default to section name
lopara(1) = LEN_TRIM(pname)
IF(operating(1:7)=='Windows' .OR. operating(1:6)=='cygwin') THEN 
   opara (2) = 'foxit'
   lopara(2) =   5
ELSE
   opara (2) = pdf_viewer(1) ! Always provide fresh default values
   lopara(2) =   8
ENDIF
owerte    =  (/  0.0    , 0.0 /)
!
CALL get_params (line, ianz, cpara, lpara, maxw, length)
IF (ier_num.ne.0) THEN
   RETURN
ENDIF
!
!  Did user provide a section / viewer ?
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
manual  = opara(1)(1:MIN(LEN(manual),lopara(1)))    ! defaults to section name
ier_num = -6
ier_typ = ER_COMM
man: DO i=1,MAX_MAN
   IF(str_comp (opara(1), c_manual(i), 3, lopara(1), LEN_TRIM(c_manual(i)))) THEN
      manual  = opara(1)(1:MIN(LEN(manual),lopara(1)))
      ier_num = 0
      ier_typ = ER_NONE
      EXIT man
   ENDIF
ENDDO man
IF(ier_num /=0) THEN
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Selected package is not present, choose:'
   ier_msg(2) = 'suite, discus, diffev, kuplot, package,'
   ier_msg(3) = 'mixscat'
   RETURN
ENDIF
!
IF(operating(1:7)=='Windows' .OR. operating(1:6)=='cygwin') THEN 
   ierror = -6
   win_opt: DO i=1,W_VIEWER
      string = (win_viewer(i))
      zeile  = (opara(2))
      CALL do_cap(string)
      CALL do_cap(zeile)
      IF(INDEX(string,zeile(2:4))>0) THEN
!     IF(INDEX(win_viewer(i),opara(2)(2:4))>0) THEN
         string = win_test_v(i)
         INQUIRE(FILE=string,EXIST=lexist)
         IF(lexist) THEN
            viewer = win_viewer(i)
            ierror = 0
            EXIT win_opt
         ENDIF
      ENDIF
   ENDDO win_opt
!
   IF(ierror == -6) THEN
      win_search: DO i=1,W_VIEWER
         string = win_test_v(i)
         INQUIRE(FILE=string,EXIST=lexist)
         IF(lexist) THEN
            viewer = win_viewer(i)
            ierror = 0
            EXIT win_search
         ENDIF
      ENDDO win_search
   ENDIF
!
   i=LEN_TRIM(man_dir)
   IF(operating(1:7)=='Windows') THEN
      IF(man_dir(i:i) /='\') THEN
         man_dir(i+1:i+1) = '\'
      ENDIF
   ELSEIF(operating(1:8)=='cygwin64' .AND. &
          .NOT.( man_dir(1:2)=='c:' .OR. man_dir(1:2)=='C:' )) THEN
      man_dir = 'c:/cygwin64/'//man_dir(1:i)
   ENDIF
!
   IF(ierror==0) THEN     ! Finally, everything is fine let's do it
      command = viewer(1:LEN_TRIM(viewer))//' '// &
                '"'// &
                man_dir(1:LEN_TRIM(man_dir))//manual(1:LEN_TRIM(manual))//'_man.pdf"  &'
      laenge=LEN_TRIM(command)
      CALL EXECUTE_COMMAND_LINE (command(1:laenge),  CMDSTAT=ierror, CMDMSG=message)
   ENDIF
ELSEIF(operating(1:6)=='darwin') THEN 
!
! MAC OS use 'open' command
!
   viewer = 'open'
   command = viewer(1:LEN_TRIM(viewer))//' '// &
             man_dir(1:LEN_TRIM(man_dir))//manual(1:LEN_TRIM(manual))//'_man.pdf &'
   laenge=LEN_TRIM(command)
   CALL EXECUTE_COMMAND_LINE (command(1:laenge),  CMDSTAT=ierror, CMDMSG=message)
ELSE
!
! LINUX Choose viewer
!
   command = 'which '//opara(2)(1:lopara(2))   ! Try opional parameter first
   CALL EXECUTE_COMMAND_LINE (command(1:LEN_TRIM(command)),  CMDSTAT=ierror, CMDMSG=message)
!
   IF(ierror == 0 ) THEN  ! Default / User option did not work try list
      viewer = opara(2)(1:lopara(2))
   ELSE
      search: DO i=1, N_VIEWER
         command = 'which '//pdf_viewer(i)(1:LEN_TRIM(pdf_viewer(i)))
         CALL EXECUTE_COMMAND_LINE (command(1:LEN_TRIM(command)),  CMDSTAT=ierror, CMDMSG=message)
         IF(ierror ==0)  THEN
            viewer = pdf_viewer(i)
            EXIT search
         ENDIF
      ENDDO search
   ENDIF
!
   IF(ierror==0) THEN     ! Finally, everything is fine let's do it
      command = viewer(1:LEN_TRIM(viewer))//' '// &
                man_dir(1:LEN_TRIM(man_dir))//manual(1:LEN_TRIM(manual))//'_man.pdf &'
      laenge=LEN_TRIM(command)
      CALL EXECUTE_COMMAND_LINE (command(1:laenge),  CMDSTAT=ierror, CMDMSG=message)
   ENDIF
ENDIF
IF(ierror/=0) THEN     ! Error ??
   WRITE(output_io,1000) man_dir
ENDIF
!
1000 FORMAT(' Did not find a PDF viewer. '/ &
' You will find the Manual in the folder:'/' ',a)
!
END SUBROUTINE do_manual
!
!
END MODULE lib_help
