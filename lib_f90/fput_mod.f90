MODULE fput_mod
!
CONTAINS
!*****7***********************************************************      
      SUBROUTINE do_fopen (zeile, lp) 
!                                                                       
USE ber_params_mod
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE macro_mod 
USE precision_mod
      USE prompt_mod 
USE sys_compiler
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) ::  zeile 
      INTEGER          , INTENT(INOUT) :: lp 
!
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(LEN=PREC_STRING) :: message
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
      INTEGER lpara (maxw)
      INTEGER :: ios
      INTEGER ianz, ianzz 
      INTEGER ii 
      LOGICAL lappend 
      LOGICAL one_open 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      ianzz = 1 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (ianz.ge.1) THEN 
         CALL ber_params (ianzz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) RETURN 
         IF (ianz.eq.1) THEN 
            ianz = 0 
         ELSE 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
         ENDIF 
         IF (ier_num.ne.0) RETURN 
         ii = nint (werte (1) ) 
         IF (ii.lt.0.or.MAC_MAX_IO.lt.ii) THEN 
            ier_num = - 13 
            ier_typ = ER_IO 
            RETURN 
         ENDIF 
!                                                                       
         IF (ianz.ge.1) THEN 
            IF (io_open (ii) ) THEN 
               ier_num = - 10 
               ier_typ = ER_IO 
               RETURN 
            ENDIF 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ianz.eq.2) THEN 
               lappend = str_comp (cpara (2) , 'append', 1, lpara (2) , &
               6)                                                       
            ELSE 
               lappend = .false. 
            ENDIF 
            IF (lappend) THEN 
               CALL oeffne_append (io_unit (ii) , cpara (1) , 'unknown')
               IF (ier_num.ne.0) RETURN 
!DBG            open (unit=io_unit(ii),file=cpara(1),status='unknown',  
!DBG     &                       access='append',err=999)               
!DBG95     &                       access='sequential',err=999)         
               WRITE (output_io, 1000) cpara (1) (1:lpara (1) ) 
            ELSE 
               OPEN (unit = io_unit (ii) , file = cpara (1) , &
                     status = 'unknown', IOSTAT=ios,IOMSG=message)                                    
               IF(ios/=0) THEN
                  ier_num = -2
                  ier_typ = ER_IO
                  ier_msg(1) ='Could not open the NULL file'
                  ier_msg(2) = message(1:MAX(1,MIN(80,LEN_TRIM(message))))
                  RETURN
               ENDIF
               WRITE (output_io, 1100) cpara (1) (1:lpara (1) ) 
            ENDIF 
            io_open (ii) = .true. 
            io_file (ii) = cpara (1) (1:lpara (1) ) 
         ELSE 
            IF (io_open (ii) ) THEN 
               WRITE (output_io, 2000) ii, io_file (ii) 
            ELSE 
               WRITE (output_io, 2010) 
            ENDIF 
         ENDIF 
      ELSE 
         one_open = .false. 
         DO ii = 1, MAC_MAX_IO 
         IF (io_open (ii) ) THEN 
            WRITE (output_io, 2000) ii, io_file (ii) 
            one_open = .true. 
         ENDIF 
         ENDDO 
         IF (.not.one_open) THEN 
            WRITE (output_io, 2010) 
         ENDIF 
      ENDIF 
      RETURN 
!                                                                       
! 999 CONTINUE 
      ier_num = - 2 
      ier_typ = ER_IO 
!                                                                       
 1000 FORMAT     (' ------ > File ',a,' opened (appending) ...') 
 1100 FORMAT     (' ------ > File ',a,' opened (overwriting) ...') 
 2000 FORMAT     (' ------ > IO stream Nr. ',i2,' open, file : ',a) 
 2010 FORMAT     (' ------ > No IO stream open') 
      END SUBROUTINE do_fopen                       
!*****7***********************************************************      
      SUBROUTINE do_fclose (zeile, lp) 
!                                                                       
USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE macro_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 25) 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) ::  zeile 
      INTEGER          , INTENT(INOUT) :: lp 
!
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw)
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
      INTEGER lpara (maxw)
      INTEGER ianz, ii 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (ianz.eq.1) THEN 
         IF (str_comp (cpara (1) , 'all', 2, lpara (1) , 3) ) THEN 
            DO ii = 1, MAC_MAX_IO 
            IF (io_open (ii) ) THEN 
               CLOSE (io_unit (ii) ) 
               io_open (ii) = .false. 
            ENDIF 
            ENDDO 
         ELSE 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) RETURN 
!                                                                       
            ii = nint (werte (1) ) 
            IF (ii.lt.0.or.MAC_MAX_IO.lt.ii) THEN 
               ier_num = - 13 
               ier_typ = ER_IO 
               RETURN 
            ENDIF 
            IF (io_open (ii) ) THEN 
               CLOSE (io_unit (ii) ) 
               io_open (ii) = .false. 
            ELSE 
               ier_num = - 11 
               ier_typ = ER_IO 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      END SUBROUTINE do_fclose                      
!*****7***********************************************************      
      SUBROUTINE do_fend (zeile, lp) 
!                                                                       
USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE macro_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, PARAMETER :: maxw = 10 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) ::  zeile 
      INTEGER          , INTENT(INOUT) :: lp 
!
      CHARACTER(LEN=PREC_STRING), DIMENSION(maxw) :: cpara
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
      INTEGER lpara (maxw)
      INTEGER ianz, ianzz 
      INTEGER ii 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      ianzz = 1 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (ianz.ge.1) THEN 
         CALL ber_params (ianzz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) RETURN 
         IF (ianz.eq.1) THEN 
            ianz = 0 
         ELSE 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
         ENDIF 
         IF (ier_num.ne.0) RETURN 
         ii = nint (werte (1) ) 
         IF (ii.lt.0.or.MAC_MAX_IO.lt.ii) THEN 
            ier_num = - 13 
            ier_typ = ER_IO 
            RETURN 
         ENDIF 
!                                                                       
         IF (ianz.ge.1) THEN 
            IF (str_comp (cpara (1) , 'error', 2, lpara (1) , 5) ) THEN
               io_eof (ii) = .false. 
            ELSEIF (str_comp(cpara(1), 'continue', 2, lpara(1), 8)) THEN
               io_eof (ii) = .true. 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fend                        
!*****7***********************************************************      
!
SUBROUTINE do_fget (zeile, lp) 
!                                                                       
USE ber_params_mod
USE calc_expr_mod
USE charact_mod
USE errlist_mod 
USE get_params_mod
USE macro_mod 
USE param_mod 
USE precision_mod
USE take_param_mod
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER :: MAXW = 25
!                                                                       
CHARACTER (LEN=*), INTENT(IN) ::  zeile 
INTEGER          , INTENT(INOUT) :: lp 
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw), line, cstr , bstr
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: fpara
CHARACTER(LEN=PREC_LSTRING) :: string 
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
INTEGER lpara (maxw), lstr
INTEGER ianz, i, igl, ii, ianzz 
INTEGER ia, ie, itab , ll
INTEGER :: k1, k2, length, idot
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER len_str 
!
DATA oname  / 'form'   /
DATA loname /  4       /
opara  =  (/ '*'       /)   ! Always provide fresh default values
lopara =  (/  1        /)
owerte =  (/  0.0      /)
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) RETURN 
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(opara(1)/='*') THEN
   line(1:lopara(1)-2) = opara(1)(2:lopara(1)-1)
   lp = lopara(1) - 2
   CALL get_params (line, ianzz, fpara, lpara, MAXW, lp)
   DO i=1,ianz
      lpara(i) = LEN_TRIM(cpara(i))   ! Restore parameter length
   ENDDO
ENDIF
!
ianzz = 1 
!                                                                       
CALL ber_params (ianzz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) RETURN 
IF (ianz.eq.1) THEN 
   ianz = 0 
ELSE 
   CALL del_params (1, ianz, cpara, lpara, maxw) 
ENDIF 
ii = nint (werte (1) ) 
IF (ii.lt.1.or.MAC_MAX_IO.lt.ii) THEN 
   ier_num = - 13 
   ier_typ = ER_IO 
   RETURN 
ENDIF 
!                                                                       
IF (.not.io_open (ii) ) THEN 
   ier_num = - 11 
   ier_typ = ER_IO 
   RETURN 
ENDIF 
!                                                                       
IF (ianz.gt.0) THEN 
   READ (io_unit (ii) , '(a)', err = 998, end = 999) string 
   DO while (string (1:1) .eq.'#') 
      READ (io_unit (ii) , '(a)', err = 998, end = 999) string 
   ENDDO 
   lstr = len_str (string) 
   itab = index (string, TAB) 
   DO while (itab.gt.0) 
      string (itab:itab) = ' ' 
      itab = index (string, TAB) 
   ENDDO 
   IF(opara(1)=='*') THEN
      ia = max (1, io_get_sub (ii, 1) ) 
      IF (io_get_sub (ii, 2) .eq. - 1) THEN 
         ie = lstr 
      ELSE 
         ie = min (lstr, io_get_sub (ii, 2) ) 
      ENDIF 
      loop_para: DO i = 1, ianz 
         line = cpara (i) (1:lpara (i) ) //'=' 
         lp = lpara (i) + 1 
         igl = index (line, '=') 
         cstr (1:1) = string (ia:ia) 
!                                                                       
!     ---Remove leading blanks or leading ","                           
!                                                                       
         DO while (cstr (1:1) .eq.' '.or.cstr (1:1) .eq.',') 
            ia = ia + 1 
            IF (ia.gt.ie) THEN 
               GOTO 998 
            ENDIF 
            cstr (1:1) = string (ia:ia) 
         ENDDO 
!                                                                       
!     ---Copy the characters into the line, until a blank or ","        
!                                                                       
         ll = 0
         bstr = ' '
         DO while (.not. (cstr (1:1) .eq.' '.or.cstr (1:1) .eq.',') ) 
            ll = ll + 1
            lp = lp + 1 
            line (lp:lp) = cstr (1:1) 
            bstr (ll:ll) = cstr (1:1)   ! make copy for character evaluation
            ia = ia + 1 
            IF (ia.gt.ie) THEN 
               GOTO 997 
            ENDIF 
            cstr (1:1) = string (ia:ia) 
         ENDDO 
997      CONTINUE 
         CALL do_math (line, igl, lp) 
         IF(ier_num /= 0) THEN    ! Error in math assum string 
            line = cpara(i)(1:lpara(i)) //' = '''// bstr(1:LEN_TRIM(bstr))//''''  
            igl = lpara(i) + 2
            lp  = LEN_TRIM(line)
            CALL do_math (line, igl, lp)
         ENDIF
      ENDDO loop_para
   ELSE
      k1 = 0
      k2 = 0
      DO i = 1, ianz
         k1 = k2 + 1
         IF(fpara(i)(1:1) == 'i') THEN
            READ(fpara(i)(2:LEN_TRIM(fpara(i))),*) length
            k2 = k2 + length
            line = cpara(i)(1:lpara(i)) //' = '//string(k1:k2)
            igl = lpara(i) + 2
            lp  = LEN_TRIM(line)
         ELSEIF(fpara(i)(1:1) == 'f') THEN
            idot = INDEX(fpara(i),'.')
            IF(idot==0) THEN
               READ(fpara(i)(2:LEN_TRIM(fpara(i))),*) length
            ELSE
               READ(fpara(i)(2:idot-1            ),*) length
            ENDIF
            k2 = k2 + length
            line = cpara(i)(1:lpara(i)) //' = '//string(k1:k2)
            igl = lpara(i) + 2
            lp  = LEN_TRIM(line)
         ELSEIF(fpara(i)(1:1) == 'a') THEN
            READ(fpara(i)(2:LEN_TRIM(fpara(i))),*) length
            k2 = k2 + length
            line = cpara(i)(1:lpara(i)) //' = '''//string(k1:k2)//''''
            igl = lpara(i) + 2
            lp  = LEN_TRIM(line)
         ENDIF
         CALL do_math (line, igl, lp)
      ENDDO
   ENDIF
   res_para (0) = 0 
ELSE 
   READ (io_unit (ii), *, err = 998, end = 999) 
   res_para (0) = 0 
ENDIF 
!                                                                       
RETURN 
!                                                                       
998 CONTINUE 
ier_num = - 3 
ier_typ = ER_IO 
RETURN 
!                                                                       
999 CONTINUE 
IF (io_eof (ii) ) THEN 
   res_para (0) = - 1 
ELSE 
   ier_num = - 6 
   ier_typ = ER_IO 
ENDIF 
!                                                                       
END SUBROUTINE do_fget                        
!*****7***********************************************************      
      SUBROUTINE do_fformat (zeile, lp) 
!                                                                       
USE ber_params_mod
      USE debug_mod 
      USE errlist_mod 
      USE get_params_mod
      USE macro_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER, DIMENSION(MAXW) :: lpara
INTEGER                  :: ii, iii, ianz 
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
!                                                                       
      INTEGER len_str 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (ianz.eq.0) THEN 
         WRITE (output_io, 1000) 
         DO ii = 1, MAC_MAX_FORM, 6 
            WRITE (output_io, 1100) (ii + iii, iii = 0, MIN (5, MAC_MAX_FORM - ii) )
            WRITE (output_io, 1200) (io_out_format (ii + iii)          &
                  (2:MIN(len_str(io_out_format(ii+iii))-1, 8)), iii=0, &
                  MIN (5, MAC_MAX_FORM - ii) )                                           
         ENDDO 
      ELSEIF (ianz.eq.2) THEN 
         CALL ber_params (1, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) RETURN 
         ii = NINT (werte (1) ) 
         IF (ii.gt.0.and.ii.le.MAC_MAX_FORM) THEN 
            io_out_format (ii) = '('//cpara (2) (1:lpara (2) ) //')' 
         ELSE 
            ier_num = - 28
            ier_typ = ER_IO 
            WRITE(ier_msg(3),9000) MAC_MAX_FORM
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Current format setting for file output :',/) 
 1100 FORMAT     ('   Column : ',6(i8,2x)) 
 1200 FORMAT     ('   Format : ',6(a8,2x),/) 
 9000 FORMAT ( 'The first parameter must be in range 0 to',i3)

      END SUBROUTINE do_fformat                     
!*****7***********************************************************      
      SUBROUTINE do_fgetsub (zeile, lp) 
!                                                                       
USE ber_params_mod
      USE debug_mod 
      USE errlist_mod 
      USE get_params_mod
      USE macro_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER, DIMENSION(MAXW) :: lpara
INTEGER                  :: ii, iii, ianz 
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      iii = 1 
      CALL ber_params (iii, cpara, lpara, werte, maxw) 
      IF (ianz.eq.1) THEN 
         ianz = 0 
      ELSE 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
      ENDIF 
      ii = nint (werte (1) ) 
      IF (ii.lt.0.or.MAC_MAX_IO.lt.ii) THEN 
         ier_num = - 13 
         ier_typ = ER_IO 
         RETURN 
      ENDIF 
!                                                                       
      IF (ianz.eq.0) THEN 
         io_get_sub (ii, 1) = 1 
         io_get_sub (ii, 2) = - 1 
      ELSEIF (ianz.eq.2) THEN 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) RETURN 
         io_get_sub (ii, 1) = nint (werte (1) ) 
         io_get_sub (ii, 2) = nint (werte (2) ) 
         IF (io_get_sub (ii, 2) .ne. - 1.and.io_get_sub (ii, 1)         &
         .gt.io_get_sub (ii, 2) ) THEN                                  
            ier_num = - 26 
            ier_typ = ER_IO 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!DBG                                                                    
!DBG      write(*,*) 'cpara      ',cpara                                
!DBG      write(*,*) 'werte      ',werte                                
!DBG      write(*,*) 'io_get_sub ',io_get_sub                           
!                                                                       
      END SUBROUTINE do_fgetsub                     
!*****7***********************************************************      
      SUBROUTINE do_fput (zeile, lp) 
!                                                                       
USE ber_params_mod
USE berechne_mod
USE blanks_mod
USE build_name_mod
USE calc_expr_mod
USE debug_mod 
USE errlist_mod 
USE get_params_mod
USE macro_mod 
USE precision_mod
USE take_param_mod
!
IMPLICIT none 
!
INTEGER, PARAMETER :: maxw = MAC_MAX_FORM 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: fpara
CHARACTER(LEN=PREC_STRING)                  :: cstr, line 
CHARACTER(LEN=1)                     :: quote 
REAL                     :: wert 
INTEGER, DIMENSION(MAXW) :: lpara
INTEGER                  :: ianz, lstr, i, ie, ianzz, ii 
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER len_str 
!
DATA oname  / 'form'   /
DATA loname /  4       /
opara  =  (/ '*'       /)   ! Always provide fresh default values
lopara =  (/  1        /)
owerte =  (/  0.0      /)
!                                                                       
quote = achar (39) 
ianzz = 1 
!                                                                       
lp = -lp
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) RETURN 
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(opara(1)/='*') THEN
   line = ' '
   line(1:lopara(1)-2) = opara(1)(2:lopara(1)-1)
   lp = lopara(1) - 2
   CALL get_params (line, ianzz, fpara, lpara, MAXW, lp)
   DO i=1,ianz
      lpara(i) = LEN_TRIM(cpara(i))   ! Restore parameter length
   ENDDO
   IF(ianz-1<= MAC_MAX_FORM) THEN
      DO i=1,ianz-1
         io_out_format (i) = '('//fpara(i)(1:MIN(20,LEN_TRIM(fpara(i))))//')'
      ENDDO
   ENDIF
ENDIF
!                                                                       
CALL ber_params (ianzz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) RETURN 
IF (ianz.eq.1) THEN 
   ianz = 0 
ELSE 
   CALL del_params (1, ianz, cpara, lpara, maxw) 
ENDIF 
ii = nint (werte (1) ) 
IF (ii.lt.0.or.MAC_MAX_IO.lt.ii) THEN 
   ier_num = - 13 
   ier_typ = ER_IO 
   RETURN 
ENDIF 
!                                                                       
IF (.not.io_open (ii) ) THEN 
   ier_num = - 11 
   ier_typ = ER_IO 
   RETURN 
ENDIF 
!                                                                       
IF (ianz.ge.1) THEN 
   ie = 0 
   CALL rem_leading_bl(cpara(1),lpara(1))
   IF (cpara (1) (1:1) .eq.'"') THEN 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      line = cpara (1) 
      ie = lpara (1) 
   ELSE 
      DO i = 1, ianz 
         IF (index (cpara (i) (1:lpara (i) ), quote) .ne.0) THEN 
            line(ie+1:ie+lpara (i) - 2) = cpara(i)(2:lpara(i) - 1)
            ie = ie+lpara (i) - 2 
         ELSE 
            cstr = '('//cpara (i) (1:lpara (i) ) //')' 
            lstr = lpara (i) + 2 
            CALL rem_bl (cstr, lstr) 
            wert = berechne (cstr, lstr) 
            IF (ier_num.ne.0) RETURN 
            IF (io_out_format (i) .eq.'(*)') THEN 
               WRITE (cstr, *, err = 999) wert 
            ELSE 
               IF (index(io_out_format (i) , 'i') .ne.0 .OR.  &
                   index(io_out_format (i) , 'I') .ne.0 .OR.  &
                   index(io_out_format (i) , 'z') .ne.0 .OR.  &
                   index(io_out_format (i) , 'Z') .ne.0)    THEN                  
                  WRITE(cstr, io_out_format(i), err = 999) NINT( wert)
               ELSE 
                  WRITE(cstr, io_out_format(i), err = 999) wert 
               ENDIF 
            ENDIF 
            lstr = len_str (cstr) 
            line (ie+1:ie+lstr) = cstr (1:lstr) 
            ie = ie+lstr 
         ENDIF 
      ENDDO 
   ENDIF 
   WRITE (io_unit (ii), 9999, err = 999) line (1:ie) 
   IF (dbg) WRITE ( *, 1000) line (1:ie) 
ELSE 
   WRITE (io_unit (ii), *, err = 999) 
ENDIF 
RETURN 
!                                                                       
  999 ier_num = - 12 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 1000 FORMAT    ( ' debug  > Written: ',a) 
 9999 FORMAT    (a) 
!
END SUBROUTINE do_fput                        
!
END MODULE fput_mod
