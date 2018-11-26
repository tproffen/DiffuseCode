MODULE do_variable_mod
!
CONTAINS
!*****7**************************************************************** 
SUBROUTINE ersetz_variable (line, length) 
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate user defined variable.                              
!     This is needed if the parameter is read.                          
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@fau.de)    
!+                                                                      
USE blanks_mod
USE constants_mod
USE errlist_mod 
USE precision_mod
USE string_extract_mod
USE variable_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(INOUT) :: length 
CHARACTER(LEN=1024) :: string 
CHARACTER(LEN=1024) :: substring 
CHARACTER(LEN=1024) :: zeile 
CHARACTER(LEN=1024) :: dummy 
!                                                                       
INTEGER :: i, ianf, iend, ll 
INTEGER :: linsert 
INTEGER :: laenge
INTEGER :: istart, istop
INTEGER :: s1,s2,s3
LOGICAL :: success = .FALSE.
!                                                                       
INTEGER len_str 
!                                                                       
string = line
istart = 1
istop  = length
s1 = 0
s2 = 0
s3 = 0
line = ' '
!
main: DO WHILE(s2<istop)     ! Loop over all non-quoted section of string
   CALL string_extract(string,istart, istop, substring, s1,s2,s3)
!
   linsert = 0
   laenge  = s2 - s1 + 1     ! Length of current string
   ll = laenge 
!                                                                       
!     Loop over all defined variables and search for coincidences       
!     of the variable name. Works since names like "cos" etc are        
!     forbidden, and the names have been sorted in alphabetically        
!     descending order.                                                 
!                                                                       
!                                                                       
!     --If an apostrophe is found, ignore the string                    
!                                                                       
!IF (max (INDEX (string, '''') , INDEX (string, '"') ) .gt.0) THEN 
!   RETURN 
!ENDIF 
   names:DO i = 1, var_num 
   success = .FALSE.
   ianf = INDEX (substring, var_name (i) (1:var_l (i) ) ) 
   DO while (ianf.ne.0) 
      IF(var_entry(i)>0) CYCLE names        ! This is a variable field
      zeile = ' ' 
      iend = ianf + var_l (i) - 1 
      IF (ianf.gt.1) zeile (1:ianf - 1) = substring (1:ianf - 1) 
      IF (var_type (i) .eq.      IS_REAL) THEN 
                                                                        
         WRITE (dummy (1:PREC_WIDTH) , PREC_F_REAL) var_val (i) 
         dummy (PREC_MANTIS+1:PREC_MANTIS+1) = 'e' 
         ll = PREC_WIDTH
         CALL rem_bl (dummy, ll) 
         zeile (ianf:ianf + ll - 1) = dummy (1:ll) 
         linsert = ll 
      ELSEIF (var_type (i) .eq.      IS_INTE) THEN 
         WRITE (dummy (1:PREC_WIDTH) , PREC_F_INTE) nint (var_val (i) ) 
         ll = PREC_WIDTH 
         CALL rem_bl (dummy, ll) 
         zeile (ianf:ianf + ll - 1) = dummy (1:ll) 
         linsert = ll 
      ELSEIF (var_type (i) .eq.      IS_CHAR) THEN 
         ll = len_str (var_char (i) ) 
         zeile (ianf:ianf) = '''' 
         zeile (ianf + 1:ianf + ll) = var_char (i) (1:ll) 
         zeile (ianf + ll + 1:ianf + ll + 1) = '''' 
         linsert = ll + 2 
      ENDIF 
      ll = laenge+linsert - (iend-ianf + 1) 
      IF (iend.lt.laenge) zeile (ianf + linsert:ll) = substring (iend+1:   &
      laenge)                                                           
      substring = zeile 
      istart = ianf+linsert                  ! Start further search after insert
      laenge = ll 
      success = .TRUE.
      ianf = INDEX (substring(istart:len_trim(substring)), var_name (i) (1:var_l (i) ) ) 
      IF(ianf>0) ianf = ianf + istart - 1    ! if found, correct the offset
   ENDDO 
   ENDDO names
!
   line = line(1:LEN_TRIM(line))//substring(1:LEN_TRIM(substring))  &
          //string(s2+1:s3-1)
   istart = s3
!
ENDDO main
length = LEN_TRIM(line)
!
END SUBROUTINE ersetz_variable                
!
!*****7**************************************************************** 
!
SUBROUTINE upd_variable (string, laenge, wert, dummy, length) 
!-                                                                      
!       updates the user defined variable by the value wert             
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)    
!+                                                                      
      USE constants_mod
      USE lib_f90_allocate_mod
      USE errlist_mod 
      USE param_mod
      USE variable_mod
      IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER(LEN=*), INTENT(IN) :: string 
INTEGER         , INTENT(IN) :: laenge 
REAL            , INTENT(IN) :: wert 
CHARACTER(LEN=*), INTENT(IN) :: dummy 
INTEGER         , INTENT(IN) :: length 
!                                                                       
      INTEGER i 
!                                                                       
      ier_num = - 24 
      ier_typ = ER_FORT 
!                                                                       
!     loop over all variable names. Strict equality is required         
!                                                                       
      DO i = 1, var_num 
      IF (string (1:laenge) .eq.var_name (i) (1:var_l (i) ) ) THEN 
         IF(var_entry(i)>0) RETURN     ! This is a variable array
         IF (var_type (i) .eq.      IS_REAL) THEN 
            var_val (i) = wert 
         ELSEIF (var_type (i) .eq.      IS_INTE) THEN 
            var_val (i) = nint (wert) 
         ELSEIF (var_type (i) .eq.      IS_CHAR) THEN 
            var_char (i) = dummy (1:length) 
         ENDIF 
         IF(string (1:laenge) == 'REF_DIMENSION') THEN
            IF(var_val (i)>UBOUND(ref_para,1)) THEN
               CALL alloc_ref_para(NINT(wert))
            ENDIF
         ENDIF
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
      ENDDO 
END SUBROUTINE upd_variable                   
!
!*****7**************************************************************** 
!
SUBROUTINE rese_variables(is_diffev)
!-
!     Removes all user variables
!+
USE variable_mod
IMPLICIT none 
!
LOGICAL, INTENT(IN)    :: is_diffev
!
INTEGER :: i, j
INTEGER :: ndel
!
DO i=var_sys+1, var_num
   IF(var_diff(i) .EQV. is_diffev) THEN
      var_name( i) = ' '
      var_char( i) = ' '
      var_l   ( i) = 1
      var_type( i) = 0
      var_diff( i) = .FALSE.
      var_val ( i) = 0
      IF(var_entry(i)/=0) THEN      ! GOT An array
         IF(ALLOCATED(var_field(var_entry(i))%var_value)) &
          DEALLOCATE(var_field(var_entry(i))%var_value)
         var_field(var_entry(i))%var_shape(:) = 0
         var_entry(i) = 0
      ENDIF
   ENDIF
ENDDO
ndel = 0
i = var_sys+1
remove: DO 
   IF(var_name(i) == ' ') THEN
      ndel = ndel +1
      IF(i>var_num-ndel) EXIT remove
      DO j=i, var_num-ndel
         var_name (j) = var_name (j+1)
         var_char (j) = var_char (j+1)
         var_l    (j) = var_l    (j+1)
         var_type (j) = var_type (j+1)
         var_diff (j) = var_diff (j+1)
         var_val  (j) = var_val  (j+1)
         var_entry(j) = var_entry(j+1)
      ENDDO
   ELSE
      i=i+1
   ENDIF
ENDDO remove
var_num = var_num - ndel
!
END SUBROUTINE rese_variables 
!
!*****7**************************************************************** 
 SUBROUTINE del_variables (MAXW, ianz, cpara, lpara, is_diffev)
!-
!     Removes all or specified user variables
!+
USE errlist_mod
USE variable_mod
!
IMPLICIT none 
!
INTEGER                               , INTENT(IN)    :: MAXW
INTEGER                               , INTENT(IN)    :: ianz
CHARACTER(LEN=1024), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
INTEGER            , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
LOGICAL                               , INTENT(IN)    :: is_diffev
!
INTEGER :: i, j, k, upper
!
LOGICAL :: str_comp
!
IF(ianz==1) THEN
   CALL rese_variables (is_diffev)
ELSEIF(str_comp (cpara (2) , 'all', 3, lpara(1), 3) ) THEN 
   CALL rese_variables (is_diffev)
ELSE
   main: DO i=2, ianz
      upper = var_num
      loop: DO j=var_sys+1, upper
         IF(cpara(i) == var_name(j)) THEN 
            IF(var_diff(j) .EQV. is_diffev) THEN
               IF(var_entry(j)/=0) THEN      ! GOT An array
                  IF(ALLOCATED(var_field(var_entry(j))%var_value)) &
                    DEALLOCATE(var_field(var_entry(j))%var_value)
                  var_field(var_entry(j))%var_shape(:) = 0
                  var_entry(j) = 0
               ENDIF
               DO k=j+1, var_num
                  var_name (k-1) = var_name(k)
                  var_char (k-1) = var_char(k)
                  var_l    (k-1) = var_l(k)
                  var_type (k-1) = var_type(k)
                  var_diff (k-1) = var_diff(k)
                  var_val  (k-1) = var_val(k)
                  var_entry(k-1) = var_entry(k)
               ENDDO
               var_num = var_num - 1
               CYCLE main
            ELSE
               ier_num =  -39
               ier_typ = ER_FORT
               ier_msg(1) = 'Offending variable: '//cpara(i)(1:LEN_TRIM(cpara(i)))
               RETURN
            ENDIF
         ENDIF
      ENDDO loop
      ier_num = -24
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending variable: '//cpara(i)(1:LEN_TRIM(cpara(i)))
      RETURN
   ENDDO main
ENDIF
!
END SUBROUTINE del_variables 
!
!*****7**************************************************************** 
!
SUBROUTINE show_variables 
!-                                                                      
!     shows the variable definitions                                    
!+                                                                      
USE constants_mod
USE blanks_mod
USE prompt_mod 
USE variable_mod
IMPLICIT none 
!
!
CHARACTER(LEN=6), DIMENSION(0:1) :: cdiffev
CHARACTER(LEN=1024) :: string
INTEGER :: i , j
!
INTEGER :: length
INTEGER :: len_str 
DATA cdiffev/'DIFFEV','      '/
!
IF (var_num.gt.0) THEN 
   WRITE (output_io, 2000) var_num, VAR_MAX 
   DO i = 1, var_num 
      j=1
      IF(var_diff(i)) j=0
      string = ' '
      length = 0
      IF(var_entry(i)>0) THEN
         IF(var_field(var_entry(i))%var_shape(2)>1) THEN
            WRITE(string,'(a1,i8,a1,i8,a1)') '[',var_field(var_entry(i))%var_shape(1), &
            ',',var_field(var_entry(i))%var_shape(2),']'
            length = 19
            CALL rem_bl(string, length)
         ELSE
            WRITE(string,'(a1,i8,a1)') '[',var_field(var_entry(i))%var_shape(1),']'
            length = 10
            CALL rem_bl(string, length)
         ENDIF
      ENDIF
      IF (var_type (i) .eq.      IS_REAL) THEN 
         IF (ABS(var_val(i)) .lt.1e-5.or. ABS(var_val(i)).ge.1e6) THEN
            WRITE(output_io, 2100) var_name (i), var_val (i) , cdiffev(j), string(1:length)
         ELSE 
            WRITE(output_io, 2150) var_name (i), var_val (i) , cdiffev(j), string(1:length)
         ENDIF 
      ELSEIF (var_type (i) .eq.      IS_INTE) THEN 
         WRITE(output_io, 2200) var_name (i), nint (var_val (i) ) , cdiffev(j), string(1:length)
      ELSEIF (var_type (i) .eq.      IS_CHAR) THEN 
         WRITE(output_io, 2300) var_name(i), var_char(i)(1:len_str(var_char(i))), cdiffev(j)
      ENDIF 
   ENDDO 
ELSE 
   WRITE (output_io, 2050) VAR_MAX 
ENDIF 
WRITE (output_io, * ) 
 2000 FORMAT    (' User defined variables',i4,                          &
     &                  ' Maximum number: ',i3/)                        
 2050 FORMAT    (' No user defined variables',                          &
     &                  ' Maximum number: ',i3/)                        
 2100 FORMAT    (' Name: ',a30,' = ',8x,e20.8e2,' Real      ',a,a) 
 2150 FORMAT    (' Name: ',a30,' = ',8x,f16.8  ,'     Real      ',a,a) 
 2200 FORMAT    (' Name: ',a30,' = ',i15,13x,   ' Integer   ',a,a) 
 2300 FORMAT    (' Name: ',a30,' = ''',a,     ''' Character ',a) 
END SUBROUTINE show_variables                 
!
!*****7**************************************************************** 
!
SUBROUTINE validate_variable (zeile, lp) 
!-                                                                      
!       Checks whether the variable name is legal                       
!     All intrinsic function names, variable names and parts of         
!     these names are illegal                                           
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@mail.uni-wuerzburg.de)    
!+                                                                      
      USE reserved_mod
      USE charact_mod
      USE errlist_mod 
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: zeile 
      INTEGER,           INTENT(IN) :: lp 
!                                                                       
      INTEGER i, ii 
      INTEGER :: length
      LOGICAL lok 
!                                                                       
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
      main:DO i = 1, lib_reserved_n 
!        IF (INDEX (lib_reserved (i), zeile (1:lp) ) .NE.0) THEN 
         length = MAX(LEN_TRIM(lib_reserved(i)), LEN_TRIM(zeile(1:lp)))
         length = MIN(length, LEN(lib_reserved), LEN(zeile))
         IF ((lib_reserved (i)(1:length)== zeile (1:length) ) ) THEN 
            ier_num = - 25 
            ier_typ = ER_FORT 
            EXIT main
         ENDIF 
      ENDDO main 
!                                                                       
!     Check against variables/functions of the main program             
!                                                                       
      IF (ier_num.eq.0) THEN 
         CALL p_validate_var_spec (zeile, lp) 
      ENDIF 
!                                                                       
!     Check that the name contains only letters, Numbers and the "_"    
!                                                                       
      IF (ier_num.eq.0) THEN 
         lok = .true. 
         DO i = 1, lp 
            ii = iachar (zeile (i:i) ) 
            lok = lok.and. (zero.le.ii.and.ii.le.nine.or.aa.le.ii.and. & 
                            ii.le.zz .or.u.eq.ii.or.a.le.ii.and.ii.le.z)                               
         ENDDO 
         IF (.not.lok) THEN 
            ier_num = - 26 
            ier_typ = ER_FORT 
         ENDIF 
      ENDIF 
!                                                                       
END SUBROUTINE validate_variable              
!
END MODULE do_variable_mod
