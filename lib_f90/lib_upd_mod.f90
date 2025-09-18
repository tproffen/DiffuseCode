MODULE lib_upd_mod
!
interface lib_check_dim_para
   module procedure lib_check_dim_para_D1, lib_check_dim_para_D2
end interface lib_check_dim_para
!
interface lib_set_para
   module procedure lib_set_para_inte_D1, &
                    lib_set_para_inte_D2, &
                    lib_set_para_real_D1, &
                    lib_set_para_real_D2, &
                    lib_set_para_real_D2_22, &
                    lib_set_para_char_D2
end interface lib_set_para
!
private
!
public lib_ersetz_para
public lib_upd_para
public lib_set_para
public lib_check_dim_para

CONTAINS
!
SUBROUTINE lib_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
USE constants_mod
USE blanks_mod
USE errlist_mod 
USE lib_length
USE param_mod 
USE precision_mod
USE precision_command_mod
USE random_mod 
USE variable_mod
IMPLICIT none 
!                                                                       
INTEGER,                    INTENT(IN   ) :: ikl
INTEGER,                    INTENT(IN   ) :: iklz
CHARACTER (LEN=*),          INTENT(INOUT) :: string 
INTEGER,                    INTENT(INOUT) :: ll
INTEGER,                    INTENT(IN   ) :: maxw
REAL(KIND=PREC_DP), DIMENSION(1:maxw), INTENT(IN   ) :: ww
INTEGER,                    INTENT(IN   ) :: ianz
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) :: zeile 
!                                                                       
INTEGER :: laenge, ltyp, kpara, kpara2
INTEGER :: i
INTEGER :: lstr
INTEGER :: lcomm 
LOGICAL :: success = .FALSE.
!                                                                       
laenge = ll 
ltyp = 1 
zeile = ' ' 
kpara = NINT (ww (1) ) 
kpara2 = 0
IF (ianz.ge.2) THEN 
   kpara2 = NINT (ww (2) ) 
ENDIF 
!                                                                       
lcomm = length_com (string, ikl) 
success = .FALSE.
!
!  Test for Expressions stored in var_exp
!
if(ikl-var_l(VAR_EXPRESSION)>=1) then      ! Enough space to search for "EXPR"
   if(var_name(VAR_EXPRESSION) == string(ikl -  var_l(VAR_EXPRESSION):ikl-1)) then ! found "EXPR"
      lcomm = var_l(VAR_EXPRESSION)
      zeile(1:ikl - lcomm-1) = string(1:ikl - lcomm-1)
      zeile(ikl-lcomm:ikl+len_trim(var_expr(kpara))) = var_expr(kpara)(1:len_trim(var_expr(kpara)))
      success = .TRUE.
   endif
endif
!
!  Test User defined variable arrays
!
search_var: DO i=var_sys+1, var_num
!write(*,*) 'search ', string(ikl - lcomm:ikl - 1),ikl, lcomm, iklz,ianz,'|',string(1:50)
!  IF(var_name(i) == string(ikl - lcomm:ikl - 1)) THEN
   IF(ikl-var_l(i)<1) CYCLE search_var
   IF(var_name(i) == string(ikl - var_l(i):ikl - 1)) THEN
      lcomm = var_l(i)
      IF(var_entry(i)>0) THEN
!write(*,*) ' kpara ', kpara, kpara2,var_field(var_entry(i))%var_shape(:), maxw
         IF(0<kpara .AND. kpara<=var_field(var_entry(i))%var_shape(1) ) THEN
            IF((ianz==2 .AND. 0<kpara2 .AND. kpara2<=var_field(var_entry(i))%var_shape(2)) .OR.  &
               (ianz==1 .AND. kpara2==0 .AND. var_field(var_entry(i))%var_shape(2)==1 ) )  THEN
               kpara2 = MAX(1, kpara2)
               zeile(1:ikl - lcomm-1) = string(1:ikl - lcomm-1)
               IF(var_type(i)==      IS_INTE) THEN
                  WRITE(zeile(ikl - lcomm:ikl + PREC_WIDTH-2),PREC_F_INTE)             &
                  NINT(var_field(var_entry(i))%var_value(kpara,kpara2))
               ELSEIF(var_type(i)==      IS_REAL) THEN
                  WRITE(zeile(ikl - lcomm:ikl + PREC_WIDTH-2),PREC_F_REAL)             &
                  var_field(var_entry(i))%var_value(kpara,kpara2)
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF(var_type(i)==      IS_CHAR) THEN
                  lstr = LEN_TRIM(var_field(var_entry(i))%var_char(kpara,kpara2))
                  zeile(ikl - lcomm:ikl + lstr                                                    ) = &
                        ''''//var_field(var_entry(i))%var_char(kpara,kpara2)(1:lstr)//''''
               ENDIF
               success = .TRUE.
!write(*,*) 'PLACED ', zeile(1:50)
               EXIT search_var
            ELSE
               ier_num = -40
               ier_typ = ER_FORT
               RETURN
            ENDIF
         ELSE
            ier_num = -40
            ier_typ = ER_FORT
            RETURN
         ENDIF
      ENDIF
   ENDIF
ENDDO search_var
!
!                                                                       
IF(.NOT.success) THEN
lcomm = length_com (string, ikl) 
   IF (lcomm.eq.1) THEN 
!                                                                       
      IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1: ikl - lcomm - 1)
         IF (string (ikl - 1:ikl - 1) .eq.'i') THEN 
            IF (ianz.eq.1) THEN 
               IF (0.le.kpara.and.kpara.le.MAXPAR) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) inpara (   &
                  kpara)                                                
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'r') THEN 
            IF (ianz.eq.1) THEN 
               IF (0.le.kpara.and.kpara.le.MAXPAR) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) rpara (&
                  kpara)                                                
                  zeile (ikl + PREC_MANTIS-1:ikl + PREC_MANTIS-1) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!
   ELSEIF (lcomm.eq.3) THEN 
!                                                                       
         IF (string (ikl - 3:ikl - 1) .eq.'res') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) THEN 
                  WRITE (zeile (ikl - 3:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  res_para (kpara)                                      
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
   ELSEIF (lcomm.eq.4) THEN 
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'seed') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
      WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE) idum
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!
   ELSEIF(lcomm==6) THEN
!                                                                       
         IF (string (ikl - 6:ikl - 1) .eq.'F_PARA') THEN
            IF (ianz.eq.1) THEN
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)
               IF (0.lt.kpara.and.kpara.le.MAXPAR) THEN
                  WRITE(zeile(ikl - 6:ikl + PREC_WIDTH-2), PREC_F_REAL) kupl_para(kpara)
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd'
               ELSE
                  ier_num = -133
                  ier_typ = ER_APPL
                  RETURN
               ENDIF
            ELSE
               ier_num = - 13
               ier_typ = ER_FORT
               RETURN
            ENDIF
         ELSE
            ier_num = - 2
            ier_typ = ER_FORT
         ENDIF
!
   ELSEIF(lcomm==7) THEN
!                                                                       
         IF (string (ikl - 7:ikl - 1) .eq.'F_DERIV' .OR.   &
             string (ikl - 7:ikl - 1) .eq.'F_DeRIV' .OR.   &
             string (ikl - 7:ikl - 1) .eq.'f_deriv') THEN
            IF (ianz.eq.1) THEN
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)
               IF (0.lt.kpara.and.kpara.le.MAXPAR) THEN
                  WRITE(zeile(ikl - 7:ikl + PREC_WIDTH-2), PREC_F_REAL) kupl_deriv(kpara)
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd'
               ELSE
                  ier_num = -133
                  ier_typ = ER_APPL
                  RETURN
               ENDIF
            ELSE
               ier_num = - 13
               ier_typ = ER_FORT
               RETURN
            ENDIF
         ELSE
            ier_num = - 2
            ier_typ = ER_FORT
         ENDIF
!
   ELSEIF (lcomm.eq.8) THEN
!
        IF (string (ikl - 8:ikl - 1) .eq.'ref_para') THEN
            IF (ianz.eq.1) THEN
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)
               IF (0.lt.kpara.and.kpara.le.MAXPAR_REF   ) THEN
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  ref_para (kpara)
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd'
               ELSE
                  ier_num = -133
                  ier_typ = ER_APPL
                  RETURN
               ENDIF
            ELSE
               ier_num = - 13
               ier_typ = ER_FORT
               RETURN
            ENDIF
         ELSE
            ier_num = - 2
            ier_typ = ER_FORT
         ENDIF
!
   ELSE 
!
         ier_num = - 2 
         ier_typ = ER_FORT 
!
   ENDIF 
ENDIF 
!
      IF (ier_num.eq.0) THEN 
         ll = laenge+PREC_WIDTH - ltyp - (iklz - ikl + 1) 
         IF (iklz + 1.le.laenge) zeile (ikl + PREC_WIDTH-1:ll) = string (iklz + 1:&
         laenge)                                                        
         string = zeile 
      ELSE 
         ll = min (40, laenge) 
         WRITE (ier_msg (1), '(a)') string (1:ll) 
      ENDIF 

      ll  = LEN_TRIM(string)
!                                                                       
END SUBROUTINE    lib_ersetz_para                    
!*****7*************************************************************************
SUBROUTINE lib_upd_para (ctype, lower_limit, upper_limit, maxw, lrange, wert, ianz,  &
           cstring, substr, lexpr, line_expression) 
!-                                                                      
!       updates the parameter specified by ctype, index lower_limit  to the      
!       new value of wert                                               
!+                                                                      
USE constants_mod
USE errlist_mod 
USE param_mod 
USE variable_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*),          INTENT(IN) :: ctype 
INTEGER,                    INTENT(IN) :: maxw
INTEGER,                    INTENT(IN) :: ianz 
INTEGER, DIMENSION(1:MAXW), INTENT(inout) :: lower_limit
integer, dimension(1:MAXW), intent(inout) :: upper_limit
logical                   , intent(in) :: lrange
REAL(KIND=PREC_DP)        , INTENT(IN) :: wert 
CHARACTER (LEN=*),          INTENT(IN) :: cstring 
INTEGER, DIMENSION(2)     , INTENT(IN) :: substr ! Indices of substring
logical                   , intent(in) :: lexpr
CHARACTER (LEN=*),          INTENT(inout) :: line_expression 
!
INTEGER :: i, ww2
integer, dimension(1:2) :: user_lower_limit
integer, dimension(1:2) :: user_upper_limit
logical :: llimited                ! Parameter value "wert" needs to be limited
real(kind=PREC_DP), dimension(2) :: r_wert_limit     ! Low and high value limit
integer           , dimension(2) :: i_wert_limit     ! Low and high value limit
!
llimited = .false.
r_wert_limit(1) = -huge(1.0_PREC_DP)   ! Default to no limits
r_wert_limit(2) =  huge(1.0_PREC_DP)
i_wert_limit(1) = -huge(1          )
i_wert_limit(2) =  huge(1          )
!
search_var: DO i=var_sys+1, var_num
   IF(var_name(i) == ctype(1:len_trim(ctype))   ) THEN
      IF(var_entry(i)>0) THEN
!write(*,*) ' VARIABLE ARRAY ', i, ' ', var_name(i)(1:len_trim(var_name(i))), var_field(var_entry(i))%var_shape(1), var_field(var_entry(i))%var_shape(2),' MAXW ', maxw, ianz
!write(*,*) ' LOWER_LIMIT    ', lower_limit(1:MAXW), lrange, lexpr
!write(*,*) ' UPPER_LIMIT    ', upper_limit(1:MAXW)
!write(*,'(a,15f8.3)') ' CURRENT', var_field(var_entry(i))%var_value(1:10,1) !var_field(var_entry(i))%var_shape(1),1)
         user_lower_limit = 1
         user_upper_limit = 1
         user_upper_limit(1:MAXW) = var_field(var_entry(i))%var_shape(1:MAXW)
!write(*,*) ' USER  LOWER    ', user_lower_limit(1:MAXW)
!write(*,*) ' USER  UPPER    ', user_upper_limit(1:MAXW)
         if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, user_lower_limit, user_upper_limit)         ) then
!write(*,*) 'Correct shape ', ianz, MAXW
           if(var_type(i)==  IS_CHAR) then
!write(*,*) ' CHARACTER VARIABLE ', shape(var_field(var_entry(i))%var_char)
              call lib_set_para(user_lower_limit,user_upper_limit, &
              var_field(var_entry(i))%var_char(user_lower_limit(1):user_upper_limit(1),user_lower_limit(2):user_upper_limit(2)), &
              1,1,  2, lower_limit, upper_limit, lrange, lexpr, substr, cstring, line_expression)
!write(*,*) ' IERROR ', ier_num, ier_typ
!do j = 1, 10
!write(*,'(a,a20   )') ' RESULT ', var_field(var_entry(i))%var_char(j   ,1) !var_field(var_entry(i))%var_shape(1),1)
!enddo
              return !cycle search_var
           else
           if(MAXW==1) then    ! 1D field
!write(*,*) ' SETTING A 1 D ARRAY ', MAXW
              call lib_set_para(user_lower_limit,user_upper_limit, &
              var_field(var_entry(i))%var_value(user_lower_limit(1):user_upper_limit(1),user_lower_limit(2):user_upper_limit(2)), 1,1,  2, &
              lower_limit, upper_limit, lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
!write(*,*) ' IERROR ', ier_num, ier_typ
!write(*,'(a,15f8.3)') ' RESULT ', var_field(var_entry(i))%var_value(1:10,1) !var_field(var_entry(i))%var_shape(1),1)
              return !cycle search_var
           elseif(MAXW==2) then
!write(*,*) ' SETTING A 2 D ARRAY ', MAXW
              call lib_set_para(user_lower_limit,user_upper_limit, &
              var_field(var_entry(i))%var_value(user_lower_limit(1):user_upper_limit(1),user_lower_limit(2):user_upper_limit(2)), 2, &
              lower_limit, upper_limit, lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
!write(*,*) ' IERROR ', ier_num, ier_typ
!do j=1, 3
!write(*,'(a,15f8.3)') ' RESULT ', var_field(var_entry(i))%var_value(j,1:4) !var_field(var_entry(i))%var_shape(1),1)
!enddo
              return !cycle search_var
           endif
         endif
         else
!write(*,*) 'FALSE shape ', ianz, MAXW
            return
         endif
         return !cycle search_var
!OBSOLETE ?
         IF(0<lower_limit(1) .AND. lower_limit(1)<=var_field(var_entry(i))%var_shape(1) ) THEN
            IF(maxw==1) THEN
               IF(var_field(var_entry(i))%var_shape(2)>1) THEN
                  ier_num = -40
                  ier_typ = ER_FORT
                  RETURN
               ELSE
                  ww2 = 1
               ENDIF
            ELSEIF(maxw==2) THEN
               IF(0>=lower_limit(2) .OR. lower_limit(2)> var_field(var_entry(i))%var_shape(2)) THEN
                  ier_num = -40
                  ier_typ = ER_FORT
                  RETURN
               ELSE
                  ww2 = lower_limit(2)
               ENDIF
            ELSE
               ier_num = -40
               ier_typ = ER_FORT
               RETURN
            ENDIF
               IF(var_type(i)==      IS_INTE) THEN
                  var_field(var_entry(i))%var_value(lower_limit(1),ww2) = INT(wert)
               ELSEIF(var_type(i)==      IS_REAL) THEN
                  var_field(var_entry(i))%var_value(lower_limit(1),ww2) = wert
               ELSEIF(var_type(i)==  IS_CHAR) THEN
                  var_field(var_entry(i))%var_char(lower_limit(1),ww2)(substr(1):substr(2)) = cstring(1:len_trim(cstring))
!write(*,*) ' SET ENTRY ', ww(1),ww2,' in ', var_name(i)(1:len_trim(var_name(i))) ,' TO ', cstring(1:len_trim(cstring))
               ENDIF
            RETURN
         ELSE
            ier_num = -40
            ier_typ = ER_FORT
            RETURN
         ENDIF
      ENDIF
   ENDIF
ENDDO search_var
!                                                                       
IF (ctype.eq.'i') THEN 
!
   if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, 0, MAXPAR)) then
      call lib_set_para(0,MAXPAR, inpara, MAXW, lower_limit, upper_limit, lrange, &
           lexpr, wert, line_expression, llimited, i_wert_limit)
   else
      return
   endif
!
ELSEIF (ctype.eq.'r') THEN 
!
   if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, 0, MAXPAR)) then
      call lib_set_para(0,MAXPAR, rpara, MAXW, lower_limit, upper_limit, lrange, &
           lexpr, wert, line_expression, llimited, r_wert_limit)
   else
      return
   endif
!                                                                       
ELSEIF (ctype.eq.'res') THEN 
!
   if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, 0, MAXPAR_RES)) then
      call lib_set_para(0,MAXPAR_RES, res_para, MAXW, lower_limit, upper_limit, lrange, &
           lexpr, wert, line_expression, llimited, r_wert_limit)
   else
      return
   endif
!                                                                       
ELSEIF (ctype.eq.'F_DERIV') THEN
!
   if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, 0, MAXPAR)) then
      call lib_set_para(0,MAXPAR, kupl_deriv, MAXW, lower_limit, upper_limit, lrange, &
           lexpr, wert, line_expression, llimited, r_wert_limit)
   else
      return
   endif
!                                                                       
!        IF (ianz.eq.1) THEN 
!           IF (0.le.lower_limit (1) .and.lower_limit (1) .le.MAXPAR) THEN 
!              kupl_deriv (lower_limit (1) ) = wert 
!           ELSE 
!              ier_num = - 8 
!              ier_typ = ER_FORT 
!           ENDIF 
!        ELSE 
!           ier_num = - 13 
!           ier_typ = ER_FORT 
!           RETURN 
!        ENDIF 
ELSEIF (ctype.eq.'F_PARA') THEN
!
   if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, 0, MAXPAR)) then
      call lib_set_para(0,MAXPAR, kupl_para, MAXW, lower_limit, upper_limit, lrange, &
           lexpr, wert, line_expression, llimited, r_wert_limit)
   else
      return
   endif
!                                                                       
!        IF (ianz.eq.1) THEN 
!           IF (0.le.lower_limit (1) .and.lower_limit (1) .le.MAXPAR) THEN 
!              kupl_para (lower_limit (1) ) = wert 
!           ELSE 
!              ier_num = - 8 
!              ier_typ = ER_FORT 
!           ENDIF 
!        ELSE 
!           ier_num = - 13 
!           ier_typ = ER_FORT 
!           RETURN 
!        ENDIF 
ELSEIF (ctype.eq.'ref_para') THEN
!
   if(lib_check_dim_para(ianz, lower_limit, upper_limit, lrange, 0, MAXPAR_REF)) then
      call lib_set_para(0,MAXPAR_REF,  ref_para, MAXW, lower_limit, upper_limit, lrange, &
           lexpr, wert, line_expression, llimited, r_wert_limit)
   else
      return
   endif
!                                                                       
!        IF (ianz.eq.1) THEN 
!           IF (0.le.lower_limit (1) .and.lower_limit (1) .le.MAXPAR_REF) THEN 
!              ref_para (lower_limit (1) ) = wert 
!           ELSE 
!              ier_num = - 8 
!              ier_typ = ER_FORT 
!           ENDIF 
!        ELSE 
!           ier_num = - 13 
!           ier_typ = ER_FORT 
!           RETURN 
!        ENDIF 
ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
         WRITE (ier_msg (1), '(a)') ctype 
ENDIF 
!
END SUBROUTINE lib_upd_para                       
!
!*******************************************************************************
!*******************************************************************************
!
subroutine lib_set_para_inte_D1(LOW, HIGH, user_variable, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, i_wert_limit)
!-
! sets a variable to the individual value or a range
! Version for 1D inte valued user variable
! user_variable[ fixed, insert_here] = wert
!+
!
use berechne_mod
use errlist_mod
use do_replace_expr_mod
use precision_mod
!
implicit none
!
integer                                        , intent(in   ) :: LOW             ! Dimensions for user variable
integer                                        , intent(in   ) :: HIGH            ! Dimensions for user variable
integer           , dimension(LOW:HIGH        ), intent(inout) :: user_variable
integer                                        , intent(in   ) :: MAXW            ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit     ! Limits for replacement range
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit     ! Limits for replacement range
logical                                        , intent(in   ) :: lrange          ! User provided a range
logical                                        , intent(in   ) :: lexpr           ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert            ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited        ! Parameter value "wert" needs to be limited
integer           , dimension(2)               , intent(in   ) :: i_wert_limit    ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer :: i
integer :: length
real(kind=PREC_DP) :: wwert
!
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         wwert  = berechne(line, length)                             ! Calculate current value of EXPR
         if(ier_num/=0) return
         if(llimited .and. (i_wert_limit(1)> int(wwert) .or. int(wwert)> i_wert_limit(2))) then
            ier_num = -50
            ier_typ = ER_FORT
         else
            user_variable(i) = int(wwert)
         endif
      enddo
   else
      if(llimited .and. (i_wert_limit(1)> int(wert) .or. int(wert)> i_wert_limit(2))) then
         ier_num = -50
         ier_typ = ER_FORT
      else
         user_variable(lower_limit(1):upper_limit(1) ) = int(wert)
      endif
   endif
else 
   if(llimited .and. (i_wert_limit(1)> int(wert) .or. int(wert)> i_wert_limit(2))) then
      ier_num = -50
      ier_typ = ER_FORT
   else
      user_variable(lower_limit(1) ) = int(wert )
   endif
ENDIF 
!
end subroutine lib_set_para_inte_D1
!
!*******************************************************************************
!
subroutine lib_set_para_real_D1(LOW, HIGH, user_variable, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
!-
! sets a variable to the individual value or a range
! Version for 1D real valued user variable
! user_variable[ insert_here] = wert
!+
!
use berechne_mod
use errlist_mod
use do_replace_expr_mod
use precision_mod
!
implicit none
!
integer                                        , intent(in   ) :: LOW             ! Dimensions for user variable
integer                                        , intent(in   ) :: HIGH            ! Dimensions for user variable
real(kind=PREC_DP), dimension(LOW:HIGH        ), intent(inout) :: user_variable
integer                                        , intent(in   ) :: MAXW            ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange          ! User provided a range
logical                                        , intent(in   ) :: lexpr           ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert            ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited        ! Parameter value "wert" needs to be limited
real(kind=PREC_DP), dimension(2)               , intent(in   ) :: r_wert_limit    ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer :: i
integer :: length
real(kind=PREC_DP) :: wwert
!
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         wwert  = berechne(line, length)                             ! Calculate current value of EXPR
         if(ier_num/=0) return
         if(llimited .and. (r_wert_limit(1)> int(wwert) .or. int(wwert)> r_wert_limit(2))) then
            ier_num = -50
            ier_typ = ER_FORT
         else
            user_variable(i) = wwert 
         endif
      enddo
   else
      if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
         ier_num = -50
         ier_typ = ER_FORT
      else
         user_variable(lower_limit(1):upper_limit(1) ) = wert 
      endif
   endif
else 
   if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
      ier_num = -50
      ier_typ = ER_FORT
   else
      user_variable(lower_limit(1) ) = wert 
   endif
ENDIF 
!
end subroutine lib_set_para_real_D1
!
!*******************************************************************************
!
subroutine lib_set_para_inte_D2(LOW, HIGH, user_variable, index_var, index1, MAXW, &
           lower_limit, upper_limit, lrange, lexpr, wert, line_expression, llimited, &
           i_wert_limit)
!-
! Interface to specific 2D real valued "set" subroutines 
!+
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
integer           , dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index_var      ! Variable Index at dimension index_var
integer                                        , intent(in   ) :: index1         ! Index at dimension 1
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited        ! Parameter value "wert" needs to be limited
integer           , dimension(2)               , intent(in   ) :: i_wert_limit    ! Low and high value limit
!
if(ubound(HIGH,1)==2) then
   if(index_var==2) then
   call lib_set_para_inte_D2_1(LOW, HIGH, user_variable, index1, MAXW, lower_limit, upper_limit, lrange, lexpr, wert, line_expression, llimited, i_wert_limit)
   endif
endif
!
end subroutine lib_set_para_inte_D2
!
!*******************************************************************************
!
subroutine lib_set_para_inte_D2_1(LOW, HIGH, user_variable, index1, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, i_wert_limit)
!-
! sets a variable to the individual value or a range
! Version for 2D inte valued user variable, insert at second index
! user_variable[ fixed, insert_here] = wert
!+
!
use berechne_mod
use errlist_mod
use do_replace_expr_mod
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
integer           , dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index1         ! Index at dimension 1
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
integer           , dimension(2)               , intent(in   ) :: i_wert_limit     ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer :: i
integer :: length
real(kind=PREC_DP) :: wwert
!
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         wwert  = berechne(line, length)                             ! Calculate current value of EXPR
         if(ier_num/=0) return
         if(llimited .and. (i_wert_limit(1)> int(wwert) .or. int(wwert)> i_wert_limit(2))) then
            ier_num = -50
            ier_typ = ER_FORT
         else
            user_variable(INDEX1, i) = int(wwert)
         endif
      enddo
   else
      if(llimited .and. (i_wert_limit(1)> int(wert) .or. int(wert)> i_wert_limit(2))) then
         ier_num = -50
         ier_typ = ER_FORT
      else
         user_variable(INDEX1, lower_limit(1):upper_limit(1) ) = int(wert )
      endif
   endif
else 
   if(llimited .and. (i_wert_limit(1)> int(wert) .or. int(wert)> i_wert_limit(2))) then
      ier_num = -50
      ier_typ = ER_FORT
   else
      user_variable(INDEX1, lower_limit(1) ) = int(wert)
   endif
ENDIF 
!
end subroutine lib_set_para_inte_D2_1
!
!*******************************************************************************
!
subroutine lib_set_para_real_D2(LOW, HIGH, user_variable, index_var, index1, MAXW, &
           lower_limit, upper_limit, lrange, lexpr, wert, line_expression, llimited, &
           r_wert_limit)
!-
! Interface to specific 2D real valued "set" subroutines 
!+
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
real(kind=PREC_DP), dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index_var      ! Variable Index at dimension index_var
integer                                        , intent(in   ) :: index1         ! Index at dimension 1
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
real(kind=PREC_DP), dimension(2)               , intent(in   ) :: r_wert_limit     ! Low and high value limit
!
if(ubound(HIGH,1)==2) then
   if(index_var==1) then
     call lib_set_para_real_D2_1(LOW, HIGH, user_variable, index1, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
   elseif(index_var==2) then
     call lib_set_para_real_D2_2(LOW, HIGH, user_variable, index1, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
   endif
endif
!
end subroutine lib_set_para_real_D2
!
!*******************************************************************************
!
subroutine lib_set_para_real_D2_1(LOW, HIGH, user_variable, index1, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
!-
! sets a variable to the individual value or a range
! Version for 2D real valued user variable, insert at first index
! user_variable[ fixed, insert_here] = wert
!+
!
use berechne_mod
use errlist_mod
use do_replace_expr_mod
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
real(kind=PREC_DP), dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index1         ! Index at dimension 1
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
real(kind=PREC_DP), dimension(2)               , intent(in   ) :: r_wert_limit     ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer :: i
integer :: length
real(kind=PREC_DP) :: wwert
!
!write(*,*) ' IN lib_set_para_real_D2_1', MAXW
!write(*,*) ' VARIABLE DIMENSION', LOW, ' <> ', HIGH
!write(*,*) '             lower ', lower_limit, ' | ', maxw, wert, ' > ', index1
!write(*,*) '             upper ', upper_limit, ' | ', lrange, lexpr
!write(*,*) ' limited           ', llimited, r_wert_limit(1:2), r_wert_limit(1)> (wert), (wert)> r_wert_limit(2)
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         wwert  = berechne(line, length)                             ! Calculate current value of EXPR
         if(ier_num/=0) return
         if(llimited .and. (r_wert_limit(1)> (wwert) .or. (wwert)> r_wert_limit(2))) then
            ier_num = -50
            ier_typ = ER_FORT
         else
!write(*,*) ' SETTING FIXED to EXPRESSION  ',wwert
            user_variable(i, INDEX1) = wwert 
         endif
      enddo
   else
      if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
         ier_num = -50
         ier_typ = ER_FORT
      else
!write(*,*) ' SETTING RANGE to fixed value ', wert
         user_variable(lower_limit(1):upper_limit(1), INDEX1 ) = wert 
      endif
   endif
else 
   if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
      ier_num = -50
      ier_typ = ER_FORT
   else
!write(*,*) ' SETTING FIXED to fixed value ', lower_limit(1), INDEX1, wert
      user_variable(lower_limit(1), INDEX1 ) = wert 
   endif
ENDIF 
!
end subroutine lib_set_para_real_D2_1
!
!*******************************************************************************
!
subroutine lib_set_para_real_D2_2(LOW, HIGH, user_variable, index1, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
!-
! sets a variable to the individual value or a range
! Version for 2D real valued user variable, insert at second index
! user_variable[ fixed, insert_here] = wert
!+
!
use berechne_mod
use errlist_mod
use do_replace_expr_mod
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
real(kind=PREC_DP), dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index1         ! Index at dimension 1
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
real(kind=PREC_DP), dimension(2)               , intent(in   ) :: r_wert_limit     ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer :: i
integer :: length
real(kind=PREC_DP) :: wwert
!
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         wwert  = berechne(line, length)                             ! Calculate current value of EXPR
         if(ier_num/=0) return
         if(llimited .and. (r_wert_limit(1)> (wwert) .or. (wwert)> r_wert_limit(2))) then
            ier_num = -50
            ier_typ = ER_FORT
         else
            user_variable(INDEX1, i) = wwert 
         endif
      enddo
   else
      if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
         ier_num = -50
         ier_typ = ER_FORT
      else
         user_variable(INDEX1, lower_limit(1):upper_limit(1) ) = wert 
      endif
   endif
else 
   if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
      ier_num = -50
      ier_typ = ER_FORT
   else
      user_variable(INDEX1, lower_limit(1) ) = wert 
   endif
ENDIF 
!
end subroutine lib_set_para_real_D2_2
!
!*******************************************************************************
!
subroutine lib_set_para_real_D2_22(LOW, HIGH, user_variable, MAXW, lower_limit, upper_limit, &
           lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
!-
! sets a variable to the individual value or a range
! Version for 2D real valued user variable, insert at both indices
! user_variable[ insert_here, insert_here] = wert
!+
!
use berechne_mod
use errlist_mod
use do_replace_expr_mod
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
real(kind=PREC_DP), dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
!integer           , dimension(2)               , intent(in   ) :: index12        ! Index at dimension 1
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
real(kind=PREC_DP), dimension(2)               , intent(in   ) :: r_wert_limit     ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer :: i,j
integer :: length
real(kind=PREC_DP) :: wwert
!
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do j=lower_limit(2), upper_limit(2)
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         wwert  = berechne(line, length)                             ! Calculate current value of EXPR
         if(ier_num/=0) return
         if(llimited .and. (r_wert_limit(1)> (wwert) .or. (wwert)> r_wert_limit(2))) then
            ier_num = -50
            ier_typ = ER_FORT
         else
            user_variable(i, j) = wwert 
         endif
      enddo
      enddo
   else
      if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
         ier_num = -50
         ier_typ = ER_FORT
      else
         user_variable(lower_limit(1):upper_limit(1), lower_limit(2):upper_limit(2) ) = wert 
      endif
   endif
else 
   if(llimited .and. (r_wert_limit(1)> (wert) .or. (wert)> r_wert_limit(2))) then
      ier_num = -50
      ier_typ = ER_FORT
   else
      user_variable(lower_limit(1), lower_limit(2) ) = wert 
   endif
ENDIF 
!
end subroutine lib_set_para_real_D2_22
!
!*******************************************************************************
!
subroutine lib_set_para_char_D2(LOW, HIGH, user_variable, index_var, index_fix, &
           MAXW, lower_limit, upper_limit, lrange, lexpr, substr, cstring, line_expression)
!-
! Interface to specific char valued "lib_set_para" subroutines 
! Version for 2D character arrays, substitution along a single index
!+
!
use errlist_mod
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
character(len=*)  , dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index_var      ! Variable Index at dimension index_var
integer                                        , intent(in   ) :: index_fix      ! Index at fixed dimension
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
integer           , dimension(2)               , intent(in   ) :: substr         ! Indices for substring substitution
character(len=*)                               , intent(in   ) :: cstring        ! Fixed sting value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
!real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
!logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
!real(kind=PREC_DP), dimension(2)               , intent(in   ) :: r_wert_limit     ! Low and high value limit
!
if(substr(1) < 1 .or. substr(1)>substr(2) .or. substr(2)>len(user_variable)) then
   ier_num = -29
   ier_typ = ER_FORT
   return
endif
if(ubound(HIGH,1)==2) then
   if(index_var==1) then
continue
     call lib_set_para_char_D2_1(LOW, HIGH, user_variable, index_fix, MAXW, &
           lower_limit, upper_limit, lrange, lexpr, substr, cstring, line_expression)
!  elseif(index_var==2) then
!    call lib_set_para_real_D2_2(LOW, HIGH, user_variable, index_fix, MAXW, &
!          lower_limit, upper_limit, lrange, lexpr, wert, line_expression, llimited, r_wert_limit)
   endif
endif
!
end subroutine lib_set_para_char_D2
!
!*******************************************************************************
!
subroutine lib_set_para_char_D2_1(LOW, HIGH, user_variable, index_fix, MAXW, lower_limit, upper_limit, lrange, lexpr, substr, cstring, line_expression)
!-
! sets a variable to the individual value or a range
! Version for 2D char valued user variable, insert at first index
! user_variable[ insert_here, fixed ] = cstring
!+
!
use berechne_mod
use build_name_mod
use do_replace_expr_mod
use errlist_mod
use get_params_mod
use precision_mod
!
implicit none
!
integer           , dimension(2)               , intent(in   ) :: LOW            ! Dimensions for user variable
integer           , dimension(2)               , intent(in   ) :: HIGH           ! Dimensions for user variable
character(len=*)  , dimension(LOW(1):HIGH(1), LOW(2):HIGH(2)), intent(inout) :: user_variable
integer                                        , intent(in   ) :: index_fix      ! Index at fixed dimension 2
integer                                        , intent(in   ) :: MAXW           ! Dimension for limits
integer           , dimension(1:MAXW)          , intent(inout) :: lower_limit
integer           , dimension(1:MAXW)          , intent(inout) :: upper_limit
logical                                        , intent(in   ) :: lrange         ! User provided a range
logical                                        , intent(in   ) :: lexpr          ! User provided an expr as EXPR[]
integer           , dimension(2)               , intent(in   ) :: substr         ! Indices for substring substitution
character(len=*)                               , intent(in   ) :: cstring        ! Fixed sting value
character(len=*)                               , intent(in   ) :: line_expression ! line with expression
!real(kind=PREC_DP)                             , intent(in   ) :: wert           ! Fixed value
!logical                                        , intent(in   ) :: llimited                ! Parameter value "wert" needs to be limited
!integer           , dimension(2)               , intent(in   ) :: i_wert_limit     ! Low and high value limit
!
character(len=PREC_STRING) :: line
integer, parameter :: MAXWW = 20
character(len=PREC_STRING), dimension(MAXWW) :: cpara
integer                   , dimension(MAXWW) :: lpara
real(kind=PREC_DP)        , dimension(MAXWW) :: werte
integer :: i,ianz
integer :: length
!real(kind=PREC_DP) :: wwert
!
!write(*,*) ' IN lib_set_para_char_D2_1', MAXW
!write(*,*) ' VARIABLE DIMENSION', LOW, ' <> ', HIGH
!write(*,*) '             lower ', lower_limit, ' | ', maxw,  ' > ', index_fix
!write(*,*) '             upper ', upper_limit, ' | ', lrange, lexpr
!
!write(*,*) ' CHARACTER SUBST ', substr(1),substr(2), '>',cstring(1:len_trim(cstring)),'<'
!write(*,*) ' LINE_EXPRESSION ', line_expression(1:len_trim(line_expression))
if(lrange) then                                                      ! A range was specified
   if(lexpr) then                                                    ! An expression was used
      do i=lower_limit(1), upper_limit(1)
         line   = line_expression(1:len_trim(line_expression))
         length = len_trim(line)
         call do_replace_expr(line, length)                          ! needed for EXPR[1] + EXPR[2] or similar
         length = len_trim(line)
         call get_params (line, ianz, cpara, lpara, MAXWW, length)
         call do_build_name (ianz, cpara, lpara, werte, MAXWW, 1)
         line = cpara(1)
         user_variable(i, index_fix)(substr(1):substr(2)) = line(1:len_trim(line))
      enddo
   else
      user_variable(lower_limit(1):upper_limit(1), index_fix)(substr(1):substr(2)) = cstring(1:len_trim(cstring))
   endif
else 
   user_variable(lower_limit(1), index_fix)(substr(1):substr(2)) = cstring(1:len_trim(cstring))
ENDIF 
!
end subroutine lib_set_para_char_D2_1
!
!*******************************************************************************
!
logical function lib_check_dim_para_D1(MAXW, lower_limit, upper_limit, lrange, user_low, user_up)
!-
! Check if the limits: lower_limit upper_limit are inside user_* range
! Version for 1D arrays
!+
!
use errlist_mod
use lib_errlist_func,  only: no_error
!
implicit none
!
integer                   , intent(in   ) :: MAXW
integer, dimension(1:MAXW), INTENT(inout) :: lower_limit
integer, dimension(1:MAXW), intent(inout) :: upper_limit
logical                   , intent(in   ) :: lrange
integer                   , INTENT(in   ) :: user_low
integer                   , intent(in   ) :: user_up
!
lib_check_dim_para_D1 = .true.
call no_error
!
cond_dim: if(MAXW == 1) then
   if(lower_limit(1) == -huge(lower_limit)) lower_limit(1) = user_low   ! [:upper] => [user_low:upper]
   if(upper_limit(1) ==  huge(upper_limit)) upper_limit(1) = user_up    ! [lower:] => [lower:user_up ]
   cond_lower: if(lower_limit(1)>=user_low .and. lower_limit(1)<=user_up) then
      cond_range: if(lrange) then
         cond_upper: if(upper_limit(1)>=user_low .and. upper_limit(1)<=user_up .and. &
            lower_limit(1)<=upper_limit(1)                                      ) then
            lib_check_dim_para_D1 = .true.
         else cond_upper
            ier_msg(1) = 'Upper limit is outside'
            ier_num = -8 
            ier_typ = ER_FORT
            lib_check_dim_para_D1 = .false.
            exit cond_dim
         endif cond_upper
      else cond_range 
         lib_check_dim_para_D1 = .true.
      endif cond_range 
   else cond_lower
      ier_msg(1) = 'Lower limit is outside'
      ier_num = -8 
      ier_typ = ER_FORT
      lib_check_dim_para_D1 = .false.
      exit cond_dim
   endif cond_lower
else cond_dim
   ier_num = -13
   ier_typ = ER_FORT
   lib_check_dim_para_D1 = .false.
   exit cond_dim
endif cond_dim
!
end function lib_check_dim_para_D1
!
!*******************************************************************************
!
logical function lib_check_dim_para_D2(MAXW, lower_limit, upper_limit, lrange, user_low, user_up)
!-
! Check if the limits: lower_limit upper_limit are inside user_* range
! Version for 2D arrays
!+
!
use errlist_mod
!
implicit none
!
integer                   , intent(in   ) :: MAXW
integer, dimension(1:MAXW), INTENT(inout) :: lower_limit
integer, dimension(1:MAXW), intent(inout) :: upper_limit
logical                   , intent(in) :: lrange
integer, dimension(1:MAXW), INTENT(inout) :: user_low
integer, dimension(1:MAXW), intent(inout) :: user_up
!
integer :: i
!
lib_check_dim_para_D2 = .true.
!
cond_dim: if(MAXW >= 1) then
   loop_dim: do i=1, MAXW
   if(lower_limit(i) == -huge(lower_limit)) lower_limit(i) = user_low(i)! [:upper] => [user_low:upper]
   if(upper_limit(i) ==  huge(upper_limit)) upper_limit(i) = user_up(i) ! [lower:] => [lower:user_up ]
   cond_lower: if(lower_limit(i)>=user_low(i) .and. lower_limit(i)<=user_up(i)) then
      cond_range: if(lrange) then
         cond_upper: if(upper_limit(i)>=user_low(i) .and. upper_limit(i)<=user_up(i) .and. &
            lower_limit(i)<=upper_limit(i)                                      ) then
            lib_check_dim_para_D2 = .true.
         else cond_upper
            lib_check_dim_para_D2 = .false.
            exit cond_dim
         endif cond_upper
      else cond_range 
         lib_check_dim_para_D2 = .true.
      endif cond_range 
   else cond_lower
      lib_check_dim_para_D2 = .false.
      exit cond_dim
   endif cond_lower
   enddo loop_dim
else cond_dim
      ier_num = - 13
      ier_typ = ER_FORT
   lib_check_dim_para_D2 = .false.
   exit cond_dim
endif cond_dim
!
end function lib_check_dim_para_D2
!
!*******************************************************************************
!
END MODULE lib_upd_mod
