MODULE lib_upd_mod
!
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
!*****7*****************************************************************
SUBROUTINE    lib_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr) 
!-                                                                      
!       updates the parameter specified by ctype, index ww  to the      
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
INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
REAL(KIND=PREC_DP)        , INTENT(IN) :: wert 
CHARACTER (LEN=*),          INTENT(IN) :: cstring 
INTEGER, DIMENSION(2)     , INTENT(IN) :: substr ! Indices of substring
!
INTEGER :: i, ww2
!
search_var: DO i=var_sys+1, var_num
   IF(var_name(i) == ctype(1:len_trim(ctype))   ) THEN
      IF(var_entry(i)>0) THEN
         IF(0<ww(1) .AND. ww(1)<=var_field(var_entry(i))%var_shape(1) ) THEN
            IF(maxw==1) THEN
               IF(var_field(var_entry(i))%var_shape(2)>1) THEN
                  ier_num = -40
                  ier_typ = ER_FORT
                  RETURN
               ELSE
                  ww2 = 1
               ENDIF
            ELSEIF(maxw==2) THEN
               IF(0>=ww(2) .OR. ww(2)> var_field(var_entry(i))%var_shape(2)) THEN
                  ier_num = -40
                  ier_typ = ER_FORT
                  RETURN
               ELSE
                  ww2 = ww(2)
               ENDIF
            ELSE
               ier_num = -40
               ier_typ = ER_FORT
               RETURN
            ENDIF
               IF(var_type(i)==      IS_INTE) THEN
                  var_field(var_entry(i))%var_value(ww(1),ww2) = INT(wert)
               ELSEIF(var_type(i)==      IS_REAL) THEN
                  var_field(var_entry(i))%var_value(ww(1),ww2) = wert
               ELSEIF(var_type(i)==  IS_CHAR) THEN
                  var_field(var_entry(i))%var_char(ww(1),ww2)(substr(1):substr(2)) = cstring(1:len_trim(cstring))
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
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) THEN 
               inpara (ww (1) ) = int (wert) 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
ELSEIF (ctype.eq.'r') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) THEN 
               rpara (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
!                                                                       
ELSEIF (ctype.eq.'res') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_RES) THEN 
               res_para (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
ELSEIF (ctype.eq.'F_DERIV') THEN
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) THEN 
               kupl_deriv (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
ELSEIF (ctype.eq.'F_PARA') THEN
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) THEN 
               kupl_para (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
ELSEIF (ctype.eq.'ref_para') THEN
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_REF) THEN 
               ref_para (ww (1) ) = wert 
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
         WRITE (ier_msg (1), '(a)') ctype 
ENDIF 
!
END SUBROUTINE lib_upd_para                       
!
END MODULE lib_upd_mod
