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
USE param_mod 
USE random_mod 
USE variable_mod
IMPLICIT none 
!                                                                       
INTEGER,                    INTENT(IN   ) :: ikl
INTEGER,                    INTENT(IN   ) :: iklz
CHARACTER (LEN=*),          INTENT(INOUT) :: string 
INTEGER,                    INTENT(INOUT) :: ll
INTEGER,                    INTENT(IN   ) :: maxw
REAL   , DIMENSION(1:maxw), INTENT(IN   ) :: ww
INTEGER,                    INTENT(IN   ) :: ianz
!                                                                       
CHARACTER(LEN=1024) :: zeile 
!                                                                       
INTEGER :: laenge, ltyp, kpara, kpara2
INTEGER :: i
INTEGER :: lcomm 
INTEGER :: length_com 
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
!
!  Test User defined variable arrays
!
success = .FALSE.
search_var: DO i=var_sys+1, var_num
!write(*,*) 'search ', string(ikl - lcomm:ikl - 1),ikl, lcomm, iklz,ianz,'|',string(1:50)
   IF(var_name(i) == string(ikl - lcomm:ikl - 1)) THEN
      IF(var_entry(i)>0) THEN
!write(*,*) ' kpara ', kpara, kpara2,var_field(var_entry(i))%var_shape(:), maxw
         IF(0<kpara .AND. kpara<=var_field(var_entry(i))%var_shape(1) ) THEN
            IF((ianz==2 .AND. 0<kpara2 .AND. kpara2<=var_field(var_entry(i))%var_shape(2)) .OR.  &
               (ianz==1 .AND. kpara2==0 .AND. var_field(var_entry(i))%var_shape(2)==1 ) )  THEN
               kpara2 = MAX(1, kpara2)
               zeile(1:ikl - lcomm-1) = string(1:ikl - lcomm-1)
               IF(var_type(i)==      IS_INTE) THEN
                  WRITE(zeile(ikl - lcomm:ikl + 13),'(i15)')                 &
                  NINT(var_field(var_entry(i))%var_value(kpara,kpara2))
               ELSEIF(var_type(i)==      IS_REAL) THEN
                  WRITE(zeile(ikl - lcomm:ikl + 13),'(e15.8e2)')             &
                  var_field(var_entry(i))%var_value(kpara,kpara2)
                  zeile (ikl + 11-lcomm:ikl + 11-lcomm) = 'e' 
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
   IF (lcomm.eq.1) THEN 
!                                                                       
      IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1: ikl - lcomm - 1)
         IF (string (ikl - 1:ikl - 1) .eq.'i') THEN 
            IF (ianz.eq.1) THEN 
               IF (0.le.kpara.and.kpara.le.MAXPAR) THEN 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') inpara (   &
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
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') rpara (&
                  kpara)                                                
                  zeile (ikl + 10:ikl + 10) = 'e' 
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
                  WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)')        &
                  res_para (kpara)                                      
                  zeile (ikl + 8:ikl + 8) = 'e' 
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
      WRITE (zeile (ikl - 4:ikl + 13) , '(i15    )') idum
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
                  WRITE(zeile(ikl - 6:ikl + 13), '(e15.8e2)') kupl_para(kpara)
                  zeile (ikl + 5:ikl + 5) = 'e'
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
                  WRITE(zeile(ikl - 7:ikl + 13), '(e15.8e2)') kupl_deriv(kpara)
                  zeile (ikl + 4:ikl + 4) = 'e'
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
                  WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)')        &
                  ref_para (kpara)
                  zeile (ikl + 3:ikl + 3) = 'e'
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
         ll = laenge+15 - ltyp - (iklz - ikl + 1) 
         IF (iklz + 1.le.laenge) zeile (ikl + 14:ll) = string (iklz + 1:&
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
SUBROUTINE    lib_upd_para (ctype, ww, maxw, wert, ianz) 
!-                                                                      
!       updates the parameter specified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
USE constants_mod
USE errlist_mod 
USE param_mod 
USE variable_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*),          INTENT(IN) :: ctype 
INTEGER,                    INTENT(IN) :: maxw
INTEGER,                    INTENT(IN) :: ianz 
INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
REAL   ,                    INTENT(IN) :: wert 
!
INTEGER :: i, ww2
!
search_var: DO i=var_sys+1, var_num
   IF(var_name(i) == ctype(1:len_trim(ctype))   ) THEN
      IF(var_entry(i)>0) THEN
         IF(0<ww(1) .AND. ww(1)<=var_field(var_entry(i))%var_shape(1) ) THEN
            IF(maxw==1 .AND. var_field(var_entry(i))%var_shape(2)>1) THEN
               ier_num = -40
               ier_typ = ER_FORT
               RETURN
            ELSE
               ww2 = 1
            ENDIF
            IF(maxw==2) THEN
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
