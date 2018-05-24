MODULE calc_expr_mod
!
CONTAINS
!
!****7***************************************************************** 
!
SUBROUTINE do_math (line, indxg, length) 
!-                                                                      
!     Calculates the value of an expression and stores the result       
!     in the proper variable                                            
!                                                                       
      USE berechne_mod
USE ber_params_mod
USE constants_mod
USE do_variable_mod
      USE errlist_mod 
      USE do_string_alloc_mod
      USE get_params_mod
      USE set_sub_generic_mod
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 10
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(IN   ) :: indxg
INTEGER          , INTENT(INOUT) :: length
!
CHARACTER(LEN=1024)   :: zeile, cpara (maxw) 
INTEGER, DIMENSION(3) :: var_is_type
!                                                                       
      INTEGER lpara (maxw) 
      INTEGER i, ikk, iii (maxw), ianz, lll 
INTEGER :: indxb
!                                                                       
      REAL wert, werte (maxw) 
!                                                                       
!     String substitution???                                            
!                                                                       
lll = indxg -1
indxb = INDEX(line(1:lll),'[')
IF(indxb >0) lll = indxb - 1

CALL p_get_var_type(line, lll, var_is_type)
IF(var_is_type(1)== IS_UNKNOWN) THEN     ! unknown variable name on left side
   ier_num = -2
   ier_typ = ER_FORT
   ier_msg(1) = 'Offending name is '//line(1:lll)
   RETURN
ENDIF
IF(var_is_type(3)==IS_READ) THEN
   ier_num = -41
   ier_typ = ER_FORT
   ier_msg(1) = 'Offending name is '//line(1:lll)
   RETURN
ENDIF
!
IF (INDEX (line, '"') .gt.0.or.INDEX (line, '''') .gt.0) then 
   IF(var_is_type(1)==IS_CHAR) THEN
      CALL do_string_alloc (line, indxg, length) 
      RETURN 
   ELSE
      ier_num = -42
      ier_typ = ER_FORT
      ier_msg(1) = 'Offending name is '//line(1:lll)
      RETURN 
   ENDIF 
ENDIF 
!                                                                       
!     Get the expression                                                
!                                                                       
lll = length - (indxg + 1) + 1 
CALL get_params (line (indxg + 1:length), ianz, cpara, lpara, maxw, lll)
IF (ier_num.ne.0) then 
   RETURN 
ELSEIF (ianz.eq.0) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
zeile = '('//cpara(1)(1:lpara(1)) //')' 
i     = lpara(1) + 2 
!                                                                       
!     Calculate the expression                                          
!                                                                       
IF(var_is_type(1)/=IS_CHAR) THEN     ! Variable on left side is numeric
      wert = berechne (zeile, i) 
      IF (ier_num.eq.0) then 
!                                                                       
!-----evaluate the index of the variable                                
!                                                                       
         lll = indxg - 1 
         CALL get_params (line (1:indxg - 1), ianz, cpara, lpara, maxw, lll)
         IF (ier_num.eq.0) then 
            line = cpara (1) 
            i = lpara (1) 
            ikk = INDEX (line, '[') 
            IF (ikk.lt.i.and.ikk.gt.0) then 
               IF (line (i:i) .eq.']') then 
                  IF (i.gt.ikk + 1) then 
                     zeile = ' ' 
                     zeile (1:i - ikk - 1) = line (ikk + 1:i - 1) 
                     lll = i - ikk - 1 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw, lll)
                     IF (ier_num.eq.0) then 
                        IF (ianz.ge.1.or.ianz.le.3) then 
                           CALL ber_params (ianz, cpara, lpara, werte, maxw)
                           IF (ier_num.eq.0) then 
                              DO i = 1, ianz 
                              iii (i) = nint (werte (i) ) 
                              ENDDO 
!                                                                       
!     ------------Store result in the variable                          
!                                                                       
                              CALL p_upd_para (line (1:ikk - 1), iii,ianz, wert, ianz)
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_FORT 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 10 
                     ier_typ = ER_FORT 
                  ENDIF 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSEIF (ikk.eq.0) then 
               CALL upd_variable (line (1:i), i, wert, cpara (1), lpara (1) )
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
      ENDIF 
ELSE             ! Variable on left side is CHARACTER
         zeile = line (1:indxg) //'"%c",'//line (indxg + 1:length) 
         length = length + 5 
         CALL do_string_alloc (zeile, indxg, length) 
ENDIF 
!                                                                       
END SUBROUTINE do_math                        
!
!*****7**************************************************************** 
!
END MODULE calc_expr_mod
