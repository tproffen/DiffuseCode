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
      CHARACTER(LEN=1024) :: zeile, cpara (maxw) 
!                                                                       
      INTEGER lpara (maxw) 
      INTEGER i, ikk, iii (maxw), ianz, lll 
!                                                                       
      REAL wert, werte (maxw) 
!                                                                       
!     String substitution???                                            
!                                                                       
      IF (INDEX (line, '"') .gt.0.or.INDEX (line, '''') .gt.0) then 
         CALL do_string_alloc (line, indxg, length) 
         RETURN 
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
      i = lpara (1) 
      zeile = '('//cpara (1) (1:i) //')' 
      i = i + 2 
!                                                                       
!     Calculate the expression                                          
!                                                                       
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
      ELSE 
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
