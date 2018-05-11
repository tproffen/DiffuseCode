MODULE do_eval_mod
!
CONTAINS
!
!****7***************************************************************** 
!
SUBROUTINE do_eval (line, i) 
!-                                                                      
!     evaluates the expression stored in line                           
!                                                                       
USE ber_params_mod
      USE errlist_mod 
      USE do_string_alloc_mod
      USE get_params_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, PARAMETER :: maxw = 10
!                                                                       
      CHARACTER(LEN=1024), INTENT(INOUT) :: line
      INTEGER            , INTENT(INOUT) :: i
!
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) cstr 
      INTEGER lpara (maxw) 
      INTEGER ianz, il 
      INTEGER length 
      INTEGER :: indxg = 0
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (line.eq.' ') then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ELSE 
!                                                                       
!     String substitution???                                            
!                                                                       
      IF (INDEX (line, '"') .gt.0.or.INDEX (line, '''') .gt.0 ) THEN
         CALL do_string_alloc (line, indxg, i) 
         WRITE (output_io, 3000) line(1:LEN_TRIM(line))
         RETURN 
      ENDIF 
!                                                                       
!     --Non blank line. Make length negative to avoid removing blanks   
!                                                                       
         i = - i 
         CALL get_params (line, ianz, cpara, lpara, maxw, i) 
         IF (ier_num.eq.0) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               length = lpara (1) 
               DO i = 1, ianz 
               length = max (length, lpara (i) ) 
               ENDDO 
               DO i = 1, ianz 
!              WRITE ( *, 2222) cpara (i) (1:length), werte (i) 
!              IF (output_status.eq.OUTPUT_FILE) then 
                  WRITE (output_io, 2222) cpara (i) (1:length), werte ( i)
!              ENDIF 
               IF (lconn.and.lsocket.and.i.eq.1) then 
                  WRITE (cstr, 2222) cpara (i) (1:lpara (i) ), werte (i) 
                  il = len_str (cstr) 
                  CALL socket_send (s_conid, cstr, il) 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
 2222 FORMAT    (' Value of ',a,' = ',g15.8) 
3000  FORMAT(' Value is ',a)
!
END SUBROUTINE do_eval                        
!
!****7***************************************************************** 
!
END MODULE do_eval_mod
