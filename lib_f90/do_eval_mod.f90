MODULE do_eval_mod
!
CONTAINS
!
!****7***************************************************************** 
!
SUBROUTINE do_eval (line, i, lout) 
!-                                                                      
!     evaluates the expression stored in line                           
!                                                                       
USE blanks_mod
USE ber_params_mod
      USE errlist_mod 
      USE do_string_alloc_mod
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, PARAMETER :: MAXW = 10
!                                                                       
      CHARACTER(LEN=1024), INTENT(INOUT) :: line
      INTEGER(KIND=PREC_INT_WORD) , INTENT(INOUT) :: i
      LOGICAL            , INTENT(IN   ) :: lout
!
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) cstr , zeile
CHARACTER(LEN=32) :: form_s
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
WRITE(form_s,'(A,I2.2,a,I2.2,A)')  '('' Value of '',a,'' = '',g',PREC_WIDTH,'.',PREC_DIGIT,')'
!                                                                       
!     String substitution???                                            
!                                                                       
      CALL rem_insig_bl(line,i)
!     IF (INDEX (line, '"') .gt.0.or.INDEX (line, '''') .gt.0 ) THEN
      IF (INDEX (line, '"')  == 1.or.INDEX (line, '''') == 1  ) THEN
         CALL do_string_alloc (line, indxg, i) 
         IF(lout) WRITE (output_io, 3000) line(1:LEN_TRIM(line))
         i = LEN_TRIM(line)
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
                  IF(lout) WRITE (output_io, form_s) cpara (i) (1:length), werte ( i)
!              ENDIF 
               IF (lconn.and.lsocket.and.i.eq.1) then 
                  IF(lout) WRITE (cstr, form_s) cpara (i) (1:lpara (i) ), werte (i) 
                  il = len_str (cstr) 
                  CALL socket_send (s_conid, cstr, il) 
               ENDIF 
               ENDDO 
               WRITE(line,'(G18.10E4)')  werte(1)
               i = LEN_TRIM(line)
            ELSEIF(ier_num==-1 .AND. ier_typ == ER_FORT) THEN   ! Try string substitution
               length = LEN_TRIM(line)
               indxg  = 0
               zeile  = line (1:indxg) //'"%c",'//line (indxg + 1:length)
               length = LEN_TRIM(zeile)
               CALL do_string_alloc (zeile, indxg, length)
               IF(lout) WRITE (output_io, 3000) zeile(1:LEN_TRIM(zeile))
               line = zeile
               i = LEN_TRIM(line)
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
3000  FORMAT(' Value is ',a)
!
END SUBROUTINE do_eval                        
!
!****7***************************************************************** 
!
END MODULE do_eval_mod
