!*****7*****************************************************************
SUBROUTINE echo (zeile, lp) 
!-                                                                      
!     writes an echo to the screen                                      
!+                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxp 
      PARAMETER (maxp = 12) 
      CHARACTER ( * ) zeile 
      CHARACTER(1024) string 
      CHARACTER(1024) cpara (maxp) 
      CHARACTER(1024) cstr 
      INTEGER lpara (maxp) 
      INTEGER ianz 
      INTEGER lp, i, il 
      INTEGER iko, iqo, iqo2, lstring 
      REAL(KIND=PREC_DP), DIMENSION(MAXP) :: werte
!                                                                       
      INTEGER len_str 
      INTEGER socket_send 
!                                                                       
!     Find any "" which would mean that we need parameter substitution  
!                                                                       
      iqo = index (zeile (1:lp) , '"') 
      IF (iqo.eq.0) THEN 
!                                                                       
!     --None foud, harmless echo                                        
!                                                                       
         WRITE (output_io, 2010) zeile (1:lp) 
         WRITE (cstr, 2010) zeile (1:lp) 
         IF (lconn.and.lsocket) THEN 
            il = len_str (cstr) 
            ier_num = socket_send (s_conid, cstr, il) 
            IF(ier_num < 0) THEN
               ier_num = -19
               RETURN
            ELSE
               ier_num = 0
            ENDIF
         ENDIF 
!                                                                       
!        IF (output_status.eq.OUTPUT_FILE) THEN 
!           WRITE (output_io, 2010) zeile (1:lp) 
!        ENDIF 
      ELSE 
!                                                                       
!     --Look for matching ""                                            
!                                                                       
         iqo2 = index (zeile (iqo + 1:lp) , '"') 
         IF (iqo2.eq.0) THEN 
            ier_num = -  4 
            ier_typ = ER_COMM 
            RETURN 
         ENDIF 
!                                                                       
!     --save the string in quotation marks for later use as             
!       first parameter                                                 
!                                                                       
         iqo2 = iqo2 + iqo 
         iko = index (zeile (iqo2 + 1:lp) , ',') 
         IF (iko.ne.0) THEN 
            string = zeile (iqo2 + iko + 1:lp) 
            lstring = - (lp - iqo2 - iko) 
!                                                                       
!     --get all other parameters                                        
!                                                                       
            CALL get_params (string, ianz, cpara, lpara, maxp, lstring) 
         ELSE 
            string = ' ' 
            lstring = 1 
            ianz = 0 
         ENDIF 
!                                                                       
         IF (ier_num.eq.0) THEN 
            DO i = ianz, 1, - 1 
            cpara (i + 1) = cpara (i) 
            lpara (i + 1) = lpara (i) 
            ENDDO 
            cpara (1) = zeile (1:iqo2) 
            lpara (1) = iqo2 
            ianz = ianz + 1 
            CALL do_build_name (ianz, cpara, lpara, werte, maxp, 1) 
            IF (ier_num.eq.0) THEN 
!                                                                       
!     ------The first comma that follows the quotation marks must be    
!           omitted, all others retained as part of the echo            
!                                                                       
               IF (ianz.eq.1) THEN 
                  WRITE (output_io, 2000) cpara (ianz) (1:lpara (ianz) ) 
                  WRITE (cstr, 2000) cpara (ianz) (1:lpara (ianz) ) 
!                 IF (output_status.eq.OUTPUT_FILE) THEN 
!                    WRITE (output_io, 2000) cpara (ianz) (1:lpara (    &
!                    ianz) )                                            
!                 ENDIF 
               ELSEIF (ianz.eq.2) THEN 
                  WRITE (output_io, 2000) cpara (1) (1:lpara (1) ), cpara (    &
                  ianz) (1:lpara (ianz) )                               
                  WRITE (cstr, 2000) cpara (1) (1:lpara (1) ), cpara (  &
                  ianz) (1:lpara (ianz) )                               
!                 IF (output_status.eq.OUTPUT_FILE) THEN 
!                    WRITE (output_io, 2000) cpara (1) (1:lpara (1) ),  &
!                    cpara (ianz) (1:lpara (ianz) )                     
!                 ENDIF 
               ELSE 
                  WRITE (output_io, 2000) cpara (1) (1:lpara (1) ) , (cpara ( &
                  i) (1:lpara (i) ) , ',', i = 2, ianz - 1) , cpara (   &
                  ianz) (1:lpara (ianz) )                               
                  WRITE (cstr, 2000) cpara (1) (1:lpara (1) ) , (cpara (&
                  i) (1:lpara (i) ) , ',', i = 2, ianz - 1) , cpara (   &
                  ianz) (1:lpara (ianz) )                               
!                 IF (output_status.eq.OUTPUT_FILE) THEN 
!                    WRITE (output_io, 2000) cpara (1) (1:lpara (1) ) , &
!                    (cpara (i) (1:lpara (i) ) , ',', i = 2, ianz - 1) ,&
!                    cpara (ianz) (1:lpara (ianz) )                     
!                 ENDIF 
               ENDIF 
               IF (lconn.and.lsocket) THEN 
                  il = len_str (cstr) 
                  ier_num = socket_send (s_conid, cstr, il) 
                  IF(ier_num < 0) THEN
                     ier_num = -19
                     RETURN
                  ELSE
                     ier_num = 0
                  ENDIF
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
 2000 FORMAT    (1x,21a) 
 2010 FORMAT    (1x,  a) 
!                                                                       
      END SUBROUTINE echo                           
!
SUBROUTINE do_flush (zeile, lp) 
!-
! print simple message to * if error code is /= 0
!
USE errlist_mod
!
CHARACTER(LEN=*), INTENT(IN) :: zeile
INTEGER         , INTENT(IN) :: lp
!
WRITE(*,'(a)') zeile(1:MIN(lp,LEN_TRIM(zeile)))
IF(ier_num/=0) THEN
WRITE(*,'(a, 2i4)' ) ' ERROR Number /TYP ', ier_num, ier_typ
ENDIF
!
END SUBROUTINE do_flush
!*****7**************************************************************** 
