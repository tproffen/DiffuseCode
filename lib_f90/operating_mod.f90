!*****7*****************************************************************
SUBROUTINE do_operating (zeile, lp) 
!-                                                                      
!     writes an echo to the screen                                      
!+                                                                      
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE prompt_mod 
USE sys_compiler
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!                                                                       
INTEGER, PARAMETER :: maxp = 12
!
CHARACTER(LEN=PREC_STRING) :: string , line
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXP) :: werte
!
INTEGER :: ianz 
INTEGER :: lline 
INTEGER :: i 
INTEGER :: iko, iqo, iqo2, lstring 
!                                                                       
!     Find any "" which would mean that we need parameter substitution  
!                                                                       
iqo = index (zeile (1:lp) , '"') 
IF (iqo == 0) THEN 
!                                                                       
!     --None foud, harmless echo                                        
!                                                                       
   CALL do_operating_comm (zeile (1:lp) ) 
ELSE 
!                                                                       
!     --Look for matching ""                                            
!                                                                       
   iqo2 = index (zeile (iqo + 1:lp) , '"') 
   IF (iqo2 == 0) THEN 
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
   IF (iko /= 0) THEN 
      string = zeile (iqo2 + iko + 1:lp) 
      lstring = - (lp - iqo2 - iko) 
   ELSE 
      string = ' ' 
      lstring = 1 
   ENDIF 
!
   line = ' '
   lline = 0
   IF(iqo>1) THEN
      line = zeile(1:iqo-1)
      lline = iqo-1
   ENDIF
!                                                                       
!     --get all other parameters                                        
!                                                                       
   CALL get_params (string, ianz, cpara, lpara, maxp, lstring) 
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
         IF (ianz == 1) THEN 
            string = cpara (ianz) (1:lpara (ianz) ) 
            lstring = lpara (ianz) 
         ELSEIF (ianz == 2) THEN 
            string = cpara (1) (1:lpara (1) ) //cpara (ianz)      &
            (1:lpara (ianz) )                                     
            lstring = lpara (1) + lpara (2) 
         ELSE 
            string = cpara (1) (1:lpara (1) ) 
            lstring = lpara (1) 
            DO i = 1, ianz 
               zeile = string (1:lstring) //cpara (i) (1:lpara (i) ) 
               lstring = lstring + lpara (i) 
               string = zeile 
            ENDDO 
         ENDIF 
         IF(lline>0) string=line(1:lline)//string(1:LEN_TRIM(string))
         CALL do_operating_comm (string) 
      ENDIF 
   ENDIF 
ENDIF 
!                                                                       
END SUBROUTINE do_operating                   
!                                                                       
!*****7***********************************************************      
