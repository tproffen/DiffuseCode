MODULE ber_params_mod
!
CONTAINS
!
!*****7***************************************************************  
!
SUBROUTINE ber_params (ianz, cpara, lpara, werte, maxpara) 
!-                                                                      
!     Calculated the value of all expressions stored in cpara           
!+                                                                      
USE berechne_mod
!USE calc_expr_mod
USE errlist_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ianz
INTEGER, INTENT(IN) :: maxpara 
CHARACTER(LEN=*   ), DIMENSION(MAXPARA), INTENT(IN)  :: cpara
INTEGER            , DIMENSION(MAXPARA), INTENT(IN)  :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXPARA), INTENT(OUT) :: werte
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: line 
INTEGER             :: ll, i ,j
REAL(KIND=PREC_DP)  :: wert 
!                                                                       
!                                                                       
main: DO i = 1, ianz 
   ll   = lpara (i) 
   line = ' ' 
   line = '('//cpara (i) (1:ll) //')' 
   ll   = ll + 2 
   DO j=1, ll
      IF(IACHAR(line(j:j))==9) line(j:j) = ' '
   ENDDO
   wert = berechne (line, ll) 
   IF (ier_num /= 0) EXIT main
   werte (i) = wert 
ENDDO  main
!
END SUBROUTINE ber_params                     
!
!*****7***************************************************************  
!
SUBROUTINE ber_param  (ipara, cpara, lpara, werte, maxpara) 
!-                                                                      
!     Calculated the value of the expressions stored in cpara ipara
!+                                                                      
USE berechne_mod
!USE calc_expr_mod
USE errlist_mod 
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: ipara
INTEGER, INTENT(IN) :: maxpara 
CHARACTER(LEN=*   ), DIMENSION(MAXPARA), INTENT(IN)  :: cpara
INTEGER            , DIMENSION(MAXPARA), INTENT(IN)  :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXPARA), INTENT(OUT) :: werte
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(CPARA))) :: line 
INTEGER             :: ll, i ,j
REAL(KIND=PREC_DP)  :: wert 
!                                                                       
!                                                                       
i = ipara 
   ll   = lpara (i) 
   line = ' ' 
   line = '('//cpara (i) (1:ll) //')' 
   ll   = ll + 2 
   DO j=1, ll
      IF(IACHAR(line(j:j))==9) line(j:j) = ' '
   ENDDO
   wert = berechne (line, ll) 
   IF (ier_num /= 0) RETURN
   werte (i) = wert 
!
END SUBROUTINE ber_param                     
!
!*****7***************************************************************  
!
END MODULE ber_params_mod
