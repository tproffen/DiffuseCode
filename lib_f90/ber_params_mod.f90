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
CHARACTER(LEN=1024), DIMENSION(MAXPARA), INTENT(IN)  :: cpara
INTEGER            , DIMENSION(MAXPARA), INTENT(IN)  :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXPARA), INTENT(OUT) :: werte
!                                                                       
CHARACTER(LEN=1024) :: line 
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
END MODULE ber_params_mod
