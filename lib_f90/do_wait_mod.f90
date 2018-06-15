MODULE do_wait_mod
!
LOGICAL :: wait_active = .TRUE.
CONTAINS
!                                                                       
! Process the 'wait' command
!
!*****7***********************************************************      
SUBROUTINE do_input (zeile, lp) 
!                                                                       
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 12
!                                                                       
CHARACTER ( LEN=* ), INTENT(INOUT)   ::  zeile 
INTEGER            , INTENT(INOUT)   ::  lp 
!
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) line 
      CHARACTER(1) cdummy 
      REAL werte (maxw) 
      INTEGER i 
      INTEGER lpara (maxw)
      INTEGER ianz 
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
IF(wait_active) THEN
      lp = - lp 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (ianz.eq.0) THEN 
         WRITE ( *, 1000,advance='no') 
         READ ( *, 5000, err = 50, end = 50) cdummy 
      ELSEIF (ianz.eq.1.or.ianz.eq.2) THEN 
         IF (str_comp (cpara (1) , 'return', 1, lpara (1) , 6) ) THEN 
            WRITE ( *, 1000,advance='no') 
            READ ( *, 5000, err = 50, end = 50) cdummy 
         ELSEIF (str_comp (cpara (1) , 'input', 1, lpara (1) , 5) )     &
         THEN                                                           
            IF (ianz.eq.1) THEN 
               WRITE ( *, 1500,advance='no') 
            ELSEIF (ianz.eq.2) THEN 
               WRITE ( *, 1600,advance='no') cpara (2) (1:lpara (2) ) 
            ENDIF 
            READ ( *, 5000, err = 50, end = 50) line 
            lp = len_str (line) 
            CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.le.maxw) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     res_para (0) = float (ianz) 
                     DO i = 1, ianz 
                     res_para (i) = werte (i) 
                     ENDDO 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      RETURN 
!                                                                       
   50 CONTINUE 
      ier_num = - 9 
      ier_typ = ER_IO 
ENDIF
!                                                                       
 1000 FORMAT     (' ------ > Waiting for <RETURN> : ') 
 1500 FORMAT     (' ------ > Waiting for input : ') 
 1600 FORMAT     (a,' ') 
 5000 FORMAT     (a) 
!
END SUBROUTINE do_input                       
!
!*****7***********************************************************      
!
END MODULE do_wait_mod
