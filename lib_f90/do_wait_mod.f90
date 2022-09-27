MODULE do_wait_mod
!
LOGICAL :: wait_active = .TRUE.
!
private
!
public wait_active
public do_input
!
CONTAINS
!                                                                       
! Process the 'wait' command
!
!*****7***********************************************************      
SUBROUTINE do_input (zeile, lp) 
!                                                                       
USE ber_params_mod
use doact_mod
USE errlist_mod 
USE get_params_mod
USE lib_length
USE param_mod 
USE precision_mod
USE str_comp_mod
!
USE class_macro_internal
USE lib_macro_func
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 12
!                                                                       
CHARACTER ( LEN=* ), INTENT(INOUT)   ::  zeile 
INTEGER            , INTENT(INOUT)   ::  lp 
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
CHARACTER(LEN=PREC_STRING) :: line 
CHARACTER(LEN=4) :: cdummy 
REAL(KIND=PREC_DP), DIMENSION(maxw) ::  werte !(maxw) 
INTEGER :: i 
INTEGER :: lpara (maxw)
INTEGER :: ianz 
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
            IF(cdummy=='stop') THEN
!
               if(lblock) then
                  lblock = .false.
                  lblock_dbg = .false.
               endif
               if(lmakro) then
               WRITE( *, '(a,a1)') ' Macro stopped ', char (7)
               lmakro = .false.
               lmakro_error = .false.    ! Macro termination error off
               macro_level = 0
               CALL macro_close
               endif
!            line = '#'
!            il = 1
            ENDIF
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
                     res_para (0) = REAL(ianz) 
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
