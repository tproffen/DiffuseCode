!
! routines specific for MAC_OS
!
LOGICAL FUNCTION IS_IOSTAT_END ( status_flag )
!
! Is actually FORTRAN2003 standard, should not be needed???
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: status_flag
!
IS_IOSTAT_END = status_flag == -1
!
END FUNCTION IS_IOSTAT_END
!*****7***********************************************************      
      SUBROUTINE do_fexist (zeile, lp, lout) 
!                                                                       
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
      INTEGER , INTENT(INOUT)::lp 
      LOGICAL , INTENT(IN)   ::lout
!
      CHARACTER(1024) cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz 
      LOGICAL lexist 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.1) then 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         INQUIRE (file = cpara (1) (1:lpara (1) ), exist = lexist) 
         IF (lexist) then 
            res_para (0) = 1 
            res_para (1) = 1 
            IF(lout) WRITE (output_io, 1000) cpara (1) (1:lpara (1) ) 
         ELSE 
            res_para (0) = 1 
            res_para (1) = 0 
            IF(lout) WRITE (output_io, 1100) cpara (1) (1:lpara (1) ) 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > File ',a,' exists ...') 
 1100 FORMAT     (' ------ > File ',a,' does NOT exist ...') 
      END SUBROUTINE do_fexist                      
!*****7***********************************************************      
      LOGICAL FUNCTION is_nan(x)
!
      REAL, INTENT(IN) :: x
!
      is_nan = isnan(x)
!
      END FUNCTION is_nan
!*****7***********************************************************      
