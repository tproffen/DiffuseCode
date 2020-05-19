!*****7**************************************************************** 
SUBROUTINE start_learn (zeile, lcomm) 
!-                                                                      
!     These routines control if the commands typed in are               
!     within a learning ==ence or not ..                              
!                                                                       
!-                                                                      
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE learn_mod 
USE precision_mod
USE prompt_mod 
USE sys_compiler
!
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER :: maxw = 12 
!                                                                       
CHARACTER (LEN= * ), INTENT(INOUT) :: zeile 
INTEGER            , INTENT(INOUT) :: lcomm 
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) :: cpara
INTEGER,             DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP),  DIMENSION(MAXW) :: werte
INTEGER :: ip, ianz, len_str 
!                                                                       
IF (llearn) THEN 
   ier_num = - 7 
   ier_typ = ER_IO 
ELSE 
!                                                                       
   CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
   IF (ier_num==0) THEN 
      IF (ianz==1) THEN 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num==0) THEN 
            fname = cpara (1) 
         ENDIF 
      ELSEIF (ianz==0) THEN 
         fname = pname (1:len_str (pname) ) //'.mac' 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      IF (ier_num==0) THEN 
         ip = INDEX (fname, '.') 
         IF (ip==0) THEN 
            ip = INDEX (fname, ' ') 
            fname = fname (1:ip - 1) //'.mac' 
         ENDIF 
         CALL oeffne (33, fname, 'unknown') 
         IF (ier_num /= 0) return 
         llearn = .true. 
         WRITE (output_io, 1000) fname (1:LEN_TRIM (fname) ) 
      ENDIF 
   ENDIF 
ENDIF 
!                                                                       
1000 FORMAT    (' ------ > Learning started, makrofile: ',a) 
!
END SUBROUTINE start_learn                    
!
!*****7**************************************************************** 
!
SUBROUTINE ende_learn 
!+                                                                      
!     End of learning ==enze                                          
!-                                                                      
USE errlist_mod 
USE learn_mod 
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
IF (.NOT.llearn) THEN 
   ier_num = -8 
   ier_typ = ER_IO 
ELSE 
   llearn = .false. 
   CLOSE (33) 
   WRITE (output_io, 1000) fname (1:LEN_TRIM (fname) ) 
ENDIF 
!
1000 FORMAT    (' ------ > Learning finished, makrofile: ',a) 
!                                                                       
END SUBROUTINE ende_learn                     
