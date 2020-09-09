MODULE do_show_mod
!
CONTAINS
!*****7**************************************************************** 
!
SUBROUTINE do_show_generic (cpara, lpara, MAXW) 
!-                                                                      
!     shows something related to the general command language           
!+                                                                      
USE do_variable_mod
USE errlist_mod 
USE str_comp_mod
!
IMPLICIT none 
!                                                                       
INTEGER                           , INTENT(IN) :: MAXW
CHARACTER (LEN=*), DIMENSION(MAXW), INTENT(IN) :: cpara
INTEGER          , DIMENSION(MAXW), INTENT(IN) :: lpara
!                                                                       
!                                                                       
IF (str_comp (cpara (1) , 'error', 2, lpara (1) , 5) ) THEN 
   CALL do_show_error 
ELSEIF (str_comp (cpara (1) , 'parallel', 2, lpara (1) , 8) ) THEN 
   CALL do_show_parallel 
!                                                                       
!     ----Show result array                'result'                     
!                                                                       
ELSEIF (str_comp (cpara (1) , 'res', 1, lpara (1) , 3) ) THEN 
   CALL do_show_res 
!                                                                       
!     ----Show variables                   'variables'                  
!                                                                       
ELSEIF (str_comp (cpara (1) , 'variables', 1, lpara (1) , 9) ) THEN
   CALL show_variables 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE do_show_generic                
!
!*****7*****************************************************************
!
SUBROUTINE do_show_res 
!-                                                                      
!     Shows the result array                                            
!+                                                                      
USE errlist_mod 
USE param_mod 
USE prompt_mod 
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER i, j, k, k1, k2 
!                                                                       
IF (NINT (res_para (0) ) == 0) THEN 
   WRITE (output_io, * ) 'Result array is empty' 
ELSE 
   WRITE (output_io, 2000) NINT (res_para (0) ) 
   j = NINT (res_para (0) ) / 5 
   IF (MOD (NINT (res_para (0) ), 5) == 0) THEN 
      j = j - 1 
   ENDIF 
   DO i = 0, j 
      k1 = 5 * i + 1 
      k2 = MIN (5 * i + 5, NINT (res_para (0) ) ) 
      WRITE (output_io, 2001) (res_para (k), k = k1, k2) 
   ENDDO 
ENDIF 
!                                                                       
 2000 FORMAT    ( i8,' Values in the result array') 
 2001 FORMAT    (5(2x,g14.7)) 
END SUBROUTINE do_show_res                    
!
!*****7*****************************************************************
!
SUBROUTINE do_show_error 
!-                                                                      
!     Shows the result array                                            
!+                                                                      
USE errlist_mod 
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
!------ - Error status setting                                          
!                                                                       
IF (ier_sta.eq.ER_S_CONT) THEN 
   WRITE (output_io, 2100) 
ELSEIF (ier_sta.eq.ER_S_LIVE) THEN 
   WRITE (output_io, 2102) 
ELSE 
   WRITE (output_io, 2105) 
ENDIF 
!                                                                       
 2100 FORMAT  (' Program continues after display of error message') 
 2102 FORMAT  (' Program lives on  after display of error message') 
 2105 FORMAT  (' Program terminates after display of error message') 
!
END SUBROUTINE do_show_error                  
!
!*****7*****************************************************************
!
SUBROUTINE do_show_parallel
!
!$ USE omp_lib
USE parallel_mod
USE prompt_mod
USE param_mod
!
IMPLICIT none 
!
INTEGER :: nthreads
INTEGER :: tid
!
nthreads = 1
IF(par_omp_use) THEN
!$OMP PARALLEL PRIVATE(tid)
!$ tid = OMP_GET_THREAD_NUM()
!$ IF(tid == 0) THEN
!$    IF(par_omp_maxthreads == -1) THEN
!$       nthreads = OMP_GET_NUM_THREADS()
!$    ELSE
!$       nthreads = MAX(1,MIN(par_omp_maxthreads, OMP_GET_NUM_THREADS()))
!$    ENDIF
!$ ENDIF
!$OMP END PARALLEL
   WRITE(output_io,'(a,i6,a)') 'OpenMP is active with ',nthreads, ' threads'
   res_para(0) = 2
   res_para(1) = 1
   res_para(2) = REAL(nthreads)
ELSE
   WRITE(output_io,'(a)') 'OpenMP is inactive'
   res_para(0) = 2
   res_para(1) = 0
   res_para(2) = 1
ENDIF
!                                                                       
END SUBROUTINE do_show_parallel
!
!*****7*****************************************************************
!
END MODULE do_show_mod
