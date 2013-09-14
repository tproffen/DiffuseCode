      SUBROUTINE initarrays 
!+                                                                      
!       Startvalues for all important arrays.                           
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
!                                                                       
      USE config_mod 
      USE mixscat_mod 
!
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k 
!                                                                       
!------ Sockt defaults                                                  
!                                                                       
      s_port = 3334 
      s_ipallowed = "localhost" 
!                                                                       
!                                                                       
      exp_nd = 0 
      elem_n = 0 
      pdf_diff_np = 0 
      pdf_xq = 0.0 
!                                                                       
      elem_diff (1) = 0 
      elem_diff (2) = 0 
      elem_name (0) = 'UNKN' 
!                                                                       
      DO i = 1, MAXDSET 
      exp_scal (i) = 1.0 
      exp_np (i) = 0 
      pdf_bave2 (i) = 0.0 
!                                                                       
      DO j = 1, MAXELEM 
      pdf_scat_int (j, i) = .true. 
      ENDDO 
      ENDDO 
!                                                                       
!     /errors/                                                   
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      ier_sta = ER_S_CONT 
      DO i = 1, 3 
      ier_msg (i) = ' ' 
      ENDDO 
!                                                                       
      END SUBROUTINE initarrays                     
!*****7*****************************************************************
      SUBROUTINE autodef 
!-                                                                      
!     Tries to open a default file for the integer and real variables   
!+                                                                      
      USE envir_mod 
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
!
      IMPLICIT NONE
!
!                                                                       
      INTEGER idef 
!                                                                       
      DATA idef / 34 / 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      CALL open_def (idef) 
      IF (ier_num.eq.0) then 
         ier_num = - 3 
         ier_typ = ER_IO 
         READ (idef, *, end = 998, err = 998) inpara 
         READ (idef, *, end = 998, err = 998) rpara 
      ENDIF 
      ier_num = 0 
      ier_typ = ER_NONE 
  998 CONTINUE 
      IF (ier_num.ne.0) then 
         CALL errlist 
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
      CLOSE (idef) 
!                                                                       
      WRITE (output_io, * ) 
 1000 FORMAT    (a) 
      END SUBROUTINE autodef                        
