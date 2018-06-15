!*****7*****************************************************************
!     This file contains routines to save information.                  
!*****7*****************************************************************
      SUBROUTINE do_save (zeile, lp) 
!-                                                                      
!     Main menu for save command.                                       
!+                                                                      
      USE ber_params_mod
      USE get_params_mod
      USE build_name_mod
      USE errlist_mod 
      USE do_show_mod
      USE prompt_mod 
      USE times_mod 
!                                                                       
      USE config_mod 
      USE mixscat_mod 
      IMPLICIT none 
      include'date.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, i, ia, ib, io 
      REAL werte (maxw), r, cc 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
!---    'pdf' : save differential PDF                                   
!                                                                       
      IF (str_comp (cpara (1) , 'pdf', 2, lpara (1) , 3) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         CALL oeffne (57, cpara (1) , 'unknown') 
         IF (ier_num.ne.0) return 
         WRITE (output_io, 1500) cpara (1) (1:lpara (1) ) 
         WRITE (57, 5000) elem_name (elem_diff (1) ) (1:len_str (       &
         elem_name (elem_diff (1) ) ) ), elem_name (elem_diff (2) )     &
         (1:len_str (elem_name (elem_diff (2) ) ) )                     
         DO i = 1, pdf_diff_np 
         r = (i - 1) * pdf_diff_deltar + pdf_diff_rmin 
         WRITE (57, 5010) r, pdf_diff (i), pdf_diff_err (i) 
         ENDDO 
         CLOSE (57) 
!                                                                       
!---    'results' : save results to text file                           
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'result', 2, lpara (1) , 6) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         CALL oeffne (57, cpara (1) , 'unknown') 
         IF (ier_num.ne.0) return 
         WRITE (output_io, 1700) cpara (1) (1:lpara (1) ) 
         io = output_io 
         output_io = 57 
         CALL datum 
         WRITE (output_io, 7000) version, cdate, f_date 
         CALL do_show_all 
         WRITE (output_io, * ) 
         CALL do_show_error 
         output_io = io 
         CLOSE (57) 
!                                                                       
 7000 FORMAT      (1x,78('-'),/,1x,'Program MIXSCAT - Version ',a,/,    &
     &                1x,78('-'),/,1x,'Build date : ',a,/,              &
     &                1x,'Save date  : ',a,/,                           &
     &                1x,'Homepage   : http://discus.sourceforge.net',/,&
     &                1x,78('-'),/)                                     
!                                                                       
!---    'weights' : save weights for DISCUS PDF calculation             
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'weights', 2, lpara (1) , 7) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         CALL oeffne (57, cpara (1) , 'unknown') 
         IF (ier_num.ne.0) return 
         CALL do_calc_weights 
         WRITE (output_io, 1600) cpara (1) (1:lpara (1) ) 
         DO ia = 1, elem_n 
         DO ib = ia, elem_n 
         cc = 1.0 / elem_conc (ia) / elem_conc (ib) 
         IF (exp_w (1) .lt.exp_w (2) ) then 
            WRITE (57, 2000) elem_name (ia), elem_name (ib), (          &
            pdf_weights (ia, ib, 1) / exp_w (1) - pdf_weights (ia, ib,  &
            2) / exp_w (2) ) * cc                                       
         ELSE 
            WRITE (57, 2000) elem_name (ia), elem_name (ib), (          &
            pdf_weights (ia, ib, 2) / exp_w (2) - pdf_weights (ia, ib,  &
            1) / exp_w (1) ) * cc                                       
         ENDIF 
         ENDDO 
         ENDDO 
!                                                                       
         IF (exp_w (1) .lt.exp_w (2) ) then 
            WRITE (57, 2100) 1.0 / exp_w (1) - 1.0 / exp_w (2) 
         ELSE 
            WRITE (57, 2100) 1.0 / exp_w (2) - 1.0 / exp_w (1) 
         ENDIF 
!                                                                       
         CLOSE (57) 
!                                                                       
!-----  General show command                                            
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1500 FORMAT   (' Saving differential PDF to file : ',A,' ...') 
 1600 FORMAT   (' Saving weights of PDF to file : ',A,' ...') 
 1700 FORMAT   (' Saving results to file : ',A,' ...') 
 2000 FORMAT     ('set partial,',a,',',a,',',g20.8) 
 2100 FORMAT     ('set rden,',g20.8) 
 5000 FORMAT     ('#S 1 - Differential G(r): ',a,'-',a,' removed  ',/,  &
     &               '#L  r  diffG  s_diffG')                           
 5010 FORMAT   (F9.4,2(3X,F21.10)) 
      END SUBROUTINE do_save                        
