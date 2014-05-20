!*****7*****************************************************************
!     This file contains routines to display information.               
!*****7*****************************************************************
      SUBROUTINE mixscat_do_show (zeile, lp) 
!-                                                                      
!     Main menu for show command.                                       
!+                                                                      
      USE errlist_mod 
      USE config_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), length 
      INTEGER ianz 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
      IF (ianz.eq.0) then 
         CALL do_show_all 
         RETURN 
      ENDIF 
!                                                                       
!------ - 'config' : show various limits ...                            
!                                                                       
      IF (str_comp (cpara (1) , 'config', 4, lpara (1) , 6) ) then 
         CALL do_show_config 
!                                                                       
!------ - 'analysis' : show error analysis ...                          
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'analysis', 4, lpara (1) , 8) )     &
      then                                                              
         CALL do_show_analysis 
!                                                                       
!------ - 'scat' : show scattering curve                                
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'scat', 1, lpara (1) , 4) ) then 
         CALL do_show_scat (ianz, cpara, lpara, werte, maxw) 
!                                                                       
!-----  - General show command                                          
!                                                                       
      ELSE 
         CALL do_show_generic (cpara, lpara, maxw) 
      ENDIF 
!                                                                       
      END SUBROUTINE mixscat_do_show                        
!********************************************************************** 
      SUBROUTINE do_show_config 
!-                                                                      
!     Shows information about PDFFIT configuration                      
!+                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      WRITE (output_io, 1000) 
      WRITE (output_io, 2030) MAXDSET 
      WRITE (output_io, 2040) MAXDAT 
      WRITE (output_io, 2050) MAXELEM 
      WRITE (output_io, * ) 
!                                                                       
 1000 FORMAT     (' Current MIXSCAT configuration :',/) 
 2030 FORMAT     ('   Maximum number of data sets              : ',i8) 
 2040 FORMAT     ('   Maximum number of data points per file   : ',i8) 
 2050 FORMAT     ('   Maximum number of elements               : ',i8) 
!                                                                       
      END SUBROUTINE do_show_config                 
!*****7*****************************************************************
      SUBROUTINE do_show_analysis 
!+                                                                      
!     Show information about errors                                     
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL sum_d (MAXDSET), err_d (MAXDSET) 
      REAL sum_diff, err_diff 
      INTEGER ip, is, i 
      INTEGER len_str 
!                                                                       
      IF (pdf_diff_np.eq.0) then 
         ier_num = - 11 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      sum_diff = 0.0 
      err_diff = 0.0 
      DO is = 1, exp_nd 
      sum_d (is) = 0.0 
      err_d (is) = 0.0 
      ENDDO 
!                                                                       
      IF (exp_w (1) .lt.exp_w (2) ) then 
         WRITE (output_io, 1000) (elem_name (elem_diff (i) ) (1:len_str &
         (elem_name (elem_diff (i) ) ) ), i = 1, 2), 1, 2               
      ELSE 
         WRITE (output_io, 1000) (elem_name (elem_diff (i) ) (1:len_str &
         (elem_name (elem_diff (i) ) ) ), i = 1, 2), 2, 1               
      ENDIF 
!                                                                       
      DO ip = 1, pdf_diff_np 
      sum_diff = sum_diff + abs (pdf_diff (ip) ) 
      err_diff = err_diff + pdf_diff_err (ip) 
      DO is = 1, exp_nd 
      sum_d (is) = sum_d (is) + abs (exp_data (is, ip) ) 
      err_d (is) = err_d (is) + exp_data_err (is, ip) 
      ENDDO 
      ENDDO 
!                                                                       
      WRITE (output_io, 2000) 
      DO is = 1, exp_nd 
      WRITE (output_io, 2100) is, exp_name (is) (1:len_str (exp_name (  &
      is) ) ), sum_d (is) / float (pdf_diff_np), err_d (is) / float (   &
      pdf_diff_np), 100. * err_d (is) / sum_d (is)                      
      ENDDO 
      WRITE (output_io, 3000) 
      WRITE (output_io, 2200) 'Extracted differential', sum_diff /      &
      float (pdf_diff_np) , err_diff / float (pdf_diff_np) , 100. *     &
      err_diff / sum_diff                                               
      WRITE (output_io, 3000) 
!                                                                       
 1000 FORMAT (1x,78('-'),/,                                             &
     &        1x,'ERROR ANALYSIS - Removed contribution : ',a4,         &
     &           ' - ',a4,8x,'Difference: ',i2,' - ',i2,/,1x,78('-'))   
 2000 FORMAT (1x,'Data',23x,'File',6x,'Average data',6x,                &
     &           'Average sigma',2x,'Relative')                         
 2100 FORMAT (3x,i2,2x,a25,1x,f17.8,2x,f17.8,3x,f6.2,'%') 
 2200 FORMAT (7x,a25,1x,f17.8,2x,f17.8,3x,f6.2,'%') 
 3000 FORMAT (1x,78('-')) 
      END SUBROUTINE do_show_analysis               
!********************************************************************** 
      SUBROUTINE do_show_all 
!+                                                                      
!     Show information                                                  
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CALL output_pdf 
      CALL output_elem 
      CALL output_weights 
!                                                                       
      END SUBROUTINE do_show_all                    
!*****7*****************************************************************
      SUBROUTINE output_elem 
!+                                                                      
!     Show information about sample settings                            
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ie, id, i 
      INTEGER len_str 
!                                                                       
      CALL do_calc_weights 
!                                                                       
      WRITE (output_io, 1000) (elem_name (elem_diff (i) ) (1:len_str (  &
      elem_name (elem_diff (i) ) ) ), i = 1, 2)                         
      WRITE (output_io, 1050) (i, i = 1, 2) 
      WRITE (output_io, 1100) (exp_w (i), i = 1, 2) 
      WRITE (output_io, 1200) (pdf_bave2 (i), i = 1, 2) 
!                                                                       
      DO ie = 1, elem_n 
      WRITE (output_io, 1300) elem_name (ie), (elem_w (ie, i), i = 1, 2)&
      , elem_conc (ie)                                                  
      ENDDO 
      WRITE (output_io, 2000) 
!                                                                       
 1000 FORMAT     (1x,78('-'),/,                                         &
     &               1x,'SETUP - Removing contribution : ',a,' - ',a,/, &
     &               1x,78('-'))                                        
 1050 FORMAT     (1x,'Data set    :',8x,i18,10x,i18) 
 1100 FORMAT     (1x,'Weight w_mn :',8x,f18.8,10x,f18.8) 
 1200 FORMAT     (1x,'<b>^2       :',8x,f18.8,10x,f18.8) 
 1300 FORMAT     (1x,'b(',a4,')     :',8x,f18.8,10x,f18.8,4x,'c=',f5.3) 
 2000 FORMAT     (1x,78('-')) 
!                                                                       
      END SUBROUTINE output_elem                    
!*****7*****************************************************************
      SUBROUTINE output_weights 
!+                                                                      
!     Show information about weights in resulting differential PDF      
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ia, ib, i 
      INTEGER len_str 
!                                                                       
      CALL do_calc_weights 
      WRITE (output_io, 1000) 
!                                                                       
      DO ia = 1, elem_n 
      DO ib = ia, elem_n 
      IF (exp_w (1) .lt.exp_w (2) ) then 
         WRITE (output_io, 2000) elem_name (ia), elem_name (ib),        &
         pdf_weights (ia, ib, 1) / exp_w (1) - pdf_weights (ia, ib, 2)  &
         / exp_w (2)                                                    
      ELSE 
         WRITE (output_io, 2000) elem_name (ia), elem_name (ib),        &
         pdf_weights (ia, ib, 2) / exp_w (2) - pdf_weights (ia, ib, 1)  &
         / exp_w (1)                                                    
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      IF (exp_w (1) .lt.exp_w (2) ) then 
         WRITE (output_io, 2100) 1.0 / exp_w (1) - 1.0 / exp_w (2) 
      ELSE 
         WRITE (output_io, 2100) 1.0 / exp_w (2) - 1.0 / exp_w (1) 
      ENDIF 
      WRITE (output_io, 3000) 
!                                                                       
 1000 FORMAT     (1x,'DIFFERENTIAL WEIGHTS',/,1x,78('-')) 
 2000 FORMAT     (1x,'w(',a4,'-',a4,'):',12x,f14.8) 
 2100 FORMAT     (1x,'w(rho_0)    :',12x,f14.8) 
 3000 FORMAT     (1x,78('-')) 
!                                                                       
      END SUBROUTINE output_weights                 
!*****7*****************************************************************
      SUBROUTINE output_pdf 
!+                                                                      
!     Show information about data                                       
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(25) rad (2) 
      INTEGER i, j, k, lc, ll, id, iset 
      INTEGER len_str 
!                                                                       
      DO iset = 1, exp_nd 
      IF (exp_lxray (iset) ) then 
         WRITE (rad (iset) , 2110) 'x-rays', pdf_xq 
      ELSE 
         WRITE (rad (iset) , 2115) 'Neutrons' 
      ENDIF 
      ENDDO 
!                                                                       
      WRITE (output_io, 1000) 
      WRITE (output_io, 1100) (i, i = 1, exp_nd) 
      WRITE (output_io, 1200) (exp_name (i) (1:len_str (exp_name (i) ) )&
      , i = 1, exp_nd)                                                  
      WRITE (output_io, 1300) (rad (i), i = 1, exp_nd) 
      WRITE (output_io, 1400) (exp_rmin (i), i = 1, exp_nd) 
      WRITE (output_io, 1500) (exp_rmax (i), i = 1, exp_nd) 
      WRITE (output_io, 1600) (exp_deltar (i), i = 1, exp_nd) 
      WRITE (output_io, 1700) (exp_qmax (i), i = 1, exp_nd) 
      WRITE (output_io, 1800) (exp_scal (i), i = 1, exp_nd) 
!                                                                       
 1000 FORMAT     (1x,78('-'),/,1x,'DATA SETS',/,1x,78('-')) 
 1100 FORMAT     (1x,'Data set    : ',20x,i5,23x,i5) 
 1200 FORMAT     (1x,'Filename    : ',a25,3x,a25) 
 1300 FORMAT     (1x,'Radiation   : ',a25,3x,a25) 
 1400 FORMAT     (1x,'Rmin (A)    : ',18x,f7.3,21x,f7.3) 
 1500 FORMAT     (1x,'Rmax (A)    : ',18x,f7.3,21x,f7.3) 
 1600 FORMAT     (1x,'Dr (A)      : ',18x,f7.3,21x,f7.3) 
 1700 FORMAT     (1x,'Qmax (A^-1) : ',18x,f7.3,21x,f7.3) 
 1800 FORMAT     (1x,'Scale       : ',11x,f14.8,14x,f14.8) 
 2110 FORMAT     (a12,1x,'(',f5.1,' A^-1)') 
 2115 FORMAT     (a25) 
!                                                                       
      END SUBROUTINE output_pdf                     
!*****7*****************************************************************
      SUBROUTINE do_show_scat (ianz, cpara, lpara, werte, maxw) 
!-                                                                      
!     Shows the scattering curve of an atom type.                       
!+                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      INTEGER ianz 
      INTEGER lpara (maxw) 
      REAL werte (maxw), ris 
!                                                                       
      INTEGER i, j, is 
!                                                                       
      CHARACTER(9) at_name_d 
      CHARACTER(9) at_name 
!                                                                       
      LOGICAL latom (0:MAXELEM) 
!                                                                       
      IF (ianz.lt.3) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL ber_params (1, cpara (1), lpara (1), ris, 1) 
      IF (ier_num.ne.0) return 
!                                                                       
      is = nint (ris) 
      IF (is.le.0.or.is.gt.exp_nd) then 
         ier_num = - 106 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      CALL dlink (is) 
!                                                                       
      DO i = 0, elem_n 
      latom (i) = .false. 
      ENDDO 
      CALL get_iscat (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.eq.0) then 
         IF (werte (1) .eq. - 1) then 
            DO i = 0, elem_n 
            latom (i) = .true. 
            ENDDO 
         ELSE 
            DO i = 1, ianz 
            latom (nint (werte (i) ) ) = .true. 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      WRITE (output_io, 3000) is 
!                                                                       
      DO i = 1, elem_n 
      IF (latom (i) ) then 
         at_name_d = elem_name (i) 
         WRITE (output_io, 3010) at_name_d, (pdf_scat (2 * j, i, is),   &
         j = 1, 4), (pdf_scat (2 * j + 1, i, is), j = 1, 4), pdf_scat ( &
         1, i, is), pdf_scat_int (i, is)                                
      ENDIF 
      ENDDO 
!                                                                       
 3000 FORMAT (1x,78('-'),/,1x,'SCATTERING LENGTHS - data set ',i2,/,    &
     &        1x,78('-'))                                               
 3010 FORMAT(1x,a9,'a(i)      : ',4(f12.7,3x)/                          &
     &        10x,'b(i)      : ',4(f12.7,3x)/                           &
     &        10x,'c         : ',  f12.7,29x,'Lookup table : ',L1,/,    &
     &        1x,78('-'))                                               
!                                                                       
      END SUBROUTINE do_show_scat                   
