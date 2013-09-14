!*****7*****************************************************************
      SUBROUTINE do_scat (zeile, lp) 
!-                                                                      
!       Set scattering factors                                          
!+                                                                      
      USE errlist_mod 
!                                                                       
      USE config_mod 
      USE mixscat_mod 
!
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 11) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), length 
      INTEGER ianz, is, i, j, k 
      INTEGER jj (maxw) 
      REAL werte (maxw), ris 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
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
!                                                                       
      IF (ianz.eq.10) then 
         i = 1 
         CALL get_iscat (i, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) then 
            IF (werte (1) .gt.0) then 
               DO k = 1, i 
               jj (k) = nint (werte (1) ) 
               ENDDO 
               cpara (1) = '0.0' 
               lpara (1) = 3 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  DO k = 1, i 
                  DO i = 2, 9 
                  pdf_scat (i, jj (k), is) = werte (i) 
                  ENDDO 
                  pdf_scat (1, jj (k), is) = werte (10) 
                  pdf_scat_int (jj (k), is) = .false. 
                  ENDDO 
               ENDIF 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF (ianz.eq.2) then 
         i = 1 
         CALL get_iscat (i, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) then 
            IF (str_comp (cpara (2) , 'internal', 2, lpara (1) , 8) )   &
            then                                                        
               IF (werte (1) .gt.0) then 
                  k = nint (werte (1) ) 
                  pdf_scat_int (k, is) = .true. 
               ELSEIF (werte (1) .eq. - 1) then 
                  DO k = 0, elem_n 
                  pdf_scat_int (k, is) = .true. 
                  ENDDO 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      END SUBROUTINE do_scat                        
!*****7**************************************************************** 
      SUBROUTINE do_scale (zeile, lp) 
!-                                                                      
!       Set scale factor                                                
!+                                                                      
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), length 
      INTEGER ianz, id 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
      IF (ianz.eq.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         IF (werte (1) .gt.0.and.werte (1) .le.exp_nd) then 
            id = nint (werte (1) ) 
            exp_scal (id) = werte (2) 
         ELSE 
            ier_num = - 106 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_scale                       
!*****7**************************************************************** 
      SUBROUTINE do_match (zeile, lp) 
!-                                                                      
!       Get scales from fitting rho0 slopes                             
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
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
      INTEGER ianz, id, ip 
      REAL werte (maxw) 
      REAL rsl, drsl 
      REAL rmin, rmax, rho0 
!                                                                       
      IF (exp_nd.lt.2) then 
         ier_num = - 8 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ianz.eq.3) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         rmin = werte (1) 
         rmax = werte (2) 
         rho0 = werte (3) 
!                                                                       
         IF (rmin.ge.rmax) then 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ELSE 
            res_para (0) = 2.0 * exp_nd 
            ip = 1 
            DO id = 1, exp_nd 
            CALL get_slope (id, rmin, rmax, rsl, drsl) 
            exp_scal (id) = - 4.0 * 3.1415927 * rho0 / rsl 
            WRITE (output_io, 1000) id, rsl, drsl, exp_scal (id),       &
            r4 * 100.                                                   
            res_para (ip) = rsl 
            res_para (ip + 1) = drsl 
            ip = ip + 2 
            ENDDO 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (1x,'Data set ',i2,': Slope: ',f11.5,' +- ',f11.5,     &
     &               2x,'Scale: ',f11.5,2x,'(R=',f5.1,'%)')             
      END SUBROUTINE do_match                       
!*****7**************************************************************** 
      SUBROUTINE get_slope (id, rmin, rmax, sl, dsl) 
!-                                                                      
!       Gets slope of rho0 part                                         
!+                                                                      
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL rmin, rmax, r, sl, dsl 
      INTEGER i, id 
!                                                                       
      fit_np = 0 
      DO i = 1, exp_np (id) 
      r = exp_rmin (id) + (i - 1) * exp_deltar (id) 
      IF (r.ge.rmin.and.r.le.rmax) then 
         fit_np = fit_np + 1 
         fit_dat_x (fit_np) = r 
         fit_dat_y (fit_np) = exp_data (id, i) 
         fit_dat_w (fit_np) = exp_wic (id, i) 
      ENDIF 
      ENDDO 
!                                                                       
      npara = 2 
      ncycle = 100 
      urf = 0.5 
!                                                                       
      pinc (1) = 0. 
      pinc (2) = 1. 
      p (1) = 0.0 
      p (2) = 1.0 
!                                                                       
      CALL fit_poly 
!                                                                       
      sl = p (2) 
      dsl = dp (2) 
!                                                                       
      END SUBROUTINE get_slope                      
!*****7**************************************************************** 
      SUBROUTINE do_remove (zeile, lp) 
!-                                                                      
!       Sets pair to be removed                                         
!+                                                                      
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, i, j 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
      IF (ianz.eq.2) then 
         DO i = 1, ianz 
         elem_diff (i) = 0 
         CALL do_cap (cpara (i) ) 
         DO j = 1, elem_n 
         IF (cpara (i) .eq.elem_name (j) ) elem_diff (i) = j 
         ENDDO 
         IF (elem_diff (i) .eq.0) then 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
         ENDDO 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_remove                      
!*****7**************************************************************** 
      SUBROUTINE do_elements (zeile, lp) 
!-                                                                      
!       Sets sample composition                                         
!+                                                                      
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2 * MAXELEM) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), length 
      INTEGER ianz, ie 
      REAL werte (maxw) 
      REAL csum 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
      IF (ianz.ge.2) then 
         IF ( (ianz / 2) .gt.MAXELEM) then 
            ier_num = - 1 
            ier_typ = ER_APPL 
         ELSE 
            elem_n = 0 
            csum = 0.0 
            DO ie = 2, ianz, 2 
            elem_n = elem_n + 1 
            elem_name (elem_n) = cpara (ie-1) (1:lpara (ie-1) ) 
            CALL do_cap (elem_name (elem_n) ) 
            CALL ber_params (1, cpara (ie), lpara (ie), werte, maxw) 
            IF (ier_num.ne.0) return 
            elem_conc (elem_n) = werte (1) 
            csum = csum + werte (1) 
            ENDDO 
!                                                                       
            DO ie = 1, elem_n 
            elem_conc (ie) = elem_conc (ie) / csum 
            ENDDO 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_elements                    
!*****7*****************************************************************
      SUBROUTINE do_calc 
!-                                                                      
!       Calculates the differential PDF                                 
!+                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
      LOGICAL lgrid 
!                                                                       
      IF (elem_n.le.0) then 
         ier_num = - 5 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      IF (elem_diff (1) .le.0.or.elem_diff (2) .le.0) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL do_calc_weights 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (exp_nd.lt.2) then 
         ier_num = - 8 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      lgrid = .true. 
      DO i = 2, exp_nd 
      lgrid = lgrid.and. (exp_rmin (i) .eq.exp_rmin (1) ) 
      lgrid = lgrid.and. (exp_deltar (i) .eq.exp_deltar (1) ) 
      ENDDO 
      IF (.not.lgrid) then 
         ier_num = - 9 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      pdf_diff_rmin = exp_rmin (1) 
      pdf_diff_rmax = min (exp_rmax (1), exp_rmax (2) ) 
      pdf_diff_np = min (exp_np (1), exp_np (2) ) 
      pdf_diff_deltar = exp_deltar (1) 
!                                                                       
      IF (exp_w (1) .lt.exp_w (2) ) then 
         DO i = 1, pdf_diff_np 
         pdf_diff (i) = exp_scal (1) * exp_data (1, i) / exp_w (1)      &
         - exp_scal (2) * exp_data (2, i) / exp_w (2)                   
         ENDDO 
      ELSE 
         DO i = 1, pdf_diff_np 
         pdf_diff (i) = exp_scal (2) * exp_data (2, i) / exp_w (2)      &
         - exp_scal (1) * exp_data (1, i) / exp_w (1)                   
         ENDDO 
      ENDIF 
!                                                                       
      DO i = 1, pdf_diff_np 
      pdf_diff_err (i) = exp_scal (1) **2 * exp_data_err (1, i) **2 /   &
      exp_w (1) **2 + exp_scal (2) **2 * exp_data_err (2, i) **2 /      &
      exp_w (2) **2                                                     
      pdf_diff_err (i) = sqrt (pdf_diff_err (i) ) 
      ENDDO 
!                                                                       
      CALL do_show_error 
!                                                                       
      END SUBROUTINE do_calc                        
!*****7*****************************************************************
      SUBROUTINE do_calc_weights 
!-                                                                      
!       Calculated weights                                              
!+                                                                      
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL hh, form, bave 
      INTEGER is, ie, ia, ib 
!                                                                       
      hh = (pdf_xq / 2 / 3.1415927) **2 
      DO is = 1, exp_nd 
      CALL dlink (is) 
      bave = 0.0 
      DO ie = 1, elem_n 
      elem_w (ie, is) = form (ie, hh, is) 
      bave = bave+elem_conc (ie) * elem_w (ie, is) 
      ENDDO 
      pdf_bave2 (is) = bave**2 
      IF (elem_diff (1) .gt.0.and.elem_diff (2) .gt.0) then 
         exp_w (is) = elem_conc (elem_diff (1) ) * elem_conc (elem_diff &
         (2) ) * elem_w (elem_diff (1), is) * elem_w (elem_diff (2),    &
         is)                                                            
         exp_w (is) = exp_w (is) / pdf_bave2 (is) 
      ENDIF 
!                                                                       
      DO ia = 1, elem_n 
      DO ib = 1, elem_n 
      pdf_weights (ia, ib, is) = elem_conc (ia) * elem_conc (ib)        &
      * elem_w (ia, is) * elem_w (ib, is) / pdf_bave2 (is)              
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE do_calc_weights                
!*****7*********************************************************        
      SUBROUTINE get_iscat (ianz, cpara, lpara, werte, maxw) 
!-                                                                      
!     Determines the scattering type of the parameter                   
!+                                                                      
      USE charact_mod
      USE errlist_mod 
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(1024) zeile 
      INTEGER lpara (maxw) 
      INTEGER i, j, l, ianz, jj, jp 
      REAL werte (maxw) 
!                                                                       
      REAL berechne 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      DO i = 1, maxw 
      werte (i) = 0.0 
      ENDDO 
      jj = 1 
      jp = 0 
      j = 1 
      DO while (j.le.ianz.and.ier_num.eq.0) 
      i = ichar (cpara (j) (1:1) ) 
      IF (cpara (j) .eq.'all') then 
         werte (1) = - 1 
         RETURN 
      ELSEIF ( ( (a.le.i.and.i.le.z) .or. (aa.le.i.and.i.le.zz) ) .and. &
      (index (cpara (j) , '[') .eq.0) ) then                            
         CALL do_cap (cpara (j) ) 
         ier_num = - 27 
         ier_typ = ER_APPL 
         DO i = 0, elem_n 
         IF (cpara (j) .eq.elem_name (i) ) then 
            werte (jj) = i 
            jj = jj + 1 
            jp = jp + 1 
            ier_num = 0 
            ier_typ = ER_NONE 
         ENDIF 
         ENDDO 
      ELSE 
         zeile = ' ' 
         l = lpara (j) 
         zeile (1:1) = '(' 
         zeile (2:l + 1) = cpara (j) 
         zeile (l + 2:l + 2) = ')' 
         l = l + 2 
         werte (jj) = berechne (zeile, l) 
         IF (ier_num.eq.0) then 
            IF (0.le.nint (werte (jj) ) .and.nint (werte (jj) )         &
            .le.elem_n) then                                            
               jj = jj + 1 
               ier_num = 0 
               ier_typ = ER_NONE 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
      ENDIF 
      j = j + 1 
      ENDDO 
      ianz = max (ianz, jp) 
!                                                                       
      END SUBROUTINE get_iscat                      
