!*****7*****************************************************************
      SUBROUTINE do_read (zeile, lp) 
!-                                                                      
!       Main menu for read command.                                     
!+                                                                      
      USE errlist_mod 
!                                                                       
      USE config_mod 
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
      INTEGER ianz 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                                                                       
      IF (ianz.eq.0) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ELSE 
!                                                                       
!-----    Reading a PDF data file                                       
!                                                                       
         IF (str_comp (cpara (1) , 'data', 2, lpara (1) , 4) ) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_read_data (ianz, cpara, lpara, werte, maxw, .false.) 
!                                                                       
!-----    Unknown command                                               
!                                                                       
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_read                        
!*****7*****************************************************************
      SUBROUTINE do_read_data (ianz, cpara, lpara, werte, maxw) 
!+                                                                      
!     Reads observed PDF as xy ASCII file.                              
!-                                                                      
      USE debug_mod
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
!                                                                       
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) pfile 
      CHARACTER(80) line 
      INTEGER lpara (maxw) 
      INTEGER lpfile 
      INTEGER ianz, iianz, i, ip, iset 
      REAL werte (maxw) 
      REAL ra, re, dr, da, de 
      LOGICAL lwei 
!                                                                       
!------ Check if there is space for another data set                    
!                                                                       
      IF (exp_nd.eq.MAXDSET) then 
         ier_num = - 104 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      iset = exp_nd+1 
!                                                                       
!------ Now analyse given parameters                                    
!                                                                       
      IF (ianz.ge.2) then 
!                                                                       
!------ - get exp. method (neutron/x-ray)                               
!                                                                       
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:1) .eq.'N') then 
            exp_lxray (iset) = .false. 
         ELSEIF (cpara (1) (1:1) .eq.'X') then 
            exp_lxray (iset) = .true. 
         ELSE 
            ier_num = - 105 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
!------ - get filename                                                  
!                                                                       
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         pfile = cpara (1) 
         lpfile = lpara (1) 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!------ Read observed PDF for given plane                               
!                                                                       
      CALL oeffne (17, pfile, 'old', .false.) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Has the file a prepended History ?                              
!------ Get information from history and store in res[i]                
!                                                                       
!------ res[1] = Temperature                                            
!                                                                       
      res_para (0) = 0 
!                                                                       
      READ (17, 1000) line 
      IF (line (1:7) .eq.'History') then 
         WRITE (output_io, 3000) 
   22    CONTINUE 
         READ (17, 1000, end = 999) line 
         CALL extract_key (line, 'temp=') 
         CALL extract_key (line, 'Qmax=') 
         IF (res_para (res_para (0) ) >0) then 
            exp_qmax (iset) = res_para (res_para (0) ) 
         ENDIF 
         IF (line (1:16) .eq.'##### start data') goto 33 
         GOTO 22 
      ELSE 
         BACKSPACE (17) 
      ENDIF 
   33 CONTINUE 
!                                                                       
!------ Any other header lines starting with # ?                        
!                                                                       
   44 CONTINUE 
      READ (17, 1000, end = 999) line 
      IF (line (1:1) .eq.'#') goto 44 
      BACKSPACE (17) 
!                                                                       
!------ Finally we actually read the data                               
!                                                                       
      READ (17, 1000, end = 20, err = 999) line 
      CALL count_col (line, ianz) 
      IF (ianz.lt.2) goto 999 
      BACKSPACE (17) 
!                                                                       
      ip = 1 
      lwei = .true. 
!                                                                       
   10 CONTINUE 
      READ (17, *, end = 20, err = 999) (werte (i), i = 1, ianz) 
      IF (ip.eq.1) ra = werte (1) 
      re = werte (1) 
      exp_data (iset, ip) = werte (2) 
      IF (ianz.eq.3.and.werte (3) .gt.0.0) then 
         exp_data_err (iset, ip) = werte (3) 
         exp_wic (iset, ip) = 1.0 / werte (3) **2 
      ELSEIF (ianz.ge.4.and.werte (4) .gt.0.0) then 
         exp_data_err (iset, ip) = werte (4) 
         exp_wic (iset, ip) = 1.0 / werte (4) **2 
      ELSE 
         exp_data_err (iset, ip) = 0.0 
         exp_wic (iset, ip) = 1.0 
         lwei = .false. 
      ENDIF 
      ip = ip + 1 
      IF (ip.gt.MAXDAT) goto 9998 
      GOTO 10 
   20 CONTINUE 
      CLOSE (17) 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
      exp_nd = iset 
!                                                                       
      exp_rmin (iset) = ra 
      exp_rmax (iset) = re 
      exp_np (iset) = ip - 1 
      exp_deltar (iset) = (re-ra) / float (exp_np (iset) - 1) 
      exp_name (iset) = pfile (1:lpfile) 
!                                                                       
      WRITE (output_io, 9100) iset, ra, re, exp_np (iset) 
      IF (.not.lwei) write (output_io, 9200) 
      RETURN 
!                                                                       
  999 CONTINUE 
      CLOSE (17) 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 9998 CONTINUE 
      CLOSE (17) 
      ier_num = - 111 
      ier_typ = ER_APPL 
      RETURN 
!                                                                       
 9999 CONTINUE 
      CLOSE (17) 
      ier_num = - 115 
      ier_typ = ER_APPL 
      RETURN 
!                                                                       
 1000 FORMAT     (a) 
 3000 FORMAT     (' History information found ...') 
 9100 FORMAT     (' Read PDF data set ',I3,'  (r = ',F7.3,' to ',F7.3,  &
     &               ' A, ',I5,' points) ...')                          
 9200 FORMAT     (' No sigmas for G(r) found, using unit weights ..') 
      END SUBROUTINE do_read_data                   
!*****7*****************************************************************
      SUBROUTINE extract_key (line, key) 
!+                                                                      
!     Gets numbers from history part of data file                       
!-                                                                      
      USE param_mod 
      USE config_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) line, key 
      INTEGER is, ie, ll, lk 
!                                                                       
      INTEGER len_str 
!                                                                       
      ll = len_str (line) 
      lk = len_str (key) 
!                                                                       
      is = index (line, key (1:lk) ) 
      IF (is.ne.0) then 
         is = is + lk 
         DO while (line (is:is) .eq.' ') 
         is = is + 1 
         ENDDO 
         ie = index (line (is:ll) , ' ') 
         IF (ie.eq.0) ie = ll 
         res_para (0) = res_para (0) + 1 
         READ (line (is:is + ie), * ) res_para (res_para (0) ) 
      ENDIF 
!                                                                       
      END SUBROUTINE extract_key                    
!*****7*****************************************************************
      SUBROUTINE do_xray (zeile, lp) 
!-                                                                      
!     Setting Q to calculate f(xray) ..                                 
!+                                                                      
      USE errlist_mod 
      USE prompt_mod 
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
      INTEGER ianz 
      REAL werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         WRITE (output_io, 1000) pdf_xq 
      ELSEIF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         pdf_xq = werte (1) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (1x,'------ > Current Q for X-ray weights : ',g14.6) 
      END SUBROUTINE do_xray                        
