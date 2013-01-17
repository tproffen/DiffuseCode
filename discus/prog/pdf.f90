!*****7*****************************************************************
      SUBROUTINE pdf 
!-                                                                      
!     This sublevel contains all routines dealing with the              
!     PDF anlysis part of DISCUS. This segment uses variables           
!     of the RMC and CHEM segment which are simply overwritten !!       
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE modify_mod
      USE pdf_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'macro.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'param.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 24  ! A command requires at least these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER(1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      REAL           , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
      REAL           , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: wwerte
      REAL           , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: wwwerte
      INTEGER        , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
!
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile, cdummy 
      INTEGER lbeg (3) 
      INTEGER lp, length
      INTEGER indxg, ianz, lbef, i, ia, is, ic (3), iianz, jjanz 
      LOGICAL lout, ldummy 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
!
      IF( cr_nscat > PDF_MAXSCAT .or. MAXSCAT > PDF_MAXSCAT) THEN
                 pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
                 pdf_ndat  = MAX(pdf_ndat ,           PDF_MAXDAT)
                 pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
                 CALL alloc_pdf( pdf_nscat, pdf_ndat, pdf_nbnd )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      CALL no_error 
!                                                                       
   10 CONTINUE 
!                                                                       
      prom = prompt (1:len_str (prompt) ) //'/pdf' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
            CALL do_math (line, indxg, length) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            CALL file_kdo (line (2:length), length - 1) 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
            CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression 'eval'                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
            CALL do_eval (zeile, lp) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.  &
              str_comp (befehl, '?   ', 1, lbef, 4) ) then                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 11 
               CALL do_hel ('discus pdf '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
            cdummy = ' ' 
            IF (zeile.ne.' ') then 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ Waiting for user input                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
            CALL do_input (zeile, lp) 
!                                                                       
!------ ------------------------------------------------------------    
!------ Here start the PDF specific commands ...                        
!------ ------------------------------------------------------------    
!                                                                       
!------ Just calculate PDF for current parameters                       
!                                                                       
         ELSEIF (str_comp (befehl, 'calc', 2, lbef, 4) ) then 
            CALL pdf_setup 
            IF (ier_num.eq.0) then 
               pdf_skal = 1.0 / rmc_skal (1) 
               CALL pdf_determine (.true.) 
            ENDIF 
!                                                                       
!------ Read observed PDF from XY file (ASCII)                          
!                                                                       
         ELSEIF (str_comp (befehl, 'data', 2, lbef, 4) ) then 
            CALL pdf_readdata (zeile, lp) 
!                                                                       
!------ Runs PDF fit                                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
            CALL pdf_setup 
            IF (ier_num.eq.0) then 
               CALL pdf_run 
            ENDIF 
!                                                                       
!------ Save structure or PDF                                           
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
            CALL pdf_save (zeile, lp) 
!                                                                       
!------ Set PDF parameters                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
            CALL pdf_set (zeile, lp) 
!                                                                       
!------ Reset PDF segement                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'rese', 2, lbef, 4) ) then 
            pdf_ldata = .false. 
!                                                                       
!------ Show current PDF parameters                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.0) then 
                  CALL pdf_show ('PDF') 
               ELSEIF (ianz.eq.1) then 
                  CALL do_cap (cpara (1) ) 
                  CALL pdf_show (cpara (1) (1:3) ) 
               ELSE 
                  ier_typ = ER_COMM 
                  ier_num = - 6 
               ENDIF 
            ENDIF 
!                                                                       
!------ selecting/deselecting atoms for PDF calculation                 
!                                                                       
         ELSEIF (str_comp (befehl, 'isele', 3, lbef, 5) .or.        &
                 str_comp (befehl, 'idese', 3, lbef, 5) ) then                            
!                                                                       
            CALL atom_select (zeile, lp, 0, MAXSCAT, pdf_allowed_i,    &
            ldummy, .false., str_comp (befehl,             &
            'isele', 3, lbef, 5) )                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'jsele', 3, lbef, 5) .or.        &
                 str_comp (befehl, 'jdese', 3, lbef, 5) ) then                            
!                                                                       
            CALL atom_select (zeile, lp, 0, MAXSCAT, pdf_allowed_j, &
            ldummy, .false., str_comp (befehl,             &
            'jsele', 3, lbef, 5) )                                      
!                                                                       
!                                                                       
!------ selecting/deselecting atoms                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'sele', 3, lbef, 4) .or.         &
                 str_comp (befehl, 'dese', 2, lbef, 4) ) then                             
!                                                                       
            CALL atom_select (zeile, lp, 0, MAXSCAT, rmc_allowed, &
            rmc_sel_atom, .false., str_comp (              &
            befehl, 'sele', 3, lbef, 4) )                               
!                                                                       
!------ selecting/deselecting of molecules                              
!                                                                       
         ELSEIF (str_comp (befehl, 'msel', 2, lbef, 4) .or.   &
                 str_comp (befehl, 'mdes', 2, lbef, 4) ) then                             
!                                                                       
            CALL mole_select (zeile, lp, 0, MAXSCAT, rmc_allowed, &
            rmc_sel_atom, .false., str_comp (  &
            befehl, 'msel', 2, lbef, 4) )                               
!                                                                       
!------ no command found                                                
!                                                                       
         ELSE 
            ier_num = - 8 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
!------ any errors ?                                                    
!                                                                       
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
      GOTO 10 
!                                                                       
 9999 CONTINUE 
!                                                                       
      END SUBROUTINE pdf                            
!*****7*****************************************************************
      SUBROUTINE pdf_setup 
!+                                                                      
!     Setup for various arrays and functions for PDF calculation.       
!-                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'debug.inc' 
      include'errlist.inc' 
      include'wink.inc' 
!                                                                       
      REAL sincut, rcut, z, bave, hh, fac, rtot, ract 
      INTEGER :: max_bnd
      INTEGER i, j, ia, is, js, nn, nnn 
      LOGICAL ltot 
!                                                                       
      REAL form, quad 
!                                                                       
      WRITE (output_io, 1000) 
!                                                                       
!------ Do we have a crystal ?                                          
!                                                                       
      IF (cr_natoms.lt.1) then 
         ier_num = - 4 
         ier_typ = ER_PDF 
         RETURN 
      ENDIF 
!
!------ Allocate arrays
!
      sincut   = 0.025
      rcut     = 1.0 / (pdf_qmax * sincut)
      nn       = INT (rcut / pdf_deltar) + 100
      nnn      = INT ( (pdf_rfmax + rcut) / pdf_deltar)
      pdf_ndat = MAX ( nn, nnn) + 1
      max_bnd = 1
      DO i=1,3
         max_bnd = MAX(max_bnd, INT(pdf_rmax / cr_a0 (i)))
      ENDDO
      pdf_nbnd  = MAX ( cr_icc (1), cr_icc (2), cr_icc (3) ) + 1 + max_bnd
!
      IF(cr_nscat > PDF_MAXSCAT .or. pdf_ndat > PDF_MAXDAT .or.         &
         pdf_nbnd > PDF_MAXBND                                  ) THEN
        pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
        pdf_ndat  = MAX(pdf_ndat ,           PDF_MAXDAT)
        pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
        CALL alloc_pdf( pdf_nscat, pdf_ndat, pdf_nbnd )
      ENDIF
!                                                                       
!------ Setting up array for periodic boundaries                        
!                                                                       
      IF (cr_icc (1) .gt.PDF_MAXBND.or.cr_icc (2) .gt.PDF_MAXBND.or. &
          cr_icc (3) .gt.PDF_MAXBND) then                                                  
         ier_num = - 3 
         ier_typ = ER_PDF 
         RETURN 
      ENDIF 
!                                                                       
      DO j = 1, 3 
      DO i = - PDF_MAXBND, PDF_MAXBND+PDF_MAXBND 
      pdf_bnd (j, i) = i 
      DO while (pdf_bnd (j, i) .le.0) 
      pdf_bnd (j, i) = pdf_bnd (j, i) + cr_icc (j) 
      ENDDO 
      DO while (pdf_bnd (j, i) .gt.cr_icc (j) ) 
      pdf_bnd (j, i) = pdf_bnd (j, i) - cr_icc (j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!------ Setting up weighting (b(i)b(j)/<b**2>)                          
!                                                                       
      CALL dlink (pdf_lxray, ano, lambda, rlambda) 
      bave = 0.0 
      hh = pdf_xq**2 
!                                                                       
      cr_n_real_atoms = 0 
      DO ia = 1, cr_natoms 
      IF (cr_at_lis (cr_iscat (ia) ) .ne.'VOID') then 
         bave = bave+form (cr_iscat (ia), cr_scat, pdf_lxray, hh) 
         cr_n_real_atoms = cr_n_real_atoms + 1 
      ENDIF 
      ENDDO 
      bave = bave / float (cr_n_real_atoms) 
!                                                                       
      IF (.not.pdf_lweights) then 
         DO is = 0, cr_nscat 
         DO js = 0, cr_nscat 
         pdf_weight (is, js) = form (is, cr_scat, pdf_lxray, hh)        &
         * form (js, cr_scat, pdf_lxray, hh) / bave**2                  
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
!------ Get the ratio total pairs/selected weight in structure          
!                                                                       
      ltot = .true. 
      DO is = 1, cr_nscat 
      ltot = ltot.and.pdf_allowed_i (is) .and.pdf_allowed_j (is) 
      ENDDO 
!                                                                       
      IF (.not.ltot) then 
         rtot = 0.0 
         ract = 0.0 
!                                                                       
         DO i = 1, cr_natoms 
         DO j = 1, cr_natoms 
         is = cr_iscat (i) 
         js = cr_iscat (j) 
         rtot = rtot + pdf_weight (is, js) 
         IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or.  &
              (pdf_allowed_j (is) .and.pdf_allowed_i (js) ) )     &
                 ract = ract + pdf_weight (is, js)                                            
         ENDDO 
         ENDDO 
         pdf_dnorm = ract / rtot 
      ELSE 
         pdf_dnorm = 1.0 
      ENDIF 
!                                                                       
!------ Setting up SINC function for convolution with PDF               
!                                                                       
      IF (pdf_qmax.gt.0.0) then 
         sincut = 0.025 
         rcut = 1.0 / (pdf_qmax * sincut) 
         z = pdf_deltar * pdf_qmax 
         nn = int (rcut / pdf_deltar) + 100 
!                                                                       
         IF ( (pdf_rfmax + rcut) .gt.pdf_rmax) then 
            nnn = int ( (pdf_rfmax + rcut) / pdf_deltar) 
            IF (nnn.le.PDF_MAXDAT) then 
               pdf_rmax = pdf_rfmax + rcut 
               pdf_bin = int (pdf_rmax / pdf_deltar) 
               WRITE (output_io, 1100) pdf_rmax 
            ELSE 
               ier_num = - 2 
               ier_typ = ER_PDF 
               RETURN 
            ENDIF 
         ENDIF 
!                                                                       
         DO i = 1, nn 
         pdf_sinc (i) = sin (z * float (i) ) / (pdf_deltar * float (i) ) 
         ENDDO 
         DO i = nn + 1, 2 * PDF_MAXDAT 
         pdf_sinc (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Setting up PDF segment ...') 
 1100 FORMAT     (' Extending PDF search distance to ',F8.4,' A ...') 
      END SUBROUTINE pdf_setup                      
!*****7*****************************************************************
      SUBROUTINE pdf_show (cmd) 
!+                                                                      
!     Shows current parameters                                          
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(9) at_name, at_lis (MAXSCAT+1) 
      CHARACTER ( * ) cmd 
      CHARACTER(9) cpoly (5) 
      INTEGER i, j, k 
!                                                                       
      DATA cpoly / 'linear  ', 'square  ', '3. order', '4. order', &
                   '5. order' /                                                           
!                                                                       
      IF (cmd.eq.'ALL'.or.cmd.eq.'PDF') then 
         CALL pdf_setup 
         IF (ier_num.ne.0) return 
!                                                                       
         pdf_skal = 1.0 / rmc_skal (1) 
!                                                                       
         WRITE (output_io, 1000) 
         WRITE (output_io, 2000) pdf_rmax, pdf_deltar, pdf_bin 
!                                                                       
         IF (pdf_lxray) then 
            WRITE (output_io, 2010) 'X-Rays', pdf_xq 
         ELSE 
            WRITE (output_io, 2015) 'Neutrons' 
         ENDIF 
!                                                                       
         WRITE (output_io, 2100) 
         IF (pdf_qmax.eq.0.0) then 
            WRITE (output_io, 2110) 
         ELSE 
            WRITE (output_io, 2120) pdf_qmax 
         ENDIF 
         IF (pdf_sigmaq.eq.0.0) then 
            WRITE (output_io, 2210) 
         ELSE 
            WRITE (output_io, 2220) pdf_sigmaq 
         ENDIF 
         IF (chem_period (1) .and.chem_period (2) .and.        &
             chem_period (3) )                            then
            IF (pdf_2d) then 
               WRITE (output_io, 2300) 'applied - 2D' 
            ELSE 
               WRITE (output_io, 2300) 'applied - 3D' 
            ENDIF 
            IF (pdf_finite.eq.PDF_BACK_PERIOD) then 
               WRITE (output_io, 2310) 
            ELSEIF (pdf_finite.eq.PDF_BACK_SPHERE) then 
               WRITE (output_io, 2320) pdf_sphere 
            ELSE 
               WRITE (output_io, 2300) 'illegal' 
            ENDIF 
         ELSE 
            WRITE (output_io, 2300) 'not applied' 
            IF (pdf_finite.eq.PDF_BACK_PERIOD) then 
               WRITE (output_io, 2350) 
            ELSEIF (pdf_finite.eq.PDF_BACK_SPHERE) then 
               WRITE (output_io, 2360) pdf_sphere 
            ELSEIF (pdf_finite.eq.PDF_BACK_POLY) then 
               WRITE (output_io, 2370) pdf_diam_poly
               DO i = 1, pdf_poly_n 
               WRITE (output_io, 2371) cpoly (i), pdf_poly (i) 
               ENDDO 
            ELSEIF (pdf_finite.eq.PDF_BACK_TANH) then 
               WRITE (output_io, 2380) pdf_diam, pdf_shape 
            ENDIF 
         ENDIF 
         WRITE (output_io, 2390) pdf_scale 
         IF (pdf_gauss) then 
            WRITE (output_io, 2400) 'applied' 
            WRITE (output_io, 2405) pdf_qalp 
            WRITE (output_io, 2410) pdf_delta 
            WRITE (output_io, 2415) pdf_gamma 
            WRITE (output_io, 2420) pdf_srat, pdf_rcut 
         ELSE 
            WRITE (output_io, 2400) 'not applied' 
            WRITE (output_io, 2405) pdf_qalp 
         ENDIF 
         IF (pdf_lrho0) then 
            WRITE (output_io, 2443) 
         ELSE 
            IF (pdf_lrho0_rel) then 
               WRITE (output_io, 2445) pdf_rho0, 'relative' 
            ELSE 
               WRITE (output_io, 2445) pdf_rho0, 'absolute' 
            ENDIF 
         ENDIF 
         WRITE (output_io, 2450) pdf_dnorm 
         IF (pdf_lexact) then 
            WRITE (output_io, 2455) 'all atoms' 
         ELSE 
            WRITE (output_io, 2455) 'neighboring unit cells' 
         ENDIF 
!                                                                       
         WRITE (output_io, 3000) 
         WRITE (output_io, 3010) pdf_rfmin, int (pdf_rfmin / pdf_deltar) 
         WRITE (output_io, 3020) pdf_rfmax, int (pdf_rfmax / pdf_deltar) 
         WRITE (output_io, 3030) pdf_skal, rmc_doskal 
!                                                                       
         j = 0 
         DO i = 0, cr_nscat 
         IF (pdf_allowed_i (i) ) then 
            j = j + 1 
            at_lis (j) = at_name (i) 
         ENDIF 
         ENDDO 
         WRITE (output_io, 3100) (at_lis (k), k = 1, j) 
                                                                        
         j = 0 
         DO i = 0, cr_nscat 
         IF (pdf_allowed_j (i) ) then 
            j = j + 1 
            at_lis (j) = at_name (i) 
         ENDIF 
         ENDDO 
         WRITE (output_io, 3110) (at_lis (k), k = 1, j) 
      ENDIF 
!                                                                       
      IF (pdf_lweights) then 
         WRITE (output_io, 3500) 'user defined' 
      ELSE 
         WRITE (output_io, 3500) 'internal' 
      ENDIF 
!                                                                       
      DO i = 1, cr_nscat 
      DO j = i, cr_nscat 
      WRITE (output_io, 3510) at_name (i), at_name (j), pdf_weight (i, j)
      ENDDO 
      ENDDO 
!                                                                       
      IF (cmd.ne.'PDF') WRITE (output_io, 4000) 
      IF (cmd.eq.'ALL'.or.cmd.eq.'RMC') then 
         CALL rmc_show ('MOD') 
         CALL rmc_show ('ATO') 
      ELSE 
         CALL rmc_show (cmd) 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Current PDF calculation settings : ') 
 2000 FORMAT     ('   Maximum r [A]              : ',F8.4,/             &
     &                   '   Grid size DR [A]           : ',F8.4,4X,    &
     &                                                 '(',I13,' pts)') 
 2010 FORMAT     ('   Radiation                  : ',A8,4X,'(at ',F8.4, &
     &                                                 ' A**-1)')       
 2015 FORMAT     ('   Radiation                  : ',A8) 
 2100 FORMAT     (/,'   Applied corrections        : ') 
 2110 FORMAT     ('     Q termination (SINC)     : not applied') 
 2120 FORMAT     ('     Q termination (SINC)     : applied, Qmax = ',   &
     &                                                    F8.4,' A**-1')
 2210 FORMAT     ('     Instrument resolution    : not applied') 
 2220 FORMAT     ('     Instrument resolution    : applied, sigQ = ',   &
     &                                                    F8.4,' A**-1')
 2300 FORMAT     ('     Periodic boundaries      : ',a) 
 2310 FORMAT     ('     Particle size            : ','infinite') 
 2320 FORMAT ('     Particle size is sphere  : ',F8.4,' A diameter') 
 2330 FORMAT ('     Particle size is         : treated by polynomial'   &
     &                   '              parameters      : ',5(F8.4,2x)) 
 2350 FORMAT ('     4 Pi Rho r correction    : none') 
 2360 FORMAT ('     4 Pi Rho r correction    : ',F8.4,              &
     &                                         ' A diameter sphere')    
 2370 FORMAT ('     4 Pi Rho r correction    : ',F8.4,              &
     &                                         ' A diameter polynomial')
 2371 FORMAT     ('              ',a9,   '       : ',F12.8) 
 2380 FORMAT     ('     4 Pi Rho r correction    : treated by tanh',/,  &
     &                   '        Particle diameter     : ',F8.4,/,     &
     &                   '        Particle shape param. : ',F8.4)       
 2390 FORMAT     ('     Weight correction        : ',F12.8) 
 2400 FORMAT     ('     Convolution therm. Gauss.: ',a) 
 2405 FORMAT     ('     Resolution broadening    : ',F8.4) 
 2410 FORMAT     ('     Quad. correlation fac.   : ',F8.4) 
 2415 FORMAT     ('     Linear correlation fac.  : ',F8.4) 
 2420 FORMAT     ('     PDF peak width ratio     : ',F8.4, ' below ',   &
     &                                                      F8.4,' A')  
 2443 FORMAT     ('     Number density           : automatic') 
 2445 FORMAT     ('     Number density           : ',F8.4,3x,           &
     &                   '(Mode: ',a,')')                               
 2450 FORMAT     ('     Correction for RHO0      : ',F8.4) 
 2455 FORMAT     ('     PDF calculation mode     : ',A) 
 3000 FORMAT     (/,'   Refinement settings        : ') 
 3010 FORMAT ('     Fit minimum r [A]        : ',F8.4,3X,'(pt.',I5,')') 
 3020 FORMAT ('     Fit maximum r [A]        : ',F8.4,3X,'(pt.',I5,')') 
 3030 FORMAT ('     Current scale factor     : ',F8.4,3x,'(refined = ', &
     &                   L1,')')                                        
 3100 FORMAT     (/,'   Selected atoms for PDF calculation :',/         &
     &                   '     Atoms (i)                : ',100(A9,1X)) 
 3110 FORMAT     ('     Atoms (j)                : ',100(A9,1X)) 
 3500 FORMAT     (/,'   Weights for PDF calculation (',a,') :') 
 3510 FORMAT     ('     ',a9,' - ',a9,'    : ',g18.8) 
 4000 FORMAT     (/,' Current RMC settings (not all used here) :') 
!                                                                       
      END SUBROUTINE pdf_show                       
!*****7*****************************************************************
      SUBROUTINE pdf_save (zeile, lp) 
!+                                                                      
!     Save PDF or current structure                                     
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
      include'debug.inc' 
       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER i, il, lp, nmi, nma 
      REAL r 
!                                                                       
      CHARACTER(1024) cdummy, cpara (maxw) 
      INTEGER ianz, lpara (maxw), ip, is 
      REAL werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ianz.ge.2) then 
         CALL do_cap (cpara (1) ) 
         cdummy = cpara (1) 
         IF (ier_num.ne.0) return 
!                                                                       
!------ - Save structure                                                
!                                                                       
         IF (cdummy (1:3) .eq.'STR') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            WRITE (output_io, 1500) cpara (1) (1:lpara (1) ) 
            CALL save_struc (cpara (1), lpara (1) ) 
!                                                                       
!------ - Save PDF (G(r))                                               
!                                                                       
         ELSEIF (cdummy (1:3) .eq.'PDF') then 
                                                                        
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            WRITE (output_io, 1505) cpara (1) (1:lpara (1) ) 
            cdummy = cpara (1) (1:lpara (1) ) 
            CALL oeffne (57, cdummy, 'unknown', .false.) 
            IF (ier_num.eq.0) then 
               nmi = int (pdf_rfmin / pdf_deltar) 
               nma = int (pdf_rfmax / pdf_deltar) 
               DO i = nmi, nma 
               r = float (i) * pdf_deltar 
               WRITE (57, 5000) r, pdf_skal * pdf_calc (i), 0.0, 1.0 
               ENDDO 
               CLOSE (57) 
            ENDIF 
!                                                                       
!------ - Save markers                                                  
!                                                                       
         ELSEIF (cdummy (1:3) .eq.'MAR') then 
            IF (pdf_gauss) then 
               ier_num = - 7 
               ier_typ = ER_PDF 
            ELSE 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               WRITE (output_io, 1508) cpara (1) (1:lpara (1) ) 
               cdummy = cpara (1) (1:lpara (1) ) 
               CALL oeffne (57, cdummy, 'unknown', .false.) 
               IF (ier_num.eq.0) then 
                  nmi = int (pdf_rfmin / pdf_deltar) 
                  nma = int (pdf_rfmax / pdf_deltar) 
                  DO i = nmi, nma 
                  r = float (i) * pdf_deltar 
                  IF (pdf_calc (i) .gt.0.0) then 
                     WRITE (57, 5100) r, 0.0 
                  ENDIF 
                  ENDDO 
                  CLOSE (57) 
               ENDIF 
            ENDIF 
!                                                                       
!------ --unknown subcommand given                                      
!                                                                       
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1500 FORMAT     (' Saving structure to file : ',A,' ...') 
 1505 FORMAT     (' Saving PDF to file : ',A,' ...') 
 1508 FORMAT     (' Saving distance markers to file : ',A,' ...') 
 5000 FORMAT     (F9.4,3X,F21.10,5X,2(F6.2,1X)) 
 5100 FORMAT     (F9.4,3X,F6.2) 
      END SUBROUTINE pdf_save                       
!*****7*****************************************************************
      SUBROUTINE pdf_readdata (zeile, lp) 
!+                                                                      
!     Reads observed PDF as xy ASCII file.                              
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
      include'debug.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) datafile 
      INTEGER lpara (maxw) 
      INTEGER ianz, ip 
      REAL ra, re, dr 
      REAL werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
      datafile = cpara (1) (1:lpara (1) ) 
!                                                                       
!------ Read observed PDF for given plane                               
!                                                                       
      IF (ianz.eq.1) then 
         CALL oeffne (17, datafile, 'old', .false.) 
         IF (ier_num.ne.0) return 
         CALL extract_hist (17) 
         CALL skip_spec (17) 
         ip = 1 
         READ (17, *, end = 20, err = 999) ra, pdf_obs (ip), dr, pdf_wic (ip)
         ip = ip + 1 
   10    CONTINUE 
         READ (17, *, end = 20, err = 999) re, pdf_obs (ip), dr, pdf_wic (ip)
         ip = ip + 1 
         IF (ip.gt.PDF_MAXDAT) goto 9999 
         GOTO 10 
   20    CONTINUE 
         CLOSE (17) 
         IF (ier_num.ne.0) return 
!                                                                       
         pdf_rmax = re 
         pdf_bin = ip - 1 
         pdf_deltar = (re-ra) / float (pdf_bin - 1) 
         pdf_rfmin = ra 
         pdf_rfmax = re 
!                                                                       
         IF (abs (ra - pdf_deltar) .gt.1e-6) then 
            ier_num = - 5 
            ier_typ = ER_PDF 
            RETURN 
         ENDIF 
!                                                                       
         WRITE (output_io, 9000) ra, re, pdf_bin 
         pdf_ldata = .true. 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      RETURN 
!                                                                       
  999 CONTINUE 
      CLOSE (17) 
      ier_num = - 3 
      ier_typ = ER_IO 
!                                                                       
 9999 CONTINUE 
      CLOSE (17) 
      ier_num = - 1 
      ier_typ = ER_PDF 
      RETURN 
!                                                                       
 9000 FORMAT     (' Read PDF data (r = ',F7.3,' to ',F7.3,' A, ',I5,    &
     &                   ' points) ...')                                
      END SUBROUTINE pdf_readdata                   
!*****7*****************************************************************
      SUBROUTINE extract_hist (ifil) 
!                                                                       
!     Extract history information if present                            
!                                                                       
      USE config_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'param.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER ifil 
!                                                                       
      CHARACTER(1024) line 
      INTEGER is, ie 
!                                                                       
      INTEGER len_str 
!                                                                       
      READ (ifil, 5000) line 
      IF (line (1:7) .eq.'History') then 
         res_para (0) = 0 
   22    CONTINUE 
         READ (ifil, 5000) line 
         CALL extract_key (line, 'temp=') 
         CALL extract_key (line, 'Qmax=') 
         IF (line (1:16) .eq.'##### start data') goto 33 
         GOTO 22 
   33    CONTINUE 
         pdf_qmax = res_para (2) 
         WRITE (output_io, 1000) pdf_qmax 
      ELSE 
         BACKSPACE (ifil) 
      ENDIF 
!                                                                       
 1000 FORMAT    ( ' History information found, setting Qmax=',f5.1,     &
     &                   ' A**-1 ...')                                  
 5000 FORMAT    (a) 
!                                                                       
      END SUBROUTINE extract_hist                   
!*****7*****************************************************************
      SUBROUTINE skip_spec (ifil) 
!+                                                                      
!     Reads over #xxx lines (SPEC headers)                              
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER ifil 
!                                                                       
      CHARACTER(1) cstr 
!                                                                       
   11 CONTINUE 
      READ (ifil, 1000, end = 22) cstr 
      IF (cstr (1:1) .eq.'#') goto 11 
   22 CONTINUE 
      BACKSPACE (ifil) 
!                                                                       
 1000 FORMAT     (a) 
      END SUBROUTINE skip_spec                      
!*****7*****************************************************************
      SUBROUTINE extract_key (line, key) 
!+                                                                      
!     Gets numbers from history part of data file                       
!-                                                                      
      USE config_mod 
      IMPLICIT none 
!                                                                       
      include'param.inc' 
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
         READ (line (is:is + ie), * ) res_para (nint (res_para (0) ) ) 
      ENDIF 
!                                                                       
      END SUBROUTINE extract_key                    
!*****7**************************************************************** 
      SUBROUTINE pdf_set (zeile, lp) 
!+                                                                      
!     Sets most parameters for 'pdf' section                            
!-                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE modify_mod
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw), wa (maxw), wb (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, nn 
      INTEGER i, j, ia, ib 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.ge.2) then 
         CALL do_cap (cpara (1) ) 
!                                                                       
!------ - set boundary: toggle periodic boundaries                      
!                                                                       
         IF (cpara (1) (1:3) .eq.'BOU') then 
            IF (ianz.eq.2.or.ianz.eq.3) then 
               chem_period (1) = str_comp (cpara (2) , 'period', 3,     &
               lpara (2) , 6)                                           
               chem_period (2) = chem_period (1) 
               chem_period (3) = chem_period (1) 
               IF (chem_period (1) ) then 
                  IF (ianz.eq.3) then 
                     pdf_2d = str_comp (cpara (3),'2D',1,lpara(3), 2)
                     pdf_lexact = .false. 
                     IF (pdf_2d) then 
      WRITE (output_io, 1000) 'periodic bound. 2D, unit cell' 
                     ELSE 
      WRITE (output_io, 1000) 'periodic bound. 3D, unit cell' 
                     ENDIF 
                  ENDIF 
               ELSE 
                  IF (ianz.eq.3) then 
                     pdf_lexact = str_comp (cpara(3),'exact',2,lpara (3),5)
                  ENDIF 
                  IF (pdf_lexact) then 
      WRITE (output_io, 1000) 'no periodic bound., exact' 
                  ELSE 
      WRITE (output_io, 1000) 'no periodic bound., unit cell' 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set mode: toggle calculation mode                             
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'CAL') then 
            IF (ianz.eq.2) then 
               pdf_lexact = str_comp (cpara (2) , 'exact', 2, lpara (2),5)
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set finite: finite size correction model                      
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'FIN') then 
            IF (ianz.ge.2) then 
               IF (str_comp (cpara (2) , 'peri', 3, lpara (2) , 4) ) then
                  pdf_finite = PDF_BACK_PERIOD 
               ELSEIF (str_comp (cpara (2) ,'sphe',3,lpara(2),4)) then
                  pdf_finite = PDF_BACK_SPHERE 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) return 
                  IF (ianz.eq.1) then 
                     pdf_sphere = werte (1) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSEIF (str_comp (cpara (2),'poly',3,lpara(2),4)) then
                  pdf_finite = PDF_BACK_POLY 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) return 
                  IF (ianz.eq.1) then 
                     pdf_diam_poly = werte (1) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSEIF (str_comp (cpara (2),'tanh',3,lpara (2),4)) then
                  pdf_finite = PDF_BACK_TANH 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set qalpha: sets peak broadening factor                       
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'QAL') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_qalp = werte (1) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set gamma: sets linear correlation factor                     
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'GAM') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_gamma = werte (1) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set delta: sets quadratic correlation factor                  
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'DEL') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_delta = werte (1) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set density: sets number density                              
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'DEN'.or.             &
                 cpara (1) (1:3) .eq.'RDE')       then
            pdf_lrho0_rel = (cpara (1) (1:3) .eq.'RDE') 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ianz.eq.1) then 
               IF (str_comp (cpara (1) , 'auto', 2, lpara (1) , 4) ) then
                  pdf_lrho0 = .true. 
               ELSE 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) return 
                  pdf_rho0 = werte (1) 
                  pdf_lrho0 = .false. 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set diam: sets diameter for finite particles                  
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'DIA') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_diam = abs (werte (1) ) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set frange: sets range in r to be used for fit                
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'FRA') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.2) then 
               pdf_rfmin = werte (1) 
               pdf_rfmax = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set qmax: sets Qmax for termination correction                
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'QMA') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_qmax = werte (1) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set qsig: sets SIGMA Q for resolution correction              
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'QSI') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_sigmaq = werte (1) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set partial: weights for PDF partials                         
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'PAR') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (str_comp (cpara(1),'internal',3,lpara(1),8)) then
               pdf_lweights = .false. 
            ELSE 
               IF (ianz.ge.3) then 
                  ia = 1 
                  CALL get_iscat (ia, cpara, lpara, wa, maxw, .false.) 
                  IF (wa (1) .lt.0) then 
                     ier_num = - 8 
                     ier_typ = ER_PDF 
                  ENDIF 
                  IF (ier_num.ne.0) return 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
!                                                                       
                  ib = 1 
                  CALL get_iscat (ib, cpara, lpara, wb, maxw, .false.) 
                  IF (wb (1) .lt.0) then 
                     ier_num = - 8 
                     ier_typ = ER_PDF 
                  ENDIF 
                  IF (ier_num.ne.0) return 
!                                                                       
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) return 
!                                                                       
                  DO i = 1, ia 
                  DO j = 1, ib 
                  pdf_weight (NINT(wa (i)), NINT(wb (j)) ) = werte (1) 
                  pdf_weight (NINT(wb (j)), NINT(wa (i)) ) = werte (1) 
                  ENDDO 
                  ENDDO 
!                                                                       
                  pdf_lweights = .true. 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!------ - set poly: background polynomial                               
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'POL') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.ge.1) then 
               pdf_poly_n = ianz 
               DO i = 1, ianz 
               pdf_poly (i) = werte (i) 
               ENDDO 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set radiation: sets radiation uses                            
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'RAD') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1.or.ianz.eq.2) then 
               CALL do_cap (cpara (1) ) 
               IF (cpara (1) (1:1) .eq.'N') then 
                  pdf_lxray = .false. 
               ELSEIF (cpara (1) (1:1) .eq.'X') then 
                  pdf_lxray = .true. 
                  IF (ianz.eq.2) then 
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.ne.0) return 
                     pdf_xq = werte (1) 
                  ELSE 
                     pdf_xq = 0.0 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set range: set range in r for PDF calculation                 
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'RAN') then 
            IF (pdf_ldata) then 
               ier_num = - 6 
               ier_typ = ER_PDF 
               RETURN 
            ENDIF 
!                                                                       
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.2) then 
               nn = int (werte (1) / werte (2) )  + 1
               IF (nn.lt.1) then 
                  ier_num = - 1 
                  ier_typ = ER_PDF 
                  RETURN 
               ELSEIF (nn.gt.PDF_MAXDAT) then
                 pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
                 pdf_ndat  = MAX(pdf_ndat , nn      , PDF_MAXDAT)
                 pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
                 CALL alloc_pdf( pdf_nscat, pdf_ndat, pdf_nbnd )
               ENDIF 
               pdf_rmax = werte (1) 
               pdf_deltar = werte (2) 
               pdf_rfmin = werte (2) 
               pdf_rfmax = pdf_rmax 
               pdf_bin = nn 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set weight: sets scale parameter to correct the weighting     
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'WEI') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_scale = abs (werte (1) ) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set shap: sets shape parameter for finite particles           
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'SHA') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_shape = abs (werte (1) ) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set srat: sets peak sharpening SRAT                           
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'SRA') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.2) then 
               pdf_srat = werte (1) 
               pdf_rcut = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set therm: toggle convolution with thermal Gaussian           
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'THE') then 
            IF (ianz.eq.2) then 
               pdf_gauss = str_comp (cpara(2),'gaus',3,lpara(2), 4)                                                       
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - if not found here try rmc_set (that is where the rest is)     
!                                                                       
         ELSEIF (cpara (1) (1:3) .eq.'CYC'.or.      &
                 cpara (1) (1:3) .eq.'DIS'.or.      &
                 cpara (1) (1:3) .eq.'SIG'.or.      &
                 cpara (1) (1:3) .eq.'MOD'.or.      &
                 cpara (1) (1:3) .eq.'MDI'.or.      &
                 cpara (1) (1:3) .eq.'SCA'.or.      &
                 cpara (1) (1:3) .eq.'MOV'     ) then
            CALL rmc_set (zeile, lp) 
!                                                                       
!------ - invalid command entered                                       
!                                                                       
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (1x,'Setting PDF calculation to ',a,' mode ..') 
      END SUBROUTINE pdf_set                        
!*****7*****************************************************************
      SUBROUTINE pdf_run 
!+                                                                      
!     Main PDF fit loop - called by run command                         
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE pdf_mod 
      USE rmc_mod 
      USE structur, ONLY: update_cr_dim
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
      include'errlist.inc' 
      include'param.inc' 
!                                                                       
      REAL(8) cc, c, ce, e, ee, wtot 
!     REAL pdf_old (MAXDAT) 
      REAL, DIMENSION(PDF_MAXDAT) ::  pdf_old !  (MAXDAT) 
      REAL cnew, cold, sig2, sum 
      REAL prob, psum, p2sum, pave, psig, pmax, pn 
      REAL start, zeit, seknds 
      REAL p_new (3, rmc_max_atom) 
      REAL p_old (3, rmc_max_atom) 
      INTEGER i_new (rmc_max_atom) 
      INTEGER i_old (rmc_max_atom) 
      INTEGER isel (rmc_max_atom), natoms 
      INTEGER imol (rmc_max_atom) 
      INTEGER zh, zm, zs, nmi, nma 
      INTEGER i, j, ip, iq, is, iii 
      INTEGER igen, itry, iacc_good, iacc_bad 
      LOGICAL loop, laccept 
!                                                                       
      REAL ran1 
!                                                                       
      IF (pdf_obs (1) .eq. - 9999.) then 
         ier_num = - 4 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
!                                                                       
      igen = 0 
      itry = 0 
      iacc_good = 0 
      iacc_bad = 0 
      loop = .true. 
      laccept = .true. 
!                                                                       
      psum = 0.0 
      p2sum = 0.0 
      pmax = 0.0 
      pn = 0.0 
!                                                                       
      nmi = int (pdf_rfmin / pdf_deltar) 
      nma = int (pdf_rfmax / pdf_deltar) 
!                                                                       
      CALL chem_aver (.false.) 
!                                                                       
!------ Calculate maximal value in rmc_mindist array                    
!                                                                       
      rmc_mindist_max = rmc_mindist (1, 1) 
      DO i = 1, cr_nscat 
      DO j = 1, cr_nscat 
      rmc_mindist_max = max (rmc_mindist_max, rmc_mindist (i, j) ) 
      ENDDO 
      ENDDO 
!                                                                       
!------ Write some start information                                    
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), cr_natoms 
!                                                                       
!------ calculate sums from exp. data needed and initial chi2           
!                                                                       
      CALL pdf_determine (.false.) 
!                                                                       
      cold = 0.0 
      wtot = 0.0 
      e = 0.0 
      ee = 0.0 
      c = 0.0 
      cc = 0.0 
      ce = 0.0 
!                                                                       
      DO ip = nmi, nma 
      wtot = wtot + pdf_wic (ip) 
      e = e+pdf_wic (ip) * pdf_obs (ip) 
      ee = ee+pdf_wic (ip) * pdf_obs (ip) **2 
      c = c + pdf_wic (ip) * pdf_calc (ip) 
      cc = cc + pdf_wic (ip) * pdf_calc (ip) **2 
      ce = ce+pdf_wic (ip) * pdf_calc (ip) * pdf_obs (ip) 
      ENDDO 
      IF (rmc_doskal) then 
         pdf_skal = ce / cc 
      ELSE 
         pdf_skal = 1.0 / rmc_skal (1) 
      ENDIF 
      cold = ee+pdf_skal**2 * cc - 2.0 * pdf_skal * ce 
      cold = cold / wtot 
!                                                                       
!------ Here is the main loop                                           
!                                                                       
      start = seknds (0.0) 
      IF (rmc_mode.eq.rmc_mode_swchem) then 
         WRITE (output_io, 2000) 'Mode SWCHEM' 
      ELSEIF (rmc_mode.eq.rmc_mode_swdisp) then 
         WRITE (output_io, 2000) 'Mode SWDISP' 
      ELSEIF (rmc_mode.eq.rmc_mode_shift) then 
         WRITE (output_io, 2000) 'Mode SHIFT' 
      ELSEIF (rmc_mode.eq.rmc_mode_external) then 
         WRITE (output_io, 2000) 'external defined mode' 
      ENDIF 
!                                                                       
      sig2 = rmc_sigma**2 / 2.0 
!                                                                       
      DO while (loop) 
      laccept = .true. 
      igen = igen + 1 
!                                                                       
!-------- generate move and check for limits                            
!                                                                       
      IF (rmc_sel_atom) then 
         CALL rmc_genmove (laccept, natoms, p_new, i_new, isel) 
      ELSE 
         CALL rmc_genmove_mol (laccept, natoms, p_new, i_new, isel,     &
         imol)                                                          
      ENDIF 
      IF (ier_num.ne.0) return 
      IF (laccept) then 
!                                                                       
!-------- - save old positions ...                                      
!                                                                       
         DO i = 1, natoms 
         i_old (i) = cr_iscat (isel (i) ) 
         p_old (1, i) = cr_pos (1, isel (i) ) 
         p_old (2, i) = cr_pos (2, isel (i) ) 
         p_old (3, i) = cr_pos (3, isel (i) ) 
         ENDDO 
!                                                                       
         DO i = 1, pdf_bin 
         pdf_old (i) = pdf_corr (i) 
         ENDDO 
!                                                                       
!-------- - Calc new PDF and chi2                                       
!                                                                       
         DO i = 1, natoms 
         CALL pdf_addcorr (isel (i), - 2.0, sum) 
         ENDDO 
         CALL pdf_makemove (natoms, i_new, p_new, isel, imol) 
         DO i = 1, natoms 
         CALL pdf_addcorr (isel (i), 2.0, sum) 
         ENDDO 
         CALL pdf_convert 
!                                                                       
         itry = itry + 1 
         cnew = 0.0 
         c = 0.0 
         cc = 0.0 
         ce = 0.0 
!                                                                       
         DO ip = nmi, nma 
         c = c + pdf_wic (ip) * pdf_calc (ip) 
         cc = cc + pdf_wic (ip) * pdf_calc (ip) **2 
         ce = ce+pdf_wic (ip) * pdf_calc (ip) * pdf_obs (ip) 
         ENDDO 
         IF (rmc_doskal) then 
            pdf_skal = ce / cc 
         ELSE 
            pdf_skal = 1.0 / rmc_skal (1) 
         ENDIF 
         cnew = ee+pdf_skal**2 * cc - 2.0 * pdf_skal * ce 
         cnew = cnew / wtot 
!                                                                       
!     ----Accept move ?                                                 
!                                                                       
         prob = cnew - cold 
!                                                                       
         IF (prob.lt.0) then 
            laccept = .true. 
         ELSE 
            IF (sig2.gt.0.0) then 
               psum = psum + prob 
               p2sum = p2sum + prob**2 
               pmax = max (pmax, prob) 
               pn = pn + 1 
               prob = exp ( - prob / sig2) 
               laccept = (prob.gt.ran1 (idum) ) 
            ELSE 
               laccept = .false. 
            ENDIF 
         ENDIF 
!                                                                       
         IF (rmc_sigma.eq. - 9999.) laccept = .true. 
!                                                                       
!------ ----if accepted make move                                       
!                                                                       
         IF (laccept) then 
            cold = cnew 
            IF (prob.lt.0) then 
               iacc_good = iacc_good+1 
            ELSE 
               iacc_bad = iacc_bad+1 
            ENDIF 
         ELSE 
            CALL pdf_makemove (natoms, i_old, p_old, isel, imol) 
            DO i = 1, pdf_bin 
            pdf_corr (i) = pdf_old (i) 
            ENDDO 
            CALL pdf_convert 
         ENDIF 
!                                                                       
!------ --WRITE info and terminate or loop again                        
!                                                                       
      ENDIF 
!                                                                       
 2222 CONTINUE 
!                                                                       
      loop = (itry.lt.rmc_maxcyc) 
!                                                                       
      IF (igen.gt.1000 * rmc_display.and.itry.eq.0) then 
         ier_num = - 20 
         ier_typ = ER_RMC 
         loop = .false. 
      ENDIF 
!                                                                       
      IF (mod (igen, rmc_display) .eq.0.or..not.loop) then 
         IF (.not.loop) WRITE (output_io, 1250) pdf_skal 
         WRITE (output_io, 1300) igen, itry, iacc_good, iacc_bad, cold 
      ENDIF 
      ENDDO 
!                                                                       
!------ WRITE timing summary and "prob" statistics                      
!                                                                       
      IF (pn.gt.0) then 
         pave = psum / pn 
         psig = sqrt (p2sum / pn - (psum / pn) **2) 
         WRITE (output_io, 3000) pave, psig, pmax 
      ENDIF 
!                                                                       
      zeit = seknds (start) 
      zh = int (zeit / 3600.) 
      zm = int ( (zeit - zh * 3600.) / 60.) 
      zs = int (zeit - zh * 3600 - zm * 60.) 
      WRITE (output_io, 4000) zh, zm, zs, zeit / itry 
!                                                                       
      rmc_skal (1) = 1.0 / pdf_skal 
!                                                                       
!------ save some results to res[i] blo                                 
!                                                                       
      res_para (0) = 8 
!                                                                       
      res_para (1) = cold 
      res_para (2) = float (itry) 
      res_para (3) = float (iacc_good) 
      res_para (4) = float (iacc_bad) 
      res_para (5) = pave 
      res_para (6) = psig 
      res_para (7) = pmax 
      res_para (8) = zeit / float (itry) 
!                                                                       
!------ Update crystal dimensions                                       
!                                                                       
      CALL update_cr_dim 
!                                                                       
!------ formats                                                         
!                                                                       
 1000 FORMAT (' Running PDF-Fit ...',//                                 &
     &        ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,     &
     &        ' containing ',I9,' atoms')                               
 1250 FORMAT (/,' ---- Final configuration (s= ',F8.4,') ---- ',/) 
 1300 FORMAT (  ' Gen: ',I6,' try: ',I6,' acc: (good/bad): ',I6,        &
     &          ' / ',I6,' s2x2: ',G15.8)                               
 2000 FORMAT (/,' Starting main RMC loop ...',/,' (',A,') ',/) 
 3000 FORMAT (/,' Delta Chi    : ave: ',G15.4,' sig: ',G15.4,           &
     &          ' max: ',G15.4)                                         
 4000 FORMAT (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/     &
     &          ' Time/cycle   : ',F9.3,' sec',/)                       
      END SUBROUTINE pdf_run                        
!*****7*****************************************************************
      SUBROUTINE pdf_makemove (natoms, i_new, p_new, isel, imol) 
!+                                                                      
!     make accepted move                                                
!                                                                       
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL p_new (3, rmc_max_atom) 
      INTEGER i_new (rmc_max_atom) 
      INTEGER isel (rmc_max_atom) 
      INTEGER imol (rmc_max_atom) 
      INTEGER natoms 
!                                                                       
      INTEGER i, j, is 
!                                                                       
      DO i = 1, natoms 
      cr_iscat (isel (i) ) = i_new (i) 
      DO j = 1, 3 
      cr_pos (j, isel (i) ) = p_new (j, i) 
      ENDDO 
      ENDDO 
!                                                                       
      IF (.not.rmc_sel_atom.and.rmc_mode.eq.rmc_mode_swchem) then 
         is = mole_type (imol (1) ) 
         mole_type (imol (1) ) = mole_type (imol (2) ) 
         mole_type (imol (2) ) = is 
      ENDIF 
!                                                                       
      END SUBROUTINE pdf_makemove                   
!*****7*****************************************************************
      SUBROUTINE pdf_determine (lout) 
!+                                                                      
!     Calculate PDF of current structure                                
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
      include'debug.inc' 
      include'errlist.inc' 
       
!                                                                       
      INTEGER i, ia, id 
      INTEGER is, js 
      REAL done, sum 
      INTEGER nmi, nma 
      REAL r 
      LOGICAL lout 
      REAL seknds, ss 
!                                                                       
      ss = seknds (0.0) 
!                                                                       
      IF(.NOT. pdf_lexact) THEN
         IF (cr_ncatoms.eq.0) THEN
            ier_num = -9
            ier_typ = ER_PDF
            ier_msg(1) = 'Quick mode and periodic boundaries'
            ier_msg(2) = 'must be disabled. Use exact mode.'
            ier_msg(3) = 'Crystal created in stack menue? '
            RETURN
         ENDIF
      ENDIF
      IF (lout) then 
         IF (pdf_lexact) then 
            WRITE (output_io, 500) 'exact mode' 
         ELSE 
            WRITE (output_io, 500) 'unit cell mode' 
         ENDIF 
      ENDIF 
!                                                                       
!------ Reset arrays                                                    
!                                                                       
      DO i = 1, PDF_MAXDAT 
      pdf_corr (i) = 0.0 
      ENDDO 
      DO i = 1, PDF_MAXDAT 
      DO is = 1, cr_nscat 
      DO js = 1, cr_nscat 
      pdf_temp (i, is, js) = 0 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      sum = 0.0 
!                                                                       
!------ Start the calculation                                           
!                                                                       
      id = max (1, cr_natoms / 5) 
!                                                                       
      IF (pdf_lexact) then 
         DO ia = 1, cr_natoms 
         CALL pdf_addcorr_e (ia) 
         IF (lout.and. (mod (ia, id) .eq.0) ) then 
            done = 100.0 * float (ia) / float (cr_natoms) 
            WRITE (output_io, 1000) done 
         ENDIF 
         ENDDO 
      ELSE 
         DO ia = 1, cr_natoms 
         CALL pdf_addcorr_n (ia) 
         IF (lout.and. (mod (ia, id) .eq.0) ) then 
            done = 100.0 * float (ia) / float (cr_natoms) 
            WRITE (output_io, 1000) done 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      CALL pdf_convtherm (1.0, sum) 
!                                                                       
!------ Convert to proper G(r)                                          
!                                                                       
      CALL pdf_convert 
!                                                                       
      ss = seknds (ss) 
      IF (lout) WRITE (output_io, 3000) ss 
!                                                                       
      IF (lout) WRITE (output_io, 2000) sum 
!                                                                       
  500 FORMAT    (  ' Calculating PDF (',a,') ...') 
 1000 FORMAT     (  '   ',f6.2,' % done  ...') 
 2000 FORMAT     (/,' Total sum of weights : ',g18.10) 
 3000 FORMAT     (  ' Required time        : ',f8.2,' s') 
      END SUBROUTINE pdf_determine                  
!*****7*****************************************************************
      SUBROUTINE pdf_convert 
!+                                                                      
!     Convert to G(r) and do convolution                                
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE atom_env_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'debug.inc' 
       
      include'errlist.inc' 
      include'param.inc' 
      include'wink.inc' 
!                                                                       
      INTEGER i, k, ia, ncc 
!     REAL ppp (MAXDAT) 
      REAL, DIMENSION(PDF_MAXDAT   ) :: ppp ! (MAXDAT) 
      REAL norm, fac, r, r0 
      REAL rr 
      INTEGER nmi, nma 
!                                                                       
      ncc = cr_icc (1) * cr_icc (2) * cr_icc (3) 
      IF (.not.pdf_lrho0) then 
         IF (pdf_lrho0_rel) then 
            r0 = pdf_rho0 * cr_ncatoms / cr_v 
         ELSE 
            r0 = pdf_rho0 
         ENDIF 
      ELSE 
         r0 = cr_ncatoms / cr_v 
      ENDIF 
      norm = 1.0 / float (cr_n_real_atoms) * pdf_scale 
      IF (.not.pdf_gauss.and.pdf_qalp.eq.0.0) then 
         norm = norm / pdf_deltar 
      ENDIF 
!                                                                       
      DO i = 1, pdf_bin 
      r = float (i) * pdf_deltar 
      IF (pdf_finite.eq.PDF_BACK_PERIOD) then 
         rr = 2.0 * zpi * r * r0 * pdf_dnorm 
      ELSEIF (pdf_finite.eq.PDF_BACK_POLY) then 
         rr = 2.0 * zpi * r * r0 * pdf_dnorm 
         IF (r.lt.pdf_diam_poly) then 
            DO k = 1, pdf_poly_n 
               rr = rr - pdf_poly (k) * r**k 
            ENDDO 
            rr = max (0.0, rr) 
         ELSE 
            rr = 0.0 
         ENDIF 
      ELSEIF (pdf_finite.eq.PDF_BACK_SPHERE) then 
         rr = 2.0 * zpi * r * r0 * pdf_dnorm 
         IF (r.lt.pdf_sphere) then 
            rr = rr * (1. - 1.5 * (r / pdf_sphere) + .5 * (r /          &
            pdf_sphere) **3)                                            
         ELSE 
            rr = 0.0 
         ENDIF 
      ELSEIF (pdf_finite.eq.PDF_BACK_TANH) then 
         rr = max (0.0, - 2.0 * zpi * r * r0 * pdf_dnorm * tanh (       &
         pdf_shape * (r - pdf_diam) ) )                                 
      ENDIF 
      IF (chem_period (1) ) then 
         IF (pdf_finite.eq.PDF_BACK_PERIOD) then 
            pdf_calc (i) = norm * pdf_corr (i) / r - rr 
         ELSEIF (pdf_finite.eq.PDF_BACK_SPHERE) then 
            IF (r.lt.pdf_sphere) then 
                                                                        
               pdf_calc (i) = norm * pdf_corr (i) / r * (1. - 1.5 *     &
               (r / pdf_sphere) + .5 * (r / pdf_sphere) **3) - rr       
            ELSE 
               pdf_calc (i) = 0.0 
            ENDIF 
         ENDIF 
      ELSE 
         pdf_calc (i) = norm * pdf_corr (i) / r - rr 
      ENDIF 
      ENDDO 
!                                                                       
!------ Apply instrument resolution correction                          
!                                                                       
      IF (pdf_sigmaq.gt.0.0) then 
         DO i = 1, pdf_bin 
         r = float (i) * pdf_deltar 
         pdf_calc (i) = pdf_calc (i) * exp ( - (r * pdf_sigmaq) **2 /   &
         2.0)                                                           
         ENDDO 
      ENDIF 
!                                                                       
!------ Convolute with SINC function                                    
!                                                                       
      IF (pdf_qmax.gt.0.0) then 
         DO i = 1, pdf_bin 
         ppp (i) = pdf_calc (i) * (pdf_qmax - pdf_sinc (2 * i) ) 
         DO k = 1, i - 1 
         ppp (i) = ppp (i) + pdf_calc (k) * (pdf_sinc (i - k) -         &
         pdf_sinc (i + k) )                                             
         ENDDO 
         DO k = i + 1, pdf_bin 
         ppp (i) = ppp (i) + pdf_calc (k) * (pdf_sinc (k - i) -         &
         pdf_sinc (k + i) )                                             
         ENDDO 
         ENDDO 
!                                                                       
         DO i = 1, pdf_bin 
         pdf_calc (i) = ppp (i) * pdf_deltar / zpi * 2.0 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE pdf_convert                    
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr (ia, rsign, sum) 
!+                                                                      
!     Calculate correlation for given atom ia                           
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE atom_env_mod 
      USE modify_mod
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'param.inc' 
      include'wink.inc' 
!                                                                       
      INTEGER ig, igaus, ib, ie 
      INTEGER i, j, k, ii, jj, kk, is, js, ks, ia, iatom, ibin, nn, ncc 
      INTEGER istart (3), iend (3), iii (3), cell (3) 
!     REAL(8) ppp (MAXDAT), gaus ( - MAXDAT:MAXDAT) 
      REAL(8), DIMENSION( PDF_MAXDAT)            :: ppp   !(MAXDAT)
      REAL(8), DIMENSION(-PDF_MAXDAT:PDF_MAXDAT) :: gaus  ! ( - MAXDAT:MAXDAT) 
      REAL rsign, sum 
      REAL asym, gnorm, norm, dist, dist2, r, r0, rg 
      REAL sigma, fac 
      REAL dd (3), d (3), offset (3) 
      LOGICAL lout 
!                                                                       
      fac = 1.0 / (2.0 * zpi**2) 
!                                                                       
      is = cr_iscat (ia) 
      IF (pdf_allowed_i (is) .or.pdf_allowed_j (is) ) then 
         CALL indextocell (ia, iii, ks) 
         DO i = 1, 3 
         istart (i) = iii (i) - 1 - int (pdf_rmax / cr_a0 (i) ) 
         iend (i) = iii (i) + 1 + int (pdf_rmax / cr_a0 (i) ) 
         IF (pdf_2d.and.cr_icc (i) .le.1) then 
            istart (i) = iii (i) 
            iend (i) = iii (i) 
         ENDIF 
         ENDDO 
!                                                                       
!------ - In case we do not want to loop around cut range               
!                                                                       
         IF ( (.not.chem_period (1) ) .and. (.not.chem_period (2) )     &
         .and. (.not.chem_period (3) ) ) then                           
            DO i = 1, 3 
            istart (i) = max (1, istart (i) ) 
            iend (i) = min (cr_icc (i), iend (i) ) 
            ENDDO 
         ENDIF 
!                                                                       
!------ - Here starts the inner loop                                    
!                                                                       
         DO k = istart (3), iend (3) 
         DO j = istart (2), iend (2) 
         DO i = istart (1), iend (1) 
         cell (1) = i 
         cell (2) = j 
         cell (3) = k 
!                                                                       
!---------- Here we apply periodic boundaries                           
!---------- Modified to allow more than a single cycle !                
!                                                                       
         DO ii = 1, 3 
         iii (ii) = cell (ii) 
         cell (ii) = pdf_bnd (ii, cell (ii) ) 
         offset (ii) = float (iii (ii) - cell (ii) ) 
         ENDDO 
!                                                                       
!------ --- Now look for neighbours in surrounding unit cells only      
!                                                                       
         DO ii = 1, cr_ncatoms 
         CALL celltoindex (cell, ii, iatom) 
         js = cr_iscat (iatom) 
         IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or. (      &
         pdf_allowed_j (is) .and.pdf_allowed_i (js) ) ) then            
            DO jj = 1, 3 
            dd (jj) = cr_pos (jj, ia) - cr_pos (jj, iatom) - offset (jj) 
            d (jj) = abs (dd (jj) ) * cr_a0 (jj) 
            ENDDO 
!             dist2 = skalpro(dd,dd,cr_gten)                            
            dist2 = dd (1) * dd (1) * cr_gten (1, 1) + dd (2) * dd (2)  &
            * cr_gten (2, 2) + dd (3) * dd (3) * cr_gten (3, 3) + 2. *  &
            dd (1) * dd (2) * cr_gten (1, 2) + 2. * dd (1) * dd (3)     &
            * cr_gten (1, 3) + 2. * dd (2) * dd (3) * cr_gten (2, 3)    
            dist = sqrt (dist2) 
!                                                                       
            IF (dist.le.pdf_rmax.and.dist.gt.pdf_deltar) then 
               ibin = nint (dist / pdf_deltar) 
!                                                                       
!------ --------- Convolute with Gaussian                               
!                                                                       
               IF (pdf_gauss.or.pdf_qalp.gt.0.0) then 
                  sigma = 0.0 
                  IF (pdf_gauss) then 
                     sigma = fac * (cr_dw (is) + cr_dw (js) ) 
                     sigma = sigma - pdf_delta / dist2 
                     sigma = sigma - pdf_gamma / dist 
                     sigma = max (0.0, sigma) 
                  ENDIF 
                  sigma = sigma + pdf_qalp**2 * dist2 
                  sigma = sqrt (sigma) 
!                                                                       
                  IF (dist.le.pdf_rcut) then 
                     sigma = sigma * pdf_srat 
                  ENDIF 
!                                                                       
                  igaus = 1 + nint (5.0 * sigma / pdf_deltar) 
                  ib = max (1, ibin - igaus + 1) 
                  ie = min (pdf_bin, ibin + igaus - 1) 
!                                                                       
                  IF (sigma.le.0.0.or.igaus.lt.2) then 
                     pdf_corr (ibin) = pdf_corr (ibin) + rsign *        &
                     pdf_weight (is, js) / pdf_deltar                   
                  ELSE 
                     gnorm = 1.0 / (sqrt (zpi) * sigma) 
!                                                                       
                     DO ig = - igaus, igaus 
                     rg = (ig - 1) * pdf_deltar 
                     asym = 1.0 + rg / dist 
                     gaus (ig) = gnorm * asym * exp ( - 0.5 * (rg /     &
                     sigma) **2)                                        
                     ENDDO 
!                                                                       
                     DO ig = ib, ie 
                     kk = ig - ibin + 1 
                     pdf_corr (ig) = pdf_corr (ig) + rsign * pdf_weight &
                     (is, js) * gaus (kk)                               
                     ENDDO 
                  ENDIF 
!                                                                       
!------ --------- Just take the distance                                
!                                                                       
               ELSE 
                  pdf_corr (ibin) = pdf_corr (ibin) + rsign *           &
                  pdf_weight (is, js)                                   
               ENDIF 
               sum = sum + pdf_weight (is, js) 
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
         ENDDO 
 1111    CONTINUE 
      ENDIF 
!                                                                       
      END SUBROUTINE pdf_addcorr                    
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_n (ia) 
!+                                                                      
!     Calculate correlation for given atom ia, fast version             
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE modify_mod
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER i, j, k, ii, jj, is, js, ks, ia, iatom, ibin 
      INTEGER istart (3), iend (3), iii (3), cell (3) 
      REAL dist, dist2 
      REAL dd (3), d (3), offset (3) 
!                                                                       
      is = cr_iscat (ia) 
      IF (pdf_allowed_i (is) .or.pdf_allowed_j (is) ) then 
         CALL indextocell (ia, iii, ks) 
         DO i = 1, 3 
         istart (i) = iii (i) - 1 - int (pdf_rmax / cr_a0 (i) ) 
         iend (i) = iii (i) + 1 + int (pdf_rmax / cr_a0 (i) ) 
         IF (pdf_2d.and.cr_icc (i) .le.1) then 
            istart (i) = iii (i) 
            iend (i) = iii (i) 
         ENDIF 
         ENDDO 
!                                                                       
!------ - In case we do not want to loop around cut range               
!                                                                       
         IF ( (.not.chem_period (1) ) .and. (.not.chem_period (2) )     &
         .and. (.not.chem_period (3) ) ) then                           
            DO i = 1, 3 
            istart (i) = max (1, istart (i) ) 
            iend (i) = min (cr_icc (i), iend (i) ) 
            ENDDO 
         ENDIF 
!                                                                       
!------ - Here starts the inner loop                                    
!                                                                       
         DO k = istart (3), iend (3) 
         DO j = istart (2), iend (2) 
         DO i = istart (1), iend (1) 
         cell (1) = i 
         cell (2) = j 
         cell (3) = k 
!                                                                       
!---------- Here we apply periodic boundaries                           
!---------- Modified to allow more than a single cycle !                
!                                                                       
         DO ii = 1, 3 
         iii (ii) = cell (ii) 
         cell (ii) = pdf_bnd (ii, cell (ii) ) 
         offset (ii) = float (iii (ii) - cell (ii) ) 
         ENDDO 
!                                                                       
!------ --- Now look for neighbours in surrounding unit cells only      
!                                                                       
         DO ii = 1, cr_ncatoms 
         CALL celltoindex (cell, ii, iatom) 
         js = cr_iscat (iatom) 
         IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or. (      &
         pdf_allowed_j (is) .and.pdf_allowed_i (js) ) ) then            
            DO jj = 1, 3 
            dd (jj) = cr_pos (jj, ia) - cr_pos (jj, iatom) - offset (jj) 
            d (jj) = abs (dd (jj) ) * cr_a0 (jj) 
            ENDDO 
!             dist2 = skalpro(dd,dd,cr_gten)                            
            dist2 = dd (1) * dd (1) * cr_gten (1, 1) + dd (2) * dd (2)  &
            * cr_gten (2, 2) + dd (3) * dd (3) * cr_gten (3, 3) + 2. *  &
            dd (1) * dd (2) * cr_gten (1, 2) + 2. * dd (1) * dd (3)     &
            * cr_gten (1, 3) + 2. * dd (2) * dd (3) * cr_gten (2, 3)    
            dist = sqrt (dist2) 
!                                                                       
            IF (dist.le.pdf_rmax.and.dist.gt.pdf_deltar) then 
               ibin = nint (dist / pdf_deltar) 
               pdf_temp (ibin, is, js) = pdf_temp (ibin, is, js)        &
               + 1                                                      
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE pdf_addcorr_n                  
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_e (ia) 
!+                                                                      
!     Calculate correlation for given atom ia, exact version            
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER is, js, ia, jj, iatom, ibin 
      REAL dist, dist2 
      REAL dd (3), d (3) 
!                                                                       
      is = cr_iscat (ia) 
      IF (pdf_allowed_i (is) .or.pdf_allowed_j (is) ) then 
!                                                                       
!------ - Here starts the inner loop over all atoms                     
!                                                                       
         DO iatom = 1, cr_natoms 
         js = cr_iscat (iatom) 
         IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or. (      &
         pdf_allowed_j (is) .and.pdf_allowed_i (js) ) ) then            
            DO jj = 1, 3 
            dd (jj) = cr_pos (jj, ia) - cr_pos (jj, iatom) 
            d (jj) = abs (dd (jj) ) * cr_a0 (jj) 
            ENDDO 
            dist2 = dd (1) * dd (1) * cr_gten (1, 1) + dd (2) * dd (2)  &
            * cr_gten (2, 2) + dd (3) * dd (3) * cr_gten (3, 3) + 2. *  &
            dd (1) * dd (2) * cr_gten (1, 2) + 2. * dd (1) * dd (3)     &
            * cr_gten (1, 3) + 2. * dd (2) * dd (3) * cr_gten (2, 3)    
            dist = sqrt (dist2) 
!                                                                       
            IF (dist.le.pdf_rmax.and.dist.gt.pdf_deltar) then 
               ibin = nint (dist / pdf_deltar) 
               pdf_temp (ibin, is, js) = pdf_temp (ibin, is, js)        &
               + 1                                                      
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE pdf_addcorr_e                  
!*****7*****************************************************************
      SUBROUTINE pdf_convtherm (rsign, sum) 
!+                                                                      
!     Convolute the pair correlation histograms with the                
!     thermal Gaussian                                                  
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'wink.inc' 
!                                                                       
      INTEGER ig, igaus, ib, ie 
      INTEGER ii, kk, is, js, ibin 
!     REAL(8) gaus ( - MAXDAT:MAXDAT) 
      REAL(8), DIMENSION(- PDF_MAXDAT:PDF_MAXDAT) :: gaus ! ( - MAXDAT:MAXDAT) 
      REAL rsign, sum 
      REAL asym, gnorm, dist, dist2, rg 
      REAL sigma, fac 
!                                                                       
      fac = 1.0 / (2.0 * zpi**2) 
      DO is = 1, cr_nscat 
      DO js = 1, cr_nscat 
      IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or. (         &
      pdf_allowed_j (is) .and.pdf_allowed_i (js) ) ) then               
                                                                        
         ii = int (pdf_rmax / pdf_deltar) + 1 
         DO ibin = 1, ii 
         dist = ibin * pdf_deltar 
         dist2 = dist**2 
!                                                                       
!------ --------- Convolute with Gaussian                               
!                                                                       
         IF (pdf_gauss.or.pdf_qalp.gt.0) then 
            sigma = 0.0 
            IF (pdf_gauss) then 
               sigma = fac * (cr_dw (is) + cr_dw (js) ) 
               sigma = sigma - pdf_delta / dist2 
               sigma = sigma - pdf_gamma / dist 
               sigma = max (0.0, sigma) 
            ENDIF 
            sigma = sigma + pdf_qalp**2 * dist2 
            sigma = sqrt (sigma) 
!                                                                       
            IF (dist.le.pdf_rcut) then 
               sigma = sigma * pdf_srat 
            ENDIF 
!                                                                       
            igaus = 1 + nint (5.0 * sigma / pdf_deltar) 
            ib = max (1, ibin - igaus + 1) 
            ie = min (pdf_bin, ibin + igaus - 1) 
!                                                                       
            IF (sigma.le.0.0.or.igaus.lt.2) then 
               pdf_corr (ibin) = pdf_corr (ibin) + pdf_temp (ibin, is,  &
               js) * rsign * pdf_weight (is, js) / pdf_deltar           
            ELSE 
               gnorm = 1.0 / (sqrt (zpi) * sigma) 
!                                                                       
               DO ig = - igaus, igaus 
               rg = (ig - 1) * pdf_deltar 
               asym = 1.0 + rg / dist 
               gaus (ig) = gnorm * asym * exp (-0.5*(rg/sigma)** 2)
               ENDDO 
!                                                                       
               DO ig = ib, ie 
               kk = ig - ibin + 1 
               pdf_corr (ig) = pdf_corr (ig) + pdf_temp (ibin, is, js)  &
               * rsign * pdf_weight (is, js) * gaus (kk)                
               ENDDO 
            ENDIF 
!                                                                       
!------ --------- Just take the distance                                
!                                                                       
         ELSE 
            pdf_corr (ibin) = pdf_corr (ibin) + pdf_temp (ibin, is, js) &
            * rsign * pdf_weight (is, js)                               
         ENDIF 
         sum = sum + pdf_weight (is, js) 
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE pdf_convtherm                  
!*****7*****************************************************************
      SUBROUTINE errlist_pdf 
!-                                                                      
!     Displays error Messages for the error type PDF                    
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = - 9, IO = 0) 
!                                                                       
      CHARACTER(41) ERROR (IU:IO) 
!                                                                       
      DATA ERROR ( -9:  0) /                              &
          'No atoms in asymmetric unit',                  & !  -r98
          'Atom type ALL not allowed              ',      & !  -8
          'Disable Gaussian mode and recalculate  ',      & !  -7
          'PDF range fixed with data loaded ',            & !  -6
          'PDF data must start with r=Dr          ',      & !  -5
          'No structure defined yet (>= 1 atoms)  ',      & !  -4
          'Crystal too large for peridic bound .   ',     & !  -3
          'Cannot extend r-range for corr. convol.',      & !  -2
          'Too many points in PDF                 ',      & !  -1
          ' '                                             & !   0
          /                                  
!                                                                       
      CALL disp_error ('PDF ', error, iu, io) 
      END SUBROUTINE errlist_pdf                    
