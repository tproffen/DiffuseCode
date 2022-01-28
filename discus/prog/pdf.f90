MODULE pdf_menu
!
CONTAINS
!*****7*****************************************************************
SUBROUTINE pdf 
!-                                                                      
!     This sublevel contains all routines dealing with the              
!     PDF anlysis part of DISCUS. This segment uses variables           
!     of the RMC and CHEM segment which are simply overwritten !!       
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE modify_mod
      USE pdf_mod 
      USE rmc_mod 
!
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
      USE class_macro_internal
      USE get_params_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      USE string_convert_mod
      USE sup_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 24  ! A command requires at least these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER(LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      INTEGER        , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
!
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(LEN=PREC_STRING) :: line, zeile, cdummy 
      INTEGER lp, length
      INTEGER indxg, ianz, lbef
      INTEGER :: n_nscat = 1  ! Dummy for RMC allocation
      INTEGER :: n_nsite = 1  ! Dummy for RMC allocation
      LOGICAL ldummy 
!                                                                       
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
!
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/pdf' 
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
             line == char(13) .or. line(1:1) == '!'  ) THEN
            IF(linteractive .OR. lmakro) THEN
               GOTO 10
            ELSE
               RETURN
            ENDIF
         ENDIF
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0.and.                                      &
          .not. (str_comp (befehl, 'echo', 2, lbef, 4) ) .and. &
          .not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and. &
          .not. (str_comp (befehl, 'help', 2, lbef, 4)   .or.  &
                 str_comp (befehl, '?   ', 2, lbef, 4) ) ) then                                              
            CALL do_math (line, indxg, length) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            line(1:length-1) = line(2:length)
            line(length:length) = ' '
            length = length - 1
            CALL file_kdo(line, length)
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
            CALL do_eval (zeile, lp, .TRUE.) 
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
!------ Reset PDF segement                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'rese', 2, lbef, 4) ) THEN 
            pdf_ldata = .false. 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  IF(str_comp(cpara(1), 'all', 2, lpara(1), 3) ) THEN
                     CALL pdf_reset
                  ELSE 
                     ier_typ = ER_COMM 
                     ier_num = - 6 
                  ENDIF 
               ELSE 
                  ier_typ = ER_COMM 
                  ier_num = - 6 
               ENDIF 
            ENDIF 
        ELSE
!                                                                       
!------ ------------------------------------------------------------    
!------ Here start the PDF specific commands ...                        
!------ ------------------------------------------------------------    
!
           IF( cr_nscat > PDF_MAXSCAT .or. MAXSCAT > PDF_MAXSCAT) THEN
                      pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
                      pdf_nsite = MAX(pdf_nsite, cr_ncatoms, PDF_MAXSITE, PDF_MAXSCAT, MAXSCAT)
                      pdf_ndat  = MAX(pdf_ndat ,           PDF_MAXDAT)
                      pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
                      CALL alloc_pdf( pdf_nscat, pdf_nsite, pdf_ndat, pdf_nbnd )
              IF ( ier_num < 0 ) THEN
                 RETURN
              ENDIF
           ENDIF
!                                                                       
!------ Just calculate PDF for current parameters                       
!                                                                       
             IF (str_comp (befehl, 'calc', 2, lbef, 4) ) then 
            pdf_mode = PDF_DO_CALC
            pdf_success = .FALSE.
            CALL pdf_setup (pdf_mode)
            IF (ier_num.eq.0) then 
               pdf_skal = 1.0 / rmc_skal (1) 
               CALL pdf_determine (.true.) 
            ENDIF 
            IF(ier_num==0) pdf_success = .TRUE.
!                                                                       
!------ Read observed PDF from XY file (ASCII)                          
!                                                                       
         ELSEIF (str_comp (befehl, 'data', 2, lbef, 4) ) then 
            pdf_mode = PDF_DO_FIT   ! Data were read, assume fit mode
            CALL pdf_readdata (zeile, lp) 
!                                                                       
!------ Runs PDF fit                                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
            pdf_mode = PDF_DO_FIT
            CALL pdf_setup (pdf_mode )
            IF (ier_num.eq.0) then 
               CALL pdf_run 
               IF(ier_num==0) pdf_success = .TRUE.
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
            pdf_lsite_i, 0, PDF_MAXSITE,                   &
            ldummy, .false., str_comp (befehl,             &
            'isele', 3, lbef, 5) )                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'jsele', 3, lbef, 5) .or.        &
                 str_comp (befehl, 'jdese', 3, lbef, 5) ) then                            
!                                                                       
            CALL atom_select (zeile, lp, 0, MAXSCAT, pdf_allowed_j, &
            pdf_lsite_j, 0, PDF_MAXSITE,                   &
            ldummy, .false., str_comp (befehl,             &
            'jsele', 3, lbef, 5) )                                      
!                                                                       
!                                                                       
!------ selecting/deselecting atoms                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'sele', 3, lbef, 4) .or.         &
                 str_comp (befehl, 'dese', 2, lbef, 4) ) then                             
!                                                                       
!           This might be the first time RMC arrays are referenced
            IF(cr_nscat > RMC_MAXSCAT .or. MAXSCAT > RMC_MAXSCAT) THEN
               n_nscat = MAX(cr_nscat, MAXSCAT, RMC_MAXSCAT)
               n_nsite = MAX(cr_ncatoms, MAXSCAT, RMC_MAXSITE)
               CALL alloc_rmc ( n_nscat, n_nsite )
               IF ( ier_num < 0 ) THEN
                  RETURN
               ENDIF
            ENDIF
            CALL atom_select (zeile, lp, 0, MAXSCAT, rmc_allowed, &
            rmc_lsite  , 0, RMC_MAXSITE,                   &
            rmc_sel_atom, .false., str_comp (              &
            befehl, 'sele', 3, lbef, 4) )                               
!                                                                       
!------ selecting/deselecting of molecules                              
!                                                                       
         ELSEIF (str_comp (befehl, 'msel', 2, lbef, 4) .or.   &
                 str_comp (befehl, 'mdes', 2, lbef, 4) ) then                             
!                                                                       
            CALL mole_select (zeile, lp, 0, MAXSCAT, rmc_allowed, &
            rmc_sel_atom, str_comp (  &
            befehl, 'msel', 2, lbef, 4) )                               
!                                                                       
!------ no command found                                                
!                                                                       
         ELSE 
            ier_num = - 8 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
      ENDIF 
!                                                                       
!------ any errors ?                                                    
!                                                                       
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in pdf menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN
               ELSE
                  IF(lmacro_close) THEN
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (lblock) THEN 
               ier_num = - 11 
               ier_typ = ER_COMM 
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN 
            ENDIF 
            CALL no_error 
            lmakro_error = .FALSE.
            sprompt = ' '
         ENDIF 
      ENDIF 
      IF(linteractive .OR. lmakro) THEN
         GOTO 10
      ELSE
         RETURN
      ENDIF
!                                                                       
 9999 CONTINUE 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE pdf                            
!*****7*****************************************************************
      SUBROUTINE pdf_setup (mode)
!+                                                                      
!     Setup for various arrays and functions for PDF calculation.       
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE pdf_mod 
      USE wink_mod
!
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE precision_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: mode
       
!                                                                       
      REAL sincut, rcut, z, bave, hh, rtot, ract 
      REAL(PREC_DP) :: factor
      INTEGER :: max_bnd
      INTEGER i, j, ia, is, js, nn, nnn 
      INTEGER, DIMENSION(0:MAXSCAT) :: pdf_natoms
      LOGICAL ltot 
!                                                                       
!     REAL form, quad 
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
      IF(mode == PDF_DO_CALC) THEN ! Prepare for simple calculation
         pdf_rmax  = pdf_rmaxu
         pdf_rmin  = pdf_deltaru
         pdf_rfmin = pdf_deltaru
         pdf_rfmax = pdf_rmaxu
      ELSEIF(mode == PDF_DO_FIT) THEN ! Prepare for refinement
         IF(pdf_rfminu == 0.01) pdf_rfminu = pdf_rfminf  ! No user values, use file
         IF(pdf_rfmaxu == 0.01) pdf_rfmaxu = pdf_rfmaxf  ! No user values, use file
         IF(pdf_rfminu > pdf_rfmaxf .or. &               ! User values outside data range
            pdf_rfminu < pdf_rfminf      ) THEN
            ier_num =  -9
            ier_typ = ER_PDF
            RETURN
         ENDIF
         IF(pdf_rfmaxu < pdf_rfminf .or. &
            pdf_rfmaxu > pdf_rfmaxf      ) THEN         ! User values outside data range
            ier_num =  -10
            ier_typ = ER_PDF
            RETURN
         ENDIF
         pdf_rfmin = MAX(pdf_rfminu, pdf_rfminf)
         pdf_rfmax = MAX(pdf_rfmaxu, pdf_rfmaxf)
         pdf_rmin  = pdf_deltar
         pdf_rmax  = pdf_rfmax
      ENDIF
!
!------ Allocate arrays
!
      sincut   = 0.025
      rcut     = 1.0
      IF(pdf_qmax> 0) THEN
         rcut     = 1.0 / (pdf_qmax * sincut)
      ENDIF
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
        pdf_nsite = MAX(pdf_nsite, cr_ncatoms, PDF_MAXSITE, PDF_MAXSCAT, MAXSCAT)
        pdf_ndat  = MAX(pdf_ndat ,           PDF_MAXDAT)
        pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
        CALL alloc_pdf( pdf_nscat, pdf_nsite, pdf_ndat, pdf_nbnd )
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
      CALL dlink (ano, lambda, rlambda,  renergy, l_energy, &
                  pdf_radiation, RAD_WAAS, pdf_power) 
      bave = 0.0 
      hh = pdf_xq**2 
!                                                                       
!RBN  cr_n_real_atoms = 0 
!RBN  DO ia = 1, cr_natoms 
!RBN  IF (cr_at_lis (cr_iscat (ia) ) .ne.'VOID') then 
!RBN     bave = bave+form (cr_iscat (ia), cr_scat, pdf_lxray, hh,  pdf_power) *cr_occ(cr_iscat(ia))
!RBN     cr_n_real_atoms = cr_n_real_atoms + 1 
!RBN  ENDIF 
!RBN  ENDDO 
!RBN  bave = bave / REAL(cr_n_real_atoms) 
!
      pdf_natoms(:) = 0
      cr_n_real_atoms = 0 
      DO ia = 1, cr_natoms
         IF (cr_at_lis (cr_iscat (ia) ) .ne.'VOID') then 
            pdf_natoms(cr_iscat(ia)) = pdf_natoms(cr_iscat(ia)) + 1
         ENDIF 
      ENDDO 
!
      bave = 0.0
      DO is = 0, cr_nscat 
         pdf_natoms(is) = NINT( REAL(pdf_natoms(is)) * cr_occ(is))
         cr_n_real_atoms = cr_n_real_atoms + pdf_natoms(is)
         bave = bave+form (is, cr_scat, pdf_lxray, hh,  pdf_power) *pdf_natoms(is)
      ENDDO
      bave = bave / REAL(cr_n_real_atoms) 
!                                                                       
      IF (.not.pdf_lweights) then 
         DO is = 0, cr_nscat 
         DO js = 0, cr_nscat 
         pdf_weight (is, js) = form (is, cr_scat, pdf_lxray, hh, pdf_power)        &
         * form (js, cr_scat, pdf_lxray, hh, pdf_power) / bave**2                  &
         *cr_occ(is)*cr_occ(js)
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
!        pdf_sinc = 0.0
         j = SIZE(pdf_sincc)+1
!         DO i = 1, SIZE(pdf_sinc)/2 ! INT(nn *1.5)
!         pdf_sinc (i+1) = sin (z * REAL(i) ) / (pdf_deltar * REAL(i) ) 
!         ENDDO 
         DO i = 1,j/2
!        pdf_sincc(i)   = sin (z * REAL(i)*0.91 ) / (pdf_deltar * REAL(i)*0.91 ) 
!        pdf_sincc(j-i) = sin (z * REAL(i)*0.91 ) / (pdf_deltar * REAL(i)*0.91 ) 
         pdf_sincc(i+1) = sin (z * REAL(i)      ) / (pdf_deltar * REAL(i)      ) 
         pdf_sincc(j-i) = sin (z * REAL(i)      ) / (pdf_deltar * REAL(i)      ) 
         ENDDO 
         pdf_sincc(1) = pdf_qmax
!        pdf_sinc (1) = pdf_qmax
!        DO i = nn + 1, 2 * PDF_MAXDAT 
!        pdf_sinc (i) = 0.0 
!        ENDDO 
      ENDIF 
!
      IF(pdf_gauss.or.pdf_qalp>0.0) THEN
         IF(pdf_gauss_init) THEN
         factor = -0.5D0*pdf_gauss_step**2/1.D0/1.D0    ! 1/2 *r^2 /sigma^2
         j = UBOUND(pdf_exp,1)
         DO i=0,j !/2
            pdf_exp(i)      = exp(factor * REAL(i*i))
!           pdf_exp(j-i) = pdf_exp(i)
         ENDDO
           pdf_gauss_init = .false.
        ENDIF
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
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE chem_mod 
      USE pdf_mod 
      USE rmc_mod 
      USE rmc_sup_mod
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: cmd 
!                                                                       
      CHARACTER(9) at_lis (0:MAXSCAT+1) 
      CHARACTER(9) cpoly (5) 
      INTEGER i, j, k 
!                                                                       
      DATA cpoly / 'linear  ', 'square  ', '3. order', '4. order', &
                   '5. order' /                                                           
!                                                                       
      IF (cmd.eq.'ALL'.or.cmd.eq.'PDF') then 
         CALL pdf_setup (pdf_mode)
         IF (ier_num.ne.0) return 
!                                                                       
         pdf_skal = 1.0 / rmc_skal (1) 
!                                                                       
         WRITE (output_io, 1000) 
         WRITE (output_io, 2000) pdf_rmax, pdf_deltaru, pdf_bin 
!                                                                       
         SELECTCASE(pdf_radiation)
            CASE(PDF_RAD_XRAY)
               WRITE (output_io, 2010) 'X-Rays', pdf_xq 
            CASE(PDF_RAD_NEUT)
               WRITE (output_io, 2015) 'Neutrons' 
            CASE(PDF_RAD_ELEC)
               WRITE (output_io, 2010) 'Electrons', pdf_xq 
         END SELECT
!        IF (pdf_lxray) then 
!           WRITE (output_io, 2010) 'X-Rays', pdf_xq 
!        ELSE 
!           WRITE (output_io, 2015) 'Neutrons' 
!        ENDIF 
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
            WRITE (output_io, 2410) pdf_cquad_a 
            WRITE (output_io, 2415) pdf_clin_a 
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
         WRITE (output_io, 3010) pdf_rfmin, int (pdf_rfmin / pdf_deltaru) 
         WRITE (output_io, 3020) pdf_rfmax, int (pdf_rfmax / pdf_deltaru) 
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
      USE discus_config_mod 
      USE crystal_mod 
      USE pdf_mod 
      USE save_menu
!
      USE debug_mod 
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
USE precision_mod
      USE string_convert_mod
USE support_mod
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: zeile 
      INTEGER          , INTENT(INOUT) :: lp
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      INTEGER i, nmi, nma, nmd 
      INTEGER ::  pdf_calc_l, pdf_calc_u
      REAL r 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))) :: cdummy, cpara (maxw) 
      INTEGER ianz, lpara (maxw)
      REAL(KIND=PREC_DP) :: werte (maxw) 
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
            IF(pdf_success) THEN                                        
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            WRITE (output_io, 1505) cpara (1) (1:lpara (1) ) 
            cdummy = cpara (1) (1:lpara (1) ) 
            pdf_calc_l = LBOUND(pdf_calc,1)
            pdf_calc_u = UBOUND(pdf_calc,1)
            CALL pdf_save_file(cdummy, pdf_rfmin, pdf_rfmax, pdf_deltar, &
                 pdf_us_int, pdf_calc_l, pdf_calc_u, pdf_skal, pdf_calc)
            ELSE
               ier_num = -12
               ier_typ = ER_PDF
               RETURN
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
               CALL oeffne (57, cdummy, 'unknown') 
               IF (ier_num.eq.0) then 
                  nmi = int (pdf_rfmin / pdf_deltar) 
                  nma = int (pdf_rfmax / pdf_deltar) 
                  nmd = pdf_us_int   ! step width = (delta r user)/(deltar internal)
                  DO i = nmi, nma , nmd
                  r = REAL(i) * pdf_deltar 
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
! 5000 FORMAT     (F9.4,3X,F21.10,5X,2(F6.2,1X)) 
 5100 FORMAT     (F9.4,3X,F6.2) 
      END SUBROUTINE pdf_save                       
!*****7*****************************************************************
      SUBROUTINE pdf_readdata (zeile, lp) 
!+                                                                      
!     Reads observed PDF as xy ASCII file.                              
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE pdf_mod 
!
      USE debug_mod 
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
USE support_mod
      USE ISO_FORTRAN_ENV
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: zeile 
      INTEGER          , INTENT(INOUT) :: lp
!                                                                       
      INTEGER, PARAMETER :: maxw = 5 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW) :: cpara
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile)))                    :: datafile 
      INTEGER            , DIMENSION(1:MAXW) :: lpara (maxw) 
      INTEGER ianz, ip 
      INTEGER  :: iostatus
      INTEGER  :: n_dat
      REAL ra, re, dr 
      REAL                                   :: r_dummy1, r_dummy2
      REAL(KIND=PREC_DP) , DIMENSION(1:MAXW) :: werte
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
      datafile = cpara(1)(1:lpara(1)) 
!                                                                       
!------ Read observed PDF for given plane                               
!                                                                       
      IF (ianz.eq.1) then 
         CALL oeffne (17, datafile, 'old') 
         IF (ier_num.ne.0) return 
         pdf_obs(:) = 0.0    ! reset data file
         pdf_wic(:) = 0.0    ! reset weighting scheme
         n_dat      = 0
         READ (17, *, IOSTAT = iostatus  ) ra, r_dummy1    , dr, r_dummy2
preread: DO
            IF(IS_IOSTAT_END(iostatus)) EXIT preread
            n_dat = n_dat + 1
            READ (17, *, IOSTAT = iostatus  ) ra, r_dummy1   , dr, r_dummy2
         ENDDO preread
         REWIND(17)
         IF(n_dat > PDF_MAXDAT) THEN
            pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
            pdf_nsite = MAX(pdf_nsite, cr_ncatoms, PDF_MAXSITE, PDF_MAXSCAT, MAXSCAT)
            pdf_ndat  = MAX(pdf_ndat , n_dat   , PDF_MAXDAT)
            pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
            CALL alloc_pdf( pdf_nscat, pdf_nsite, pdf_ndat, pdf_nbnd )
            IF ( ier_num < 0 ) THEN
               RETURN
            ENDIF
         ENDIF
         CALL extract_hist (17) 
         CALL skip_spec (17) 
         ip = 1 
!
!        READ (17, *, end = 20, err = 999) ra, pdf_obs (ip), dr, pdf_wic (ip)
         READ (17, *, IOSTAT = iostatus  ) ra, pdf_obs (ip), dr, pdf_wic (ip)
         IF(IS_IOSTAT_END(iostatus)) THEN
            CLOSE(17)
         ELSE IF(IS_IOSTAT_EOR(iostatus)) THEN
            CLOSE(17)
            ier_num = - 3
            ier_typ = ER_IO                                                           
            RETURN
         ENDIF
         IF(pdf_wic(ip) /= 0.0) THEN
            pdf_wic(ip) = 1./pdf_wic(ip)**2
         ELSE
            pdf_wic(ip) = 1.0
         ENDIF
!!!         ip = ip + 1 
!!!   10    CONTINUE 
!!!         READ (17, *, end = 20, err = 999) re, pdf_obs (ip), dr, pdf_wic (ip)
         ier_num = 0
         ier_typ = ER_NONE
main:    DO
            ip = ip + 1 
            IF (ip.gt.PDF_MAXDAT) goto 9999 
            READ (17, *, IOSTAT = iostatus  ) re, pdf_obs (ip), dr, pdf_wic (ip)
            IF(IS_IOSTAT_END(iostatus)) THEN
               EXIT main
            ELSE IF(IS_IOSTAT_EOR(iostatus)) THEN
               ier_num = - 3
               ier_typ = ER_IO                                                           
               EXIT main
            ENDIF
            IF(pdf_wic(ip) /= 0.0) THEN
               pdf_wic(ip) = 1./pdf_wic(ip)**2
            ELSE
               pdf_wic(ip) = 1.0
            ENDIF
         END DO main
!!!         ip = ip + 1 
!!!         IF (ip.gt.PDF_MAXDAT) goto 9999 
!!!         GOTO 10 
!!!   20    CONTINUE 
         CLOSE (17) 
         IF (ier_num.ne.0) return 
!!!         pdf_rmax = re 
!!!         pdf_bin = ip - 1 
!!!         pdf_deltar = (re-ra) / REAL(pdf_bin - 1) 
!!!         pdf_deltaru= pdf_deltar    ! copy delta r from file into user 
!!!         pdf_us_int = 1
         pdf_deltaru =  (re-ra)/REAL(ip-2)
         IF(pdf_deltaru > pdf_deltari) THEN
            pdf_us_int = NINT(MAX(pdf_deltaru/pdf_deltari,1.))
            pdf_deltar = pdf_deltaru/pdf_us_int ! internal delta R ~ pdf_deltari always
         ELSE
            pdf_deltar = pdf_deltaru
            pdf_us_int = 1
         ENDIF
         pdf_bin  = NINT((re-ra)/pdf_deltar)
         pdf_rmax = re
         pdf_deltars=pdf_deltar/2.
         pdf_rfminf = ra      ! Set limits from file
         pdf_rfmaxf = re      ! final value will be set in pdf_setup
!                                                                       
         IF (abs (ra - pdf_deltaru) .gt.1e-6) then 
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
!!!  999 CONTINUE 
!!!      CLOSE (17) 
!!!      ier_num = - 3 
!!!      ier_typ = ER_IO 
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
      USE discus_config_mod 
      USE pdf_mod 
!
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER , INTENT(IN) :: ifil 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line 
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
      INTEGER , INTENT(IN) :: ifil 
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
      USE discus_config_mod 
USE lib_length
      USE param_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: line
      CHARACTER (LEN=*), INTENT(IN) :: key 
!
      INTEGER is, ie, ll, lk 
!                                                                       
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
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE get_iscat_mod
      USE modify_mod
      USE pdf_mod 
      USE rmc_mod 
      USE rmc_sup_mod
!
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
      INTEGER          , INTENT(INOUT) :: lp 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))) :: cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw)
REAL(KIND=PREC_DP) ::  wa (maxw), wb (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, nn 
      INTEGER i, j, ia, ib 
      INTEGER  :: n_nscat = 1! dummy for rmc_allocation
      INTEGER  :: n_nsite = 1! dummy for rmc_allocation
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.ge.2) then 
         CALL do_cap (cpara (1) ) 
!                                                                       
!------ - set boundary: toggle periodic boundaries                      
!                                                                       
         IF (str_comp(cpara (1),'BOUNDARY',3,lpara(1),8)) THEN
            IF (ianz==2 .OR. ianz==3 .OR. ianz==4) THEN 
               IF(str_comp(cpara(2), 'periodic', 3, lpara(2), 8) .AND. chem_purge) THEN
                  ier_num = -31
                  ier_typ = ER_CHEM
                  ier_msg(1) = "Use >set crystal< in chem to define "
                  ier_msg(2) = "Number of unit cells and atoms per unit cell"
                  ier_msg(3) = "Or read a new cell/structure"
                  RETURN
               ELSE
               chem_period (1) = str_comp(cpara(2), 'periodic', 3, lpara(2), 8)
               chem_period (2) = chem_period (1) 
               chem_period (3) = chem_period (1) 
               IF (ianz.eq.3) THEN 
                  cpara(4) = ' '
                  lpara(4) =  1
               ENDIF
               pdf_2d = str_comp (cpara (3),'2D',1,lpara(3), 2) .OR. &
                        str_comp (cpara (4),'2D',1,lpara(4), 2)
               pdf_lexact = str_comp (cpara (3),'exact',1,lpara(3), 2) .OR. &
                            str_comp (cpara (4),'exact',1,lpara(4), 2)
               IF (chem_period(1) ) THEN
                  IF (pdf_2d) THEN 
                     IF (pdf_lexact) THEN 
                       WRITE (output_io, 1000) 'periodic bound. 2D, exact ' 
                     ELSE 
                       WRITE (output_io, 1000) 'periodic bound. 2D, unit cell' 
                     ENDIF 
                  ELSE 
                     IF (pdf_lexact) THEN 
                       WRITE (output_io, 1000) 'periodic bound. 3D, exact ' 
                     ELSE 
                       WRITE (output_io, 1000) 'periodic bound. 3D, unit cell' 
                     ENDIF 
                  ENDIF 
               ELSE 
                  IF (pdf_2d) THEN 
                     IF (pdf_lexact) THEN 
                        WRITE (output_io, 1000) 'no periodic bound., 2D, exact' 
                     ELSE 
                        WRITE (output_io, 1000) 'no periodic bound., 2D, unit cell' 
                     ENDIF 
                  ELSE 
                     IF (pdf_lexact) THEN 
                        WRITE (output_io, 1000) 'no periodic bound., 3D, exact' 
                     ELSE 
                        WRITE (output_io, 1000) 'no periodic bound., 3D, unit cell' 
                     ENDIF 
                  ENDIF 
               ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set mode: toggle calculation mode                             
!                                                                       
         ELSEIF (str_comp(cpara (1),'CALC',3,lpara(1),4)) THEN
            IF (ianz.eq.2) then 
               pdf_lexact = str_comp (cpara (2) , 'exact', 2, lpara (2),5)
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set finite: finite size correction model                      
!                                                                       
         ELSEIF (str_comp(cpara (1),'FINITE',3,lpara(1),6)) THEN
            IF (ianz.ge.2) then 
               IF (str_comp (cpara (2) , 'periodic', 3, lpara (2) , 8) ) then
                  pdf_finite = PDF_BACK_PERIOD 
               ELSEIF (str_comp (cpara (2) ,'sphere',3,lpara(2),6)) then
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
               ELSEIF (str_comp (cpara (2),'polygon',3,lpara(2),7)) then
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
         ELSEIF (str_comp(cpara (1),'QALPHA',3,lpara(1),6) .OR. &
                 str_comp(cpara (1),'QBROAD',3,lpara(1),6)) then 
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
         ELSEIF (str_comp(cpara (1),'CORRLINEAR',5,lpara(1),10) .OR. &
                 str_comp(cpara (1),'GAMMA',3,lpara(1),5)) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_clin_a = werte (1) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set delta: sets quadratic correlation factor                  
!                                                                       
         ELSEIF (str_comp(cpara (1),'CORRQUADRATIC',5,lpara(1),13) .OR. &
                 str_comp(cpara (1),'DELTA',3,lpara(1),5)) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1) then 
               pdf_cquad_a = werte (1) 
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
         ELSEIF (str_comp(cpara (1),'DIAMETER',3,lpara(1),8)) THEN
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
         ELSEIF (str_comp(cpara (1),'FRANGE',3,lpara(1),6)) THEN
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.2) then 
               pdf_rfminu = werte (1) ! user value internal to be set 
               pdf_rfmaxu = werte (2) ! with actual data
!              pdf_rfmin  = werte (1) 
!              pdf_rfmax  = werte (2) 
               pdf_mode   = PDF_DO_FIT ! assume fit mode
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ - set qmax: sets Qmax for termination correction                
!                                                                       
         ELSEIF (str_comp(cpara (1),'QMAX',3,lpara(1),4)) THEN
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
         ELSEIF (str_comp(cpara (1),'QSIGMA',3,lpara(1),6) .OR. &
                 str_comp(cpara (1),'QDAMP' ,3,lpara(1),5)) then 
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
         ELSEIF (str_comp(cpara (1),'PARTIAL',3,lpara(1),7)) THEN
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
         ELSEIF (str_comp(cpara (1),'POLYGON',3,lpara(1),7)) THEN
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
         ELSEIF (str_comp(cpara (1),'RADIATION',3,lpara(1),9)) THEN
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) return 
            IF (ianz.eq.1.or.ianz.eq.2) then 
               CALL do_cap (cpara (1) ) 
               IF (cpara (1) (1:1) .eq.'N') then 
                  pdf_lxray = .false. 
                  pdf_radiation = PDF_RAD_NEUT
               ELSEIF (cpara (1) (1:1) .eq.'X') then 
                  pdf_lxray = .true. 
                  pdf_radiation = PDF_RAD_XRAY
                  IF (ianz.eq.2) then 
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.ne.0) return 
                     pdf_xq = werte (1) 
                  ELSE 
                     pdf_xq = 0.0 
                  ENDIF 
               ELSEIF (cpara (1) (1:1) .eq.'E') then 
                  pdf_lxray = .true. 
                  pdf_radiation = PDF_RAD_ELEC
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
         ELSEIF (str_comp(cpara (1),'RANGE',3,lpara(1),5)) THEN
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
               ENDIF 
               pdf_mode   = PDF_DO_CALC ! assume simple calculation mode
               pdf_rmaxu  = werte (1)   ! user supplied Rmax
               pdf_rmax   = werte (1) 
!              pdf_deltar = werte (2) 
               pdf_deltaru= werte (2)   ! user supplied delta R
               IF(pdf_deltaru>pdf_deltari) then
                  pdf_us_int = NINT(MAX(pdf_deltaru/pdf_deltari,1.))
                  pdf_deltar = pdf_deltaru/pdf_us_int ! internal delta R ~ pdf_deltari always
               ELSE
                  pdf_deltar = pdf_deltaru
                  pdf_us_int = 1
               ENDIF
               pdf_deltars= pdf_deltar/2.
               pdf_rfmin = pdf_deltaru  ! werte (2) 
               pdf_rfmax = pdf_rmax 
!               ELSEIF (nn.gt.PDF_MAXDAT) then
               nn = int (werte (1) / pdf_deltar) + 1
                 pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
                 pdf_ndat  = MAX(pdf_ndat , nn      , PDF_MAXDAT)
                 pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
!                There is no need to call alloc at this point, is done in  setup
!                CALL alloc_pdf( pdf_nscat, pdf_ndat, pdf_nbnd )
               pdf_bin = nn 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!
!------ - set refinement flags
!
         ELSEIF (str_comp(cpara (1),'REFINE',3,lpara(1),6)) THEN
!                                                                       
!------ - set weight: sets scale parameter to correct the weighting     
!                                                                       
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (str_comp(cpara (1),'none',3,lpara(1),4)) THEN
               pdf_refine_scale   = .false.
               pdf_refine_density = .false.
               pdf_refine_lattice = .false.
            ELSEIF (str_comp(cpara (1),'scale',3,lpara(1),5)) THEN
               pdf_refine_scale   = .true.
            ELSEIF (str_comp(cpara (1),'density',3,lpara(1),7)) THEN
               pdf_refine_density = .true.
            ELSEIF (str_comp(cpara (1),'lattice',3,lpara(1),7)) THEN
               IF (ianz == 1) then 
                  IF(cr_syst == CR_CUBIC) THEN
                     pdf_refine_lattice    = .false.
                     pdf_refine_lattice(1) = .true.
                  ELSEIF(cr_syst == CR_HEXAGONAL .or. &
                         cr_syst == CR_TRIGONAL  .or. &
                         cr_syst == CR_TETRAGONAL    ) THEN
                     pdf_refine_lattice    = .false.
                     pdf_refine_lattice(1) = .true.
                     pdf_refine_lattice(3) = .true.
                  ELSEIF(cr_syst == CR_RHOMBOHED) THEN
                     pdf_refine_lattice    = .false.
                     pdf_refine_lattice(1) = .true.
                     pdf_refine_lattice(4) = .true.
                  ELSEIF(cr_syst == CR_ORTHO) THEN
                     pdf_refine_lattice(4:6) = .false.
                     pdf_refine_lattice(1:3) = .true.
                  ELSEIF(cr_syst == CR_MONOCLINICC) THEN
                     pdf_refine_lattice(1:3) = .true.
                     pdf_refine_lattice(4:5) = .false.
                     pdf_refine_lattice(6)   = .true.
                  ELSEIF(cr_syst == CR_MONOCLINICB) THEN
                     pdf_refine_lattice(1:3) = .true.
                     pdf_refine_lattice(4:6) = .false.
                     pdf_refine_lattice(5)   = .true.
                  ELSEIF(cr_syst == CR_TRICLINIC  ) THEN
                     pdf_refine_lattice      = .true.
                  ENDIF
               ELSE
                  IF (str_comp(cpara (2),'a',1,lpara(2),1) .and.    &
                      lpara(2)==1                             ) THEN
                     pdf_refine_lattice(1)   = .true.
                  ELSEIF (str_comp(cpara (2),'b',1,lpara(2),1).and. &
                      lpara(2)==1                             ) THEN
                     IF(cr_syst < CR_TETRAGONAL) THEN
                        pdf_refine_lattice(2)   = .true.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF
                  ELSEIF (str_comp(cpara (2),'c',1,lpara(2),1)) THEN
                     IF(cr_syst < CR_TETRAGONAL) THEN
                        pdf_refine_lattice(3)   = .true.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF
                  ELSEIF (str_comp(cpara (2),'alpha',2,lpara(2),5)) THEN
                     IF(cr_syst == CR_TRICLINIC .or. &
                        cr_syst == CR_RHOMBOHED      ) THEN
                        pdf_refine_lattice(4)   = .true.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF
                  ELSEIF (str_comp(cpara (2),'beta',2,lpara(2),4)) THEN
                     IF(cr_syst == CR_TRICLINIC .or. &
                        cr_syst == CR_MONOCLINICB    ) THEN
                        pdf_refine_lattice(5)   = .true.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF
                  ELSEIF (str_comp(cpara (2),'gamma',2,lpara(2),5)) THEN
                     IF(cr_syst == CR_TRICLINIC .or. &
                        cr_syst == CR_MONOCLINICC    ) THEN
                        pdf_refine_lattice(6)   = .true.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
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
         ELSEIF (str_comp(cpara (1),'WEIGHT',3,lpara(1),6)) THEN
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
         ELSEIF (str_comp(cpara (1),'SHAPE',3,lpara(1),5)) THEN
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
         ELSEIF (str_comp(cpara (1),'SRATIO',3,lpara(1),6)) THEN
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
         ELSEIF (str_comp(cpara (1),'THERMAL',3,lpara(1),7)) THEN
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
!           This might be the first time RMC arrays are referenced
            IF(cr_nscat > RMC_MAXSCAT .or. MAXSCAT > RMC_MAXSCAT) THEN
               n_nscat = MAX(cr_nscat, MAXSCAT, RMC_MAXSCAT)
               n_nsite = MAX(cr_ncatoms, MAXSCAT, RMC_MAXSITE)
               CALL alloc_rmc ( n_nscat, n_nsite )
               IF ( ier_num < 0 ) THEN
                  RETURN
               ENDIF
            ENDIF
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
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_aver_mod
      USE diffuse_mod 
      USE pdf_mod 
      USE refine_mod
      USE rmc_mod 
      USE rmc_sup_mod
      USE update_cr_dim_mod
      USE random_mod
!                                                                       
      USE errlist_mod 
USE lib_random_func
      USE param_mod 
      USE prompt_mod 
      USE precision_mod
USE support_mod
      IMPLICIT none 
       
!                                                                       
      REAL(PREC_DP) cc, c, ce, e, ee, wtot, cold, cnew
!     REAL pdf_old (MAXDAT) 
      REAL(PREC_DP), DIMENSION(PDF_MAXDAT) ::  pdf_old !  (MAXDAT) 
      REAL sig2, sumbad 
      REAL prob, psum, p2sum, pave, psig, pmax, pn 
      REAL start, zeit
      REAL p_new (3, rmc_max_atom) 
      REAL p_old (3, rmc_max_atom) 
      INTEGER i_new (rmc_max_atom) 
      INTEGER i_old (rmc_max_atom) 
      INTEGER isel (rmc_max_atom), natoms 
      INTEGER imol (rmc_max_atom) 
      INTEGER zh, zm, zs, nmi, nma 
      INTEGER i, j, ip, ipc
      INTEGER igen, itry, iacc_good, iacc_bad 
      LOGICAL loop, laccept 
!
!                                                                       
!open(89,file='shift.log', status='unknown')
      IF (pdf_obs (1) .eq. - 9999.) then 
         ier_num = - 4 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
!
      CALL refine_alloc     ! Allocate generic refinement settings
!                                                                       
      igen      = 0 
      itry      = 0 
      iacc_good = 0 
      iacc_bad  = 0 
      loop      = .true. 
      laccept   = .true. 
!                                                                       
      psum  = 0.0 
      p2sum = 0.0 
      pmax  = 0.0 
      pn    = 0.0 
!                                                                       
      nmi = nint (pdf_rfmin / pdf_deltar) 
      nma = nint (pdf_rfmax / pdf_deltar) 
!                                                                       
      IF(rmc_move_prob(RMC_MODE_SWDISP)> 0.0) THEN
         CALL chem_aver (.false., .true.) 
      ENDIF
!                                                                       
!------ Calculate maximal value in rmc_mindist array                    
!                                                                       
      rmc_mindist_max = rmc_mindist (1, 1) 
      DO i = 1, cr_nscat 
         DO j = 1, cr_nscat 
            rmc_mindist_max = MAX(rmc_mindist_max,rmc_mindist(i,j)) 
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
      cold = 0.0d0
      wtot = 0.0d0 
      e    = 0.0d0
      ee   = 0.0d0 
      c    = 0.0d0 
      cc   = 0.0d0 
      ce   = 0.0d0 
!                                                                       
      DO ip = nmi, nma / pdf_us_int
         ipc = (ip-1)*pdf_us_int + 1
         wtot = wtot + pdf_wic (ip)
         e    = e    + pdf_wic (ip) * pdf_obs (ip) 
         ee   = ee   + pdf_wic (ip) * pdf_obs (ip) **2 
         c    = c    + pdf_wic (ip) * pdf_calc (ipc) 
         cc   = cc   + pdf_wic (ip) * pdf_calc (ipc) **2 
         ce   = ce   + pdf_wic (ip) * pdf_calc (ipc) * pdf_obs (ip) 
      ENDDO 
      IF (rmc_doskal) then 
         pdf_skal = REAL(ce / cc )
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
      elseIF (rmc_mode.eq.rmc_mode_rotate) then 
         WRITE (output_io, 2000) 'Mode ROTATE' 
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
main: DO while (loop) 
   IF(ier_ctrlc) THEN
      ier_num = -14
      ier_typ = ER_COMM
      RETURN
   ENDIF
   IF(ier_num/=0) RETURN      ! An error occured or CTRL-C
   laccept = .true. 
   igen = igen + 1 
   IF((igen<0.05*rmc_maxcyc.and.MOD(igen,10)==2) .or. MOD(igen, 10)==2) THEN
         CALL pdf_rmc_scale( nmi, nma, cold, wtot, ee, psum, p2sum, &
                          sig2, pmax, pn )
   ELSEIF((igen<0.05*rmc_maxcyc.and.MOD(igen,10)==5) .or. MOD(igen,200)==5) THEN
         CALL pdf_rmc_lattice( nmi, nma, cold, wtot, ee, psum, p2sum, &
                               sig2, pmax, pn )
   ELSE
!                                                                       
!-------- generate move and check for limits                            
!                                                                       
   IF (rmc_sel_atom) then 
      CALL rmc_genmove (laccept, natoms, p_new, i_new, isel) 
   ELSE 
      CALL rmc_genmove_mol(laccept, natoms, p_new, i_new, isel, imol)
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
               CALL pdf_addcorr (isel (i), - 2.0, sumbad) 
            ENDDO 
            CALL pdf_makemove (natoms, i_new, p_new, isel, imol) 
            DO i = 1, natoms 
               CALL pdf_addcorr (isel (i), 2.0, sumbad) 
            ENDDO 
            CALL pdf_convert 
!                                                                       
            itry = itry + 1 
            cnew = 0.0d0 
            c    = 0.0 
            cc   = 0.0 
            ce   = 0.0 
!                                                                       
            DO ip = nmi, nma / pdf_us_int
               ipc = (ip-1)*pdf_us_int + 1
               c   = c  + pdf_wic (ip) * pdf_calc (ipc) 
               cc  = cc + pdf_wic (ip) * pdf_calc (ipc) **2 
               ce  = ce + pdf_wic (ip) * pdf_calc (ipc) * pdf_obs (ip) 
            ENDDO 
            IF (rmc_doskal) then 
               pdf_skal = REAL(ce / cc )
            ELSE 
               pdf_skal = 1.0 / rmc_skal (1) 
            ENDIF 
            cnew = ee+pdf_skal**2 * cc - 2.0 * pdf_skal * ce 
            cnew = cnew / wtot 
!                                                                       
!     ----Accept move ?                                                 
!                                                                       
            prob = REAL( cnew - cold )
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
!
!           CALL refine_adapt_move (natoms, isel, 1)
            ELSE 
               CALL pdf_makemove (natoms, i_old, p_old, isel, imol) 
               DO i = 1, pdf_bin 
                  pdf_corr (i) = pdf_old (i) 
               ENDDO 
               CALL pdf_convert 
!              CALL refine_adapt_move (natoms, isel, 0)
!
            ENDIF 
!write(89,7777) itry, rmc_maxmove(1,1), cold
!7777 format(i5,2(1x,G20.6e3))
         ENDIF
!                                                                       
!------ --WRITE info and terminate or loop again                        
!                                                                       
         ENDIF 
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
      ENDDO  main    ! end of main loop
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
      rmc_skal (1) = REAL(1.0d0 / pdf_skal )
!                                                                       
!------ save some results to res[i] blo                                 
!                                                                       
      res_para (0) = 8 
!                                                                       
      res_para (1) = REAL(cold) 
      res_para (2) = REAL(itry) 
      res_para (3) = REAL(iacc_good) 
      res_para (4) = REAL(iacc_bad) 
      res_para (5) = pave 
      res_para (6) = psig 
      res_para (7) = pmax 
      res_para (8) = zeit / REAL(itry) 
!                                                                       
!------ Update crystal dimensions                                       
!                                                                       
      CALL update_cr_dim 
!
      CALL refine_dealloc     ! Deallocate generic refinement settings
close(89)
!                                                                       
!------ formats                                                         
!                                                                       
 1000 FORMAT (' Running PDF-Fit ...',//                                 &
     &        ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,     &
     &        ' containing ',I9,' atoms')                               
 1250 FORMAT (/,' ---- Final configuration (s= ',F8.4,') ---- ',/) 
 1300 FORMAT (  ' Gen: ',I8,' try: ',I8,' acc: (good/bad): ',I8,        &
     &          ' / ',I8,' s2x2: ',G15.8)                               
 2000 FORMAT (/,' Starting main RMC loop ...',/,' (',A,') ',/) 
 3000 FORMAT (/,' Delta Chi    : ave: ',G15.4,' sig: ',G15.4,           &
     &          ' max: ',G15.4)                                         
 4000 FORMAT (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/     &
     &          ' Time/cycle   : ',F9.3,' sec',/)                       
      END SUBROUTINE pdf_run                        
!*****7*****************************************************************
      SUBROUTINE  pdf_rmc_scale( nmi, nma, cold, wtot, ee, psum, p2sum, &
                                    sig2, pmax, pn )

!
      USE pdf_mod
      USE refine_mod
      USE rmc_mod
USE lib_random_func
      USE random_mod
      USE precision_mod
!
      IMPLICIT NONE
!
      REAL(PREC_DP), INTENT(INOUT) :: cold
      REAL(PREC_DP), INTENT(IN)    :: wtot
      REAL(PREC_DP), INTENT(IN)    :: ee
      REAL    , INTENT(INOUT) :: psum
      REAL    , INTENT(INOUT) :: p2sum
      REAL    , INTENT(INOUT) :: pmax
      REAL    , INTENT(INOUT) :: pn 
      REAL    , INTENT(IN)    :: sig2
!
      INTEGER  :: i, ip, nmi, nma, ipc
      LOGICAL  :: laccept
      REAL(PREC_DP) :: c, cc, ce, cnew
      REAL     :: pdf_old_scale
      REAL     :: pdf_old_rho0 
      REAL     :: prob
      REAL(PREC_DP), DIMENSION(PDF_MAXDAT) ::  pdf_old !  (MAXDAT) 
!
!
      laccept = .false.
         pdf_old_scale = pdf_scale
         pdf_old_rho0  = pdf_rho0 
         IF(pdf_refine_scale )  THEN
            pdf_scale     = pdf_scale + gasdev(DBLE(ref_maxpdfsd(1)))
            laccept       = .true.
         ENDIF
         IF(pdf_refine_density) THEN
            pdf_rho0      = pdf_rho0  + gasdev(DBLE(ref_maxpdfsd(2)))
            laccept       = .true.
         ENDIF
      IF(laccept) THEN
         DO i = 1, pdf_bin 
            pdf_old (i) = pdf_corr (i) 
         ENDDO 
         CALL pdf_convert
         cnew = 0.0d0 
         c    = 0.0 
         cc   = 0.0 
         ce   = 0.0 
!                                                                       
         DO ip = nmi, nma/pdf_us_int
            ipc = (ip-1)*pdf_us_int + 1
            c   = c  + pdf_wic (ip) * pdf_calc (ipc) 
            cc  = cc + pdf_wic (ip) * pdf_calc (ipc) **2 
            ce  = ce + pdf_wic (ip) * pdf_calc (ipc) * pdf_obs (ip) 
         ENDDO 
         IF (rmc_doskal) then 
            pdf_skal = REAL(ce / cc )
         ELSE 
            pdf_skal = 1.0 / rmc_skal (1) 
         ENDIF 
         cnew = ee+pdf_skal**2 * cc - 2.0 * pdf_skal * ce 
         cnew = cnew / wtot 
!                                                                       
!     ----Accept move ?                                                 
!                                                                       
         prob = REAL( cnew - cold )
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
               DO i = 1, pdf_bin 
                  pdf_corr (i) = pdf_old (i) 
               ENDDO 
               CALL pdf_convert
            ENDIF 
         ENDIF 
!write(*,*) ' TESTED sca/den ACCEPT : ', laccept, pdf_scale, pdf_rho0,cnew, cold
!                                                                       
         IF (rmc_sigma.eq. - 9999.) laccept = .true. 
         IF (laccept) then 
            cold = cnew 
            CALL refine_adapt_pdf_sd (1)
         ELSE
            pdf_scale = pdf_old_scale
            pdf_rho0  = pdf_old_rho0 
            CALL refine_adapt_pdf_sd (0)
         ENDIF
      ENDIF
      END SUBROUTINE  pdf_rmc_scale
!*****7*****************************************************************
      SUBROUTINE  pdf_rmc_lattice( nmi, nma, cold, wtot, ee, psum, p2sum, &
                                    sig2, pmax, pn )

!
      USE crystal_mod
      USE pdf_mod
      USE rmc_mod
      USE refine_mod
USE lib_random_func
      USE random_mod
      USE spcgr_apply
      USE precision_mod
!
      IMPLICIT NONE
!
      REAL(PREC_DP), INTENT(INOUT) :: cold
      REAL(PREC_DP), INTENT(IN)    :: wtot
      REAL(PREC_DP), INTENT(IN)    :: ee
      REAL    , INTENT(INOUT) :: psum
      REAL    , INTENT(INOUT) :: p2sum
      REAL    , INTENT(INOUT) :: pmax
      REAL    , INTENT(INOUT) :: pn 
      REAL    , INTENT(IN)    :: sig2
!
      INTEGER  :: i, ip, nmi, nma, ipc
      LOGICAL  :: laccept
      REAL(PREC_DP) :: c, cc, ce, cnew
      REAL, DIMENSION(6) :: pdf_old_lattice
      REAL     :: sum
      REAL     :: prob
      REAL(PREC_DP), DIMENSION(PDF_MAXDAT) ::  pdf_old !  (MAXDAT) 
!
!
      pdf_old_lattice(1:3) = cr_a0     ! Backup old values
      pdf_old_lattice(4:6) = cr_win
laccept = .false.
      IF(pdf_refine_lattice(1)) THEN   ! Change lattice params, check crystal system
         cr_a0(1) = cr_a0(1) + gasdev(DBLE(ref_maxlatt(1)))
         laccept = .true.
         IF(cr_syst == cr_cubic) THEN
            cr_a0(2) = cr_a0(1)
            cr_a0(3) = cr_a0(1)
         ELSEIF(cr_syst == cr_hexagonal) THEN
            cr_a0(2) = cr_a0(1)
         ELSEIF(cr_syst == cr_trigonal ) THEN
            cr_a0(2) = cr_a0(1)
         ELSEIF(cr_syst == cr_tetragonal ) THEN
            cr_a0(2) = cr_a0(1)
         ENDIF
      ENDIF
      IF(pdf_refine_lattice(2)) THEN
         cr_a0(2) = cr_a0(2) + gasdev(DBLE(ref_maxlatt(1)))
         laccept = .true.
      ENDIF
      IF(pdf_refine_lattice(3)) THEN
         cr_a0(3) = cr_a0(3) + gasdev(DBLE(ref_maxlatt(1)))
         laccept = .true.
      ENDIF
      IF(pdf_refine_lattice(4)) THEN
         cr_win(1) = cr_win(1) + gasdev(DBLE(ref_maxlatt(2)))
         laccept = .true.
         IF(cr_syst == cr_rhombohed) THEN
            cr_win(2) = cr_win(1)
            cr_win(3) = cr_win(1)
         ENDIF
      ENDIF
      IF(pdf_refine_lattice(5)) THEN
         cr_win(2) = cr_win(2) + gasdev(DBLE(ref_maxlatt(2)))
         laccept = .true.
      ENDIF
      IF(pdf_refine_lattice(6)) THEN
         cr_win(3) = cr_win(3) + gasdev(DBLE(ref_maxlatt(2)))
         laccept = .true.
      ENDIF
!
!     Any lattice parameter to refine ?
!
      IF(laccept) THEN  ! LABEL laccept
!
!     Define new lattice parameters
!
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,       &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., cr_gmat, cr_fmat, &
         cr_cartesian,                                                     &
              cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
!
         DO i = 1, pdf_bin 
            pdf_old (i) = pdf_corr (i) 
         ENDDO 
         pdf_corr = 0.0  ! Clear pdf_corr, as we loop over all atoms
!
!        Calculate new PDF and chi2
!
         DO i = 1, cr_natoms 
            CALL pdf_addcorr (i, 1.0, sum) 
         ENDDO 
         CALL pdf_convert
         cnew = 0.0d0 
         c    = 0.0 
         cc   = 0.0 
         ce   = 0.0 
!                                                                       
         DO ip = nmi, nma/pdf_us_int
            ipc = (ip-1)*pdf_us_int + 1
            c   = c  + pdf_wic (ip) * pdf_calc (ipc) 
            cc  = cc + pdf_wic (ip) * pdf_calc (ipc) **2 
            ce  = ce + pdf_wic (ip) * pdf_calc (ipc) * pdf_obs (ip) 
         ENDDO 
         IF (rmc_doskal) then 
            pdf_skal = REAL(ce / cc )
         ELSE 
            pdf_skal = 1.0 / rmc_skal (1) 
         ENDIF 
         cnew = ee+pdf_skal**2 * cc - 2.0 * pdf_skal * ce 
         cnew = cnew / wtot 
!                                                                       
!     ----Accept move ?                                                 
!                                                                       
         prob = REAL( cnew - cold )
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
!write(*,*) 'Tested Lattice Acc ', laccept, cr_a0(1),pdf_old_lattice(1), psum,cold, cnew
         IF (rmc_sigma.eq. - 9999.) laccept = .true. 
         IF (laccept) then 
            cold = cnew 
            CALL refine_adapt_lattice (1)
         ELSE
!
!           Restore old lattice parameters
!
            cr_a0  = pdf_old_lattice(1:3)
            cr_win = pdf_old_lattice(4:6)
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,       &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., cr_gmat, cr_fmat, &
            cr_cartesian,                                                     &
            cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
            pdf_corr = 0.0  ! Clear pdf_corr, as we loop over all atoms
            DO i = 1, cr_natoms 
               CALL pdf_addcorr (i, 1.0, sum) 
            ENDDO 
            CALL pdf_convert
            CALL refine_adapt_lattice (0)
         ENDIF
      ENDIF ! LABEL laccept
      END SUBROUTINE  pdf_rmc_lattice
!*****7*****************************************************************
      SUBROUTINE pdf_makemove (natoms, i_new, p_new, isel, imol) 
!+                                                                      
!     make accepted move                                                
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: natoms 
      REAL    , INTENT(IN) :: p_new (3, rmc_max_atom) 
      INTEGER , INTENT(IN) :: i_new (rmc_max_atom) 
      INTEGER , INTENT(IN) :: isel (rmc_max_atom) 
      INTEGER , INTENT(IN) :: imol (rmc_max_atom) 
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
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE discus_plot_mod
      USE chem_mod
      USE crystal_mod 
      USE molecule_mod 
      USE pdf_mod 
      USE powder_pdf_hist_mod, ONLY : powder_trans_atoms_tocart, powder_trans_atoms_fromcart
      USE discus_plot_init_mod
!
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
      LOGICAL , INTENT(IN) :: lout 
!                                                                       
      INTEGER i, j, ia, id 
      REAL done, sum 
!
      INTEGER              :: npoint   !Number of points for histogram in exact mode
      INTEGER              :: nlook    !Number of look up dimensions 
      LOGICAL              :: all_atoms ! Are all atoms included in PDF?
      LOGICAL              :: do_mol   ! Take moleculoar B-values into account
      REAL, DIMENSION(1:3) :: u        ! Crystal diagonal
      REAL ss 
!                                                                       
      ss = seknds (0.0) 
      u  = 0.00
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
!------ Any molecules with b-value /= zero ?
!
      do_mol   = .false.
      pdf_nmol = 0
      nlook    = 0
      search_mol: DO i=1, mole_num_type
         IF(mole_biso(i) > 0.0) THEN
            do_mol   = .true.
            pdf_nmol = mole_num_type + mole_num_type*(mole_num_type+1)/2
            EXIT search_mol
         ENDIF
      ENDDO search_mol
!
!     Lay out look_up table for molecule entries
!
      IF(ALLOCATED(pdf_look_mol)) DEALLOCATE(pdf_look_mol)
      ALLOCATE(pdf_look_mol(0:mole_num_type,0:mole_num_type))
      IF(ALLOCATED(pdf_bvalue_mole)) DEALLOCATE(pdf_bvalue_mole)
      ALLOCATE(pdf_bvalue_mole(0:pdf_nmol))
      IF(ALLOCATED(pdf_clin_mole)) DEALLOCATE(pdf_clin_mole)
      ALLOCATE(pdf_clin_mole(0:pdf_nmol))
      IF(ALLOCATED(pdf_cqua_mole)) DEALLOCATE(pdf_cqua_mole)
      ALLOCATE(pdf_cqua_mole(0:pdf_nmol))
      pdf_look_mol    = 0
      pdf_bvalue_mole = 0.0
      pdf_clin_mole   = 0.0
      pdf_cqua_mole   = 0.0
      IF(pdf_nmol>0) THEN    ! Non-zero molecular bvalues
         nlook = 0
         DO i=1,mole_num_type
            IF(ABS(mole_biso(i))>0.0) THEN
               nlook = nlook + 1
               pdf_look_mol(0,i) = nlook
               pdf_look_mol(i,0) = nlook
               pdf_bvalue_mole(nlook) = mole_biso(i)
               pdf_clin_mole(nlook)   = mole_clin(i)
               pdf_cqua_mole(nlook)   = mole_cqua(i)
            ELSE
               pdf_look_mol(0,i) = 0
               pdf_look_mol(i,0) = 0
            ENDIF
         ENDDO
         nlook = mole_num_type
         DO i=1,mole_num_type
            DO j = i,mole_num_type
               IF(ABS(mole_biso(i))>0.0 .AND. ABS(mole_biso(j))>0.0) THEN
               nlook = nlook + 1
               pdf_look_mol(i,j) = nlook
               pdf_look_mol(j,i) = nlook
               pdf_bvalue_mole(nlook) = mole_biso(i) + mole_biso(j)
               pdf_clin_mole(nlook)   = mole_clin(i) + mole_clin(j)
               pdf_cqua_mole(nlook)   = mole_cqua(i) + mole_cqua(j)
               ELSEIF(ABS(mole_biso(i))==0.0 .AND. ABS(mole_biso(j))>0.0) THEN
                  pdf_look_mol(i,j) = pdf_look_mol(0,j)
                  pdf_look_mol(j,i) = pdf_look_mol(j,0)
               ELSEIF(ABS(mole_biso(i))>0.0 .AND. ABS(mole_biso(j))==0.0) THEN
                  pdf_look_mol(i,j) = pdf_look_mol(0,i)
                  pdf_look_mol(j,i) = pdf_look_mol(i,0)
               ELSE
                  pdf_look_mol(i,j) = 0
                  pdf_look_mol(j,i) = 0
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
      npoint = 1
      IF (pdf_lexact .AND. .NOT.chem_period(1)) then 
!
!------ Convert to cartesian
!
      CALL plot_ini_trans (1.0D0,                              &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
      CALL powder_trans_atoms_tocart (u)
      npoint    = INT(SQRT(u(1)**2+u(2)**2+u(3)**2)/pdf_deltar)
      npoint = MIN(npoint, pdf_ndat)
         npoint = pdf_ndat
!      IF(npoint > UBOUND(pdf_temp,1)) THEN
!        pdf_nscat = MAX(pdf_nscat, cr_nscat, PDF_MAXSCAT, MAXSCAT)
!        pdf_ndat  = MAX(pdf_ndat , npoint  , PDF_MAXDAT)
!        pdf_nbnd  = MAX(pdf_nbnd ,           PDF_MAXBND)
!        CALL alloc_pdf( pdf_nscat, pdf_ndat, pdf_nbnd )
!        IF ( ier_num < 0 ) THEN
!           RETURN
!        ENDIF
      ELSE
!
!------ Unit cell mode, size of pdf_temp is equal to pdf_ndat
!
         npoint = pdf_ndat
      ENDIF
!
!--------Reallocate pdf_temp only if dimensions have changed. 
!        An automatic change collides with RMC, which needs 
!        pdf_temp preserved accross its cycles!!!!!!!!!!!!!
         pdf_has_atom(:) = 0
         DO i=1, cr_natoms
            IF(pdf_allowed_i(cr_iscat(i)) .OR. pdf_allowed_j(cr_iscat(i)) ) THEN
               pdf_has_atom(cr_iscat(i)) = pdf_has_atom(cr_iscat(i)) + 1
            ENDIF
         ENDDO
         DO i=1, cr_nscat
            IF(pdf_has_atom(i)>0) pdf_nscat = i
         ENDDO
!      pdf_nscat = MAXVAL(cr_iscat(1:cr_natoms))
!     IF(npoint    > UBOUND(pdf_temp,1) .OR.       &
!        pdf_nscat > UBOUND(pdf_temp,2) .OR.       &
!        nlook     > UBOUND(pdf_temp,4)      ) THEN
!         pdf_ntemp  = INT(MAX(npoint,    PDF_MAXTEMP))
         pdf_ntemp = npoint
         IF(ALLOCATED(pdf_temp)) DEALLOCATE(pdf_temp)
         ALLOCATE(pdf_temp(0:pdf_ntemp,0:pdf_nscat,0:pdf_nscat,0:nlook))
!     ENDIF
!                                                                       
!------ Reset arrays                                                    
!                                                                       
!      DO i = 1, PDF_MAXDAT 
!      pdf_corr (i) = 0.0 
!      ENDDO 
!      DO i = 1, PDF_MAXDAT 
!      DO is = 1, cr_nscat 
!      DO js = 1, cr_nscat 
!      pdf_temp (i, is, js, 0) = 0
!      ENDDO 
!      ENDDO 
!      ENDDO 
      pdf_corr(:)  = 0.0
      pdf_temp(:,:,:,:)  = 0
!
      all_atoms = .true.
      do i=1,cr_nscat
         all_atoms = all_atoms .and. pdf_allowed_i(i) .and. pdf_allowed_j(i)
      enddo
!                                                                       
      sum = 0.0 
!                                                                       
!------ Start the calculation                                           
!                                                                       
      id = max (1, cr_natoms / 5) 
!                                                                       
      IF (pdf_lexact) THEN     ! Use exact loop over all atoms
         IF(all_atoms) THEN    ! All atom types are selected
            IF(do_mol) THEN
               CALL pdf_addcorr_e_all_mol(lout)
            ELSE
               IF(chem_period(1)) THEN
                  CALL pdf_addcorr_ep_all(lout)
               ELSE
                  CALL pdf_addcorr_e_all(lout)
               ENDIF
            ENDIF
         ELSE                  ! Partial PDF
            CALL pdf_addcorr_e (lout) 
         ENDIF 
      ELSE                     ! Use unit cell indexing for large periodic objects 
         IF(all_atoms) THEN    ! All atom types are selected
            loop1: DO ia = 1, cr_natoms 
               CALL pdf_addcorr_n_fast (ia) 
               IF(ier_ctrlc) THEN
                  ier_num = -14
                  ier_typ = ER_COMM
                  EXIT loop1
               ENDIF
               IF(ier_num/=0) EXIT loop1      ! An error occured or CTRL-C
               IF (lout.and. (mod (ia, id) .eq.0) ) THEN 
                  done = 100.0 * REAL(ia) / REAL(cr_natoms) 
                  WRITE (output_io, 1000) done 
               ENDIF 
            ENDDO loop1
        ELSE
            loop2: DO ia = 1, cr_natoms 
               CALL pdf_addcorr_n (ia) 
               IF(ier_ctrlc) THEN
                  ier_num = -14
                  ier_typ = ER_COMM
                  EXIT loop2
               ENDIF
               IF(ier_num/=0) EXIT loop2      ! An error occured or CTRL-C
               IF (lout.and. (mod (ia, id) .eq.0) ) THEN 
                  done = 100.0 * REAL(ia) / REAL(cr_natoms) 
                  WRITE (output_io, 1000) done 
               ENDIF 
            ENDDO loop2
         ENDIF 
      ENDIF 
!do is=1,cr_nscat
!   do js= 1, cr_nscat
!do ia=0,int(pdf_rmax/pdf_deltar)+1
!   write(10*is+js,*) ia*pdf_deltar,pdf_temp(ia,is,js,0)
!enddo
!enddo
!enddo
!
      IF (pdf_lexact .AND. .NOT.chem_period(1) .AND.(ier_num==0 .AND. .NOT.ier_ctrlc)) then 
!
!------ Convert back to crystal metric
!
         CALL powder_trans_atoms_fromcart
      ELSEIF (pdf_lexact .AND. .NOT.chem_period(1) .AND.(ier_num/=0 .OR. ier_ctrlc)) then 
!
!------ Convert back to crystal metric
!
         CALL powder_trans_atoms_fromcart
         RETURN     !ERRR, skip convtherm and convert
      ENDIF
!                                                                       
      CALL pdf_convtherm (1.0, sum) 
!                                                                       
!------ Convert to proper G(r)                                          
!                                                                       
      CALL pdf_convert 
!
         IF(ALLOCATED(pdf_temp)) DEALLOCATE(pdf_temp)
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
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE atom_env_mod 
      USE pdf_mod 
      USE wink_mod
!
      USE debug_mod 
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, k, ncc 
      INTEGER :: jpdf_bin
!     REAL ppp (MAXDAT) 
!      REAL, DIMENSION(PDF_MAXDAT   ) :: ppp ! (MAXDAT) 
!     REAL(PREC_DP), DIMENSION(:), ALLOCATABLE :: ppp ! (MAXDAT) 
      REAL norm, r, r0 
      REAL rr 
      REAL :: factor,fac4
      REAL :: c_sphere
!     REAL(PREC_DP) :: convlv
!     INTEGER (SELECTED_INT_KIND(9)) :: isign = 1
!                                                                       
      rr = 0.0
!     ALLOCATE(ppp(1:SIZE(pdf_calc)))
      ncc = cr_icc (1) * cr_icc (2) * cr_icc (3) 
      IF (.not.pdf_lrho0) then 
         IF (pdf_lrho0_rel) then 
            r0 = pdf_rho0 * cr_ncreal  / cr_v 
            r0 = pdf_rho0 * cr_n_real_atoms  / cr_v / ncc
         ELSE 
            r0 = pdf_rho0 
         ENDIF 
      ELSE 
         r0       = cr_n_real_atoms  / cr_v  /ncc
         pdf_rho0 = r0
      ENDIF 
      norm = 1.0 / REAL(cr_n_real_atoms) * pdf_scale 
      IF (.not.pdf_gauss.and.pdf_qalp.eq.0.0) then 
         norm = norm / pdf_deltar 
      ENDIF 
!                                                                       
      IF(r0>0.0) THEN
         c_sphere = cr_n_real_atoms/(4./3.*pi*(pdf_sphere/2)**3*r0)
      ELSE
         c_sphere = 1.0
      ENDIF
! c_sphere=1.0
      DO i = 1, pdf_bin 
      r = REAL(i) * pdf_deltar 
      IF (pdf_finite.eq.PDF_BACK_PERIOD) then 
         rr = 2.0 * REAL(zpi) * r * r0 * pdf_dnorm 
      ELSEIF (pdf_finite.eq.PDF_BACK_POLY) then 
         rr = 2.0 * REAL(zpi) * r * r0 * pdf_dnorm 
         IF (r.lt.pdf_diam_poly) then 
            DO k = 1, pdf_poly_n 
               rr = rr - pdf_poly (k) * r**k 
            ENDDO 
            rr = max (0.0, rr) 
         ELSE 
            rr = 0.0 
         ENDIF 
      ELSEIF (pdf_finite.eq.PDF_BACK_SPHERE) then 
         rr = 2.0 * REAL(zpi) * r * r0 * pdf_dnorm 
         IF (r.lt.pdf_sphere) then 
            rr = rr * (1. - 1.5 * (r / pdf_sphere) + .5 * (r /          &
            pdf_sphere) **3)*c_sphere
         ELSE 
            rr = 0.0 
         ENDIF 
      ELSEIF (pdf_finite.eq.PDF_BACK_TANH) then 
         rr = max (0.0, - 2.0 * REAL(zpi) * r * r0 * pdf_dnorm * tanh (       &
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
         factor = (pdf_deltar*pdf_sigmaq) * (pdf_deltar*pdf_sigmaq) /2.0
         fac4   = REAL(pdf_deltar/pdf_gauss_step*pdf_sigmaq )
         jpdf_bin = MIN(pdf_bin, IABS(INT(UBOUND(pdf_exp,1)/fac4-1)))
         DO i = 1, jpdf_bin 
!         r = REAL(i) * pdf_deltar 
!         pdf_calc (i) = pdf_calc (i) * exp ( - (r * pdf_sigmaq) **2 /   &
!         2.0)                                                           
            pdf_calc (i) = pdf_calc (i) * pdf_exp (NINT((i+1)*fac4))
         ENDDO 
      ENDIF 
!                                                                       
!------ Convolute with SINC function                                    
!                                                                       
!
      IF (pdf_qmax.gt.0.0) then 
!         DO i = 1, pdf_bin 
!         ppp (i) = pdf_calc (i) * (pdf_qmax - pdf_sinc (2 * i) ) 
!         DO k = 1, i - 1 
!         ppp (i) = ppp (i) + pdf_calc (k) * (pdf_sinc (i - k) - pdf_sinc (i + k) )
!         ENDDO 
!         DO k = i + 1, pdf_bin 
!         ppp (i) = ppp (i) + pdf_calc (k) * (pdf_sinc (k - i) - pdf_sinc (k + i) )
!         ENDDO 
!         ENDDO 
!                                                                       
         pdf_ppp = 0.0d0
         CALL CONVLV_SUB(SIZE(pdf_calc), SIZE(pdf_sincc),pdf_ppp,pdf_calc, pdf_sincc, 1)
         factor = pdf_deltar / REAL(zpi) * 2.
         DO i = 1, pdf_bin 
            pdf_calc (i) = pdf_ppp (i) * factor
!           pdf_calc (i) = pdf_ppp (i) * pdf_deltar / zpi * 2.0 
         ENDDO 
      ENDIF 
!     DEALLOCATE(ppp)
!                                                                       
      END SUBROUTINE pdf_convert                    
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr (ia, rsign, sum) 
!+                                                                      
!     Calculate correlation for given atom ia                           
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE atom_env_mod 
      USE celltoindex_mod
      USE molecule_mod
      USE pdf_mod 
      USE param_mod 
      USE wink_mod
      USE errlist_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN)  :: ia 
      REAL    , INTENT(IN)  :: rsign
      REAL    , INTENT(OUT) :: sum 
!                                                                       
      INTEGER ig, igaus, ib, ie 
      INTEGER :: jgaus  ! Limit checked igaus
      INTEGER i, j, k, ii, jj, kk, is, js, ks, iatom, ibin
      INTEGER istart (3), iend (3), iii (3), cell (3) 
!     REAL(PREC_DP) ppp (MAXDAT), gaus ( - MAXDAT:MAXDAT) 
!      REAL(PREC_DP), DIMENSION( PDF_MAXDAT)            :: ppp   !(MAXDAT)
      REAL(PREC_DP), DIMENSION(-PDF_MAXDAT:PDF_MAXDAT) :: gaus  ! ( - MAXDAT:MAXDAT) 
      REAL asym, gnorm, dist, dist2 !, rg 
      REAL sigma, fac , factor, fac4
      REAL dd (3), d (3), offset (3) 
!                                                                       
      fac = 1.0 / (2.0 * REAL(zpi)**2) 
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
         offset (ii) = REAL(iii (ii) - cell (ii) ) 
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
!write(*,*) ia, iatom,is,js, cr_mole(ia),cr_mole(iatom), mole_type(cr_mole(ia)), mole_type(cr_mole(iatom)),&
!           mole_biso(mole_type(cr_mole(ia))), mole_biso(mole_type(cr_mole(iatom)))
!                     IF(cr_mole(ia)/=cr_mole(iatom)) THEN
!                        sigma = sigma + fac*(mole_biso(mole_type(cr_mole(ia))) + &
!                                             mole_biso(mole_type(cr_mole(iatom))) )
!                     ENDIF
                     sigma = sigma - pdf_cquad_a / dist2 
                     sigma = sigma - pdf_clin_a / dist 
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
                  factor = REAL(pdf_deltar/pdf_gauss_step/sigma)
                  IF (sigma.le.0.0.or.igaus.lt.2.or. factor>UBOUND(pdf_exp,1)/2) then 
                     pdf_corr (ibin) = pdf_corr (ibin) + rsign *        &
                     pdf_weight (is, js) / pdf_deltar                   
                  ELSE 
                     factor = REAL(pdf_deltar/pdf_gauss_step/sigma)
                     fac4   = pdf_deltar/dist
                     gnorm = 1.0 / (sqrt (REAL(zpi)) * sigma) 
!                                                                       
                     jgaus = MIN(igaus, IABS(INT(UBOUND(pdf_exp,1)/factor+1)))
                     DO ig = - jgaus, jgaus 
!                    rg = (ig - 1) * pdf_deltar 
!                    asym = 1.0 + rg / dist 
                     asym = 1.0 + (ig-0)*fac4  ! Corrected shift in Gaussian
!                    gaus (ig) = gnorm * asym * exp ( - 0.5 * (rg /     &
!                    sigma) **2)                                        
                     gaus (ig) = gnorm * asym * pdf_exp(IABS(INT((IABS(ig)-0)*factor)))
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
      ENDIF 
!
!                                                                       
      END SUBROUTINE pdf_addcorr                    
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_n_fast (ia) 
!+                                                                      
!     Calculate correlation for given atom ia, fast version             
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE celltoindex_mod
!     USE modify_mod
      USE molecule_mod
      USE pdf_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: ia
!                                                                       
      INTEGER i, j, k, ii, is, js, ks, iatom, ibin , islook
      INTEGER  :: ipdf_rmax
      INTEGER istart (3), iend (3), iii (3), cell (3) 
      INTEGER  :: offzero
      REAL dist
      REAL dd (3), offset (3)
!
      ipdf_rmax = int(pdf_rmax/pdf_deltar)+1
!
      is = cr_iscat (ia) 
      IF (pdf_allowed_i (is) .or.pdf_allowed_j (is) ) then 
         pdf_temp(0,is,is,0) = pdf_temp(0,is,is,0) + 1       ! Add self correlation peak
         CALL indextocell (ia, iii, ks) 
         DO i = 1, 3 
         istart (i) = iii (i) - 1 - int (pdf_rmax / cr_a0 (i) ) 
         iend (i)   = iii (i) + 1 + int (pdf_rmax / cr_a0 (i) ) 
         IF (pdf_2d.and.cr_icc (i) .le.1) then 
            istart (i) = iii (i) 
            iend (i)   = iii (i) 
         ENDIF 
         ENDDO 
!                                                                       
!------ - In case we do not want to loop around cut range               
!                                                                       
         IF ( (.not.chem_period (1) ) .and. (.not.chem_period (2) )     &
                                      .and. (.not.chem_period (3) ) ) then                           
            DO i = 1, 3 
            istart (i) = max (1, istart (i) ) 
            iend (i)   = min (cr_icc (i), iend (i) ) 
            ENDDO 
         ENDIF 
!                                                                       
!------ - Here starts the inner loop                                    
!                                                                       
         DO k = istart (3), iend (3) 
         DO j = istart (2), iend (2) 
         DO i = istart (1), iend (1) 
         cell (3) = k 
         cell (2) = j 
         cell (1) = i 
!                                                                       
!---------- Here we apply periodic boundaries                           
!---------- Modified to allow more than a single cycle !                
!                                                                       
         offzero = 0
         DO ii = 1, 3 
            iii (ii)    = cell (ii) 
            cell (ii)   = pdf_bnd (ii, cell (ii) ) 
            offset (ii) = REAL(iii (ii) - cell (ii) ) 
            offzero     = MAX(offzero,NINT(ABS(offset(ii))))
         ENDDO 
!                                                                       
!------ --- Now look for neighbours in surrounding unit cells only      
!                                                                       
         DO ii = 1, cr_ncatoms 
         CALL celltoindex (cell, ii, iatom) 
         js = cr_iscat (iatom) 
         IF ( (pdf_allowed_i (js) .and.pdf_allowed_j (js) ) )THEN
               IF(offzero==0 .AND. cr_mole(ia)==cr_mole(iatom)) THEN
                  islook = 0   ! Atoms are within the same molecule
               ELSE
                  islook = pdf_look_mol(mole_type(cr_mole(ia)),mole_type(cr_mole(iatom)))
               ENDIF
!           DO jj = 1, 3 
!              dd (jj) = cr_pos (jj, ia) - cr_pos (jj, iatom) - offset (jj) 
!           ENDDO 
               dd ( 1) = cr_pos ( 1, ia) - cr_pos ( 1, iatom) - offset ( 1) 
               dd ( 2) = cr_pos ( 2, ia) - cr_pos ( 2, iatom) - offset ( 2) 
               dd ( 3) = cr_pos ( 3, ia) - cr_pos ( 3, iatom) - offset ( 3) 
!           ibin = int((SQRT(                                             &
!                   dd(1)*dd(1)*cr_gten(1,1) + dd(2)*dd(2)*cr_gten(2,2) + &
!                   dd(3)*dd(3)*cr_gten(3,3) + 2. * (                     &
!                   dd(1)*dd(2)*cr_gten(1,2) + dd(1)*dd(3)*cr_gten(1,3) + &
!                   dd(2)*dd(3)*cr_gten(2,3)))+pdf_deltars)/pdf_deltar)
!           ibin = int((SQRT(                                             &
            dist  = SQRT(dd(1)*dd(1)*cr_gten(1,1) + dd(2)*dd(2)*cr_gten(2,2) + &
                    dd(3)*dd(3)*cr_gten(3,3) + 2. * (                     &
                    dd(1)*dd(2)*cr_gten(1,2) + dd(1)*dd(3)*cr_gten(1,3) + &
                    dd(2)*dd(3)*cr_gten(2,3)))
!            dist = sqrt (dist2) 
!                                                                       
             IF (dist.le.pdf_rmax                       ) then 
!            IF (ibin.le.ipdf_rmax) then 
!              ibin = nint (dist / pdf_deltar) 
               ibin =  int((dist+pdf_deltars) / pdf_deltar) 
               pdf_temp (ibin, is, js,islook) = pdf_temp (ibin, is, js,islook) + 1
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE pdf_addcorr_n_fast                  
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_n (ia) 
!+                                                                      
!     Calculate correlation for given atom ia, fast version             
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE celltoindex_mod
!     USE modify_mod
      USE molecule_mod
      USE pdf_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: ia
!                                                                       
      INTEGER i, j, k, ii, jj, is, js, ks, iatom, ibin 
      INTEGER istart (3), iend (3), iii (3), cell (3) 
      INTEGER  :: offzero, islook
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
         offzero = 0
         DO ii = 1, 3 
         iii (ii) = cell (ii) 
         cell (ii) = pdf_bnd (ii, cell (ii) ) 
         offset (ii) = REAL(iii (ii) - cell (ii) ) 
            offzero     = MAX(offzero,NINT(ABS(offset(ii))))
         ENDDO 
!                                                                       
!------ --- Now look for neighbours in surrounding unit cells only      
!                                                                       
         DO ii = 1, cr_ncatoms 
         CALL celltoindex (cell, ii, iatom) 
         js = cr_iscat (iatom) 
         IF ( (pdf_allowed_i (js) .and.pdf_allowed_j (js) ) )THEN
               IF(offzero==0 .AND. cr_mole(ia)==cr_mole(iatom)) THEN
                  islook = 0   ! Atoms are within the same molecule
               ELSE
                  islook = pdf_look_mol(mole_type(cr_mole(ia)),mole_type(cr_mole(iatom)))
               ENDIF
            DO jj = 1, 3 
            dd (jj) = cr_pos (jj, ia) - cr_pos (jj, iatom) - offset (jj) 
            d (jj) = abs (dd (jj) ) * cr_a0 (jj) 
            ENDDO 
            dist2 = dd (1) * dd (1) * cr_gten (1, 1) + dd (2) * dd (2)  &
            * cr_gten (2, 2) + dd (3) * dd (3) * cr_gten (3, 3) + 2. *  &
            dd (1) * dd (2) * cr_gten (1, 2) + 2. * dd (1) * dd (3)     &
            * cr_gten (1, 3) + 2. * dd (2) * dd (3) * cr_gten (2, 3)    
            dist = sqrt (dist2) 
!                                                                       
            IF (dist.le.pdf_rmax                       ) then 
               ibin = nint (dist / pdf_deltar) 
               pdf_temp (ibin, is, js,islook) = pdf_temp (ibin, is, js,islook) + 1
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
      SUBROUTINE pdf_addcorr_e (lout) 
!+                                                                      
!     Calculate correlation for given atom ia, exact version            
!     Version for partial PDF, i.e. check on pdf-allowed is required
!     As this is not super fast, the check on molecules is included here
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE molecule_mod 
      USE pdf_mod 
      USE errlist_mod 
      USE prompt_mod
      IMPLICIT none 
!
      LOGICAL, INTENT(IN) :: lout
!
!
      INTEGER   :: id
      INTEGER   :: is, js, ia, iatom, ibin , islook
      INTEGER   :: ipdf_rmax
      REAL      :: done
      REAL dd (3)
!                                                                       
      id = MAX(100, cr_natoms/5)    ! Progress report 20% or every 100 atoms
      ipdf_rmax = int(pdf_rmax/pdf_deltar)+1
main: DO ia=1,cr_natoms    ! Outer loop over all atoms
         IF (lout.and. (mod (ia, id) .eq.0) ) THEN 
            done = 100.0 * REAL(ia) / REAL(cr_natoms) 
            WRITE (output_io, 1000) done 
         ENDIF 
         is = cr_iscat (ia) 
         IF (pdf_allowed_i (is) .or.pdf_allowed_j (is) ) THEN 
            pdf_temp (0, is, is,0) = pdf_temp (0, is, is,0) +1  ! Element is,is is 
!                                                           ! excluded in the inner loop
!                                                                       
!------ - Here starts the inner loop over all atoms                     
!                                                                       
inner:      DO iatom = ia+1, cr_natoms 
               js = cr_iscat (iatom) 
               IF(cr_mole(ia)==cr_mole(iatom)) THEN
                  islook = 0   ! Atoms are within the same molecule
               ELSE
                  islook = pdf_look_mol(mole_type(cr_mole(ia)),mole_type(cr_mole(iatom)))
               ENDIF
               IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or.  &
                    (pdf_allowed_j (is) .and.pdf_allowed_i (js) ) ) THEN
                  dd (1) = cr_pos (1, ia) - cr_pos (1, iatom) 
                  dd (2) = cr_pos (2, ia) - cr_pos (2, iatom) 
                  dd (3) = cr_pos (3, ia) - cr_pos (3, iatom) 
                  ibin=int((SQRT(dd (1) * dd (1) + &
                                 dd (2) * dd (2) + &
                                 dd (3) * dd (3)  )+pdf_deltars)/pdf_deltar)
!                                                                       
!                  IF (dist.le.pdf_rmax) THEN
                  IF (ibin.le.ipdf_rmax) THEN
                     pdf_temp (ibin, is, js,islook) = pdf_temp (ibin, is, js,islook) +1
                     pdf_temp (ibin, js, is,islook) = pdf_temp (ibin, js, is,islook) +1
                  ENDIF 
               ENDIF 
            ENDDO inner
         ENDIF 
      ENDDO  main
!                                                                       
1000  FORMAT     (  '   ',f6.2,' % done  ...') 
!
      END SUBROUTINE pdf_addcorr_e                  
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_e_all(lout)
!+                                                                      
!     Calculate correlation for given atom ia, exact version, all atoms
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      USE errlist_mod 
      USE prompt_mod
USE support_mod
      IMPLICIT none 
!
      LOGICAL, INTENT(IN) :: lout
!                                                                       
      INTEGER :: id
      INTEGER :: is, js, ia, iatom, ibin 
      INTEGER   :: ipdf_rmax
      REAL                  :: done
      REAL   , DIMENSION(3) :: dd
      REAL :: ss
!
      id = MAX(100, cr_natoms/5)    ! Progress report 20% or every 100 atoms
      ipdf_rmax = int(pdf_rmax/pdf_deltar)+1
      ss = seknds (0.0)
main: DO ia=1,cr_natoms    ! Outer loop over all atoms
         IF (lout.and. (mod (ia, id) .eq.0) ) then 
            done = 100.0 * REAL(ia) / REAL(cr_natoms) 
            WRITE (output_io, 1000) done 
         ENDIF 
         is = cr_iscat (ia) 
         IF (is /= 0 ) THEN   ! disregard VOIDs
            pdf_temp (0, is, is,0) = pdf_temp (0, is, is,0) +1  ! Element is,is is 
!                                                           ! excluded in the inner loop
!                                                                       
!------ - Here starts the inner loop over all atoms                     
!                                                                       
inner:      DO iatom = ia+1, cr_natoms 
               js = cr_iscat (iatom) 
               IF ( js /= 0) THEN 
                  dd (1) = cr_pos (1, ia) - cr_pos (1, iatom) 
                  dd (2) = cr_pos (2, ia) - cr_pos (2, iatom) 
                  dd (3) = cr_pos (3, ia) - cr_pos (3, iatom) 
                  ibin=int((SQRT(dd (1) * dd (1) + &
                                 dd (2) * dd (2) + &
                                 dd (3) * dd (3)  )+pdf_deltars)/pdf_deltar)
!                                                                       
!                  IF (dist.le.pdf_rmax) THEN 
                  IF (ibin.le.ipdf_rmax) THEN
                     pdf_temp (ibin, is, js,0) = pdf_temp (ibin, is, js,0) +1
                     pdf_temp (ibin, js, is,0) = pdf_temp (ibin, js, is,0) +1
                  ENDIF 
               ENDIF 
            ENDDO  inner
         ENDIF 
      ENDDO  main
      ss = seknds (ss )
      WRITE (output_io, 4000) ss
4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec')
!                                                                       
1000  FORMAT     (  '   ',f6.2,' % done  ...') 
!
      END SUBROUTINE pdf_addcorr_e_all
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_ep_all(lout)
!+                                                                      
!     Calculate correlation for given atom ia, exact version, all atoms
!     Periodic version
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      USE errlist_mod 
      USE prompt_mod
USE support_mod
      IMPLICIT none 
!
      LOGICAL, INTENT(IN) :: lout
!                                                                       
      INTEGER :: id
      INTEGER :: is, js, ia, iatom, ibin 
      INTEGER   :: ipdf_rmax
      INTEGER, DIMENSION(3) :: ncell
      INTEGER, DIMENSION(3) :: icell
      INTEGER               :: ix,iy,iz   !Loop over periodic boundaries
      REAL                  :: done, dist
      REAL   , DIMENSION(3) :: dd
      REAL :: ss
!
      icell(1) = MAX(cr_icc(1),NINT(cr_dim0(1,2)-cr_dim0(1,1)))
      icell(2) = MAX(cr_icc(1),NINT(cr_dim0(2,2)-cr_dim0(2,1)))
      icell(3) = MAX(cr_icc(1),NINT(cr_dim0(3,2)-cr_dim0(3,1)))
      ncell(1) = INT(pdf_rmax/cr_a0(1)/icell(1)) + 1
      ncell(2) = INT(pdf_rmax/cr_a0(2)/icell(2)) + 1
      ncell(3) = INT(pdf_rmax/cr_a0(3)/icell(3)) + 1
!
      id = MAX(100, cr_natoms/5)    ! Progress report 20% or every 100 atoms
      ipdf_rmax = int(pdf_rmax/pdf_deltar)+1
      ss = seknds (0.0)
main: DO ia=1,cr_natoms    ! Outer loop over all atoms
         IF (lout.and. (mod (ia, id) .eq.0) ) then 
            done = 100.0 * REAL(ia) / REAL(cr_natoms) 
            WRITE (output_io, 1000) done 
         ENDIF 
         is = cr_iscat (ia) 
         IF (is /= 0 ) THEN   ! disregard VOIDs
!           pdf_temp (0, is, is,0) = pdf_temp (0, is, is,0) +1  ! Element (is,is) is 
!                                                           ! excluded in the inner loop
!                                                                       
!------ - Here starts the inner loop over all atoms                     
!                                                                       
!inner:      DO iatom = ia+1, cr_natoms 
inner:      DO iatom =    1, cr_natoms 
               js = cr_iscat (iatom) 
               IF ( js /= 0) THEN 
                  DO ix=-ncell(1),ncell(1)
                  DO iy=-ncell(2),ncell(2)
                  DO iz=-ncell(3),ncell(3)
                  dd (1) = cr_pos (1, ia) - cr_pos (1, iatom) + REAL(ix*icell(1))
                  dd (2) = cr_pos (2, ia) - cr_pos (2, iatom) + REAL(iy*icell(2))
                  dd (3) = cr_pos (3, ia) - cr_pos (3, iatom) + REAL(iz*icell(3))
            dist  = SQRT(dd(1)*dd(1)*cr_gten(1,1) + dd(2)*dd(2)*cr_gten(2,2) + &
                    dd(3)*dd(3)*cr_gten(3,3) + 2. * (                     &
                    dd(1)*dd(2)*cr_gten(1,2) + dd(1)*dd(3)*cr_gten(1,3) + &
                    dd(2)*dd(3)*cr_gten(2,3)))
            ibin=INT((dist+pdf_deltars)/pdf_deltar)
!                 ibin=int((SQRT(dd (1) * dd (1) + &
!                                dd (2) * dd (2) + &
!                                dd (3) * dd (3)  )+pdf_deltars)/pdf_deltar)
!                                                                       
!                  IF (dist.le.pdf_rmax) THEN 
                  IF (ibin.le.ipdf_rmax) THEN
                     pdf_temp (ibin, is, js,0) = pdf_temp (ibin, is, js,0) +1
!                    pdf_temp (ibin, js, is,0) = pdf_temp (ibin, js, is,0) +1
                  ENDIF 
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF 
            ENDDO  inner
         ENDIF 
      ENDDO  main
      ss = seknds (ss )
      WRITE (output_io, 4000) ss
4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec')
!                                                                       
1000  FORMAT     (  '   ',f6.2,' % done  ...') 
!
      END SUBROUTINE pdf_addcorr_ep_all
!*****7*****************************************************************
      SUBROUTINE pdf_addcorr_e_all_mol(lout)
!+                                                                      
!     Calculate correlation for given atom ia, exact version, all atoms
!     Version for molecular b-value
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      USE molecule_mod 
      USE errlist_mod 
      USE prompt_mod
USE support_mod
      IMPLICIT none 
!
      LOGICAL, INTENT(IN) :: lout
!                                                                       
      INTEGER :: id
      INTEGER :: is, js, ia, iatom, ibin , islook
      INTEGER   :: ipdf_rmax
      REAL                  :: done
      REAL   , DIMENSION(3) :: dd
      REAL :: ss
!
      id = MAX(100, cr_natoms/5)    ! Progress report 20% or every 100 atoms
      ipdf_rmax = int(pdf_rmax/pdf_deltar)+1
      ss = seknds (0.0)
main: DO ia=1,cr_natoms    ! Outer loop over all atoms
         IF (lout.and. (mod (ia, id) .eq.0) ) then 
            done = 100.0 * REAL(ia) / REAL(cr_natoms) 
            WRITE (output_io, 1000) done 
         ENDIF 
         is = cr_iscat (ia) 
         IF (is /= 0 ) THEN   ! disregard VOIDs
            pdf_temp (0, is, is,0) = pdf_temp (0, is, is,0) +1  ! Element is,is is 
!                                                           ! excluded in the inner loop
!                                                                       
!------ - Here starts the inner loop over all atoms                     
!                                                                       
inner:      DO iatom = ia+1, cr_natoms 
               js = cr_iscat (iatom) 
               IF(cr_mole(ia)==cr_mole(iatom)) THEN
                  islook = 0   ! Atoms are within the same molecule
               ELSE
                  islook = pdf_look_mol(mole_type(cr_mole(ia)),mole_type(cr_mole(iatom)))
               ENDIF
               IF ( js /= 0) THEN 
                  dd (1) = cr_pos (1, ia) - cr_pos (1, iatom) 
                  dd (2) = cr_pos (2, ia) - cr_pos (2, iatom) 
                  dd (3) = cr_pos (3, ia) - cr_pos (3, iatom) 
                  ibin=int((SQRT(dd (1) * dd (1) + &
                                 dd (2) * dd (2) + &
                                 dd (3) * dd (3)  )+pdf_deltars)/pdf_deltar)
!                                                                       
!                  IF (dist.le.pdf_rmax) THEN 
                  IF (ibin.le.ipdf_rmax) THEN
                     pdf_temp (ibin, is, js,islook) = pdf_temp (ibin, is, js,islook) +1
                     pdf_temp (ibin, js, is,islook) = pdf_temp (ibin, js, is,islook) +1
                  ENDIF 
               ENDIF 
            ENDDO  inner
         ENDIF 
      ENDDO  main
      ss = seknds (ss )
      WRITE (output_io, 4000) ss
4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec')
!                                                                       
1000  FORMAT     (  '   ',f6.2,' % done  ...') 
!
      END SUBROUTINE pdf_addcorr_e_all_mol
!*****7*****************************************************************
      SUBROUTINE pdf_convtherm (rsign, sum) 
!+                                                                      
!     Convolute the pair correlation histograms with the                
!     thermal Gaussian                                                  
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE pdf_mod 
      USE errlist_mod 
      USE wink_mod
      USE precision_mod
      IMPLICIT none 
!                                                                       
      REAL , INTENT(IN)  :: rsign
      REAL , INTENT(OUT) :: sum 
!                                                                       
      INTEGER ig, igaus, ib, ie , jgaus
      INTEGER ii, is, js, ibin , ibin1
      INTEGER :: il   ! Index for mole B-values
!     REAL(PREC_DP) gaus ( - MAXDAT:MAXDAT) 
      REAL(PREC_DP), DIMENSION(- PDF_MAXDAT:PDF_MAXDAT) :: gaus ! ( - MAXDAT:MAXDAT) 
      REAL asym, gnorm, dist, dist2 !, rg 
      REAL sigma, fac , factor, fac4
      REAL :: sqrt_zpi
!                                                                       
!   open(45,file='POWDER/hist.pdf_init',status='unknown')
!   do ii=1,UBOUND(pdf_temp,1)
!   write(45, '(i7,4(1x,F18.6))') ii,REAL(pdf_temp(ii,1,1,0))
!   enddo
!   close (45)

      fac = 1.0 / (2.0 * REAL(zpi)**2) 
      sqrt_zpi =1.0/sqrt(REAL(zpi))
      loop_is: DO is = 1, cr_nscat 
      IF(pdf_has_atom(is)>0) THEN
      loop_js: DO js = 1, cr_nscat 
      IF(pdf_has_atom(js)>0) THEN
      IF ( (pdf_allowed_i (is) .and.pdf_allowed_j (js) ) .or. &
           (pdf_allowed_j (is) .and.pdf_allowed_i (js) ) ) then               
         mole_loop: DO il=0,UBOUND(pdf_temp,4)
         ii = int (pdf_rmax / pdf_deltar) + 1 
         DO ibin = 1, ii 
         zero: IF(pdf_temp(ibin,is,js,il)>0) THEN
         dist = ibin * pdf_deltar 
         dist2 = dist**2 
!                                                                       
!------ --------- Convolute with Gaussian                               
!                                                                       
         IF (pdf_gauss.or.pdf_qalp.gt.0) then 
            sigma = 0.0 
            IF (pdf_gauss) then 
               sigma = fac * (cr_dw (is) + cr_dw (js) + pdf_bvalue_mole(il)) 
               sigma = sigma - pdf_cquad_a / dist2  - pdf_cqua_mole(il)/dist2
               sigma = sigma - pdf_clin_a / dist   - pdf_clin_mole(il)/dist
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
               js,il) * rsign * pdf_weight (is, js) / pdf_deltar           
            ELSE 
!              gnorm = 1.0 / (sqrt (zpi) * sigma) 
               gnorm =        sqrt_zpi / sigma 
!                                                                       
               factor = REAL(pdf_deltar/pdf_gauss_step/(sigma))
               fac4   = REAL(pdf_deltar/dist)
               jgaus = MIN(igaus, IABS(INT(UBOUND(pdf_exp,1)/factor+0)), &
                                  ABS(UBOUND(gaus,1)))
!if(jgaus>UBOUND(gaus,1) ) then
!   write(*,*) ' ELEMENTS ', is, js
!   write(*,*) ' sigma    ', sigma
!   write(*,*) ' igaus    ', igaus
!   write(*,*) ' jgaus    ', jgaus
!   write(*,*) ' factor   ', factor
!   write(*,*) ' PDF_EXP  ', UBOUND(pdf_exp,1)
!endif
               DO ig = - jgaus, jgaus 
!                 rg = (ig - 1) * pdf_deltar 
!                 asym = 1.0 + rg / dist 
!                 gaus (ig) = gnorm * asym * exp (-0.5*(rg/sigma)** 2)
                  asym = 1.0 + (ig + 0)*fac4   ! Corrected Shift in Gauss and in next line
                  gaus (ig) = gnorm * asym * pdf_exp(IABS(INT((IABS(ig)-0)*factor)))
               ENDDO 
!                                                                       
               fac4 = pdf_temp (ibin, is, js,il) * rsign * pdf_weight (is, js)
               ibin1 = ibin !-1   ! Corrected Shift in Gauss
               DO ig = ib, ie 
!                 kk = ig - ibin + 1 
!                 pdf_corr (ig) = pdf_corr (ig) + pdf_temp (ibin, is, js,il)  &
!                 * rsign * pdf_weight (is, js) * gaus (kk)                
                  pdf_corr (ig) = pdf_corr (ig) + fac4 *gaus(ig - ibin1  )
               ENDDO 
            ENDIF 
!                                                                       
!------ --------- Just take the distance                                
!                                                                       
         ELSE 
            pdf_corr (ibin) = pdf_corr (ibin) + pdf_temp (ibin, is, js,il) &
            * rsign * pdf_weight (is, js)                               
         ENDIF 
         sum = sum + pdf_weight (is, js) 
            ENDIF zero
         ENDDO 
         ENDDO  mole_loop
      ENDIF 
      ENDIF 
      ENDDO loop_js
      ENDIF 
      ENDDO loop_is
!                                                                       
!   open(45,file='POWDER/hist.pdf_conv',status='unknown')
!   do ii=1,UBOUND(pdf_corr,1)
!   write(45, '(2F18.8)') ii*pdf_deltar,pdf_corr(ii)
!   write(45, '(i7,4(1x,F18.6))') ii,pdf_temp(ii, 1, 1, 0)
!   enddo
!   close (45)
      END SUBROUTINE pdf_convtherm                  
!
!*****7*****************************************************************
!
SUBROUTINE pdf_reset
!
USE discus_allocate_appl_mod
USE pdf_mod
USE precision_mod
!
IMPLICIT NONE
!
CALL alloc_pdf(1, 1,  1, 1)
PDF_MAXSCAT      = 1
PDF_MAXDAT       = 1
PDF_MAXBND       = 1
PDF_MAXTEMP      = 1
PDF_MAXSINCC     = 2**12+1
!
pdf_nscat = 1
pdf_ndat  = 1
pdf_nbnd  = 1
pdf_ntemp = 1
!
IF(ALLOCATED(pdf_calc))   pdf_calc(:)       = 0.0D0  ! (MAXDAT)
IF(ALLOCATED(pdf_ppp))    pdf_ppp (:)       = 0.0D0  ! (MAXDAT)
IF(ALLOCATED(pdf_corr))   pdf_corr(:)       = 0.0D0  ! (MAXDAT)
IF(ALLOCATED(pdf_temp))   pdf_temp(:,:,:,:) = 0      ! (MAXTEMP,0:MAXSCAT,0:MAXSCAT)
IF(ALLOCATED(pdf_obs))    pdf_obs(:)        = 0.0    ! (MAXDAT)
IF(ALLOCATED(pdf_wic))    pdf_wic(:)        = 0.0    ! (MAXDAT)
!
IF(ALLOCATED(pdf_sincc))  pdf_sincc(:)      = 0.0D0  ! (2*MAXDAT)
IF(ALLOCATED(pdf_weight)) pdf_weight(:,:)   = 0.0    ! (0:PDF_MAXSCAT,0:PDF_MAXSCAT)
IF(ALLOCATED(pdf_exp))    pdf_exp(:)        = 0.0D0  ! (4000)
!
pdf_qmax   = 30.00
pdf_deltari=  0.001   ! internal delta r
pdf_deltar =  0.001   ! internal delta r
pdf_deltars=  0.0005
pdf_deltaru=  0.01    ! User supplied delta r
pdf_rmin   = pdf_deltari ! Minimum distance to calculate, internal value
pdf_rminu  =  0.01    ! Minimum distance to calculate, user value
pdf_rmax   = 50.00    ! Maximum distance to calculate, internal value
pdf_rmaxu  = 50.00    ! Maximum distance to calculate, user value
pdf_us_int =  1       ! Ratio user steps to internal steps
pdf_mode   =  PDF_DO_CALC ! PDF mode, default to 'calc'
pdf_skal   =  1.00
pdf_sigmaq =  0.00
pdf_xq     =  0.00
pdf_rfmin  =  0.01    ! distance range for refinement, internal value
pdf_rfmax  =  0.01    ! distance range for refinement, internal value
pdf_rfminu =  0.01    ! distance range for refinement, user value
pdf_rfmaxu =  0.01    ! distance range for refinement, user value
pdf_rfminf =  0.01    ! distance range for refinement, file value
pdf_rfmaxf =  0.01    ! distance range for refinement, file value
pdf_cquad_a=  0.00
pdf_rcut   =  0.00
pdf_srat   =  1.00
pdf_clin_a  =  0.00
pdf_qalp   =  0.00
pdf_dnorm  =  1.00
pdf_rho0   =  0.00
pdf_sphere =  0.00
pdf_diam_poly =  0.00
pdf_diam   =  0.00
pdf_shape  =  0.00
pdf_scale  =  1.00
pdf_poly(5)=  0.00
pdf_success= .FALSE.
!
IF(ALLOCATED(pdf_bnd))  pdf_bnd(:,:) = 0     ! (3,-MAXBND:2*MAXBND)
!
pdf_bin    = 1
pdf_finite = PDF_BACK_PERIOD
pdf_poly_n = 0
pdf_sel_prop(0:1) = 0
!                                                
pdf_radiation = PDF_RAD_XRAY
pdf_power     = 4
pdf_nmol      = 0   ! pdf_temp dimension if molecules are relevant
pdf_lxray       = .false.
pdf_gauss       = .false.
pdf_gauss_init  = .true.
pdf_2d          = .false.
IF(ALLOCATED(pdf_allowed_i))   pdf_allowed_i(:)   = .TRUE. ! (0:PDF_MAXSCAT)
IF(ALLOCATED(pdf_allowed_j))   pdf_allowed_j(:)   = .TRUE. ! (0:PDF_MAXSCAT)
IF(ALLOCATED(pdf_lsite_i))     pdf_lsite_i  (:)   = .TRUE. ! (0:PDF_MAXSCAT)
IF(ALLOCATED(pdf_lsite_j))     pdf_lsite_j  (:)   = .TRUE. ! (0:PDF_MAXSCAT)
IF(ALLOCATED(pdf_look_mol))    pdf_look_mol(:,:)  = 0      ! (0:PDF_MAXSCAT)
IF(ALLOCATED(pdf_bvalue_mole)) pdf_bvalue_mole(:) = 0.0    ! effective mol bvalues
IF(ALLOCATED(pdf_clin_mole))   pdf_clin_mole(:)   = 0.0    ! linear correction mol
IF(ALLOCATED(pdf_cqua_mole))   pdf_cqua_mole(:)   = 0.0    ! quadratic correction mol
!
pdf_ldata     = .false.
pdf_lweights  = .false.
pdf_lrho0     = .true.
pdf_lexact    = .false.
pdf_lrho0_rel = .false.
pdf_size_of   = 0
!
pdf_refine_scale   = .false.
pdf_refine_density = .false.
pdf_refine_lattice(:) = .false.
!
END SUBROUTINE pdf_reset
!
!*****7*****************************************************************
!
END MODULE pdf_menu
