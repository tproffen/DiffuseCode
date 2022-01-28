MODULE powder
!
USE errlist_mod 
!
IMPLICIT NONE
!
PUBLIC
!
CONTAINS
!+                                                                      
!     Calculation of powder diffraction pattern.                        
!                                                                       
!*****7*****************************************************************
SUBROUTINE do_powder 
!-                                                                      
!     Main menu for powder diffraction pattern                          
!+                                                                      
USE discus_config_mod 
USE diffuse_mod 
USE crystal_mod 
USE diffuse_mod 
USE fourier_sup
USE phases_mod
USE phases_set_mod
USE powder_mod 
USE powder_write_mod 
USE powder_pdf_hist_mod
USE pdf_mod
USE discus_show_menu
!
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
USE get_params_mod
USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE class_macro_internal
USE precision_mod
USE prompt_mod 
USE str_comp_mod
USE sup_mod
use take_param_mod
USE wink_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MIN_PARA= 5
CHARACTER(LEN=5) :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=PREC_STRING) :: line, zeile
CHARACTER (LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: cpara ! (MIN(10,MAXSCAT)) 
INTEGER                    , DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: lpara
REAL(KIND=PREC_DP)         , DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: werte
integer :: MAXW
INTEGER :: ianz
INTEGER :: lp, length, lbef 
INTEGER :: indxg
LOGICAL :: lend
real(kind=PREC_DP) :: fwhm  ! FWHM at Qmax
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_TABLE   = 1
character(LEN=   5), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'table'  /
data loname /  2       /
opara  =  (/ 'waas  ' /)   ! Always provide fresh default values
lopara =  (/  4       /)
owerte =  (/  0.0     /)
!
MAXW = MAX(MIN_PARA,MAXSCAT+1)
!                                                                       
lend = .false. 
CALL no_error 
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/powd' 
diff_lsingle = .FALSE.
!                                                                       
main: DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) THEN 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') THEN 
                  IF (length.ge.2) THEN 
                     line(1:length-1) = line(2:length)
                     length = 1
                     CALL file_kdo(line, length)
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----list asymmetric unit 'asym'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) THEN 
                  CALL show_asym 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) THEN 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----list atoms present in the crystal 'chem'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'chem', 2, lbef, 4) ) THEN 
                  CALL show_chem 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
                  CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) THEN                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus powd '//zeile, lp) 
                  ENDIF 
!                                                                       
!     switch to electrons diffraction 'electron'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'electron', 1, lbef, 8) ) THEN 
                  lxray = .true. 
                  diff_radiation = RAD_ELEC
                  diff_table = RAD_INTER
!                                                                       
!     switch to neutron diffraction 'neut'                              
!                                                                       
               ELSEIF (str_comp (befehl, 'neut', 1, lbef, 4) ) THEN 
                  lxray = .false. 
                  diff_radiation = RAD_NEUT
                  diff_table = RAD_INTER
!
!     ----rese powder patter settings 'rese'
!
               ELSEIF (str_comp (befehl, 'rese', 2, lbef, 4) ) THEN 
                  CALL powder_reset
!                                                                       
!     ----run transformation 'run'                                      
!                                                                       
               ELSEIF(str_comp(befehl, 'run ', 2, lbef, 4)) THEN 
                  CALL dlink(ano, lambda, rlambda, renergy, l_energy, &
                             diff_radiation, diff_table, diff_power) 
                  IF(ier_num.eq.0) THEN 
                     CALL phases_set(zeile, lp)
                     CALL pow_conv_limits
                     pow_qmin_u = pow_qmin   !! save user values
                     pow_qmax_u = pow_qmax
                     pow_deltaq_u=pow_deltaq
                     fwhm = SQRT(MAX(ABS(pow_u*pow_qmax**2 + pow_v*pow_qmax + pow_w), 0.00001D0) ) 
                     pow_qmax = pow_qmax + 2.0D0*pow_width*fwhm
                     pow_ds_max = (pow_qmax+pow_deltaq)/(zpi)
                     pow_ds_min = pow_qmin/(zpi)
                     check_dq: DO 
                        IF(NINT( (pow_qmax-pow_qmin)/pow_deltaq+1)>2**18) THEN
                           pow_deltaq = pow_deltaq*2.0D0
                        ELSE
                           EXIT check_dq
                        ENDIF
                     ENDDO check_dq
                     CALL powder_qcheck
                     IF(ier_num==0) THEN
                        IF(.NOT.pha_multi) pha_frac(1) = 1.0E0
!
                     IF(pow_four_type.eq.POW_DEBYE) THEN 
                        CALL pow_pdf_hist
                     ELSE
                        CALL powder_complete 
                     ENDIF
                     IF(ier_num == 0) THEN
                        four_was_run = .true.
                        CALL powder_convolute   ! convolute with profile
                        CALL phases_place       ! Copy current powder pattern into proper phase entry
                        IF (pow_four_type.eq.POW_DEBYE) THEN 
                           four_last = POWD_DY
                           if(pow_lperiod) call pow_pdf_hist_prep_period
                        ELSE
                           four_last = POWD_CO
                        ENDIF 
                     ENDIF 
                     ENDIF 
                     pow_qmax = pow_qmax_u ! Restore user settings
                  ENDIF 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
                  CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                              diff_radiation, diff_table, diff_power) 
                  IF(str_comp(cpara(1), 'scat', 4, lpara(1), 4)) THEN
                     CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                                 diff_radiation, diff_table, diff_power)
                     IF (ier_num.ne.0) THEN
                       RETURN
                     ENDIF
                     CALL get_params(zeile, ianz, cpara, lpara, UBOUND(cpara,1), lp)
                     IF(ier_num==0) &
                     CALL do_show_scat(ianz, cpara, lpara, werte, UBOUND(cpara,1))
                  ELSE
                     CALL pow_conv_limits
                     CALL pow_show 
                  ENDIF
!                                                                       
!------- -Set values 'set'                                              
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) THEN 
                  CALL do_pow_set (zeile, lp) 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
                  IF (zeile.ne.' ') THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
                  CALL do_input (zeile, lp) 
!                                                                       
!     switch to x-ray diffraction 'xray'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'xray', 1, lbef, 4) ) THEN 
                  lxray = .true. 
                  diff_radiation = RAD_XRAY
                  diff_table = RAD_INTER
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp)
                  if(ianz>0) then
                  CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                       oname, loname, opara, lopara, lpresent, owerte)
                  endif
                  if(opara(O_TABLE)=='waas') diff_table= RAD_WAAS
!                                                                       
!     ------unknown command                                             
!                                                                       
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in powder menu'
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
               CYCLE main
         ELSE
               EXIT main
         ENDIF 
      ENDDO  main
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE do_powder                      
!*****7*****************************************************************
      SUBROUTINE pow_show 
!-                                                                      
!     Prints summary of powder diffraction settings                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE metric_mod
      USE powder_mod 
      USE trig_degree_mod
      USE wink_mod
!
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(8) radiation 
      CHARACTER (LEN=8), DIMENSION(3), PARAMETER :: c_rad = (/ &
         'X-ray   ', 'neutron ', 'electron' /)
      CHARACTER(14) cfour (0:1) 
      CHARACTER(28) ccalc (0:5) 
      CHARACTER(21) cpref (1:2) 
      CHARACTER(29) cprofile (0:2) 
!                                                                       
!     REAL hkl (3), h1, h2, dstar, ttheta 
!     REAL del_tth_min (3) 
!     REAL del_tth_max (3) 
!                                                                       
!     REAL skalpro 
!     REAL asind 
!     REAL sind 
!                                                                       
      DATA cfour / 'normal Fourier', 'Stacking fault' / 
      DATA ccalc / 'rez. space integration      ', 'Debye formula       &
     &        ', 'Debye formula long          ', 'Debye formula fast    &
     &      ', 'Debye formula via histogram ', 'new integration       ' &
     &/                                                                 
      DATA cpref / 'Rietveld Toraya model', 'Modified March model ' / 
      DATA cprofile / 'Profile function switched off', 'Gaussian        &
     &             ', 'Pseudo-Voigt                 ' /                 
!                                                                       
CALL pow_conv_limits
!
      WRITE (output_io, 1000) 
!     radiation = 'neutron' 
!     IF (lxray) radiation = 'x-ray' 
      radiation = c_rad(diff_radiation)
      IF (lambda.eq.' ') THEN 
         IF(diff_radiation==2) THEN
            WRITE (output_io, 1201) radiation, rlambda , renergy
         ELSE 
            WRITE (output_io, 1200) radiation, rlambda , renergy
         ENDIF 
      ELSE 
         WRITE (output_io, 1210) radiation, lambda, rlambda 
      ENDIF 
!                                                                       
!     IF (pow_axis.eq.POW_AXIS_DSTAR) THEN 
!        WRITE (output_io, 1211) 'dstar=2 sin(Theta)/lambda' 
!     ELSEIF (pow_axis.eq.POW_AXIS_TTH) THEN 
!        WRITE (output_io, 1211) '2-Theta' 
!     ELSEIF (pow_axis.eq.POW_AXIS_Q) THEN 
         WRITE (output_io, 1211) 'Q=4pi sin(Theta)/lambda' 
!     ELSE 
!        WRITE (output_io, 1211) 'Has not been defined!' 
!     ENDIF 
!     IF (pow_axis.eq.POW_AXIS_DSTAR.or.pow_axis.eq.POW_AXIS_TTH) THEN 
!        IF (rlambda.ne.0.0) THEN 
!           pow_ds_max = 2. * sind (pow_tthmax * 0.5) / rlambda 
!           pow_ds_min = 2. * sind (pow_tthmin * 0.5) / rlambda 
!        ENDIF 
!        WRITE (output_io, 1220) pow_tthmin, pow_tthmax 
!        WRITE (output_io, 1221) pow_ds_min, pow_ds_max 
!        WRITE (output_io, 1230) pow_deltatth 
!        WRITE (output_io, 1240) pow_hkl_del 
!                                                                       
!        hkl (1) = 1.0 
!        hkl (2) = 0.0 
!        hkl (3) = 0.0 
!        h1 = pow_ds_min / sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        h2 = h1 + pow_hkl_del (1) 
!        hkl (1) = h2 
!        dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!        del_tth_min (1) = abs (ttheta - pow_tthmin) 
!        hkl (1) = 1.0 
!        hkl (2) = 0.0 
!        hkl (3) = 0.0 
!        h1 = pow_ds_max / sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        h2 = h1 - pow_hkl_del (1) 
!        hkl (1) = h2 
!        dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!        del_tth_max (1) = abs (ttheta - pow_tthmax) 
!                                                                       
!        hkl (2) = 1.0 
!        hkl (1) = 0.0 
!        hkl (3) = 0.0 
!        h1 = pow_ds_min / sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        h2 = h1 + pow_hkl_del (2) 
!        hkl (2) = h2 
!        dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!        del_tth_min (2) = abs (ttheta - pow_tthmin) 
!        hkl (2) = 1.0 
!        hkl (1) = 0.0 
!        hkl (3) = 0.0 
!        h1 = pow_ds_max / sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        h2 = h1 - pow_hkl_del (2) 
!        hkl (2) = h2 
!        dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!        del_tth_max (2) = abs (ttheta - pow_tthmax) 
!                                                                       
!        hkl (3) = 1.0 
!        hkl (2) = 0.0 
!        hkl (1) = 0.0 
!        h1 = pow_ds_min / sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        h2 = h1 + pow_hkl_del (3) 
!        hkl (3) = h2 
!        dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!        del_tth_min (3) = abs (ttheta - pow_tthmin) 
!        hkl (3) = 1.0 
!        hkl (2) = 0.0 
!        hkl (1) = 0.0 
!        h1 = pow_ds_max / sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        h2 = h1 - pow_hkl_del (3) 
!        hkl (3) = h2 
!        dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!        ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!        del_tth_max (3) = abs (ttheta - pow_tthmax) 
!                                                                       
!        WRITE (output_io, 1245) 
!        WRITE (output_io, 1250) del_tth_min 
!        WRITE (output_io, 1260) del_tth_max 
!        WRITE (output_io, 1270) pow_hkl_shift 
!        WRITE (output_io, 1240) pow_hkl_del 
!                                                                       
!        WRITE (output_io, 2100) cprofile (pow_profile) 
!        IF (pow_profile.eq.0) THEN 
!           CONTINUE 
!        ELSEIF (pow_profile.eq.POW_PROFILE_GAUSS) THEN 
!           WRITE (output_io, 1235) pow_delta 
!           WRITE (output_io, 2123) pow_width 
!        ELSEIF (pow_profile.eq.POW_PROFILE_PSVGT) THEN 
!           WRITE (output_io, 2120) pow_u, pow_v, pow_w 
!           WRITE (output_io, 2121) pow_eta, pow_eta_l , pow_eta_q
!           WRITE (output_io, 2122) pow_p1, pow_p2, pow_p3, pow_p4 
!           WRITE (output_io, 2123) pow_width 
!        ENDIF 
!                                                                       
!     ELSEIF (pow_axis.eq.POW_AXIS_Q) THEN 
         WRITE (output_io, 1290) pow_qmin, pow_qmax 
         WRITE (output_io, 1291) pow_deltaq 
         WRITE (output_io, 1292) pow_qzero 
         WRITE (output_io, 1220) pow_tthmin, pow_tthmax 
         WRITE (output_io, 1230) pow_deltatth 
         WRITE (output_io, 1231) pow_tthzero 
         WRITE (output_io, 2100) cprofile (pow_profile) 
         IF (pow_profile.eq.0) THEN 
            CONTINUE 
         ELSEIF (pow_profile.eq.POW_PROFILE_GAUSS) THEN 
            WRITE (output_io, 1235) pow_delta 
            WRITE (output_io, 2123) pow_width 
         ELSEIF (pow_profile.eq.POW_PROFILE_PSVGT) THEN 
            WRITE (output_io, 2120) pow_u, pow_v, pow_w 
            WRITE (output_io, 2121) pow_eta, pow_eta_l, pow_eta_q
            WRITE (output_io, 2122) pow_asym(:,1)
            WRITE (output_io, 2123) pow_asym(:,2)
            WRITE (output_io, 2124) pow_asym(:,3)
            WRITE (output_io, 2125) pow_width 
         ENDIF 
!     ENDIF 
!                                                                       
      IF (ldbw) THEN 
         WRITE (output_io, 1300) 'used' 
      ELSE 
         WRITE (output_io, 1300) 'ignored' 
      ENDIF 
!                                                                       
      IF (ano) THEN 
         WRITE (output_io, 1310) 'used' 
      ELSE 
         WRITE (output_io, 1310) 'ignored' 
      ENDIF 
!                                                                       
      IF (pow_pref) THEN 
         WRITE (output_io, 1340) cpref (pow_pref_type) 
         WRITE (output_io, 1341) pow_pref_g1, pow_pref_g2 
         WRITE (output_io, 1342) pow_pref_hkl 
      ELSE 
         WRITE (output_io, 1340) 'off' 
      ENDIF 
!                                                                       
      WRITE (output_io, 1400) cfour (pow_four_mode) 
      IF (pow_four_mode.eq.0.or.pow_four_mode.eq.5) THEN 
         IF (pow_l_all) THEN 
            WRITE (output_io, 1450) 
         ELSE 
            WRITE (output_io, 1455) 
         ENDIF 
      ELSEIF(pow_four_mode==POW_DEBYE) THEN
         IF(pow_lperiod) THEN
            WRITE(output_io,'(A)' ) '   Create periodic powder/PDF '
            WRITE(output_io, 1460) pow_period
         ELSE
            WRITE(output_io,'(A)' ) '   Create finite   powder/PDF '
         ENDIF 
      ENDIF 
      IF (pow_lp.eq.POW_LP_NONE) THEN 
         WRITE (output_io, 1500) 
      ELSEIF (pow_lp.eq.POW_LP_BRAGG) THEN 
         WRITE (output_io, 1510) pow_lp_fac, pow_lp_ang 
      ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
         WRITE (output_io, 1520) 
      ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
         WRITE (output_io, 1530) pow_lp_fac, pow_lp_ang 
      ENDIF 
      WRITE (output_io, 1600) ccalc (pow_four_type) 
!                                                                       
 1000 FORMAT    ( ' Settings for Powder Diffraction segment :') 
 1200 FORMAT    ( '   Radiation               : ',A,', wavelength = ',  &
     &          F7.4,' A == ', F9.4,' keV')                               
 1201 FORMAT    ( '   Radiation               : ',A,', wavelength = ',  &
     &          F7.4,' A == ', F9.4,' meV')
 1210 FORMAT    ( '   Radiation               : ',A,', wavelength = ',A4,   &
     &                    ' = ',F8.4,' A')                              
 1211 FORMAT    ( '   Calculations for axis   : ',a) 
 1220 FORMAT    ( '   TTHmin, TTHmax          : ',f10.5,2x,f10.5) 
 1221 FORMAT    ( '   d* min, d* max          : ',f10.5,2x,f10.5) 
 1230 FORMAT    ( '   DELTA TTH               : ',f10.5) 
 1231 FORMAT    ( '   TTH zero; Obs0 at true  : ',f10.5) 
 2100 FORMAT    ( '   Instrument resolution   : ',A) 
 1235 FORMAT( '   Instrument resolution   : ',f10.5,2x,'(0.0 = off)') 
 2120 FORMAT    ( '       Profile U,V,W       : ',f10.5,2x,f10.5,2x,    &
     &                                                   f10.5)         
 2121 FORMAT    ( '       Profile Eta, q, l   : ',f10.5,2x,f10.5, 2x,f10.5) 
 2122 FORMAT    ( '       Profile asymmetry   : ',4(f10.5,2x)) 
 2123 FORMAT    ( '       Profile asym linear : ',4(f10.5,2x)) 
 2124 FORMAT    ( '       Profile asym square : ',4(f10.5,2x)) 
 2125 FORMAT    ( '       Profile width *FWHM : ',1(f10.5,2x)) 
 1240 FORMAT    ( '   dH, dK, dL              : ',3(f10.5,2x)) 
 1245 FORMAT    ( '   Corr. steps in TTH') 
 1250 FORMAT    ( '   at TTHmin               : ',3(f10.5,2x)) 
 1260 FORMAT    ( '   at TTHmax               : ',3(f10.5,2x)) 
 1270 FORMAT    ( '   shift for dH,dK,dL      : ',3(f10.5,2x)) 
 1290 FORMAT    ( '   Q  min, Q  max          : ',f10.5,2x,f10.5) 
 1291 FORMAT    ( '   DELTA Q                 : ',f10.5) 
 1292 FORMAT    ( '   Q zero (Obs0 at true Q) : ',f10.5) 
 1300 FORMAT    ( '   Temp. factors           : ',A) 
 1310 FORMAT    ( '   Anomalous scat.         : ',A) 
 1340 FORMAT    ( '   Preferred orientation   : ',A) 
 1341 FORMAT    ( '   Preferred damp,portion  : ',2(f10.5,2x)) 
 1342 FORMAT    ( '   Preferred axis hkl      : ',3(f10.5,2x)) 
 1400 FORMAT    ( '   Fourier Mode            : ',A) 
 1450 FORMAT    ( '       Bragg Reflections   : ','included') 
 1455 FORMAT    ( '       Bragg Reflections   : ','excluded') 
 1460 FORMAT    ( '   Radius for periodicity  : ', F10.5, 'A')
 1500 FORMAT    ( '   Powder diffractometer   : ','none specified',     &
     &                   ' no Lorentz/Polarisation effect calculated')  
 1510 FORMAT    ( '   Powder diffractometer   : ','Bragg-Brentano',/    &
     &                   '   Polarisation fraction   : ',f10.5,/        &
     &                   '   2-Theta Monochromator   : ',f10.5)         
 1520 FORMAT    ( '   Powder diffractometer   : ',                      &
     &                   'Neutron Debye-Scherrer')                      
 1530 FORMAT    ( '   Powder diffractometer   : ','3-axis Synchrotron',/&
     &                   '   Polarisation fraction   : ',f10.5,/        &
     &                   '   2-Theta Monochromator   : ',f10.5)         
 1600 FORMAT    ( '   Fourier calculation via : ',A) 
      END SUBROUTINE pow_show                       
!
!*****7*****************************************************************
!
SUBROUTINE do_pow_set(zeile, lcomm) 
!-                                                                      
!     Set various paramters for the powder diffraction                  
!+                                                                      
USE discus_config_mod 
USE debye_mod 
USE diffuse_mod 
USE phases_mod
USE pdf_mod
USE powder_mod 
USE ber_params_mod
USE get_params_mod
USE lib_errlist_func
USE precision_mod
USE trig_degree_mod
USE str_comp_mod
USE string_convert_mod
USE take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 7 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile 
INTEGER         , INTENT(INOUT) :: lcomm 
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_RCUT    = 1
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
!
DATA oname  / 'rcut' /
DATA loname /  4     /
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile)))  :: cpara (MAXW) 
CHARACTER(LEN=PREC_STRING) :: symbol
INTEGER :: lpara (MAXW) 
INTEGER :: lsymbol
INTEGER :: ianz 
INTEGER :: i 
REAL(KIND=PREC_DP) :: werte (MAXW) 
!                                                                       
!     REAL cosd 
!                                                                       
opara  =  (/ '0.0000'/)   ! Always provide fresh default values
lopara =  (/  6      /)
owerte =  (/  0.0    /)
CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
err_para: IF (ier_num.eq.0) THEN 
!
      IF(str_comp(cpara(1), 'axis', 2, lpara(1), 4)) THEN 
         pow_axis = POW_AXIS_Q 
      ELSEIF (str_comp (cpara (1) , 'back', 2, lpara (1) , 4) ) THEN 
         IF (ianz.ge.2) THEN 
            cpara (1) = '0' 
            lpara (1) = 1 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) THEN 
               DO i = 2, ianz 
               pow_back (i - 2) = werte (i) 
               ENDDO 
               pow_nback = ianz - 2
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSEIF(str_comp(cpara(1), 'bragg', 2, lpara(1), 5)) THEN
         IF (ianz.eq.2) THEN 
            IF (str_comp(cpara(2), 'incl', 1, lpara(2), 4) ) THEN
               pow_l_all = .true. 
            ELSEIF(str_comp(cpara(2), 'excl', 1, lpara(2), 4) ) THEN
               pow_l_all = .false. 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSEIF(str_comp (cpara (1) , 'disp', 2, lpara (1) , 4) ) THEN 
         IF(ianz.eq.2) THEN 
            IF(str_comp(cpara(2), 'anom', 1, lpara(2), 4) ) THEN
               ano = .true. 
            ELSEIF(str_comp (cpara (2) , 'off', 1, lpara (2) , 3) ) THEN
               ano = .false. 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
   ELSEIF(str_comp(cpara(1), 'corrlin', 5, lpara(1), 7)) THEN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      err_opt1: IF (ier_num.eq.0) THEN 
         IF (ianz == 2) THEN 
            cpara (1) = '0' 
            lpara (1) = 1 
            CALL ber_params(ianz, cpara, lpara, werte, maxw) 
            IF(ier_num == 0) THEN 
               pdf_clin_a = werte(2) 
               pdf_rcut   = owerte(O_RCUT)
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF  err_opt1
   ELSEIF(str_comp(cpara(1), 'corrquad', 5, lpara(1), 8)) THEN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      err_opt2: IF (ier_num.eq.0) THEN 
         IF (ianz == 2) THEN 
            cpara (1) = '0' 
            lpara (1) = 1 
            CALL ber_params(ianz, cpara, lpara, werte, maxw) 
            IF(ier_num == 0) THEN 
               pdf_cquad_a = werte(2) 
               pdf_rcut   = owerte(O_RCUT)
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF  err_opt2
   ELSEIF(str_comp(cpara(1), 'delta', 2, lpara(1), 5)) THEN
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_delta = werte (2) 
                  pow_profile = POW_PROFILE_GAUSS 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'dtth', 2, lpara (1) , 4) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  IF(werte(2)<0.0) THEN
                     ier_num = -107
                     ier_typ = ER_APPL
                     ier_msg(1) = '2Theta step must be positive!'
                  ELSEIF(werte(2)>180.0) THEN
                     ier_num = -107
                     ier_typ = ER_APPL
                     ier_msg(1) = '2Theta step must be less than 180degrees'
                  ELSE
                     pow_deltatth = werte (2) 
                     pow_deltaqtth = .FALSE.
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'dq', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_deltaq = NINT(werte(2)*1.D5)/1.D5
                  pow_deltaqtth = .TRUE.
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'dh', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_hkl_del (1) = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'dk', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_hkl_del (2) = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'dl', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_hkl_del (3) = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     set the energy of the radiation to be used 'energy'                             
!                                                                       
         ELSEIF (str_comp (cpara(1), 'energy', 2, lpara(1), 6) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               renergy  = werte (2) 
               lambda   = ' ' 
               l_energy = .true.
            ELSE 
               ier_num = -6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     set the Ka2/Ka1 ratio 'ka21' ~ 0.50
!                                                                       
         ELSEIF (str_comp (cpara(1), 'ka21', 4, lpara(1), 4) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               pow_ka21   = werte (2) 
               pow_ka21_u = .TRUE.
            ELSE 
               ier_num = -6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     set the periodic radius 'period' 
!                                                                       
         ELSEIF(str_comp(cpara(1), 'period', 4, lpara(1), 6)) THEN 
            IF(ianz.eq.2) THEN 
               cpara(1) = '0' 
               lpara(1) = 1 
               IF(str_comp(cpara(2), 'off', 3, lpara(2), 3)) THEN
                  pow_period  = 0.00
                  pow_lperiod = .FALSE.
               ELSE
                  CALL ber_params(ianz, cpara, lpara, werte, maxw) 
                  IF(ier_num == 0) THEN
                     pow_period  = werte(2) 
                     pow_lperiod = .TRUE.
                  ELSE 
                     ier_num = -6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = -6 
               ier_typ = ER_COMM 
            ENDIF 
!
         ELSEIF (str_comp (cpara (1) , 'pref', 2, lpara (1) , 4) ) THEN 
            IF (ianz.ge.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               IF (str_comp (cpara (2) , 'off', 2, lpara (2) , 3) )     &
               THEN                                                     
                  pow_pref = .false. 
               ELSE 
                  IF (str_comp (cpara (2) , 'riet', 2, lpara (2) , 4) ) &
                  THEN                                                  
                     pow_pref_type = POW_PREF_RIET 
                     pow_pref = .true. 
                  ELSEIF (str_comp (cpara (2) , 'march', 2, lpara (2) , &
                  4) ) THEN                                             
                     pow_pref_type = POW_PREF_MARCH 
                     pow_pref = .true. 
                  ELSEIF (str_comp (cpara (2) , 'damp', 2, lpara (2) ,  &
                  4) .or.str_comp (cpara (2) , 'g1', 2, lpara (2) , 2) )&
                  THEN                                                  
                     cpara (2) = '0' 
                     lpara (2) = 1 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) THEN 
                        pow_pref_g1 = werte (3) 
                     ENDIF 
                  ELSEIF (str_comp (cpara (2) , 'portion', 2, lpara (2) &
                  , 7) .or.str_comp (cpara (2) , 'g2', 2, lpara (2) , 2)&
                  ) THEN                                                
                     cpara (2) = '0' 
                     lpara (2) = 1 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) THEN 
                        pow_pref_g2 = werte (3) 
                     ENDIF 
                  ELSEIF (str_comp (cpara (2) , 'hkl', 2, lpara (2) , 7)&
                  ) THEN                                                
                     cpara (2) = '0' 
                     lpara (2) = 1 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) THEN 
                        pow_pref_hkl (1) = werte (3) 
                        pow_pref_hkl (2) = werte (4) 
                        pow_pref_hkl (3) = werte (5) 
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
         ELSEIF (str_comp (cpara (1) , 'profile', 2, lpara (1) , 7) )   &
         THEN                                                           
            IF (ianz.ge.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               IF(str_comp(cpara(2), 'off', 2, lpara(2), 3)) THEN
                  pow_profile = 0 
                  pow_delta = 0.0 
                  pow_eta = 0.5 
               ELSEIF(str_comp(cpara(2), 'gauss', 2, lpara(2), 5)) THEN
                  pow_profile = POW_PROFILE_GAUSS 
               ELSEIF(str_comp(cpara(2), 'pseudo', 2, lpara(2) , 6)) THEN
                  pow_profile = POW_PROFILE_PSVGT 
               ELSEIF(str_comp(cpara(2), 'eta', 2, lpara(2), 3)) THEN
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  pow_eta   = 0.5 
                  pow_eta_l = 0.0 
                  pow_eta_q = 0.0 
                  IF (ier_num == 0) THEN 
                     pow_eta = werte(3) 
                     IF(ianz == 4) THEN 
                        pow_eta_l = werte(4) 
                        pow_eta_q = 0.0 
                     ELSEIF(ianz == 5) THEN 
                        pow_eta_l = werte(4) 
                        pow_eta_q = werte(5) 
                     ELSE 
                        pow_eta_l = 0.0 
                        pow_eta_q = 0.0 
                     ENDIF 
                  ENDIF 
               ELSEIF(str_comp(cpara(2), 'uvw', 2, lpara(2), 3)) THEN
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     pow_u = werte (3) 
                     pow_v = werte (4) 
                     pow_w = werte (5) 
                  ENDIF 
               ELSEIF(str_comp(cpara(2), 'asym', 2, lpara(2), 4)  &
                      .AND.    lpara(2)<=4) THEN                                                     
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     pow_asym(:,1) = werte (3:6) 
                  ENDIF 
               ELSEIF(str_comp(cpara(2), 'asym_l', 2, lpara(2), 6)) THEN
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     pow_asym(:,2) = werte (3:6) 
                  ENDIF 
               ELSEIF(str_comp(cpara(2), 'asym_q', 2, lpara(2), 6)) THEN
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     pow_asym(:,3) = werte (3:6) 
                  ENDIF 
               ELSEIF (str_comp (cpara (2) , 'width', 2, lpara (2) , 5) &
               ) THEN                                                   
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     pow_width = werte (3) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF(str_comp (cpara (1) , 'scale', 2, lpara (1) , 5) ) THEN
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_scale = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'sh', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_hkl_shift (1) = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'sk', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_hkl_shift (2) = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'stepr', 2, lpara (1) , 5) )     &
         THEN                                                           
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_del_hist = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'sl', 2, lpara (1) , 2) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_hkl_shift (3) = werte (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     Switch Fourier type between normal Fourier and DEBYE calculation  
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'calc', 1, lpara (1) , 4) ) THEN 
            IF (ianz.ge.2) THEN 
               IF (str_comp (cpara (2) , 'comp', 1, lpara (2) , 4) )    &
               THEN                                                     
                  pow_four_type = POW_COMPL 
               ELSEIF (str_comp (cpara (2) , 'debye', 1, lpara (2) , 5) ) THEN                                                   
                     pow_four_type = POW_DEBYE 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ Switch Fourier mode between normal Fourier and Stacking         
!       Fault 'four'                                                    
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'four', 1, lpara (1) , 4) ) THEN 
            IF (ianz.eq.2) THEN 
               IF (str_comp (cpara (2) , 'four', 1, lpara (2) , 4) )    &
               THEN                                                     
                  pow_four_mode = POW_FOURIER 
               ELSEIF (str_comp (cpara (2) , 'stack', 1, lpara (2) , 5) &
               ) THEN                                                   
                  pow_four_mode = POW_STACK 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     Switch application of LP correction to proper instrument 'lpcor'  
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'lpcor', 1, lpara (1) , 5) ) THEN                                                           
            IF (ianz.ge.2) THEN 
               IF(str_comp(cpara(2), 'bragg', 4, lpara(2), 5)) THEN
                  pow_lp = POW_LP_BRAGG 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  IF(ianz == 1) THEN 
                     CALL ber_params(ianz, cpara, lpara, werte, maxw) 
                     IF(ier_num == 0) THEN 
                        pow_lp_ang = werte (1) 
                        pow_lp_fac = (cosd(pow_lp_ang))**2 
                     ELSE 
                        RETURN 
                     ENDIF 
                  ELSEIF(ianz == 0) THEN 
                     pow_lp_fac = 0.5          ! Not polarised
                     pow_lp_ang = 0.0
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSEIF (str_comp (cpara (2) , 'corre', 4, lpara (2) , 5) )   &
               THEN                                                     
                  pow_lp = POW_LP_CORRE 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  IF (ianz.eq.1) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) THEN 
                        pow_lp_ang = werte (1) 
                        pow_lp_fac = 0.5
                     ELSE 
                        RETURN 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSEIF (str_comp (cpara (2) , 'neutron', 4, lpara (2) ,  &
               7) ) THEN                                                
                  pow_lp = POW_LP_NEUT 
               ELSEIF (str_comp (cpara (2) , 'none', 4, lpara (2) , 4) )&
               THEN                                                     
                  pow_lp = POW_LP_NONE 
               ELSEIF (str_comp (cpara (2) , 'synchrotron', 4, lpara (2)&
               , 11) ) THEN                                             
                  pow_lp = POW_LP_SYNC 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  IF (ianz.eq.1.or.ianz.eq.2) THEN 
                     werte (2) = 0.0 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) THEN 
                        pow_lp_fac = werte (1) 
                        pow_lp_ang = werte (2) 
                        pow_lp_cos = abs (cosd (pow_lp_ang) ) 
                     ELSE 
                        RETURN 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 90 
                  ier_TYP = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     Switch usage of temperature coefficients on/off 'temp'            
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'temp', 1, lpara (1) , 4) ) THEN 
            IF (ianz.eq.2) THEN 
               IF (str_comp (cpara (2) , 'igno', 1, lpara (2) , 4) )    &
               THEN                                                     
                  ldbw = .false. 
               ELSEIF (str_comp (cpara (2) , 'use', 1, lpara (2) , 3) ) &
               THEN                                                     
                  ldbw = .true. 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF(str_comp(cpara(1), 'tthmax', 5, lpara(1), 6)) THEN
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  IF(werte(2)<0.0) THEN
                     ier_num = -107
                     ier_typ = ER_APPL
                     ier_msg(1) = '2Theta max must be positive!'
                  ELSEIF(werte(2)>180.0) THEN
                     ier_num = -107
                     ier_typ = ER_APPL
                     ier_msg(1) = '2Theta max must be less than 180degrees'
                  ELSE
                     pow_tthmax = werte (2) 
                     pow_qtthmax = .FALSE.
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF(str_comp(cpara(1), 'tthmin', 5, lpara(1), 6)) THEN
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  IF(werte(2)<0.0) THEN
                     ier_num = -107
                     ier_typ = ER_APPL
                     ier_msg(1) = '2Theta min must be positive!'
                  ELSEIF(werte(2)>180.0) THEN
                     ier_num = -107
                     ier_typ = ER_APPL
                     ier_msg(1) = '2Theta min must be less than 180degrees'
                  ELSE
                     pow_tthmin = werte (2) 
                     pow_qtthmin = .FALSE.
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF(str_comp(cpara(1), 'tthzero', 4, lpara(1), 4)) THEN 
            IF(ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_tthzero    = werte(2) 
                  pow_qtthzero = .TRUE.
               ENDIF 
            ELSE 
               ier_num = -6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'qmax', 4, lpara (1) , 4) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_qmax = werte (2) 
                  pow_qtthmax = .TRUE.
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'qmin', 4, lpara (1) , 4) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_qmin = werte (2) 
                  pow_qtthmin = .TRUE.
               ENDIF 
            ELSE 
               ier_num = -6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF(str_comp(cpara(1), 'qzero', 4, lpara(1), 4)) THEN 
            IF(ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pow_qzero    = werte(2) 
                  pow_qtthzero = .TRUE.
               ENDIF 
            ELSE 
               ier_num = -6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     set the wave length to be used 'wvle'                             
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'wvle', 1, lpara (1) , 4) ) THEN 
            IF (ianz.eq.2) THEN 
               cpara (1) = '0' 
               lpara (1) = 1 
               symbol    = cpara(2)
               lsymbol   = lpara(2)
               CALL do_cap (symbol) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  rlambda = werte (2) 
                  lambda = ' ' 
                  l_energy = .false.
               ELSEIF (ichar ('A') <=  ichar (symbol    (1:1) ) .AND.&
                       ichar (symbol    (1:1) ) <= ichar ('Z') ) THEN                 
                  lambda = symbol(1:lsymbol)  
                  l_energy = .false.
                  CALL no_error
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
   ELSE 
      ier_num = - 8 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF err_para
!                                                                       
END SUBROUTINE do_pow_set                     
!
!*****7*****************************************************************
!
SUBROUTINE powder_qcheck 
!-                                                                      
!     Perform common error checking for all powder modules
!+                                                                      
USE diffuse_mod 
USE powder_mod 
USE wink_mod
!
IMPLICIT none 
!                                                                       
IF (rlambda < TINY(0.0)) THEN
   ier_num = -99
   ier_typ = ER_APPL
   ier_msg(1) = 'Wave length is zero; check powder settings'
   RETURN
ENDIF 
!                                                                       
!      Perform error checking                                           
!                                                                       
IF (pow_hkl_del (1) .eq.0.and.pow_hkl_del (1)               &
      .eq.0.and.pow_hkl_del (1) .eq.0) THEN                       
   ier_num = - 106 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!
IF (pow_qmax.le.pow_qmin.or.pow_deltaq.le.0.0) THEN 
   ier_num = - 108 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!
IF(pow_qmax*rlambda/2./zpi > 1.0) THEN
   ier_num = -108
   ier_typ = ER_APPL
   ier_msg(1) = 'Qmax is too large for current wave length'
   ier_msg(2) = 'Qmax*lambda/(4pi) is greater than one!'
   ier_msg(3) = 'Reduce Qmax or the wave length'
   RETURN
ENDIF
!                                                                       
END SUBROUTINE powder_qcheck                     
!
!*****7*****************************************************************
!
SUBROUTINE powder_complete ()
!-                                                                      
!     Calculate global parameters and start the individual modes        
!+                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE diffuse_mod 
USE fourier_sup
USE metric_mod
USE phases_mod
USE phases_set_mod
USE powder_mod 
USE powder_scat_mod 
USE powder_tables_mod 
USE stack_menu
USE stack_mod
USE wink_mod
!                                                                       
USE param_mod 
USE prompt_mod 
USE precision_mod 
USE trig_degree_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=PREC_STRING) :: line 
INTEGER :: laenge 
INTEGER :: i 
INTEGER :: h_start, h_end 
INTEGER :: k_start, k_end 
INTEGER :: l_start=0, l_end=1
INTEGER :: ih, ik 
INTEGER                    :: n_qxy   = 1
INTEGER                    :: n_nscat = 1
INTEGER                    :: n_natom = 1
INTEGER                    :: n_pkt   = 1
INTEGER                    :: n_pha   = 1
INTEGER :: itth 
!
LOGICAL :: l_twoparts
LOGICAL :: l_ano 
LOGICAL :: l_hh_real 
LOGICAL :: l_kk_real 
LOGICAL :: l_ll_real 
LOGICAL                   :: calc_f2aver, rept_f2aver
REAL(kind=PREC_DP) :: llstart, llend 
REAL(kind=PREC_DP) :: llstart2=0.0, llend2=1.0
REAL(kind=PREC_DP):: xstart, xdelta   ! start/step in dstar for sinthea/lambda table
REAL(kind=PREC_DP) :: hh, kk, ll 
REAL(kind=PREC_DP) :: rr, rrr, rtm 
REAL(kind=PREC_DP), dimension(3) :: hkl (3) 
REAL(kind=PREC_DP) :: dstar , q !, ttheta
REAL(KIND=PREC_DP)  :: inten 
REAL(kind=PREC_DP) :: u (3), v (3), w_min (3), w_max (3) 
REAL(kind=PREC_DP) :: u2, vv, ww 
REAL(kind=PREC_DP) :: aaa, bbb, ccc 
REAL(kind=PREC_DP) :: llstartmini 
REAL(kind=PREC_DP) :: llendmini 
REAL(kind=PREC_DP) :: ss 
!                                                                       
!
st_new_form = .TRUE.    ! We need new form factors from stack to be placed into phases
n_qxy   = 1
n_nscat = 1
n_natom = 1
n_pkt   = 1
l_twoparts = .false.
calc_f2aver = .true.    ! Assume that we need form factors
rept_f2aver = .true.    ! Assume that we need to repeat them
!
pow_hkl_max (1) = cr_a0 (1) * pow_ds_max 
pow_hkl_max (2) = cr_a0 (2) * pow_ds_max 
pow_hkl_max (3) = cr_a0 (3) * pow_ds_max 
!                                                                       
!
ss = seknds (0.0D0) 
!                                                                       
!     Set Fourier definitions                                           
!                                                                       
inc (2) = 1 
DO i = 1, 3 
   vi (i, 1) = 0.0 
   vi (i, 2) = 0.0 
ENDDO 
vi (3, 1) = pow_hkl_del (3) 
four_log = .false. 
!
!
n_pkt = NINT((pow_qmax+pow_deltaq  -pow_qmin  )/pow_deltaq  ) + 2
n_nscat = MAX(UBOUND(pow_f2,2), MAXSCAT, DIF_MAXSCAT)
IF(n_pkt .gt. POW_MAXPKT .OR. n_nscat>POW_MAXSCAT) THEN
   CALL alloc_powder ( n_pkt, n_nscat )
ENDIF
!     reset powder diagramm                                             
!                                                                       
pow_qsp(:)    = 0.0D0   ! 0:POW_MAXPKT
pow_f2aver(:) = 0.0D0   ! 0:POW_MAXPKT
pow_faver2(:) = 0.0D0   ! 0:POW_MAXPKT
pow_nreal     = 0
pow_u2aver    = 0.0
!     DO i = 1, POW_MAXPKT 
!     pow_qsp (i) = 0.0 
!     ENDDO 
n_qxy   = MAX(n_qxy,inc(1) * inc(2),n_pkt,MAXQXY)
n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
n_natom = MAX(n_natom,cr_natoms,DIF_MAXAT)
IF (inc (1) * inc (2) .gt. MAXQXY  .OR.          &
    n_pkt             .gt. MAXQXY  .OR.          &
    cr_nscat>DIF_MAXSCAT              ) THEN
  CALL alloc_diffuse (n_qxy,  n_nscat,  n_natom )
  IF (ier_num /= 0) THEN
    RETURN
  ENDIF
ENDIF
!
IF(n_qxy > PHA_MAXPTS .OR. cr_nscat> PHA_MAXSCAT) THEN
   n_pha   = PHA_MAXPHA
   n_qxy   = MAX(PHA_MAXPTS,  n_qxy)
   n_nscat = MAX(PHA_MAXSCAT, cr_nscat)
   CALL alloc_phases(n_pha, n_qxy, n_nscat)
ENDIF
!                                                                       
!------ calculate complex exponent table, form factor table             
!                                                                       
pow_npkt = n_pkt            ! set actual number of powder data points
!
CALL four_cexpt 
CALL four_formtab
!
xstart = pow_qmin  /zpi
xdelta = pow_deltaq/zpi
CALL powder_stltab(n_qxy, xstart  ,xdelta    )   ! Really only needed for <f^2> and <f>^2 for F(Q) and S(Q)
!
CALL powder_getatoms 
!                                                                       
!     calculate global limits along h                                   
!                                                                       
!
l_ano = ano 
!
IF (l_ano) THEN 
   h_start = - int (pow_hkl_max (1) / pow_hkl_del (1) ) 
ELSE 
   h_start = 0 
ENDIF 
h_end = int (pow_hkl_max (1) / pow_hkl_del (1) ) 
!                                                                       
!     loop over h                                                       
!                                                                       
loop_h: DO ih = h_start, h_end 

      hh = ih * pow_hkl_del (1) + pow_hkl_shift (1) 
      hkl (1) = hh 
      l_hh_real = abs ( (nint (hh) - hh) / pow_hkl_del (1) ) .gt.0.51 
!                                                                       
!     --Produce output to keep the user informed                        
!                                                                       
      WRITE (output_io, 5000) hh, pow_hkl_del (1), pow_hkl_max (1) 
!                                                                       
!     --NEW LIMITS                                                      
!                                                                       
      WRITE (line, 6000) hh 
 6000 FORMAT      (f12.6,',0,0, 1,0,0, rdrr') 
      laenge = 29 
      CALL do_proj (line, laenge) 
      u (1) = res_para (1) 
      u (2) = res_para (2) 
      u (3) = res_para (3) 
      u2 = skalpro (u, u, cr_rten) 
      line = '0,1,0, 1,0,0, ddrr' 
      laenge = 18 
      CALL do_proj (line, laenge) 
      v (1) = res_para (4) 
      v (2) = res_para (5) 
      v (3) = res_para (6) 
      vv = sqrt (skalpro (v, v, cr_rten) ) 
      v (1) = v (1) / vv 
      v (2) = v (2) / vv 
      v (3) = v (3) / vv 
      rr = pow_ds_max**2 - u2 
!DBG_RBN                                                                
!DBG      write(*,*) ' Vector u ',u                                     
!DBG      write(*,*) ' Vector v ',v                                     
!DBG      write(*,*) ' x*       ',rr                                    
      IF (rr.ge.0) THEN 
         w_min (1) = u (1) - sqrt ( (rr) ) * v (1) 
         w_min (2) = u (2) - sqrt ( (rr) ) * v (2) 
         w_min (3) = u (3) - sqrt ( (rr) ) * v (3) 
         ww = sqrt (skalpro (w_min, w_min, cr_rten) ) 
!DBG        write(*,'(a,5f10.4)') ' k minimum ',w_min, ww,pow_ds_max    
         w_max (1) = u (1) + sqrt ( (rr) ) * v (1) 
         w_max (2) = u (2) + sqrt ( (rr) ) * v (2) 
         w_max (3) = u (3) + sqrt ( (rr) ) * v (3) 
         ww = sqrt (skalpro (w_max, w_max, cr_rten) ) 
!DBG        write(*,'(a,5f10.4)') ' k maximum ',w_max, ww,pow_ds_max    
!DBG        write(line,'(a4,f6.2)') 'hkl.',hh                           
!DBG        do i=5,14                                                   
!DBG          if(line(i:i).eq. ' ') THEN                                
!DBG            line(i:i) = '0'                                         
!DBG          endif                                                     
!DBG        ENDDO                                                       
!DBG        open(15,file=line,status='unknown')                         
!DBG        line(1:3) = 'edg'                                           
!DBG        open(16,file=line,status='unknown')                         
!DBG        line(1:3) = 'HHH'                                           
!DBG        open(17,file=line,status='unknown')                         
!DBG        line(1:3) = 'EEE'                                           
!DBG        open(18,file=line,status='unknown')                         
      ENDIF 
      IF (rr.gt.0) THEN 
         eck (1, 1) = hh 
         eck (1, 2) = hh 
         eck (1, 3) = hh 
         IF (.not.l_ano.and.ih.eq.0) THEN 
            k_start = 0 
         ELSE 
            k_start = int( w_min (2) / pow_hkl_del (2) )
         ENDIF 
         k_end = int( w_max (2) / pow_hkl_del (2) )
!DBG_RBN                                                                
!DBG      k_start = -6.00000/pow_hkl_del(2)                             
!DBG      k_end   =  6.00000/pow_hkl_del(2)                             
!DBG      k_start = w_min(2)/pow_hkl_del(2)                             
!DBG      k_end   = w_max(2)/pow_hkl_del(2)                             
!DBGXXX      if(hkl(1).eq.-6.) THEN                                     
!DBG        write(*,*) 'k_start k_end',k_start,k_end                    
!DBG        write(*,*)                                                  
!DBGXXX      endif                                                      
!                                                                       
!     ----Start loop over k                                             
!                                                                       
         DO ik = k_start, k_end 
!DBGXXX      if(hkl(1).eq.0.0) THEN                                     
!DBGXXX        write(*,*) ' ik*pow_hkl_del(2) ',ik*pow_hkl_del(2)       
!DBGXXX        write(*,*) ' pow_hkl_shift(2)  ',pow_hkl_shift(2)        
!DBGXXX      endif                                                      
         kk = ik * pow_hkl_del (2) + pow_hkl_shift (2) 
         hkl (2) = kk 
         l_kk_real = abs ( (nint (kk) - kk) / pow_hkl_del (2) )         &
         .gt.0.51                                                       
!                                                                       
!     --Produce output to keep the user informed                        
!                                                                       
!       write (output_io,5010) kk,pow_hkl_del(2),pow_hkl_max(2)         
!                                                                       
         aaa = cr_rten (3, 3) 
         bbb = 2 * (cr_rten (1, 3) * hkl (1) + cr_rten (2, 3) * hkl (2) &
         )                                                              
         ccc = cr_rten (1, 1) * hkl (1) **2 + 2. * cr_rten (1, 2)       &
         * hkl (1) * hkl (2) + cr_rten (2, 2) * hkl (2) **2 -           &
         pow_ds_max**2                                                  
         rrr = bbb**2 - 4. * aaa * ccc 
!DBG          if(.not.l_ano .and. ih.eq.0 .and. ik.eq.0) THEN           
!DBG            write(*,*) ' rrr ',rrr                                  
!DBG          endif                                                     
         IF (rrr.ge.0) THEN 
!DBGXXX        write(*,*) ' hkl ',hkl                                   
!DBGXXX        write(*,'(a,3f10.4)') ' aaa,bbb,ccc ',aaa,bbb,ccc        
            llstart = ( - bbb - sqrt (rrr) ) / 2. / aaa 
            hkl (3) = llstart 
!DBGXXX        ww   = sqrt(skalpro(hkl,hkl,cr_rten))                    
!DBGXXX        write(*,'(a,5f10.4)') ' l minimum ',hkl, ww,pow_ds_max   
            llend = ( - bbb + sqrt (rrr) ) / 2. / aaa 
            hkl (3) = llend 
!DBGXXX        ww   = sqrt(skalpro(hkl,hkl,cr_rten))                    
!DBGXXX        write(*,'(a,5f10.4)') ' l maximum ',hkl, ww,pow_ds_max   
!                                                                       
            ccc = cr_rten (1, 1) * hkl (1) **2 + 2. * cr_rten (1, 2)    &
            * hkl (1) * hkl (2) + cr_rten (2, 2) * hkl (2) **2 -        &
            pow_ds_min**2                                               
            rtm = bbb**2 - 4. * aaa * ccc 
!DBGXXX          if(ih.eq.0 .and. ik.eq.0) THEN                         
!DBGXXX            write(*,*) ' llstart ',llstart                       
!DBGXXX            write(*,*) ' llend   ',llend                         
!DBGXXX            write(*,*) ' rtm ',rtm                               
!DBGXXX          endif                                                  
            IF (rtm.gt.0) THEN 
!                                                                       
!     --------- Intersection with 2Theta minimum sphere                 
!                                                                       
               llendmini = ( - bbb - sqrt (rtm) ) / 2. / aaa 
               hkl (3) = llendmini 
!DBGXXX  ww   = sqrt(skalpro(hkl,hkl,cr_rten))                          
!DBGXXX  write(*,'(a,5f10.4)')'2THETA min: l minimum ',hkl,ww,pow_ds_min
               llstartmini = ( - bbb + sqrt (rtm) ) / 2. / aaa 
               hkl (3) = llstartmini 
!DBGXXX        ww   = sqrt(skalpro(hkl,hkl,cr_rten))                    
               IF (.not.l_ano.and.ih.eq.0.and.ik.eq.0) THEN 
!                                                                       
!     ----------- Save computational time, start at 2th minimum,        
!     ----------- only one line calculated                              
!                                                                       
!DBGXXX  write(*,'(a,5f10.4)')'2THETA min: l maximum ',hkl,ww,pow_ds_min
!DBGXXX  write(*,*) ' l_ano ',l_ano                                     
!DBGXXX  write(*,*) ' ih    ',ih                                        
!DBGXXX  write(*,*) ' ik    ',ik                                        
!DBGXXX  write(*,*) ' TEST  ',(.not.l_ano .and. ih.eq.0 .and. ik.eq.0)  
!DBGXXX        write(*,*)                                               
!DBGXXX                write(* ,'(2f12.4)') hkl(2),llstart              
!DBGXXX                write( *,'(2f12.4)') hkl(2),llendmini            
!DBGXXX                write( *,'(2f12.4)') hkl(2),llstartmini          
!DBGXXX                write( *,'(2f12.4)') hkl(2),llend                
                  llstart = llstartmini 
                  l_start = int( llstart / pow_hkl_del (3) )
                  l_end = int( llend / pow_hkl_del (3) )
                  l_twoparts = .false. 
!DBG                  write(15,'(2f12.4)') hkl(2),llstart               
!DBG                  write(15,'(2f12.4)') hkl(2),llend                 
               ELSE 
!                                                                       
!     ----------- Calculate all, two lines calculated                   
!                                                                       
                  llstart2 = llstartmini 
                  llend2 = llend 
                  llstart = llstart 
                  llend = llendmini 
                  l_start = int( llstart / pow_hkl_del (3) )
                  l_end = int( llend / pow_hkl_del (3) )
                  l_twoparts = .true. 
!DBG                  write(15,'(2f12.4)') hkl(2),llstart               
!DBG                  write(15,'(2f12.4)') hkl(2),llend                 
!DBGXXX      if(ih.eq.0 .and. ik.eq.0) THEN                             
!DBGXXX                write(* ,'(2f12.4)') hkl(2),llstart2             
!DBGXXX                write( *,'(2f12.4)') hkl(2),llend2               
!DBGXXX                write( *,'(2f12.4)') hkl(2),llstart              
!DBGXXX                write( *,'(2f12.4)') hkl(2),llendmini            
!DBGXXX      endif                                                      
               ENDIF 
            ELSE 
!                                                                       
!     --------- No intersection with 2-Theta minimum                    
!                                                                       
               IF (.not.l_ano.and.ih.eq.0.and.ik.eq.0) THEN 
                  l_start = 0 
                  l_end = int( llend / pow_hkl_del (3) )
!DBG                  write(15,'(2f12.4)') hkl(2),0.00                  
!DBG                  write(15,'(2f12.4)') hkl(2),llend                 
               ELSE 
                  l_start = int( llstart / pow_hkl_del (3) )
                  l_end = int( llend / pow_hkl_del (3) )
!DBG                  write(15,'(2f12.4)') hkl(2),llstart               
!DBG                  write(15,'(2f12.4)') hkl(2),llend                 
               ENDIF 
               l_twoparts = .false. 
            ENDIF 
         ENDIF 
!DBG_RBN      END OF NEW L LIMIT                                        
!                                                                       
         IF (rrr.gt.0) THEN 
            eck (2, 1) = kk 
            eck (2, 2) = kk 
            eck (2, 3) = kk 
!                                                                       
!     --------Distinguish between intersection with tthmin and no       
!     --------New code, first segment is always run                     
!                                                                       
!                                                                       
!     ----------No intersection with 2Theta min sphere                  
!               or first of the two segments                            
!                                                                       
!DBG          if(.not.l_ano .and. ih.eq.0 .and. ik.eq.0) THEN           
!DBG                write(* ,'(3f12.4)') hkl(1),hkl(2),llstart          
!DBG                write(* ,'(3f12.4)') hkl(1),hkl(2),llendmini        
!DBG                write(* ,'(3f12.4)') hkl(1),hkl(2),llstartmini      
!DBG                write(* ,'(3f12.4)') hkl(1),hkl(2),llend            
!DBG          endif                                                     
!DBG_RBN                                                                
!DBG      l_start = -6.0/pow_hkl_del(3)                                 
!DBG      l_end   =  6.0/pow_hkl_del(3)                                 
!DBG      l_start =  llstart/pow_hkl_del(3)                             
!DBG      l_end   =  llend  /pow_hkl_del(3)                             
            eck (3, 1) = l_start * pow_hkl_del (3) + pow_hkl_shift (3) 
            eck (3, 2) = l_end * pow_hkl_del (3) + pow_hkl_shift (3) 
            eck (3, 3) = l_start * pow_hkl_del (3) + pow_hkl_shift (3) 
            inc (1) = nint ( (eck (3, 2) - eck (3, 1) ) / pow_hkl_del ( &
            3) ) + 1                                                    
!DBG_RBN                                                                
!DBGXXX          if(.not.l_ano .and. ih.eq.0 .and. ik.eq.0) THEN        
!DBGXXX      write(*,*) 'l_start , l_end ',l_start , l_end              
!DBGXXX      write(*,*) 'eck(3,*)        ',eck(3,1),eck(3,2),eck(3,3)   
!DBGXXX      write(*,*) 'inc(1)          ',inc(1)                       
!DBGXXX      endif                                                      
            IF (inc (1) * inc (2) .gt. MAXQXY  .OR.          &
                cr_nscat>DIF_MAXSCAT              ) THEN
              n_qxy   = MAX(n_qxy,inc(1) * inc(2),MAXQXY)
              n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
              n_natom = MAX(n_natom,cr_natoms,DIF_MAXAT)
              CALL alloc_diffuse (n_qxy,  n_nscat,  n_natom )
              IF (ier_num.ne.0) THEN
                RETURN
              ENDIF
            ENDIF
            IF (inc (1) .gt.MAXQXY) THEN 
               ier_num = - 8 
               ier_typ = ER_APPL 
               WRITE (ier_msg (1), 8888) inc (1) 
               WRITE (ier_msg (2), 8889) MAXQXY 
               ier_msg (3) = 'Increase dl or reduce TTHMAX' 
               RETURN 
            ENDIF 
            IF (pow_four_mode.eq.POW_FOURIER) THEN 
!              IF (pow_four_type.eq.POW_COMPL) THEN 
!DBG_RBN                                                                
!                 CALL four_run 
!              ELSE 
                  CALL four_run_powder 
!              ENDIF 
            ELSEIF (pow_four_mode.eq.POW_STACK) THEN 
               CALL st_fourier(rept_f2aver)
               rept_f2aver = .false.   ! No further calculations needed
            ENDIF 
            IF(ier_ctrlc) THEN
               ier_num = -14
               ier_typ = ER_COMM
               RETURN
            ENDIF
            IF(ier_num/=0) RETURN      ! An error occured or CTRL-C
            DO i = 1, inc (1) 
            hkl (3) = eck (3, 1) + (i - 1) * vi (3, 1) 
            ll = hkl (3) 
            l_ll_real = abs ( (nint (ll) - ll) / pow_hkl_del (3) )      &
            .gt.0.51                                                    
            IF (pow_l_all.or..not.pow_l_all.and. (                      &
            l_hh_real.or.l_kk_real.or.l_ll_real) ) THEN                 
               dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!DBG_RBN                                                                
!                 IF(pow_axis==POW_AXIS_TTH) THEN
!              IF (rlambda * 0.5 * dstar.le.1.0) THEN 
!                 ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!                 IF (pow_tthmin.le.ttheta.and.ttheta.le.(pow_tthmax+pow_deltatth))    &
!                 THEN                                                  
!                    itth = nint( (ttheta - pow_tthmin) / pow_deltatth )
!                    inten = DBLE (csf (i) * conjg (csf (i) ) ) 
!                    IF (pow_pref) THEN 
!write(*,'(a,3(f4.0,1x),1x,f5.2,1x,f10.2,1x,f10.2)') 'hkl',hkl,ttheta,   &
!                        inten , inten * calc_preferred (hkl,            &
!                        pow_pref_type, pow_pref_hkl, pow_pref_g1,       &
!                        pow_pref_g2, POW_PREF_RIET, POW_PREF_MARCH)     
!                       inten = inten * DBLE(calc_preferred (hkl,            &
!                       pow_pref_type, pow_pref_hkl, pow_pref_g1,       &
!                       pow_pref_g2, POW_PREF_RIET, POW_PREF_MARCH))     
!                    ENDIF 
!                    pow_qsp (itth) = pow_qsp (itth) + inten 
!DBG_RBN                                                                
!DBG                write(16,'(2f12.4)') hkl(2),hkl(3)                  
!                 ENDIF 
!              ENDIF 
!                 ELSEIF(pow_axis==POW_AXIS_Q  ) THEN
                     q = (zpi) * dstar
                     IF( pow_qmin <= q .AND. q <= (pow_qmax+pow_deltaq) ) THEN
                        itth = nint( (q - pow_qmin) / pow_deltaq )
                        inten = DBLE (csf (i) * conjg (csf (i) ) ) 
                        IF (pow_pref) THEN 
                           inten = inten * DBLE(calc_preferred (hkl,         &
                           pow_pref_type, pow_pref_hkl, pow_pref_g1,    &
                           pow_pref_g2, POW_PREF_RIET, POW_PREF_MARCH))  
                        ENDIF 
                        pow_qsp (itth) = pow_qsp (itth) + inten 
                     ENDIF 
!                 ENDIF 
            ENDIF 
!DBG_RBN      write(13,4444) hkl,dstar,ttheta,csf(i),                   
!DBG_RBN     &                 real(csf(i)*conjg(csf(i))),              
!DBG_RBN     &               pow_l_all,l_hh_real,l_kk_real,l_ll_real,   
!DBG_RBN     &               pow_l_all .or. .not.pow_l_all .and.        
!DBG_RBN     &               (l_hh_real .or. l_kk_real .or. l_ll_real)  
! 4444 FORMAT    (3(2x,f5.1),f8.5,f8.3,2x,f12.6,f12.6,2x,f12.6,5(2x,l1)) 
            ENDDO 
            IF (l_twoparts) THEN 
!                                                                       
!     ----------Intersection with 2Theta min sphere, calculate          
!               second section along line                               
!                                                                       
               l_start = int( llstart2 / pow_hkl_del (3) )
               l_end = int( llend2 / pow_hkl_del (3) )
!DBG                write(17,'(2f12.4)') hkl(2),llstart2                
!DBG                write(17,'(2f12.4)') hkl(2),llend2                  
!DBG_RBN                                                                
!DBG      l_start = -6.0/pow_hkl_del(3)                                 
!DBG      l_end   =  6.0/pow_hkl_del(3)                                 
!DBG      l_start =  llstart2 /pow_hkl_del(3)                           
!DBG      l_end   =  llend2   /pow_hkl_del(3)                           
               eck (3, 1) = l_start * pow_hkl_del (3) 
               eck (3, 2) = l_end * pow_hkl_del (3) 
               eck (3, 3) = l_start * pow_hkl_del (3) 
               inc (1) = nint ( (eck (3, 2) - eck (3, 1) ) /            &
               pow_hkl_del (3) ) + 1                                    
               IF (inc (1) * inc (2) .gt. MAXQXY  .OR.          &
                   cr_nscat>DIF_MAXSCAT              ) THEN
                 n_qxy   = MAX(n_qxy,inc(1) * inc(2),MAXQXY)
                 n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
                 n_natom = MAX(n_natom,cr_natoms,DIF_MAXAT)
                 CALL alloc_diffuse (n_qxy,  n_nscat,  n_natom )
                 IF (ier_num.ne.0) THEN
                   RETURN
                 ENDIF
               ENDIF
               IF (inc (1) .gt.MAXQXY) THEN 
                  ier_num = - 8 
                  ier_typ = ER_APPL 
                  WRITE (ier_msg (1), 8888) inc (1) 
                  WRITE (ier_msg (2), 8889) MAXQXY 
                  ier_msg (3) = 'Increase dl or adjust TTHMIN / TTHMAX' 
                  RETURN 
               ENDIF 
               IF (pow_four_mode.eq.POW_FOURIER) THEN 
!                 IF (pow_four_type.eq.POW_COMPL) THEN 
!DBG_RBN                                                                
!                    CALL four_run 
!                 ELSE 
                     CALL four_run_powder 
!                 ENDIF 
               ELSEIF (pow_four_mode.eq.POW_STACK) THEN 
                  CALL st_fourier(rept_f2aver)
                  rept_f2aver = .false.   ! no further calculations needed
               ENDIF 
               IF(ier_ctrlc) THEN
                  ier_num = -14
                  ier_typ = ER_COMM
                  RETURN
               ENDIF
               IF(ier_num/=0) RETURN      ! An error occured or CTRL-C
               DO i = 1, inc (1) 
               hkl (3) = eck (3, 1) + (i - 1) * vi (3, 1) 
               ll = hkl (3) 
               l_ll_real = abs ( (nint (ll) - ll) / pow_hkl_del (3) )   &
               .lt.0.51                                                 
               IF (pow_l_all.or..not.pow_l_all.and. (                   &
               l_hh_real.or.l_kk_real.or.l_ll_real) ) THEN              
                  dstar = sqrt (skalpro (hkl, hkl, cr_rten) ) 
!                 IF(pow_axis==POW_AXIS_TTH) THEN
!                 IF (rlambda * 0.5 * dstar.le.1.0) THEN 
!                    ttheta = 2.0 * asind (rlambda * 0.5 * dstar) 
!                    IF (pow_tthmin.le.ttheta.and.ttheta.le.(pow_tthmax+pow_deltatth))    &
!                    THEN                                               
!                       itth = nint( (ttheta - pow_tthmin) / pow_deltatth )
!                       inten = DBLE (csf (i) * conjg (csf (i) ) ) 
!                       IF (pow_pref) THEN 
!                          inten = inten * DBLE(calc_preferred (hkl,         &
!                          pow_pref_type, pow_pref_hkl, pow_pref_g1,    &
!                          pow_pref_g2, POW_PREF_RIET, POW_PREF_MARCH))  
!                       ENDIF 
!                       pow_qsp (itth) = pow_qsp (itth) + inten 
!DBG      write(18,'(2f12.4)') hkl(2),hkl(3)                            
!                    ENDIF 
!                 ENDIF 
!                 ELSEIF(pow_axis==POW_AXIS_Q  ) THEN
                     q = (zpi) * dstar
                     IF( pow_qmin <= q .AND. q <= (pow_qmax+pow_deltaq) ) THEN
                        itth = nint( (q - pow_qmin) / pow_deltaq )
                        inten = DBLE (csf (i) * conjg (csf (i) ) ) 
                        IF (pow_pref) THEN 
                           inten = inten * DBLE(calc_preferred (hkl,         &
                           pow_pref_type, pow_pref_hkl, pow_pref_g1,    &
                           pow_pref_g2, POW_PREF_RIET, POW_PREF_MARCH))  
                        ENDIF 
                        pow_qsp (itth) = pow_qsp (itth) + inten 
                     ENDIF 
!                 ENDIF 
               ENDIF 
!DBG_RBN      write(13,4444) hkl,dstar,ttheta,csf(i),                   
!DBG_RBN     &        real(csf(i)*conjg(csf(i))),                       
!DBG_RBN     &        pow_l_all,l_hh_real,l_kk_real,l_ll_real,          
!DBG_RBN     &        pow_l_all .or. .not.pow_l_all .and.               
!DBG_RBN     &        (l_hh_real .or. l_kk_real .or. l_ll_real)         
               ENDDO 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!
ENDDO  loop_h
!
!     CALCULATE normalized average squared atomic form factor
!
      IF(calc_f2aver) THEN
         DO i=0,pow_npkt
            pow_f2aver(i) = pow_f2aver(i) / DBLE(pow_nreal)
            pow_faver2(i) = pow_faver2(i) / DBLE(pow_nreal)
         ENDDO
         pow_faver2(:) = pow_faver2(:)**2
         pow_u2aver = pow_u2aver / pow_nreal /8./(pi**2)
      ENDIF
!                                                                       
CALL dealloc_powder_nmax ! was allocated in powder_getatoms
!
!     Prepare and calculate average form factors
!
!      natom = 0
!      DO i=1,cr_natoms
!         natom(cr_iscat(i)) = natom(cr_iscat(i)) + 1
!      ENDDO
!      pow_nreal = SUM(natom)  ! Add real atom numbers 
!      CALL powder_f2aver ( cr_nscat , natom , cr_dw)
!
!
ss = seknds (ss) 
WRITE (output_io, 4000) ss 
!
!                                                                       
 4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec') 
 5000 FORMAT     (' Currently at H = ',f9.4,'   (dH = ',f9.4,           &
     &                   ', maxH = ',f9.4,')')                          
 8888 FORMAT    ('Current number = ',i10) 
 8889 FORMAT    ('Maximum number = ',i10) 
      END SUBROUTINE powder_complete                
!
!*****7*****************************************************************
!
      SUBROUTINE powder_strucf (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(kind=PREC_DP) xarg0, xincu, twopi 
      INTEGER iscat 
      INTEGER i, ii, j, k 
      INTEGER(KIND=PREC_INT_LARGE) :: iarg, iarg0, iincu, iadd 
      LOGICAL lform 
!                                                                       
      INTEGER IAND, ISHFT 
!                                                                       
      twopi = 8.0d0 * datan (1.0d0) 
!                                                                       
!------ zero fourier array                                              
!                                                                       
      DO i = 1, num (1) * num (2) 
         tcsf (i) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
      ENDDO 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
      DO k = 1, nxat 
      xarg0 = xm (1) * xat (k, 1) + xm (2) * xat (k, 2) + xm (3)        &
      * xat (k, 3)                                                      
      xincu = uin (1) * xat (k, 1) + uin (2) * xat (k, 2) + uin (3)     &
      * xat (k, 3)                                                      
!DBG        xincv = vin(1)*xat(k,1)+vin(2)*xat(k,2)+vin(3)*xat(k,3)     
!                                                                       
      iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
      iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
!DBG        iincv = nint( 64*I2PI*( xincv-int(xincv)+1.0d0))            
      iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                                                                       
      ii = 0 
!                                                                       
      DO j = 1, num (1) 
!DBG          do i=1,num(2)                                             
      iadd = ISHFT (iarg, - 6) 
      iadd = IAND (iadd, MASK) 
      ii = ii + 1 
      tcsf (ii) = tcsf (ii) + cex (iadd) / DBLE ( (xarg0 + REAL(j - 1)&
      * xincu) * twopi)                                                 
!DBG            iarg = iarg + iincv                                     
!DBG          ENDDO                                                     
      iarg = iarg0 + iincu * j 
      ENDDO 
      ENDDO 
!                                                                       
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) THEN 
         DO i = 1, num (1) * num (2) 
         tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE powder_strucf                  
!
!*****7*****************************************************************
!
      SUBROUTINE powder_cexpt 
!+                                                                      
!     This routine initialises the complex exponent table and           
!     is called only at the first Powder run.                           
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
!                                                                       
      USE prompt_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
      REAL(kind=PREC_DP) twopi, xmult, xarg, xt 
      INTEGER i 
!                                                                       
      WRITE (output_io, 1000) 
!                                                                       
      xt = 1.0d0 / REAL (I2PI, KIND=KIND(0.0D0)) 
      twopi = 8.0d0 * datan (1.0d0) 
!                                                                       
!DBG      open(9,file='CEX.DAT',status='unknown')                       
      DO i = 0, MASK 
      xmult = REAL(i, kind=PREC_DP) * xt 
      xarg = twopi * xmult 
      cex (i) = cmplx (int( sin (xarg) ), 0.0) 
!DBG      write(9,*) xarg,real(cex(i))                                  
      ENDDO 
      ffour = .false. 
!DBG      close(9)                                                      
!                                                                       
 1000 FORMAT     (' Computing complex exponent table ...') 
      END SUBROUTINE powder_cexpt                   
!*****7*****************************************************************
!!      SUBROUTINE powder_sine_f (iscat, jscat) 
!+                                                                      
!     Here the real structure factor of 'nxat' identical atoms          
!     from array 'xat' is computed.                                     
!-                                                                      
!!      USE discus_config_mod 
!!      USE debye_mod 
!!      USE diffuse_mod 
!!      IMPLICIT none 
!                                                                       
!                                                                       
!!      REAL(PREC_DP) xarg0, xincu, twopi 
!!      INTEGER iscat, jscat 
!!      INTEGER i, ii, j, k, iarg, iarg0, iincu, iadd 
!                                                                       
!!      INTEGER IAND, ISHFT 
!                                                                       
!!      twopi = 8.0d0 * datan (1.0d0) 
!                                                                       
!------ zero fourier array                                              
!                                                                       
!!      DO i = 1, num (1) * num (2) 
!!      tcsf (i) = cmplx (0.0, 0.0) 
!!      ENDDO 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
!!      DO k = 1, nxat 
!!      xarg0 = xm (1) * xat (k, 1) + xm (2) * xat (k, 2) + xm (3)        &
!!      * xat (k, 3)                                                      
!!      xincu = uin (1) * xat (k, 1) + uin (2) * xat (k, 2) + uin (3)     &
!!      * xat (k, 3)                                                      
!DBG        xincv = vin(1)*xat(k,1)+vin(2)*xat(k,2)+vin(3)*xat(k,3)     
!                                                                       
!!      iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
!!      iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
!DBG        iincv = nint( 64*I2PI*( xincv-int(xincv)+1.0d0))            
!!      iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                                                                       
!!      ii = 0 
!                                                                       
!!      DO j = 1, num (1) 
!DBG          do i=1,num(2)                                             
!!      iadd = ISHFT (iarg, - 6) 
!!      iadd = IAND (iadd, MASK) 
!!      ii = ii + 1 
!!      partial (ii, look (iscat, jscat),0 ) = partial (ii, look (iscat,    &
!!      jscat),0 ) + sinetab (iadd) / real ( (xarg0 + REAL(j - 1) * xincu)&
!!      * twopi)                                                          
!DBG            iarg = iarg + iincv                                     
!DBG          ENDDO                                                     
!!      iarg = iarg0 + iincu * j 
!!      ENDDO 
!!      ENDDO 
!                                                                       
!!      END SUBROUTINE powder_sine_f                  
!*****7*****************************************************************
      SUBROUTINE powder_getatm (iscat, i_start) 
!+                                                                      
!     This routine creates an atom list of atoms of type 'iscat'        
!     which are within the current lot.                                 
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE metric_mod
use precision_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER iscat, i_start 
      INTEGER i, j 
      REAL(kind=PREC_DP), dimension(3) :: u (3), v (3) 
!                                                                       
!     REAL do_blen 
!                                                                       
      nxat = 0 
!                                                                       
      DO i = i_start + 1, cr_natoms 
      IF (cr_iscat (i) .eq.iscat) THEN 
         nxat = nxat + 1 
         DO j = 1, 3 
         u (j) = cr_pos (j, i_start) 
         v (j) = cr_pos (j, i) 
         xat (nxat, j) = 0.0 
         ENDDO 
         xat (nxat, 1) = do_blen (.true., u, v) 
      ENDIF 
      ENDDO 
      END SUBROUTINE powder_getatm                  
!*****7*****************************************************************
      SUBROUTINE four_run_powder 
!+                                                                      
!     claculates the Fourier, complete mode                             
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE four_strucf_mod
      USE fourier_sup
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER lbeg (3), csize (3) 
      INTEGER iscat, i 
!                                                                       
      ier_num = 0 
      csize (1) = cr_icc (1) 
      csize (2) = cr_icc (2) 
      csize (3) = cr_icc (3) 
!                                                                       
!------ preset some values                                              
!                                                                       
      CALL four_layer 
!                                                                       
!------ zero some arrays                                                
!                                                                       
      DO i = 1, num (1) * num (2) 
         csf (i) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
!DBG        acsf(i) = cmplx(0.0d0,0.0d0)                                
!DBG         dsi(i) = 0.0d0                                             
      ENDDO 
!                                                                       
!------ preset some tables, calculate average structure                 
!                                                                       
      CALL four_stltab 
      IF (ier_num.ne.0) return 
      lbeg (1) = 1 
      lbeg (2) = 1 
      lbeg (3) = 1 
!                                                                       
!------ - loop over all different atom types                            
!                                                                       
      DO iscat = 1, cr_nscat 
!DBG        call four_getatm (iscat,ilots,lbeg,csize,ncell)             
!DBG        call four_strucf (iscat,.true.)                             
      CALL powder_strucfactor (iscat, .true.) 
!                                                                       
!------ --- Add this part of the structur factor to the total           
!                                                                       
      DO i = 1, num (1) * num (2) 
      csf (i) = csf (i) + tcsf (i) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE four_run_powder                
!*****7*****************************************************************
      SUBROUTINE powder_getatoms 
!+                                                                      
!     sorts all atoms by scattering type                                
!                                                                       
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE powder_scat_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j 
      INTEGER iscat 
!
      CALL alloc_powder_nmax(MAXSCAT, cr_natoms) ! will automatically be deallocated 
!
      DO j = 1, cr_nscat 
         pow_nscat (j) = 0 
      ENDDO 
      DO i = 1, cr_natoms 
         iscat = cr_iscat (i) 
         IF (iscat.gt.0) THEN 
            pow_nscat (iscat) = pow_nscat (iscat) + 1 
            pow_iatom (iscat, pow_nscat (iscat) ) = i 
         ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE powder_getatoms                
!*****7*****************************************************************
      SUBROUTINE powder_strucfactor (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE powder_scat_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP) xarg0, xincu, xincv 
      INTEGER iscat 
      INTEGER i, ii, j, k
      INTEGER(KIND=PREC_INT_LARGE) :: iarg, iarg0, iincu, iincv, iadd 
      LOGICAL lform 
!                                                                       
      INTEGER IAND, ISHFT 
!                                                                       
!------ zero fourier array                                              
!                                                                       
!DBGXXX      do i=1,num(1)*num(2)                                       
      DO i = 1, num (1) 
         tcsf (i) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
      ENDDO 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
      DO k = 1, pow_nscat (iscat) 
      xarg0 = xm (1) * cr_pos (1, pow_iatom (iscat, k) ) + xm (2)       &
      * cr_pos (2, pow_iatom (iscat, k) ) + xm (3) * cr_pos (3,         &
      pow_iatom (iscat, k) )                                            
      xincu = uin (1) * cr_pos (1, pow_iatom (iscat, k) ) + uin (2)     &
      * cr_pos (2, pow_iatom (iscat, k) ) + uin (3) * cr_pos (3,        &
      pow_iatom (iscat, k) )                                            
      xincv = vin (1) * cr_pos (1, pow_iatom (iscat, k) ) + vin (2)     &
      * cr_pos (2, pow_iatom (iscat, k) ) + vin (3) * cr_pos (3,        &
      pow_iatom (iscat, k) )                                            
!                                                                       
      iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
      iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
      iincv = nint (64 * I2PI * (xincv - int (xincv) + 1.0d0) ) 
      iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                                                                       
      ii = 0 
!                                                                       
      DO j = 1, num (1) 
!DBGXXX          do i=1,num(2)                                          
      iadd = ISHFT (iarg, - 6) 
      iadd = IAND (iadd, MASK) 
      ii = ii + 1 
      tcsf (ii) = tcsf (ii) + cex (iadd) 
      iarg = iarg + iincv 
!DBGXXX          ENDDO                                                  
      iarg = iarg0 + iincu * j 
      ENDDO 
      ENDDO 
!                                                                       
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) THEN 
!DBGXXX        do i=1,num(1)*num(2)                                     
         DO i = 1, num (1) 
         tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE powder_strucfactor             
!*****7*****************************************************************
      REAL(kind=PREC_DP) function calc_preferred (w, pow_pref_type, pow_pref_hkl,     &
      pow_pref_g1, pow_pref_g2, POW_PREF_RIET, POW_PREF_MARCH)          
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
use precision_mod
USE metric_mod
      USE trig_degree_mod
      IMPLICIT none 
!                                                                       
      REAL(kind=PREC_DP), dimension(3) ::  w (3) 
      INTEGER pow_pref_type 
      REAL(kind=PREC_DP), dimension(3) :: pow_pref_hkl (3) 
      REAL(kind=PREC_DP) :: pow_pref_g1 
      REAL(kind=PREC_DP) :: pow_pref_g2 
      INTEGER POW_PREF_RIET 
      INTEGER POW_PREF_MARCH 
!                                                                       
      LOGICAL lspace 
      REAL(kind=PREC_DP), dimension(3), parameter :: nullv = (/ 0.0D0, 0.0D0, 0.0D0 /)
      REAL(KIND=PREC_DP) :: alpha   = 0.0
      REAL(KIND=PREC_DP) :: alpha2  = 0.0
!                                                                       
      lspace = .true. 
!                                                                       
      calc_preferred = 1.0 
!                                                                       
      alpha = do_bang (lspace, w, NULLV, pow_pref_hkl) 
      IF (pow_pref_type.eq.POW_PREF_RIET) THEN 
         IF (alpha.le.90.) THEN 
            alpha2 = alpha**2 
         ELSE 
            alpha2 = (180 - alpha) **2 
         ENDIF 
         calc_preferred = pow_pref_g2 + (1. - pow_pref_g2) * exp (      &
         pow_pref_g1 * alpha2)                                          
      ELSEIF (pow_pref_type.eq.POW_PREF_MARCH) THEN 
         calc_preferred = pow_pref_g2 + (1. - pow_pref_g2) * ( (        &
         pow_pref_g1 * cosd (alpha) ) **2 + (sind (alpha) ) **2 /       &
         pow_pref_g1) ** ( - 1.5)                                       
      ENDIF 
!     write(*,'(3f4.0,2x,f10.2,2x,f10.2,2x,f8.5)') w,alpha,alpha2,  &
!                calc_preferred                                     
      END FUNCTION calc_preferred                   
!*****7*****************************************************************
      SUBROUTINE proj_preferred (w, pow_pref_hkl) 
!+                                                                      
!                                                                       
!-                                                                      
      USE metric_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE param_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL(kind=PREC_DP), DIMENSION(3) :: w (3) 
      REAL(kind=PREC_DP), dimension(3) :: pow_pref_hkl (3) 
!                                                                       
      INTEGER i 
      LOGICAL lspace 
      REAL(KIND=PREC_DP) :: uv, uu, vv 
!                                                                       
      lspace = .true. 
!                                                                       
!                                                                       
!     ------Calculate projection onto second vector, always in          
!             direct space                                              
!                                                                       
      uv = skalpro (w, pow_pref_hkl, cr_gten) 
      uu = skalpro (w, pow_pref_hkl, cr_gten) 
      vv = skalpro (w, pow_pref_hkl, cr_gten) 
      IF (vv.eq.0.0) THEN 
         w (1) = 0.0 
         w (2) = 0.0 
         w (3) = 0.0 
         RETURN 
      ENDIF 
      DO i = 1, 3 
      w (i) = pow_pref_hkl (i) * uv / vv 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE proj_preferred                 
!
!*******************************************************************************
!
SUBROUTINE pow_conv_limits
!-
! convert limits and steps
!+
USE powder_mod
USE diffuse_mod
!
USE trig_degree_mod
USE wink_mod
!
IMPLICIT NONE
!
IF(pow_qtthmin) THEN
   pow_tthmin = 2.0D0*asind(pow_qmin*rlambda/fpi)
ELSE
   pow_qmin = fpi*sind(pow_tthmin*0.5D0)/rlambda
ENDIF
!
IF(pow_qtthmax) THEN
   pow_tthmax = 2.0D0*asind(pow_qmax*rlambda/fpi)
ELSE
   pow_qmax = fpi*sind(pow_tthmax*0.5D0)/rlambda
ENDIF
!
IF(pow_deltaqtth) THEN
   pow_deltatth = 2.0D0*asind((pow_qmin+pow_deltaq)*rlambda/fpi) - &
                  2.0D0*asind((pow_qmin           )*rlambda/fpi)
ELSE
   pow_deltaq   = fpi*sind(0.5D0*(pow_tthmax             ))/rlambda     - &
                  fpi*sind(0.5D0*(pow_tthmax-pow_deltatth))/rlambda
ENDIF
!
IF(pow_qtthzero) THEN
  pow_tthzero = 2.0D0*asind(pow_qzero*rlambda/fpi)
ELSE
   pow_qzero  = fpi*sind(pow_tthzero*0.5D0)/rlambda
ENDIF
!
END SUBROUTINE pow_conv_limits                 
!
!*******************************************************************************
!
SUBROUTINE powder_reset

USE discus_allocate_appl_mod
USE powder_mod
USE powder_scat_mod
!
CALL alloc_powder(1, 1)
CALL alloc_powder_nmax(1, 1)
!
pow_axis       = POW_AXIS_Q
pow_npkt       = 1           ! Actual number of powder data points
!
pow_four_mode  = 0
pow_four_type  = POW_COMPL
!
pow_lp         = POW_LP_BRAGG
!
pow_l_all      = .true.
!
pow_qtthmin    = .TRUE.    ! User  provided Qmin(==true) TTHmin(=false)
pow_qtthmax    = .TRUE.    ! User  provided Qmax(==true) TTHmax(=false)
pow_deltaqtth  = .TRUE.    ! User  provided Qstp(==true) TTHstp(=false)
pow_qtthzero   = .TRUE.    ! User  provided Qzero(==true) TTHzero(=false)

pow_tthzero    =  0.0
pow_tthmin     =  0.1
pow_tthmax     = 40.0
pow_deltatth   =  0.05
pow_qzero      =  0.0
pow_qmin       =  0.2
pow_qmax       =  7.0
pow_deltaq     =  0.001
pow_ds_max     =  0.001
pow_ds_min     =  0.001
pow_delta      =  0.0
pow_lp_fac     =  0.88
pow_lp_ang     = 20.0
pow_lp_cos     =  0.936
!pow_bvalue     =  0.000
!
pow_nback      = 0
pow_back(:)    = 0.0
pow_scale      = 1.0
!
pow_hkl_max(:) = 4.0
pow_hkl_del(:) = 0.05
pow_hkl_shift(:)= 0.00
!
pow_pref       = .false.
pow_pref_type  = POW_PREF_RIET
pow_pref_g1    = 0.0
pow_pref_g2    = 0.0
pow_pref_hkl(:)   = (/0., 0., 1./)
!
pow_profile    = POW_PROFILE_PSVGT
pow_pr_par     =  0
pow_fwhm       =  0.01
pow_eta        =  0.5
pow_eta_l      =  0.0
pow_eta_q      =  0.0
pow_u          =  0.0
pow_v          =  0.0
pow_w          =  0.05
pow_asym       =  0.0
pow_width      = 20.0
!
pow_ka21       =  0.0
pow_ka21_u     =  .FALSE.
!
IF(ALLOCATED(pow_qsp))    pow_qsp       = 0.0D0  !  (0:POW_MAXPKT)
IF(ALLOCATED(pow_f2aver)) pow_f2aver    = 0.0D0  !  (0:POW_MAXPKT)
IF(ALLOCATED(pow_faver2)) pow_faver2    = 0.0D0  !  (0:POW_MAXPKT)
IF(ALLOCATED(pow_f2    )) pow_f2        = 0.0D0  !  (0:POW_MAXPKT)
IF(ALLOCATED(pow_conv  )) pow_conv      = 0.0D0  !  (0:POW_MAXPKT)
IF(ALLOCATED(pow_sq    )) pow_sq        = 0.0D0  !  (0:POW_MAXPKT)
pow_nreal   = 0
pow_ncreal  = 0
pow_u2aver  = 0.0
!
! powder_scat
!
POW_NMAX  = 1
POW_MAXSCAT = 1
IF(ALLOCATED(pow_nscat)) pow_nscat(:)   = 0
IF(ALLOCATED(pow_iatom)) pow_iatom(:,:) = 0
!
!
END SUBROUTINE powder_reset
!
!*******************************************************************************
!
END MODULE powder
