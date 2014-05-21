MODULE discus_init_mod
!
CONTAINS
!
      SUBROUTINE discus_initarrays 
!                                                                       
!+                                                                      
!       Startvalues for all important arrays.                           
!-                                                                      
      USE crystal_mod 
      USE config_mod 
      USE atom_env_mod 
      USE chem_mod 
      USE debye_mod 
      USE diffuse_mod 
      USE domaindis_mod 
      USE gen_add_mod 
      USE mc_mod 
      USE molecule_mod 
      USE output_mod 
      USE patters_mod 
      USE pdf_mod 
      USE plot_mod 
      USE powder_mod 
      USE prop_para_mod 
      USE rmc_mod 
      USE save_mod 
      USE shear_mod 
      USE spcgr_setup_mod
      USE stack_mod 
      USE surface_mod 
      USE sym_add_mod 
      USE symm_mod 
      USE transfrm_mod 
      USE unitcell_mod 
      USE waves_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j, k
      INTEGER stift (0:4) 
!                                                                       
      DATA stift / 5, 3, 1, 2, 4 / 
!                                                                       
!     Sockt defaults                                                    
!                                                                       
      s_port = 3330 
      s_ipallowed = "localhost" 
!                                                                       
!     /errors/                                                   
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      ier_sta = ER_S_CONT 
      ier_msg (:) = ' ' 
!
!     /space_group/
!
      CALL spcgr_setup 
!                                                                       
!     /chem/                                                
!                                                                       
      DO i = 1, CHEM_MAX_COR 
      chem_ldall (i) = .false. 
      chem_cang (i) = .false. 
      chem_ctyp (i) = CHEM_NONE 
      chem_nvec (i) = 0 
      chem_freq_sigma (i) = 0.05 
      chem_wink_sigma (i) = 0.05 
      chem_dir (1, 1, i) = - 9999 
      ENDDO 
!                                                                       
!     /mcbl/                                                
!                                                                       
      mo_energy = MC_NONE 
      mo_mode = RMC_MODE_SWCHEM 
      mo_local = RMC_LOCAL_ALL 
      mo_cyc = 100 
      mo_feed = 10 
      mo_kt = 1.0 
      mo_atom (1) = 'ALL ' 
      mo_atom (2) = 'ALL ' 
      mo_atom (3)  = '    ' 
      mo_sel_atom = .true. 
!                                                                       
      DO k = 1, CHEM_MAX_COR 
      DO i = 0, MAXSCAT 
      DO j = 0, MAXSCAT 
      mo_disp (k, i, j) = 0.0 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      mo_maxmove (1, 0) = 0.0 
      mo_maxmove (2, 0) = 0.0 
      mo_maxmove (3, 0) = 0.0 
      DO i = 1, MAXSCAT 
      mo_maxmove (1, i) = 0.2 
      mo_maxmove (2, i) = 0.2 
      mo_maxmove (3, i) = 0.2 
      ENDDO 
!                                                                       
      DO i = 1, CHEM_MAX_COR 
      mo_const (i) = 0.0 
      mo_cfac (i) = 1.0 
      ENDDO 
      mo_const (0) = 1.0 
!                                                                       
!     /mmcbl/                                                   
!                                                                       
      CALL mmc_init 
!                                                                       
!     /rmbl/                                                
!                                                                       
      rmc_display = 500 
      rmc_maxcyc = 5000 
      rmc_nplane = 0 
      rmc_calc_f = .TRUE. 
      rmc_doskal = .FALSE. 
      rmc_doback = .FALSE. 
      rmc_dosym = .FALSE. 
      rmc_nosym = .FALSE. 
      rmc_log = .FALSE. 
      rmc_ranloc = .TRUE. 
      rmc_sel_atom = .TRUE. 
      rmc_sigma = 0.01 
      rmc_ave = 0.0 
      rmc_ilots = LOT_OFF 
      rmc_nlots = 1 
      rmc_qmin = 9999.0 
      rmc_qmax = - 9999.0 
      rmc_data = RMC_DATA_NIPL 
      rmc_mode = RMC_MODE_SHIFT 
      rmc_local = RMC_LOCAL_ALL 
      rmc_sel_prop (0) = 0 
      rmc_sel_prop (1) = 0 
!                                                                       
!     rmc_allowed (0) = .TRUE. 
      rmc_maxmove (1, 0) = 0.0 
      rmc_maxmove (2, 0) = 0.0 
      rmc_maxmove (3, 0) = 0.0 
!     DO i = 1, RMC_MAXSCAT 
!        rmc_repl    (i) = 0
!        rmc_allowed (i) = .TRUE. 
!     ENDDO
!     DO i = 1, MAXSCAT 
!     rmc_maxmove (1, i) = 0.2 
!     rmc_maxmove (2, i) = 0.2 
!     rmc_maxmove (3, i) = 0.2 
!     DO j = 1, MAXSCAT 
!     rmc_mindist (i, j) = 0.5 
!     ENDDO 
!     ENDDO 
!                                                                       
!     DO i = 1, RMC_MAX_PLANES 
!     rmc_fname (i) = ' ' 
!     rmc_num (1, i) = 0 
!     rmc_num (2, i) = 0 
!     rmc_ldbw (i) = .FALSE. 
!     rmc_ano (i) = .FALSE. 
!     rmc_skal (i) = 1.0 
!     rmc_back (i) = 0.0 
!     rmc_wic_typ (i) = RMC_WIC_EINS 
!     rmc_constrain (i) = i 
!     offq (i) = 0 
!     DO j = 1, RMC_MAX_SYM 
!     offsq (i, j) = 0 
!     ENDDO 
!     ENDDO 
!
      END SUBROUTINE discus_initarrays                     
!*****7*****************************************************************
      SUBROUTINE mmc_init 
!+                                                                      
!     initializes all mmc variables                                     
!-                                                                      
      USE config_mod 
      USE mmc_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
!      INTEGER i, j, k, l, m 
!                                                                       
!     /mmcbl/                                                   
!                                                                       
      mmc_move = 0 
!     DO i = 1, MC_N_MOVE 
      mmc_local       = rmc_local_all 
      mmc_move_prob   = 0.0 
      mmc_move_cprob  = 0.0 
!     ENDDO 
!     DO i = 0, CHEM_MAX_COR 
!     DO j = 0, MC_N_ENERGY 
      mmc_cor_energy        = .false. 
      mmc_const             = 1.0 
      mmc_cfac              = 1.0 
!     ENDDO 
!     ENDDO 
!     DO i = 1, CHEM_MAX_COR * MMC_MAX_ANGLES 
!     mmc_angles (i) = 0 
!     mmc_target_angl (i) = 0.0 
!     mmc_ang_sigm (i) = 0.0 
!     ENDDO 
      mmc_n_angles = 0 
!     DO m = 0, MAXSCAT 
      mmc_allowed = .false. 
!     DO l = 0, MAXSCAT 
!     DO k = 1, CHEM_MAX_COR 
!     DO j = 1, 12 
!     DO i = 1, 4 
!     mmc_vec (i, j, k, l, m) = 0.0 
!     ENDDO 
!     ENDDO 
!     mmc_nvec (k, l, m) = 0 
!     DO j = 0, MC_N_ENERGY 
      mmc_target_corr = 0.0 
      mmc_ach_corr    = 0.0 
      mmc_pair        = 0
!     ENDDO 
!     ENDDO 
!     ENDDO 
!     ENDDO 
!     DO i = 1, MC_N_MOVE 
      mmc_local      = rmc_local_all 
      mmc_move_prob  = 0.0 
      mmc_move_cprob = 0.0 
!     ENDDO 
!                                                                       
      mmc_l_constrains = .false. 
!     DO i = 1, 3 
      mmc_c_min  = -1.e10 
      mmc_c_max  =  1.e10 
!     ENDDO 
!                                                                       
!     /mmcav/                                                   
!                                                                       
!     DO i = 0, MC_N_ENERGY 
      n_e_av_p     = 0 
      n_e_av_m     = 0 
      n_e_av_z     = 0 
      e_aver_p     = 0.0 
      e_aver_m     = 0.0 
!     ENDDO 
!                                                                       
      mmc_sel_prop = 0 
!     mmc_sel_prop = 0 
!                                                                       
      END SUBROUTINE mmc_init                       
!*****7*****************************************************************
      SUBROUTINE discus_autodef 
!-                                                                      
!     Tries to open a default file for the integer and real variables   
!+                                                                      
      USE envir_mod 
      USE param_mod 
      USE errlist_mod 
      USE prompt_mod 
!
      IMPLICIT NONE
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
      END SUBROUTINE discus_autodef                        
END MODULE discus_init_mod
