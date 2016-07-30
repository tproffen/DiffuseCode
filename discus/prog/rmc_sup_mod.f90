MODULE rmc_sup_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE rmc_set (zeile, lp) 
!+                                                                      
!     sets most parameters for 'set' section                            
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE modify_mod
      USE rmc_mod 
      USE quad_mod
      USE errlist_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER, PARAMETER :: MIN_PARA = 23  ! A command requires at leaset these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER ( LEN=* ), INTENT(INOUT) :: zeile 
      INTEGER            , INTENT(INOUT) :: lp 
!                                                                       
!     CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
!     REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
!     REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: wwerte
!     INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
      CHARACTER(LEN=1024), DIMENSION(25                     ) :: cpara
      REAL               , DIMENSION(25                     ) :: werte
      REAL               , DIMENSION(25                     ) :: wwerte
      INTEGER            , DIMENSION(25                     ) :: lpara
!
      REAL mmdis, hklmin (3), hklmax (3)
!, quad 
      INTEGER ianz, iianz, jjanz, is, js, i, j, ii, ie1, ie2, ip 
      LOGICAL flag 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      maxw = 25
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.2) then 
            CALL do_cap (cpara (1) ) 
!                                                                       
!------ --- 'set aver': sets the portion of the crystal used to calc <F>
!                                                                       
            IF (cpara (1) (1:2) .eq.'AV') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
                     IF (werte (1) .ge.0.0.and.werte (1) .le.100.0)     &
                     then                                               
                        rmc_ave = werte (1) * 0.01 
                        rmc_calc_f = .true. 
                     ELSE 
                        ier_num = - 1 
                        ier_typ = ER_FOUR 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- 'set back': background correction for all planes on/off     
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'BA') then 
               CALL do_cap (cpara (2) ) 
               rmc_doback = (cpara (2) (1:2) .eq.'ON') 
               IF (ianz.eq.3.or.ianz.eq.4) then 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  IF (ier_num.eq.0) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ianz.eq.2) then 
                        ip = nint (werte (2) ) 
                        IF (ip.gt.0.and.ip.le.rmc_max_planes) then 
                           rmc_back (ip) = werte (1) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSE 
                        DO i = 1, rmc_max_planes 
                        rmc_back (i) = werte (1) 
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  DO i = 1, rmc_max_planes 
                  rmc_back (i) = 0.0 
                  ENDDO 
               ENDIF 
!                                                                       
!------ --- 'set constrain': setting possible skaling/background        
!                            constrains                                 
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'CO') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     ie1 = nint (werte (1) ) 
                     ie2 = nint (werte (2) ) 
                     IF (ie1.gt.0.and.ie1.le.rmc_max_planes.and. &
                         ie2.gt.0.and.ie2.le.rmc_max_planes.and. &
                         ie2.le.ie1                        ) then                                  
                        rmc_constrain (ie1) = ie2 
                     ELSE 
                        ier_num = - 19 
                        ier_typ = ER_RMC 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set cycl': setting of number of RMC cycles                 
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'CY') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.ne.0) return 
                     rmc_maxcyc = int (werte (1) ) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set data' : setting data type for input data               
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'DA') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.ge.1) then 
                     CALL do_cap (cpara (1) ) 
                     IF (cpara (1) (1:1) .eq.'N') then 
                        rmc_data = rmc_data_nipl 
                     ELSEIF (cpara (1) (1:1) .eq.'P') then 
                        rmc_data = rmc_data_pgm 
                     ELSE 
                        ier_num = - 15 
                        ier_typ = ER_RMC 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set dbw'  : debye-waller factor control flag setting       
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'DB') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.ge.1) then 
                     CALL do_cap (cpara (1) ) 
                     flag = (cpara (1) (1:2) .eq.'ON') 
                     IF (ianz.gt.1) then 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        DO i = 1, ianz 
                        ii = nint (werte (i) ) 
                        IF (ii.gt.0.and.ii.le.rmc_nplane) then 
                           rmc_ldbw (ii) = flag 
                        ELSE 
                           ier_num = - 7 
                           ier_typ = ER_RMC 
                        ENDIF 
                        ENDDO 
                     ELSE 
                        DO i = 1, rmc_nplane 
                        rmc_ldbw (i) = flag 
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set disp' : setting of output intervall                    
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'DI') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     rmc_display = max (1, int (werte (1) ) ) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set log': Output RMC progress in rmc.log - on/off          
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'LOG') then 
               CALL do_cap (cpara (2) ) 
               rmc_log = (cpara (2) (1:2) .eq.'ON') 
!                                                                       
!     --- 'set lots' : Setting of RMC lots                              
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'LOT') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL do_cap (cpara (1) ) 
                     IF (cpara (1) (1:1) .eq.'O') then 
                        rmc_ilots = LOT_OFF 
                        rmc_nlots = 1 
                        rmc_calc_f = .true. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSEIF (ianz.ge.6) then 
                     CALL do_cap (cpara (1) ) 
                     IF (cpara (1) (1:1) .eq.'B') then 
                        rmc_ilots = LOT_BOX 
                        rmc_calc_f = .true. 
                     ELSEIF (cpara (1) (1:1) .eq.'E') then 
                        rmc_ilots = LOT_ELI 
                        rmc_calc_f = .true. 
                     ELSE 
                        ier_num = - 2 
                        ier_typ = ER_FOUR 
                     ENDIF 
                     IF (ier_num.eq.0) then 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (4, cpara, lpara, werte, maxw) 
                        IF (ier_num.eq.0) then 
!                          IF (nint (werte (4) ) .le.RMC_MAX_LOTS) then 
                              ls_xyz (1) = nint (werte (1) ) 
                              ls_xyz (2) = nint (werte (2) ) 
                              ls_xyz (3) = nint (werte (3) ) 
                              rmc_nlots = nint (werte (4) ) 
                              CALL do_cap (cpara (5) ) 
                              lperiod = (cpara (5) (1:1) .eq.'Y') 
                              IF (ianz.gt.5) then 
                                 CALL del_params (5, ianz, cpara, lpara,&
                                 maxw)                                  
                                 CALL do_build_name (ianz, cpara, lpara,&
                                 werte, maxw, 1)                        
                                 IF (ier_num.eq.0) then 
                                    rmc_ranloc = .FALSE. 
                                    rmc_lname = cpara (1) (1:lpara(1))
                                 ENDIF 
                              ELSE 
                                 rmc_ranloc = .TRUE. 
                              ENDIF 
!                          ELSE 
!                             ier_num = - 21 
!                             ier_typ = ER_RMC 
!                          ENDIF 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set mdis' : sets minimal allowed atom distances            
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'MD') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     iianz = 1 
                     CALL get_iscat (iianz, cpara (1), lpara (1),       &
                     werte, maxw, .false.)                              
                     IF (ier_num.eq.0) then 
                        jjanz = 1 
                        CALL get_iscat (jjanz, cpara (2), lpara (2),    &
                        wwerte, maxw, .false.)                          
                        IF (ier_num.eq.0) then 
                           CALL del_params (2, ianz, cpara, lpara, maxw) 
                           CALL ber_params (1, cpara, lpara, mmdis, 1) 
                           IF (ier_num.eq.0) then 
                              DO i = 1, iianz 
                              DO j = 1, jjanz 
                              is = int (werte (i) ) 
                              js = int (wwerte (j) ) 
                              CALL rmc_set_mindist (is, js, mmdis) 
                              ENDDO 
                              ENDDO 
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set mode': sets operation mode for RMC                     
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'MOD') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               j = 1 
               CALL ber_params (j, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL rmc_set_mode (                     ianz, cpara,lpara, werte,maxw)
!              CALL rmc_set_mode (rmc_mode, rmc_local, ianz, cpara,lpara, werte,maxw)
!                                                                       
!------ --- 'set move': sets factor for generated moves                 
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'MOV') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL rmc_set_move (rmc_maxmove, rmc_sel_atom, ianz,      &
               cpara, werte, lpara, maxw, 3)                               
!                                                                       
!------ --- 'set range': set q-range to be used in RMC runs             
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'RAN') then 
               IF (ianz.eq.2) then 
                  CALL do_cap (cpara (2) ) 
                  IF (cpara (2) (1:3) .eq.'ALL') then 
                     rmc_llim = rmc_qmin 
                     rmc_ulim = rmc_qmax 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSEIF (ianz.eq.3) then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  IF (ier_num.ne.0) return 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  rmc_llim = werte (1) / 2.0 
                  rmc_ulim = werte (2) / 2.0 
               ELSEIF (ianz.eq.7) then 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  IF (ier_num.ne.0) return 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  hklmin (1) = werte (1) 
                  hklmin (2) = werte (2) 
                  hklmin (3) = werte (3) 
                  hklmax (1) = werte (4) 
                  hklmax (2) = werte (5) 
                  hklmax (3) = werte (6) 
                  rmc_llim = sqrt (quad (hklmin, hklmin, cr_rten)       &
                  / 4.0)                                                
                  rmc_ulim = sqrt (quad (hklmax, hklmax, cr_rten)       &
                  / 4.0)                                                
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ --- 'set scal': scalingfactor for all planes on/off             
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'SCA') then 
               CALL do_cap (cpara (2) ) 
               rmc_doskal = (cpara (2) (1:2) .eq.'ON') 
               IF (ianz.eq.3.or.ianz.eq.4) then 
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  IF (ier_num.eq.0) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ianz.eq.2) then 
                        ip = nint (werte (2) ) 
                        IF (ip.gt.0.and.ip.le.rmc_max_planes) then 
                           rmc_skal (ip) = 1.0 / werte (1) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSE 
                        DO i = 1, rmc_max_planes 
                        rmc_skal (i) = 1.0 / werte (1) 
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  DO i = 1, rmc_max_planes 
                  rmc_skal (i) = 1.0 
                  ENDDO 
               ENDIF 
!                                                                       
!------ --- 'set sigma': sets SIGMA for CHI2 calculation                
!                                                                       
            ELSEIF (cpara (1) (1:2) .eq.'SI') then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     rmc_sigma = werte (1) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ --- 'set sym': include symmetry of measured data on/off         
!                                                                       
            ELSEIF (cpara (1) (1:3) .eq.'SYM') then 
               CALL do_cap (cpara (2) ) 
               rmc_dosym = (cpara (2) (1:2) .eq.'ON') 
               rmc_calc_f = .true. 
!                                                                       
               IF (rmc_dosym.and.rmc_nosym) then 
                  rmc_dosym = .false. 
                  ier_num = - 13 
                  ier_typ = ER_RMC 
               ENDIF 
!                                                                       
!------ --- unknown subcommand entered                                  
!                                                                       
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_set                        
!*****7*****************************************************************
      SUBROUTINE rmc_set_move (move, sel_atom, ianz, cpara, werte,      &
      lpara, maxw,dimen)                                                      
!+                                                                      
!     Sets RMC/MC maximum shift in 'shift' mode                         
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: maxww = 4 
!
      INTEGER                             , INTENT(IN)    :: dimen
      INTEGER                             , INTENT(IN )   :: maxw
      REAL    , DIMENSION(dimen,0:MAXSCAT), INTENT(OUT)   :: move ! (3, 0:MAXSCAT) 
      LOGICAL                             , INTENT(IN )   :: sel_atom 
      INTEGER                             , INTENT(INOUT) :: ianz
!                                                                       
      CHARACTER (LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara !(maxw) 
      REAL             , DIMENSION(1:MAXW), INTENT(INOUT) :: werte !(maxw) 
      INTEGER          , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara !(maxw) 
!
      REAL wwerte (maxww) 
      INTEGER ii, is, i, j 
!                                                                       
!------ Atoms                                                           
!                                                                       
      wwerte(:) = 0.0
      IF (sel_atom) then 
         IF ((dimen==3 .AND.ianz.eq.4) .OR.           &
             (dimen==4 .AND. (ianz==4 .OR. ianz==5))) then 
            IF(dimen==3) THEN
               ii = ianz - 3
            ELSEIF(dimen==4) THEN
               ii = 1       ! MMC allows only ONE atom name
            ENDIF
            CALL get_iscat (ii, cpara, lpara, werte, maxw, .false.) 
            IF (ier_num.eq.0) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, wwerte, maxw) 
               is = int (werte (1) ) 
               IF (is.ne. - 1) then 
                  DO i = 1, ii 
                  is = int (werte (i) ) 
                  DO j = 1, ianz 
                  move (j, is) = wwerte (j) 
                  ENDDO 
                  ENDDO 
               ELSE 
                  DO i = 1, cr_nscat 
                  DO j = 1, ianz 
                  move (j, i) = wwerte (j) 
                  ENDDO 
                  ENDDO 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ Molecules                                                       
!                                                                       
      ELSE 
         IF ((dimen==3 .AND.ianz.eq.4) .OR.           &
             (dimen==4 .AND. (ianz==4 .OR. ianz==5))) then 
            IF (cpara (1) (1:1) .eq.'a') then 
               is = - 1 
            ELSE 
               CALL ber_params (1, cpara, lpara, werte, maxw) 
               is = nint (werte (1) ) 
               IF (is.le.0.or.is.gt.mole_num_type) then 
                  ier_num = - 64 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
!                                                                       
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, wwerte, maxw) 
               IF (is.ne. - 1) then 
                  DO j = 1, ianz 
                  move (j, is) = wwerte (j) 
                  ENDDO 
               ELSE 
                  DO i = 1, mole_num_type 
                  DO j = 1, ianz 
                  move (j, i) = wwerte (j) 
                  ENDDO 
                  ENDDO 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_set_move                   
!*****7*****************************************************************
!     SUBROUTINE rmc_set_mode (imode, ilocal, ianz, cpara, lpara, werte,maxw) 
      SUBROUTINE rmc_set_mode (               ianz, cpara, lpara, werte,maxw) 
!+                                                                      
!     Sets RMC/MC mode                                                  
!-                                                                      
      USE discus_config_mod 
      USE rmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN) :: maxw 
!                                                                       
!     INTEGER           , INTENT(OUT) :: imode
!     INTEGER           , INTENT(OUT) :: ilocal 
      INTEGER           , INTENT(IN ) :: ianz
      CHARACTER (LEN=* ), DIMENSION(1:MAXW), INTENT(IN ) :: cpara  !(maxw) 
      INTEGER           , DIMENSION(1:MAXW), INTENT(IN ) :: lpara  !(maxw) 
      REAL              , DIMENSION(1:MAXW), INTENT(IN ) :: werte  !(maxw) 
!
      INTEGER :: imode
      INTEGER :: i
      REAL    :: sump
!                                                                       
      IF (ianz.ge.1) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3) .eq.'SHI') then 
            imode = rmc_mode_shift 
         ELSEIF (cpara (1) (1:3) .eq.'SWD') then 
            imode = rmc_mode_swdisp 
         ELSEIF (cpara (1) (1:3) .eq.'SWC') then 
            imode = rmc_mode_swchem 
         ELSEIF (cpara (1) (1:3) .eq.'EXT') then 
            imode = rmc_mode_external 
         ELSE 
            ier_typ = ER_RMC 
            ier_num = - 9 
         ENDIF 
!                                                                       
         IF (ianz.eq.2) then 
            CALL do_cap (cpara (2) ) 
            IF (cpara (2) (1:1) .eq.'A') then 
               rmc_move_local(imode) = rmc_local_all 
            ELSEIF (cpara (2) (1:1) .eq.'L') then 
               rmc_move_local(imode) = rmc_local_loc 
            ELSEIF (cpara (2) (1:2) .eq.'SL') then 
               rmc_move_local(imode) = rmc_local_locsite 
            ELSEIF (cpara (2) (1:2) .eq.'SI') then 
               rmc_move_local(imode) = rmc_local_site 
            ELSEIF (cpara (2) (1:2) .eq.'CO') then 
               rmc_move_local(imode) = rmc_local_conn 
            ELSE 
               ier_typ = ER_RMC 
               ier_num = - 9 
            ENDIF 
         ELSE 
            rmc_move_local(imode) = rmc_local_all 
         ENDIF 
!                                                                       
!     --Set probabilities for the different moves                       
!                                                                       
         rmc_move_prob (imode) = werte (1) 
         sump = 0.0 
         DO i = 1, RMC_N_MOVE 
         sump = sump + rmc_move_prob (i) 
         ENDDO 
         IF (sump.gt.0) then 
            rmc_move_cprob (1) = rmc_move_prob (1) / sump 
            DO i = 2, RMC_N_MOVE 
            rmc_move_cprob (i) = rmc_move_cprob (i - 1) + &
                                 rmc_move_prob  (i) / sump                                                  
            ENDDO 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_set_mode                   
!*****7*****************************************************************
      SUBROUTINE rmc_set_mindist (is, js, dist) 
!+                                                                      
!     Set minimal allowed atom distances                                
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE rmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: is
      INTEGER , INTENT(IN) :: js
      REAL    , INTENT(IN) :: dist 
!
      INTEGER :: ii, jj 
!                                                                       
                                                                        
      IF (is.ne. - 1.and.js.ne. - 1) then 
         rmc_mindist (is, js) = dist 
         rmc_mindist (js, is) = dist 
      ELSEIF (is.eq. - 1.and.js.ne. - 1) then 
         DO ii = 1, cr_nscat 
         rmc_mindist (ii, js) = dist 
         rmc_mindist (js, ii) = dist 
         ENDDO 
      ELSEIF (is.ne. - 1.and.js.eq. - 1) then 
         DO ii = 1, cr_nscat 
         rmc_mindist (ii, is) = dist 
         rmc_mindist (is, ii) = dist 
         ENDDO 
      ELSE 
         DO ii = 1, cr_nscat 
         DO jj = 1, cr_nscat 
         rmc_mindist (ii, jj) = dist 
         rmc_mindist (jj, ii) = dist 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_set_mindist                
!*****7*****************************************************************
      SUBROUTINE rmc_rvalue 
!+                                                                      
!     Calculate R-values                                                
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE rmc_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      REAL sumd, sumo, tsumd, tsumo 
      REAL cint, dint, oo, dd, r4 
      INTEGER isym (rmc_max_planes) 
      INTEGER ip, is, iq 
!                                                                       
      tsumd = 0.0 
      tsumo = 0.0 
!                                                                       
      DO ip = 1, rmc_nplane 
      IF (rmc_dosym) then 
         isym (ip) = rmc_nsym (ip) 
      ELSE 
         isym (ip) = 1 
      ENDIF 
      ENDDO 
!                                                                       
      WRITE (output_io, 1000) 
      DO ip = 1, rmc_nplane 
      DO is = 1, isym (ip) 
      sumo = 0.0 
      sumd = 0.0 
      CALL rmc_inten (ip, is, .false.) 
      DO iq = 1, rmc_num (1, ip) * rmc_num (2, ip) 
      cint = rmc_back (ip) + rmc_skal (ip) * dsi (iq) 
      dint = rmc_int (offq (ip) + iq) - cint 
      oo = rmc_wic (offq (ip) + iq) * rmc_int (offq (ip) + iq) **2 
      dd = rmc_wic (offq (ip) + iq) * dint**2 
      sumo = sumo + oo 
      sumd = sumd+dd 
      tsumo = tsumo + oo 
      tsumd = tsumd+dd 
      ENDDO 
      r4 = 0.0 
      IF (sumo.ne.0.0) r4 = sqrt (sumd / sumo) 
      WRITE (output_io, 1500) ip, rmc_nplane, is, isym (ip), r4 * 100. 
      ENDDO 
      ENDDO 
      r4 = 0.0 
      IF (tsumo.ne.0.0) r4 = sqrt (tsumd / tsumo) 
      WRITE (output_io, 1600) r4 * 100. 
!                                                                       
 1000 FORMAT     (' R-values for current RMC results : ') 
 1500 FORMAT     ('    Plane ',i3,'/',i3,', sym. op. ',i2,'/',i2,' : ', &
     &                   g12.6,' %')                                    
 1600 FORMAT     ('                    Total R-value : ',g12.6,' %') 
      END SUBROUTINE rmc_rvalue                     
!*****7*****************************************************************
      SUBROUTINE rmc_show (cmd) 
!+                                                                      
!     Show current RMC settings                                         
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE diffuse_mod 
      USE molecule_mod 
      USE rmc_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: cmd
! 
      CHARACTER(25) wic (8) 
      CHARACTER(9) at_lis (0:maxscat+1) 
      CHARACTER(8) ra 
      INTEGER mol_lis (mole_max_type+1) 
      INTEGER i, j, k, ie, len_str 
!                                                                       
      CHARACTER(9) at_name_i, at_name_j 
!     CHARACTER(9) at_name 
!                                                                       
      DATA wic / 'w(Int) = 1.0           ', 'w(Int) = sqrt(Int)     ', '&
     &w(Int) = log(Int)      ', 'w(Int) = Int           ', 'w(Int) = Int&
     &**2        ', 'w(Int) = 1/Int         ', 'w(Int) = 1/sqrt(Int)   '&
     &, 'w(Int) = read from file' /                                     
!                                                                       
      IF (cmd (1:3) .eq.'ALL'.or.cmd (1:3) .eq.'MOD') then 
         IF (rmc_mode.eq.rmc_mode_shift) then 
            WRITE (output_io, 1100) 'move atoms/molecules' 
         ELSEIF (rmc_mode.eq.rmc_mode_swchem) then 
            WRITE (output_io, 1100) 'switch atoms/molecules' 
         ELSEIF (rmc_mode.eq.rmc_mode_swdisp) then 
            WRITE (output_io, 1100) 'switch displacements' 
         ELSEIF (rmc_mode.eq.rmc_mode_external) then 
      WRITE (output_io, 1100) 'defined in user supplied subroutine' 
         ELSE 
            WRITE (output_io, 1100) '??????????' 
         ENDIF 
         IF (rmc_local.eq.rmc_local_all) then 
            WRITE (output_io, 1120) 'all' 
         ELSEIF (rmc_local.eq.rmc_local_loc) then 
            WRITE (output_io, 1120) 'only local (+-1 unit cell)' 
         ELSEIF (rmc_local.eq.rmc_local_locsite) then 
            WRITE (output_io, 1120) 'only local, same site' 
         ELSEIF (rmc_local.eq.rmc_local_site) then 
            WRITE (output_io, 1120) 'all, same site' 
         ELSE 
            WRITE (output_io, 1120) '??????????' 
         ENDIF 
         IF (rmc_doskal.and.rmc_doback) then 
            WRITE (output_io, 1150) 'scaling factor and background' 
         ELSEIF (rmc_doskal.and..not.rmc_doback) then 
            WRITE (output_io, 1150) 'scaling factor only' 
         ELSEIF (rmc_doback.and..not.rmc_doskal) then 
            WRITE (output_io, 1150) 'background only' 
         ELSE 
            WRITE (output_io, 1150) 'none' 
         ENDIF 
!                                                                       
         IF (rmc_ilots.eq.LOT_OFF) then 
            WRITE (output_io, 1200) 
         ELSEIF (rmc_ilots.eq.LOT_BOX) then 
            WRITE (output_io, 1210) rmc_nlots 
            WRITE (output_io, 1230) (ls_xyz (i), i = 1, 3), lperiod 
            IF (rmc_ranloc) then 
               WRITE (output_io, 1235) 'random lots' 
            ELSE 
               WRITE (output_io, 1235) rmc_lname (1:len_str (rmc_lname) &
               )                                                        
            ENDIF 
         ELSEIF (rmc_ilots.eq.LOT_ELI) then 
            WRITE (output_io, 1220) rmc_nlots 
            WRITE (output_io, 1230) (ls_xyz (i), i = 1, 3), lperiod 
            IF (rmc_ranloc) then 
               WRITE (output_io, 1235) 'random lots' 
            ELSE 
               WRITE (output_io, 1235) rmc_lname (1:len_str (rmc_lname) &
               )                                                        
            ENDIF 
         ENDIF 
!                                                                       
         IF (rmc_ave.gt.0.0) WRITE (output_io, 1240) rmc_ave * 100.0 
         WRITE (output_io, 1250) rmc_maxcyc 
         WRITE (output_io, 1300) rmc_display, rmc_log 
         WRITE (output_io, 1400) rmc_sigma 
      ENDIF 
!                                                                       
      IF (cmd (1:3) .eq.'ALL'.or.cmd (1:3) .eq.'DAT') then 
         IF (rmc_data.eq.rmc_data_nipl) then 
            WRITE (output_io, 1500) 'NIPL' 
         ELSEIF (rmc_data.eq.rmc_data_pgm) then 
            WRITE (output_io, 1500) 'PGM ' 
         ELSE 
            WRITE (output_io, 1500) '????' 
         ENDIF 
         DO i = 1, rmc_nplane 
         WRITE (output_io, 1800) i, rmc_num (1, i), rmc_num (2, i),     &
         rmc_fname (i) (1:len_str (rmc_fname (i) ) )                    
!        ra = 'Neutrons' 
!        IF (rmc_lxray (i) ) ra = 'X-rays' 
         SELECTCASE(rmc_radiation (i))
            CASE(RMC_RAD_XRAY)
               ra = 'X-rays'
            CASE(RMC_RAD_NEUT)
               ra = 'Neutrons'
            CASE(RMC_RAD_ELEC)
               ra = 'Electron'
         END SELECT
      IF (rmc_lambda (i) .eq.'    ') then 
            WRITE (output_io, 1850) rmc_rlambda (i), ra, rmc_ldbw (i),  &
            rmc_ano (i)                                                 
         ELSE 
            WRITE (output_io, 1860) rmc_lambda (i), ra, rmc_ldbw (i),   &
            rmc_ano (i)                                                 
         ENDIF 
         WRITE (output_io, 1900) wic (rmc_wic_typ (i) ) 
         IF (rmc_constrain (i) .ne.i) then 
            WRITE (output_io, 1950) rmc_constrain (i) 
         ENDIF 
         ENDDO 
         WRITE (output_io, 1970) 2.0 * rmc_qmin, 2.0 * rmc_qmax, 2.0 *  &
         rmc_llim, 2.0 * rmc_ulim                                       
      ENDIF 
!                                                                       
      IF (cmd (1:3) .eq.'ALL'.or.cmd (1:3) .eq.'SYM') then 
         WRITE (output_io, 2000) 
         DO i = 1, rmc_nplane 
         ie = rmc_nsym (i) 
         IF (.not.rmc_dosym) ie = 1 
         DO j = 1, ie 
         WRITE (output_io, 2020) i, rmc_nplane, j, rmc_nsym (i) 
         WRITE (output_io, 2030) (rmc_eck (k, 1, j, i), k = 1, 3) 
         WRITE (output_io, 2040) (rmc_eck (k, 2, j, i), k = 1, 3) 
         WRITE (output_io, 2050) (rmc_eck (k, 3, j, i), k = 1, 3) 
         WRITE (output_io, 2060) (rmc_vi (k, 1, j, i), k = 1, 3) 
         WRITE (output_io, 2070) (rmc_vi (k, 2, j, i), k = 1, 3) 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      IF (cmd (1:3) .eq.'ALL'.or.cmd (1:3) .eq.'ATO') then 
         WRITE (output_io, 2200) 
         DO i = 1, cr_nscat 
         DO j = i, cr_nscat 
         at_name_i = at_name (i) 
         at_name_j = at_name (j) 
         WRITE (output_io, 2300) at_name_i, at_name_j, rmc_mindist (i,  &
         j)                                                             
         ENDDO 
         ENDDO 
!                                                                       
         WRITE (output_io, 2400) 
         IF (rmc_sel_atom) then 
            DO i = 1, cr_nscat 
            at_name_i = at_name (i) 
            WRITE (output_io, 2500) at_name_i, (rmc_maxmove (j, i),     &
            j = 1, 3)                                                   
            ENDDO 
!                                                                       
            j = 0 
            DO i = 0, cr_nscat 
            IF (rmc_allowed (i) ) then 
               j = j + 1 
               at_lis (j) = at_name (i) 
            ENDIF 
            ENDDO 
            WRITE (output_io, 2700) (at_lis (k), k = 1, j) 
!                                                                       
         ELSE 
            DO i = 1, mole_num_type 
            WRITE (output_io, 2600) i, (rmc_maxmove (j, i), j = 1, 3) 
            ENDDO 
!                                                                       
            j = 0 
            DO i = 1, mole_num_type 
            IF (rmc_allowed (i) ) then 
               j = j + 1 
               mol_lis (j) = i 
            ENDIF 
            ENDDO 
            WRITE (output_io, 2800) (mol_lis (k), k = 1, j) 
!                                                                       
         ENDIF 
      ENDIF 
!                                                                       
 1100 FORMAT     (  ' RMC running mode            : ',A) 
 1120 FORMAT     (  '   Valid neighbours          : ',A) 
 1150 FORMAT     (  '   Calc. scaling/backgr.     : ',A) 
 1200 FORMAT     (  '   Fourier volume            : complete crystal') 
 1210 FORMAT     (  '   Fourier volume            : ',I4,               &
     &                       ' box shaped lots')                        
 1220 FORMAT     (  '   Fourier volume            : ',I4,               &
     &                       ' ellipsoid shaped lots')                  
 1230 FORMAT (  '   Lot size                  : ',I3,' x ',I3,' x ',I3, &
     &                       ' unit cells (per. bound. = ',L1,')')      
 1235 FORMAT     (  '   Lot mode (random / file)  : ',A) 
 1240 FORMAT (  '   Calc. of <F> based on     : ',F5.1,' % of crystal') 
 1250 FORMAT     (  '   Max. number of RMC cycles : ',I8) 
 1300 FORMAT     (  '   Output update intervall   : ',I8,               &
     &                     '   (logging to rmc.log=',L1,')')            
 1400 FORMAT     (  '   SIGMA for CHI calculation : ',F8.4) 
 1500 FORMAT     (/,' Experimental data (',A4,')    : ') 
 1800 FORMAT     (/,' Plane ',I2,' (',I3,' x ',I3,' points) : ',A) 
 1850 FORMAT     (  '    Method (wavelen.: ',F5.3,'A): ',A8,            &
     &              ' (Debye-Waller=',L1,',anomalous scat.=',L1,')')    
 1860 FORMAT     (  '    Method (wavelen.: ',A6,'): ',A8,               &
     &              ' (Debye-Waller=',L1,',anomalous scat.=',L1,')')    
 1900 FORMAT     (  '       Weighting scheme used : ',A) 
 1950 FORMAT     (  '    Scal. / back. from plane : ',I2) 
 1970 FORMAT     (/,' Q-range of input data       : ',F6.3,' - ',F6.3,  &
     &                     ' used : ',F6.3,' - ',F6.3,' [1/A]')         
 2000 FORMAT     (/,' Symmetry of input data      : ') 
 2020 FORMAT     (  '       Reciprocal data layer : plane ',I2,'/',I2,  &
     &                     ', sym.op. ',I2,'/',I2)                      
 2030 FORMAT     (  '           lower left corner : ',3(F7.3,1X)) 
 2040 FORMAT     (  '          lower right corner : ',3(F7.3,1X)) 
 2050 FORMAT     (  '           upper left corner : ',3(F7.3,1X)) 
 2060 FORMAT     (  '              increment hor. : ',3(F7.3,1X)) 
 2070 FORMAT     (  '             increment vert. : ',3(F7.3,1X),/) 
 2200 FORMAT     (  ' Minimal allowed distances   : ') 
 2300 FORMAT     (  '     ',A9,' and ',A9,' : ',F7.3) 
 2400 FORMAT     (/,' Sigma for created RMC moves : ') 
 2500 FORMAT     (  '              Atom ',A9,' : ',3(F7.3,1X)) 
 2600 FORMAT     (  '     Molecule type ',I9,' : ',3(F7.3,1X)) 
 2700 FORMAT     (/,' Selected atoms for RMC      : ',100(A9,1X)) 
 2800 FORMAT     (/,' Selected molecules for RMC  : ',100(I4,1X)) 
      END SUBROUTINE rmc_show                       
!*****7*****************************************************************
      SUBROUTINE rmc_save (zeile, lp) 
!+                                                                      
!     save scattering amplitudes / structure                            
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE rmc_mod 
      USE save_menu
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: maxw = 10 
!
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
      INTEGER          , INTENT(INOUT) :: lp 
!                                                                       
      INTEGER i, il
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
!------ --save structure                                                
!                                                                       
         IF (cdummy (1:3) .eq.'STR') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            WRITE (output_io, 1500) cpara (1) 
            CALL save_struc (cpara (1), lpara (1) ) 
!                                                                       
!------ --save lot origins                                              
!                                                                       
         ELSEIF (cdummy (1:3) .eq.'LOT') then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            WRITE (output_io, 1505) cpara (1) 
            CALL oeffne (57, cpara (1) , 'unknown') 
            IF (ier_num.eq.0) then 
               DO il = 1, rmc_nlots 
               WRITE (57, '(3(I4,1X))') (rmc_lots_orig (i, il) , i = 1, &
               3)                                                       
               ENDDO 
               CLOSE (57) 
            ENDIF 
!                                                                       
!------ --save scattering intensities                                   
!                                                                       
         ELSEIF (cdummy (1:3) .eq.'SCA') then 
            IF (ianz.ge.3) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (1, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               ip = nint (werte (1) ) 
               IF (ip.le.0.or.ip.gt.rmc_nplane) then 
                  ier_num = - 7 
                  ier_typ = ER_RMC 
               ELSE 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL do_build_name (ianz, cpara, lpara, werte, maxw,  &
                  1)                                                    
                  IF (ier_num.ne.0) return 
                  WRITE (output_io, 1510) ip, cpara (1) 
                  CALL rmc_writedat (ip, 1, cpara (1) ) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!------ --save scattering intensities of sym. equivalent planes         
!                                                                       
         ELSEIF (cdummy (1:3) .eq.'SYM') then 
            IF (ianz.ge.4) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (2, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               ip = nint (werte (1) ) 
               is = nint (werte (2) ) 
               IF (ip.le.0.or.ip.gt.rmc_nplane) then 
                  ier_num = - 7 
                  ier_typ = ER_RMC 
               ELSE 
                  IF (is.le.0.or.is.gt.rmc_nsym (ip) ) then 
                     ier_num = - 8 
                     ier_typ = ER_RMC 
                  ELSE 
                     CALL del_params (2, ianz, cpara, lpara, maxw) 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
                     IF (ier_num.ne.0) return 
                     WRITE (output_io, 1520) ip, is, cpara (1) 
                     CALL rmc_writedat (ip, is, cpara (1) ) 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
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
 1500 FORMAT     (' Saving structure to file : ',A30) 
 1505 FORMAT     (' Saving lot origins to file : ',A30) 
 1510 FORMAT     (' Saving intensities of plane ',I2,' to file : ',A30) 
 1520 FORMAT     (' Saving intensities (ip:',I2,'/is:',I2,              &
     &                   ') to file: ',A30)                             
      END SUBROUTINE rmc_save                       
!*****7*****************************************************************
      SUBROUTINE rmc_writedat (ip, is, fname) 
!+                                                                      
!     Save scattering data to disk                                      
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE rmc_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: fname 
      INTEGER          , INTENT(IN) :: ip
      INTEGER          , INTENT(IN) :: is
!
      CHARACTER(LEN=1024) :: message
      INTEGER             :: ios
      INTEGER pgmmax, i, j 
!                                                                       
      pgmmax = 255 
!                                                                       
!------ first convert to intensities and check for points w=0           
!                                                                       
      CALL rmc_inten (ip, is, .false.) 
!                                                                       
      DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
      dsi (i) = rmc_back (ip) + rmc_skal (ip) * dsi (i) 
      IF (rmc_wic (offq (ip) + i) .eq.0.0) then 
         IF (rmc_data.eq.rmc_data_nipl) then 
            dsi (i) = - 9999.0D0
         ELSEIF (rmc_data.eq.rmc_data_pgm) then 
            dsi (i) = 0.0D0 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      OPEN (unit = 44, file = fname, status = 'unknown', &
            IOSTAT=IOS, IOMSG=message) 
      IF(ios/=0) THEN
         ier_num = -2
         ier_typ = ER_IO
         ier_msg(1)(1:60) = fname(1:60)
         ier_msg(3) = message(1:80)
         RETURN
      ENDIF
!                                                                       
!------ NIPL file                                                       
!                                                                       
      IF (rmc_data.eq.rmc_data_nipl) then 
         WRITE (44, * ) (rmc_num (i, ip), i = 1, 2) 
         WRITE (44, * ) (rmc_xy (i, ip), i = 1, 4) 
         DO j = 1, rmc_num (2, ip) 
         WRITE (44, * ) (dsi ( (i - 1) * rmc_num (2, ip) + j), i = 1,   &
         rmc_num (1, ip) )                                              
         ENDDO 
!                                                                       
!------ PGM file                                                        
!                                                                       
      ELSEIF (rmc_data.eq.rmc_data_pgm) then 
         WRITE (44, 2000) 
         WRITE (44, 2100) 'lower left ', (rmc_eck (i, 1, is, ip) , i =  &
         1, 3)                                                          
         WRITE (44, 2100) 'lower right', (rmc_eck (i, 2, is, ip) , i =  &
         1, 3)                                                          
         WRITE (44, 2100) 'upper left ', (rmc_eck (i, 3, is, ip) , i =  &
         1, 3)                                                          
         WRITE (44, * ) (rmc_num (i, ip), i = 1, 2) 
         WRITE (44, * ) pgmmax 
         DO j = rmc_num (2, ip), 1, - 1 
         WRITE (44, * ) (nint (dsi ( (i - 1) * rmc_num (2, ip) + j) ),  &
         i = 1, rmc_num (1, ip) )                                       
         ENDDO 
      ENDIF 
!                                                                       
      CLOSE (44) 
!                                                                       
 2000 FORMAT     ('P2'/'# DISCUS RMC resulting plane',/,'#') 
 2100 FORMAT     ('# ',A11,' : ',3(F7.3,1X)) 
      END SUBROUTINE rmc_writedat                   
!*****7*****************************************************************
      SUBROUTINE rmc_readdata (zeile, lp) 
!+                                                                      
!     read experimental data for rmc fit                                
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod
      USE metric_mod
      USE rmc_mod 
      USE rmc_symm_mod
!
      USE debug_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER, PARAMETER :: maxw = 13 
      INTEGER, PARAMETER :: max_sym = 48 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
      INTEGER          , INTENT(INOUT) :: lp 
!                                                                       
      CHARACTER(1024) cpara (maxw), cfile, cwic 
      CHARACTER(4) cdummy 
      INTEGER lpara (maxw) , cfile_l
      INTEGER ianz, nsym, rsym 
      INTEGER ip, nx, ny, wx, wy 
      INTEGER i, j, k 
      REAL werte (maxw), d 
      REAL e1 (3), e2 (3), e3 (3), vi1 (3), vi2 (3), z (3) 
      REAL qmin, qmax 
      REAL ee1 (4, max_sym) 
      REAL ee2 (4, max_sym) 
      REAL ee3 (4, max_sym) 
      REAL zone (4, max_sym) 
      REAL mat (4, 4, max_sym) 
!
      INTEGER  :: n_planes=1 ! Number of planes
      INTEGER  :: n_qxy   =1 ! Data points in reciprocal space
      INTEGER  :: n_sq    =1 ! Data points in reciprocal space*planes
      INTEGER  :: n_natoms=1 ! Maximum number of atoms for DIFFUSE allocation
      INTEGER  :: n_nscat =1 ! Maximum number of atoms for DIFFUSE allocation
      INTEGER  :: n_sym   =48! Maximum number of symmetry operations in reciprocal space
!
      LOGICAL lexist 
!                                                                       
!     REAL rmc_dowic 
!
!
!------ Allocate RMC PLANES
!
      IF(rmc_nplane == RMC_MAX_PLANES ) THEN
         n_planes = RMC_MAX_PLANES + 5
         n_sym    = 48       ! These arrays are small, use maximum symmetry number
         CALL alloc_rmc_planes(n_planes, n_sym)
      ENDIF
!                                                                       
!------ check if there is space for another plane                       
!                                                                       
      IF (rmc_nplane.eq.rmc_max_planes) then 
         ier_num = - 1 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
      ip = rmc_nplane+1 
!                                                                       
!------ work through all the parameters given                           
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.ge.13) then 
!                                                                       
!------ - get exp. method (neutron/x-ray)                               
!                                                                       
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:1) .eq.'N') then 
            rmc_lxray (ip) = .false. 
            rmc_radiation (ip) = RMC_RAD_NEUT
         ELSEIF (cpara (1) (1:1) .eq.'X') then 
            rmc_lxray (ip) = .true. 
            rmc_radiation (ip) = RMC_RAD_XRAY
         ELSE 
            ier_num = - 5 
            ier_typ = ER_RMC 
            RETURN 
         ENDIF 
!                                                                       
!------ - get wavelength (cua1 or real number)                          
!                                                                       
         CALL do_cap (cpara (2) ) 
         IF (ichar ('A') .le.ichar (cpara (2) (1:1) ) .and.ichar (cpara &
         (2) (1:1) ) .le.ichar ('Z') ) then                             
            rmc_lambda (ip) = cpara (2) (1:lpara(2))
         ELSE 
            CALL ber_params (1, cpara (2), lpara (2), werte, maxw) 
            rmc_rlambda (ip) = werte (1) 
            rmc_lambda (ip) = ' ' 
         ENDIF 
!                                                                       
!     - build file name                                                 
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 3) 
         IF (ier_num.eq.0) then 
            cfile = cpara (3) 
            cfile_l = lpara(3)
         ELSE 
            RETURN 
         ENDIF 
         cwic = cpara (4) 
         CALL del_params (4, ianz, cpara, lpara, maxw) 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         DO i = 1, 3 
         e1 (i) = werte (i) 
         e2 (i) = werte (i + 3) 
         e3 (i) = werte (i + 6) 
         ENDDO 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
!------ check if specified file is present                              
!                                                                       
      INQUIRE (file = cfile, exist = lexist) 
      IF (.not.lexist) then 
         ier_num = - 1 
         ier_typ = ER_IO 
         RETURN 
      ENDIF 
!                                                                       
!------ get weighting scheme                                            
!                                                                       
      cdummy = cwic (1:4) 
      CALL do_cap (cdummy) 
      IF (cdummy.eq.'ONE ') then 
         rmc_wic_typ (ip) = rmc_wic_eins 
      ELSEIF (cdummy.eq.'SQRT') then 
         rmc_wic_typ (ip) = rmc_wic_sqrt 
      ELSEIF (cdummy.eq.'LOG ') then 
         rmc_wic_typ (ip) = rmc_wic_log 
      ELSEIF (cdummy.eq.'LIN ') then 
         rmc_wic_typ (ip) = rmc_wic_lin 
      ELSEIF (cdummy.eq.'SQUA') then 
         rmc_wic_typ (ip) = rmc_wic_qua 
      ELSEIF (cdummy.eq.'INV ') then 
         rmc_wic_typ (ip) = rmc_wic_inv 
      ELSEIF (cdummy.eq.'ISQ ') then 
         rmc_wic_typ (ip) = rmc_wic_isq 
      ELSE 
         INQUIRE (file = cwic, exist = lexist) 
         IF (lexist) then 
            rmc_wic_typ (ip) = rmc_wic_dat 
            CALL oeffne (18, cwic, 'old') 
         ELSE 
            ier_num = - 16 
            ier_typ = ER_RMC 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!------ read data for given plane                                       
!                                                                       
      rmc_fname (ip) = cfile (1:cfile_l)
      CALL oeffne (17, cfile, 'old') 
      IF (ier_num.ne.0) return 
!                                                                       
!------ NIPL file                                                       
!                                                                       
      IF (rmc_data.eq.rmc_data_nipl) then 
         READ (17, *, end = 99, err = 999) nx, ny 
         READ (17, *, end = 99, err = 999) (rmc_xy (j, ip), j = 1, 4) 
         IF (rmc_wic_typ (ip) .eq.rmc_wic_dat) then 
            READ (18, *, end = 99, err = 999) wx, wy 
            READ (18, *, end = 99, err = 999) d, d, d, d 
            IF (wx.ne.nx.or.wy.ne.ny) then 
               ier_num = - 17 
               ier_typ = ER_RMC 
            ENDIF 
         ENDIF 
!                                                                       
         n_qxy  = MAX( n_qxy , rmc_n_qxy, RMC_MAX_Q , offq(ip)+nx*ny)
         CALL alloc_rmc_data ( n_qxy)
!
         IF (nx * ny.gt. (RMC_MAX_Q - offq (ip) ) ) then 
            ier_num = - 2 
            ier_typ = ER_RMC 
         ENDIF 
!                                                                       
         IF (ier_num.ne.0) goto 98 
!
         DO j = 1, ny 
         READ (17, *, end = 99, err = 999) (rmc_int (offq (ip) +        &
         (i - 1) * ny + j), i = 1, nx)                                  
         IF (rmc_wic_typ (ip) .eq.rmc_wic_dat) then 
            READ (18, *, end = 99, err = 999) (rmc_wic (offq (ip)       &
            + (i - 1) * ny + j), i = 1, nx)                             
         ELSE 
            DO i = 1, nx 
            k = offq (ip) + (i - 1) * ny + j 
            rmc_wic (k) = rmc_dowic (rmc_wic_typ (ip), rmc_int (k) ) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!------ PGM file                                                        
!                                                                       
      ELSEIF (rmc_data.eq.rmc_data_pgm) then 
         CALL rmc_pgmheader (17, nx, ny, rmc_xy (1, ip), rmc_xy (2, ip),&
         rmc_xy (3, ip), rmc_xy (4, ip) )                               
         IF (rmc_wic_typ (ip) .eq.rmc_wic_dat) then 
            CALL rmc_pgmheader (18, wx, wy, d, d, d, d) 
            IF (wx.ne.nx.or.wy.ne.ny) then 
               ier_num = - 17 
               ier_typ = ER_RMC 
            ENDIF 
         ENDIF 
!                                                                       
         n_qxy  = MAX( n_qxy , rmc_n_qxy, RMC_MAX_Q , offq(ip)+nx*ny)
         CALL alloc_rmc_data ( n_qxy)
!
         IF (nx * ny.gt. (RMC_MAX_Q - offq (ip) ) ) then 
            ier_num = - 2 
            ier_typ = ER_RMC 
         ENDIF 
!                                                                       
         IF (ier_num.ne.0) goto 98 
!                                                                       
         DO j = ny, 1, - 1 
         READ (17, *, end = 99, err = 999) (rmc_int (offq (ip) +        &
         (i - 1) * ny + j), i = 1, nx)                                  
         IF (rmc_wic_typ (ip) .eq.rmc_wic_dat) then 
            READ (18, *, end = 99, err = 999) (rmc_wic (offq (ip)       &
            + (i - 1) * ny + j), i = 1, nx)                             
         ELSE 
            DO i = 1, nx 
            k = offq (ip) + (i - 1) * ny + j 
            rmc_wic (k) = rmc_dowic (rmc_wic_typ (ip), rmc_int (k) ) 
            ENDDO 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
   98 CONTINUE 
      CLOSE (17) 
      IF (rmc_wic_typ (ip) .eq.rmc_wic_dat) close (18) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ calculate increments and zone axis                              
!                                                                       
      DO i = 1, 3 
      vi1 (i) = (e2 (i) - e1 (i) ) / float (max (1, nx - 1) ) 
      vi2 (i) = (e3 (i) - e1 (i) ) / float (max (1, ny - 1) ) 
      ENDDO 
      CALL vekprod (vi1, vi2, z, cr_reps, cr_gten) 
!                                                                       
!------ perform symmetry operations                                     
!                                                                       
      DO i = 1, 3 
      zone (i, 1) = z (i) 
      ee1 (i, 1) = e1 (i) 
      ee2 (i, 1) = e2 (i) 
      ee3 (i, 1) = e3 (i) 
      ENDDO 
      zone (4, 1) = 0.0 
      ee1 (4, 1) = 0.0 
      ee2 (4, 1) = 0.0 
      ee3 (4, 1) = 0.0 
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j, 1) = 0.0 
      ENDDO 
      mat (i, i, 1) = 1.0 
      ENDDO 
!                                                                       
      CALL rmc_symmetry (nsym, zone, mat, max_sym, .true., cr_acentric) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (nsym.gt.RMC_MAX_SYM) then 
         rmc_dosym = .false. 
         rmc_nosym = .true. 
         rsym = nsym 
         nsym = 1 
      ELSE 
         rsym = nsym 
         rmc_nosym = .false. 
      ENDIF 
!                                                                       
!------ do transform q=gSg* and transform corners of input plane        
!                                                                       
      CALL rmc_gsg (mat, nsym, max_sym) 
      CALL rmc_trans (ee1, mat, nsym, max_sym) 
      CALL rmc_trans (ee2, mat, nsym, max_sym) 
      CALL rmc_trans (ee3, mat, nsym, max_sym) 
!                                                                       
!------ data accepted                                                   
!                                                                       
      rmc_nplane = ip 
      rmc_num (1, ip) = nx 
      rmc_num (2, ip) = ny 
      rmc_nsym (ip) = nsym 
!                                                                       
!------ set offq and offsq values                                       
!                                                                       
      DO i = 2, nsym 
      offsq (ip, i) = offsq (ip, 1) + (i - 1) * nx * ny 
      ENDDO 
      offsq (ip + 1, 1) = offsq (ip, nsym) + nx * ny 
      offq (ip + 1) = offq (ip) + nx * ny 
!                                                                       
      IF (dbg) then 
         DO i = 1, rmc_max_planes 
         DO j = 1, RMC_MAX_SYM    ! Here RMC_MAX_SYM=48 is OK
         WRITE (output_io, 9999) i, j, offq (i), offsq (i, j) 
         ENDDO 
         ENDDO 
      ENDIF 
 9999 FORMAT     (' DBG : ip:',i3,' is:',i3,' offq:',i8,' offsq:',i8) 
!                                                                       
!------ store plane increments and corners                              
!                                                                       
      DO i = 1, nsym 
      DO j = 1, 3 
      rmc_vi (j, 1, i, ip) = (ee2 (j, i) - ee1 (j, i) ) / float (max (1,&
      nx - 1) )                                                         
      rmc_vi (j, 2, i, ip) = (ee3 (j, i) - ee1 (j, i) ) / float (max (1,&
      ny - 1) )                                                         
      ENDDO 
      ENDDO 
      DO i = 1, nsym 
      DO j = 1, 3 
      rmc_eck (j, 1, i, ip) = ee1 (j, i) 
      rmc_eck (j, 2, i, ip) = ee2 (j, i) 
      rmc_eck (j, 3, i, ip) = ee3 (j, i) 
      ENDDO 
      ENDDO 
!
!------ Allocate initial Diffuse
!
      rmc_n_qxy = MAX(nx*ny,RMC_MAX_Q)   ! Save RMC required size for allocation prior to 'run'
      n_qxy     = MAX(n_qxy, nx*ny, RMC_MAX_Q, MAXQXY)
      n_natoms  = MAX(n_natoms, cr_natoms, DIF_MAXAT)
      n_nscat   = MAX(n_nscat, cr_nscat, DIF_MAXSCAT)
      call alloc_diffuse (n_qxy, n_nscat, n_natoms)
      rmc_n_sym = MAX(rmc_n_sym, nsym)    ! Save RMC actual number of symmetry operations
      n_sq      = n_qxy*rmc_n_sym
      rmc_n_sq  = n_sq                    ! Save RMC_actual number of data points * nsym
      CALL alloc_rmc_istl ( n_sq, n_nscat, rmc_nplane )
!                                                                       
!------ initial q range of data to be used is ALL data                  
!                                                                       
      CALL rmc_layer (1, ip) 
      CALL rmc_stltab (1, ip, .true.) 
!                                                                       
      qmin = ristl (offsq (ip, 1) + 1) * CFINC 
      qmax = ristl (offsq (ip, 1) + 1) * CFINC 
      DO i = 2, rmc_num (1, ip) * rmc_num (2, ip) 
      qmin = min (qmin, ristl (offsq (ip, 1) + i) * CFINC) 
      qmax = max (qmax, ristl (offsq (ip, 1) + i) * CFINC) 
      ENDDO 
!                                                                       
      rmc_qmin = min (rmc_qmin, qmin) 
      rmc_qmax = max (rmc_qmax, qmax) 
      rmc_llim = rmc_qmin 
      rmc_ulim = rmc_qmax 
!                                                                       
!     write comment                                                     
!                                                                       
      IF (rmc_lxray (ip) ) then 
         WRITE (output_io, 9000) 'x-ray', ip, rmc_num (1, ip) , rmc_num &
         (2, ip) , rsym                                                 
      ELSE 
         WRITE (output_io, 9000) 'neutron', ip, rmc_num (1, ip) ,       &
         rmc_num (2, ip) , rsym                                         
      ENDIF 
      RETURN 
!                                                                       
   99 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
  999 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 9000 FORMAT     (' Read ',A,'-data for plane ',I2,' (',I3,' x ',I3,    &
     &                   ' points / ',I2,' sym. equivalent planes)')    
      END SUBROUTINE rmc_readdata                   
!*****7*****************************************************************
      SUBROUTINE rmc_pgmheader (ifile, nx, ny, xmin, xmax, ymin, ymax) 
!+                                                                      
!     Reading PGM Header incl. possible comments                        
!-                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN)    :: ifile
      INTEGER, INTENT(INOUT) :: nx
      INTEGER, INTENT(INOUT) :: ny
      REAL   , INTENT(OUT)   :: xmin
      REAL   , INTENT(OUT)   :: xmax
      REAL   , INTENT(OUT)   :: ymin
      REAL   , INTENT(OUT)   :: ymax 
!
      CHARACTER(1024) zeile 
      CHARACTER(2) cmagic 
      INTEGER id 
!                                                                       
      READ (ifile, 1000, end = 99, err = 999) cmagic 
      IF (cmagic (1:2) .ne.'P2') then 
         ier_num = - 18 
         ier_typ = ER_RMC 
      ELSE 
    5    CONTINUE 
         READ (ifile, 1000, end = 99, err = 999) zeile 
         IF (zeile (1:1) .eq.'#') goto 5 
         BACKSPACE (ifile) 
         READ (ifile, *, end = 99, err = 999) nx, ny, id 
!                                                                       
         xmin = 0.0 
         xmax = float (nx - 1) 
         ymin = 0.0 
         ymax = float (ny - 1) 
      ENDIF 
      RETURN 
!                                                                       
   99 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
  999 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
 1000 FORMAT    (a) 
      END SUBROUTINE rmc_pgmheader                  
!*****7*****************************************************************
      REAL function rmc_dowic (typ, inte) 
!+                                                                      
!     Calulates the weight for intensity int                            
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: typ 
      REAL   , INTENT(IN) :: inte 
      rmc_dowic = 0.0
!                                                                       
      IF (inte.le.0.0) then 
         rmc_dowic = 0.0 
      ELSE 
         IF (typ.eq.rmc_wic_eins) then 
            rmc_dowic = 1.0 
         ELSEIF (typ.eq.rmc_wic_sqrt) then 
            rmc_dowic = sqrt (inte) 
         ELSEIF (typ.eq.rmc_wic_log) then 
            rmc_dowic = log (inte) 
         ELSEIF (typ.eq.rmc_wic_lin) then 
            rmc_dowic = inte 
         ELSEIF (typ.eq.rmc_wic_qua) then 
            rmc_dowic = inte**2 
         ELSEIF (typ.eq.rmc_wic_inv) then 
            rmc_dowic = 1.0 / inte 
         ELSEIF (typ.eq.rmc_wic_isq) then 
            rmc_dowic = 1.0 / sqrt (inte) 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION rmc_dowic                        
!*****7*****************************************************************
      SUBROUTINE rmc_run 
!+                                                                      
!     main rmc loop - called by run command                             
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE rmc_mod
      USE update_cr_dim_mod
      USE random_mod
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE precision_mod
      IMPLICIT none 
       
!                                                                       
      REAL(PREC_DP) rmc_cc, rmc_c, rmc_ce 
      REAL(PREC_DP) rmc_e (rmc_max_planes) 
      REAL(PREC_DP) rmc_ee (rmc_max_planes) 
      REAL chi2_new, chi2_old, sig2 
      REAL prob, psum, p2sum, pave, psig, pmax, pn 
      REAL chi2 (rmc_max_planes) 
      REAL skal (rmc_max_planes) 
      REAL back (rmc_max_planes) 
      REAL start, zeit, seknds 
      REAL p_new (3, rmc_max_atom) 
      INTEGER i_new (rmc_max_atom) 
      INTEGER isel (rmc_max_atom), natoms 
      INTEGER imol (rmc_max_atom) 
      INTEGER isym (rmc_max_planes) 
      INTEGER zh, zm, zs 
      INTEGER i, j, ip, iq, is, iii 
      INTEGER igen, itry, iacc_good, iacc_bad 
      LOGICAL loop, laccept 
      REAL    ::   rmc_energy = 0.00
      LOGICAL :: l_rmc_energy = .false.
!
!      INTEGER  :: n_qxy   =1 ! Data points in reciprocal space
      INTEGER  :: n_sq    =1 ! Data points in reciprocal space*planes
      INTEGER  :: n_lots  =1 ! DANumber of RMC Lots for Fourier
!      INTEGER  :: n_natoms=1 ! Maximum number of atoms for DIFFUSE allocation
!      INTEGER  :: n_nscat =1 ! Maximum number of atoms for DIFFUSE allocation
!                                                                       
      REAL ran1 
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
!------ check consistency of input data                                 
!                                                                       
      CALL rmc_check_input 
      IF (ier_num.ne.0) return 
!
!------	Allocate arrays related to SQ, LOTS
!
      n_sq   = rmc_n_qxy*rmc_n_sym   ! Could like wise be = rmc_n_sq
      n_lots = MAX( n_lots, rmc_nlots, RMC_MAX_LOTS)
      CALL alloc_rmc_q ( n_sq, n_lots)
!                                                                       
!------ sym ?                                                           
!                                                                       
      DO ip = 1, rmc_nplane 
      IF (rmc_dosym) then 
         isym (ip) = rmc_nsym (ip) 
      ELSE 
         isym (ip) = 1 
      ENDIF 
      ENDDO 
!                                                                       
!------ Write some start information                                    
!                                                                       
      WRITE (output_io, 1000) (cr_icc (i), i = 1, 3), cr_natoms 
!                                                                       
!------ calculate scattering amplitudes                                 
!                                                                       
      IF (rmc_calc_f) then 
         WRITE (output_io, 1100) 
         CALL rmc_calc_scat (isym) 
      ELSE 
         WRITE (output_io, 1150) 
      ENDIF 
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
      IF (rmc_log) THEN 
        call oeffne_append(34,'rmc.log','unknown',.false.) 
!DBG        open (34,file='rmc.log',status='unknown',access='append')   
        IF (rmc_mode.eq.rmc_mode_swchem) THEN 
          WRITE (34,2000) 'Mode SWCHEM' 
        ELSEIF (rmc_mode.eq.rmc_mode_swdisp) THEN 
          WRITE (34,2000) 'Mode SWDISP' 
        ELSEIF (rmc_mode.eq.rmc_mode_shift) THEN 
          WRITE (34,2000) 'Mode SHIFT' 
        ENDIF 
        CLOSE (34) 
      ENDIF 
!                                                                       
!------ calculate sums from exp. data needed and initial chi2           
!                                                                       
      chi2_old = 0.0 
      DO ip=1,rmc_nplane 
        rmc_wtot (ip)=0.0 
        rmc_e    (ip)=0.0 
        rmc_ee   (ip)=0.0 
        rmc_c  = 0.0 
        rmc_cc = 0.0 
        rmc_ce = 0.0 
!                                                                       
        DO is=1,isym(ip) 
!                                                                       
          call rmc_inten(ip,is,.false.) 
!                                                                       
          DO iq=1,rmc_num(1,ip)*rmc_num(2,ip) 
            IF ((ristl(offsq(ip,1)+iq)*CFINC).ge.rmc_llim .and.         &
     &          (ristl(offsq(ip,1)+iq)*CFINC).le.rmc_ulim    ) THEN
               iii=offq(ip)+iq 
               rmc_wtot(ip)=rmc_wtot(ip)+rmc_wic(iii) 
               rmc_e (ip)=rmc_e (ip)+rmc_wic(iii)*rmc_int(iii) 
               rmc_ee(ip)=rmc_ee(ip)+rmc_wic(iii)*rmc_int(iii)**2 
               rmc_c =rmc_c +rmc_wic(iii)*dsi(iq) 
               rmc_cc=rmc_cc+rmc_wic(iii)*dsi(iq)**2 
               rmc_ce=rmc_ce+rmc_wic(iii)*dsi(iq)*rmc_int(iii) 
            ENDIF 
          ENDDO 
        ENDDO 
        IF (rmc_e(ip).gt.1E-06) THEN 
          call rmc_calcskal(ip,rmc_wtot(ip),rmc_c,rmc_cc,rmc_ce,        &
                                 rmc_e(ip),skal,back)        
          chi2(ip)=real( rmc_ee(ip)+skal(ip)**2*rmc_cc+back(ip)**2*     &
                   rmc_wtot(ip)+2.0*(skal(ip)*back(ip)*rmc_c-back(ip)*  &
                   rmc_e(ip)-skal(ip)*rmc_ce) )
          chi2(ip)=chi2(ip)/rmc_wtot(ip) 
          chi2_old=chi2_old + chi2(ip) 
        ENDIF 
      ENDDO 
!                                                                       
      IF (chi2_old.lt.1E-06) THEN 
        ier_num = -14 
        ier_typ = ER_RMC 
        return 
      ENDIF 
!                                                                       
!------ Here is the main loop                                           
!                                                                       
      start = seknds(0.0) 
      IF (rmc_mode.eq.rmc_mode_swchem) THEN 
        WRITE (output_io,2000) 'Mode SWCHEM' 
      ELSEIF (rmc_mode.eq.rmc_mode_swdisp) THEN 
        WRITE (output_io,2000) 'Mode SWDISP' 
      ELSEIF (rmc_mode.eq.rmc_mode_shift) THEN 
        WRITE (output_io,2000) 'Mode SHIFT' 
      ELSEIF (rmc_mode.eq.rmc_mode_external) THEN 
        WRITE (output_io,2000) 'external defined mode' 
      ENDIF 
!                                                                       
      sig2 = rmc_sigma**2/2.0 
!                                                                       
      DO while (loop) 
        laccept = .true. 
        igen = igen + 1 
!                                                                       
!-------- generate move and check for limits                            
!                                                                       
        IF (rmc_sel_atom) THEN 
          call rmc_genmove (laccept,natoms,p_new,i_new,isel) 
        ELSE 
          call rmc_genmove_mol (laccept,natoms,p_new,i_new,isel,imol) 
        ENDIF 
        IF (ier_num.ne.0) return 
!                                                                       
!-------- calc new scattering amplitudes and chi2                       
!                                                                       
        IF (laccept) THEN 
          itry = itry + 1 
          chi2_new = 0.0 
          DO ip=1,rmc_nplane 
            rmc_c    = 0.0 
            rmc_cc   = 0.0 
            rmc_ce   = 0.0 
            call dlink(rmc_ano(ip),           &
     &                     rmc_lambda(ip),rmc_rlambda(ip),   &
                           rmc_energy, l_rmc_energy,                & 
                           rmc_radiation(ip), rmc_power(ip))
            DO is=1,isym(ip) 
              call rmc_fcalc (ip,is,natoms,i_new,p_new,isel) 
              call rmc_inten (ip,is,.true.) 
              DO iq=1,rmc_num(1,ip)*rmc_num(2,ip) 
                IF ((ristl(offsq(ip,is)+iq)*CFINC).ge.rmc_llim .and.    &
     &              (ristl(offsq(ip,is)+iq)*CFINC).le.rmc_ulim ) THEN
                  iii=offq(ip)+iq 
                  rmc_c =rmc_c +rmc_wic(iii)*dsi(iq) 
                  rmc_cc=rmc_cc+rmc_wic(iii)*dsi(iq)**2 
                  rmc_ce=rmc_ce+rmc_wic(iii)*dsi(iq)*rmc_int(iii) 
                ENDIF 
              ENDDO 
            ENDDO 
!                                                                       
!------ ----- if plane contains no valid q points at all                
!                                                                       
            IF (rmc_e(ip).lt.1E-06) goto 1111 
!                                                                       
            call rmc_calcskal(ip,rmc_wtot(ip),rmc_c,rmc_cc,rmc_ce,      &
     &                             rmc_e(ip),skal,back)      
            chi2(ip)=real( rmc_ee(ip)+skal(ip)**2*rmc_cc+back(ip)**2*   &
     &               rmc_wtot(ip)+2.0*(skal(ip)*back(ip)*rmc_c-         &
     &               back(ip)*rmc_e(ip)-skal(ip)*rmc_ce) )
            chi2(ip)=chi2(ip)/rmc_wtot(ip) 
            chi2_new=chi2_new + chi2(ip) 
!                                                                       
 1111           continue 
          ENDDO 
!                                                                       
!     ----Accept move ?                                                 
!                                                                       
          prob = chi2_new - chi2_old 
!                                                                       
          IF (prob.lt.0) THEN 
            laccept = .true. 
          ELSE 
            IF (sig2.gt.0.0) THEN 
              psum    = psum  + prob 
              p2sum   = p2sum  + prob**2 
              pmax    = max (pmax,prob) 
              pn      = pn + 1 
              prob    = exp(-prob/sig2) 
              laccept = (prob.gt.ran1(idum)) 
            ELSE 
              laccept = .false. 
            ENDIF 
          ENDIF 
!                                                                       
          IF (rmc_sigma.eq.-9999.) laccept = .true. 
!                                                                       
!------ ----if accepted make move                                       
!                                                                       
          IF (laccept) THEN 
            call  rmc_makemove (isym,natoms,i_new,p_new,isel,imol) 
            chi2_old = chi2_new 
            DO i=1,rmc_nplane 
              rmc_chi2(i)=chi2(i) 
              rmc_back(i)=back(i) 
              rmc_skal(i)=skal(i) 
            ENDDO 
            IF (prob.lt.0) THEN 
              iacc_good = iacc_good + 1 
            ELSE 
              iacc_bad  = iacc_bad + 1 
            ENDIF 
          ENDIF 
!                                                                       
!------ --WRITE info and terminate or loop again                        
!                                                                       
        ENDIF 
!                                                                       
        loop = (itry.lt.rmc_maxcyc) 
!                                                                       
        IF (igen.gt.1000*rmc_display .and. itry.eq.0) THEN 
          ier_num = -20 
          ier_typ = ER_RMC 
          loop    = .false. 
        ENDIF 
!                                                                       
        IF (mod(igen,rmc_display).eq.0 .or. .not.loop) THEN 
          IF (.not.loop) WRITE (output_io,1250) 
          WRITE (output_io,1300) igen,itry,iacc_good,iacc_bad,          &
     &                        chi2_old                                  
          DO i=1,rmc_nplane 
            WRITE (output_io,1400) i,1.0/rmc_skal(i),                   &
     &                                  rmc_back(i),rmc_chi2(i)         
          ENDDO 
          IF (rmc_log) THEN 
            call oeffne_append(34,'rmc.log','unknown',.false.) 
!DBG            open (34,file='rmc.log',status='unknown',access='append')
            IF (.not.loop) WRITE (34,1250) 
            WRITE (34,1300) igen,itry,iacc_good,iacc_bad,chi2_old 
            DO i=1,rmc_nplane 
              WRITE (34,1400) i,1.0/rmc_skal(i),rmc_back(i),rmc_chi2(i) 
            ENDDO 
            CLOSE (34) 
          ENDIF 
        ENDIF 
      ENDDO 
!                                                                       
!------ WRITE timing summary and "prob" statistics                      
!                                                                       
      IF (pn.gt.0) THEN 
        pave = psum / pn 
        psig = sqrt (p2sum/pn - (psum/pn)**2) 
        WRITE (output_io,3000) pave,psig,pmax 
      ENDIF 
!                                                                       
      zeit = seknds(start) 
      zh = int (zeit / 3600.) 
      zm = int ((zeit - zh*3600.) / 60.) 
      zs = int ( zeit - zh*3600 - zm*60.) 
!                                                                       
      IF (rmc_log) then 
         CALL oeffne_append (34, 'rmc.log', 'unknown', .false.) 
!DBG        open (34,file='rmc.log',status='unknown',access='append')   
         WRITE (34, 3000) pave, psig, pmax 
         WRITE (34, 4000) zh, zm, zs, zeit / float (itry) 
         CLOSE (34) 
      ENDIF 
!                                                                       
!------ save some results to res[i] block                               
!                                                                       
      res_para (0) = 8 
!                                                                       
      res_para (1) = chi2_old 
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
 1000 FORMAT (' Running RMC-Fit ...',//                                 &
     &          ' Size of model crystal     : ',I3,' x ',I3,' x ',I3,   &
     &          ' containing ',I9,' atoms')                             
 1100 FORMAT (/,' Calculating scattering amplitudes ... ') 
 1150 FORMAT (/,' Using old scattering amplitudes ...',/                &
     &          ' (Use command rese to force recalculation)'/)          
 1250 FORMAT (/,' ---- Final configuration ---- ') 
 1300 FORMAT (/,' Gen: ',I6,' try: ',I6,' acc: (good/bad): ',I6,        &
     &          ' / ',I6,' s2x2: ',G15.8)                               
 1400 FORMAT ('   Plane ',I2,': scal: ',G10.4,' / back: ',G10.4,        &
     &          ' / s2x2: ',G22.8)                                      
 2000 FORMAT (/,' Starting main RMC loop ...',/,                        &
     &          ' (',A,') ',/)                                          
 3000 FORMAT (/,' Delta Chi    : ave: ',G15.4,' sig: ',G15.4,           &
     &          ' max: ',G15.4)                                         
 4000 FORMAT (/,' Elapsed time : ',I4,' h ',I2,' min ',I2,' sec ',/     &
     &          ' Time/cycle   : ',F9.3,' sec',/)                       
      END SUBROUTINE rmc_run                        
!*****7*****************************************************************
      SUBROUTINE rmc_calc_scat (isym) 
!+                                                                      
!     This routine calculates initial Fourier transform                 
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE four_strucf_mod
      USE fourier_lmn_mod
      USE rmc_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER isym (rmc_max_planes) 
      INTEGER lbeg (3), ncell 
      INTEGER ip, iscat, nlot, i, k, iii 
      REAL    ::   rmc_energy = 0.00
      LOGICAL :: l_rmc_energy = .false.
!
      CALL four_cexpt 
      CALL rmc_zero 
      CALL four_csize (cr_icc, rmc_csize, lperiod, ls_xyz) 
!                                                                       
!------ Create random lots or read file                                 
!                                                                       
      IF (rmc_ranloc) then 
         DO nlot = 1, rmc_nlots 
         CALL four_ranloc (rmc_csize, lbeg) 
         DO i = 1, 3 
         rmc_lots_orig (i, nlot) = lbeg (i) 
         ENDDO 
         ENDDO 
      ELSE 
         CALL oeffne (57, rmc_lname, 'old') 
         IF (ier_num.ne.0) return 
         DO nlot = 1, rmc_nlots 
         READ (57, *, end = 99, err = 999) (rmc_lots_orig (i, nlot),    &
         i = 1, 3)                                                      
         ENDDO 
         CLOSE (57) 
      ENDIF 
!                                                                       
!------ Loop over all exp. data planes                                  
!                                                                       
loop_plane: DO ip = 1, rmc_nplane 
         CALL dlink (rmc_ano (ip), rmc_lambda (ip),        &
                     rmc_rlambda (ip),  rmc_energy, l_rmc_energy, &
                     rmc_radiation(ip), rmc_power(ip) )
         CALL rmc_formtab (ip, .true.) 
!                                                                       
!------ - Loop over all sym. equivalent planes                          
!                                                                       
         loop_sym: DO k = 1, isym (ip) 
            DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
               acsf (i) = cmplx (0.0D0, 0.0D0) 
            ENDDO 
!                                                                       
            CALL rmc_layer (k, ip) 
            CALL fourier_lmn(eck,vi,num,lmn,off_shift)
            CALL rmc_stltab (k, ip, .true.) 
            CALL four_aver (rmc_ilots, rmc_ave) 
!                                                                       
            iii = offsq (ip, k) 
!                                                                       
!------ --- Loop over all 'lots'                                        
!                                                                       
            DO nlot = 1, rmc_nlots 
               DO i = 1, 3 
                  lbeg (i) = rmc_lots_orig (i, nlot) 
               ENDDO 
!                                                                       
!------ ----- Loop over all atom types                                  
!                                                                       
            DO iscat = 1, cr_nscat 
               CALL four_getatm (iscat, rmc_ilots, lbeg, ncell) 
               CALL four_strucf (iscat, .true.) 
               DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
                  rmc_csf (iii + i, nlot) = rmc_csf (iii + i, nlot) + tcsf (i) 
               ENDDO 
            ENDDO 
!                                                                       
            DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
               rmc_csf (iii + i, nlot) = rmc_csf (iii + i, nlot) - acsf (i) 
            ENDDO 
         ENDDO 
         WRITE (output_io, 1200) ip, k 
         ENDDO  loop_sym
      ENDDO loop_plane
!                                                                       
!------ we don not need to calculate the initial Fourier again          
!                                                                       
      rmc_calc_f = .false. 
      RETURN 
!                                                                       
   99 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      CLOSE (57) 
      RETURN 
!                                                                       
  999 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      CLOSE (57) 
      RETURN 
!                                                                       
!                                                                       
 1200 FORMAT   (' Done plane ',I2,', sym.op. ',I2,' ...') 
      END SUBROUTINE rmc_calc_scat                  
!*****7*****************************************************************
      SUBROUTINE rmc_inten (ip, is, lnew) 
!+                                                                      
!     This routine copies average intensity for given plane and sym     
!     in array 'dsi' ...                                                
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN) :: ip
      INTEGER, INTENT(IN) :: is
      LOGICAL, INTENT(IN) :: lnew 
!                                                                       
      INTEGER iq, il 
!------ We want to average rmc_csf_new ..                               
!                                                                       
      IF (lnew) then 
!                                                                       
!------ - First copy lot 1 to dsi ..                                    
!                                                                       
         DO iq = 1, rmc_num (1, ip) * rmc_num (2, ip) 
         dsi (iq) = DBLE (rmc_csf_new (offsq (ip, is) + iq, 1) * conjg (&
         rmc_csf_new (offsq (ip, is) + iq, 1) ) )                       
         ENDDO 
!                                                                       
!------ - If there is more than 1 lot, add the rest ..                  
!                                                                       
         IF (rmc_nlots.gt.1) then 
            DO il = 2, rmc_nlots 
            DO iq = 1, rmc_num (1, ip) * rmc_num (2, ip) 
            dsi (iq) = dsi (iq) + DBLE (rmc_csf_new (offsq (ip, is)     &
            + iq, il) * conjg (rmc_csf_new (offsq (ip, is) + iq, il) ) )
            ENDDO 
            ENDDO 
         ENDIF 
!                                                                       
!------ We want to average rmc_csf                                      
!                                                                       
      ELSE 
!                                                                       
         DO iq = 1, rmc_num (1, ip) * rmc_num (2, ip) 
         dsi (iq) = DBLE (rmc_csf (offsq (ip, is) + iq, 1) * conjg (    &
         rmc_csf (offsq (ip, is) + iq, 1) ) )                           
         ENDDO 
!                                                                       
         IF (rmc_nlots.gt.1) then 
            DO il = 2, rmc_nlots 
            DO iq = 1, rmc_num (1, ip) * rmc_num (2, ip) 
            dsi (iq) = dsi (iq) + DBLE (rmc_csf (offsq (ip, is) + iq,   &
            il) * conjg (rmc_csf (offsq (ip, is) + iq, il) ) )          
            ENDDO 
            ENDDO 
         ENDIF 
!                                                                       
      ENDIF 
      END SUBROUTINE rmc_inten                      
!*****7*****************************************************************
      SUBROUTINE rmc_check_input 
!+                                                                      
!     This routine checks if all input parameters are valid             
!     for a RMC run.                                                    
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE chem_aver_mod
      USE molecule_mod 
      USE rmc_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER i, j, ia 
      LOGICAL laccept 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
!------ atoms present ?                                                 
!                                                                       
      IF (cr_natoms.eq.0) then 
         ier_typ = ER_RMC 
         ier_num = - 3 
         RETURN 
      ENDIF 
!                                                                       
!------ experimental data present ?                                     
!                                                                       
      IF (rmc_nplane.eq.0) then 
         ier_num = - 4 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
!                                                                       
!------ operational mode ok ?                                           
!------ first count selected atoms                                      
!                                                                       
      ia = 0 
      DO i = 0, cr_nscat 
      IF (rmc_allowed (i) ) ia = ia + 1 
      ENDDO 
!                                                                       
!------ ALL modes need at least one selected atom typ                   
!                                                                       
      IF (ia.lt.1) then 
         ier_num = - 10 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
!                                                                       
!------ SWCHEM mode needs at least two selected atoms                   
!                                                                       
      IF (rmc_mode.eq.rmc_mode_swchem.and.ia.lt.2) then 
         ier_num = - 11 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
!                                                                       
!------ SWDISP needs initial displacements > 0.000001                   
!                                                                       
      IF (rmc_mode.eq.rmc_mode_swdisp) then 
         WRITE (output_io, 1050) 
         CALL chem_aver (.false.,.true.) 
         laccept = .false. 
         DO i = 1, cr_nscat 
         laccept = laccept.or. (chem_ave_sig (1, i) .gt.0.000001) .or. (&
         chem_ave_sig (2, i) .gt.0.000001) .or. (chem_ave_sig (3, i)    &
         .gt.0.000001)                                                  
         ENDDO 
         IF (.not.laccept) then 
            ier_num = - 12 
            ier_typ = ER_RMC 
         ENDIF 
      ENDIF 
!                                                                       
!------ Some checks in case we are using molecules                      
!                                                                       
      IF (.not.rmc_sel_atom) then 
         IF (rmc_mode.ne.rmc_mode_shift) then 
            laccept = .true. 
            DO i = 1, mole_num_mole 
            DO j = i, mole_num_mole 
            IF (rmc_allowed (mole_type (i) ) .and.rmc_allowed (         &
            mole_type (j) ) ) then                                      
               laccept = laccept.and. (mole_len (i) .eq.mole_len (j) ) 
            ENDIF 
            ENDDO 
            ENDDO 
            IF (.not.laccept) then 
               ier_num = - 67 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
!                                                                       
         ia = 0 
         DO i = 1, mole_num_mole 
         IF (mole_len (i) .gt.ia) ia = mole_len (i) 
         ENDDO 
         IF (rmc_mode.ne.rmc_mode_shift) ia = 2 * ia 
         IF (ia.gt.rmc_max_atom) then 
            ier_num = - 6 
            ier_typ = ER_RMC 
         ENDIF 
      ENDIF 
!                                                                       
!------ EXTERNAL will get the average structure too ..                  
!                                                                       
      IF (rmc_mode.eq.rmc_mode_external) then 
         WRITE (output_io, 1050) 
         CALL chem_aver (.false.,.true.) 
      ENDIF 
!                                                                       
 1050 FORMAT     (/,' Calculating average structure ... ') 
      END SUBROUTINE rmc_check_input                
!*****7*****************************************************************
      SUBROUTINE rmc_calcskal (ip, wtot, c, cc, ce, se, sk, ba) 
!+                                                                      
!     Calculate scaling factor / background                             
!-                                                                      
      USE discus_config_mod 
      USE rmc_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: ip 
      REAL   , INTENT(IN) :: wtot 
!                                                                       
      REAL(PREC_DP), INTENT(IN) :: c
      REAL(PREC_DP), INTENT(IN) :: cc
      REAL(PREC_DP), INTENT(IN) :: ce
      REAL(PREC_DP), INTENT(IN) :: se
      REAL    , DIMENSION(RMC_MAX_PLANES), INTENT(OUT) :: sk !(rmc_max_planes) 
      REAL    , DIMENSION(RMC_MAX_PLANES), INTENT(OUT) :: ba !(rmc_max_planes) 
!                                                                       
!------ calculate the values                                            
!                                                                       
      IF (rmc_doskal.and.rmc_doback) then 
         sk (ip) = real( (wtot * ce-se * c) / (wtot * cc - c * c) )
         ba (ip) = real( (se-sk (ip) * c) / wtot )
      ELSEIF (rmc_doskal.and..not.rmc_doback) then 
         sk (ip) = real( (ce-rmc_back (ip) * c) / cc )
         ba (ip) = rmc_back (ip) 
      ELSEIF (.not.rmc_doskal.and.rmc_doback) then 
         sk (ip) = rmc_skal (ip) 
         ba (ip) = real( (se-sk (ip) * c) / wtot )
      ELSE 
         sk (ip) = rmc_skal (ip) 
         ba (ip) = rmc_back (ip) 
      ENDIF 
!                                                                       
!------ possible constrains                                             
!                                                                       
      IF (rmc_constrain (ip) .ne.ip) then 
         sk (ip) = sk (rmc_constrain (ip) ) 
         ba (ip) = ba (rmc_constrain (ip) ) 
      ENDIF 
!                                                                       
!------ force the values to be positive                                 
!                                                                       
      IF (sk (ip) .le.0.0) sk (ip) = rmc_skal (ip) 
      IF (ba (ip) .le.0.0) ba (ip) = rmc_back (ip) 
!                                                                       
      END SUBROUTINE rmc_calcskal                   
!*****7*****************************************************************
      SUBROUTINE rmc_stltab (is, ip, lsave) 
!+                                                                      
!     Save or restore SIN(THETA)/LAMBDA for diffuse routines            
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER , INTENT(IN) :: is
      INTEGER , INTENT(IN) :: ip
      LOGICAL , INTENT(IN) :: lsave 
!
      INTEGER i 
!                                                                       
!------ Calculate sin(theta)/lambda table for plan 'ip' and             
!------ symmetry peration 'is' and store it                             
!                                                                       
      IF (lsave) then 
         CALL four_stltab 
         DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
         ristl (offsq (ip, is) + i) = istl (i) 
         ENDDO 
!                                                                       
!------ just restore sin(theta)/lambda table                            
!                                                                       
      ELSE 
         DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
         istl (i) = ristl (offsq (ip, is) + i) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_stltab                     
!*****7*****************************************************************
      SUBROUTINE rmc_formtab (ip, lsave) 
!+                                                                      
!     Save or restore FORMTAB for diffuse routines                      
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: ip
      LOGICAL , INTENT(IN) :: lsave 
!                                                                       
      INTEGER i, j 
!     INTEGER    :: all_status
      LOGICAL bano, blxray, bldbw 
!
!     ALLOCATE(rcfact(0:CFPKT,1:cr_nscat,1),STAT=all_status)
!     IF(all_status /= 0 ) THEN
!        ier_num = -9999
!        RETURN
!     ENDIF
!                                                                       
!------ Calculate formfactor table for plan 'ip' and store it           
!                                                                       
      IF (lsave) then 
         bano = ano 
         blxray = lxray 
         bldbw = ldbw 
         ano = rmc_ano (ip) 
         lxray = rmc_lxray (ip) 
         ldbw = rmc_ldbw (ip) 
!                                                                       
         CALL four_formtab 
         DO i = 1, cr_nscat 
         DO j = 0, CFPKT 
         rcfact (j, i, ip) = cfact (j, i) 
         ENDDO 
         ENDDO 
!                                                                       
         ano = bano 
         lxray = blxray 
         ldbw = bldbw 
!                                                                       
!------ Restore diffuse cfact array for given plane                     
!                                                                       
      ELSE 
         DO i = 1, cr_nscat 
         DO j = 0, CFPKT 
         cfact (j, i) = rcfact (j, i, ip) 
         ENDDO 
         ENDDO 
      ENDIF 
!
!     DEALLOCATE(rcfact)
!                                                                       
      END SUBROUTINE rmc_formtab                    
!*****7*****************************************************************
      SUBROUTINE rmc_layer (is, ip) 
!+                                                                      
!     Set diffuse layer variables for current plane/isym                
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: is
      INTEGER , INTENT(IN) :: ip
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, 3 
      xm (i) = rmc_eck (i, 1, is, ip) 
      uin (i) = rmc_vi (i, 1, is, ip) 
      vin (i) = rmc_vi (i, 2, is, ip) 
      ENDDO 
!                                                                       
      DO i = 1, 2 
      num (i) = rmc_num (i, ip) 
      ENDDO 
      eck(1:3,1:3) = rmc_eck(1:3,1:3,is, ip)
      vi (1:3,1:2) = rmc_vi (1:3,1:2,is, ip)
      vi (1:3,3)   = 0.0
      num(3)       = 1
!                                                                       
      END SUBROUTINE rmc_layer                      
!*****7*****************************************************************
      SUBROUTINE rmc_zero 
!+                                                                      
!     Reset rmc_csf array                                               
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER ip, il 
!                                                                       
      DO ip = 1, rmc_max_sq 
      DO il = 1, rmc_nlots 
      rmc_csf (ip, il) = cmplx (0.0D0, 0.0D0) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE rmc_zero                       
!*****7*****************************************************************
      SUBROUTINE rmc_gsg (mat, nsym, n) 
!+                                                                      
!     Do transformation q = gSg*                                        
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE rmc_mod 
      USE tensors_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: n
      REAL   , DIMENSION(4,4,n), INTENT(INOUT) :: mat
      INTEGER, INTENT(IN) :: nsym
!                                                                       
      INTEGER i, j, k
      REAL a (3, 3), b (3, 3) 
!                                                                       
      DO i = 1, nsym 
      DO j = 1, 3 
      DO k = 1, 3 
      a (j, k) = mat (j, k, i) 
      ENDDO 
      ENDDO 
!                                                                       
      CALL matmulx (b, a, cr_rten) 
      CALL matmulx (a, cr_gten, b) 
      DO j = 1, 3 
      DO k = 1, 3 
      mat (j, k, i) = a (j, k) 
      ENDDO 
      mat (j, 4, i) = mat (j, 4, i) 
      mat (4, j, i) = 0.0 
      ENDDO 
      mat (4, 4, i) = 1.0 
      ENDDO 
!                                                                       
      END SUBROUTINE rmc_gsg                        
!*****7*****************************************************************
      SUBROUTINE rmc_trans (v, mat, nsym, n) 
!+                                                                      
!     Compute symmetrically equivalent vectors of 'v' using 'mat'       
!-                                                                      
      USE discus_config_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: n
      REAL   , DIMENSION(4,  n), INTENT(INOUT) :: v
      REAL   , DIMENSION(4,4,n), INTENT(INOUT) :: mat
      INTEGER, INTENT(IN) :: nsym
!                                                                       
!      REAL v (4, n), mat (4, 4, n) 
!                                                                       
      INTEGER i, j, isym 
!                                                                       
      DO isym = 2, nsym 
      DO i = 1, 3 
      v (i, isym) = 0.0 
      DO j = 1, 3 
      v (i, isym) = v (i, isym) + mat (i, j, isym) * v (j, 1) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE rmc_trans                      
!*****7*****************************************************************
      SUBROUTINE rmc_fcalc (ip, is, natoms, i_new, p_new, isel) 
!+                                                                      
!     This routine calculates the change in structure factor            
!     cause by the RMC move.                                            
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN)  :: ip
      INTEGER, INTENT(IN)  :: is
      INTEGER, INTENT(IN)  :: natoms
      REAL   , DIMENSION(3, RMC_MAX_ATOM), INTENT(IN) :: p_new !(3, rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(IN) :: i_new !(rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(IN) :: isel  !(rmc_max_atom) 
       
!                                                                       
      REAL p_old (3, rmc_max_atom) 
      REAL off (3) 
      INTEGER i_old (rmc_max_atom) 
      INTEGER i, il
!                                                                       
!     LOGICAL rmc_inlot 
!                                                                       
      DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
      DO il = 1, rmc_nlots 
      rmc_csf_new (offsq (ip, is) + i, il) = rmc_csf (offsq (ip, is)    &
      + i, il)                                                          
      ENDDO 
      ENDDO 
!                                                                       
      DO i = 1, natoms 
      p_old (1, i) = cr_pos (1, isel (i) ) 
      p_old (2, i) = cr_pos (2, isel (i) ) 
      p_old (3, i) = cr_pos (3, isel (i) ) 
      i_old (i) = cr_iscat (isel (i) ) 
      ENDDO 
!                                                                       
      CALL rmc_formtab (ip, .false.) 
      CALL rmc_layer (is, ip) 
      CALL rmc_stltab (is, ip, .false.) 
!                                                                       
      DO il = 1, rmc_nlots 
      DO i = 1, natoms 
      IF (rmc_inlot (isel (i), il, off) ) then 
         CALL rmc_strucf (i_old (i), p_old, i, is, ip, il, off, .false.) 
         CALL rmc_strucf (i_new (i), p_new, i, is, ip, il, off, .true.) 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE rmc_fcalc                      
!*****7*****************************************************************
      SUBROUTINE rmc_strucf (iscat, pos, isite, is, ip, il, off, lplus) 
!+                                                                      
!     Here is the real calculation :-)                                  
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE four_strucf_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat
      REAL   , DIMENSION(3, RMC_MAX_ATOM), INTENT(IN) :: pos !(3, rmc_max_atom)
      REAL   , DIMENSION(3)              , INTENT(IN) :: off !(3) 
      INTEGER, INTENT(IN) :: isite
      INTEGER, INTENT(IN) :: is
      INTEGER, INTENT(IN) :: ip
      INTEGER, INTENT(IN) :: il
      LOGICAL, INTENT(IN) :: lplus 
!                                                                       
      INTEGER i 
!                                                                       
      IF (iscat.eq.0) return 
!                                                                       
      nxat = 1 
      xat (nxat, 1) = pos (1, isite) - off (1) 
      xat (nxat, 2) = pos (2, isite) - off (2) 
      xat (nxat, 3) = pos (3, isite) - off (3) 
!                                                                       
      CALL four_strucf (iscat, .true.) 
!                                                                       
      IF (lplus) then 
         DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
         rmc_csf_new (offsq (ip, is) + i, il) = rmc_csf_new (offsq (ip, &
         is) + i, il) + tcsf (i)                                        
         ENDDO 
      ELSE 
         DO i = 1, rmc_num (1, ip) * rmc_num (2, ip) 
         rmc_csf_new (offsq (ip, is) + i, il) = rmc_csf_new (offsq (ip, &
         is) + i, il) - tcsf (i)                                        
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_strucf                     
!*****7*****************************************************************
      LOGICAL function rmc_inlot (ia, il, off) 
!+                                                                      
!     Checks if given atom 'ia' is in given lot 'il'                    
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE diffuse_mod 
!     USE modify_mod
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
      INTEGER              , INTENT(IN)  :: ia 
      INTEGER              , INTENT(IN)  :: il 
      REAL   , DIMENSION(3), INTENT(OUT) :: off !(3) 
!                                                                       
      REAL xtest (3), x0 
      INTEGER cr_end 
      INTEGER iz (3), izmin, izmax 
      INTEGER is, i 
!     
      rmc_inlot = .false.
      cr_end = cr_ncatoms * cr_icc (1) * cr_icc (2) * cr_icc (3) + 1
!                                                                       
!------ We are not using lots                                           
!                                                                       
      IF (rmc_ilots.eq.LOT_OFF) then 
         rmc_inlot = .true. 
         DO i = 1, 3 
         off (i) = 0.0 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (ilots.eq.LOT_BOX) then 
         IF (ia.lt.cr_end) then 
            CALL indextocell (ia, iz, is) 
         ELSE 
            DO i = 1, 3 
            iz (i) = int (cr_pos (i, ia) + cr_dim0 (i, 1) ) - 1 
            ENDDO 
         ENDIF 
!                                                                       
         rmc_inlot = .true. 
!                                                                       
         DO i = 1, 3 
         izmin = rmc_lots_orig (i, il) 
         izmax = rmc_lots_orig (i, il) + ls_xyz (i) - 1 
         off (i) = cr_dim0 (i, 1) + rmc_lots_orig (i, il) - 1 
         rmc_inlot = rmc_inlot.and. (iz (i) .ge.izmin) .and. (iz (i)    &
         .le.izmax)                                                     
!                                                                       
         IF (.not.rmc_inlot.and. (izmax.gt.cr_icc (i) ) ) then 
            izmin = 1 
            izmax = izmax - cr_icc (i) 
            off (i) = off (i) - float (cr_icc (i) ) 
            rmc_inlot = rmc_inlot.and. (iz (i) .ge.izmin) .and. (iz (i) &
            .le.izmax)                                                  
         ENDIF 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (rmc_ilots.eq.LOT_ELI) then 
         IF (ia.lt.cr_end) then 
            CALL indextocell (ia, iz, is) 
         ELSE 
            DO i = 1, 3 
            iz (i) = int (cr_pos (i, ia) + cr_dim0 (i, 1) ) - 1 
            ENDDO 
         ENDIF 
!                                                                       
         rmc_inlot = .true. 
!                                                                       
         DO i = 1, 3 
         x0 = float (ls_xyz (i) ) / 2.0 
         izmin = rmc_lots_orig (i, il) 
         izmax = rmc_lots_orig (i, il) + ls_xyz (i) - 1 
         off (i) = cr_dim0 (i, 1) + rmc_lots_orig (i, il) - 1 
         xtest (i) = (iz (i) - izmin - x0 + 0.5) **2 / x0**2 
         rmc_inlot = rmc_inlot.and. (iz (i) .ge.izmin) .and. (iz (i)    &
         .le.izmax)                                                     
!                                                                       
         IF (.not.rmc_inlot.and. (izmax.gt.cr_icc (i) ) ) then 
            izmin = 1 
            izmax = izmax - cr_icc (i) 
            off (i) = off (i) - float (cr_icc (i) ) 
            rmc_inlot = rmc_inlot.and. (iz (i) .ge.izmin) .and. (iz (i) &
            .le.izmax)                                                  
         ENDIF 
         ENDDO 
!                                                                       
         rmc_inlot = rmc_inlot.and. ( (xtest (1) + xtest (2) + xtest (3)&
         ) .le.1.0)                                                     
!                                                                       
      ENDIF 
      END FUNCTION rmc_inlot                        
!*****7*****************************************************************
      SUBROUTINE rmc_makemove (isym, natoms, i_new, p_new, isel, imol) 
!+                                                                      
!     make accepted move                                                
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE rmc_mod 
      IMPLICIT none 
!                                                                       
       
      INTEGER, DIMENSION(RMC_MAX_PLANES), INTENT(IN) :: isym !(rmc_max_planes) 
      INTEGER                           , INTENT(IN) :: natoms 
!                                                                       
      INTEGER, DIMENSION(  RMC_MAX_ATOM), INTENT(IN) :: i_new !(rmc_max_atom) 
      REAL   , DIMENSION(3,RMC_MAX_ATOM), INTENT(IN) :: p_new !(3, rmc_max_atom) 
      INTEGER, DIMENSION(  RMC_MAX_ATOM), INTENT(IN) :: isel  !(rmc_max_atom) 
      INTEGER, DIMENSION(  RMC_MAX_ATOM), INTENT(IN) :: imol  !(rmc_max_atom) 
!                                                                       
      INTEGER i, j, ip, iq, is, il 
!                                                                       
      DO ip = 1, rmc_nplane 
      DO is = 1, isym (ip) 
      DO iq = 1, rmc_num (1, ip) * rmc_num (2, ip) 
      DO il = 1, rmc_nlots 
      rmc_csf (offsq (ip, is) + iq, il) = rmc_csf_new (offsq (ip, is)   &
      + iq, il)                                                         
      ENDDO 
      ENDDO 
      ENDDO 
      ENDDO 
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
      END SUBROUTINE rmc_makemove                   
!*****7*****************************************************************
      SUBROUTINE rmc_check_dist (laccept, iatom, inumber, i_new, p_new) 
!+                                                                      
!     check distances                                                   
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE chem_mod 
      USE modify_mod
      USE rmc_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
      LOGICAL, INTENT(OUT) :: laccept 
      INTEGER, INTENT(IN)  :: iatom 
      INTEGER, INTENT(IN)  :: inumber 
      INTEGER, DIMENSION(  RMC_MAX_ATOM), INTENT(IN)  :: i_new !(rmc_max_atom)
      REAL   , DIMENSION(3,RMC_MAX_ATOM), INTENT(IN)  :: p_new !(3, rmc_max_atom) 
!                                                                       
      INTEGER i
      REAL pos (3), werte(1) 
!                                                                       
      werte = - 1 
      DO i = 1, 3 
      pos (i) = p_new (i, iatom) 
      ENDDO 
!                                                                       
!------ check distances to neighbouring atoms                           
!                                                                       
      IF (i_new (iatom) .eq.0) return 
!                                                                       
      CALL do_find_env (1, werte, 1, pos, 0.1, rmc_mindist_max,         &
      chem_quick, chem_period)                                          
      DO i = 1, atom_env (0) 
      IF (cr_iscat (atom_env (i) ) .ne.0.and.laccept) then 
         IF(inumber /= atom_env(i)) THEN ! Do not check distance to original place!
         laccept = (res_para (i) .ge.rmc_mindist (i_new (iatom),        &
         cr_iscat (atom_env (i) ) ) )                                   
            ENDIF
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE rmc_check_dist                 
!*****7*****************************************************************
      SUBROUTINE rmc_select (imode, isel, iz1, iz2, is1, is2, &
                             MAXALLOWED, local_rmc_allowed) 
!+                                                                      
!     Select two atom sites (global, local, ...)                        
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE conn_mod
      USE rmc_mod 
      USE random_mod
      IMPLICIT none 
!                                                                       
      INTEGER                         , INTENT(IN)    :: imode
      INTEGER, DIMENSION(RMC_MAX_ATOM), INTENT(INOUT) :: isel !(rmc_max_atom) 
      INTEGER, DIMENSION(3)           , INTENT(INOUT) :: iz1 !(3)
      INTEGER, DIMENSION(3)           , INTENT(INOUT) :: iz2 !(3)
      INTEGER                         , INTENT(INOUT) :: is1
      INTEGER                         , INTENT(INOUT) :: is2
      INTEGER                         , INTENT(IN   ) :: MAXALLOWED
      LOGICAL, DIMENSION(0:MAXALLOWED), INTENT(IN   ) :: local_rmc_allowed
!
      INTEGER, PARAMETER :: MAXW = 25
      INTEGER, DIMENSION(:), ALLOCATABLE   :: c_list ! Result of connectivity search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs ! Result of connectivity search
      INTEGER            :: natoms ! Number of actual actoms in connectivity list
!                                                                       
      INTEGER i,j
      REAL ran1 
!                                                                       
      j = 0
   10 CONTINUE 
      j = j + 1
      isel (1) = int (ran1 (idum) * cr_natoms) + 1 
      IF (isel (1) .gt.cr_natoms.or.isel (1) .lt.1) goto 10 
      IF(.NOT.local_rmc_allowed(cr_iscat(isel(1)))) GOTO 10
!
      IF (imode == rmc_local_conn) THEN   ! Choose second atom from connectivity
         CALL get_connectivity_list(isel(1), cr_iscat(isel(1)), 1, MAXW, c_list, c_offs, natoms)
         IF(natoms == 0) GOTO 10
         isel(2) = c_list(INT(ran1(idum)*natoms) + 1)
      ELSE
!                                                                       
   20 CONTINUE 
      isel (2) = int (ran1 (idum) * cr_natoms) + 1 
      IF (isel (2) .gt.cr_natoms.or.isel (2) .lt.1) goto 20 
!
      ENDIF
!                                                                       
      CALL indextocell (isel (1), iz1, is1) 
      CALL indextocell (isel (2), iz2, is2) 
!                                                                       
!------ Force to be same site                                           
!                                                                       
      IF (imode.eq.rmc_local_site) then 
         is2 = is1 
         CALL celltoindex (iz2, is2, isel (2) ) 
!                                                                       
!------ only local +- 1 unit cell for second atom                       
!                                                                       
      ELSEIF (imode.eq.rmc_local_loc) then 
         DO i = 1, 3 
            IF( cr_icc(i) > 1 ) THEN
               iz1 (i) = max (2, iz1 (i) ) 
               iz1 (i) = min (cr_icc (i) - 1, iz1 (i) ) 
               iz2 (i) = iz1 (i) + nint (1.0 - 2.0 * ran1 (idum) ) 
            ELSE
               iz1 (i) = 1
               iz2 (i) = 1
            ENDIF
         ENDDO 
         is2 = MIN(int (ran1 (idum) * cr_ncatoms) + 1 ,cr_ncatoms)
         CALL celltoindex (iz1, is1, isel (1) ) 
         CALL celltoindex (iz2, is2, isel (2) ) 
!                                                                       
!------ only local +- 1 unit cell for second atom, same site            
!                                                                       
      ELSEIF (imode.eq.rmc_local_locsite) then 
         DO i = 1, 3 
            IF( cr_icc(i) > 1 ) THEN
               iz1 (i) = max (2, iz1 (i) ) 
               iz1 (i) = min (cr_icc (i) - 1, iz1 (i) ) 
               iz2 (i) = iz1 (i) + nint (1.0 - 2.0 * ran1 (idum) ) 
            ELSE
               iz1 (i) = 1
               iz2 (i) = 1
            ENDIF
         ENDDO 
         is2 = is1 
         CALL celltoindex (iz1, is1, isel (1) ) 
         CALL celltoindex (iz2, is2, isel (2) ) 
      ENDIF 
!                                                                       
      END SUBROUTINE rmc_select                     
!*****7*****************************************************************
      SUBROUTINE rmc_genmove (laccept, natoms, p_new, i_new, isel) 
!+                                                                      
!     generate RMC move (for switch and relax. mode)                    
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE diffuse_mod 
      USE extrmc
      USE rmc_mod 
      USE random_mod
!
      USE debug_mod 
      USE errlist_mod
!                                                                       
      USE prompt_mod 
      USE random_mod
      IMPLICIT none 
!
      LOGICAL, INTENT(OUT) :: laccept 
      INTEGER, INTENT(OUT) :: natoms 
      REAL   , DIMENSION(3, RMC_MAX_ATOM), INTENT(OUT) :: p_new !(3, rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(OUT) :: i_new !(rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(OUT) :: isel  !(rmc_max_atom) 
!                                                                       
      INTEGER iz1 (3), iz2 (3) 
      INTEGER i, j, is1, is2, il 
      REAL  :: value
      REAL disp1, disp2 
      REAL ran1, gasdev 
      REAL dummy (3) 
      LOGICAL lflag  !, rmc_inlot 
!
!------ Randomly select a move
!
      value = ran1(idum)
      i     = 1
      DO WHILE (value>rmc_move_cprob(i) .and. i<RMC_N_MOVE)
         i = i + 1
      ENDDO
      rmc_mode  = i
      rmc_local = rmc_move_local(i)
!                                                                       
!------ Mode 'swchem': only switch atoms => pnew=pold                   
!                                                                       
      IF (rmc_mode.eq.rmc_mode_swchem) then 
         natoms = 2 
         CALL rmc_select (rmc_local, isel, iz1, iz2, is1, is2,  &
                          cr_nscat,rmc_allowed) 
         laccept = rmc_allowed (cr_iscat (isel (1) ) ) .and.rmc_allowed &
         (cr_iscat (isel (2) ) ) .and.cr_iscat (isel (1) ) .ne.cr_iscat &
         (isel (2) )                                                    
         IF (laccept) then 
            i_new (1) = cr_iscat (isel (2) ) 
            i_new (2) = cr_iscat (isel (1) ) 
            DO j = 1, 3 
            p_new (j, 1) = cr_pos (j, isel (1) ) 
            p_new (j, 2) = cr_pos (j, isel (2) ) 
            ENDDO 
         ENDIF 
!                                                                       
!------ Mode 'shift': shift atoms (gaussian distribution)               
!                                                                       
      ELSEIF (rmc_mode.eq.rmc_mode_shift) then 
         natoms = 1 
   10    CONTINUE 
         isel (1) = int (ran1 (idum) * cr_natoms) + 1 
         IF (isel (1) .gt.cr_natoms.or.isel (1) .lt.1) goto 10 
         laccept = rmc_allowed (cr_iscat (isel (1) ) ) 
         IF (laccept) then 
            i_new (1) = cr_iscat (isel (1) ) 
            DO i = 1, 3 
            p_new (i, 1) = cr_pos (i, isel (1) ) + gasdev (rmc_maxmove (&
            i, i_new (1) ) )                                            
!              IF(chem_period(i)) THEN ! Apply modulo of crystal dimensions
!                 IF(p_new(i,1) > cr_dim0(i,2) + 0.75 .or.  &
!                    p_new(i,1) < cr_dim0(i,1) - 0.75     ) THEN
!                    p_new(i,1) = MOD( (p_new(i,1)-cr_dim0(i,1))+    &
!                                       (cr_dim0(i,2)-cr_dim0(i,1)), &
!                                 (cr_dim0(i,2)-cr_dim0(i,1))+cr_dim0(i,1))
!                 ENDIF
!              ENDIF
            ENDDO 
         ENDIF 
!                                                                       
!------ Mode 'swdisp': switch displacements of two atoms                
!                                                                       
      ELSEIF (rmc_mode.eq.rmc_mode_swdisp) then 
         natoms = 2 
         CALL rmc_select (rmc_local, isel, iz1, iz2, is1, is2, &
                          cr_nscat,rmc_allowed) 
         laccept = rmc_allowed (cr_iscat (isel (1) ) ) .and.rmc_allowed &
         (cr_iscat (isel (2) ) )                                        
         IF (laccept) then 
            i_new (1) = cr_iscat (isel (1) ) 
            i_new (2) = cr_iscat (isel (2) ) 
            DO j = 1, 3 
            disp1 = cr_pos (j, isel (1) ) - chem_ave_pos (j, is1)       &
            - float (iz1 (j) - 1) - cr_dim0 (j, 1)                      
            disp2 = cr_pos (j, isel (2) ) - chem_ave_pos (j, is2)       &
            - float (iz2 (j) - 1) - cr_dim0 (j, 1)                      
            p_new (j, 1) = cr_pos (j, isel (1) ) - disp1 + disp2 
            p_new (j, 2) = cr_pos (j, isel (2) ) - disp2 + disp1 
            ENDDO 
         ENDIF 
!                                                                       
!------ Mode 'external': Moves defined by external subroutine           
!                                                                       
      ELSEIF (rmc_mode.eq.rmc_mode_external) then 
         CALL rmc_genmove_ext (laccept, natoms, p_new, i_new, isel) 
      ENDIF 
!                                                                       
!------ Check if atoms are in any defined 'lots'                        
!                                                                       
      IF (laccept.and. (rmc_ilots.ne.LOT_OFF) ) then 
         DO i = 1, natoms 
         lflag = .false. 
         DO il = 1, rmc_nlots 
         lflag = lflag.or. (rmc_inlot (isel (i), il, dummy) ) 
         ENDDO 
         laccept = laccept.and.lflag 
         ENDDO 
      ENDIF 
!                                                                       
!------ Check distances                                                 
!                                                                       
      DO i = 1, natoms 
      IF (laccept) then 
         CALL rmc_check_dist (laccept, i, isel(i),i_new, p_new) 
      ENDIF 
      ENDDO 
!                                                                       
!------ Uncomment these lines to debug generated moves                  
!                                                                       
      IF (laccept.and.dbg) then 
         DO j = 1, natoms 
         WRITE (output_io, 11) j, isel (j), cr_at_lis (cr_iscat (isel ( &
         j) ) ), cr_iscat (isel (j) ), (cr_pos (i, isel (j) ), i = 1, 3)&
         , isel (j), cr_at_lis (i_new (j) ), i_new (j), (p_new (i, j),  &
         i = 1, 3), ( (p_new (i, j) - cr_pos (i, isel (j) ) ), i = 1, 3)
         ENDDO 
      ENDIF 
!                                                                       
   11 FORMAT(' DBG: # ',I3,' old: atom ',I6,' (',A4,',#',I3,            &
     &       ') : ',3(F9.3,1X),/,12X,'new: atom ',I6,' (',A4,           &
     &       ',#',I3,') : ',3(F9.3,1X),/,35X,'shift : ',                &
     &       3(F9.3,1X),/,72('-'))                                      
!                                                                       
      END SUBROUTINE rmc_genmove                    
!********************************************************************** 
      SUBROUTINE rmc_genmove_mol (laccept, natoms, p_new, i_new, isel,  &
      imol)                                                             
!+                                                                      
!     Generate RMC move (for switch and relax. mode)                    
!     Version for molecules ..                                          
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_mod 
      USE celltoindex_mod
      USE diffuse_mod 
!     USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      USE extrmc
      USE random_mod
!
      USE debug_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!
      LOGICAL, INTENT(OUT) :: laccept 
      INTEGER, INTENT(OUT) :: natoms 
      REAL   , DIMENSION(3, RMC_MAX_ATOM), INTENT(OUT) :: p_new !(3, rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(OUT) :: i_new !(rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(OUT) :: isel  !(rmc_max_atom) 
      INTEGER, DIMENSION(   RMC_MAX_ATOM), INTENT(OUT) :: imol !(rmc_max_atom) 
!                                                                       
      INTEGER iz1 (3), iz2 (3) 
      INTEGER i, j, k, is1, is2, il, i0, j0 
      REAL disp1 (3), disp2 (3) 
      REAL ran1, gasdev 
      REAL dummy (3) 
      LOGICAL lflag  !, rmc_inlot 
!                                                                       
      natoms = 0 
!                                                                       
!------ Mode 'swchem': only switch atoms => pnew=pold                   
!                                                                       
      IF (rmc_mode.eq.rmc_mode_swchem) then 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         imol (2) = int (ran1 (idum) * mole_num_mole) + 1 
         laccept = rmc_allowed (mole_type (imol (1) ) )                 &
         .and.rmc_allowed (mole_type (imol (2) ) ) .and.mole_type (imol &
         (1) ) .ne.mole_type (imol (2) )                                
!                                                                       
         IF (laccept) then 
            natoms = 2 * mole_len (imol (1) ) 
            i0 = mole_cont (mole_off (imol (1) ) + 1) 
            j0 = mole_cont (mole_off (imol (2) ) + 1) 
            DO i = 1, mole_len (imol (1) ) 
            j = i + mole_len (imol (1) ) 
            isel (i) = mole_cont (mole_off (imol (1) ) + i) 
            isel (j) = mole_cont (mole_off (imol (2) ) + i) 
            i_new (i) = cr_iscat (isel (j) ) 
            i_new (j) = cr_iscat (isel (i) ) 
            DO k = 1, 3 
            p_new (k, i) = cr_pos (k, isel (j) ) - cr_pos (k, j0)       &
            + cr_pos (k, i0)                                            
            p_new (k, j) = cr_pos (k, isel (i) ) - cr_pos (k, i0)       &
            + cr_pos (k, j0)                                            
            ENDDO 
            ENDDO 
         ENDIF 
!                                                                       
!------ Mode 'shift': shift atoms (gaussian distribution)               
!                                                                       
      ELSEIF (rmc_mode.eq.rmc_mode_shift) then 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         laccept = rmc_allowed (mole_type (imol (1) ) ) 
!                                                                       
         IF (laccept) then 
            natoms = mole_len (imol (1) ) 
            DO j = 1, 3 
            disp1 (j) = gasdev (rmc_maxmove (j, mole_type (imol (1) ) ) &
            )                                                           
            ENDDO 
!                                                                       
            DO i = 1, mole_len (imol (1) ) 
            isel (i) = mole_cont (mole_off (imol (1) ) + i) 
            i_new (i) = cr_iscat (isel (i) ) 
            DO j = 1, 3 
            p_new (j, i) = cr_pos (j, isel (i) ) + disp1 (j) 
            ENDDO 
            ENDDO 
         ENDIF 
!                                                                       
!------ Mode 'swdisp': switch displacements of two atoms                
!                                                                       
      ELSEIF (rmc_mode.eq.rmc_mode_swdisp) then 
         imol (1) = int (ran1 (idum) * mole_num_mole) + 1 
         imol (2) = int (ran1 (idum) * mole_num_mole) + 1 
         laccept = rmc_allowed (mole_type (imol (1) ) )                 &
         .and.rmc_allowed (mole_type (imol (2) ) ) .and.mole_type (imol &
         (1) ) .ne.mole_type (imol (2) )                                
!                                                                       
         IF (laccept) then 
            natoms = 2 * mole_len (imol (1) ) 
            CALL indextocell (mole_cont (mole_off (imol (1) ) + 1),     &
            iz1, is1)                                                   
            CALL indextocell (mole_cont (mole_off (imol (2) ) + 1),     &
            iz2, is2)                                                   
            DO j = 1, 3 
            disp1 (j) = cr_pos (j, mole_cont (mole_off (imol (1) )      &
            + 1) ) - chem_ave_pos (j, is1) - float (iz1 (j) - 1)        &
            - cr_dim0 (j, 1)                                            
            disp2 (j) = cr_pos (j, mole_cont (mole_off (imol (2) )      &
            + 1) ) - chem_ave_pos (j, is2) - float (iz2 (j) - 1)        &
            - cr_dim0 (j, 1)                                            
            ENDDO 
!                                                                       
            DO i = 1, mole_len (imol (1) ) 
            j = i + mole_len (imol (1) ) 
            isel (i) = mole_cont (mole_off (imol (1) ) + i) 
            isel (j) = mole_cont (mole_off (imol (2) ) + i) 
            i_new (i) = cr_iscat (isel (i) ) 
            i_new (j) = cr_iscat (isel (j) ) 
            DO k = 1, 3 
            p_new (k, i) = cr_pos (k, isel (i) ) - disp1 (k) + disp2 (k) 
            p_new (k, j) = cr_pos (k, isel (j) ) - disp2 (k) + disp1 (k) 
            ENDDO 
            ENDDO 
         ENDIF 
!                                                                       
!------ Mode 'external': Moves defined by external subroutine           
!                                                                       
      ELSEIF (rmc_mode.eq.rmc_mode_external) then 
         CALL rmc_genmove_ext (laccept, natoms, p_new, i_new, isel) 
      ENDIF 
!                                                                       
!------ Check if atoms are in any defined 'lots'                        
!                                                                       
      IF (laccept.and. (rmc_ilots.ne.LOT_OFF) ) then 
         DO i = 1, natoms 
         lflag = .false. 
         DO il = 1, rmc_nlots 
         lflag = lflag.or. (rmc_inlot (isel (i), il, dummy) ) 
         ENDDO 
         laccept = laccept.and.lflag 
         ENDDO 
      ENDIF 
!                                                                       
!------ Check distances                                                 
!                                                                       
      DO i = 1, natoms 
      IF (laccept) then 
         CALL rmc_check_dist (laccept, i, isel(i),i_new, p_new) 
      ENDIF 
      ENDDO 
!                                                                       
!------ Uncomment these lines to debug generated moves                  
!                                                                       
      IF (laccept.and.dbg) then 
         WRITE (output_io, 11) imol 
         DO j = 1, natoms 
         WRITE (output_io, 12) j, isel (j), cr_at_lis (cr_iscat (isel ( &
         j) ) ), cr_iscat (isel (j) ), (cr_pos (i, isel (j) ), i = 1, 3)&
         , isel (j), cr_at_lis (i_new (j) ), i_new (j), (p_new (i, j),  &
         i = 1, 3), ( (p_new (i, j) - cr_pos (i, isel (j) ) ), i = 1, 3)
         ENDDO 
      ENDIF 
!                                                                       
   11 FORMAT(/,' DBG: molecules # ',10I3) 
   12 FORMAT(' DBG: # ',I3,' old: atom ',I6,' (',A4,',#',I3,            &
     &       ') : ',3(F9.3,1X),/,12X,'new: atom ',I6,' (',A4,           &
     &       ',#',I3,') : ',3(F9.3,1X),/,35X,'shift : ',                &
     &       3(F9.3,1X),/,72('-'))                                      
!                                                                       
      END SUBROUTINE rmc_genmove_mol                
!********************************************************************** 
END MODULE rmc_sup_mod
