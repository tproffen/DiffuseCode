MODULE discus_allocate_appl_mod
!-
!     Contains data and routines to allocate application dependent arrays
!
!     SUBROUTINE discus_show_config    ! Shows the current configuration
!     SUBROUTINE discus_alloc_default  ! Allocate default array sizes, mostly used at program start
!     SUBROUTINE alloc_constraint ( n_constr)  ! Allocate the number of constraints
!     SUBROUTINE alloc_population ( n_pop, n_dimx)  ! Allocate the number of members, parameters
!+
  USE allocate_generic
  USE discus_config_mod
  USE errlist_mod 
!
  CONTAINS
!
    SUBROUTINE discus_do_allocate_appl(zeile,lcomm)
!
    USE ber_params_mod
    USE get_params_mod
USE precision_mod
USE str_comp_mod
    IMPLICIT NONE
!
!
    CHARACTER (LEN=*), INTENT(IN)            :: zeile
    INTEGER          , INTENT(INOUT)         :: lcomm
!
    INTEGER , PARAMETER                      :: MAXW=10
    CHARACTER (LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW)  :: cpara
    INTEGER             , DIMENSION(1:MAXW)  :: lpara
    REAL(KIND=PREC_DP)  , DIMENSION(1:MAXW)  :: werte
    INTEGER                                  :: ianz
!
!
    CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
    IF (ier_num.eq.0) then 
       IF ( ianz == 0 ) THEN
          CALL discus_show_config
       ELSE IF (str_comp (cpara (1) , 'default'   , 1, lpara (1) , 7) )  THEN
          CALL discus_alloc_default
       ELSE IF (str_comp (cpara (1) , 'crystal'   , 1, lpara (1) , 5) )  THEN
          IF ( ianz == 3 ) THEN
             CALL del_params (1, ianz, cpara, lpara, maxw) 
             CALL ber_params (ianz, cpara, lpara, werte, maxw) 
             IF (ier_num.eq.0) then 
               CALL alloc_crystal_scat ( NINT (werte(1)))
               CALL alloc_crystal_nmax ( NINT (werte(2)))
             ELSE 
                ier_num = - 6 
                ier_typ = ER_COMM 
             ENDIF 
          ELSE 
             ier_num = - 6 
             ier_typ = ER_COMM 
          ENDIF 
       ELSE IF (str_comp (cpara (1) , 'show', 1, lpara (1) , 4) )  THEN
          CALL discus_show_config
       ELSE 
          ier_num = - 6 
          ier_typ = ER_COMM 
       ENDIF 
    ENDIF 
!
    END SUBROUTINE discus_do_allocate_appl
!
    SUBROUTINE discus_do_deallocate_appl(zeile,lcomm)
!
       USE get_params_mod
USE precision_mod
USE str_comp_mod
       IMPLICIT NONE
!
!
       CHARACTER (LEN=*), INTENT(IN)            :: zeile
       INTEGER          , INTENT(INOUT)         :: lcomm
!
       INTEGER , PARAMETER                      :: MAXW=10
       CHARACTER (LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW)  :: cpara
       INTEGER             , DIMENSION(1:MAXW)  :: lpara
       INTEGER                                  :: ianz
!
!
       CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
       IF (ier_num.eq.0) then
            IF ( ianz == 0 ) THEN
               ier_num = - 6
               ier_typ = ER_COMM
            ELSE IF (str_comp (cpara (1) , 'all'       , 3, lpara (1) , 3) )  THEN
               CALL discus_alloc_default
!
            ELSE IF (str_comp (cpara (1) , 'chemistry' , 3, lpara (1) , 9) )  THEN
               CALL dealloc_chemistry
            ELSE IF (str_comp (cpara (1) , 'domain'    , 3, lpara (1) , 6) )  THEN
               CALL dealloc_domain
            ELSE IF (str_comp (cpara (1) , 'mmc'       , 3, lpara (1) , 3) )  THEN
               CALL dealloc_mmc
            ELSE IF (str_comp (cpara (1) , 'stack'     , 3, lpara (1) , 5) )  THEN
               CALL dealloc_stack
!
            ELSE IF (str_comp (cpara (1) , 'crystal'   , 3, lpara (1) , 7) )  THEN
               CALL dealloc_crystal
            ELSE IF (str_comp (cpara (1) , 'debye'     , 3, lpara (1) , 5) )  THEN
               CALL dealloc_debye
            ELSE IF (str_comp (cpara (1) , 'diffuse'   , 3, lpara (1) , 7) )  THEN
               CALL dealloc_diffuse
            ELSE IF (str_comp (cpara (1) , 'molecule'  , 3, lpara (1) , 8) )  THEN
               CALL dealloc_molecule
            ELSE IF (str_comp (cpara (1) , 'phase'     , 3, lpara (1) , 5) )  THEN
               CALL dealloc_phases
            ELSE IF (str_comp (cpara (1) , 'pdf'       , 3, lpara (1) , 3) )  THEN
               CALL dealloc_pdf
            ELSE IF (str_comp (cpara (1) , 'plot'      , 3, lpara (1) , 4) )  THEN
               CALL dealloc_plot
            ELSE IF (str_comp (cpara (1) , 'powder'    , 3, lpara (1) , 6) )  THEN
               CALL dealloc_powder
               CALL dealloc_powder_nmax
            ELSE IF (str_comp (cpara (1) , 'rmc'       , 3, lpara (1) , 3) )  THEN
               CALL dealloc_rmc    
            ELSE IF (str_comp (cpara (1) , 'save'      , 3, lpara (1) , 4) )  THEN
               CALL dealloc_save    
            ELSE IF (str_comp (cpara (1) , 'shear'     , 3, lpara (1) , 5) )  THEN
               CALL dealloc_shear
            ELSE IF (str_comp (cpara (1) , 'super'     , 3, lpara (1) , 5) )  THEN
               CALL dealloc_super
            ELSE IF (str_comp (cpara (1) , 'surface'   , 3, lpara (1) , 7) )  THEN
               CALL dealloc_surf
            ELSE IF (str_comp (cpara (1) , 'symmetry'  , 3, lpara (1) , 8) )  THEN
               CALL dealloc_symmetry
            ELSE IF (str_comp (cpara (1) , 'transform' , 3, lpara (1) , 9) )  THEN
               CALL dealloc_transfrm
            ELSE IF (str_comp (cpara (1) , 'waves'     , 3, lpara (1) , 5) )  THEN
               CALL dealloc_waves
            ELSE
               ier_num = - 6
               ier_typ = ER_COMM
            ENDIF
       ENDIF
!
    END SUBROUTINE discus_do_deallocate_appl
!
!
    SUBROUTINE discus_show_config
!
       USE chem_mod
       USE atom_env_mod
       USE debye_mod
       USE diffuse_mod
       USE domain_mod
       USE micro_mod
       USE mmc_mod
       USE pdf_mod
       USE discus_plot_mod
       USE powder_mod
       USE rmc_mod
       USE discus_save_mod
       USE shear_mod
       USE surface_mod
       USE stack_mod
       USE symm_mod
       USE transfrm_mod
       USE waves_mod
!
       USE prompt_mod
       IMPLICIT NONE
!
!
       WRITE (output_io, 1000)
!
       WRITE (output_io, 1005) MAXSCAT
       WRITE (output_io, 1010) NMAX
       WRITE (output_io, 1020) MAXQXY
!
       WRITE (output_io, 1030) RMC_MAX_Q
       WRITE (output_io, 1040) RMC_MAX_PLANES
       WRITE (output_io, 1050) RMC_MAX_SYM
       WRITE (output_io, 1060) RMC_MAX_ATOM
       WRITE (output_io, 1070) RMC_MAX_LOTS
!
       WRITE (output_io, 1080) CHEM_MAX_VEC
       WRITE (output_io, 1090) CHEM_MAX_COR
       WRITE (output_io, 1100) CHEM_MAX_BIN
       WRITE (output_io, 1110)  MMC_MAX_ATOM
       WRITE (output_io, 1120) MAX_ATOM_ENV
!
       WRITE (output_io, 1130) MK_MAX_ATOM
       WRITE (output_io, 1140) MK_MAX_SCAT
!
       WRITE (output_io, 1160) ST_MAXQXY
       WRITE (output_io, 1170) ST_MAXLAYER
       WRITE (output_io, 1180) ST_MAXTYPE
!
       WRITE (output_io, 1190) PDF_MAXDAT
!
       WRITE (output_io, 1200) MAXHIST
!
       WRITE (output_io, 6400) PL_MAXSCAT
!
       WRITE (output_io, 6600) SAV_MAXSCAT
!
       WRITE (output_io, 6700) SHEAR_MAXSCAT
!
       WRITE (output_io, 6800) SYM_MAXSCAT
!
       WRITE (output_io, 6900) TRAN_MAXSCAT
!
       WRITE (output_io, 7000) WV_MAXSCAT
!
!
 1000 FORMAT(' Current configuration of DISCUS : ')
 1005 FORMAT('   Maximum number of different atomtypes         : ',i8)
 1010 FORMAT('   Maximum number of atoms                       : ',i8)
 1020 FORMAT('   Maximum number of points in Q                 : ',i8)
 1030 FORMAT('   Maximum number of RMC input data points       : ',i8)
 1040 FORMAT('   Maximum number of RMC input data planes       : ',i8)
 1050 FORMAT('   Maximum number of RMC sym. equivalent planes  : ',i8)
 1060 FORMAT('   Maximum number of atoms moved in one RMC move : ',i8)
 1070 FORMAT('   Maximum number of RMC lots allowed            : ',i8)
 1080 FORMAT('   Maximum number of neighbour vec. definitions  : ',i8)
 1090 FORMAT('   Maximum number of correlations                : ',i8)
 1100 FORMAT('   Maximum number of points for histograms       : ',i8)
 1110 FORMAT('   Maximum number of different atoms             : ',i8)
 1120 FORMAT('   Maximum number of neighbouring atoms/mol.     : ',i8)
 1130 FORMAT('   Maximum number of atoms in microdomain        : ',i8)
 1140 FORMAT('   Maximum number of atom types in microdomain   : ',i8)
 1160 FORMAT('   Maximum number of points of Fourier           : ',i8)
 1170 FORMAT('   Maximum number of layers                      : ',i8)
 1180 FORMAT('   Maximum number of layer types                 : ',i8)
 1190 FORMAT('   Maximum number of points in PDF               : ',i8)
 1200 FORMAT('   Maximum number of points in Debye Histogram   : ',i8)
 6400 FORMAT('   Current number of atomtypes for plot          : ',i8)
 6600 FORMAT('   Current number of atomtypes for save          : ',i8)
 6700 FORMAT('   Current number of atomtypes for shear         : ',i8)
 6800 FORMAT('   Current number of atomtypes for symmetry      : ',i8)
 6900 FORMAT('   Current number of atomtypes for transform     : ',i8)
 7000 FORMAT('   Current number of atomtypes for waves         : ',i8)
!
    END SUBROUTINE discus_show_config
!
    SUBROUTINE discus_alloc_default
!
      USE discus_config_mod
      USE precision_mod
      IMPLICIT NONE
      INTEGER(KIND=PREC_INT_LARGE), PARAMETER :: ONE=1
      integer, dimension(3) :: n_qxy
!
      n_qxy = 1
!
      CALL alloc_chem_correlation( 1     )
      CALL alloc_chem_ang ( 1,  1        )
      CALL alloc_chem_aver( 1,  1        )
      CALL alloc_chem_disp( 1,  1        )
      CALL alloc_chem_env ( 1,  1,  1    )
      CALL alloc_chem_vec ( 1,  1        )
      CALL alloc_chem_con ( 1,  1        )
      CALL alloc_chem_dir ( 1 )
      CALL alloc_chem_ran ( 1,  1,  1    )
      CALL alloc_chem_dist( 1 )
      CALL alloc_chem_hist( 1            )
      CALL alloc_conn_vect( 1            )
      CALL alloc_crystal_scat( 1)
      CALL alloc_crystal_nmax( 1)
      call alloc_unitcell ( 1 )
      call alloc_anis     ( 1 )
      CALL alloc_deco     ( 1,  1,  4,  3,   3, 2 , 3)
      CALL alloc_debye    ( 1,  1,  1, ONE )
      CALL alloc_demol    ( 1,  1        ) 
      CALL alloc_diffuse_four  ( n_qxy)
      CALL alloc_diffuse_scat  ( 1    )
      call alloc_diffuse_dbw   ( n_qxy, 1)
      CALL alloc_diffuse_atom  ( 1    )
      CALL alloc_domain   ( 1            )
      CALL alloc_micro    ( 1,  1        )
      CALL alloc_mmc      ( 1,  8,  1,  1)
      CALL alloc_mmc_angle( 1,  1        )
      CALL alloc_mmc_buck ( 1,  1        )
      CALL alloc_mmc_lenn ( 1,  1        )
      CALL alloc_mmc_move ( 1,  1,  1    )
      CALL alloc_molecule ( 1,  1,  1,  1,  1)
      CALL alloc_phases   ( 1,  n_qxy,  1        )
      CALL alloc_pdf      ( 1,  1,  1,  1    )
      CALL alloc_plot     ( 1,  1,  1    )
      CALL alloc_powder   ( 1,  1        )
      CALL alloc_powder_partial( 1       )
      CALL alloc_powder_nmax ( 1,1          )
      CALL alloc_rmc      ( 1,  1        )
      CALL alloc_rmc_data ( n_qxy        )
      CALL alloc_rmc_istl ( (/1,1,1/),  1, 1     )
      CALL alloc_rmc_q    ( (/1,1,1/),  1, 1        )
      CALL alloc_rmc_planes(1, 48        )
      CALL alloc_save     ( 1, 1         )
      CALL alloc_shear    ( 1,  1        )
      CALL alloc_stack    ( 1,  1,  .TRUE.)
      CALL alloc_stack_four(n_qxy )
      CALL alloc_stack_crystal    ( 1,  1        )
      CALL alloc_surf     ( MAXSCAT      )
      CALL alloc_super    ( 1,  1, 1, 1  )
      CALL alloc_symmetry ( 1,  1        )
      CALL alloc_transfrm ( 1,  1        )
      CALL alloc_waves    ( 1,  1        )
!
    END SUBROUTINE discus_alloc_default
!
SUBROUTINE alloc_chem_correlation( n_cor)
!-
!     Allocate the arrays needed by CHEM correlations, for all types
!+
USE chem_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: n_cor
!
INTEGER :: all_status
!
CALL alloc_arr(chem_ctyp , 1,n_cor,  all_status, 0)
CALL alloc_arr(chem_nnei , 1,n_cor,  all_status, 0)
CALL alloc_arr(chem_neig , 1,3,  1, 48,  1,n_cor,  all_status, 0.0D0)
CALL alloc_arr(chem_ldall, 1,n_cor,  all_status, .TRUE.)
!
CHEM_MAX_COR = n_cor
!
END SUBROUTINE alloc_chem_correlation
!
    SUBROUTINE alloc_chem_aver ( n_atom_cell, n_max_atom )
!-
!     Allocate the arrays needed by CHEM average structure
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_atom_cell
      INTEGER, INTENT(IN)  :: n_max_atom
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
      CALL alloc_arr ( chem_ave_n    ,1,n_atom_cell ,   &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_ave_iscat,1,n_atom_cell , 1,n_max_atom  , &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_ave_pos  ,1,3 , 1,n_atom_cell ,   &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_ave_sig  ,1,3           ,  1,n_atom_cell , &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_ave_bese ,1,n_atom_cell , 1,n_max_atom  , &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         CHEM_MAXAT_CELL   = n_atom_cell
         CHEM_MAX_AVE_ATOM = n_max_atom
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         CHEM_MAXAT_CELL   = n_atom_cell
         CHEM_MAX_AVE_ATOM = n_max_atom
         NMAX          =  1
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_aver
!
!
    SUBROUTINE alloc_chem_ang  ( n_ang, n_cor )
!-
!     Allocate the arrays needed by CHEM angle definitions
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_ang
      INTEGER, INTENT(IN)  :: n_cor
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
      CALL alloc_arr ( chem_cwin     ,1,9, 1,n_ang,      &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_use_win  ,1,n_ang, 1,n_cor,  &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_nwin             , 1,n_cor,  &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         CHEM_MAX_ANG  = n_ang
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         CHEM_MAX_ANG  = n_ang
         NMAX          =  1
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_ang
!
!
    SUBROUTINE alloc_chem_disp ( n_cor, n_scat )
!-
!     Allocate the arrays needed by CHEM displacement structure
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_cor
      INTEGER, INTENT(IN)  :: n_scat
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
      CALL alloc_arr ( chem_disp_ave    ,1,n_cor ,   &
                       0, n_scat, 0, n_scat,         &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_disp_sig    ,1,n_cor ,   &
                       0, n_scat, 0, n_scat,         &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_disp
!
!
    SUBROUTINE alloc_chem_env ( n_atom, n_env, n_cor )
!-
!     Allocate the arrays needed by CHEM displacement structure
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_atom
      INTEGER, INTENT(IN)  :: n_env
      INTEGER, INTENT(IN)  :: n_cor
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
      CALL alloc_arr ( chem_use_env     ,1,n_env  , 1, n_cor,  &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_nenv                  , 1, n_cor,  &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cenv        ,0,n_atom , 1, n_env,  &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_env_neig    ,1,n_env ,   &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_rmax_env    ,1,n_env ,   &
                       all_status, 3.5D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_rmin_env    ,1,n_env ,   &
                       all_status, 0.1D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         CHEM_MAX_ENV  = n_env
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_env
!
!
    SUBROUTINE alloc_chem_ran  ( n_ran, n_cor, n_max_atom )
!-
!     Allocate the arrays needed by CHEM range definitions
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_ran
      INTEGER, INTENT(IN)  :: n_cor
      INTEGER, INTENT(IN)  :: n_max_atom
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
      CALL alloc_arr ( chem_use_ran  ,1,n_ran, 1,n_cor,      &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_nran     ,1,n_cor,      &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_cent,0,n_max_atom, 1,n_ran, &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
 
      CALL alloc_arr ( chem_cran_neig,0,n_max_atom, 1,n_ran, &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
 
      CALL alloc_arr ( chem_cran_nuvw ,1,n_ran    ,          &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
 
      CALL alloc_arr ( chem_cran_nshort ,1,n_ran  ,          &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_uvw ,1,3         , 1,48,    &
                                      1,n_ran     ,          &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      chem_cran_uvw(1,1,CHEM_MAX_RAN+1:n_ran) = -9999
!
      CALL alloc_arr ( chem_cran_sig  ,1,n_ran       ,      &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_wsig ,1,n_ran       ,      &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_rmax ,1,n_ran       ,      &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_rmin ,1,n_ran       ,      &
                       all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_cang ,1,n_ran       ,      &
                       all_status, .false.)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_lsym ,1,n_ran       ,      &
                       all_status, .false.)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cran_short,1,n_ran       ,      &
                       all_status, .false.)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
!
      IF( lstat ) THEN                        ! Success
         CHEM_MAX_RAN  = n_ran
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         CHEM_MAX_RAN  = n_ran
         NMAX          =  1
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_ran
!
!
    SUBROUTINE alloc_chem_vec  ( n_vec, n_cor )
!-
!     Allocate the arrays needed by CHEM vector definitions
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_vec
      INTEGER, INTENT(IN)  :: n_cor
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
!
      CALL alloc_arr ( chem_cvec     ,1,5,     1,n_vec,      &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_use_vec  ,1,n_vec, 1,n_cor,      &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_nvec     ,1,n_cor,      &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         CHEM_MAX_VEC  = n_vec
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         CHEM_MAX_VEC  = n_vec
         NMAX          =  1
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_vec
!
!
    SUBROUTINE alloc_chem_con  ( n_con, n_cor )
!-
!     Allocate the arrays needed by CHEM contor definitions
!+
      USE crystal_mod
      USE chem_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_con
      INTEGER, INTENT(IN)  :: n_cor
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
!
      CALL alloc_arr ( chem_ccon     ,1,2,     1,n_con,      &
                       all_status, -9999)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cname    ,         1,n_con,      &
                       all_status, ' '  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_cname_l  ,         1,n_con,      &
                       all_status,     1)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_use_con  ,1,n_con, 1,n_cor,      &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( chem_ncon  ,1,n_cor,      &
                       all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         CHEM_MAX_CON  = n_con
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'CHEMISTRY'
         ENDIF
      ELSE                                    ! Failure
         CHEM_MAX_CON  = n_con
         NMAX          =  1
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'CHEMISTRY'
         RETURN
      END IF
    END SUBROUTINE alloc_chem_con
!
SUBROUTINE alloc_chem_dir( n_cor)
!-
!     Allocate the arrays needed by CHEM correlations, for all types
!+
USE chem_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: n_cor
!
INTEGER :: all_status
!
CALL alloc_arr(chem_dir , 1, 3, 1,2,  1,n_cor,  all_status, 0.0D0)
!
END SUBROUTINE alloc_chem_dir
!
SUBROUTINE alloc_chem_dist( n_cor)
!-
!     Allocate the arrays needed by CHEM correlations, for all types
!+
USE chem_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: n_cor
!
INTEGER :: all_status
!
CALL alloc_arr(chem_rmax       , 1,n_cor,  all_status, 0.0D0)
CALL alloc_arr(chem_rmin       , 1,n_cor,  all_status, 0.0D0)
CALL alloc_arr(chem_freq_sigma , 1,n_cor,  all_status, 0.0D0)
CALL alloc_arr(chem_wink_sigma , 1,n_cor,  all_status, 0.0D0)
CALL alloc_arr(chem_cang       , 1,n_cor,  all_status, .FALSE.)
!
END SUBROUTINE alloc_chem_dist
!
SUBROUTINE alloc_chem_hist( n_hist)
!-
!     Allocate the arrays needed by CHEM correlations, for all types
!+
USE chem_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: n_hist
!
INTEGER :: all_status
!
CALL alloc_arr(chem_hist       , 1,n_hist,  all_status, 0)
!
CHEM_MAX_BIN = n_hist
!
END SUBROUTINE alloc_chem_hist
!
!*******************************************************************************
!
SUBROUTINE alloc_conn_vect(n_vec)
!-
!     Allocate the arrays needed by CONN vector definitions
!+
USE conn_def_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_vec
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat  = .TRUE.
!
!
CALL alloc_arr(conn_vectors ,1,5, 1,n_vec, all_status, -9999)
!
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
IF( lstat ) THEN                        ! Success
   CONN_MAX_VECT = n_vec
   ier_typ       = 0
   ier_num       = 0
   IF(all_status == 1) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'CONNECTIVITY'
   ENDIF
ELSE                                    ! Failure
   CONN_MAX_VECT = n_vec
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'CONNECTIVITY'
   RETURN
END IF
!
end SUBROUTINE alloc_conn_vect
!
!*******************************************************************************
!
subroutine alloc_crystal_nmax ( n_max )
!-
!     Allocate the arrays needed by CRYSTAL
!  Arrays related to number of atoms
!+
USE crystal_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_max
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat  = .TRUE.
!
CALL alloc_arr ( cr_iscat      , 1, 3 , 1,n_max , all_status, 0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_ianis      ,1,n_max ,  all_status, 0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_prop       ,1,n_max ,  all_status, 0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_mole       ,1,n_max ,  all_status, 0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_pos,1,3,1,n_max ,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_surf,0,3   ,1,n_max ,  all_status, 0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_magn,0,3   ,1,n_max ,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
IF( lstat ) THEN                        ! Success
   NMAX          = n_max
   cr_dim_natoms = n_max                ! Array dimension
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Crystal'
   ENDIF
ELSE                                    ! Failure
   NMAX          =  1
   cr_dim_natoms =  1                   ! Array dimension
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Crystal'
   RETURN
END IF
!
end subroutine alloc_crystal_nmax
!
!*******************************************************************************
!
subroutine alloc_crystal_scat ( n_scat)
!-
!     Allocate the arrays needed by CRYSTAL
!  Arrays related to number of scattering types
!+
USE crystal_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_scat
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat  = .TRUE.
!
!
CALL alloc_arr ( cr_scat,1,11,0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_delfr,       0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_delfi,       0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_delfr_u,     0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_delfi_u,     0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_scat_int,    0,n_scat,  all_status, .true.)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_scat_equ,    0,n_scat,  all_status, .false.)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_delf_int,    0,n_scat,  all_status, .true.)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_at_equ,      0,n_scat,  all_status, ' ')
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_ndiffer,     0,n_scat,  all_status, 1  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_at_lis,      0,n_scat,  all_status, ' ')
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
cr_at_lis(0) = 'VOID'
!
CALL alloc_arr ( as_at_lis,      0,n_scat,  all_status, ' ')
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
as_at_lis(0) = 'VOID'
!
CALL alloc_arr ( as_mole ,       0,n_scat,  all_status,  0 )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( as_prop ,       0,n_scat,  all_status,  0 )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_niscat,      0,n_scat,  all_status,  0 )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( as_iscat,       0,n_scat,  all_status,  0 )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_dw,          0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_anis, 1, 6,  0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_occ,         0,n_scat,  all_status, 1.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( as_occ,         0,n_scat,  all_status, 1.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( as_pos, 1,3,0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( as_dw,          0,n_scat,  all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!     CALL alloc_arr ( as_occ,         0,n_scat,  all_status, 1.0D0)
!     lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( cr_amount,      0,n_scat,  all_status, 0  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
IF( lstat ) THEN                        ! Success
   MAXSCAT       = n_scat
   cr_dim_nscat  = n_scat               ! Array dimension
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Crystal'
   ENDIF
ELSE                                    ! Failure
   MAXSCAT       =  1
   cr_dim_nscat  =  1                   ! Array dimension
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Crystal'
   RETURN
END IF
!
end subroutine alloc_crystal_scat
!
!*******************************************************************************
!
subroutine alloc_unitcell(n_atoms)
!-
!  Allocate arrays specific fo a single unit cell
!+
use crystal_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: n_atoms  ! Number atoms per unit cell
!
integer :: all_status
!
call alloc_arr(cr_is_sym, 1, n_atoms, all_status, 0)
!
if(all_status==0) then
  cr_dim_ncatoms = n_atoms
else
  cr_dim_ncatoms =  1
endif
!
end subroutine alloc_unitcell
!
!*******************************************************************************
!
subroutine alloc_anis(nanis)
!-
!  Allocate arrays for ADPs, actual crystal
!+
use crystal_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: nanis  ! Number different ADP types
!
integer :: all_status
!
call alloc_anis_generic(nanis, cr_anis_full, cr_prin, all_status)
!
if(all_status==0) then
  cr_dim_ncatoms = nanis
else
  cr_dim_ncatoms =  1
endif
!
end subroutine alloc_anis
!
!*******************************************************************************
!
subroutine alloc_anis_generic(nanis, anis_full, prin, all_status)
!-
!  Allocate arrays specific fo a single unit cell
!+
use precision_mod
!
implicit none
!
integer, intent(in) :: nanis  ! Number atoms per unit cell
real(kind=PREC_DP), dimension(:,:)  , allocatable, intent(inout) :: anis_full
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(inout) :: prin
!
integer, intent(out) :: all_status
!
call alloc_arr(anis_full, 1, 6, 1, nanis, all_status, 1.0D0)
call alloc_arr(prin, 1, 4, 1, 3,  1, nanis, all_status, 1.0D0)
!call alloc_arr(cr_u2  , 1, 3, 1, nanis, all_status, 1.0D0)
!
end subroutine alloc_anis_generic
!
!*******************************************************************************
!
    SUBROUTINE alloc_debye ( n_scat, n_hist, n_qxy, MASK )
!-
!     Allocate the arrays needed by DEBYE
!+
      USE debye_mod
      USE precision_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_hist
      INTEGER, INTENT(IN)  :: n_qxy
      INTEGER(KIND=PREC_INT_LARGE), INTENT(IN)  :: MASK
!
      INTEGER              :: n_look
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      n_look = n_scat*(n_scat+1)/2
      lstat  = .TRUE.
!
      CALL alloc_arr ( rsf     ,1,n_qxy  ,  all_status, 0.0D0 )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sinetab ,0,INT(MASK)   ,  all_status, 0.0D0 )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         MAXLOOK       = n_look
         MAXHIST       = n_hist
         MAXDSCAT      = n_scat
         MAXDQXY       = n_qxy
         DEB_MAXMASK   = INT(MASK)
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Debye'
         ENDIF
      ELSE                                    ! Failure
         MAXLOOK       =  1
         MAXHIST       =  1
         MAXDSCAT      =  1
         MAXDQXY       = n_qxy
         DEB_MAXMASK   =  0
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Debye'
         RETURN
      END IF
    END SUBROUTINE alloc_debye
!
!
    SUBROUTINE alloc_deco ( n_scat, n_site, n_deco, n_anch, n_hkl, n_new, m_scat)
!-
!     Allocate the arrays needed by DECORATE
!+
      USE deco_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
      INTEGER, INTENT(IN)  :: n_deco
      INTEGER, INTENT(IN)  :: n_anch
      INTEGER, INTENT(IN)  :: n_hkl
      INTEGER, INTENT(IN)  :: n_new
      INTEGER, INTENT(IN)  :: m_scat
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( dc_latom       ,0,n_scat  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( dc_lsite       ,0,n_site  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_name,1,n_deco ,  all_status, ' '   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_file,1,n_deco ,  all_status, ' '   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_natoms,1,n_deco, all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_atom_name,0,m_scat, 1,n_deco, all_status, ' '   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_adp      ,0,m_scat, 1,n_deco, all_status, 0.0D0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_biso,1,n_deco ,  all_status, 0.0D0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_clin,1,n_deco ,  all_status, 0.0D0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_cqua,1,n_deco ,  all_status, 0.0D0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_mole_type,1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_type,1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_lname,1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_lfile,1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_nanch ,1,2,         1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_surf  ,0,2,0,n_anch,1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_neig,0,2,           1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_secnd,              1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_axis,0,2,           1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_lrestrict,          1,n_deco ,  all_status, .FALSE. )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_spread,             1,n_deco ,  all_status, .TRUE.  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_lform,              1,n_deco ,  all_status, .FALSE. )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_hkl ,1,3, 0,n_hkl,  1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_surfnew,0,n_new,    1,n_deco ,  all_status, 0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_dens,               1,n_deco ,  all_status, 0.0D0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_dist,1,2,           1,n_deco ,  all_status, 0.0D0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_angle,              1,n_deco ,  all_status, 170.0D0 )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_tilt,               1,n_deco ,  all_status,   0.0D0 )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_tilt_hkl, 1, 3,     1,n_deco ,  all_status,   0.0D0 )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_tilt_atom, 1, 4,    1,n_deco ,  all_status,   0   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_tilt_is_atom,       1,n_deco ,  all_status,.TRUE. )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( dcc_tilt_is_auto,       1,n_deco ,  all_status,.TRUE. )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         DCC_MAXNUM    = n_deco
          DC_MAXSITE   = n_site
         DCC_MAXANCH   = n_anch
         DCC_MAXHKL    = n_hkl
         DCC_MAXNEW    = n_new
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Save'
         ENDIF
      ELSE                                    ! Failure
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Deco'
         RETURN
      END IF
    END SUBROUTINE alloc_deco
!
!
SUBROUTINE alloc_demol ( n_atomtype, n_moletype )
!-
!     Allocate the arrays needed by demolecularize
!+
USE crystal_mod
USE demolec_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_atomtype
INTEGER, INTENT(IN)  :: n_moletype
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat     = .TRUE.
!
CALL alloc_arr ( dem_latomtype   ,0,n_atomtype,  all_status, .FALSE. )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( dem_lmoletype   ,1,n_moletype,  all_status, .FALSE. )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
IF( lstat ) THEN                        ! Success
   DEM_MAX_MOLETYPE  = n_moletype
   DEM_MAX_ATOMTYPE  = n_atomtype
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'demolec'
   ENDIF
ELSE                                    ! Failure
   DEM_MAX_MOLETYPE  = n_moletype
   DEM_MAX_ATOMTYPE  = n_atomtype
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'demolec'
   RETURN
END IF
!
END SUBROUTINE alloc_demol
!
subroutine alloc_diffuse_four ( n_qxy )
!-
!  Allocate the arrays needed by DIFFUSE
!  Arrays related to extend of reciprocal space
!+
use diffuse_mod
use precision_mod
!
implicit none
!
!      
integer, dimension(3), intent(in)  :: n_qxy
!
complex (kind=kind(0.0D0)) :: def_value
integer              :: all_status
logical              :: lstat
!
lstat     = .TRUE.
def_value = CMPLX(0.0D0,0.0D0,kind=PREC_DP)
!
call alloc_arr ( csf     ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3),  all_status, def_value)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
call alloc_arr ( tcsf    ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3), all_status, def_value)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
call alloc_arr ( acsf    ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3), all_status, def_value)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
call alloc_arr ( istl    ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3),  all_status, 0        )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
call alloc_arr ( dsi     ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3),  all_status, 0.0D0    )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
if( lstat ) then                        ! Success
   MAXQXY        = n_qxy
   ier_typ       = 0
   ier_num       = 0
   if ( all_status == 1 ) then
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Diffuse'
   endif
else                                    ! Failure
   MAXQXY        =  1
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Diffuse'
   return
end if
!
end subroutine alloc_diffuse_four
!
!
subroutine alloc_diffuse_scat (  n_scat )
!-
!  Allocate the arrays needed by DIFFUSE
!  Arrays related to extend of reciprocal space
!+
use diffuse_mod
use precision_mod
!
implicit none
!
!      
integer, intent(in)  :: n_scat
!
complex (kind=kind(0.0D0)) :: def_value
integer              :: all_status
logical              :: lstat
!
lstat     = .TRUE.
def_value = CMPLX(0.0D0,0.0D0,kind=PREC_DP)
!
!
call alloc_arr ( cfact   ,0,CFPKT  ,  1, n_scat, &
                                 all_status, def_value)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
call alloc_arr ( cfact_pure,0,CFPKT  ,  1, n_scat, &
                                 all_status, def_value)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
if( lstat ) then                        ! Success
   DIF_MAXSCAT   = n_scat
   ier_typ       = 0
   ier_num       = 0
   if ( all_status == 1 ) then
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Diffuse'
   endif
else                                    ! Failure
   DIF_MAXSCAT   = n_scat
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Diffuse'
   return
end if
!
end subroutine alloc_diffuse_scat
!
!
subroutine alloc_diffuse_dbw(n_qxy, n_anis)
!-
!  Allocate Dbewy-Waller Map
!+
!
use diffuse_mod
!
use precision_mod
!      
integer, dimension(3), intent(in)  :: n_qxy
integer              , intent(in)  :: n_anis
!
real(kind=PREC_DP)   :: def_value
integer              :: all_status
logical              :: lstat
!
lstat     = .TRUE.
def_value = 0.0D0
!
call alloc_arr(four_dbw, 1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3), 1,n_anis, &
                                 all_status, def_value)
!
end subroutine alloc_diffuse_dbw
!
!
subroutine alloc_diffuse_atom ( n_atoms )
!-
!  Allocate the arrays needed by DIFFUSE
!  Arrays related to extend of reciprocal space
!+
use diffuse_mod
use precision_mod
!
implicit none
!
!      
integer, intent(in)  :: n_atoms
!
complex (kind=kind(0.0D0)) :: def_value
integer              :: all_status
logical              :: lstat
!
lstat     = .TRUE.
def_value = CMPLX(0.0D0,0.0D0,kind=PREC_DP)
!
!
call alloc_arr ( xat     ,1,n_atoms,  1, 3, all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
if( lstat ) then                        ! Success
   DIF_MAXAT     = n_atoms
   ier_typ       = 0
   ier_num       = 0
   if ( all_status == 1 ) then
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Diffuse'
   endif
else                                    ! Failure
   DIF_MAXAT     =  1
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Diffuse'
   return
end if
!
end subroutine alloc_diffuse_atom
!
!
    SUBROUTINE alloc_domain ( n_clu )
!-
!     Allocate the arrays needed by domain
!+
      USE crystal_mod
      USE domain_mod
      USE domaindis_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_clu
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( clu_content   ,1,n_clu,  all_status, ' ' )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_name      ,1,n_clu,  all_status, ' ' )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_use_expr ,1,n_clu, 1, 2,  all_status, .FALSE.  )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_character ,1,n_clu,  all_status, -4  )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_fuzzy     ,1,n_clu,  all_status,2.55D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_orient    ,1,n_clu, 1,3,1,4, all_status, 0.0D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_or_expr   ,1,n_clu, 1,3,1,4, all_status, ' ' )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_shape     ,1,n_clu, 1,3,1,4, all_status, 0.0D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( clu_sigma     ,1,n_clu,    1,3,  all_status, 0.0D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         CLU_MAX_TYPE  = n_clu
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'domain'
         ENDIF
      ELSE                                    ! Failure
         CLU_MAX_TYPE  = n_clu
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'domain'
         RETURN
      END IF
    END SUBROUTINE alloc_domain
!
!
    SUBROUTINE alloc_micro ( n_scat, n_atoms )
!-
!     Allocate the arrays needed by micro
!+
      USE crystal_mod
      USE micro_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_atoms
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( mk_at_lis    ,0,n_scat,  all_status, ' ' )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( mk_dw         ,0,n_scat,  all_status, 0.0D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( mk_anis, 1, 6 ,0,n_scat,  all_status, 0.0D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( mk_occ        ,0,n_scat,  all_status, 1.00D0 )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         MK_MAX_SCAT   = n_scat
         MK_MAX_ATOM   = n_atoms
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'micro'
         ENDIF
      ELSE                                    ! Failure
         MK_MAX_SCAT   = n_scat
         MK_MAX_ATOM   = n_atoms
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'micro'
         RETURN
      END IF
    END SUBROUTINE alloc_micro
!
!
    SUBROUTINE alloc_mmc ( n_corr, n_ener, n_scat, n_site )
!-
!     Allocate the arrays needed by MMC 
!+
      USE chem_mod
      USE mmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_corr
      INTEGER, INTENT(IN)  :: n_ener
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
integer :: i_scat
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
!
      lstat       = .TRUE.
!
!
       CALL alloc_arr ( mmc_latom     ,0,n_scat  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!      CALL alloc_arr ( mmc_lsite     ,0,n_site  ,  all_status, .false.   )
!      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_allowed   ,0,n_scat  ,  all_status, .false.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_nvec    ,1,n_corr ,  &
                                    0,n_scat ,  &
                                    0,n_scat ,  &
                                      all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_target_corr ,1,n_corr , &
                                        0,n_ener , &
                                        0,n_scat , &
                                        0,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_depth       ,1,n_corr , &
                                        0,n_ener , &
                                        0,n_scat , &
                                        0,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_depth_def   ,1,n_corr , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ach_corr    ,1,n_corr , &
                                        0,n_ener , &
                                       -1,n_scat , &
                                       -1,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ach_pairs   ,1,n_corr , &
                                        0,n_ener , &
                                       -1,n_scat , &
                                       -1,n_scat , &
                                      all_status, 0      )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ini_corr    ,1,n_corr , &
                                        0,n_ener , &
                                       -1,n_scat , &
                                       -1,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ach_sigm    ,1,n_corr ,  0,n_ener , &
                                       -1,n_scat , -1,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ini_sigm    ,1,n_corr ,  0,n_ener , &
                                       -1,n_scat , -1,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out

!     CALL alloc_arr ( mmc_vec         ,1,4, 1,12,  1,n_corr , &
!                                      -1,n_scat , -1,n_scat , &
!                                     all_status, 0.0  )
!     lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_const       ,0,n_corr ,  &
                                        0,n_ener ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_cfac        ,0,n_corr ,  &
                                        0,n_ener ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_lfeed       ,0,n_corr ,  &
                                        0,n_ener ,  &
                                      all_status, .true.  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_pre_corr    ,1,n_corr , &
                                        0,n_ener , &
                                       -1,n_scat , &
                                       -1,n_scat , &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
 i_scat = max(3, n_scat)
      CALL alloc_arr ( mmc_pneig       ,0,i_scat , &
                                        0,i_scat , &
                                        1,n_corr , &
                                      all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_cor_energy  ,0,n_corr ,  &
                                        0,n_ener ,  &
                                      all_status,.false.)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_pair        ,1,n_corr , &
                                        0,n_ener , &
                                       -1,n_scat , &
                                       -1,n_scat , &
                                      all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_left        ,1,n_corr , &
                                        0,n_ener , &
                                        0,n_scat , &
                                      all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_right       ,1,n_corr , &
                                        0,n_ener , &
                                        0,n_scat , &
                                      all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         MMC_MAX_CORR     = n_corr
         MMC_MAX_SCAT     = n_scat
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'MMC'
         ENDIF
      ELSE                                    ! Failure
         MMC_MAX_CORR     = n_corr
         MMC_MAX_SCAT     = n_scat
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'MMC'
         RETURN
      END IF
    END SUBROUTINE alloc_mmc
!
!*******************************************************************************
!
subroutine alloc_mmc_pid ( n_corr, n_ener, n_scat, n_site )
!-
!     Allocate the arrays needed by MMC 
!+
USE chem_mod
USE mmc_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_corr
INTEGER, INTENT(IN)  :: n_ener
INTEGER, INTENT(IN)  :: n_scat
INTEGER, INTENT(IN)  :: n_site
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat       = .TRUE.
!
!
CALL alloc_arr ( mmc_pid_diff ,1,n_corr , &
                               0,n_ener , &
                              -1,n_scat , &
                              -1,n_scat , &
                              all_status, 0.0D0  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( mmc_pid_inte ,1,n_corr , &
                               0,n_ener , &
                              -1,n_scat , &
                              -1,n_scat , 0,2 ,&
                              all_status, 0.0D0  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( mmc_pid_deri ,1,n_corr , &
                               0,n_ener , &
                              -1,n_scat , &
                              -1,n_scat , 0,2, &
                              all_status, 0.0D0  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
end subroutine alloc_mmc_pid
!
!*******************************************************************************
!
SUBROUTINE alloc_mo_ach  ( n_corr)
!-
!     Allocate the arrays needed by MMC  for moves in mc_mod
!+
USE mc_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_corr    ! Number of targets
!
INTEGER              :: all_status
!
CALL alloc_arr(mo_ach_corr   ,1,n_corr                    ,  all_status, 0.0D0     )
!
END SUBROUTINE alloc_mo_ach
!
!*******************************************************************************
!
SUBROUTINE alloc_mmc_move( n_corr, n_scat, n_mole)
!-
!     Allocate the arrays needed by MMC  for moves in mc_mod
!+
USE mc_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_corr    ! Number of targets
INTEGER, INTENT(IN)  :: n_scat    ! Number of atom types
INTEGER, INTENT(IN)  :: n_mole    ! Number of molecule types
!
INTEGER              :: all_status
!
!CALL alloc_arr(mo_target_corr,1,n_corr                    ,  all_status, 0.0       )
!CALL alloc_arr(mo_ach_corr   ,1,n_corr                    ,  all_status, 0.0       )
!CALL alloc_arr(mo_const      ,0,n_corr                    ,  all_status, 0.0       )
!CALL alloc_arr(mo_cfac       ,0,n_corr                    ,  all_status, 0.0       )
!CALL alloc_arr(mo_disp       ,1,n_corr, 0,n_scat, 0,n_scat,  all_status, 0.0       )
CALL alloc_arr(mo_maxmove     ,1,4     , 0,n_scat          ,  all_status, 0.0D0     )
CALL alloc_arr(mo_maxmove_mole,1,4     , 0,n_mole          ,  all_status, 0.0D0     )
CALL alloc_arr(mo_maxrota_mole,1,4     , 0,n_mole          ,  all_status, 0.0D0     )
!
END SUBROUTINE alloc_mmc_move
!
    SUBROUTINE alloc_mmc_angle ( n_corr, n_angle )
!-
!     Allocate the arrays needed by MMC ANGLE potential
!+
      USE chem_mod
      USE mmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_corr
      INTEGER, INTENT(IN)  :: n_angle
!
      INTEGER              :: n_size
      INTEGER              :: n_old
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
!
      n_size      = n_corr       * n_angle
      n_old       = CHEM_MAX_COR * MMC_MAX_ANGLES
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_angles    ,1,n_size ,  &
                                      all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_target_angl,1,n_size ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_depth_angl,1,n_size ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ach_angl  ,1,n_size ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ang_sigm  ,1,n_size ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ini_angl  ,1,n_size ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_ini_sang  ,1,n_size ,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         MMC_MAX_ANGLES   = n_angle
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'MMC_ANGLES'
         ENDIF
      ELSE                                    ! Failure
         MMC_MAX_ANGLES   = n_angle
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'MMC_ANGLES'
         RETURN
      END IF
    END SUBROUTINE alloc_mmc_angle
!
!
    SUBROUTINE alloc_mmc_buck ( n_corr, n_scat )
!-
!     Allocate the arrays needed by MMC Buckingham potential
!+
      USE mmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_corr
      INTEGER, INTENT(IN)  :: n_scat
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_buck_a    ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_buck_rho  ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_buck_b    ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_buck_rmin ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_buck_atmin,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         MMC_BUCK_CORR    = n_corr
         MMC_BUCK_SCAT    = n_scat
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'MMC_BUCKINGHAM'
         ENDIF
      ELSE                                    ! Failure
         MMC_BUCK_CORR    = n_corr
         MMC_BUCK_SCAT    = n_scat
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'MMC_BUCKINGHAM'
         RETURN
      END IF
    END SUBROUTINE alloc_mmc_buck
!
!
    SUBROUTINE alloc_mmc_lenn ( n_corr, n_scat )
!-
!     Allocate the arrays needed by MMC Lennard potential
!+
      USE mmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_corr
      INTEGER, INTENT(IN)  :: n_scat
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_len_a     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_len_b     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_len_m     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_len_n     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         MMC_LENN_CORR    = n_corr
         MMC_LENN_SCAT    = n_scat
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'MMC_LENNARD'
         ENDIF
      ELSE                                    ! Failure
         MMC_LENN_CORR    = n_corr
         MMC_LENN_SCAT    = n_scat
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'MMC_Lennard'
         RETURN
      END IF
    END SUBROUTINE alloc_mmc_lenn
!
!
    SUBROUTINE alloc_mmc_rep ( n_corr, n_scat )
!-
!     Allocate the arrays needed by MMC Repulsive potential
!+
      USE mmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_corr
      INTEGER, INTENT(IN)  :: n_scat
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mmc_rep_a     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_rep_b     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_rep_c     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mmc_rep_m     ,0,n_corr,  &
                                      0,n_scat,  &
                                      0,n_scat,  &
                                      all_status, 0.0D0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         MMC_REP_CORR    = n_corr
         MMC_REP_SCAT    = n_scat
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'MMC_Repulsive'
         ENDIF
      ELSE                                    ! Failure
         MMC_REP_CORR    = n_corr
         MMC_REP_SCAT    = n_scat
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'MMC_Repulsive'
         RETURN
      END IF
    END SUBROUTINE alloc_mmc_rep
!
!
    SUBROUTINE alloc_molecule ( n_gene, n_symm, n_mole, n_type, n_atom )
!-
!     Allocate the arrays needed by molecules
!+
      USE molecule_mod
      USE external_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_gene
      INTEGER, INTENT(IN)  :: n_symm
      INTEGER, INTENT(IN)  :: n_mole
      INTEGER, INTENT(IN)  :: n_type
      INTEGER, INTENT(IN)  :: n_atom
!
      INTEGER              :: i
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat       = .TRUE.
!
      CALL alloc_arr ( mole_gene_power,1,n_gene,  all_status, 1  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_gene,1,4, 1,4,1,n_gene,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_symm_power,1,n_symm,  all_status, 1  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_symm,1,4, 1,4,1,n_symm,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_len       ,0,n_mole,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_off       ,0,n_mole,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_type      ,0,n_mole,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_char      ,0,n_mole,  all_status, MOLE_ATOM  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_file      ,0,n_mole,  all_status, ' '  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_biso      ,0,n_type,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_clin      ,0,n_type,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_cqua      ,0,n_type,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_dens      ,0,n_mole,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_fuzzy     ,0,n_mole,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( mole_cont      ,0,n_atom,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( exte_names     ,1,n_type,  all_status, ' '  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( exte_length    ,1,n_type,  all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( exte_type      ,1,n_type,  all_status, 0    )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      FORALL (i=1:4)
         mole_gene(i,i,:) = 1.0
         mole_symm(i,i,:) = 1.0
      END FORALL
!
      IF( lstat ) THEN                        ! Success
         MOLE_MAX_GENE    = n_gene
         MOLE_MAX_SYMM    = n_symm
         MOLE_MAX_MOLE    = n_mole
         MOLE_MAX_TYPE    = n_type
         MOLE_MAX_ATOM    = n_atom
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'MOLECULE'
         ENDIF
      ELSE                                    ! Failure
         MOLE_MAX_GENE    = n_gene
         MOLE_MAX_SYMM    = n_symm
         MOLE_MAX_MOLE    = n_mole
         MOLE_MAX_TYPE    = n_type
         MOLE_MAX_ATOM    = n_atom
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'MOLECULE'

         RETURN
      END IF
    END SUBROUTINE alloc_molecule
!
!*******************************************************************************
!
SUBROUTINE alloc_phases ( n_pha, n_qxy, n_scat )
!-
!     Allocate the arrays needed by PHASES
!+
USE phases_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_pha    ! Number of phases
INTEGER, dimension(3), INTENT(IN)  :: n_qxy    ! Number of point in powder pattern
INTEGER, INTENT(IN)  :: n_scat   ! Number of atom types
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat     = .TRUE.
!
CALL alloc_arr(pha_nscat           , 1,n_pha, all_status, 0         )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_calc            , 1,n_pha, all_status, 0         )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_frac            , 1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_occ   , 0,n_scat, 1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_weight,           1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_nreal ,           1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_ncreal,           1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_powder, 0,product(n_qxy) , 1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_form  , 0,product(n_qxy) , 0, n_scat, 1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_adp   , 0,n_scat, 1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_scale ,           1,n_pha, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr(pha_niscat, 0,n_scat, 1,n_pha, all_status, 0         )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
IF( lstat ) THEN                        ! Success
   PHA_MAXPHA    = n_pha
   PHA_MAXPTS    = product(n_qxy)
   PHA_MAXSCAT   = n_scat
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Phases'
   ENDIF
ELSE                                    ! Failure
   PHA_MAXPHA    = n_pha
   PHA_MAXPTS    = product(n_qxy)
   PHA_MAXSCAT   = n_scat
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Phases'
END IF
!
END SUBROUTINE alloc_phases
!
!*******************************************************************************
!
    SUBROUTINE alloc_pdf ( n_scat, n_site, n_dat, n_bnd )
!-
!     Allocate the arrays needed by PDF
!+
      USE pdf_mod
      USE precision_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
      INTEGER, INTENT(IN)  :: n_dat
      INTEGER, INTENT(IN)  :: n_bnd
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      REAL(PREC_DP) , PARAMETER   :: def_dbl = 0.0D0
      INTEGER              :: n_dat2    !  Size rounded up to power of 2
!
      lstat       = .TRUE.
!
      n_dat2 = 2**(INT(LOG(REAL(n_dat))/LOG(2.))+2)
      CALL alloc_arr ( pdf_calc      ,1,n_dat2,  all_status, def_dbl  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_ppp       ,1,n_dat2,  all_status, def_dbl  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_corr      ,1,n_dat,  all_status, def_dbl  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_temp      ,0,n_dat,0,n_scat, &
                                              0,n_scat, 0,0, &
                                              all_status, 0      )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      CALL alloc_arr ( pdf_lsite_i   ,0,n_site,  all_status, .true.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
      CALL alloc_arr ( pdf_weight    ,0,n_scat,                      &
                                      0,n_scat, all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_obs       ,1,n_dat,  all_status, 0.0D0      )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_wic       ,1,n_dat,  all_status, 0.0D0      )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!     CALL alloc_arr ( pdf_sinc      ,1,2*n_dat+1, all_status, def_dbl  )
!     lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_sincc     ,1,PDF_MAXSINCC, all_status, def_dbl  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_exp       ,0,20000 ,  all_status, def_dbl  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_allowed_i ,0,n_scat,  all_status, .true.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_allowed_j ,0,n_scat,  all_status, .true.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_lsite_i   ,0,n_site,  all_status, .true.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_lsite_j   ,0,n_site,  all_status, .true.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_has_atom  ,0,n_scat,  all_status, 0        )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( pdf_bnd       ,1,3, -n_bnd,2*n_bnd, &
                                                              all_status, 0        )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         PDF_MAXSCAT      = n_scat
         PDF_MAXSITE      = n_site
         PDF_MAXDAT       = n_dat
         PDF_MAXBND       = n_bnd
         pdf_allowed_i(0) = .false.
         pdf_allowed_j(0) = .false.
         ier_typ          = 0
         ier_num          = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'PDF'
         ENDIF
      ELSE                                    ! Failure
         PDF_MAXSCAT   = n_scat
         PDF_MAXSITE   = n_site
         PDF_MAXDAT    = n_dat
         PDF_MAXBND    = n_bnd
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'PDF'

         RETURN
      END IF
    END SUBROUTINE alloc_pdf
!
!
    SUBROUTINE alloc_plot ( n_scat, n_site, n_line )
!-
!     Allocate the arrays needed by PLOT
!+
      USE discus_plot_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
      INTEGER, INTENT(IN)  :: n_line
!
      INTEGER              :: i
      INTEGER, DIMENSION(0:4), PARAMETER :: stift = (/ 5, 3, 1, 2, 4 /)

      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( pl_siz       ,0,n_scat,  all_status, 1.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_rgb       ,1,3,  0,n_scat,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_bond_len  ,1,2,  0, n_scat, &
                                            0, n_scat, all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_bond_rad  ,           0, n_scat,  &
                                                 0, n_scat,  all_status, 0.2D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_bond_col  ,1,3,  0, n_scat, &
                                            0, n_scat, all_status, 0.8D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_typ       ,0,n_scat,  all_status, 3         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_color     ,0,n_scat,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_latom     ,0,n_scat,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_lsite     ,0,n_site,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_batom_a   ,0,n_scat,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_batom_e   ,0,n_scat,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_bond      ,           0, n_scat,  &
                                                 0, n_scat,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_poly_c    ,          -1, n_scat,  all_status, .FALSE.  )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( pl_poly_o    ,          -1, n_scat,  all_status, .FALSE.  )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
call alloc_arr ( pl_lines, 1,3, 1, 2, 1, n_line, all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         PL_MAXSCAT    = n_scat
         PL_MAXSITE    = n_site
         PL_MAXLINE    = n_line
         pl_rgb(1, :)  = 1
         FORALL (i=0:PL_MAXSCAT) pl_color(i) =  stift (mod (i, 5) )
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Plot'
         ENDIF
      ELSE                                    ! Failure
         PL_MAXSCAT    = n_scat
         PL_MAXSITE    = n_site
         PL_MAXLINE    = n_line
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Plot'
         RETURN
      END IF
    END SUBROUTINE alloc_plot
!
!*******************************************************************************
!
subroutine alloc_plot_line(n_line)
!-
!  Alloc the number of JMOL lines
!+
USE discus_plot_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_line
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat     = .TRUE.
!
call alloc_arr ( pl_lines, 1,3, 1, 2, 1, n_line, all_status, 0.0D0)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
IF( lstat ) THEN                        ! Success
   PL_MAXLINE    = n_line
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Plot'
   ENDIF
ELSE                                    ! Failure
   PL_MAXLINE    = n_line
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Plot'
   RETURN
END IF
!
end subroutine alloc_plot_line
!
!*******************************************************************************
!
SUBROUTINE alloc_powder ( n_qxy, n_scat )
!-
!     Allocate the arrays needed by POWDER
!+
USE powder_mod
USE powder_scat_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_qxy
INTEGER, INTENT(IN)  :: n_scat
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat     = .TRUE.
!
CALL alloc_arr ( pow_qsp  ,0,n_qxy  ,  all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_f2aver  ,0,n_qxy  ,  all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_faver2  ,0,n_qxy  ,  all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_f2      ,0,n_qxy  ,  0, n_scat, all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_fu      ,0,n_qxy  ,  all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_conv    ,0,n_qxy  ,  all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_sq      ,0,n_qxy  ,  all_status, 0.0D0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
IF( lstat ) THEN                        ! Success
   POW_MAXPKT    = n_qxy
   POW_MAXSCAT   = n_scat
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Powder'
   ENDIF
ELSE                                    ! Failure
   POW_MAXPKT    = n_qxy
   POW_MAXSCAT   = n_scat
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Powder'
   RETURN
END IF
!
END SUBROUTINE alloc_powder
!
!*******************************************************************************
!
SUBROUTINE alloc_powder_partial ( n_scat )
!-
!     Allocate the arrays needed by POWDER
!+
USE powder_mod
!USE powder_scat_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_scat
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat     = .TRUE.
!
!
call alloc_arr(pow_do_partial, 0, n_scat, 0, n_scat, all_status, .TRUE.)
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
call alloc_arr(pow_nn_partial, 0, n_scat,            all_status, 0     )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
IF( lstat ) THEN                        ! Success
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Powder partial'
   ENDIF
ELSE                                    ! Failure
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Powder partial'
   RETURN
END IF
!
END SUBROUTINE alloc_powder_partial
!
!*******************************************************************************
!
SUBROUTINE alloc_powder_nmax ( n_scat, n_max )
!-
!     Allocate the arrays needed by POWDER complete mode
!+
USE powder_scat_mod
!
IMPLICIT NONE
!
!      
INTEGER, INTENT(IN)  :: n_scat
INTEGER, INTENT(IN)  :: n_max
!
INTEGER              :: all_status
LOGICAL              :: lstat
!
lstat     = .TRUE.
!
CALL alloc_arr ( pow_nscat, 1,n_scat, all_status, 0  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
CALL alloc_arr ( pow_iatom, 1,n_scat, &
                 0,n_max ,  all_status, 0  )
lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
IF( lstat ) THEN                        ! Success
   POW_MAXSCAT   = n_scat
   POW_NMAX      = n_max
   ier_typ       = 0
   ier_num       = 0
   IF ( all_status == 1 ) THEN
      ier_typ       = 1
      ier_num       = ER_COMM
      ier_msg(1)    = 'Powder'
   ENDIF
ELSE                                    ! Failure
   POW_MAXSCAT   = n_scat
   POW_NMAX      = n_max
   ier_num       = -3
   ier_typ       = ER_COMM
   ier_msg(1)    = 'Powder'
   RETURN
END IF
END SUBROUTINE alloc_powder_nmax
!
!*******************************************************************************
!
    SUBROUTINE alloc_rmc ( n_scat, n_site )
!-
!     Allocate the arrays needed by RMC
!+
      USE rmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( rmc_allowed    ,0,n_scat  ,  all_status, .true.    )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_lsite      ,0,n_site  ,  all_status, .true.    )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_maxmove,1,3,0,n_scat  ,  all_status, 0.2D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_mindist,1,n_scat,1,n_scat,  all_status, 0.5D0  )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         RMC_MAXSCAT   = n_scat
         RMC_MAXSITE   = n_site
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'RMC'
         ENDIF
      ELSE                                    ! Failure
         RMC_MAXSCAT   = n_scat
         RMC_MAXSITE   = n_site
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'RMC'
         RETURN
      END IF
    END SUBROUTINE alloc_rmc
!
!
    SUBROUTINE alloc_rmc_data ( n_qxy )
!-
!     Allocate the data intensity arrays needed by RMC
!+
      USE rmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, dimension(3), INTENT(IN)  :: n_qxy
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( rmc_int        ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3), all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_wic        ,1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3), all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         RMC_MAX_Q     = n_qxy
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'RMC'
         ENDIF
      ELSE                                    ! Failure
         RMC_MAX_Q     = n_qxy
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'RMC'
         RETURN
      END IF
    END SUBROUTINE alloc_rmc_data
!
!
    SUBROUTINE alloc_rmc_istl ( n_sq, n_scat, n_planes )
!-
!     Allocate the arrays needed by RMC
!+
      USE rmc_mod
      USE precision_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, dimension(3), INTENT(IN)  :: n_sq
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_planes
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      COMPLEX (KIND=PREC_DP) :: def_value
!
      def_value = CMPLX(0.0D0,0.0D0,KIND=PREC_DP)
      lstat     = .TRUE.
!
       CALL alloc_arr ( ristl, 1,n_sq(1), 1, n_sq(2), 1, n_sq(3)  , all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rcfact, 0, CFPKT,1,n_scat,1, n_planes , all_status, def_value   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
!        RMC_MAXSCAT   = n_scat
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'RMC'
         ENDIF
      ELSE                                    ! Failure
!        RMC_MAXSCAT   = n_scat
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'RMC'
         RETURN
      END IF
    END SUBROUTINE alloc_rmc_istl
!
!
    SUBROUTINE alloc_rmc_q ( n_sq, n_sym, n_lots )
!-
!     Allocate the arrays needed by RMC
!+
      USE rmc_mod
      USE precision_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, dimension(3), INTENT(IN)  :: n_sq
      INTEGER, INTENT(IN)  :: n_sym 
      INTEGER, INTENT(IN)  :: n_lots
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      COMPLEX (KIND=PREC_DP) :: def_value
!
      def_value = CMPLX(0.0D0,0.0D0,KIND=PREC_DP)
      lstat     = .TRUE.
!
       CALL alloc_arr ( rmc_lots_orig,1,3, 1,n_lots, all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_csf    , 1, n_sq(1), 1, n_sq(2), 1, n_sq(3), 1, n_sym, 1, n_lots, all_status, def_value )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_csf_new, 1, n_sq(1), 1, n_sq(2), 1, n_sq(3), 1, n_sym, 1, n_lots, all_status, def_value )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
!        RMC_MAXSCAT   = n_scat
         RMC_MAX_SQ    = n_sq
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'RMC'
         ENDIF
      ELSE                                    ! Failure
!        RMC_MAXSCAT   = n_scat
         RMC_MAX_SQ    = n_sq
         ier_typ       = 0
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'RMC'
         RETURN
      END IF
    END SUBROUTINE alloc_rmc_q
!
!
    SUBROUTINE alloc_rmc_planes ( n_planes, n_sym )
!-
!     Allocate the arrays needed by RMC
!+
      USE rmc_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_planes
      INTEGER, INTENT(IN)  :: n_sym
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      INTEGER              :: i
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( rmc_fname      ,1,n_planes  ,  all_status, ' '       )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_lambda     ,1,n_planes  ,  all_status, ' '       )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_xy, 1,4    ,1,n_planes  ,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_rlambda    ,1,n_planes  ,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_skal       ,1,n_planes  ,  all_status, 1.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_back       ,1,n_planes  ,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_chi2       ,1,n_planes  ,  all_status, 0.0D00    )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_wtot       ,1,n_planes  ,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_wic_typ    ,1,n_planes  ,  all_status, RMC_WIC_EINS)
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_nsym       ,1,n_planes  ,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_constrain  ,1,n_planes  ,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
       DO i=1,n_planes
          rmc_constrain(i) = i
       ENDDO
!
       CALL alloc_arr ( rmc_num, 1,2   ,1,n_planes  ,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_radiation  ,1,n_planes  ,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_power      ,1,n_planes  ,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr (     offq       ,1,n_planes+1,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
       CALL alloc_arr ( rmc_lxray      ,1,n_planes  ,  all_status, .true.    )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_ano        ,1,n_planes  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_ldbw       ,1,n_planes  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_eck, 1,3,1,3, 1, n_sym  ,1,n_planes  ,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( rmc_vi , 1,3,1,2, 1, n_sym  ,1,n_planes  ,  all_status, 0.0D0     )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr (     offsq      ,1,n_planes+1, 1,n_sym,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         RMC_MAX_PLANES = n_planes
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'RMC'
         ENDIF
      ELSE                                    ! Failure
         RMC_MAX_PLANES = n_planes
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'RMC'
         RETURN
      END IF
    END SUBROUTINE alloc_rmc_planes
!
!
    SUBROUTINE alloc_save ( n_scat, n_site )
!-
!     Allocate the arrays needed by SAVE
!+
      USE discus_save_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( sav_latom      ,0,n_scat  ,  all_status, .true.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( sav_t_latom    ,0,n_scat  ,  all_status, .true.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( sav_lsite      ,0,n_site  ,  all_status, .true.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         SAV_MAXSCAT   = n_scat
         SAV_MAXSITE   = n_site
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Save'
         ENDIF
      ELSE                                    ! Failure
         SAV_MAXSCAT   = n_scat
         SAV_MAXSITE   = n_site
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Save'
         RETURN
      END IF
    END SUBROUTINE alloc_save
!
!
    SUBROUTINE alloc_shear ( n_scat, n_site )
!-
!     Allocate the arrays needed by Shear
!+
      USE shear_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( shear_latom ,0,n_scat,  all_status, .false. )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( shear_lsite ,0,n_site,  all_status, .false. )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         SHEAR_MAXSCAT = n_scat
         SHEAR_MAXSITE = n_site
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Shear'
         ENDIF
      ELSE                                    ! Failure
         SHEAR_MAXSCAT = n_scat
         SHEAR_MAXSITE = n_site
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Shear'
         RETURN
      END IF
    END SUBROUTINE alloc_shear
!
!
    SUBROUTINE alloc_stack ( n_types, n_layers, l_rot )
!-
!     Allocate the arrays needed by POWDER
!+
      USE stack_mod
      USE precision_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_types
      INTEGER, INTENT(IN)  :: n_layers
      LOGICAL, INTENT(IN)  :: l_rot
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      COMPLEX (KIND=PREC_DP) :: def_value
!
      lstat      = .TRUE.
      def_value  = CMPLX(0.0D0,0.0D0,KIND=PREC_DP)
!
      CALL alloc_arr (  st_layer    ,0,n_types  ,  all_status, ' '       )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_layer_c  ,0,n_types  ,  all_status, ' '       )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_llayer   ,0,n_types  ,  all_status, 0         )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_number   ,0,n_types  ,  all_status, 0         )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_ndisp    ,0,n_types  ,  all_status, 0         )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_chem     ,0,n_types  ,  all_status, 0         )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_disp     ,1,3      ,                                           &
                                     0,n_types  ,  all_status, 0.0D0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_corr     ,0,n_types  ,                                  &
                                     0,n_types  ,  all_status, 0.0D0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_sigma    ,0,n_types  ,                                  &
                                     0,n_types  ,                                  &
                                     1,3      , all_status, 0.0D0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_trans    ,0,n_types  , &
                                     0,n_types  , &
                                     1,3      ,  all_status, 0.0D0      )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_type     ,0,n_layers,  all_status, 0         )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_internal ,0,n_types ,  all_status, .false.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr (  st_origin   ,1,3      ,                                          &
                                     0,n_layers,  all_status, 0.0D0     )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF ( l_rot ) THEN   ! Allocate rotational stacking faults, this is rarely used
!
         CALL alloc_arr (  st_rot_ang_no, 0,n_layers,  all_status, 0.0D0 )
         lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
         CALL alloc_arr (  st_rot_ang_m1, 0,n_layers,  all_status, 0.0D0 )
         lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
         CALL alloc_arr (  st_rot_ang_m2, 0,n_layers,  all_status, 0.0D0 )
         lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      ENDIF
!
      IF( lstat ) THEN                        ! Success
         ST_MAXTYPE    = n_types
         ST_MAXLAYER   = n_layers
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Stack'
         ENDIF
      ELSE                                    ! Failure
         ST_MAXTYPE    = n_types
         ST_MAXLAYER   = n_layers
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Stack'
         RETURN
      END IF
    END SUBROUTINE alloc_stack
!
!
    SUBROUTINE alloc_stack_four ( n_qxy)
!-
!     Allocate the arrays needed by POWDER
!+
      USE stack_mod
      USE precision_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, dimension(3), INTENT(IN)  :: n_qxy
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
      COMPLEX (KIND=PREC_DP) :: def_value
!
      lstat      = .TRUE.
      def_value  = CMPLX(0.0D0,0.0D0,KIND=PREC_DP)
!
!
      CALL alloc_arr (  st_csf, 1,n_qxy(1), 1,n_qxy(2), 1,n_qxy(3),  all_status, def_value )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         ST_MAXQXY     = n_qxy
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Stack_Fourier'
         ENDIF
      ELSE                                    ! Failure
         ST_MAXQXY     = n_qxy
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Stack_Fourier'
         RETURN
      END IF
    END SUBROUTINE alloc_stack_four
!
!
    SUBROUTINE alloc_stack_crystal ( n_scat, n_max )
!-
!     Allocate the arrays needed by CRYSTAL
!+
      USE stack_cr_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_max
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat  = .TRUE.
!
      CALL alloc_arr ( st_at_lis     ,0,n_scat,  all_status, ' ')
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sa_at_lis     ,0,n_scat,  all_status, ' ')
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_scat,1,11,0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_dw,          0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_anis,    1,6,0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sa_dw,          0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_occ,         0,n_scat,  all_status, 1.00D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sa_occ,         0,n_scat,  all_status, 1.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_delfr       ,0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_delfi       ,0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_mole ,       1,n_max,   all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_surf ,0, 3,  1,n_max ,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_magn ,0, 3,  1,n_max ,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sa_iscat,       1,n_scat,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sa_prop,        1,n_scat,  all_status, 0  )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_delfi       ,0,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_iscat,1,3  ,1,n_max,         all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_prop       ,1,n_max ,  all_status, 0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( st_pos,1,3,1,n_max ,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sa_pos,1,3,1,n_scat,  all_status, 0.0D0)
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
!
      IF( lstat ) THEN                        ! Success
         ST_MAX_SCAT   = n_scat
         ST_MMAX       = n_max
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Crystal'
         ENDIF
      ELSE                                    ! Failure
         ST_MAX_SCAT   =  1
         ST_MMAX       =  1
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Crystal'
         RETURN
      END IF
    END SUBROUTINE alloc_stack_crystal
!
!*******************************************************************************
!
subroutine alloc_super(n_site, n_scat, n_waves, n_groups )
!-
!  Allocate super space
!+
use superspace_mod
!
implicit none
!
integer, intent(in) :: n_site
integer, intent(in) :: n_scat
integer, intent(in) :: n_waves
integer, intent(in) :: n_groups
!
integer             :: all_status
!
call alloc_arr(sup_group,                      1, n_waves, all_status, 1     )
call alloc_arr(sup_atom, 1, 2,      1, n_site, 1, n_waves, 1, n_groups, all_status, ' '   )
!call alloc_arr(sup_repl,            1, n_site, all_status, 'VOID')
call alloc_arr(sup_irepl,           1, n_site, 1, n_waves, all_status, 0     )
call alloc_arr(sup_char,            1, n_site, 1, n_waves, all_status, 0     )
call alloc_arr(sup_func,            1, n_site, 1, n_waves, all_status, 0     )
call alloc_arr(sup_func_p,          1, n_site, 1, n_waves, all_status, 0.0D0 )
call alloc_arr(sup_old,  0, n_scat, 1, n_site, 1, n_waves, all_status, .FALSE.)
call alloc_arr(sup_new,  0, n_scat, 1, n_site, 1, n_waves, all_status, -1    )
call alloc_arr(sup_ampl, 1, 4, 1,3, 1, n_site, 1, n_waves, all_status, 0.0D0)
call alloc_arr(sup_phase,      1,3, 1, n_site, 1, n_waves, all_status, 0.0D0)
call alloc_arr(sup_prob, 1, 3,      1, n_site, 1, n_waves, all_status, 0.0D0)
call alloc_arr(sup_qvec, 1, 3,                 1, n_waves, all_status, 0.0D0)
!
sup_nwaves  = n_waves
sup_ngroups = n_groups
!
end subroutine alloc_super
!
!*******************************************************************************
!
    SUBROUTINE alloc_surf ( n_scat )
!-
!     Allocate the arrays needed by Shear
!+
      USE surface_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( surf_ex_dist    ,0,n_scat  , &
                                      all_status, SURF_DIST_DEF )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( surf_in_dist    ,0,n_scat  ,  &
                                     all_status, SURF_DIST_DEF )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         SURF_MAXSCAT  = n_scat
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Shear'
         ENDIF
      ELSE                                    ! Failure
         SURF_MAXSCAT  = n_scat
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Shear'
         RETURN
      END IF
    END SUBROUTINE alloc_surf
!
    SUBROUTINE alloc_symmetry ( n_scat, n_site )
!-
!     Allocate the arrays needed by SYMMETRY
!+
      USE symm_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
      CALL alloc_arr ( sym_latom      ,0,n_scat  ,  all_status, .false.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      CALL alloc_arr ( sym_lsite      ,0,n_site  ,  all_status, .false.   )
      lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         SYM_MAXSCAT   = n_scat
         SYM_MAXSITE   = n_site
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Symmetry'
         ENDIF
      ELSE                                    ! Failure
         SYM_MAXSCAT   = n_scat
         SYM_MAXSITE   = n_site
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Symmetry'
         RETURN
      END IF
    END SUBROUTINE alloc_symmetry
!
    SUBROUTINE alloc_transfrm ( n_scat, n_site )
!-
!     Allocate the arrays needed by TRANSFRM
!+
      USE transfrm_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( tran_latom     ,0,n_scat  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( tran_lsite     ,0,n_site  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         TRAN_MAXSCAT  = n_scat
         TRAN_MAXSITE  = n_site
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Transform'
         ENDIF
      ELSE                                    ! Failure
         TRAN_MAXSCAT  = n_scat
         TRAN_MAXSITE  = n_site
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Transform'
         RETURN
      END IF
    END SUBROUTINE alloc_transfrm
!
!
    SUBROUTINE alloc_waves ( n_scat, n_site )
!-
!     Allocate the arrays needed by WAVES
!+
      USE waves_mod
!
      IMPLICIT NONE
!
!      
      INTEGER, INTENT(IN)  :: n_scat
      INTEGER, INTENT(IN)  :: n_site
!
      INTEGER              :: all_status
      LOGICAL              :: lstat
!
      lstat     = .TRUE.
!
       CALL alloc_arr ( wv_repl      ,0,n_scat  ,  all_status, 0         )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( wv_latom     ,0,n_scat  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( wv_latom_rot ,0,n_scat  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
       CALL alloc_arr ( wv_lsite     ,0,n_site  ,  all_status, .false.   )
       lstat = lstat .and. all_status >= 0     ! This will be true if all worked out
!
      IF( lstat ) THEN                        ! Success
         WV_MAXSCAT    = n_scat
         WV_MAXSITE    = n_site
         ier_typ       = 0
         ier_num       = 0
         IF ( all_status == 1 ) THEN
            ier_typ       = 1
            ier_num       = ER_COMM
            ier_msg(1)    = 'Waves'
         ENDIF
      ELSE                                    ! Failure
         WV_MAXSCAT    = n_scat
         WV_MAXSITE    = n_site
         ier_num       = -3
         ier_typ       = ER_COMM
         ier_msg(1)    = 'Waves'
         RETURN
      END IF
    END SUBROUTINE alloc_waves
!
    SUBROUTINE dealloc_chemistry
!-
!     Deallocate the arrays for CHEMISTRY
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_chem_correlation( 1     )
      CALL alloc_chem_ang ( 1,  1        )
      CALL alloc_chem_aver( 1,  1        )
      CALL alloc_chem_disp( 1,  1        )
      CALL alloc_chem_vec ( 1,  1        )
      CALL alloc_chem_con ( 1,  1        )
      CALL alloc_chem_ran ( 1,  1,  1    )
    END SUBROUTINE dealloc_chemistry
!
    SUBROUTINE dealloc_crystal
!-
!     Deallocate the arrays for CRYSTAL
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_crystal_scat( 1)
      CALL alloc_crystal_nmax( 1)
!
    END SUBROUTINE dealloc_crystal
!
    SUBROUTINE dealloc_debye
!-
!     Deallocate the arrays for DEBYE
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      USE precision_mod
      IMPLICIT NONE
!
      INTEGER(KIND=PREC_INT_LARGE), PARAMETER :: ONE=1
      CALL alloc_debye ( 1, 1, 1, ONE )
!
    END SUBROUTINE dealloc_debye
!
    SUBROUTINE dealloc_diffuse
!-
!     Deallocate the arrays for DIFFUSE (FOURIER; PATTERSON ETC)
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_diffuse_four ( (/1,1,1/))
      CALL alloc_diffuse_scat ( 1)
      call alloc_diffuse_dbw  ( (/1,1,1/), 1)
      CALL alloc_diffuse_atom ( 1)
!
    END SUBROUTINE dealloc_diffuse
!
    SUBROUTINE dealloc_domain
!-
!     Deallocate the arrays for DIFFUSE (FOURIER; PATTERSON ETC)
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_domain   ( 1            )
!
    END SUBROUTINE dealloc_domain
!
    SUBROUTINE dealloc_micro
!-
!     Deallocate the arrays for DIFFUSE (FOURIER; PATTERSON ETC)
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_micro   ( 1, 1         )
!
    END SUBROUTINE dealloc_micro
!
    SUBROUTINE dealloc_mmc
!-
!     Deallocate the arrays for MMC
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_mmc      ( 1,  8,  1,  1)
      CALL alloc_mmc_angle( 1,  1        )
      CALL alloc_mmc_buck ( 1,  1        )
      CALL alloc_mmc_lenn ( 1,  1        )
      CALL alloc_mmc_rep  ( 1,  1        )
!
    END SUBROUTINE dealloc_mmc
!
    SUBROUTINE dealloc_molecule
!-
!     Deallocate the arrays for PDF
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_molecule ( 1, 1, 1, 1, 1 )
!
    END SUBROUTINE dealloc_molecule
!
!*******************************************************************************
!
SUBROUTINE dealloc_phases
!-
!     Deallocate the arrays for PHASES
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
IMPLICIT NONE
!
CALL alloc_phases ( 1, (/1,1,1/),  1    )
!
END SUBROUTINE dealloc_phases
!
!*******************************************************************************
!
    SUBROUTINE dealloc_pdf
!-
!     Deallocate the arrays for PDF
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_pdf ( 1, 1,  1, 1 )
!
    END SUBROUTINE dealloc_pdf
!
    SUBROUTINE dealloc_plot
!-
!     Deallocate the arrays for PLOT
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_plot ( 1, 1, 1 )
!
    END SUBROUTINE dealloc_plot
!
    SUBROUTINE dealloc_powder
!-
!     Deallocate the arrays for POWDER
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_powder ( 1, 1  )
      CALL alloc_powder_partial ( 1     )
!
    END SUBROUTINE dealloc_powder
!
    SUBROUTINE dealloc_powder_nmax
!-
!     Deallocate the arrays for POWDER
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_powder_nmax ( 1, 1 )
!
    END SUBROUTINE dealloc_powder_nmax
!
    SUBROUTINE dealloc_rmc
!-
!     Deallocate the arrays for RMC
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_rmc ( 1, 1 )
!
    END SUBROUTINE dealloc_rmc
!
    SUBROUTINE dealloc_shear
!-
!     Deallocate the arrays for SHEAR
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_shear ( 1, 1 )
!
    END SUBROUTINE dealloc_shear
!
    SUBROUTINE dealloc_stack
!-
!     Deallocate the arrays for STACK
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_stack         ( 1, 1, .TRUE.)
      CALL alloc_stack_four    ( (/1,1,1/)  )
      CALL alloc_stack_crystal ( 1,  1        )
!
    END SUBROUTINE dealloc_stack
!
!*******************************************************************************
!
SUBROUTINE dealloc_super
!-
!     Deallocate the arrays for SUPER SPACE
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_super ( 1, 1, 1, 1 )
!
END SUBROUTINE dealloc_super
!
!*******************************************************************************
!
    SUBROUTINE dealloc_surf
!-
!     Deallocate the arrays for SURFACE
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_surf ( 1 )
!
    END SUBROUTINE dealloc_surf
!
    SUBROUTINE dealloc_save
!-
!     Deallocate the arrays for SAVE
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_save ( 1, 1 )
!
    END SUBROUTINE dealloc_save
!
    SUBROUTINE dealloc_symmetry
!-
!     Deallocate the arrays for SYMMETRY
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_symmetry ( 1, 1 )
!
    END SUBROUTINE dealloc_symmetry
!
    SUBROUTINE dealloc_transfrm
!-
!     Deallocate the arrays for TRANSFRM
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_transfrm ( 1, 1 )
!
    END SUBROUTINE dealloc_transfrm
!
    SUBROUTINE dealloc_waves
!-
!     Deallocate the arrays for WAVES
!     To avoid possible pitfals with old code, the arrays are simply
!     reallocated to a size of 1.
!+
      IMPLICIT NONE
!
      CALL alloc_waves ( 1,  1 )
!
    END SUBROUTINE dealloc_waves
!
!
END MODULE discus_allocate_appl_mod
