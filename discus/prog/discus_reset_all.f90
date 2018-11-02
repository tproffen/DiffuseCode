MODULE discus_reset_all_mod
!
CONTAINS
!
SUBROUTINE discus_reset_all
!
! Put all arrays and all menues into program startup mode
!
USE chem_menu
USE class_internal
USE crystal_mod
USE conn_mod
USE discus_allocate_appl_mod
USE discus_init_mod
USE discus_plot_menu
USE domain_menu
USE fourier_reset_mod
USE insert_menu
USE mole_surf_mod
USE output_menu
USE pdf_menu
USE powder
USE prop_para_func
USE rmc_menu
USE save_menu
USE shear
USE stack_menu
USE symm_sup_mod
USE structur
USE transform_menu
USE waves_do_menu
!
IMPLICIT NONE
!
INTEGER, PARAMETER  :: code_res   = -2
!
CHARACTER(LEN=1024) :: zeile
INTEGER             :: lp
!
CALL discus_alloc_default               ! Shrink all arrays to minimum size
CALL discus_initarrays                  ! Reset most arrays, does mmc_init
CALL rese_cr                            ! Clean crystal
CALL deco_reset                         ! Surface decoration
CALL chem_reset                         ! Complete chemisty menu
zeile = ' '
lp    = 1
CALL conn_do_set(code_res,zeile, lp)    ! Connectivity
CALL domain_reset                       ! Domain menu
CALL fourier_reset                      ! Fourier/diffuse variables 
CALL insert_reset                       ! Insertion of molecules 
CALL output_reset                       ! Output of Fourier results
CALL pdf_reset                          ! PDF reset
CALL plot_reset                         ! plot_reset
CALL powder_reset                       ! reset powder menu
zeile = 'ignore,all'
lp    = 10
CALL property_select (zeile, lp, cr_sel_prop)
CALL rmc_reset                          ! RMC  reset
CALL save_reset                         ! Save menu
CALL shear_reset                        ! Shear menu
CALL stack_reset                        ! Stack menu
CALL store_remove_all(store_root)       ! Internal structure storage
CALL symm_reset                         ! Symmetry menu
CALL tran_reset                         ! Transform metric menu
CALL waves_reset                        ! Modulation waves menu
!
END SUBROUTINE discus_reset_all
!
END MODULE discus_reset_all_mod
