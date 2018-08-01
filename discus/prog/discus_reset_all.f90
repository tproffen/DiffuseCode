MODULE discus_reset_all_mod
!
CONTAINS
!
SUBROUTINE discus_reset_all
!
! Put all arrays and all menues into program startup mode
!
USE class_internal
USE conn_mod
USE discus_plot_menu
USE structur
!
IMPLICIT NONE
!
INTEGER, PARAMETER  :: code_res   = -2
!
CHARACTER(LEN=1024) :: zeile
INTEGER             :: lp
!
CALL rese_cr                            ! Clean crystal
CALL store_remove_all(store_root)
zeile = ' '
lp    = 1
CALL conn_do_set(code_res,zeile, lp)    ! Connectivity
CALL plot_reset                         ! plot_reset
!
END SUBROUTINE discus_reset_all
!
END MODULE discus_reset_all_mod
