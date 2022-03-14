MODULE kuplot_setup_sub_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE kuplot_set_sub
!
! Sets the specific DIFFEV interfaces four routines that are refecenced in
! LIB_F90 by their generic names
!
use kuplot_errlist_mod
use kuplot_update_mod
USE set_sub_generic_mod
USE prompt_mod
!
INTERFACE
   SUBROUTINE kuplot_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE kuplot_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_branch(zeile, length, lreset, lloop)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
integer          , INTENT(IN) :: lloop
!
   END SUBROUTINE kuplot_branch
END INTERFACE
!
INTERFACE
   SUBROUTINE kuplot_top(zeile)
!
   CHARACTER (LEN=*), INTENT(IN) :: zeile
   END SUBROUTINE kuplot_top
END INTERFACE
!
p_mache_kdo         => kuplot_mache_kdo
p_errlist_appl      => kuplot_errlist_appl
p_ersetz_para       => kuplot_ersetz_para
p_upd_para          => kuplot_upd_para
p_calc_intr_spec    => kuplot_calc_intr_spec
p_calc_intr_log_spec=> kuplot_calc_intr_log_spec
p_validate_var_spec => kuplot_validate_var_spec
p_branch            => kuplot_branch
!p_loop_mpi          => dummy_loop_mpi
p_get_var_type      => kuplot_get_var_type
IF(lstandalone) THEN
p_top               => kuplot_top
ENDIF
!
END SUBROUTINE kuplot_set_sub
!
!*****7*****************************************************************
!
END MODULE kuplot_setup_sub_mod
