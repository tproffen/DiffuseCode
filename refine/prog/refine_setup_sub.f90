MODULE refine_setup_sub_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_set_sub
!
! Sets the specific REFINE interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
use refine_errlist_mod
!use refine_mache_kdo_mod
use refine_update_mod
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE refine_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE refine_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE refine_branch(zeile, length, lreset)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
   END SUBROUTINE refine_branch
END INTERFACE
!
!INTERFACE
!   SUBROUTINE refine_loop_mpi(prog_n, prog_l, mac_n, mac_l, out_n, out_l, repeat, nindiv)
!!
!CHARACTER (LEN=*), INTENT(IN) :: prog_n
!CHARACTER (LEN=*), INTENT(IN) :: mac_n
!CHARACTER (LEN=*), INTENT(IN) :: out_n
!INTEGER          , INTENT(IN) :: prog_l
!INTEGER          , INTENT(IN) :: mac_l
!INTEGER          , INTENT(IN) :: out_l
!LOGICAL          , INTENT(IN) :: repeat
!INTEGER          , INTENT(IN) :: nindiv
!!
!   END SUBROUTINE refine_loop_mpi
!END INTERFACE
!
p_mache_kdo         => refine_mache_kdo
p_errlist_appl      => refine_errlist_appl
p_ersetz_para       => refine_ersetz_para
p_upd_para          => refine_upd_para
p_calc_intr_spec    => refine_calc_intr_spec
p_calc_intr_log_spec=> refine_calc_intr_log_spec
p_validate_var_spec => refine_validate_var_spec
p_branch            => refine_branch
!p_loop_mpi          => refine_loop_mpi
p_get_var_type      => refine_get_var_type
!
END SUBROUTINE refine_set_sub
!
!
END MODULE refine_setup_sub_mod
