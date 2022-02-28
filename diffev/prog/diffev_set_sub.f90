MODULE diffev_setup_sub_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE diffev_set_sub
!
! Sets the specific DIFFEV interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
use diffev_errlist_mod
!use diffev_mache_kdo_mod
use diffev_update_mod
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE diffev_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE diffev_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_branch(zeile, length, lreset)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
   END SUBROUTINE diffev_branch
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_loop_mpi(prog_n, prog_l, mac_n, mac_l, out_n, out_l, repeat, nindiv)
!
CHARACTER (LEN=*), INTENT(IN) :: prog_n
CHARACTER (LEN=*), INTENT(IN) :: mac_n
CHARACTER (LEN=*), INTENT(IN) :: out_n
INTEGER          , INTENT(IN) :: prog_l
INTEGER          , INTENT(IN) :: mac_l
INTEGER          , INTENT(IN) :: out_l
LOGICAL          , INTENT(IN) :: repeat
INTEGER          , INTENT(IN) :: nindiv
!
   END SUBROUTINE diffev_loop_mpi
END INTERFACE
!
p_mache_kdo         => diffev_mache_kdo
p_errlist_appl      => diffev_errlist_appl
p_ersetz_para       => diffev_ersetz_para
p_upd_para          => diffev_upd_para
p_calc_intr_spec    => diffev_calc_intr_spec
p_calc_intr_log_spec=> diffev_calc_intr_log_spec
p_validate_var_spec => diffev_validate_var_spec
p_branch            => diffev_branch
p_loop_mpi          => diffev_loop_mpi
p_get_var_type      => diffev_get_var_type
!
END SUBROUTINE diffev_set_sub
!
!*****7*****************************************************************
!
END MODULE diffev_setup_sub_mod
