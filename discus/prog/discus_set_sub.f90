MODULE discus_setup_sub_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE discus_set_sub
!
! Sets the specific DISCUS interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
use discus_errlist_mod
!use discus_mache_kdo_mod
use discus_update_mod
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE discus_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE discus_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE discus_branch(zeile, length, lreset)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
   END SUBROUTINE discus_branch
END INTERFACE
!
p_mache_kdo         => discus_mache_kdo
p_errlist_appl      => discus_errlist_appl
p_ersetz_para       => discus_ersetz_para
p_upd_para          => discus_upd_para
p_calc_intr_spec    => discus_calc_intr_spec
p_calc_intr_log_spec=> discus_calc_intr_log_spec
p_validate_var_spec => discus_validate_var_spec
p_branch            => discus_branch
!p_loop_mpi          => dummy_loop_mpi
p_get_var_type      => discus_get_var_type
!
END SUBROUTINE discus_set_sub
!
!*******************************************************************************
!
END MODULE discus_setup_sub_mod
