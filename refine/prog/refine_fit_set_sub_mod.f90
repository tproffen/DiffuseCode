MODULE refine_fit_set_sub_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_fit_set_sub
!
USE set_sub_generic_mod
!
IMPLICIT NONE
!
INTERFACE
   SUBROUTINE    refine_fit_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE    refine_fit_mache_kdo
END INTERFACE
!
p_mache_kdo         => refine_fit_mache_kdo
!
END SUBROUTINE refine_fit_set_sub
!
!*******************************************************************************
!
END MODULE refine_fit_set_sub_mod
