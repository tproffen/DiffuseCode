MODULE fit_set_sub_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE fit_set_sub
!
USE set_sub_generic_mod
!
IMPLICIT NONE
!
INTERFACE
   SUBROUTINE    fit_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE    fit_mache_kdo
END INTERFACE
!
p_mache_kdo         => fit_mache_kdo
!
END SUBROUTINE fit_set_sub
!
!*******************************************************************************
!
END MODULE fit_set_sub_mod
