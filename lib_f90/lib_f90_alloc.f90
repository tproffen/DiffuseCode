MODULE lib_f90_allocate_mod
!
USE allocate_generic
USE errlist_mod
!
CONTAINS
!
   SUBROUTINE alloc_param(n_res)
!
   USE param_mod
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)  :: n_res
   INTEGER              :: all_status
   INTEGER              :: size_of
!
   CALL alloc_arr(res_para, 0, n_res, all_status, 0.0, size_of )
!
   END SUBROUTINE alloc_param
END MODULE lib_f90_allocate_mod
