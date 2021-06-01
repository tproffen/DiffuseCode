module lib_f90_default_mod
!
use lib_f90_allocate_mod
use param_mod
!
implicit none
!
contains
!
!*******************************************************************************
!
subroutine lib_alloc_default
!
integer :: n_res
integer :: n_para
integer :: n_expr
!
n_res = MAX(MAXPAR_RES, 6000)
call alloc_param(n_res)
MAXPAR_RES = n_res
n_para = 1
call alloc_ref_para(n_para)
n_expr = 1
call alloc_expr(n_expr)
!
end subroutine lib_alloc_default
!
!*******************************************************************************
!
end module lib_f90_default_mod
