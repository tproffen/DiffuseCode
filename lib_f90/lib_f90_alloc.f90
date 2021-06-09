module lib_f90_allocate_mod
!
use allocate_generic
use errlist_mod
!
contains
!
   subroutine alloc_param(n_res)
!
   use param_mod
!
   implicit none
!
   integer, intent(in)  :: n_res
   integer              :: all_status
   integer              :: size_of
!
   call alloc_arr(res_para, 0, n_res, all_status, 0.0D0, size_of )
!
   end subroutine alloc_param
!
!*******************************************************************************
!
   subroutine alloc_ref_para(n_para)
!
   use param_mod
!
   implicit none
!
   integer, intent(in)  :: n_para
   integer              :: all_status
   integer              :: size_of
!
   call alloc_arr(ref_para, 0, n_para, all_status, 0.0, size_of )
   MAXPAR_REF = n_para
!
   end subroutine alloc_ref_para
!
!*******************************************************************************
!
subroutine alloc_expr(n_expr)
!-
!  Allocate the array for expressions
!+
!
use variable_mod
!
implicit none
!
integer, intent(in)  :: n_expr
integer              :: all_status
integer              :: size_of
!
call alloc_arr(var_expr, 1, n_expr, all_status, ' ', size_of)
call alloc_arr(var_expr_val, 1, n_expr, all_status, 0.0D0, size_of)
var_entry(VAR_EXPRESSION) = n_expr    ! Store dimension
!
end subroutine alloc_expr
!
!*******************************************************************************
!
end module lib_f90_allocate_mod
