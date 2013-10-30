module discus_wrap
use iso_c_binding
use config_mod
use crystal_mod
use pdf_mod
implicit none

contains
  
  subroutine get_cr_pos_c(cr_pos_out,n) bind(c)
    real(c_float), intent (out), dimension(3,n) :: cr_pos_out
    integer(c_int), intent(in), value :: n
    cr_pos_out=cr_pos
  end subroutine get_cr_pos_c
  
  subroutine set_cr_pos_c(cr_pos_in,n) bind(c)
    real(c_float), intent (in), dimension(3,n) :: cr_pos_in
    integer(c_int), intent(in), value :: n
    cr_pos=cr_pos_in
  end subroutine set_cr_pos_c
  
  subroutine get_pdf_c(out,nmi,n) bind(c)
    real(c_double), intent(out), dimension(n) :: out
    integer(c_int), intent(in), value :: nmi, n
    integer :: i,j
    j=0
    do i=nmi,nmi+n-1
       j=j+1
       out(j)=pdf_skal * pdf_calc(i)
    end do
  end subroutine get_pdf_c
  
end module discus_wrap
