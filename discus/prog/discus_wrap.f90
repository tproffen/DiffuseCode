MODULE discus_wrap
USE iso_c_binding
USE config_mod
USE crystal_mod
implicit none

contains
  
  SUBROUTINE get_cr_pos_c(cr_pos_out,n) BIND(C)
    REAL(c_float), BIND(C), intent (out), dimension(3,n) :: cr_pos_out
    INTEGER(c_int), intent(in), value :: n
    cr_pos_out=cr_pos
  END SUBROUTINE get_cr_pos_c
  
  !SUBROUTINE set_cr_pos() BIND(C, name='set_cr_pos_c')
    
  !END SUBROUTINE set_cr_pos
  
  
END MODULE discus_wrap
