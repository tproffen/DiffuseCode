module lib_write_mod
!
!  A few simple write routines, useful for debug
!+
!
interface tofile
  module procedure complex_tofile, real_tofile
end interface tofile
!
private
public tofile
!
contains
!
!###############################################################################
!
subroutine complex_tofile(idims, cname, field, xmin, xstep)
!-
!  Complex version writes x, real, imag, 0.0
!+
!
use precision_mod
!
integer                                  , intent(in) :: idims
character(len=*)                         , intent(in) :: cname
complex(kind=PREC_DP), dimension(1:idims), intent(in) :: field
real(kind=PREC_DP)                       , intent(in) :: xmin 
real(kind=PREC_DP)                       , intent(in) :: xstep
!
integer :: i
!
open(unit=66, file=cname, status='unknown')
do i=1, idims
  write(66, '(4(g20.8e3,2x))') (i-1)*xstep+xmin, real(field(i)), imag(field(i)), 0.0D0
enddo
close(unit=66)
!
end subroutine complex_tofile
!
!*******************************************************************************
!
subroutine real_tofile(idims, cname, field, xmin, xstep)
!-
!  Complex version writes x, real, 0.0 , 0.0
!+
!
use precision_mod
!
integer                               , intent(in) :: idims
character(len=*)                      , intent(in) :: cname
real(kind=PREC_DP), dimension(1:idims), intent(in) :: field
real(kind=PREC_DP)                    , intent(in) :: xmin 
real(kind=PREC_DP)                    , intent(in) :: xstep
!
integer :: i
!
open(unit=66, file=cname, status='unknown')
do i=1, idims
  write(66, '(4(g20.8e3,2x))') (i-1)*xstep+xmin, field(i), 0.0D0, 0.0D0
enddo
close(unit=66)
!
end subroutine real_tofile
!
!###############################################################################
!
end module lib_write_mod
