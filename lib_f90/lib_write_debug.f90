module lib_write_mod
!
!  A few simple write routines, useful for debug
!+
!
interface tofile
  module procedure complex_tofile, real_tofile, real_tonipl
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
!*******************************************************************************
!
subroutine real_tonipl(idims, cname, field, xmin, xstep)
!-
!  Real Version 2D field to NIPL type output
!+
!
use precision_mod
!
integer           , dimension(2)      , intent(in) :: idims
character(len=*)                      , intent(in) :: cname
real(kind=PREC_DP), dimension(idims(1), idims(2)), intent(in) :: field
real(kind=PREC_DP), dimension(2)      , intent(in) :: xmin
real(kind=PREC_DP), dimension(2)      , intent(in) :: xstep
!
integer :: i, j
!
!
open(unit=66, file=cname, status='unknown')
write(66, '(2i10)') idims
write(66, '(4(f15.4,2x))') xmin(1), xmin(1) + (idims(1)-1)*xstep(1), &
                           xmin(2), xmin(2) + (idims(2)-1)*xstep(2)
do j=1, idims(2)
   write(66, '(5(g15.8,2x))') (field(i,j), i=1, idims(1))
enddo
close(unit=66)
!
end subroutine real_tonipl
!
!###############################################################################
!
end module lib_write_mod
