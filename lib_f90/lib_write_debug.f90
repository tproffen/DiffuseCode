module lib_write_mod
!
!  A few simple write routines, useful for debug
!+
!
interface tofile
  module procedure complex_tofile, real_tofile, &
                   real_tonipl, complex_tonipl, &
                   real_tohdf5, complex_tohdf5
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
!*******************************************************************************
!
subroutine complex_tonipl(idims, cname, field, xmin, xstep)
!-
!  Real Version 2D field to NIPL type output
!+
!
use precision_mod
!
integer           , dimension(2)      , intent(in) :: idims
character(len=*)                      , intent(in) :: cname
complex(kind=PREC_DP), dimension(idims(1), idims(2)), intent(in) :: field
real(kind=PREC_DP), dimension(2)      , intent(in) :: xmin
real(kind=PREC_DP), dimension(2)      , intent(in) :: xstep
!
character(len=PREC_STRING)                         :: string
integer :: i, j
!
!
string = cname(1:len_trim(cname))//'.real'
open(unit=66, file=string, status='unknown')
write(66, '(2i10)') idims
write(66, '(4(f15.4,2x))') xmin(1), xmin(1) + (idims(1)-1)*xstep(1), &
                           xmin(2), xmin(2) + (idims(2)-1)*xstep(2)
do j=1, idims(2)
   write(66, '(5(g15.8,2x))') (real(field(i,j)), i=1, idims(1))
enddo
close(unit=66)
!
string = cname(1:len_trim(cname))//'.imag'
open(unit=66, file=string, status='unknown')
write(66, '(2i10)') idims
write(66, '(4(f15.4,2x))') xmin(1), xmin(1) + (idims(1)-1)*xstep(1), &
                           xmin(2), xmin(2) + (idims(2)-1)*xstep(2)
do j=1, idims(2)
   write(66, '(5(g15.8,2x))') (imag(field(i,j)), i=1, idims(1))
enddo
close(unit=66)
!
end subroutine complex_tonipl
!
!###############################################################################
!
subroutine real_tohdf5(idims, cname, field, xmin, xstep)
!SUBROUTINE real_tohdf5(value, laver, outfile, out_inc, out_eck, out_vi, &
!                       extr_abs, extr_ord, extr_top,                       &
!                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
!                       ier_num, ier_typ, ER_IO, ER_APPL)
!
use gen_hdf_write_mod
use precision_mod
!
IMPLICIT NONE
integer           , dimension(3)      , intent(in) :: idims
character(len=*)                      , intent(in) :: cname
real(kind=PREC_DP), dimension(idims(1), idims(2),idims(3)), intent(in) :: field
real(kind=PREC_DP), dimension(3)      , intent(in) :: xmin
real(kind=PREC_DP), dimension(3)      , intent(in) :: xstep
!
!INTEGER, PARAMETER:: PREC_DP=SELECTED_REAL_KIND(p=15,r=307)  ! double precision
!
INTEGER  :: value
LOGICAL  :: laver
CHARACTER(LEN=200)  :: outfile
INTEGER, DIMENSION(3)    :: out_inc
REAL(kind=PREC_DP)   , DIMENSION(3,4)  :: out_eck ! (3,4)
REAL(kind=PREC_DP)   , DIMENSION(3,3)  :: out_vi
integer                               :: extr_abs
integer                               :: extr_ord
integer                               :: extr_top
REAL(kind=PREC_DP)   , DIMENSION(3)    :: cr_a0
REAL(kind=PREC_DP)   , DIMENSION(3)    :: cr_win
REAL(kind=PREC_DP)   , DIMENSION(:,:,:), allocatable :: qvalues
INTEGER                  :: VAL_PDF
INTEGER                  :: VAL_3DPDF
REAL(KIND=PREC_DP)        :: valmax
INTEGER                 :: ier_num
INTEGER                 :: ier_typ
INTEGER                  :: ER_IO
INTEGER                  :: ER_APPL
!
value    = 1
laver    = .false.
outfile  = cname
val_pdf  = 2
val_3dpdf = 3
out_inc  = idims
out_eck  = 0.0D0
out_eck(1,2) = real(idims(1), kind=PREC_DP)
out_eck(2,3) = real(idims(2), kind=PREC_DP)
out_eck(3,4) = real(idims(3), kind=PREC_DP)
out_vi   = 0.0D0
out_vi(1,1)    = 1.0D0
out_vi(2,2)    = 1.0D0
out_vi(3,3)    = 1.0D0
cr_a0    = 1.0D0
cr_win   =90.0D0
extr_abs = 1
extr_ord = 2
extr_top = 3
ER_IO = 1
ER_APPL = 1
valmax  = 0.0D0
!
allocate(qvalues(out_inc(1), out_inc(2), out_inc(3)))
qvalues = field
!
call gen_hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                       extr_abs, extr_ord, extr_top,                       &
                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
                       ier_num, ier_typ, ER_IO, ER_APPL)
deallocate(qvalues)
!
end subroutine real_tohdf5
!
!###############################################################################
!
subroutine complex_tohdf5(idims, cname, field, xmin, xstep)
!SUBROUTINE real_tohdf5(value, laver, outfile, out_inc, out_eck, out_vi, &
!                       extr_abs, extr_ord, extr_top,                       &
!                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
!                       ier_num, ier_typ, ER_IO, ER_APPL)
!
use gen_hdf_write_mod
use precision_mod
!
IMPLICIT NONE
integer           , dimension(3)      , intent(in) :: idims
character(len=*)                      , intent(in) :: cname
complex(kind=PREC_DP), dimension(idims(1), idims(2),idims(3)), intent(in) :: field
real(kind=PREC_DP), dimension(3)      , intent(in) :: xmin
real(kind=PREC_DP), dimension(3)      , intent(in) :: xstep
!
!INTEGER, PARAMETER:: PREC_DP=SELECTED_REAL_KIND(p=15,r=307)  ! double precision
!
INTEGER  :: value
LOGICAL  :: laver
CHARACTER(LEN=200)  :: outfile
INTEGER, DIMENSION(3)    :: out_inc
REAL(kind=PREC_DP)   , DIMENSION(3,4)  :: out_eck ! (3,4)
REAL(kind=PREC_DP)   , DIMENSION(3,3)  :: out_vi
integer                               :: extr_abs
integer                               :: extr_ord
integer                               :: extr_top
REAL(kind=PREC_DP)   , DIMENSION(3)    :: cr_a0
REAL(kind=PREC_DP)   , DIMENSION(3)    :: cr_win
REAL(kind=PREC_DP)   , DIMENSION(:,:,:), allocatable :: qvalues
INTEGER                  :: VAL_PDF
INTEGER                  :: VAL_3DPDF
REAL(KIND=PREC_DP)        :: valmax
INTEGER                 :: ier_num
INTEGER                 :: ier_typ
INTEGER                  :: ER_IO
INTEGER                  :: ER_APPL
!
value    = 1
laver    = .false.
outfile  = cname(1:len_trim(cname))//'.real'
val_pdf  = 2
val_3dpdf = 3
out_inc  = idims
out_eck  = 0.0D0
out_eck(1,2) = real(idims(1), kind=PREC_DP)
out_eck(2,3) = real(idims(2), kind=PREC_DP)
out_eck(3,4) = real(idims(3), kind=PREC_DP)
out_vi   = 0.0D0
out_vi(1,1)    = 1.0D0
out_vi(2,2)    = 1.0D0
out_vi(3,3)    = 1.0D0
cr_a0    = 1.0D0
cr_win   =90.0D0
extr_abs = 1
extr_ord = 2
extr_top = 3
ER_IO = 1
ER_APPL = 1
valmax  = 0.0D0
!
allocate(qvalues(out_inc(1), out_inc(2), out_inc(3)))
!
qvalues = real(field,kind=PREC_DP)
!
call gen_hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                       extr_abs, extr_ord, extr_top,                       &
                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
                       ier_num, ier_typ, ER_IO, ER_APPL)
!
outfile  = cname(1:len_trim(cname))//'.imag'
qvalues = imag(field)
!
call gen_hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                       extr_abs, extr_ord, extr_top,                       &
                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
                       ier_num, ier_typ, ER_IO, ER_APPL)
deallocate(qvalues)
!
end subroutine complex_tohdf5
!
!###############################################################################
!
end module lib_write_mod
