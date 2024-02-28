MODULE gen_hdf_write_mod
!
contains
!
!*****7*****************************************************************
!
SUBROUTINE gen_hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                       extr_abs, extr_ord, extr_top,                       &
                       cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax, &
                       ier_num, ier_typ, ER_IO, ER_APPL)

!USE hdf5
!use diffuse_mod
!use fourier_sup
use precision_mod
!
IMPLICIT NONE
!
!INTEGER, PARAMETER:: PREC_DP=SELECTED_REAL_KIND(p=15,r=307)  ! double precision
!
INTEGER, INTENT(IN) :: value
LOGICAL, INTENT(IN) :: laver
CHARACTER(LEN=200), INTENT(IN) :: outfile
INTEGER, DIMENSION(3)  , INTENT(IN) :: out_inc
REAL(kind=PREC_DP)   , DIMENSION(3,4), INTENT(IN) :: out_eck ! (3,4)
REAL(kind=PREC_DP)   , DIMENSION(3,3), INTENT(IN) :: out_vi
integer                              , intent(in) :: extr_abs
integer                              , intent(in) :: extr_ord
integer                              , intent(in) :: extr_top
REAL(kind=PREC_DP)   , DIMENSION(3)  , INTENT(IN) :: cr_a0
REAL(kind=PREC_DP)   , DIMENSION(3)  , INTENT(IN) :: cr_win
REAL(kind=PREC_DP)   , DIMENSION(out_inc(1), out_inc(2), out_inc(3)), INTENT(IN) :: qvalues
INTEGER                , INTENT(IN) :: VAL_PDF
INTEGER                , INTENT(IN) :: VAL_3DPDF
REAL(KIND=PREC_DP)     , INTENT(IN)  :: valmax
INTEGER                , INTENT(OUT) :: ier_num
INTEGER                , INTENT(OUT) :: ier_typ
INTEGER                , INTENT(IN) :: ER_IO
INTEGER                , INTENT(IN) :: ER_APPL
!
ier_num = -172  ! HDF5 not supported
ier_typ = ER_APPL
!
END SUBROUTINE gen_hdf5_write
!
!*****7*****************************************************************
!
end MODULE gen_hdf_write_mod
