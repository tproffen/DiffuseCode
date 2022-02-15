!
!*****7*****************************************************************
!
SUBROUTINE hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                       cr_a0, cr_win, VAL_PDF, VAL_3DPDF,               &
                       ier_num, ier_typ, ER_IO, ER_APPL)
!
use precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: value
LOGICAL, INTENT(IN) :: laver
CHARACTER(LEN=200), INTENT(IN) :: outfile
INTEGER, DIMENSION(3)  , INTENT(IN) ::  out_inc
REAL(kind=PREC_DP)   , DIMENSION(3,4), INTENT(IN) ::  out_eck ! (3,4)
REAL(kind=PREC_DP)   , DIMENSION(3,3), INTENT(IN) ::  out_vi 
REAL(kind=PREC_DP)   , DIMENSION(3)  , INTENT(IN) ::  cr_a0
REAL(kind=PREC_DP)   , DIMENSION(3)  , INTENT(IN) ::  cr_win
INTEGER                , INTENT(IN) :: VAL_PDF
INTEGER                , INTENT(IN) :: VAL_3DPDF
!
INTEGER,                 INTENT(OUT) :: ier_num
INTEGER,                 INTENT(OUT) :: ier_typ
INTEGER,                 INTENT(IN) :: ER_IO
INTEGER,                 INTENT(IN) :: ER_APPL
!
ier_num = -172  ! HDF5 not supported
ier_typ = ER_APPL
!
END SUBROUTINE hdf5_write
!
!*****7*****************************************************************
!
