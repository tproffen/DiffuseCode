MODULE kuplot_load_h5
!-
!  Contains routines to read a HDF5 file written by DISCUS
!+
!USE hdf5
USE prompt_mod
USE precision_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC hdf5_read_kuplot
!
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE hdf5_read_kuplot(infile, length, O_LAYER, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte, iz, ku_ndims,     &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!
!SUBROUTINE hdf5_read_kuplot(infile, length, O_LAYER, NOPTIONAL, opara, lopara,         &
!                     lpresent, owerte,               &
!                     MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
!                     xmin, xmax, ymin, ymax, offxy, offz, lni, lh5, lenc,       &
!                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!
use kuplot_config
use kuplot_global
use lib_data_struc_h5
!
IMPLICIT NONE
!
CHARACTER(LEN=1024), INTENT(IN) :: infile
INTEGER            , INTENT(IN) :: length
INTEGER            , INTENT(IN) :: O_LAYER
INTEGER            , INTENT(IN) :: NOPTIONAL
CHARACTER(LEN=*)   , DIMENSION(NOPTIONAL), INTENT(IN) :: opara
INTEGER            , DIMENSION(NOPTIONAL), INTENT(IN) :: lopara
LOGICAL            , DIMENSION(NOPTIONAL), INTENT(IN) :: lpresent
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL), INTENT(IN) :: owerte
INTEGER                                  , INTENT(INOUT) :: iz     ! KUPLOT data set number
INTEGER, DIMENSION(  MAXKURVTOT)         , INTENT(INOUT) :: ku_ndims
!
INTEGER,                            INTENT(OUT)   :: ier_num
INTEGER,                            INTENT(OUT)   :: ier_typ
INTEGER,                            INTENT(IN )   :: idims
CHARACTER(LEN=*), DIMENSION(idims), INTENT(INOUT) :: ier_msg    ! Error message
INTEGER,                            INTENT(IN )   :: ER_APPL
INTEGER,                            INTENT(IN )   :: ER_IO
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!
ier_num = -72   ! HDF5 not supported
ier_typ = ER_APPL !
!
END SUBROUTINE hdf5_read_kuplot
!
!*******************************************************************************
!
!
!*******************************************************************************
!
END MODULE kuplot_load_h5
