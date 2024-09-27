MODULE kuplot_load_h5
!-
!  Contains routines to read a HDF5 file written by DISCUS
!+
!USE hdf5
!
USE kuplot_config
!
use lib_hdf5_read_mod
use lib_data_struc_h5
use hdf5_def_mod
USE precision_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC hdf5_read_kuplot
PUBLIC nexus_read_kuplot
!PUBLIC hdf5_place_kuplot
!
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE hdf5_read_kuplot(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte, iz, ku_ndims,     &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, lout, output_io)
!
use kuplot_global
!
USE ber_params_mod
use lib_data_struc_h5
!
IMPLICIT NONE
!
CHARACTER(LEN=1024)                      , INTENT(INOUT) :: infile
INTEGER                                  , INTENT(IN) :: length
INTEGER                                  , INTENT(IN) :: O_LAYER
INTEGER                                  , INTENT(IN) :: O_TRANS
INTEGER                                  , INTENT(IN) :: NOPTIONAL
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
logical,                            intent(in)    :: lout
INTEGER,                            INTENT(IN )   :: ER_IO
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!
CHARACTER(LEN=14)   :: dataname    ! Dummy name for HDF5 datasets
!
integer               :: node_number = 0
integer               :: ndims = 0
integer               :: ik
integer, dimension(3) :: dims  = 1
!
dataname = ' '
!

call hdf5_read(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number,ndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
call dgl5_set_h5_is_ku(iz, node_number)
call dgl5_set_ku_is_h5(node_number, iz)
if(ier_num/=0) return
!
ku_ndims(iz) = ndims
ik = iz 
call dgl5_set_ku_is_h5(iz, node_number)
call dgl5_set_h5_is_ku(node_number, iz)
!
call data2kuplot(ik, infile, lout)
!
END SUBROUTINE hdf5_read_kuplot
!
!*******************************************************************************
!
!
!*******************************************************************************
!
SUBROUTINE nexus_read_kuplot(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte, iz, ku_ndims,     &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, lout, output_io)
!-
!  Read a Nexus file format specified as teh common develoepers format
!+
!
use kuplot_global
!
USE ber_params_mod
use lib_data_struc_h5
use lib_nx_read_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1024)                      , INTENT(INOUT) :: infile
INTEGER                                  , INTENT(IN) :: length
INTEGER                                  , INTENT(IN) :: O_LAYER
INTEGER                                  , INTENT(IN) :: O_TRANS
INTEGER                                  , INTENT(IN) :: NOPTIONAL
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
logical,                            intent(in)    :: lout
INTEGER,                            INTENT(IN )   :: ER_IO
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!
CHARACTER(LEN=14)   :: dataname    ! Dummy name for HDF5 datasets
!
integer               :: node_number = 0
integer               :: ndims = 0
integer               :: ik
integer, dimension(3) :: dims  = 1
!
dataname = ' '
!

call nx_read_scattering(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number,ndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
call dgl5_set_h5_is_ku(iz, node_number)
call dgl5_set_ku_is_h5(node_number, iz)
if(ier_num/=0) return
!
ku_ndims(iz) = ndims
ik = iz 
call dgl5_set_ku_is_h5(iz, node_number)
call dgl5_set_h5_is_ku(node_number, iz)
!
call data2kuplot(ik, infile, lout)
!
END SUBROUTINE nexus_read_kuplot
!
!*******************************************************************************
!
END MODULE kuplot_load_h5
