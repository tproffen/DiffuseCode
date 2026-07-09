module kuplot_load_rmc
!-
!  Load an ASCII *.dat file in the rmc data file format
!  np 1
!  i j k h k l inte
!+
!
private
public rmc_read_kuplot
!
contains
!
!*******************************************************************************
!
subroutine rmc_read_kuplot(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte, iz, ku_ndims,     &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, lout, output_io)
!-
! Read the data
!+
use kuplot_config
!use kuplot_mod
use kuplot_global
use lib_data_struc_h5
use lib_load_mod
use precision_mod
!
implicit none
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
character(len=14)   :: dataname    ! Dummy name for HDF5 datasets
integer             :: ik
integer             :: node_number = 0
!
dataname = ' '
!
call gen_load_rmc(infile, node_number, lout)
call dgl5_set_ku_is_h5(iz, node_number)
call dgl5_set_h5_is_ku(node_number, iz)
ku_ndims(iz) = 3
ik = iz
!
call data2kuplot(ik, infile  , lout  )
!
end subroutine rmc_read_kuplot
!
!*******************************************************************************
!
end module kuplot_load_rmc
