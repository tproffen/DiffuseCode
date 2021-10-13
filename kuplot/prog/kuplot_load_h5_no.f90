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
PUBLIC hdf5_read
PUBLIC hdf5_place_kuplot
PUBLIC hdf5_get_layer
PUBLIC hdf5_get_height
PUBLIC hdf5_get_direct
PUBLIC hdf5_reset
PUBLIC hdf5_copy_node
PUBLIC hdf5_get_dims
PUBLIC hdf5_get_map
PUBLIC hdf5_set_map
PUBLIC hdf5_set_pointer
!
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE hdf5_read(infile, length, O_LAYER, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
                     xmin, xmax, ymin, ymax, offxy, offz, lni, lh5, lenc,       &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!
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
INTEGER, INTENT(IN)    :: MAXARRAY     ! KUPLOT array size
INTEGER, INTENT(IN)    :: MAXKURVTOT   ! KUPLOT array size
CHARACTER(LEN=200), DIMENSION(MAXKURVTOT), INTENT(INOUT) :: fname
INTEGER, INTENT(INOUT) :: iz     ! KUPLOT data set number
REAL   , DIMENSION(MAXARRAY)  , INTENT(INOUT) :: x
REAL   , DIMENSION(MAXARRAY)  , INTENT(INOUT) :: y
REAL   , DIMENSION(MAXARRAY)  , INTENT(INOUT) :: z
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: nx
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ny
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: xmax ! (maxkurvtot)
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: xmin ! (maxkurvtot)
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ymax ! (maxkurvtot)
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ymin
INTEGER, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offxy
INTEGER, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offz
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lni
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lh5
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: lenc
!
INTEGER,                            INTENT(OUT)   :: ier_num
INTEGER,                            INTENT(OUT)   :: ier_typ
INTEGER,                            INTENT(IN )   :: idims
CHARACTER(LEN=*), DIMENSION(idims), INTENT(INOUT) :: ier_msg    ! Error message
INTEGER,                            INTENT(IN )   :: ER_APPL
INTEGER,                            INTENT(IN )   :: ER_IO
INTEGER, INTENT(IN)    :: output_io
!
!
ier_num = -72   ! HDF5 not supported
ier_typ = ER_APPL !
!
END SUBROUTINE hdf5_read
!
!*******************************************************************************
!
SUBROUTINE hdf5_place_kuplot(nlayer, lset, lnew, lshow,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
   xmin, xmax, ymin, ymax, &
   offxy, offz, lni, lh5, lenc, ier_num, ier_typ, output_io)
!-
! PLace a curve into the kuplot section, 
! IF lset==TRUE set absolute layer , else increment
! IF lnew==TRUE, make new curve, 
! IF lshow = TRUE display data
!+
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nlayer    ! Cut this layer from the data
LOGICAL, INTENT(IN) :: lset      ! absolute layer setting
LOGICAL, INTENT(IN) :: lnew      ! make new curve
LOGICAL, INTENT(IN) :: lshow     ! show data
INTEGER, INTENT(IN)    :: MAXARRAY     ! KUPLOT array size
INTEGER, INTENT(IN)    :: MAXKURVTOT   ! KUPLOT array size
CHARACTER(LEN=200), DIMENSION(MAXKURVTOT), INTENT(INOUT) :: fname
INTEGER, INTENT(INOUT) :: iz     ! KUPLOT data set number
REAL   , DIMENSION(MAXARRAY)  , INTENT(INOUT) :: x
REAL   , DIMENSION(MAXARRAY)  , INTENT(INOUT) :: y
REAL   , DIMENSION(MAXARRAY)  , INTENT(INOUT) :: z
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: nx
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ny
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: xmax ! (maxkurvtot)
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: xmin ! (maxkurvtot)
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ymax ! (maxkurvtot)
REAL   , DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ymin
INTEGER, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offxy
INTEGER, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offz
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lni
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lh5
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: lenc
INTEGER,                 INTENT(OUT) :: ier_num
INTEGER,                 INTENT(OUT) :: ier_typ
INTEGER, INTENT(IN)    :: output_io 
!
ier_num = -72   ! HDF5 not supported
ier_typ = 6     ! ER_APPL
!
END SUBROUTINE hdf5_place_kuplot
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_layer()
!
IMPLICIT NONE
!
hdf5_get_layer = 1
!
END FUNCTION hdf5_get_layer
!
!*******************************************************************************
!
LOGICAL FUNCTION hdf5_get_direct()
!
IMPLICIT NONE
!
hdf5_get_direct = .TRUE.
!
END FUNCTION hdf5_get_direct
!
!*******************************************************************************
!
REAL FUNCTION hdf5_get_height()
!
IMPLICIT NONE
!
hdf5_get_height = 0.0
!
END FUNCTION hdf5_get_height
!
!*******************************************************************************
!
SUBROUTINE hdf5_copy_node(old, new)
!-
!  Copies old node to new node
!+
INTEGER, INTENT(IN)  :: old
INTEGER, INTENT(OUT) :: new
!
end subroutine hdf5_copy_node
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_dims(idata, dims)
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN)  :: idata
INTEGER, DIMENSION(3), INTENT(OUT) :: dims
!
dims = 1
!
END SUBROUTINE hdf5_get_dims
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_map(dims, odata)
!
IMPLICIT NONE
!
INTEGER,            DIMENSION(3),                         INTENT(IN)  :: dims
REAL(KIND=PREC_DP), DIMENSION(dims(1), dims(2), dims(3)), INTENT(OUT) :: odata
!
!INTEGER :: i,j,k
!
odata = 0.0D0
!
END SUBROUTINE hdf5_get_map
!
!*******************************************************************************
!
SUBROUTINE hdf5_set_map(dims, odata)
!
IMPLICIT NONE
!
INTEGER,            DIMENSION(3),                         INTENT(IN) :: dims
REAL(KIND=PREC_DP), DIMENSION(dims(1), dims(2), dims(3)), INTENT(IN) :: odata
!
!
END SUBROUTINE hdf5_set_map
!
!*******************************************************************************
!
SUBROUTINE hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
!-
!  Find the node associated to kuplot data set number izz
!+
INTEGER, INTENT(IN ) :: izz
INTEGER, INTENT(OUT) :: ier_num
INTEGER, INTENT(OUT) :: ier_typ
INTEGER, INTENT(OUT) :: node_number
!
end subroutine hdf5_set_pointer
!
!*******************************************************************************
!
SUBROUTINE hdf5_reset
!
implicit none
!
END SUBROUTINE hdf5_reset
!
!*******************************************************************************
!
END MODULE kuplot_load_h5
