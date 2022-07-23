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
!PUBLIC hdf5_place_kuplot
!
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE hdf5_read_kuplot(infile, length, O_LAYER, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
                     xmin, xmax, ymin, ymax, offxy, offz, lni, lh5, ku_ndims,lenc,       &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!
use kuplot_place
!
USE hdf5
USE iso_c_binding
!
USE ber_params_mod
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
REAL(kind=PREC_DP), DIMENSION(MAXARRAY)  , INTENT(INOUT) :: x
REAL(kind=PREC_DP), DIMENSION(MAXARRAY)  , INTENT(INOUT) :: y
REAL(kind=PREC_DP), DIMENSION(MAXARRAY)  , INTENT(INOUT) :: z
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: nx
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ny
REAL(kind=PREC_DP), DIMENSION(MAXKURVTOT), INTENT(INOUT) :: xmax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(MAXKURVTOT), INTENT(INOUT) :: xmin ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ymax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(MAXKURVTOT), INTENT(INOUT) :: ymin
INTEGER, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offxy
INTEGER, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offz
LOGICAL, DIMENSION(  maxkurvtot), INTENT(INOUT) :: lni
LOGICAL, DIMENSION(0:maxkurvtot), INTENT(INOUT) :: lh5
INTEGER, DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ku_ndims
INTEGER, DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: lenc
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
CHARACTER(LEN=14)   :: dataname    ! Dummy name for HDF5 datasets
!
integer               :: idata = 0
integer, dimension(3) :: h5_dims
!
INTEGER, PARAMETER                          :: MAXW = 1
INTEGER                                     :: ianz
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER                   , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
INTEGER :: nlayer                   ! Layer to place into KUPLOT
!
dataname = ' '
!

call hdf5_read(infile, length, O_LAYER, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
                     xmin, xmax, ymin, ymax, offxy, offz, lni, lh5, lenc,       &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!
call hdf5_get_dims(idata, h5_dims)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy into KUPLOT array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IF(opara(O_LAYER)=='bottom') THEN
   nlayer = 1
ELSEIF(opara(O_LAYER)=='middle') THEN
   nlayer = INT((h5_dims(1)+1)/2)
ELSEIF(opara(O_LAYER)=='top') THEN
   nlayer = h5_dims(1)
ELSE
   cpara(1) = opara(O_LAYER)
   lpara(1) = lopara(O_LAYER)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, maxw)
   nlayer = NINT(werte(1))
   IF(nlayer  <=0) THEN
      ier_num = -71
      ier_typ = ER_APPL
      ier_msg(1) = 'Layer number <= 0'
      ier_msg(2) = 'FILE '//infile (1:LEN(ier_msg)-5)
      RETURN
   ELSEIF(nlayer  >  h5_dims(1)) THEN
      ier_num = -71
      ier_typ = ER_APPL
      WRITE(ier_msg(1),'(a,i4)') 'Layer number > ', nlayer
      ier_msg(2) = 'FILE '//infile (1:LEN(ier_msg)-5)
      RETURN
   ENDIF
ENDIF
!
CALL place_kuplot(h5_dims, nlayer, .TRUE.,.TRUE., .TRUE.,               &
     MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
     xmin, xmax, ymin, ymax, &
     offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
write(*,*) ' PLACED INTO KUPLOT ', lh5(:2), ' DIMS ', ku_ndims(:4)
!if(h5_temp%h5_dims(1)==1 .and. h5_temp%h5_dims(2)==1) then
!if(h5_dims(1)==1 .and. h5_dims(2)==1) then
!   CALL hdf5_place_kuplot_1d(nlayer, .TRUE.,.TRUE., .TRUE.,               &
!      MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
!      xmin, xmax, ymin, ymax, &
!      offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
!else
!   CALL hdf5_place_kuplot(nlayer, .TRUE.,.TRUE., .TRUE.,               &
!      MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
!      xmin, xmax, ymin, ymax, &
!      offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
!endif
!
!
END SUBROUTINE hdf5_read_kuplot
!
!*******************************************************************************
!
!
!*******************************************************************************
!
END MODULE kuplot_load_h5
