MODULE kuplot_place
!
!  Place a multidimensional data set into the KUPLOT data structure
!
contains
!
subroutine place_kuplot(h5_dims, nlayer, lset, lnew, lshow,                     &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
      xmin, xmax, ymin, ymax, &
      offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
!-
! Place into 1D or 2D field
!
!use errlist_mod
use precision_mod
!
implicit none
!
integer, dimension(3), intent(in) :: h5_dims
INTEGER, INTENT(IN) :: nlayer    ! Cut this layer from the data
LOGICAL, INTENT(IN) :: lset      ! absolute layer setting
LOGICAL, INTENT(IN) :: lnew      ! make new curve
LOGICAL, INTENT(IN) :: lshow     ! show data
INTEGER, INTENT(IN)    :: MAXARRAY     ! KUPLOT array size
INTEGER, INTENT(IN)    :: MAXKURVTOT   ! KUPLOT array size
CHARACTER(LEN=200), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: fname
INTEGER                                    , INTENT(INOUT) :: iz     ! KUPLOT data set number
REAL(kind=PREC_DP), DIMENSION(MAXARRAY)    , INTENT(INOUT) :: x
REAL(kind=PREC_DP), DIMENSION(MAXARRAY)    , INTENT(INOUT) :: y
REAL(kind=PREC_DP), DIMENSION(MAXARRAY)    , INTENT(INOUT) :: z
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: nx
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ny
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: xmax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: xmin ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ymax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ymin
INTEGER           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offxy
INTEGER           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offz
LOGICAL           , DIMENSION(  maxkurvtot), INTENT(INOUT) :: lni
LOGICAL           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: lh5
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ku_ndims
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: lenc
INTEGER                                    ,                 INTENT(OUT) :: ier_num
INTEGER                                    ,                 INTENT(OUT) :: ier_typ
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!if(h5_temp%h5_dims(1)==1 .and. h5_temp%h5_dims(2)==1) then
write(*,*) ' H5_DIMS place-kuplot ', h5_dims(:)
if(h5_dims(1)==1 .and. h5_dims(2)==1) then
   CALL hdf5_place_kuplot_1d(nlayer, lset  , lnew, lshow  ,               &
      MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
      xmin, xmax, ymin, ymax, &
      offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
else
write(*,*) ' CALLING Ndimensional '
   CALL hdf5_place_kuplot(nlayer,  lset, lnew, lshow,                     &
      MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
      xmin, xmax, ymin, ymax, &
      offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
endif
write(*,*) ' LOADED    ', ku_ndims(:4)
write(*,*) ' LH5 lbound', lbound(lh5), ubound(lh5)
write(*,*) ' ku_ndims  ', lbound(ku_ndims), ubound(ku_ndims)
!
!
END SUBROUTINE place_kuplot
!
!*******************************************************************************
!
SUBROUTINE hdf5_place_kuplot(nlayer, lset, lnew, lshow,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
   xmin, xmax, ymin, ymax, &
   offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
!
!-
! PLace a curve into the kuplot section, 
! IF lset==TRUE set absolute layer , else increment
! IF lnew==TRUE, make new curve, 
! IF lshow = TRUE display data
!+
!
use kuplot_show_mod
!
use lib_data_struc_h5
use precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nlayer    ! Cut this layer from the data
LOGICAL, INTENT(IN) :: lset      ! absolute layer setting
LOGICAL, INTENT(IN) :: lnew      ! make new curve
LOGICAL, INTENT(IN) :: lshow     ! show data
INTEGER, INTENT(IN)    :: MAXARRAY     ! KUPLOT array size
INTEGER, INTENT(IN)    :: MAXKURVTOT   ! KUPLOT array size
CHARACTER(LEN=200), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: fname
INTEGER                                    , INTENT(INOUT) :: iz     ! KUPLOT data set number
REAL(kind=PREC_DP), DIMENSION(  MAXARRAY)  , INTENT(INOUT) :: x
REAL(kind=PREC_DP), DIMENSION(  MAXARRAY)  , INTENT(INOUT) :: y
REAL(kind=PREC_DP), DIMENSION(  MAXARRAY)  , INTENT(INOUT) :: z
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: nx
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ny
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: xmax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: xmin ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ymax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ymin
INTEGER           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offxy
INTEGER           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offz
LOGICAL           , DIMENSION(  maxkurvtot), INTENT(INOUT) :: lni
LOGICAL           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: lh5
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ku_ndims
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: lenc
INTEGER                                    , INTENT(IN)    :: output_io   ! KUPLOT array size
!
INTEGER,                 INTENT(OUT) :: ier_num
INTEGER,                 INTENT(OUT) :: ier_typ
!
INTEGER :: i,j,k, ll             ! dummy indices
INTEGER :: izz
INTEGER :: node_number
character(len=PREC_STRING)       :: h5_infile
integer, dimension(3)            :: h5_dims
integer                          :: h5_layer
integer                          :: h5_number
real(kind=PREC_DP), dimension(3) :: h5_llims
real(kind=PREC_DP), dimension(3,3) :: h5_steps
!
IF(lnew) THEN            ! This is a new data set, from 'load' command
   h5_number = hdf5_get_number()
   izz = iz
!  h5_h5_is_ku(h5_number) = izz
!  h5_ku_is_h5(izz      ) = h5_number
   call hdf5_set_h5_is_ku(h5_number, izz)
   call hdf5_set_ku_is_h5(izz, h5_number)
   node_number = h5_number
ELSE                     ! Overwrite current KUPLOT data set
   izz = iz - 1
ENDIF
!                        ! Locate this data set in the h5 storage
CALL hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
if(ier_num /= 0) return
!
call hdf5_get_dims(node_number, h5_dims)
call hdf5_get_llims(node_number, h5_llims)
call hdf5_get_steps(node_number, h5_steps)
call hdf5_get_infile(node_number, h5_infile)
!
write(*,*) ' h5_DIMS ', h5_dims(:), ' IZZ ', izz
if(h5_dims(1)>1 .and. h5_dims(2)>1 .and. h5_dims(3)>1) then
  ku_ndims(izz) = 3
elseif(h5_dims(1)==1 .and. h5_dims(2)>1 .and. h5_dims(3)>1) then
  ku_ndims(izz) = 2
endif
!
IF(lset) THEN
  call hdf5_set_layer(nlayer)
  h5_layer = nlayer
ELSE
  h5_layer =  hdf5_get_layer()
  h5_layer = MAX(1,MIN(INT(           h5_dims(1)),         h5_layer+nlayer))
ENDIF
ll = 0
k = h5_layer
!
DO i = 1, h5_dims(3)
  DO j = 1, h5_dims(2)
     ll = ll + 1
     z(offz(izz - 1) + ll ) = hdf5_get_data(i,j,k)
  ENDDO
ENDDO
!
nx(izz) =            h5_dims(3)
ny(izz) =            h5_dims(2)
xmin(izz) = h5_llims(1)
xmax(izz) = h5_llims(1) + (nx(izz)-1)*h5_steps(1,1)
ymin(izz) = h5_llims(2)
ymax(izz) = h5_llims(2) + (ny(izz)-1)*h5_steps(2,2)
!
DO i = 1, nx(izz)
    x(offxy(izz - 1) + i) = xmin(izz) + (i - 1) * h5_steps(1,1)
ENDDO
DO i = 1, ny(izz)
   y(offxy(izz - 1) + i) = ymin(izz) + (i - 1) * h5_steps(2,2)
ENDDO
lni (izz) = .TRUE.
lh5 (izz) = .TRUE.
lenc(izz) = MAX(nx(izz), ny(izz))
offxy(izz) = offxy(izz - 1) + lenc(izz)
offz (izz) = offz (izz - 1) + nx(izz) * ny(izz)
fname(izz) = h5_infile(1:LEN_TRIM(h5_infile))
call hdf5_set_h5_is_ku(node_number, izz) ! H5 Data set 1 is stored in Kuplot as number izz
call hdf5_set_h5_is_ku(izz, node_number) ! Kuplot data set izz is stored in H5 number 1
call hdf5_set_layer(h5_layer)
IF(lnew) iz = iz + 1
!
IF(lshow) THEN
   CALL show_data(iz - 1)!
   WRITE(output_io,1000) h5_dims(3), h5_dims(2), h5_dims(1)
   WRITE(output_io,1100) nlayer
   1000 FORMAT('   Full size:', 2(i7,' x'), i7, ' points')
   1100 FORMAT('   At  layer:',   i7      ,/)
ENDIF
write(*,*) ' h5_ku_N   ', ku_ndims(:4)
!
END SUBROUTINE hdf5_place_kuplot
!
!*******************************************************************************
!
SUBROUTINE hdf5_place_kuplot_1d(nlayer, lset, lnew, lshow,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
   xmin, xmax, ymin, ymax, &
   offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
!
!-
! PLace a 1D curve into the kuplot section, 
! IF lset==TRUE set absolute layer , else increment
! IF lnew==TRUE, make new curve, 
! IF lshow = TRUE display data
!+
!
use kuplot_show_mod
use lib_data_struc_h5
use precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nlayer    ! Cut this layer from the data
LOGICAL, INTENT(IN) :: lset      ! absolute layer setting
LOGICAL, INTENT(IN) :: lnew      ! make new curve
LOGICAL, INTENT(IN) :: lshow     ! show data
INTEGER, INTENT(IN)    :: MAXARRAY     ! KUPLOT array size
INTEGER, INTENT(IN)    :: MAXKURVTOT   ! KUPLOT array size
CHARACTER(LEN=200), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: fname
INTEGER                                    , INTENT(INOUT) :: iz     ! KUPLOT data set number
REAL(kind=PREC_DP), DIMENSION(  MAXARRAY)  , INTENT(INOUT) :: x
REAL(kind=PREC_DP), DIMENSION(  MAXARRAY)  , INTENT(INOUT) :: y
REAL(kind=PREC_DP), DIMENSION(  MAXARRAY)  , INTENT(INOUT) :: z
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: nx
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ny
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: xmax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: xmin ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ymax ! (maxkurvtot)
REAL(kind=PREC_DP), DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ymin
INTEGER           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offxy
INTEGER           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: offz
LOGICAL           , DIMENSION(  maxkurvtot), INTENT(INOUT) :: lni
LOGICAL           , DIMENSION(0:maxkurvtot), INTENT(INOUT) :: lh5
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: ku_ndims
INTEGER           , DIMENSION(  MAXKURVTOT), INTENT(INOUT) :: lenc
INTEGER                                    , INTENT(OUT)   :: ier_num
INTEGER                                    , INTENT(OUT)   :: ier_typ
INTEGER                                    , INTENT(IN)    :: output_io   ! KUPLOT array size
!
!
INTEGER :: i                     ! dummy indices
INTEGER :: izz
INTEGER :: node_number
character(len=PREC_STRING)       :: h5_infile
integer, dimension(3)            :: h5_dims
integer                          :: h5_layer
integer                          :: h5_number
real(kind=PREC_DP), dimension(3) :: h5_llims
real(kind=PREC_DP), dimension(3,3) :: h5_steps
!
!
IF(lnew) THEN            ! This is a new data set, from 'load' command
   h5_number = hdf5_get_number()
   izz = iz
   call hdf5_set_h5_is_ku(h5_number, izz)
   call hdf5_set_ku_is_h5(izz, h5_number)
   node_number = h5_number
ELSE                     ! Overwrite current KUPLOT data set
   izz = iz - 1
ENDIF
!                        ! Locate this data set in the h5 storage
CALL hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
if(ier_num /= 0) return
!
call hdf5_get_dims(node_number, h5_dims)
call hdf5_get_llims(node_number, h5_llims)
call hdf5_get_steps(node_number, h5_steps)
call hdf5_get_infile(node_number, h5_infile)
!
IF(lset) THEN
  call hdf5_set_layer(nlayer)
ELSE
  h5_layer =  hdf5_get_layer()
  h5_layer = MAX(1,MIN(INT(           h5_dims(1)),         h5_layer+nlayer))
ENDIF
!
lenc(izz) = h5_dims(3)
xmin(izz) = h5_llims(1)
xmax(izz) = h5_llims(1) + (lenc(izz)-1)*h5_steps(1,1)
DO i = 1, h5_dims(3)
  x(offxy(izz - 1) + i) = xmin(izz) + (i - 1) * h5_steps(1,1)
  y(offxy(izz - 1) + i) = hdf5_get_data  (i,1,1)
ENDDO
!
ymin(izz) = minval(y(offxy(izz - 1) + 1: offxy(izz - 1) + lenc(izz)))
ymax(izz) = maxval(y(offxy(izz - 1) + 1: offxy(izz - 1) + lenc(izz)))
!
ku_ndims(izz) = 1
lni (izz) = .false.
lh5 (izz) = .true. 
offxy(izz) = offxy(izz - 1) + lenc(izz)
fname(izz) = h5_infile(1:LEN_TRIM(h5_infile))
call hdf5_set_h5_is_ku(node_number, izz) ! H5 Data set 1 is stored in Kuplot as number izz
call hdf5_set_h5_is_ku(izz, node_number) ! Kuplot data set izz is stored in H5 number 1
IF(lnew) iz = iz + 1
!
IF(lshow) THEN
   CALL show_data(iz - 1)!
ENDIF
!
END SUBROUTINE hdf5_place_kuplot_1d
!
!*******************************************************************************
!
END MODULE kuplot_place
