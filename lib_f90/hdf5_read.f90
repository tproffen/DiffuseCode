module lib_hdf5_read_mod
!
! Load an HDF5 file into the lib_f90 data structure
!
use errlist_mod
use lib_data_struc_h5
use lib_hdf5_params_mod
use hdf5_def_mod
use precision_mod
!
private
public hdf5_read     ! Read the pure Yell/discus file format
public nx_read_scattering_place   ! Read the DiffuseDevelopers format via FORPY
                                  ! and place into global data structure
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: h5_datasets       ! Names of the data set in file
INTEGER                                               :: ndims             ! Number of dimensions
INTEGER                                               :: one_ndims         ! Number of dimensions
INTEGER                                               :: H5_MAX_DATASETS   ! Current MAX data sets
INTEGER                                               :: h5_data_type      ! Data type 
INTEGER                                               :: h5_n_datasets     ! Current actual data sets
INTEGER                                               :: h5_layer=1        ! Current layer in data set
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: h5_dims           ! Actual dimensions
INTEGER                  , DIMENSION(3)               :: d5_dims           ! Actual dimensionsA in transposed sequence
!NTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: maxdims           ! Maximum dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_dims          ! Actual dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_maxdims       ! Maximum dimensions
LOGICAL                                               :: h5_direct         ! Direct space == TRUE
LOGICAL                                               :: h5_is_grid=.true. ! Data on periodic grid
LOGICAL                                               :: h5_has_dxyz=.false. ! Data on periodic grid
LOGICAL                                               :: h5_has_dval=.false. ! Data on periodic grid
REAL(KIND=PREC_DP)   , DIMENSION(3,4)                 :: h5_corners        ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_vectors        ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(6)                   :: h5_unit           ! Lattice parameters
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_x              ! Actual x-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_y              ! Actual y-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_z              ! Actual z-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_dx             ! Actual x-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_dy             ! Actual y-coordinates
REAL(KIND=PREC_DP)   , DIMENSION(:)    , ALLOCATABLE  :: h5_dz             ! Actual z-coordinates
REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_data           ! Actual diffraction data
REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_sigma          ! Actual diffraction data
REAL(KIND=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  :: d5_data           ! Actual diffraction data  in transposed sequence 
REAL(KIND=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  :: d5_sigma          ! Actual diffraction sigma in transposed sequence 
REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_llims          ! Lower limits
REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_steps          ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_steps_full     ! steps in H, K, L
logical                                               :: h5_calc_coor      ! TRUE: Need to calculate coordinates
integer              , dimension(3)                   :: h5_use_coor       ! Use this index for x, y, z, coordinates
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE hdf5_read(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number, nndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!-
!  Read the pure Yell / DISCUS file format from HDF5 in fortran directly
!+
!
use hdf5
use iso_c_binding
!
use ber_params_mod
use lib_trans_mod
use lib_data_types_mod
use lib_use_coor_mod
use precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1024), INTENT(INOUT) :: infile
INTEGER            , INTENT(IN) :: length
INTEGER            , INTENT(IN) :: O_LAYER    ! Number optional parameter "layer:"
INTEGER            , INTENT(IN) :: O_TRANS    ! Number optional parameter "trans:"
INTEGER            , INTENT(IN) :: NOPTIONAL
CHARACTER(LEN=*)   , DIMENSION(NOPTIONAL), INTENT(IN) :: opara
INTEGER            , DIMENSION(NOPTIONAL), INTENT(IN) :: lopara
LOGICAL            , DIMENSION(NOPTIONAL), INTENT(IN) :: lpresent
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL), INTENT(IN) :: owerte
!
integer, intent(out) :: node_number
integer, intent(out) :: nndims
integer, dimension(3), intent(out) :: dims
INTEGER,                            INTENT(OUT)   :: ier_num
INTEGER,                            INTENT(OUT)   :: ier_typ
INTEGER,                            INTENT(IN )   :: idims
CHARACTER(LEN=*), DIMENSION(idims), INTENT(INOUT) :: ier_msg    ! Error message
INTEGER,                            INTENT(IN )   :: ER_APPL
INTEGER,                            INTENT(IN )   :: ER_IO
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!INTEGER(KIND=SIZE_T), PARAMETER :: f_str_len = 20  ! Length F-string to receive  
!
CHARACTER(LEN=14)   :: dataname    ! Dummy name for HDF5 datasets
!CHARACTER(LEN=f_str_len)   :: progname    ! Program that created H5 FILE 
!CHARACTER(LEN=f_str_len), DIMENSION(:), ALLOCATABLE, TARGET   :: rstring    ! string that is read
!CHARACTER(LEN=9), DIMENSION(1:1), TARGET ::  rstring
INTEGER(KIND=HID_T) :: file_id     ! File identifier
!INTEGER(KIND=HID_T) :: filetype    ! Filetype identifier
!INTEGER(KIND=HID_T) ::  memtype    ! Mem type identifier
INTEGER(KIND=HID_T) :: dset_id     ! dataset identifier
INTEGER(KIND=HID_T) :: space_id    ! space identifier
!INTEGER(KIND=SIZE_T) :: c_str_len  ! Length C-string in 'PROGRAM'
!INTEGER(KIND=HSIZE_T), DIMENSION(2) :: cdims
!INTEGER(KIND= SIZE_T), DIMENSION(2) :: clens=(/8,1/)
INTEGER             :: hdferr      ! Error number
INTEGER(KIND=HSIZE_T) :: idx
INTEGER             :: ret_value   ! Error number
TYPE(C_FUNPTR)      :: funptr      ! Pointer to display function
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata   ! READ buffer
TYPE(C_PTR)         :: f_ptr        ! Pointer to data elements
TYPE(C_PTR) :: ptr
CHARACTER(LEN=1024, KIND=c_char),  POINTER ::  rstring
INTEGER(KIND=2), PARAMETER :: SHORT_TWO = 2
INTEGER(KIND=2), TARGET :: r_is_direct
!INTEGER(KIND=2), DIMENSION(10), TARGET :: rstring
!INTEGER(HSIZE_T) :: nummer
!
!INTEGER, PARAMETER                          :: MAXW = 1
!INTEGER                                     :: ianz
!CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
!INTEGER                   , DIMENSION(MAXW) :: lpara
!REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP)        , DIMENSION(3)    :: steps     ! dummy steps in H, K, L
INTEGER :: i,j,k                    ! Dummy loop indices
!
!integer, parameter :: VAL_PDF   = 14
!integer, parameter :: VAL_3DPDF = 15
!integer :: value  ! Intensity = 1; VAL_PDF or VAL_3DPDF
!
!                                                                          ! If data are transformed upon input use "new"
!character(len=PREC_STRING)                                :: new_outfile   ! New file name  from transformed HDF5 data
!integer                   , dimension(3)                  :: new_inc       ! New dimensions from transformed HDF5 data
!real(kind=PREC_DP)        , dimension(3,4)                :: new_eck       ! New corners    from transformed HDF5 data
!real(kind=PREC_DP)        , dimension(3,3)                :: new_vi        ! New vectors    from transformed HDF5 data
!real(kind=PREC_DP)        , dimension(:,:,:), allocatable :: new_data      ! New data       from transformed HDF5 data
!
real(kind=PREC_DP), dimension(3,3) :: temp_vi    ! Temporary copy
integer :: extr_abs, extr_ord, extr_top
!
! Make sure all old data are gone
!
if(allocated(h5_datasets)) deallocate(h5_datasets)
if(allocated(h5_data)) deallocate(h5_data)
if(allocated(d5_data)) deallocate(d5_data)
if(allocated(h5_sigma)) deallocate(h5_sigma)
if(allocated(d5_sigma)) deallocate(d5_sigma)
if(allocated(h5_x)) deallocate(h5_x)
if(allocated(h5_y)) deallocate(h5_y)
if(allocated(h5_z)) deallocate(h5_z)
if(allocated(h5_dx)) deallocate(h5_dx)
if(allocated(h5_dy)) deallocate(h5_dy)
if(allocated(h5_dz)) deallocate(h5_dz)
!
h5_infile = infile
dataname = ' '
h5_steps_full = 0.0D0 
!
h5_steps      = 0.0D0
!
!
H5_MAX_DATASETS = 10                                        ! Initial estimate of dataset number
ALLOCATE(h5_datasets(H5_MAX_DATASETS))
h5_n_datasets = 0                                           ! Currently no datasets found
!
h5_dims    = 1
one_maxdims = 1
CALL H5open_f(hdferr)                                       ! Open access to HDF5 stream
CALL H5Eset_auto_f(1, hdferr)                              ! Turn Error messages off
!
CALL H5Fopen_f(h5_infile, H5F_ACC_RDWR_F, file_id, hdferr)     ! Open existing file
IF(hdferr/=0) THEN
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -2, ER_IO, 'Could not open ''H5'' file')
   RETURN
ENDIF
!
idx = 0
funptr = C_FUNLOC(op_func) ! call back function
ptr    = C_NULL_PTR
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iterate across file to obtain format, and list of datasets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CALL H5Literate_f(file_id, H5_INDEX_NAME_F, H5_ITER_NATIVE_F, idx, funptr, ptr, ret_value, hdferr)
IF(hdferr/=0) THEN
   dataname = ' '
   dset_id  = 0
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -70, ER_APPL, 'Initial iteration failed')
   RETURN
ENDIF
!
IF(h5_n_datasets<1) THEN
   dataname = ' '
   dset_id  = 0
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -69, ER_APPL, 'No datasets in H5 file')
   RETURN
ENDIF
!write(*,*) 'NUMBER of data sets found ', h5_n_datasets
yd_present = .false.              ! Assume all data set missing
DO i=1,h5_n_datasets
   loop_ydnd:do j=1, YD_ND
      if(h5_datasets(i)==yd_datasets(j)) then
         yd_present(j) = .true.
         exit loop_ydnd
      endif
   enddo loop_ydnd
!   WRITE(*,*) ' DATASETS ', h5_datasets(i)(1:LEN_TRIM(h5_datasets(i)))
ENDDO
!do j=1, YD_ND
!   write(*,*) ' DATASETS ', yd_present(j),yd_datasets(j)
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the format identifier
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dataname = 'format'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the 'format'
IF(hdferr/=0) THEN
   CALL hdf5_error(h5_infile, file_id, dataname, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -69, ER_APPL, 'Dataset does not exist')
   RETURN
ENDIF
!CALL H5Dget_storage_size_f(dset_id, nummer, hdferr)          ! Get rough storage size
!
ALLOCATE(rdata(1:1024))
f_ptr = C_LOC(rdata(1))
CALL H5Dread_F(dset_id, H5T_STRING, f_ptr, hdferr)
CALL C_F_POINTER(rdata(1), rstring)
j = 0
DO
   IF(j>8) EXIT
   IF(rstring(j+1:j+1)==C_NULL_CHAR) EXIT
   j = j+1
ENDDO
!write(*,*) ' RSTRING ', rstring(1:j), ' J ', j, LEN_TRIM(rstring)
CALL H5Dclose_f(dset_id , hdferr)
!
DEALLOCATE(rdata)
!
IF(rstring(1:j)/='Yell 1.0') THEN
   CALL hdf5_error(h5_infile, file_id, rstring, dset_id, ier_num, ier_typ, idims,ier_msg,&
                   -69, ER_APPL, 'format dataset has wrong value')
   RETURN
ENDIF
!
!Qdataname  = 'PROGRAM'
!QCALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the 'PROGRAM'
!Qwrite(*,*) ' PROGRAM open  ', hdferr
!Q!
!Q! Get the datatype and its size.
!Q!
!QCALL H5Dget_type_f(dset_id, filetype, hdferr)
!Qwrite(*,*) ' GET_TYPE  ', hdferr
!QCALL H5Tget_size_f(filetype, c_str_len, hdferr)
!Qwrite(*,*) ' C_STR_LEN ', c_str_len, hdferr
!QIF(c_str_len> f_str_len+1) THEN                              ! Fortran string is too short
!Q   CALL H5Dclose_f(dset_id , hdferr)
!Q   CALL H5close_f(hdferr)                                    ! Close HDF interface
!Q   RETURN
!QENDIF
!Q!
!Q! Get dataspace.
!Q!
!QCALL H5Dget_space_f(dset_id, space_id, hdferr)
!Qwrite(*,*) ' GET_SPACE  ', hdferr
!QCALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)
!Qwrite(*,*) ' GET_EXTENT ', one_dims, hdferr
!QALLOCATE(rstring(1:one_dims(1)))
!Q!
!Q! Create the memory datatype.
!Q!
!QCALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
!QCALL H5Tset_size_f(memtype, f_str_len, hdferr)
!Q!
!Q! Read the data.
!Q!
!Qf_ptr = C_LOC(rstring(1)(1:1))
!QCALL H5Dread_f(dset_id, memtype, f_ptr, hdferr, space_id)
!Qprogname = ' '
!Qprogname = rstring(1)
!Qwrite(*,*) ' PROGNAME ', progname
!Qwrite(*,*) ' RSTRING  ', rstring
!Q!
!QCALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!QCALL H5Sclose_f(space_id, hdferr)
!QCALL H5Tclose_f(filetype, hdferr)
!QCALL H5Tclose_f(memtype, hdferr)
!Q!CALL H5Fclose_f(file_id, hdferr)
!Qwrite(*,*) ' PROGRAM close ', hdferr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
one_maxdims = 1
dataname='data'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, h5_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
IF(ALLOCATED(h5_data)) DEALLOCATE(h5_data)
ALLOCATE(h5_data(h5_dims(1), h5_dims(2), h5_dims(3)))
CALL H5Dread_f(dset_id, H5T_NATIVE_REAL, h5_data, h5_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'lower_limits'
one_dims    = 1
one_maxdims = 1
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, h5_llims, one_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the is_direct dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'is_direct'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
f_ptr = C_LOC(r_is_direct)
CALL H5Dread_F(dset_id, H5T_STD_I8BE , f_ptr, hdferr)
CALL H5Dclose_f(dset_id, hdferr)
!h5_direct = 1 == (r_is_direct-8192)
h5_direct = (abs(mod(r_is_direct,SHORT_TWO))==1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the unit cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'unit_cell'
one_dims    = 1
one_maxdims = 1
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, h5_unit , one_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the stepss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dataname = 'step_sizes'
one_dims    = 1
one_maxdims = 1
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, h5_steps, one_dims, hdferr)
CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the steps DISCUS style
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if(yd_present(YD_step_sizes_abs)) then
   dataname = 'step_sizes_abs'
   one_dims    = 1
   one_maxdims = 1
   CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
   if(hdferr==0) then
      CALL H5Dget_space_f(dset_id, space_id, hdferr)
      CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
      CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
      CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, steps, one_dims, hdferr)
      h5_steps_full(:,1) = steps
      CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
   endif
else
   h5_steps_full(1,1) = h5_steps(1)
endif
!
if(yd_present(YD_STEP_SIZES_ORD)) then
   dataname = 'step_sizes_ord'
   one_dims    = 1
   one_maxdims = 1
   CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
   if(hdferr==0) then
      CALL H5Dget_space_f(dset_id, space_id, hdferr)
      CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
      CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
      CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, steps, one_dims, hdferr)
      h5_steps_full(:,2) = steps
      CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
   endif
else
   h5_steps_full(2,2) = h5_steps(2)
endif
!
if(yd_present(YD_STEP_SIZES_TOP)) then
   dataname = 'step_sizes_top'
   one_dims    = 1
   one_maxdims = 1
   CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
   if(hdferr==0) then
      CALL H5Dget_space_f(dset_id, space_id, hdferr)
      CALL H5Sget_simple_extent_ndims_f(space_id, one_ndims, hdferr)   ! Get the number of dimensions in data set
      CALL H5Sget_simple_extent_dims_f(space_id, one_dims, one_maxdims, hdferr)   ! Get the dimensions in data set
      CALL H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, steps, one_dims, hdferr)
      h5_steps_full(:,3) = steps
      CALL H5Dclose_f(dset_id, hdferr)                             ! Close the dataset file
   endif
else
   h5_steps_full(3,3) = h5_steps(3)
endif
!
CALL h5fclose_f(file_id, hdferr)                             ! Close the input file
!
CALL H5close_f(hdferr)                                    ! Close HDF interface
!
!  If full data are present, anylyse to get good estimates for the short 
!  version of the steps
if(yd_present(YD_STEP_SIZES_ABS) .and. yd_present(YD_STEP_SIZES_ORD) .and. &
   yd_present(YD_STEP_SIZES_TOP) ) then
      temp_vi = h5_steps_full
   extr_abs = maxloc(abs(temp_vi(:,1)), dim=1)
   temp_vi(extr_abs,2) = 0.0D0
   extr_ord = maxloc(abs(temp_vi(:,2)), dim=1)
   extr_top = 3
   if(    extr_abs==1 .and. extr_ord==2) then
      extr_top = 3
   elseif(extr_abs==1 .and. extr_ord==3) then
      extr_top = 2
   elseif(extr_abs==2 .and. extr_ord==3) then
      extr_top = 1
   endif
   h5_steps(1) = h5_steps_full(extr_abs,1)
   h5_steps(2) = h5_steps_full(extr_ord,2)
   h5_steps(3) = h5_steps_full(extr_top,3)
endif

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy into H5 storage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
d5_dims(1) = int(h5_dims(3))
d5_dims(2) = int(h5_dims(2))
d5_dims(3) = int(h5_dims(1))
h5_is_grid  = .true.
h5_has_dxyz = .false.
h5_has_dval = .false.
!
allocate(d5_data(d5_dims(1), d5_dims(2), d5_dims(3)))
do i=1, d5_dims(1)
   do j=1, d5_dims(2)
      do k=1, d5_dims(3)
         d5_data(i,j,k) = real(h5_data(k,j,i),kind=PREC_DP)
!        d5_data(i,j,k) = real(h5_data(i,j,k),kind=PREC_DP)
      enddo
   enddo
enddo
if(allocated(h5_sigma)) then
   allocate(d5_sigma(d5_dims(1), d5_dims(2), d5_dims(3)))
   do i=1, d5_dims(1)
      do j=1, d5_dims(2)
         do k=1, d5_dims(3)
            d5_sigma(i,j,k) = h5_sigma(k,j,i)
         enddo
      enddo
   enddo
endif
h5_vectors      = h5_steps_full
h5_corners(:,1) = h5_llims                                          ! Lower left
h5_corners(:,2) = h5_corners(:,1) + (d5_dims(1)-1)* h5_vectors(:,1)   ! Lower right
h5_corners(:,3) = h5_corners(:,1) + (d5_dims(2)-1)* h5_vectors(:,2)   ! Upper left
h5_corners(:,4) = h5_corners(:,1) + (d5_dims(3)-1)* h5_vectors(:,3)   ! Top left
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call hdf5_trans_store &
                    (infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number, nndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!write(*,*) ' DID NEW ROUTINE ', ier_num, ier_typ
if(ier_num==0) return
!QQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!QQ! If requested transform into new orientation
!QQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!QQif(opara(O_TRANS)=='yes') then    ! Transform into different orientation
!QQ   value = 1                      ! Assume Intensity in reciprocal space
!QQ   if(h5_direct) value = VAL_PDF  ! Direct space 
!QQ   call lib_trans_menu(0, value, .false., h5_infile, d5_dims, h5_corners, h5_vectors, h5_unit(1:3),    &
!QQ           h5_unit(4:6), d5_data, VAL_PDF, VAL_3DPDF , &
!QQ        new_outfile, new_inc, new_eck, new_vi, new_data)
!QQ   deallocate(d5_data)
!QQ   if(ier_num/=0) then
!QQ!
!QQ      if(allocated(h5_datasets)) deallocate(h5_datasets)
!QQ      if(allocated(h5_data)) deallocate(h5_data)
!QQ      if(allocated(d5_data)) deallocate(d5_data)
!QQ      if(allocated(h5_sigma)) deallocate(h5_sigma)
!QQ      if(allocated(d5_sigma)) deallocate(d5_sigma)
!QQ      if(allocated(h5_x)) deallocate(h5_x)
!QQ      if(allocated(h5_y)) deallocate(h5_y)
!QQ      if(allocated(h5_z)) deallocate(h5_z)
!QQ      if(allocated(h5_dx)) deallocate(h5_dx)
!QQ      if(allocated(h5_dy)) deallocate(h5_dy)
!QQ      if(allocated(h5_dz)) deallocate(h5_dz)
!QQ      if(allocated(new_data)) deallocate(new_data)
!QQ      return
!QQ   endif
!QQ   h5_infile     = new_outfile
!QQ      infile     = new_outfile
!QQ   h5_llims(1)   = new_eck(1,1)
!QQ   h5_llims(2)   = new_eck(2,1)
!QQ   h5_llims(3)   = new_eck(3,1)
!QQ   h5_dims       = new_inc
!QQ   d5_dims       = new_inc
!QQ   h5_corners    = new_eck
!QQ   h5_vectors    = new_vi
!QQ   h5_steps_full = new_vi
!QQ   h5_steps(1)   = new_vi(1,1)
!QQ   h5_steps(2)   = new_vi(2,2)
!QQ   h5_steps(3)   = new_vi(3,3)
!QQ   allocate(d5_data(d5_dims(1), d5_dims(2), d5_dims(3)))
!QQ   d5_data       = new_data
!QQ   deallocate(new_data)
!QQ   if(allocated(d5_sigma)) deallocate(d5_sigma)
!QQ   allocate(d5_sigma(d5_dims(1), d5_dims(2), d5_dims(3)))
!QQ   d5_sigma = 0.0    ! Currently sigmas are not treated properly
!QQendif
!QQ!
!QQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!QQ! Copy into KUPLOT array
!QQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!QQ!
!QQ!write(*,*) ' D5_DIMS ', d5_dims, opara(O_LAYER)(1:6), INT((d5_dims(3)+1)/2), INT((d5_dims(3)  )/2)
!QQIF(opara(O_LAYER)=='bottom') THEN
!QQ   h5_layer = 1
!QQELSEIF(opara(O_LAYER)=='middle') THEN
!QQ   if(mod(d5_dims(3),2)==0) then
!QQ      h5_layer = INT((d5_dims(3)+1)/2)
!QQ   else
!QQ      h5_layer = INT((d5_dims(3)+1)/2)
!QQ   endif
!QQELSEIF(opara(O_LAYER)=='top') THEN
!QQ   h5_layer = d5_dims(3)
!QQELSE
!QQ   cpara(1) = opara(O_LAYER)
!QQ   lpara(1) = lopara(O_LAYER)
!QQ   ianz = 1
!QQ   CALL ber_params (ianz, cpara, lpara, werte, maxw)
!QQ   h5_layer = NINT(werte(1))
!QQ   IF(h5_layer  <=0) THEN
!QQ      ier_num = -71
!QQ      ier_typ = ER_APPL
!QQ      ier_msg(1) = 'Layer number <= 0'
!QQ      ier_msg(2) = 'FILE '//h5_infile (1:LEN(ier_msg)-5)
!QQ   ELSEIF(h5_layer  >  d5_dims(3)) THEN
!QQ      ier_num = -71
!QQ      ier_typ = ER_APPL
!QQ      WRITE(ier_msg(1),'(a,i4)') 'Layer number > ', h5_layer
!QQ      ier_msg(2) = 'FILE '//h5_infile (1:LEN(ier_msg)-5)
!QQ   ENDIF
!QQENDIF
!QQ!
!QQallocate(h5_x(1:d5_dims(1)))
!QQallocate(h5_y(1:d5_dims(2)))
!QQallocate(h5_z(1:d5_dims(3)))
!QQallocate(h5_dx(1:d5_dims(1)))
!QQallocate(h5_dy(1:d5_dims(2)))
!QQallocate(h5_dz(1:d5_dims(3)))
!QQh5_dx = 0.0D0
!QQh5_dy = 0.0D0
!QQh5_dz = 0.0D0
!QQdo i=1, d5_dims(1)
!QQ  h5_x(i) = h5_llims(1) + (i-1)*h5_steps_full(1,1)
!QQenddo
!QQdo i=1, d5_dims(2)
!QQ  h5_y(i) = h5_llims(2) + (i-1)*h5_steps_full(2,2)
!QQenddo
!QQdo i=1, d5_dims(3)
!QQ  h5_z(i) = h5_llims(3) + (i-1)*h5_steps_full(3,3)
!QQenddo
!QQ!
!QQ!write(*,'(a,l2)') ' DIRECT ', h5_direct
!QQ!write(*,'(a,i5)') ' LAYER  ', h5_layer
!QQ!write(*,'(a,3i5)') ' nndims ',    nndims
!QQ!write(*,'(a,3i5)') '   dims ', d5_dims
!QQ!write(*,'(a,3l2)') ' GRID   ', h5_is_grid, h5_has_dxyz, h5_has_dval
!QQ!write(*,'(a,3f12.5,3f8.1)') ' UNIT   ', h5_unit
!QQ!write(*,'(a,3f12.5)') ' Vector1', h5_vectors(:,1)
!QQ!write(*,'(a,3f12.5)') ' Vector2', h5_vectors(:,2)
!QQ!write(*,'(a,3f12.5)') ' Vector3', h5_vectors(:,3)
!QQ!write(*,'(a,3f12.5)') ' Steps 1', h5_steps_full(:,1)
!QQ!write(*,'(a,3f12.5)') ' Steps 2', h5_steps_full(:,2)
!QQ!write(*,'(a,3f12.5)') ' Steps 3', h5_steps_full(:,3)
!QQ!write(*,'(a,3f12.5)') ' L L B  ', h5_corners(:,1)
!QQ!write(*,'(a,3f12.5)') ' R L B  ', h5_corners(:,2)
!QQ!write(*,'(a,3f12.5)') ' L U B  ', h5_corners(:,3)
!QQ!write(*,'(a,3f12.5)') ' L L T  ', h5_corners(:,4)
!QQ!write(*,'(a,3f12.5)') ' x X    ', h5_x(1), h5_x(d5_dims(1))
!QQ!write(*,'(a,3f12.5)') ' y Y    ', h5_y(1), h5_y(d5_dims(2))
!QQ!write(*,'(a,3f12.5)') ' z Z    ', h5_z(1), h5_z(d5_dims(3))
!QQcall lib_get_use_coor(h5_vectors, h5_calc_coor, h5_use_coor)
!QQif(ier_num==0) then
!QQ!
!QQ   call dgl5_new_node
!QQ   node_number = dgl5_get_number()
!QQ   nndims = 0
!QQ   if(d5_dims(3)>1) nndims = nndims + 1
!QQ   if(d5_dims(2)>1) nndims = nndims + 1
!QQ   if(d5_dims(1)>1) nndims = nndims + 1
!QQ   dims   = d5_dims
!QQ   if(h5_direct) then
!QQ      if(nndims==3) then
!QQ         h5_data_type = H5_3D_DIRECT
!QQ      elseif(nndims==2) then
!QQ         h5_data_type = H5_2D_DIRECT
!QQ      elseif(nndims==1) then
!QQ         h5_data_type = H5_1D_DIRECT
!QQ      endif
!QQ   else
!QQ      if(nndims==3) then
!QQ         h5_data_type = H5_3D_RECI
!QQ      elseif(nndims==2) then
!QQ         h5_data_type = H5_2D_RECI
!QQ      elseif(nndims==1) then
!QQ         h5_data_type = H5_1D_RECI
!QQ      endif
!QQ   endif
!QQ   call dgl5_set_node(h5_infile, h5_data_type, h5_layer, h5_direct, nndims, d5_dims ,         &
!QQ                   h5_is_grid, h5_has_dxyz, h5_has_dval, h5_calc_coor, h5_use_coor, &
!QQ                   h5_corners, h5_vectors,&
!QQ                   h5_unit(1:3), h5_unit(4:6), h5_x, h5_y, h5_z, h5_dx, h5_dy,  &
!QQ                   h5_dz,      d5_data               , d5_sigma, h5_llims,      &
!QQ                   h5_steps, h5_steps_full)
!QQelse
!QQ   ndims = 0
!QQ   dims  = 0
!QQendif
!QQ!
!QQif(allocated(h5_datasets)) deallocate(h5_datasets)
!QQif(allocated(h5_data)) deallocate(h5_data)
!QQif(allocated(d5_data)) deallocate(d5_data)
!QQif(allocated(h5_sigma)) deallocate(h5_sigma)
!QQif(allocated(d5_sigma)) deallocate(d5_sigma)
!QQif(allocated(h5_x)) deallocate(h5_x)
!QQif(allocated(h5_y)) deallocate(h5_y)
!QQif(allocated(h5_z)) deallocate(h5_z)
!QQif(allocated(h5_dx)) deallocate(h5_dx)
!QQif(allocated(h5_dy)) deallocate(h5_dy)
!QQif(allocated(h5_dz)) deallocate(h5_dz)
!
END SUBROUTINE hdf5_read
!
!*******************************************************************************
!
subroutine nx_read_scattering_place(infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number, nndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!-
! Read a NEXUS HDF5 file in the Diffuse Developers common file format via FORPY
!+
!
use lib_nx_read_mod
!
use lib_forpython_mod
use forpy_mod
!
use iso_fortran_env, only: real64
use iso_c_binding, only:C_CHAR
!
implicit none
!
character(LEN=1024), intent(inout) :: infile
integer            , intent(in) :: length
integer            , intent(in) :: O_LAYER    ! Number optional parameter "layer:"
integer            , intent(in) :: O_TRANS    ! Number optional parameter "trans:"
integer            , intent(in) :: NOPTIONAL
character(LEN=*)   , dimension(NOPTIONAL), intent(in) :: opara
integer            , dimension(NOPTIONAL), intent(in) :: lopara
LOGICAL            , dimension(NOPTIONAL), intent(in) :: lpresent
real(KinD=PREC_DP) , dimension(NOPTIONAL), intent(in) :: owerte
!
integer, intent(out) :: node_number
integer, intent(out) :: nndims
integer, dimension(3), intent(out) :: dims
integer,                            intent(out)   :: ier_num
integer,                            intent(out)   :: ier_typ
integer,                            intent(in )   :: idims
character(LEN=*), dimension(idims), intent(inout) :: ier_msg    ! Error message
integer,                            intent(in )   :: ER_APPL
integer,                            intent(in )   :: ER_IO
integer, intent(in)    :: output_io   ! KUPLOT array size
!
!
integer :: i ,j ,k    ! Dummy loop arrays
character(len=PREC_STRING)  :: file_type   ! File type should be "Disorder scattering"
character(len=PREC_STRING)  :: file_version! File version 
character(len=PREC_STRING)  :: file_date   ! File date 
character(len=PREC_STRING)  :: file_program! File program
character(len=PREC_STRING)  :: file_author ! File author
character(len=PREC_STRING)  :: signal      ! Data contain 'measurement', 'background' ,,, ,,,
character(len=PREC_STRING)  :: radiation   ! Data were measured/calculated for this radiation
character(len=PREC_STRING)  :: is_space       ! Data are 'reciprocal or 'direct'
character(len=PREC_STRING)  :: value_type  ! Data contain this value type
character(len=PREC_STRING)  :: axes        ! string with ["h", "k","l"] or so
character(len=PREC_STRING)  :: space_group ! Hermann-Mauguin Symbol
integer                     :: symmetry_applied ! Data conform to symmetry
integer                     :: symmetry_n_mat   ! Number of Symmetry matrices
real(kind=PREC_DP)        , dimension(:,:,:), allocatable :: symmetry_mat ! Actual Symmetry matrices
!
! TEMPORARY DEBUG
!character(kind=C_CHAR, len=:), allocatable :: dname
!
call nx_read_scattering(infile, &
file_type, file_version, file_date, file_program, file_author, &
ndims, d5_dims, d5_data, &
signal, radiation, is_space, value_type, &
h5_llims, h5_steps_full, axes, &
h5_x, h5_y, h5_z, &
h5_unit(1:3), h5_unit(4:6), &
space_group, symmetry_applied, symmetry_n_mat, symmetry_mat,    &
ier_num &
)
!
!
if(    is_space=='direct' .or. is_space=='patterson') then
   h5_direct = .true.
elseif(is_space=='reciprocal') then
   h5_direct = .false.
endif
!
h5_is_grid  = .true.
h5_has_dxyz = .false.
h5_has_dval = .false.
h5_steps(1) = (h5_x(d5_dims(1))-h5_x(1))/real(d5_dims(1)-1,PREC_DP)
h5_steps(2) = (h5_y(d5_dims(2))-h5_y(1))/real(d5_dims(2)-1,PREC_DP)
h5_steps(3) = (h5_z(d5_dims(3))-h5_z(1))/real(d5_dims(3)-1,PREC_DP)
!Interim, no sigma
if(allocated(h5_sigma)) then
   allocate(d5_sigma(d5_dims(1), d5_dims(2), d5_dims(3)))
   do i=1, d5_dims(1)
      do j=1, d5_dims(2)
         do k=1, d5_dims(3)
            d5_sigma(i,j,k) = h5_sigma(k,j,i)
         enddo
      enddo
   enddo
endif
!
h5_vectors      = h5_steps_full
h5_corners(:,1) = h5_llims                                          ! Lower left
h5_corners(:,2) = h5_corners(:,1) + (d5_dims(1)-1)* h5_vectors(:,1)   ! Lower right
h5_corners(:,3) = h5_corners(:,1) + (d5_dims(2)-1)* h5_vectors(:,2)   ! Upper left
h5_corners(:,4) = h5_corners(:,1) + (d5_dims(3)-1)* h5_vectors(:,3)   ! Top left
allocate(h5_dx(ubound(h5_x,1)))
allocate(h5_dy(ubound(h5_y,1)))
allocate(h5_dz(ubound(h5_z,1)))
h5_dx = 0.0_PREC_DP
h5_dy = 0.0_PREC_DP
h5_dz = 0.0_PREC_DP
h5_infile = infile
!
call hdf5_trans_store &
                    (infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number, nndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!
deallocate(symmetry_mat)
!
end subroutine nx_read_scattering_place
!
!*******************************************************************************
!
subroutine hdf5_trans_store &
                    (infile, length, O_LAYER, O_TRANS, NOPTIONAL, opara, lopara,         &
                     lpresent, owerte,               &
                     node_number, nndims, dims, &
                     ier_num, ier_typ, idims, ier_msg, ER_APPL, ER_IO, output_io)
!+
! Transform the h5 data set into new orientation 
! place into kuplot data set
!+
!
use ber_params_mod
use lib_trans_mod
use lib_data_types_mod
use lib_use_coor_mod
use precision_mod
!
implicit none
!
!
CHARACTER(LEN=1024), INTENT(INOUT) :: infile
INTEGER            , INTENT(IN) :: length
INTEGER            , INTENT(IN) :: O_LAYER    ! Number optional parameter "layer:"
INTEGER            , INTENT(IN) :: O_TRANS    ! Number optional parameter "trans:"
INTEGER            , INTENT(IN) :: NOPTIONAL
CHARACTER(LEN=*)   , DIMENSION(NOPTIONAL), INTENT(IN) :: opara
INTEGER            , DIMENSION(NOPTIONAL), INTENT(IN) :: lopara
LOGICAL            , DIMENSION(NOPTIONAL), INTENT(IN) :: lpresent
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL), INTENT(IN) :: owerte
!
integer, intent(out) :: node_number
integer, intent(out) :: nndims
integer, dimension(3), intent(out) :: dims
INTEGER,                            INTENT(OUT)   :: ier_num
INTEGER,                            INTENT(OUT)   :: ier_typ
INTEGER,                            INTENT(IN )   :: idims
CHARACTER(LEN=*), DIMENSION(idims), INTENT(INOUT) :: ier_msg    ! Error message
INTEGER,                            INTENT(IN )   :: ER_APPL
INTEGER,                            INTENT(IN )   :: ER_IO
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!
INTEGER, PARAMETER                          :: MAXW = 1
INTEGER                                     :: ianz
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER                   , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
!REAL(KIND=PREC_DP)        , DIMENSION(3)    :: steps     ! dummy steps in H, K, L
integer :: i
!
integer, parameter :: VAL_PDF   = 14
integer, parameter :: VAL_3DPDF = 15
integer :: value  ! Intensity = 1; VAL_PDF or VAL_3DPDF
!
!                                                                          ! If data are transformed upon input use "new"
character(len=PREC_STRING)                                :: new_outfile   ! New file name  from transformed HDF5 data
integer                   , dimension(3)                  :: new_inc       ! New dimensions from transformed HDF5 data
real(kind=PREC_DP)        , dimension(3,4)                :: new_eck       ! New corners    from transformed HDF5 data
real(kind=PREC_DP)        , dimension(3,3)                :: new_vi        ! New vectors    from transformed HDF5 data
real(kind=PREC_DP)        , dimension(:,:,:), allocatable :: new_data      ! New data       from transformed HDF5 data
!
!real(kind=PREC_DP), dimension(3,3) :: temp_vi    ! Temporary copy
!integer :: extr_abs, extr_ord, extr_top
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If requested transform into new orientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(opara(O_TRANS)=='yes') then    ! Transform into different orientation
   value = 1                      ! Assume Intensity in reciprocal space
   if(h5_direct) value = VAL_PDF  ! Direct space 
   call lib_trans_menu(0, value, .false., h5_infile, d5_dims, h5_corners, h5_vectors, h5_unit(1:3),    &
           h5_unit(4:6), d5_data, VAL_PDF, VAL_3DPDF , &
        new_outfile, new_inc, new_eck, new_vi, new_data)
   deallocate(d5_data)
   if(ier_num/=0) then
!
      if(allocated(h5_datasets)) deallocate(h5_datasets)
      if(allocated(h5_data)) deallocate(h5_data)
      if(allocated(d5_data)) deallocate(d5_data)
      if(allocated(h5_sigma)) deallocate(h5_sigma)
      if(allocated(d5_sigma)) deallocate(d5_sigma)
      if(allocated(h5_x)) deallocate(h5_x)
      if(allocated(h5_y)) deallocate(h5_y)
      if(allocated(h5_z)) deallocate(h5_z)
      if(allocated(h5_dx)) deallocate(h5_dx)
      if(allocated(h5_dy)) deallocate(h5_dy)
      if(allocated(h5_dz)) deallocate(h5_dz)
      if(allocated(new_data)) deallocate(new_data)
      return
   endif
   h5_infile     = new_outfile
      infile     = new_outfile
   h5_llims(1)   = new_eck(1,1)
   h5_llims(2)   = new_eck(2,1)
   h5_llims(3)   = new_eck(3,1)
   h5_dims       = new_inc
   d5_dims       = new_inc
   h5_corners    = new_eck
   h5_vectors    = new_vi
   h5_steps_full = new_vi
   h5_steps(1)   = new_vi(1,1)
   h5_steps(2)   = new_vi(2,2)
   h5_steps(3)   = new_vi(3,3)
   allocate(d5_data(d5_dims(1), d5_dims(2), d5_dims(3)))
   d5_data       = new_data
   deallocate(new_data)
   if(allocated(d5_sigma)) deallocate(d5_sigma)
   allocate(d5_sigma(d5_dims(1), d5_dims(2), d5_dims(3)))
   d5_sigma = 0.0    ! Currently sigmas are not treated properly
endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy into KUPLOT array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IF(opara(O_LAYER)=='bottom') THEN
   h5_layer = 1
ELSEIF(opara(O_LAYER)=='middle') THEN
   if(mod(d5_dims(3),2)==0) then
      h5_layer = INT((d5_dims(3)+1)/2)
   else
      h5_layer = INT((d5_dims(3)+1)/2)
   endif
ELSEIF(opara(O_LAYER)=='top') THEN
   h5_layer = d5_dims(3)
ELSE
   cpara(1) = opara(O_LAYER)
   lpara(1) = lopara(O_LAYER)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, maxw)
   h5_layer = NINT(werte(1))
   IF(h5_layer  <=0) THEN
      ier_num = -71
      ier_typ = ER_APPL
      ier_msg(1) = 'Layer number <= 0'
      ier_msg(2) = 'FILE '//h5_infile (1:LEN(ier_msg)-5)
   ELSEIF(h5_layer  >  d5_dims(3)) THEN
      ier_num = -71
      ier_typ = ER_APPL
      WRITE(ier_msg(1),'(a,i4)') 'Layer number > ', h5_layer
      ier_msg(2) = 'FILE '//h5_infile (1:LEN(ier_msg)-5)
   ENDIF
ENDIF
!
if(.not.allocated(h5_x)) then    ! HDF5_read dowes not, NX_read does 
allocate(h5_x(1:d5_dims(1)))
allocate(h5_y(1:d5_dims(2)))
allocate(h5_z(1:d5_dims(3)))
allocate(h5_dx(1:d5_dims(1)))
allocate(h5_dy(1:d5_dims(2)))
allocate(h5_dz(1:d5_dims(3)))
h5_dx = 0.0D0
h5_dy = 0.0D0
h5_dz = 0.0D0
do i=1, d5_dims(1)
  h5_x(i) = h5_llims(1) + (i-1)*h5_steps_full(1,1)
enddo
do i=1, d5_dims(2)
  h5_y(i) = h5_llims(2) + (i-1)*h5_steps_full(2,2)
enddo
do i=1, d5_dims(3)
  h5_z(i) = h5_llims(3) + (i-1)*h5_steps_full(3,3)
enddo
endif
!
!write(*,'(a,l2)') ' DIRECT ', h5_direct
!write(*,'(a,i5)') ' LAYER  ', h5_layer
!write(*,'(a,3i5)') ' nndims ',    nndims
!write(*,'(a,3i5)') '   dims ', d5_dims
!write(*,'(a,3l2)') ' GRID   ', h5_is_grid, h5_has_dxyz, h5_has_dval
!write(*,'(a,3f12.5,3f8.1)') ' UNIT   ', h5_unit
!write(*,'(a,3f12.5)') ' Vector1', h5_vectors(:,1)
!write(*,'(a,3f12.5)') ' Vector2', h5_vectors(:,2)
!write(*,'(a,3f12.5)') ' Vector3', h5_vectors(:,3)
!write(*,'(a,3f12.5)') ' Steps 1', h5_steps_full(:,1)
!write(*,'(a,3f12.5)') ' Steps 2', h5_steps_full(:,2)
!write(*,'(a,3f12.5)') ' Steps 3', h5_steps_full(:,3)
!write(*,'(a,3f12.5)') ' L L B  ', h5_corners(:,1)
!write(*,'(a,3f12.5)') ' R L B  ', h5_corners(:,2)
!write(*,'(a,3f12.5)') ' L U B  ', h5_corners(:,3)
!write(*,'(a,3f12.5)') ' L L T  ', h5_corners(:,4)
!write(*,'(a,3f12.5)') ' x X    ', h5_x(1), h5_x(d5_dims(1))
!write(*,'(a,3f12.5)') ' y Y    ', h5_y(1), h5_y(d5_dims(2))
!write(*,'(a,3f12.5)') ' z Z    ', h5_z(1), h5_z(d5_dims(3))
call lib_get_use_coor(h5_vectors, h5_calc_coor, h5_use_coor)
if(ier_num==0) then
!
   call dgl5_new_node
   node_number = dgl5_get_number()
   nndims = 0
   if(d5_dims(3)>1) nndims = nndims + 1
   if(d5_dims(2)>1) nndims = nndims + 1
   if(d5_dims(1)>1) nndims = nndims + 1
   dims   = d5_dims
   if(h5_direct) then
      if(nndims==3) then
         h5_data_type = H5_3D_DIRECT
      elseif(nndims==2) then
         h5_data_type = H5_2D_DIRECT
      elseif(nndims==1) then
         h5_data_type = H5_1D_DIRECT
      endif
   else
      if(nndims==3) then
         h5_data_type = H5_3D_RECI
      elseif(nndims==2) then
         h5_data_type = H5_2D_RECI
      elseif(nndims==1) then
         h5_data_type = H5_1D_RECI
      endif
   endif
   call dgl5_set_node(h5_infile, h5_data_type, h5_layer, h5_direct, nndims, d5_dims ,         &
                   h5_is_grid, h5_has_dxyz, h5_has_dval, h5_calc_coor, h5_use_coor, &
                   h5_corners, h5_vectors,&
                   h5_unit(1:3), h5_unit(4:6), h5_x, h5_y, h5_z, h5_dx, h5_dy,  &
                   h5_dz,      d5_data               , d5_sigma, h5_llims,      &
                   h5_steps, h5_steps_full)
else
   ndims = 0
   dims  = 0
endif
!
if(allocated(h5_datasets)) deallocate(h5_datasets)
if(allocated(h5_data)) deallocate(h5_data)
if(allocated(d5_data)) deallocate(d5_data)
if(allocated(h5_sigma)) deallocate(h5_sigma)
if(allocated(d5_sigma)) deallocate(d5_sigma)
if(allocated(h5_x)) deallocate(h5_x)
if(allocated(h5_y)) deallocate(h5_y)
if(allocated(h5_z)) deallocate(h5_z)
if(allocated(h5_dx)) deallocate(h5_dx)
if(allocated(h5_dy)) deallocate(h5_dy)
if(allocated(h5_dz)) deallocate(h5_dz)
!
end subroutine hdf5_trans_store
!
!*******************************************************************************
!
INTEGER FUNCTION op_func(loc_id, name, info, operator_data) bind(C)
     
    USE allocate_generic
    USE HDF5
    USE ISO_C_BINDING
use trig_degree_mod
    IMPLICIT NONE
     
    INTEGER, PARAMETER ::MAXSTR = 1024
    INTEGER(HID_T), VALUE :: loc_id
    CHARACTER(LEN=1), DIMENSION(1:MAXSTR) :: name ! must have LEN=1 for bind(C) strings
    TYPE(C_PTR) :: info
    TYPE(C_PTR) :: operator_data
     
    INTEGER   :: status, i, length
 
    TYPE(H5O_info_t), TARGET :: infobuf
!   TYPE(C_PTR) :: ptr
    CHARACTER(LEN=MAXSTR) :: name_string
!
    INTEGER :: all_status ! allocation status
    INTEGER :: ndata      ! upper limit for allocation
 
    !
    ! Get type of the object and display its name and type.
    ! The name of the object is passed to this FUNCTION by
    ! the Library.
    !
 
    DO i = 1, MAXSTR
       name_string(i:i) = name(i)(1:1)
    ENDDO
 
    CALL H5Oget_info_by_name_f(loc_id, name_string, infobuf, status)
 
    ! Include the string up to the C NULL CHARACTER
    length = 0
    DO
       IF(length>=MAXSTR) EXIT
       IF(name_string(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.MAXSTR) EXIT
       length = length + 1
    ENDDO
 
    IF(infobuf%type.EQ.H5O_TYPE_GROUP_F)THEN
!      WRITE(*,*) "Group: ", name_string(1:length)
       CONTINUE
    ELSE IF(infobuf%type.EQ.H5O_TYPE_DATASET_F)THEN
!      WRITE(*,*) "Dataset: ", name_string(1:length)
       IF(h5_n_datasets==H5_MAX_DATASETS) THEN
          ndata = H5_MAX_DATASETS + 10
          CALL alloc_arr(h5_datasets, 1, ndata, all_status, ' ' )
          H5_MAX_DATASETS = H5_MAX_DATASETS + 10
       ENDIF 
       h5_n_datasets = h5_n_datasets + 1
       h5_datasets(h5_n_datasets) = name_string(1:length)
    ELSE IF(infobuf%type.EQ.H5O_TYPE_NAMED_DATATYPE_F)THEN
!      WRITE(*,*) "Datatype: ", name_string(1:length)
       CONTINUE
    ELSE
!      WRITE(*,*) "Unknown: ", name_string(1:length)
       CONTINUE
    ENDIF
 
    op_func = 0 ! return successful
 
  END FUNCTION op_func
!
!*******************************************************************************
!
SUBROUTINE hdf5_error(infile, file_id, dataname, dset_id, ier_num, ier_typ, &
                      idims,ier_msg,er_nr, er_type, message)
!
IMPLICIT NONE
!
CHARACTER(LEN=*)       , INTENT(IN) :: infile     ! File name
INTEGER(KIND=LIB_HID_T), INTENT(IN) :: file_id    ! File identifier
CHARACTER(LEN=*)       , INTENT(IN) :: dataname   ! dataset name
INTEGER(KIND=LIB_HID_T), INTENT(IN) :: dset_id    ! dataset identifier
INTEGER            , INTENT(OUT) :: ier_num     ! Error number
INTEGER            , INTENT(OUT) :: ier_typ     ! Error type
INTEGER            , INTENT(IN)  :: idims       ! Dimension of ier_msg
CHARACTER(LEN=*), DIMENSION(idims)   , INTENT(OUT) :: ier_msg    ! Error message
INTEGER            , INTENT(IN) :: er_nr      ! Error number
INTEGER            , INTENT(IN) :: er_type    ! Error type
CHARACTER(LEN=*)   , INTENT(IN) :: message    ! Error message
!
INTEGER             :: hdferr      ! Error number
!
ier_num = er_nr             ! Copy error numbers into errlist_mod
ier_typ = er_type
ier_msg(1) = message(1:LEN(ier_msg))
IF(dataname /= ' ') ier_msg(2) = 'Dataset: '//dataname
IF(infile   /= ' ') ier_msg(3) = 'H5 File: '//infile
!
IF(ALLOCATED(h5_datasets)) DEALLOCATE(h5_datasets)
IF(ALLOCATED(h5_data))     DEALLOCATE(h5_data)
!
! CLOSE up, ignore error messages
IF(dset_id /= 0) CALL H5Dclose_f(dset_id , hdferr)                         ! Close dataset
IF(file_id /= 0) CALL h5fclose_f(file_id, hdferr)                          ! Close the input file
CALL H5close_f(hdferr)                                    ! Close HDF interface
!
END SUBROUTINE hdf5_error
!
!*******************************************************************************
!
end module lib_hdf5_read_mod
