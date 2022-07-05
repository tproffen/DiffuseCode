MODULE kuplot_load_h5
!-
!  Contains routines to read a HDF5 file written by DISCUS
!+
USE hdf5
!
USE kuplot_config
!
use hdf5_def_mod
USE precision_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC hdf5_read
PUBLIC hdf5_new_node
PUBLIC hdf5_set_node
PUBLIC hdf5_place_kuplot
PUBLIC hdf5_set_pointer
PUBLIC hdf5_find_node
PUBLIC hdf5_copy_node
PUBLIC hdf5_get_layer
PUBLIC hdf5_get_height
PUBLIC hdf5_get_direct
PUBLIC hdf5_get_dims
PUBLIC hdf5_get_llims
PUBLIC hdf5_get_steps
PUBLIC hdf5_get_map
PUBLIC hdf5_set_map
PUBLIC hdf5_get_tmap
PUBLIC hdf5_reset
!
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: h5_datasets       ! Names of the data set in file
CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
INTEGER, DIMENSION(MAXKURVTOT)                        :: h5_ku_is_h5 = 0   ! Pointer from kuplot number to h5 number
INTEGER, DIMENSION(MAXKURVTOT)                        :: h5_h5_is_ku = 0   ! Pointer from h5 number to kuplot number
INTEGER                                               :: h5_number   = 0   ! Currently loaded h5 data sets
!
TYPE :: h5_data_struc
   INTEGER                                               :: h5_data_num       ! Current data set number
   CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
   INTEGER                                               :: h5_layer=1        ! Current layer in data set
   LOGICAL                                               :: h5_direct         ! Direct space == TRUE
   INTEGER                                               :: ndims             ! Number of dimensions
   INTEGER                                               :: one_ndims         ! Number of dimensions
   INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: h5_dims           ! Actual dimensions
   INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: maxdims           ! Maximum dimensions
   INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: one_dims          ! Actual dimensions
   INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: one_maxdims       ! Maximum dimensions
   REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_data           ! Actual diffraction data
   REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_llims          ! Lower limits
   REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_steps          ! steps in H, K, L
   REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_steps_full     ! steps in H, K, L
   TYPE(h5_data_struc), POINTER                          :: after
END TYPE h5_data_struc
!
TYPE(h5_data_struc), POINTER                          :: h5_root => NULL()
TYPE(h5_data_struc), POINTER                          :: h5_temp => NULL()
TYPE(h5_data_struc), POINTER                          :: h5_find => NULL()
!
INTEGER                                               :: H5_MAX_DATASETS   ! Current MAX data sets
INTEGER                                               :: h5_n_datasets     ! Current actual data sets
INTEGER                                               :: h5_layer=1        ! Current layer in data set
LOGICAL                                               :: h5_direct         ! Direct space == TRUE
INTEGER                                               :: ndims             ! Number of dimensions
INTEGER                                               :: one_ndims         ! Number of dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: h5_dims           ! Actual dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: maxdims           ! Maximum dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: one_dims          ! Actual dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3)                   :: one_maxdims       ! Maximum dimensions
REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_data           ! Actual diffraction data
REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_llims          ! Lower limits
REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_steps          ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_steps_full     ! steps in H, K, L
!                                                                          ! 1.st dim: [hkl], 2nd:[abs, ord, top]
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
INTEGER(KIND=2), TARGET :: r_is_direct
!INTEGER(KIND=2), DIMENSION(10), TARGET :: rstring
!INTEGER(HSIZE_T) :: nummer
!
INTEGER, PARAMETER                          :: MAXW = 1
INTEGER                                     :: ianz
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER                   , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP)        , DIMENSION(3)    :: steps     ! dummy steps in H, K, L
INTEGER :: i,j                      ! Dummy loop indices
INTEGER :: nlayer                   ! Layer to place into KUPLOT
!
h5_infile = infile
dataname = ' '
!dset_id  = 0
!file_id  = 0
CALL hdf5_new_node
!
!
H5_MAX_DATASETS = 10                                        ! Initial estimate of dataset number
ALLOCATE(h5_datasets(H5_MAX_DATASETS))
h5_n_datasets = 0                                           ! Currently no datasets found
!
h5_dims    = 1
maxdims = 1
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
dataname='data'
CALL H5Dopen_f(file_id, dataname, dset_id, hdferr)           ! Open the dataset
CALL H5Dget_space_f(dset_id, space_id, hdferr)
CALL H5Sget_simple_extent_ndims_f(space_id, ndims, hdferr)   ! Get the number of dimensions in data set
CALL H5Sget_simple_extent_dims_f(space_id, h5_dims, maxdims, hdferr)   ! Get the dimensions in data set
!write(*,*) ' h5_dims ', h5_dims
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
h5_direct = 1 == r_is_direct
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
!write(*,*) ' H5_steps      ', h5_steps
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
!write(*,*) ' H5_steps abs  ', h5_steps_full(:,1)
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
!write(*,*) ' H5_steps ord  ', h5_steps_full(:,2)
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
!write(*,*) ' H5_steps top  ', h5_steps_full(:,3)
endif
!
CALL h5fclose_f(file_id, hdferr)                             ! Close the input file
!
CALL H5close_f(hdferr)                                    ! Close HDF interface
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy into H5 storage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CALL hdf5_set_node(h5_infile, h5_layer, h5_direct, ndims, one_ndims, h5_dims, &
                   maxdims, one_dims, one_maxdims, h5_data, h5_llims, h5_steps,&
                   h5_steps_full)
 
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
      ier_msg(2) = 'FILE '//h5_infile (1:LEN(ier_msg)-5)
      DEALLOCATE(h5_datasets)
      RETURN
   ELSEIF(nlayer  >  h5_dims(1)) THEN
      ier_num = -71
      ier_typ = ER_APPL
      WRITE(ier_msg(1),'(a,i4)') 'Layer number > ', nlayer
      ier_msg(2) = 'FILE '//h5_infile (1:LEN(ier_msg)-5)
      DEALLOCATE(h5_datasets)
      RETURN
   ENDIF
ENDIF
!
if(h5_temp%h5_dims(1)==1 .and. h5_temp%h5_dims(2)==1) then
   CALL hdf5_place_kuplot_1d(nlayer, .TRUE.,.TRUE., .TRUE.,               &
      MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
      xmin, xmax, ymin, ymax, &
      offxy, offz, lni, lh5, lenc, ier_num, ier_typ, output_io)
else
   CALL hdf5_place_kuplot(nlayer, .TRUE.,.TRUE., .TRUE.,               &
      MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
      xmin, xmax, ymin, ymax, &
      offxy, offz, lni, lh5, lenc, ier_num, ier_typ, output_io)
endif
!
DEALLOCATE(h5_datasets)
DEALLOCATE(h5_data)
!
END SUBROUTINE hdf5_read
!
!*******************************************************************************
!
SUBROUTINE hdf5_new_node
!-
! Create a new node
!+
IF(ASSOCIATED(h5_root)) THEN                                ! A root node exists
   h5_temp => h5_root                                       ! Point to current node
   find_node: DO WHILE(ASSOCIATED(h5_temp%after))           ! Does next node exist?
      h5_temp => h5_temp%after                              ! Next node exists, point to this next node
   ENDDO find_node
   ALLOCATE(h5_temp%after)                                  ! Create next node
   h5_temp => h5_root%after                                 ! Point to Current working node
ELSE
   ALLOCATE(h5_root)                                        ! Create first node
   NULLIFY(h5_root%after)
   h5_temp => h5_root                                       ! Point to Current working node
ENDIF
! Work on current node
NULLIFY(h5_temp%after)
h5_temp%h5_data_num = h5_number + 1                         ! Increment the data number
h5_number = h5_number + 1                                   ! Increment the global data number
!
END SUBROUTINE hdf5_new_node
!
!*******************************************************************************
!
SUBROUTINE hdf5_copy_node(old, new)
!-
!  Copies old node to new node
!+
INTEGER, INTENT(IN)  :: old
INTEGER, INTENT(OUT) :: new
INTEGER :: ier_num
INTEGER :: ier_typ
!
CALL hdf5_find_node(old, ier_num, ier_typ)
CALL hdf5_new_node
h5_temp%h5_infile   = h5_find%h5_infile         ! input file
h5_temp%h5_layer    = h5_find%h5_layer          ! Current layer in data set
h5_temp%h5_direct   = h5_find%h5_direct         ! Direct space == TRUE
h5_temp%ndims       = h5_find%ndims          ! Number of dimensions
h5_temp%one_ndims   = h5_find%one_ndims      ! Number of dimensions
h5_temp%h5_dims     = h5_find%h5_dims           ! Actual dimensions
h5_temp%maxdims     = h5_find%maxdims        ! Maximum dimensions
h5_temp%one_dims    = h5_find%one_dims       ! Actual dimensions
h5_temp%one_maxdims = h5_find%one_maxdims    ! Maximum dimensions
h5_temp%h5_data     = h5_find%h5_data           ! Actual diffraction data
h5_temp%h5_llims    = h5_find%h5_llims          ! Lower limits
h5_temp%h5_steps    = h5_find%h5_steps          ! steps in H, K, L
!
new = h5_temp%h5_data_num
!
END SUBROUTINE hdf5_copy_node
!
!*******************************************************************************
!
SUBROUTINE hdf5_set_node(l_infile, l_layer, l_direct, l_ndims, l_one_ndims, l_dims, &
                   l_maxdims, l_one_dims, l_one_maxdims, l_data, l_llims, l_steps,  &
                   l_steps_full)
!-
!  Place the temporary values into the current hdf5 node
!+
CHARACTER(LEN=*)                   , INTENT(IN)       :: l_infile         ! Input file
INTEGER                            , INTENT(IN)       :: l_layer          ! Current layer in data set
LOGICAL                            , INTENT(IN)       :: l_direct         ! Direct space == TRUE
INTEGER                            , INTENT(IN)       :: l_ndims          ! Number of dimensions
INTEGER                            , INTENT(IN)       :: l_one_ndims      ! Number of dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3), INTENT(IN)       :: l_dims           ! Actual dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3), INTENT(IN)       :: l_maxdims        ! Maximum dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3), INTENT(IN)       :: l_one_dims       ! Actual dimensions
INTEGER(KIND=HSIZE_T), DIMENSION(3), INTENT(IN)       :: l_one_maxdims    ! Maximum dimensions
REAL(KIND=PREC_SP)   , DIMENSION(l_dims(1), l_dims(2), l_dims(3)), INTENT(IN):: l_data           ! Actual diffraction data
REAL(KIND=PREC_DP)   , DIMENSION(3), INTENT(IN)       :: l_llims          ! Lower limits
REAL(KIND=PREC_DP)   , DIMENSION(3), INTENT(IN)       :: l_steps          ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3),INTENT(IN)       :: l_steps_full     ! steps in H, K, L
!
h5_temp%h5_infile   = l_infile         ! input file
h5_temp%h5_layer    = l_layer          ! Current layer in data set
h5_temp%h5_direct   = l_direct         ! Direct space == TRUE
h5_temp%ndims       = l_ndims          ! Number of dimensions
h5_temp%one_ndims   = l_one_ndims      ! Number of dimensions
h5_temp%h5_dims     = l_dims           ! Actual dimensions
h5_temp%maxdims     = l_maxdims        ! Maximum dimensions
h5_temp%one_dims    = l_one_dims       ! Actual dimensions
h5_temp%one_maxdims = l_one_maxdims    ! Maximum dimensions
!ALLOCATE(h5_temp%h5_data(h5_temp%h5_dims(1), h5_temp%h5_dims(2), h5_temp%h5_dims(3)))
h5_temp%h5_data     = l_data           ! Actual diffraction data
h5_temp%h5_llims    = l_llims          ! Lower limits
h5_temp%h5_steps    = l_steps          ! steps in H, K, L
h5_temp%h5_steps_full    = l_steps_full          ! steps in H, K, L
!
END SUBROUTINE hdf5_set_node
!
!*******************************************************************************
!
INTEGER FUNCTION op_func(loc_id, name, info, operator_data) bind(C)
     
    USE allocate_generic
    USE HDF5
    USE ISO_C_BINDING
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
    INTEGER :: size_d     ! allocated size
 
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
          CALL alloc_arr(h5_datasets, 1, ndata, all_status, ' ', size_d)
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
CHARACTER(LEN=*)   , INTENT(IN) :: infile     ! File name
INTEGER(KIND=HID_T), INTENT(IN) :: file_id    ! File identifier
CHARACTER(LEN=*)   , INTENT(IN) :: dataname   ! dataset name
INTEGER(KIND=HID_T), INTENT(IN) :: dset_id    ! dataset identifier
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
SUBROUTINE hdf5_place_kuplot(nlayer, lset, lnew, lshow,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
   xmin, xmax, ymin, ymax, &
   offxy, offz, lni, lh5, lenc, ier_num, ier_typ, output_io)
!
!-
! PLace a curve into the kuplot section, 
! IF lset==TRUE set absolute layer , else increment
! IF lnew==TRUE, make new curve, 
! IF lshow = TRUE display data
!+
!
use kuplot_show_mod
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
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lni
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lh5
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: lenc
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
INTEGER,                 INTENT(OUT) :: ier_num
INTEGER,                 INTENT(OUT) :: ier_typ
!
INTEGER :: i,j,k, ll             ! dummy indices
INTEGER :: izz
INTEGER :: node_number
!
IF(lnew) THEN            ! This is a new data set, from 'load' command
   izz = iz
   h5_h5_is_ku(h5_number) = izz
   h5_ku_is_h5(izz      ) = h5_number
   node_number = h5_number
ELSE                     ! Overwrite current KUPLOT data set
   izz = iz - 1
ENDIF
!                        ! Locate this data set in the h5 storage
IF(.NOT. ASSOCIATED(h5_root)) THEN
   ier_num = -74         ! Root node does not exist !
   ier_typ =   6         ! ER_APPL
   RETURN
ENDIF
CALL hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
!
IF(lset) THEN
  h5_temp%h5_layer = nlayer
ELSE
  h5_temp%h5_layer = MAX(1,MIN(INT(h5_temp%h5_dims(1)), h5_temp%h5_layer+nlayer))
ENDIF
ll = 0
k = h5_temp%h5_layer
DO i = 1, h5_temp%h5_dims(3)
  DO j = 1, h5_temp%h5_dims(2)
     ll = ll + 1
     z(offz(izz - 1) + ll ) = h5_temp%h5_data(k,j,i)
  ENDDO
ENDDO
!
nx(izz) = h5_temp%h5_dims(3)
ny(izz) = h5_temp%h5_dims(2)
xmin(izz) = h5_temp%h5_llims(1)
xmax(izz) = h5_temp%h5_llims(1) + (nx(izz)-1)*h5_temp%h5_steps(1)
ymin(izz) = h5_temp%h5_llims(2)
ymax(izz) = h5_temp%h5_llims(2) + (ny(izz)-1)*h5_temp%h5_steps(2)
!
DO i = 1, nx(izz)
    x(offxy(izz - 1) + i) = xmin(izz) + (i - 1) * h5_temp%h5_steps(1)
ENDDO
DO i = 1, ny(izz)
   y(offxy(izz - 1) + i) = ymin(izz) + (i - 1) * h5_temp%h5_steps(2)
ENDDO
lni (izz) = .TRUE.
lh5 (izz) = .TRUE.
lenc(izz) = MAX(nx(izz), ny(izz))
offxy(izz) = offxy(izz - 1) + lenc(izz)
offz (izz) = offz (izz - 1) + nx(izz) * ny(izz)
fname(izz) = h5_temp%h5_infile(1:LEN_TRIM(h5_temp%h5_infile))
h5_h5_is_ku(node_number) = izz       ! H5 Data set 1 is stored in Kuplot as number izz
h5_ku_is_h5(izz        ) = node_number ! Kuplot data set izz is stored in H5 number 1
IF(lnew) iz = iz + 1
!
IF(lshow) THEN
   CALL show_data(iz - 1)!
   WRITE(output_io,1000) h5_temp%h5_dims(3), h5_temp%h5_dims(2), h5_temp%h5_dims(1)
   WRITE(output_io,1100) nlayer
   1000 FORMAT('   Full size:', 2(i7,' x'), i7, ' points')
   1100 FORMAT('   At  layer:',   i7      ,/)
ENDIF
!
END SUBROUTINE hdf5_place_kuplot
!
!*******************************************************************************
!
SUBROUTINE hdf5_place_kuplot_1d(nlayer, lset, lnew, lshow,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
   xmin, xmax, ymin, ymax, &
   offxy, offz, lni, lh5, lenc, ier_num, ier_typ, output_io)
!
!-
! PLace a 1D curve into the kuplot section, 
! IF lset==TRUE set absolute layer , else increment
! IF lnew==TRUE, make new curve, 
! IF lshow = TRUE display data
!+
!
use kuplot_show_mod
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
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lni
LOGICAL, DIMENSION(maxkurvtot), INTENT(INOUT) :: lh5
INTEGER, DIMENSION(MAXKURVTOT), INTENT(INOUT) :: lenc
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
INTEGER,                 INTENT(OUT) :: ier_num
INTEGER,                 INTENT(OUT) :: ier_typ
!
INTEGER :: i,j,k, ll             ! dummy indices
INTEGER :: izz
INTEGER :: node_number
!
IF(lnew) THEN            ! This is a new data set, from 'load' command
   izz = iz
   h5_h5_is_ku(h5_number) = izz
   h5_ku_is_h5(izz      ) = h5_number
   node_number = h5_number
ELSE                     ! Overwrite current KUPLOT data set
   izz = iz - 1
ENDIF
!                        ! Locate this data set in the h5 storage
IF(.NOT. ASSOCIATED(h5_root)) THEN
   ier_num = -74         ! Root node does not exist !
   ier_typ =   6         ! ER_APPL
   RETURN
ENDIF
CALL hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
!
IF(lset) THEN
  h5_temp%h5_layer = nlayer
ELSE
  h5_temp%h5_layer = MAX(1,MIN(INT(h5_temp%h5_dims(1)), h5_temp%h5_layer+nlayer))
ENDIF
!
lenc(izz) = h5_temp%h5_dims(3)
xmin(izz) = h5_temp%h5_llims(1)
xmax(izz) = h5_temp%h5_llims(1) + (lenc(izz)-1)*h5_temp%h5_steps(1)
DO i = 1, h5_temp%h5_dims(3)
  x(offxy(izz - 1) + i) = xmin(izz) + (i - 1) * h5_temp%h5_steps(1)
  y(offxy(izz - 1) + i) = h5_temp%h5_data(1,1,i)
ENDDO
!
ymin(izz) = minval(h5_temp%h5_data(1,1,1:h5_temp%h5_dims(3)))
ymax(izz) = maxval(h5_temp%h5_data(1,1,1:h5_temp%h5_dims(3)))
!
lni (izz) = .false.
lh5 (izz) = .true. 
offxy(izz) = offxy(izz - 1) + lenc(izz)
fname(izz) = h5_temp%h5_infile(1:LEN_TRIM(h5_temp%h5_infile))
h5_h5_is_ku(node_number) = izz       ! H5 Data set 1 is stored in Kuplot as number izz
h5_ku_is_h5(izz        ) = node_number ! Kuplot data set izz is stored in H5 number 1
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
SUBROUTINE hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
!-
!  Find the node associated to kuplot data set number izz
!+
INTEGER, INTENT(IN ) :: izz
INTEGER, INTENT(OUT) :: ier_num
INTEGER, INTENT(OUT) :: ier_typ
INTEGER, INTENT(OUT) :: node_number
!
h5_temp => h5_root
find_node: DO            ! Search for node
   IF(h5_ku_is_h5(izz) == h5_temp%h5_data_num) THEN
      node_number = h5_temp%h5_data_num
      EXIT find_node
   ELSE
      IF(ASSOCIATED(h5_temp%after)) THEN    ! A next node exists 
         h5_temp => h5_temp%after
      ELSE
         ier_num = -74         ! Root node does not exist !
         ier_typ =   6         ! ER_APPL
         RETURN
      ENDIF
   ENDIF
ENDDO find_node
!
END SUBROUTINE hdf5_set_pointer
!
!*******************************************************************************
!
SUBROUTINE hdf5_find_node(node_number, ier_num, ier_typ)
!-
!  Find the node with node_number
!+
INTEGER, INTENT(IN)  :: node_number
INTEGER, INTENT(OUT) :: ier_num
INTEGER, INTENT(OUT) :: ier_typ
!
h5_find => h5_root
find_node: DO            ! Search for node
   IF(node_number == h5_find%h5_data_num) THEN
      EXIT find_node
   ELSE
      IF(ASSOCIATED(h5_find%after)) THEN    ! A next node exists 
         h5_find => h5_find%after
      ELSE
         ier_num = -74         ! Root node does not exist !
         ier_typ =   6         ! ER_APPL
         RETURN
      ENDIF
   ENDIF
ENDDO find_node
!
END SUBROUTINE hdf5_find_node
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_layer()
!
IMPLICIT NONE
!
hdf5_get_layer = h5_temp%h5_layer
!
END FUNCTION hdf5_get_layer
!
!*******************************************************************************
!
LOGICAL FUNCTION hdf5_get_direct()
!
IMPLICIT NONE
!
hdf5_get_direct = h5_temp%h5_direct
!
END FUNCTION hdf5_get_direct
!
!*******************************************************************************
!
REAL FUNCTION hdf5_get_height()
!
IMPLICIT NONE
!
hdf5_get_height = h5_temp%h5_llims(3) + (h5_temp%h5_layer-1)*h5_temp%h5_steps(3)
!
END FUNCTION hdf5_get_height
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
dims = h5_temp%h5_dims
!
END SUBROUTINE hdf5_get_dims
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_llims(idata, llims)
!
use precision_mod
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN)  :: idata
real(kind=PREC_DP), DIMENSION(3), INTENT(OUT) :: llims
!
llims = h5_temp%h5_llims
!
END SUBROUTINE hdf5_get_llims
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_steps(idata, steps)
!
use hdf5_def_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN)  :: idata
real(kind=PREC_DP), DIMENSION(3,3), INTENT(OUT) :: steps
!
!write(*,*) yd_present(YD_step_sizes_abs:YD_step_sizes_TOP), &
!       ALL(yd_present(YD_step_sizes_abs:YD_step_sizes_TOP))
if(ALL(yd_present(YD_step_sizes_abs:YD_step_sizes_TOP))) then
   steps = h5_temp%h5_steps_full
else
   steps = 0.0D0
   steps(1,1) = h5_temp%h5_steps(1)
   steps(2,2) = h5_temp%h5_steps(2)
   steps(3,3) = h5_temp%h5_steps(3)
endif
!write(*,*) ' GETS ', steps(:,1)
!write(*,*) ' GETS ', steps(:,2)
!write(*,*) ' GETS ', steps(:,3)
!
END SUBROUTINE hdf5_get_steps
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
INTEGER :: i,j,k
!
DO i=1, dims(1)
   DO j=1, dims(2)
      DO k=1, dims(3)
         odata(i,j,k) = h5_temp%h5_data(i,j,k)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE hdf5_get_map
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_tmap(dims, odata)
!-
! Get the matrix in transposed form, which is the regular Fortran style
!+
!
IMPLICIT NONE
!
INTEGER,            DIMENSION(3),                         INTENT(IN)  :: dims
REAL(KIND=PREC_DP), DIMENSION(dims(1), dims(2), dims(3)), INTENT(OUT) :: odata
!
INTEGER :: i,j,k
!
DO i=1, dims(1)
   DO j=1, dims(2)
      DO k=1, dims(3)
         odata(i,j,k) = h5_temp%h5_data(k,j,i)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE hdf5_get_tmap
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
INTEGER :: i,j,k
!
DO i=1, dims(1)
   DO j=1, dims(2)
      DO k=1, dims(3)
         h5_temp%h5_data(i,j,k)= odata(i,j,k) 
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE hdf5_set_map
!
!*******************************************************************************
!
SUBROUTINE hdf5_reset
!
TYPE(h5_data_struc), POINTER :: h5_current => NULL()
!
IF(ASSOCIATED(h5_root)) THEN       ! A storage does exist
   h5_temp => h5_root
   IF(ALLOCATED(h5_temp%h5_data)) DEALLOCATE(h5_temp%h5_data)
   find_node: DO 
      IF(ASSOCIATED(h5_temp%after)) THEN   ! A next node exists
         h5_current => h5_temp             ! Point to current
         h5_temp    => h5_temp%after       ! Point to next node
         DEALLOCATE(h5_current)            ! Clean up current node
      ELSE
         h5_current => h5_temp             ! Point to current
         DEALLOCATE(h5_current)            ! Clean up current node
         EXIT find_node                    ! We are done
      ENDIF
   ENDDO find_node
ENDIF
NULLIFY(h5_temp)
NULLIFY(h5_root)
h5_number   = 0
h5_h5_is_ku = 0
h5_ku_is_h5 = 0
!
END SUBROUTINE hdf5_reset
!
!*******************************************************************************
!
END MODULE kuplot_load_h5
