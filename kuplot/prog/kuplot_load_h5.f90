MODULE kuplot_load_h5
!-
!  Contains routines to read a HDF5 file written by DISCUS
!+
USE hdf5
!
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
!
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: h5_datasets       ! Names of the data set in file
CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
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
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
!INTEGER(KIND=SIZE_T), PARAMETER :: f_str_len = 20  ! Length F-string to receive  
!
CHARACTER(LEN=12)   :: dataname    ! Dummy name for HDF5 datasets
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
INTEGER :: i,j                      ! Dummy loop indices
INTEGER :: nlayer                   ! Layer to place into KUPLOT
!
h5_infile = infile
dataname = ' '
!dset_id  = 0
!file_id  = 0
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
!DO i=1,h5_n_datasets
!   WRITE(*,*) ' DATASETS ', h5_datasets(i)(1:LEN_TRIM(h5_datasets(i)))
!ENDDO
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
CALL H5Dread_F(dset_id, H5T_STD_I32LE, f_ptr, hdferr)
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
!
CALL h5fclose_f(file_id, hdferr)                             ! Close the input file
!
CALL H5close_f(hdferr)                                    ! Close HDF interface
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
CALL hdf5_place_kuplot(nlayer, .TRUE.,.TRUE., .TRUE.,               &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, &
   xmin, xmax, ymin, ymax, &
   offxy, offz, lni, lh5, lenc, ier_num, ier_typ, output_io)
!
DEALLOCATE(h5_datasets)
!
END SUBROUTINE hdf5_read
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
INTEGER, INTENT(IN)    :: output_io   ! KUPLOT array size
!
INTEGER,                 INTENT(OUT) :: ier_num
INTEGER,                 INTENT(OUT) :: ier_typ
!
INTEGER :: i,j,k, ll             ! dummy indices
INTEGER :: izz
!
IF(lset) THEN
  h5_layer = nlayer
ELSE
  h5_layer = MAX(1,MIN(INT(h5_dims(1)), h5_layer+nlayer))
ENDIF
IF(lnew) THEN
   izz = iz
ELSE
   izz = iz - 1
ENDIF
ll = 0
k = h5_layer
DO i = 1, h5_dims(3)
  DO j = 1, h5_dims(2)
     ll = ll + 1
     z(offz(izz - 1) + ll ) = h5_data(k,j,i)
  ENDDO
ENDDO
nx(izz) = h5_dims(3)
ny(izz) = h5_dims(2)
xmin(izz) = h5_llims(1)
xmax(izz) = h5_llims(1) + (nx(izz)-1)*h5_steps(1)
ymin(izz) = h5_llims(2)
ymax(izz) = h5_llims(2) + (ny(izz)-1)*h5_steps(2)
DO i = 1, nx(izz)
   x(offxy(izz - 1) + i) = xmin(izz) + (i - 1) * h5_steps(1)
ENDDO
DO i = 1, ny(izz)
   y(offxy(izz - 1) + i) = ymin(izz) + (i - 1) * h5_steps(2)
ENDDO
lni (izz) = .TRUE.
lh5 (izz) = .TRUE.
lenc(izz) = MAX(nx(izz), ny(izz))
offxy(izz) = offxy(izz - 1) + lenc(izz)
offz (izz) = offz (izz - 1) + nx(izz) * ny(izz)
fname(izz) = h5_infile(1:LEN_TRIM(h5_infile))
IF(lnew) iz = iz + 1
!
IF(lshow) THEN
   CALL show_data(iz - 1)!
   WRITE(output_io,1000) h5_dims(3), h5_dims(2), h5_dims(1)
   WRITE(output_io,1100) nlayer
   1000 FORMAT('   Full size:', 2(i7,' x'), i7, ' points')
   1100 FORMAT('   At  layer:',   i7      ,/)
ENDIF
!
END SUBROUTINE hdf5_place_kuplot
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_layer()
!
IMPLICIT NONE
!
hdf5_get_layer = h5_layer
!
END FUNCTION hdf5_get_layer
!
!*******************************************************************************
!
LOGICAL FUNCTION hdf5_get_direct()
!
IMPLICIT NONE
!
hdf5_get_direct = h5_direct
!
END FUNCTION hdf5_get_direct
!
!*******************************************************************************
!
REAL FUNCTION hdf5_get_height()
!
IMPLICIT NONE
!
hdf5_get_height = h5_llims(3) + (h5_layer-1)*h5_steps(3)
!
END FUNCTION hdf5_get_height
!
!*******************************************************************************
!
END MODULE kuplot_load_h5
