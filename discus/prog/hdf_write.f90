!MODULE hdf_write_mod
!
!*****7*****************************************************************
!
!CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE hdf5_write (value, laver, outfile, out_inc, out_eck, out_vi, &
                       cr_a0, cr_win, qval, VAL_PDF, VAL_3DPDF, valmax, &
                       ier_num, ier_typ, ER_IO, ER_APPL)
!
USE hdf5
!
IMPLICIT NONE
!
INTEGER, PARAMETER:: PREC_DP=SELECTED_REAL_KIND(p=15,r=307)  ! double precision
!
INTEGER, INTENT(IN) :: value
LOGICAL, INTENT(IN) :: laver
CHARACTER(LEN=200), INTENT(IN) :: outfile
INTEGER, DIMENSION(3)  , INTENT(IN) :: out_inc
REAL   , DIMENSION(3,4), INTENT(IN) :: out_eck ! (3,4)
REAL   , DIMENSION(3,3), INTENT(IN) :: out_vi 
REAL   , DIMENSION(3)  , INTENT(IN) :: cr_a0
REAL   , DIMENSION(3)  , INTENT(IN) :: cr_win
INTEGER                , INTENT(IN) :: VAL_PDF
INTEGER                , INTENT(IN) :: VAL_3DPDF
REAL(KIND=PREC_DP)     , INTENT(IN)  :: valmax
INTEGER                , INTENT(OUT) :: ier_num
INTEGER                , INTENT(OUT) :: ier_typ
INTEGER                , INTENT(IN) :: ER_IO
INTEGER                , INTENT(IN) :: ER_APPL
!
INTEGER(KIND=HSIZE_T), PARAMETER       :: dim0 = 1        ! Dimension of "Yell 1.0"
INTEGER(KIND=HSIZE_T), PARAMETER       :: sdim = 8        ! String length of "Yell 1.0"
INTEGER(KIND= SIZE_T), PARAMETER       :: sdimd= 8        ! String length of "DISCUS60"
INTEGER(KIND=2         ), PARAMETER    :: SHORT_NULL = 0  ! a 32 bit 0 for "is_direct"
INTEGER(KIND=2         ), PARAMETER    :: SHORT_ONE  = 1  ! a 32 bit 1 for "is_direct"
!
CHARACTER(LEN=1024) :: line
CHARACTER(LEN=4), PARAMETER :: dataset = "data"           ! Dummy name for HDF5
CHARACTER(LEN=sdim), DIMENSION(1:dim0), TARGET ::  wdata = (/"Yell 1.0"/) ! Write buffer
LOGICAL                                :: isda            ! File foud yes/no
INTEGER                                :: i,j,k,l         ! Dummy indices
INTEGER                                :: hdferr          ! Error returned by HDF5
INTEGER, TARGET                        :: is_direct       ! Das are 3DPDF or diffraction pattern
!
INTEGER(KIND=HSIZE_T), DIMENSION(1)    :: str_dims = (/dim0 /)      ! "Dimensions of the string "DISCSU60" 
INTEGER(KIND=HSIZE_T), DIMENSION(2)    :: data_dims = (/sdim,dim0/) ! "Dimensions of the string "Yell 1.0" 
INTEGER(KIND=SIZE_T),  DIMENSION(1)    :: str_len = (/8/)
INTEGER(KIND=HID_T)                    :: dset            ! dataset Indicator to HDF5
INTEGER(KIND=HID_T)                    :: file_id         ! file Indicator to HDF5
INTEGER(KIND=HID_T)                    :: filetype        ! file Indicator to HDF5
INTEGER(KIND=HID_T)                    ::  memtype        ! memory Indicator to HDF5
INTEGER(KIND=HID_T)                    :: space           ! Communication between calls to HDF5
INTEGER(KIND=HSIZE_T), DIMENSION(1:3)  :: hdims           ! Dimensions in sequence 3,2,1
INTEGER(KIND=HSIZE_T), DIMENSION(1:3)  :: dim_one         ! Dimensions for Yell stuff
REAL(KIND=PREC_DP)                     :: max_data        ! Max data value for normalization
REAL(KIND=PREC_DP), DIMENSION(1:3), TARGET     :: llim    ! Lower Left Corner
REAL(KIND=PREC_DP), DIMENSION(1:3), TARGET     :: steps   ! Step size in H,K,L
REAL(KIND=PREC_DP), DIMENSION(1:6), TARGET     :: uc_hdf  ! Unit cell dimensions
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: values ! Actual diffraction data

TYPE(C_PTR)                            :: f_ptr           ! C-type pointer to data array               
!
INTERFACE
  REAL FUNCTION qval(i, value, ix, iy, laver)
     INTEGER, INTENT(IN) :: i
     INTEGER, INTENT(IN) :: value
     INTEGER, INTENT(IN) :: ix
     INTEGER, INTENT(IN) :: iy
     LOGICAL, INTENT(IN) :: laver
  END FUNCTION qval
END INTERFACE
!
CALL H5open_f(hdferr)                                     ! Open access to HDF5 stream
IF(hdferr/=0) THEN
   ier_num = -2
   ier_typ = ER_IO
   RETURN
ENDIF
!
INQUIRE(FILE=outfile, EXIST=isda)                         ! If file exists, remove
IF(ISDA) THEN
   line = 'rm -f '//outfile(1:LEN_TRIM(outfile))
   CALL EXECUTE_COMMAND_LINE(line, WAIT=.TRUE.)
ENDIF
!
CALL H5Fcreate_f(outfile, H5F_ACC_TRUNC_f, file_id, hdferr)    ! Create output file
IF(hdferr/=0) THEN
   ier_num = -2
   ier_typ = ER_IO
   RETURN
ENDIF
!
is_direct = SHORT_NULL
hdims(3) = out_inc(1)                                        ! Transpose dimensions
hdims(2) = out_inc(2)
hdims(1) = out_inc(3)
!
CALL H5Screate_simple_f(3, hdims, space, hdferr)            ! Create data set
!
ALLOCATE(values(hdims(1), hdims(2), hdims(3)))
!
l = 0                                                       ! Copy proper "value"
DO i = 1, out_inc(1)
   DO j = 1, out_inc(2)
      DO k = 1, out_inc(3)
         l = l + 1
         values(k,j,i) = REAL(qval(l, value, i, j, laver), KIND=PREC_DP)
      ENDDO
   ENDDO
ENDDO
!
IF(valmax>0.0D0) THEN
   max_data = MAXVAL(values)
   IF(max_data >= 0.0) THEN
      values = valmax*values/max_data                          ! Normalize data for HDF5
   ENDIF
ENDIF
!
CALL H5Dcreate_f(file_id, dataset, H5T_NATIVE_DOUBLE, space, dset, hdferr)
f_ptr = C_LOC(values(1,1,1))
CALL H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
!
CALL H5Dclose_f(dset , hdferr)
CALL H5Sclose_f(space, hdferr)
!
!Write other stuff for yell format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the format dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
wdata = 'Yell 1.0'
CALL H5Tcopy_f(H5T_STRING, filetype, hdferr)
CALL H5Tset_strpad_f(filetype, H5T_STR_NULLTERM_F, hdferr)
CALL H5Screate_f(H5S_SCALAR_F, space, hdferr)
CALL H5Dcreate_f(file_id, 'format', filetype, space, dset, hdferr)
CALL H5Dwrite_vl_f(dset, filetype, wdata, data_dims, str_len, hdferr, space)
CALL H5Dclose_f(dset , hdferr)
CALL H5Sclose_f(space, hdferr)
CALL H5Tclose_f(filetype, hdferr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the is_direct dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IF(value==val_3DPDF .OR. value==val_PDF) THEN
   is_direct = SHORT_ONE
ELSE
   is_direct = SHORT_NULL
ENDIF
f_ptr = C_LOC(is_direct)
CALL H5Screate_f(H5S_SCALAR_F, space, hdferr)
CALL H5Dcreate_f(file_id, 'is_direct', H5T_STD_I32LE, space, dset, hdferr)
CALL H5Dwrite_f(dset, H5T_STD_I32LE, f_ptr         , hdferr)
CALL H5Dclose_f(dset , hdferr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the lower_limits dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dim_one = 3
llim(1) = out_eck(1,1)                                      ! Place lower left corner
llim(2) = out_eck(2,1)
llim(3) = out_eck(3,1)
f_ptr = C_LOC(llim(1))
CALL H5Screate_simple_f(1, dim_one, space, hdferr)
CALL H5Dcreate_f(file_id, 'lower_limits',H5T_IEEE_F64BE, space, dset, hdferr)
CALL H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr         , hdferr)
CALL H5Dclose_f(dset , hdferr)
CALL H5Sclose_f(space, hdferr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the step_sizes dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dim_one = 3
steps(1) = out_vi(1,1)
steps(2) = out_vi(2,2)
steps(3) = out_vi(3,3)
f_ptr = C_LOC(steps(1))
CALL H5Screate_simple_f(1, dim_one, space, hdferr)
CALL H5Dcreate_f(file_id, 'step_sizes',H5T_IEEE_F64BE, space, dset, hdferr)
CALL H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr         , hdferr)
CALL H5Dclose_f(dset , hdferr)
CALL H5Sclose_f(space, hdferr)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the unit_cell dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
dim_one = 6
uc_hdf(1:3) = cr_a0(1:3)
uc_hdf(4:6) = cr_win(1:3)
f_ptr = C_LOC(uc_hdf(1))
CALL H5Screate_simple_f(1, dim_one, space, hdferr)
CALL H5Dcreate_f(file_id, 'unit_cell',H5T_IEEE_F64BE, space, dset, hdferr)
CALL H5Dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr         , hdferr)
CALL H5Dclose_f(dset , hdferr)
CALL H5Sclose_f(space, hdferr)
!
str_dims = 1
wdata = 'DISCUS60'
!dataset = 'PROGRAM'
CALL H5Tcopy_f(H5T_C_S1, filetype, hdferr)           ! Copy C-strign into 'filetype'
CALL H5Tset_size_f(filetype, sdimd+1, hdferr)        ! +1 include C nullstring
CALL H5Tcopy_f( H5T_FORTRAN_S1, memtype, hdferr)     ! Copy F-string into 'memtype'
CALL H5Tset_size_f(memtype, sdimd, hdferr)           ! Required size for Fortran string
CALL H5Screate_simple_f(1, str_dims, space, hdferr)     ! Create data space 
CALL H5Dcreate_f(file_id, 'PROGRAM', filetype, space, dset, hdferr) ! Create dataset; write C-string info
f_ptr = C_LOC(wdata(1)(1:1))
CALL H5Dwrite_f(dset, memtype, f_ptr, hdferr)        ! Now write Fortran string
CALL H5Dclose_f(dset , hdferr)                       ! Close data set
CALL H5Sclose_f(space, hdferr)                       ! Close data set
CALL H5Tclose_f(filetype, hdferr)                    ! Close C-string space
CALL H5Tclose_f(memtype, hdferr)                     ! Close F-string space
!
CALL H5Fclose_f(file_id , hdferr)
!
CALL h5close_f(hdferr)
DEALLOCATE(values)
!
END SUBROUTINE hdf5_write
!
!*****7*****************************************************************
!
!END MODULE hdf_write_mod
