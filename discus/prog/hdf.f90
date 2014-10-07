MODULE hdf_discus
!
USE H5LT
USE errlist_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC hdf_write
!
CHARACTER(LEN=10), DIMENSION(5), PARAMETER :: cvalue =(/        &  ! name for output value
        'I(hkl)    ', 'A(hkl)    ', 'Phase(hkl)', &
        'Real(hkl) ', 'Imag(hkl) '                &
        /)
CHARACTER(LEN=1 ), DIMENSION(3), PARAMETER :: chkl =(/'h','k','l'/) ! Names for axes
CHARACTER(LEN=8) :: caxes                         ! The actual axes variable
CHARACTER(LEN=2) :: cabs_h                        ! Short names 
CHARACTER(LEN=2) :: cord_k
CHARACTER(LEN=2) :: ctop_l
!
!  This interface block will serve to write 1D, 2D, 3D, etc data with one name 
INTERFACE hdf_write_data          ! Define a generic name 
   MODULE PROCEDURE hdf_write_3D !, hdf_write_2D, hdf_write_3D
END INTERFACE hdf_write_data
!
CONTAINS
!
! Write output in Hdf format
!
   SUBROUTINE hdf_write ( value, laver)
!
   USE output_mod 
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: value
   LOGICAL, INTENT(IN) :: laver
!
   INTEGER(HID_T)  :: fileId
   INTEGER   :: error
!
   write(*,*) ' THIS IS HDF WRITE'
   CALL no_error
!
   call H5open_f(error)
   call H5Fcreate_f(outfile, H5F_ACC_TRUNC_F, fileId, error)
   IF(error == 0) THEN
      IF( out_inc(1) > 1 .AND. out_inc(2) == 1 .AND. out_inc(3) == 1 ) THEN
         !CALL hdf_write_data( value, laver, fileId, nx_status, out_inc(1))
      ELSEIF( out_inc(1) > 1 .AND. out_inc(2) >  1 .AND. out_inc(3) == 1 ) THEN
         !CALL hdf_write_data( value, laver, fileId, nx_status, out_inc(1), out_inc(2))
      ELSE
         CALL hdf_write_data( value, laver, fileId, error, out_inc(1), out_inc(2), out_inc(3))
      ENDIF
      call H5Fclose_f(fileId,error)
      call H5close_f(error)
      IF(error /= 0) THEN
         write(*,*) ' Could not close Hdf File', error
      ENDIF
   ELSE
      write(*,*) ' Could not open Hdf File', error
   ENDIF
!
   END SUBROUTINE hdf_write
!
!
   SUBROUTINE hdf_write_3D(value, laver, fileId, error, dimx, dimy, dimz)
!
   USE crystal_mod 
   USE diffuse_mod 
   USE fourier_sup
   USE output_mod 
   USE random_mod
   USE qval_mod 
   USE iso_c_binding
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)              :: value
   LOGICAL, INTENT(IN)              :: laver
   INTEGER(HID_T)   , INTENT(INOUT) :: fileId
   INTEGER(HID_T)                   :: entry0, data
   INTEGER          , INTENT(OUT)   :: error
   INTEGER          , INTENT(IN)    :: dimx
   INTEGER          , INTENT(IN)    :: dimy
   INTEGER          , INTENT(IN)    :: dimz
   INTEGER(HSIZE_T)                 :: dimxx(1)
   INTEGER(HSIZE_T)                 :: dimyy(1)
   INTEGER(HSIZE_T)                 :: dimzz(1)
!  
   INTEGER                             :: status_al
   INTEGER                             :: i,j,l
   CHARACTER(LEN=80)                   :: title = 'Default DISCUS title'
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: sq
   INTEGER                             :: signal(1) = 1
   INTEGER(size_t)                     :: signal_size = 1
   REAL, DIMENSION(:)    , ALLOCATABLE :: qabs_h
   REAL, DIMENSION(:)    , ALLOCATABLE :: qord_k
   REAL, DIMENSION(:)    , ALLOCATABLE :: qtop_l
   INTEGER(HSIZE_T) :: dims(3)
   dims = out_inc
   dimxx = dimx
   dimyy = dimy
   dimzz = dimz
!
!  REAL :: qval
!
   ALLOCATE ( sq(dimz, dimy, dimx), STAT=status_al)
   ALLOCATE ( qabs_h(dimx), STAT=status_al)
   ALLOCATE ( qord_k(dimy), STAT=status_al)
   ALLOCATE ( qtop_l(dimz), STAT=status_al)
   sq     = 0.0
   qabs_h = 0.0
   qord_k = 0.0
   qtop_l = 0.0
!
   DO i=1,out_inc(1)
      qabs_h(i) = out_eck(extr_abs,1) + (i-1)*out_vi(extr_abs,1)
   ENDDO
!
   DO i=1,out_inc(2)
      qord_k(i) = out_eck(extr_ord,1) + (i-1)*out_vi(extr_ord,2)
   ENDDO
!
   DO i=1,out_inc(3)
      qtop_l(i) = out_eck(extr_top,1) + (i-1)*out_vi(extr_top,3)
   ENDDO
!
   DO i = 1, out_inc(1)
      DO j = 1, out_inc(2)
         DO l = 1, out_inc(3)
            sq(l,j,i) =  qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
                                (j - 1) * out_inc(3)             + l,     &
                                value,  i, j, laver)
         ENDDO
      ENDDO
   ENDDO
!
!  Prepare Character strings for axes names
   WRITE(caxes, 2000) chkl(out_extr_abs), chkl(out_extr_ord), chkl(out_extr_top)
   WRITE(cabs_h,2100) chkl(out_extr_abs)
   WRITE(cord_k,2100) chkl(out_extr_ord)
   WRITE(ctop_l,2100) chkl(out_extr_top)
write(*,*) 'CAXES ', caxes,' ', cabs_h,' ', cord_k,' ', ctop_l,' ',cvalue(value)
!
   !NXUwritegroup ( fileId, "entry", "NXentry")
   call h5gcreate_f ( fileId, "entry", entry0, error)
   !NXUwritegroup ( fileId, "data",  "NXdata")
   call h5gcreate_f ( entry0, "data",  data, error)
   !NXUwritedata  ( fileId, "title", title)
   call h5ltmake_dataset_string_f(data,"title",title,error)
   !NXUwritedata  ( fileId, "Sq", Sq)                ! "Sq" is hopefully flexible...
   call h5ltmake_dataset_float_f(data,"Sq",3,dims,sq,error)
   !NXputattr     ( fileId, "signal", 1)
   call h5ltset_attribute_int_f(data,"Sq","signal",signal,signal_size,error)
   !NXputattr     ( fileId, "axes", caxes )          ! Usually "Qh:Qk:Ql")
   call h5ltset_attribute_string_f(data,"Sq","axes",caxes,error)
   !NXputattr     ( fileId, "long_name", cvalue(value)  ) ! Name of output field "I(hkl)" etc
   call h5ltset_attribute_string_f(data,"Sq","long_name",cvalue(value),error)
   !NXUwritedata  ( fileId, cabs_h, qabs_h, "rlu")   ! Write value of abszissa usually h
   call h5ltmake_dataset_float_f(data,cabs_h,1,dimxx,qabs_h,error)
   !NXUwritedata  ( fileId, cord_k, qord_k, "rlu")   ! Write value of ordinate usually k
   call h5ltmake_dataset_float_f(data,cord_k,1,dimyy,qord_k,error)
   !NXUwritedata  ( fileId, ctop_l, qtop_l, "rlu")   ! Write value of top axis usually l
   call h5ltmake_dataset_float_f(data,ctop_l,1,dimzz,qtop_l,error)
   !NXclosegroup  ( fileId)
   call  h5gclose_f ( data, error)
   !NXclosegroup  ( fileId)
   call  h5gclose_f ( entry0, error)
!
2000 FORMAT('Q',a1,':Q',a1,':Q',a1)
2100 FORMAT('Q',a1)
!
   END SUBROUTINE hdf_write_3D
!
END MODULE hdf_discus
