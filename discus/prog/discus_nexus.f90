MODULE nexus_discus
!
USE NXUmodule
USE errlist_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC nexus_write
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
INTERFACE nexus_write_data          ! Define a generic name 
   MODULE PROCEDURE nexus_write_1D, nexus_write_2D, nexus_write_3D
END INTERFACE nexus_write_data
!
CONTAINS
!
! Write output in NeXus format
!
   SUBROUTINE nexus_write ( value, laver)
!
   USE output_mod 
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: value
   LOGICAL, INTENT(IN) :: laver
!
   TYPE(NXhandle)  :: fileId
   INTEGER   :: nx_status
!
   write(*,*) ' THIS IS NEXUS WRITE'
   CALL no_error
!
   nx_status = NXopen(outfile, NXACC_CREATE5, fileId)
   IF(nx_status == NX_OK) THEN
      IF( out_inc(1) > 1 .AND. out_inc(2) == 1 .AND. out_inc(3) == 1 ) THEN
         CALL nexus_write_data( value, laver, fileId, nx_status, out_inc(1))
      ELSEIF( out_inc(1) > 1 .AND. out_inc(2) >  1 .AND. out_inc(3) == 1 ) THEN
         CALL nexus_write_data( value, laver, fileId, nx_status, out_inc(1), out_inc(2))
      ELSE
         CALL nexus_write_data( value, laver, fileId, nx_status, out_inc(1), out_inc(2), out_inc(3))
      ENDIF
      nx_status = NXclose(fileId)
      IF(nx_status /= NX_OK) THEN
         write(*,*) ' Could not close NeXus File', nx_status
      ENDIF
   ELSE
      write(*,*) ' Could not open NeXus File', nx_status
   ENDIF
!
   END SUBROUTINE nexus_write
!
!
   SUBROUTINE nexus_write_1D(value, laver, fileId, nx_status, dimx)
!
   USE crystal_mod 
   USE diffuse_mod 
   USE fourier_sup
   USE output_mod 
   USE qval_mod
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)              :: value
   LOGICAL, INTENT(IN)              :: laver
   TYPE(NXhandle)   , INTENT(INOUT) :: fileId
   INTEGER          , INTENT(OUT)   :: nx_status
   INTEGER          , INTENT(IN)    :: dimx
!  
   INTEGER                             :: status_al
   INTEGER                             :: i,j
   CHARACTER(LEN=80)                   :: title = 'Default DISCUS title'
   REAL, DIMENSION(:)    , ALLOCATABLE :: sq
   REAL, DIMENSION(:)    , ALLOCATABLE :: qabs_h
!
!  REAL :: qval
!
   ALLOCATE ( sq(dimx),     STAT=status_al)
   ALLOCATE ( qabs_h(dimx), STAT=status_al)
   sq     = 0.0
   qabs_h = 0.0
!
   DO i=1,out_inc(1)
      qabs_h(i) = out_eck(extr_abs,1) + (i-1)*out_vi(extr_abs,1)
   ENDDO
!
      DO i = 1, out_inc(1)
         sq(i) =  qval ( i, value,  i, j, laver)
      ENDDO
!
!  Prepare Character strings for axes names
   WRITE(caxes, 2000) chkl(out_extr_abs), chkl(out_extr_ord)
   WRITE(cabs_h,2100) chkl(out_extr_abs)
write(*,*) 'CAXES ', caxes,' ', cabs_h,' ', cvalue(value)
!
   nx_status = NXUwritegroup ( fileId, "entry", "NXentry")
   nx_status = NXUwritegroup ( fileId, "data",  "NXdata")
   nx_status = NXUwritedata  ( fileId, "title", title)
   nx_status = NXUwritedata  ( fileId, "Sq", Sq)                ! "Sq" is hopefully flexible...
   nx_status = NXputattr     ( fileId, "signal", 1)
   nx_status = NXputattr     ( fileId, "axes", caxes(1:2) )     ! Usually "Qh")
   nx_status = NXputattr     ( fileId, "long_name", cvalue(value)  ) ! Name of output field "I(hkl)" etc
   nx_status = NXUwritedata  ( fileId, cabs_h, qabs_h, "rlu")   ! Write value of abszissa usually h
   nx_status = NXclosegroup  ( fileId)
   nx_status = NXclosegroup  ( fileId)
!
2000 FORMAT('Q',a1,':Q',a1)
2100 FORMAT('Q',a1)
!
   END SUBROUTINE nexus_write_1D
!
!
   SUBROUTINE nexus_write_2D(value, laver, fileId, nx_status, dimx, dimy)
!
   USE crystal_mod 
   USE diffuse_mod 
   USE fourier_sup
   USE output_mod 
   USE qval_mod
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)              :: value
   LOGICAL, INTENT(IN)              :: laver
   TYPE(NXhandle)   , INTENT(INOUT) :: fileId
   INTEGER          , INTENT(OUT)   :: nx_status
   INTEGER          , INTENT(IN)    :: dimx
   INTEGER          , INTENT(IN)    :: dimy
!  
   INTEGER                             :: status_al
   INTEGER                             :: i,j
   CHARACTER(LEN=80)                   :: title = 'Default DISCUS title'
   REAL, DIMENSION(:,:  ), ALLOCATABLE :: sq
   REAL, DIMENSION(:)    , ALLOCATABLE :: qabs_h
   REAL, DIMENSION(:)    , ALLOCATABLE :: qord_k
!
!  REAL :: qval
!
   ALLOCATE ( sq(dimy, dimx), STAT=status_al)
   ALLOCATE ( qabs_h(dimx), STAT=status_al)
   ALLOCATE ( qord_k(dimy), STAT=status_al)
   sq     = 0.0
   qabs_h = 0.0
   qord_k = 0.0
!
   DO i=1,out_inc(1)
      qabs_h(i) = out_eck(extr_abs,1) + (i-1)*out_vi(extr_abs,1)
   ENDDO
!
   DO i=1,out_inc(2)
      qord_k(i) = out_eck(extr_ord,1) + (i-1)*out_vi(extr_ord,2)
   ENDDO
!
      DO i = 1, out_inc(1)
         DO j = 1, out_inc(2)
            sq(j,i) =  qval ( (i - 1) * out_inc (2) + j,     &
                              value,  i, j, laver)
         ENDDO
      ENDDO
!
!  Prepare Character strings for axes names
   WRITE(caxes, 2000) chkl(out_extr_abs), chkl(out_extr_ord)
   WRITE(cabs_h,2100) chkl(out_extr_abs)
   WRITE(cord_k,2100) chkl(out_extr_ord)
write(*,*) 'CAXES ', caxes,' ', cabs_h,' ', cord_k,' ', cvalue(value)
!
   nx_status = NXUwritegroup ( fileId, "entry", "NXentry")
   nx_status = NXUwritegroup ( fileId, "data",  "NXdata")
   nx_status = NXUwritedata  ( fileId, "title", title)
   nx_status = NXUwritedata  ( fileId, "Sq", Sq)                ! "Sq" is hopefully flexible...
   nx_status = NXputattr     ( fileId, "signal", 1)
   nx_status = NXputattr     ( fileId, "axes", caxes(1:5) )     ! Usually "Qh:Qk")
   nx_status = NXputattr     ( fileId, "long_name", cvalue(value)  ) ! Name of output field "I(hkl)" etc
   nx_status = NXUwritedata  ( fileId, cabs_h, qabs_h, "rlu")   ! Write value of abszissa usually h
   nx_status = NXUwritedata  ( fileId, cord_k, qord_k, "rlu")   ! Write value of ordinate usually k
   nx_status = NXclosegroup  ( fileId)
   nx_status = NXclosegroup  ( fileId)
!
2000 FORMAT('Q',a1,':Q',a1)
2100 FORMAT('Q',a1)
!
   END SUBROUTINE nexus_write_2D
!
!
   SUBROUTINE nexus_write_3D(value, laver, fileId, nx_status, dimx, dimy, dimz)
!
   USE crystal_mod 
   USE diffuse_mod 
   USE fourier_sup
   USE output_mod 
   USE random_mod
   USE qval_mod 
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)              :: value
   LOGICAL, INTENT(IN)              :: laver
   TYPE(NXhandle)   , INTENT(INOUT) :: fileId
   INTEGER          , INTENT(OUT)   :: nx_status
   INTEGER          , INTENT(IN)    :: dimx
   INTEGER          , INTENT(IN)    :: dimy
   INTEGER          , INTENT(IN)    :: dimz
!  
   INTEGER                             :: status_al
   INTEGER                             :: i,j,l
   CHARACTER(LEN=80)                   :: title = 'Default DISCUS title'
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: sq
   REAL, DIMENSION(:)    , ALLOCATABLE :: qabs_h
   REAL, DIMENSION(:)    , ALLOCATABLE :: qord_k
   REAL, DIMENSION(:)    , ALLOCATABLE :: qtop_l
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
   nx_status = NXUwritegroup ( fileId, "entry", "NXentry")
   nx_status = NXUwritegroup ( fileId, "data",  "NXdata")
   nx_status = NXUwritedata  ( fileId, "title", title)
   nx_status = NXUwritedata  ( fileId, "Sq", Sq)                ! "Sq" is hopefully flexible...
   nx_status = NXputattr     ( fileId, "signal", 1)
   nx_status = NXputattr     ( fileId, "axes", caxes )          ! Usually "Qh:Qk:Ql")
   nx_status = NXputattr     ( fileId, "long_name", cvalue(value)  ) ! Name of output field "I(hkl)" etc
   nx_status = NXUwritedata  ( fileId, cabs_h, qabs_h, "rlu")   ! Write value of abszissa usually h
   nx_status = NXUwritedata  ( fileId, cord_k, qord_k, "rlu")   ! Write value of ordinate usually k
   nx_status = NXUwritedata  ( fileId, ctop_l, qtop_l, "rlu")   ! Write value of top axis usually l
   nx_status = NXclosegroup  ( fileId)
   nx_status = NXclosegroup  ( fileId)
!
2000 FORMAT('Q',a1,':Q',a1,':Q',a1)
2100 FORMAT('Q',a1)
!
   END SUBROUTINE nexus_write_3D
!
END MODULE nexus_discus
