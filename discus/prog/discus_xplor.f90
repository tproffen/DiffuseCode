MODULE discus_xplor
!
USE errlist_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC xplor_write
!
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE xplor_write( value, laver)
!
USE crystal_mod
USE diffuse_mod
USE output_mod
USE qval_mod
!
USE envir_mod
USE errlist_mod
USE lib_errlist_func
USE lib_length
USE precision_mod
USE string_convert_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: value
LOGICAL, INTENT(IN) :: laver
!
INTEGER, PARAMETER :: IMRC = 97
!
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=PREC_STRING) :: message
INTEGER            :: i, j,l, ios
INTEGER            :: l_datei
INTEGER, DIMENSION(3)   :: start      ! start point in pixels
INTEGER, DIMENSION(3)   :: fini       ! End point in pixels
INTEGER, DIMENSION(3,2) :: mp         ! End point in pixels
!
CALL no_error
!
l_datei = LEN_STR(outfile)
string  = outfile(l_datei-5:l_datei)
CALL do_cap(string)
IF(string /= '.XPLOR')  THEN
  outfile = outfile(1:l_datei)//'.xplor'
ENDIF
!
IF (outfile (1:1) .eq.'~') THEN 
   line = ' '
   line = home_dir (1:home_dir_l) //outfile (2:l_datei)
   outfile = line(1:200)
ENDIF
! 
OPEN(UNIT=IMRC,FILE=outfile,STATUS='unknown', &
     FORM='formatted', ACCESS='sequential',   &
     IOSTAT=ios,IOMSG=message)
IF(ios/=0) THEN
   ier_num = -2
   ier_typ = ER_IO
   ier_msg(3) = message(1:80)
   RETURN
ENDIF
!
start = 1                       ! Default is full field
fini  = out_inc
IF(out_lrange/=0) THEN          ! User limited output range
   IF(out_lcenter == 0) THEN    ! Center at midpoint
      out_center(1) = MAX(1,NINT(out_inc(1)/2.))
      out_center(2) = MAX(1,NINT(out_inc(2)/2.))
      out_center(3) = MAX(1,NINT(out_inc(3)/2.))
   ENDIF
   IF(out_lpixel==0) THEN
      out_pixel(1) = out_inc(1)
      out_pixel(2) = out_inc(2)
      out_pixel(3) = out_inc(3)
   ENDIF
   IF(out_lrange==2) THEN
      DO i=1,3
         IF(out_quad(i:i) == 'l') THEN
            out_pixel (i) = (out_inc(i)+MOD(out_inc(i),2))/2
            out_center(i) = (out_inc(i)+MOD(out_inc(i),2))/4 + 1
         ELSEIF(out_quad(i:i) == 'r') THEN
            out_pixel (i) = (out_inc(i)+MOD(out_inc(i),2))/2
            out_center(i) = (out_inc(i)+MOD(out_inc(i),2))/4 + out_pixel(i) + 1 - MOD(out_inc(i),2)
         ENDIF
      ENDDO
   ENDIF
   mp(:,1) = out_pixel/2
   mp(:,2) = out_pixel/2
   IF(MOD(out_pixel(1),2)==0) mp(1,2) = mp(1,2) -1
   IF(MOD(out_pixel(2),2)==0) mp(2,2) = mp(2,2) -1
   IF(MOD(out_pixel(3),2)==0) mp(3,2) = mp(3,2) -1
!
   start(1) = MIN(MAX(1,out_center(1) - mp(1,1)), out_inc(1))
   start(2) = MIN(MAX(1,out_center(2) - mp(2,1)), out_inc(2))
   start(3) = MIN(MAX(1,out_center(3) - mp(3,1)), out_inc(3))
   fini (1) = MIN(MAX(1,out_center(1) + mp(1,2)), out_inc(1))
   fini (2) = MIN(MAX(1,out_center(2) + mp(2,2)), out_inc(2))
   fini (3) = MIN(MAX(1,out_center(3) + mp(3,2)), out_inc(3))
ENDIF
!
!  Write header is identical to all 
!
WRITE(IMRC,*)                    ! Empty header line
line = 'TITLE ' //cr_name(1:LEN_TRIM(cr_name)) // ' ' // cvalue(value) // ' written by DISCUS '
!
WRITE(IMRC,'(a80)') line
WRITE(IMRC,'(9i8)') fini(1)-start(1)+1, 1, fini(1)-start(1)+1,    &
                    fini(2)-start(2)+1, 1, fini(2)-start(2)+1,    &
                    fini(3)-start(3)+1, 1, fini(3)-start(3)+1
If(value==val_3DPDF) THEN
   WRITE(IMRC,'(6e12.5)') cr_a0(1:3), cr_win(1:3)
ELSEIF(value/=val_pdf) THEN
   WRITE(IMRC,'(6e12.5)') cr_ar(1:3), cr_wrez(1:3)
ENDIF
WRITE(IMRC,'(a5)') '  ZYX'
DO l=start(3), fini(3)
   WRITE(IMRC,'(i8)') l-start(3) + 1
   WRITE(IMRC,'(6e12.5)') ((qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
                                   (j - 1) * out_inc(3)             + l,     &
                                    value,  i, j, laver),                    &
                            i=start(1), fini(1)), j = start(2), fini(2))
ENDDO
WRITE(IMRC,'(i8)') -9999
WRITE(IMRC,'(2(e12.4, 1x))') 1.00D0, 1.00D0
CLOSE(IMRC)
!
CLOSE(IMRC)
!
END SUBROUTINE xplor_write
!
!*******************************************************************************
!
END MODULE discus_xplor
