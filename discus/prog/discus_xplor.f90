MODULE discus_xplor
!
USE errlist_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC xplor_write
public grd_write
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
!   WRITE(IMRC,'(6e12.5)') ((qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                                   (j - 1) * out_inc(3)             + l,     &
!                                    value,  i, j, laver),                    &
!                            i=start(1), fini(1)), j = start(2), fini(2))
   WRITE(IMRC,'(6e12.5)') ((qval(i,j,l, value,  i, j, laver),                &
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
subroutine grd_write( value, laver, cpatt, spatt, dpatt)
!-
!  Writes volumetric data in ASCII "grd" format
!
!  Line 1: title 
!  line 2: lattice parameters (direct / reciprocal)
!  line 3: NPx  NPy  NPz
!  Lines : (((values(ix, iy, iz), iz=1,NP3), iy=1,NP2), ix=1,NP1)
!
use crystal_mod
use diffuse_mod
use output_mod
use qval_mod
!
use envir_mod
use errlist_mod
use lib_errlist_func
use lib_length
use matrix_mod
use precision_mod
use string_convert_mod
!
implicit none
!
integer, intent(in) :: value
logical, intent(in) :: laver
character(len=*), intent(in) :: cpatt        ! Patterson overlay = none, unit, full
character(len=*), intent(in) :: spatt        ! Selected atoms
character(len=*), intent(in) :: dpatt        ! Deselected atoms
!
integer, parameter :: IMRC = 97
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: file_base
character(len=PREC_STRING) :: file_ext
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: message
integer            :: i, j,l, ios
integer            :: l_datei
integer, dimension(3)   :: start      ! start point in pixels
integer, dimension(3)   :: fini       ! End point in pixels
integer, dimension(3)   :: npx        ! number of    pixels
integer, dimension(3,2) :: mp         ! End point in pixels
logical                 :: lreal          ! Real space / reciprocal T/F
logical                 :: lpatt          ! Write patterson overlay T/F
logical                 :: lpatt_full     ! Write patterson overlay T/F
integer      , dimension(2) :: isurf       ! Values for surfaces (+, -)
real(kind=PREC_DP), dimension(2) :: rsurf       ! Values for surfaces (+, -)
real(kind=PREC_DP), dimension(2) :: extrema     ! Values for surfaces (+, -)
real(kind=PREC_DP), dimension(:), allocatable :: qvals
real(kind=PREC_DP), dimension(3,3)            :: gmat
real(kind=PREC_DP), dimension(3,3)            :: rmat
!
call no_error
!
lreal = .false.                       ! Assume reciprocal space
lpatt = .false.                       ! Assume no patterson overlay
lpatt_full = .false.                  ! Assume no patterson overlay over the full image
if(cpatt=='none') then
   lpatt = .false.                       ! Assume no patterson overlay
   lpatt_full = .false.                  ! Assume no patterson overlay over the full image
else if(cpatt=='unit') then
   lpatt = .true.                        ! Use  patterson overlay
   lpatt_full = .false.                  ! No patterson overlay over the full image
else if(cpatt=='full') then
   lpatt = .true.                        ! Use patterson overlay
   lpatt_full = .true.                   ! Use patterson overlay over the full image
end if
   
!
l_datei = len_trim(outfile)
string  = outfile(l_datei-5:l_datei)
call do_cap(string)
file_base = ' '
if(string == '.VESTA')  then
  file_base(1:l_datei-6) = outfile(1:l_datei-6)
else if(string(3:6) == '.GRD')  then
  file_base(1:l_datei-4) = outfile(1:l_datei-4)
else
  file_base(1:l_datei) = outfile(1:l_datei)
end if
outfile = ' '
outfile(1:len_trim(file_base)+4) = file_base(1:len_trim(file_base)) // '.grd'
file_ext = '.grd'
!
if (outfile (1:1) .eq.'~') then 
   line = ' '
   line = home_dir (1:home_dir_l) //outfile (2:l_datei)
   outfile = line(1:200)
end if
! 
open(unit=IMRC, file=outfile, status='unknown', &
     form='formatted', access='sequential',   &
     iostat=ios, iomsg=message)
if(ios/=0) then
   ier_num = -2
   ier_typ = ER_IO
   ier_msg(3) = message(1:80)
   return
end if
!
start = 1                       ! Default is full field
fini  = out_inc
if(out_lrange/=0) then          ! User limited output range
   if(out_lcenter == 0) then    ! Center at midpoint
      out_center(1) = max(1,nint(out_inc(1)/2.))
      out_center(2) = max(1,nint(out_inc(2)/2.))
      out_center(3) = max(1,nint(out_inc(3)/2.))
   end if
   if(out_lpixel==0) then
      out_pixel(1) = out_inc(1)
      out_pixel(2) = out_inc(2)
      out_pixel(3) = out_inc(3)
   end if
   if(out_lrange==2) then
      do i=1,3
         if(out_quad(i:i) == 'l') then
            out_pixel (i) = (out_inc(i)+mod(out_inc(i),2))/2
            out_center(i) = (out_inc(i)+mod(out_inc(i),2))/4 + 1
         elseif(out_quad(i:i) == 'r') then
            out_pixel (i) = (out_inc(i)+mod(out_inc(i),2))/2
            out_center(i) = (out_inc(i)+mod(out_inc(i),2))/4 + out_pixel(i) + 1 - mod(out_inc(i),2)
         end if
      end do
   end if
   mp(:,1) = out_pixel/2
   mp(:,2) = out_pixel/2
   if(mod(out_pixel(1),2)==0) mp(1,2) = mp(1,2) -1
   if(mod(out_pixel(2),2)==0) mp(2,2) = mp(2,2) -1
   if(mod(out_pixel(3),2)==0) mp(3,2) = mp(3,2) -1
!
   start(1) = min(max(1,out_center(1) - mp(1,1)), out_inc(1))
   start(2) = min(max(1,out_center(2) - mp(2,1)), out_inc(2))
   start(3) = min(max(1,out_center(3) - mp(3,1)), out_inc(3))
   fini (1) = min(max(1,out_center(1) + mp(1,2)), out_inc(1))
   fini (2) = min(max(1,out_center(2) + mp(2,2)), out_inc(2))
   fini (3) = min(max(1,out_center(3) + mp(3,2)), out_inc(3))
end if
!
!  Write header is identical to all 
!
line = 'TITLE ' //cr_name(1:len_trim(cr_name)) // ' ' // cvalue(value) // ' written by DISCUS '
!
write(IMRC,'(a80)') line                            ! Title line
if(value==val_3DPDF) then
   write(IMRC,'(6e12.5)') cr_a0(1:3), cr_win(1:3)   ! Lattice parameters
   lreal = .true.
elseif(value/=val_pdf) then
   write(IMRC,'(6e12.5)') cr_ar(1:3), cr_wrez(1:3)  ! reciprocal lattice parameters
   lreal = .false.
end if
write(IMRC,'(3i8)') fini(1)-start(1)+1, fini(2)-start(2)+1, fini(3)-start(3)+1
!
allocate(qvals(start(3): fini(3)))
extrema = qval(1,1,1, value, 1, 1, laver)
qvals = 0.0
!
do i=start(1), fini(1)
   do j=start(2), fini(2)
      do l=start(3), fini(3)
!        qvals(l) = qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                          (j - 1) * out_inc(3)             + l,     &
!                          value,  i, j, laver)
         qvals(l) = qval(i,j,l, value,  i, j, laver)
      enddo
      extrema(1) = min(extrema(1), minval(qvals))
      extrema(2) = max(extrema(2), maxval(qvals))
      write(IMRC,'(6g16.6e3)') (qvals(l),l=start(3), fini(3))
!                               (qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
!                                       (j - 1) * out_inc(3)             + l,     &
!                                        value,  i, j, laver),                    &
!                               l=start(3), fini(3))
   end do
end do
!do i=start(1)-1, fini(1)-1
!   do j=start(2)-1, fini(2)-1
!      write(IMRC,'(6g16.6e3)') (qvals(l,j,i),l=start(3)-1, fini(3)-1)
!   end do
!end do
deallocate(qvals)
!
close(IMRC)
!
if(lreal) then
   if(fave==0) then
      rsurf(1) =  0.20*abs(extrema(2))     ! 20
      rsurf(2) =  0.70*abs(extrema(2))     ! 80
      isurf(1) = 1
      isurf(2) = 1
   else
      rsurf(1) =  0.20*abs(extrema(1))
      rsurf(2) =  0.20*abs(extrema(1))
      isurf(1) = 2
      isurf(2) = 1
   end if
else
   if(fave==0) then
      rsurf(1) = max(0.0D0, extrema(1)) + (extrema(2)-extrema(1))*0.0025D0
      rsurf(2) = max(0.0D0, extrema(1)) + (extrema(2)-extrema(1))*0.01D0
   else
      rsurf(1) = max(0.0D0, extrema(1)) + (extrema(2)-extrema(1))*0.333D0
      rsurf(2) = max(0.0D0, extrema(1)) + (extrema(2)-extrema(1))*0.500D0
   end if
      isurf(1) = 1
      isurf(2) = 1
end if
!
gmat(1,:) = out_vi(1,:)*1.0D0*      (fini(1)-start(1))
gmat(2,:) = out_vi(2,:)*1.0D0*      (fini(2)-start(2))
gmat(3,:) = out_vi(3,:)*1.0D0*      (fini(3)-start(3))
call matinv3(gmat, rmat) 
npx = fini-start
call vesta_write(file_base, file_ext, lreal, lpatt, lpatt_full, &
                 isurf, rsurf, gmat, rmat, npx, spatt, dpatt)
!
end subroutine grd_write
!
!*******************************************************************************
!
subroutine vesta_write(file_base, file_ext, lreal, lpatt, lpatt_full, &
                       isurf, rsurf, gmat, rmat, npx, spatt, dpatt)
!-
!  Write the VEST file for the current volumetric data
!+
!
use crystal_mod
!
use errlist_mod
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: file_base      ! File name base
character(len=*), intent(in) :: file_ext       ! file name extension including '.'
logical         , intent(in) :: lreal          ! Real space / reciprocal T/F
logical         , intent(in) :: lpatt          ! Write patterson overlay T/F
logical         , intent(in) :: lpatt_full     ! Write patterson overlay T/F
integer      , dimension(2), intent(in) :: isurf       ! Values for surfaces (+, -)
real(kind=PREC_DP), dimension(2), intent(in) :: rsurf       ! Values for surfaces (+, -)
real(kind=PREC_DP), dimension(3,3), intent(in)  :: gmat
real(kind=PREC_DP), dimension(3,3), intent(in)  :: rmat
integer      , dimension(3), intent(in) :: npx         ! Number of pixels - 1
character(len=*), intent(in) :: spatt        ! Selected atoms
character(len=*), intent(in) :: dpatt        ! Deselected atoms
!
integer, parameter :: IMRC = 97
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: outfile
character(len=PREC_STRING) :: message
integer :: ios                               ! I/O status
integer :: isep                              ! Directory separator

!
outfile = file_base(1:len_trim(file_base)) // '.vesta'
!
open(unit=IMRC, file=outfile, status='unknown', &
     form='formatted', access='sequential',   &
     iostat=ios, iomsg=message)
if(ios/=0) then
   ier_num = -2
   ier_typ = ER_IO
   ier_msg(3) = message(1:80)
   return
end if
isep = max(index(file_base(1:len_trim(file_base)),'/'), &
           index(file_base(1:len_trim(file_base)),'\'))
if(isep>0) then
  line = file_base(isep+1:len_trim(file_base)) // file_ext (1:len_trim(file_ext))
else
  line = file_base(     1:len_trim(file_base)) // file_ext (1:len_trim(file_ext))
endif
!
write(IMRC,'(a)') '#VESTA_FORMAT_VERSION 3.5.0'
write(IMRC,'(a)')
write(IMRC,'(a)')
write(IMRC,'(a)') 'CRYSTAL'
write(IMRC,'(a)')
write(IMRC,'(a)') 'TITLE'
write(IMRC,'(a)') cr_name(1:len_trim(cr_name))
write(IMRC,'(a)')
write(IMRC,'(a)') 'IMPORT_DENSITY 1'
write(IMRC,'(a,a,a)') '+1.000000E+00 ',line(1:len_trim(line))
write(IMRC,'(a)')
write(IMRC,'(a)') 'GROUP'
write(IMRC,'(a)') '1 1 P 1'
write(IMRC,'(a)') 'SYMOP'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1'
write(IMRC,'(a)') ' -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0'
write(IMRC,'(a)') 'TRANM 0'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1'
write(IMRC,'(a)') 'LTRANSL'
write(IMRC,'(a)') ' -1'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') 'LORIENT'
write(IMRC,'(a)') ' -1   0   0   0   0'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') 'LMATRIX'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000'
write(IMRC,'(a)') 'CELLP'
if(lreal) then
   write(IMRC,'(5(f10.6,1x),f10.6)') cr_a0, cr_win   ! Direct space lattice parameters
else
   write(IMRC,'(5(f10.6,1x),f10.6)') cr_ar, cr_wrez  ! reciprocal lattice parameters
endif
write(IMRC,'(a)') '  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000'
write(IMRC,'(a)') 'STRUC'
!
if(lpatt) then
   call vesta_real(IMRC, gmat, rmat, npx, lpatt_full, spatt, dpatt)
else
   write(IMRC,'(a)') '  0 0 0 0 0 0 0'
   write(IMRC,'(a)') 'THERI 1'
   write(IMRC,'(a)') '  0 0 0'
   write(IMRC,'(a)') 'SHAPE'
   write(IMRC,'(a)') '  0       0       0       0   0.000000  0   192   192   192   192'
   write(IMRC,'(a)') 'BOUND'
   write(IMRC,'(a)') '       0        1         0        1         0        1'
   write(IMRC,'(a)') '  0   0   0   0  0'
   write(IMRC,'(a)') 'SBOND'
   write(IMRC,'(a)') '  0 0 0 0'
   write(IMRC,'(a)') 'SITET'
   write(IMRC,'(a)') '  0 0 0 0 0 0'
   write(IMRC,'(a)') 'VECTR'
   write(IMRC,'(a)') ' 0 0 0 0 0'
   write(IMRC,'(a)') 'VECTT'
   write(IMRC,'(a)') ' 0 0 0 0 0'
   write(IMRC,'(a)') 'SPLAN'
   write(IMRC,'(a)') '  0   0   0   0'
   write(IMRC,'(a)') 'LBLAT'
   write(IMRC,'(a)') ' -1'
   write(IMRC,'(a)') 'LBLSP'
   write(IMRC,'(a)') ' -1'
   write(IMRC,'(a)') 'DLATM'
   write(IMRC,'(a)') ' -1'
   write(IMRC,'(a)') 'DLBND'
   write(IMRC,'(a)') ' -1'
   write(IMRC,'(a)') 'DLPLY'
   write(IMRC,'(a)') ' -1'
   write(IMRC,'(a)') 'PLN2D'
   write(IMRC,'(a)') '  0   0   0   0'
   write(IMRC,'(a)') 'ATOMT'
   write(IMRC,'(a)') '  0 0 0 0 0 0'
endif
!
write(IMRC,'(a)') 'SCENE'
write(IMRC,'(a)') ' 0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') '  0.000   0.000'
write(IMRC,'(a)') '  0.000'
write(IMRC,'(a)') '  1.000'
write(IMRC,'(a)') 'HBOND 0 2'
write(IMRC,'(a)') ''
write(IMRC,'(a)') 'STYLE'
write(IMRC,'(a)') 'DISPF 37753794'
write(IMRC,'(a)') 'MODEL   0  1  0'
write(IMRC,'(a)') 'SURFS   0  1  1'
write(IMRC,'(a)') 'SECTS  32  1'
write(IMRC,'(a)') 'FORMS   0  1'
write(IMRC,'(a)') 'ATOMS   0  0  1'
write(IMRC,'(a)') 'BONDS   1'
write(IMRC,'(a)') 'POLYS   1'
write(IMRC,'(a)') 'VECTS 1.000000'
write(IMRC,'(a)') 'FORMP'
write(IMRC,'(a)') '  1  1.0   0   0   0'
write(IMRC,'(a)') 'ATOMP'
write(IMRC,'(a)') ' 24  24   0  50  2.0   0'
write(IMRC,'(a)') 'BONDP'
write(IMRC,'(a)') '  1  16  0.250  2.000 127 127 127'
write(IMRC,'(a)') 'POLYP'
write(IMRC,'(a)') ' 204 1  1.000 180 180 180'
write(IMRC,'(a)') 'ISURF'
if(lreal) then
   if(isurf(1) == 2) then
      write(IMRC,'(a,i1,G16.6E3,a)') '  1   ',isurf(1),abs(rsurf(1)),'  92 178 253 127 127'
   else
      write(IMRC,'(a,i1,G16.6E3,a)') '  1   ',isurf(1),abs(rsurf(1)),' 155 155   0  60  30'
   endif
   write(IMRC,'(a,i1,G16.6E3,a)') '  2   ',isurf(2),abs(rsurf(2)),' 253  95  89 127 127'
else
   write(IMRC,'(a,i1, G16.6E3,a)') '  1   ',isurf(1),abs(rsurf(1)),' 255 255   0  60  30'
   write(IMRC,'(a,i1, G16.6E3,a)') '  2   ',isurf(2),abs(rsurf(2)), ' 255  40  40 200  30'
!  write(IMRC,'(a,G16.6E3,a)') '  1   1     172813 255 255   0 127 127'
endif
write(IMRC,'(a)') '  0   0   0   0'
write(IMRC,'(a)') 'TEX3P'
write(IMRC,'(a)') '  1  0.00000E+00  1.00000E+00'
write(IMRC,'(a)') 'SECTP'
write(IMRC,'(a)') '  1  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00'
write(IMRC,'(a)') 'CONTR'
write(IMRC,'(a)') ' 0.1 -1 1 1 10 -1 2 5'
write(IMRC,'(a)') ' 2 1 2 1'
write(IMRC,'(a)') '   0   0   0'
write(IMRC,'(a)') '   0   0   0'
write(IMRC,'(a)') '   0   0   0'
write(IMRC,'(a)') '   0   0   0'
write(IMRC,'(a)') 'HKLPP'
write(IMRC,'(a)') ' 192 1  1.000 255   0 255'
write(IMRC,'(a)') 'UCOLP'
write(IMRC,'(a)') '   0   1  1.000   0   0   0'
write(IMRC,'(a)') 'COMPS 1'
write(IMRC,'(a)') 'LABEL 1    12  1.000 0'
write(IMRC,'(a)') 'PROJT 0  0.962'
write(IMRC,'(a)') 'BKGRC'
write(IMRC,'(a)') ' 255 255 255'
write(IMRC,'(a)') 'DPTHQ 1 -0.5000  3.5000'
write(IMRC,'(a)') 'LIGHT0 1'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 20.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 -1.000000'
write(IMRC,'(a)') '  26  26  26 255'
write(IMRC,'(a)') ' 179 179 179 255'
write(IMRC,'(a)') ' 255 255 255 255'
write(IMRC,'(a)') 'LIGHT1'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 20.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 -1.000000'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') 'LIGHT2'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 20.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 -1.000000'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') 'LIGHT3'
write(IMRC,'(a)') ' 1.000000  0.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  1.000000  0.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  1.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000  0.000000  1.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 20.000000  0.000000'
write(IMRC,'(a)') ' 0.000000  0.000000 -1.000000'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') '   0   0   0   0'
write(IMRC,'(a)') 'ATOMM'
write(IMRC,'(a)') ' 204 204 204 255'
write(IMRC,'(a)') ' 126.720'
write(IMRC,'(a)') 'BONDM'
write(IMRC,'(a)') ' 255 255 255 255'
write(IMRC,'(a)') ' 128.000'
write(IMRC,'(a)') 'POLYM'
write(IMRC,'(a)') ' 255 255 255 255'
write(IMRC,'(a)') ' 128.000'
write(IMRC,'(a)') 'SURFM'
write(IMRC,'(a)') '   0   0   0 255'
write(IMRC,'(a)') ' 128.000'
write(IMRC,'(a)') 'FORMM'
write(IMRC,'(a)') ' 255 255 255 255'
write(IMRC,'(a)') ' 128.000'
write(IMRC,'(a)') 'HKLPM'
write(IMRC,'(a)') ' 255 255 255 255'
write(IMRC,'(a)') ' 128.000'
!
close(IMRC)
!
end subroutine vesta_write
!
!*******************************************************************************
!
subroutine vesta_real(IMRC, gmat, rmat, npx, lpatt_full, spatt, dpatt)
!-
!   Write the structure part for Vesta
!+
!
use crystal_mod
use chem_aver_mod
use chem_mod
use modify_mod
!
use param_mod
!
implicit none
!
integer, intent(in) :: IMRC
real(kind=PREC_DP), dimension(3,3), intent(in) :: gmat
real(kind=PREC_DP), dimension(3,3), intent(in) :: rmat
integer           , dimension(3)  , intent(in) :: npx
logical                           , intent(in) :: lpatt_full     ! Write patterson overlay T/F
character(len=*), intent(in) :: spatt        ! Selected atoms
character(len=*), intent(in) :: dpatt        ! Deselected atoms
!
real(kind=PREC_DP), parameter :: EPS = 1.0e-6
integer :: npat    ! Number of Patterson maxima
integer :: i
integer :: j
integer :: k
integer :: l
integer :: m
integer :: lp
integer, dimension(2,3) :: ind_pat              !  Indices of full overlay
logical :: lout
logical :: lsite
logical :: pat_sel_atom
logical, dimension(0:MAXSCAT)    :: pat_latom    ! selected atoms
logical, dimension(1:cr_ncatoms) :: pat_lsite    ! selected sites
real(kind=PREC_DP), dimension(:,:), allocatable :: pat_pos
real(kind=PREC_DP), dimension(3)                :: pat_vec
!
ind_pat = 0
lout  = .true.
lsite = .true.
call chem_aver (lout, lsite)
!
pat_latom    = .true.                ! Assume all atoms were selected
pat_latom(0) = .false.
pat_lsite    = .true.                ! Assume all sites were selected
if(spatt=='all') then
   pat_latom    = .true.                ! All atoms were selected
   pat_latom(0) = .false.
else if(spatt=='none') then
   pat_latom = .false.               ! All atoms were deselected
else
   pat_latom = .false.               ! Assume all atoms were deselected
   pat_sel_atom = .true.
   lp = len_trim(spatt)
   call atom_select(spatt, lp, 0, MAXSCAT, pat_latom, pat_lsite, 1, cr_ncatoms, & 
        pat_sel_atom, .true., .true.)
endif
if(dpatt=='all') then
   pat_latom = .false.               ! All atoms were deselected
else if(dpatt=='none') then
   continue
else
   pat_sel_atom = .true.
   lp = len_trim(dpatt)
   call atom_select(dpatt, lp, 0, MAXSCAT, pat_latom, pat_lsite, 1, cr_ncatoms, & 
        pat_sel_atom, .true., .false.)
endif
!
!
npat = 1
allocate(pat_pos(3, 1 + cr_ncatoms*(cr_ncatoms-1)))
pat_pos(:,:) = 0.0                ! Central maximum
outer: do i = 1, cr_ncatoms
  if(pat_latom(i)) then
     inner: do j = i+1, cr_ncatoms
        if(pat_latom(j)) then
           do k = 1, 3
              pat_vec(k) = chem_ave_pos(k,i) - chem_ave_pos(k,j)
              if(pat_vec(k) < -0.5) then
                  pat_vec(k) =  1.0 + pat_vec(k)
              else if(pat_vec(k) > 0.5) then
                  pat_vec(k) = -1.0 + pat_vec(k)
              end if
           end do
           do k = 1, npat
              if(abs(pat_vec(1)-pat_pos(1,k))<EPS .and.  &   
                 abs(pat_vec(2)-pat_pos(2,k))<EPS .and.  &
                 abs(pat_vec(3)-pat_pos(3,k))<EPS      ) cycle inner   ! found identical patterson vector
           end do
           npat = npat + 1
           pat_pos(:,npat) =  pat_vec
           npat = npat + 1
           pat_pos(:,npat) = -pat_vec
        end if
     end do inner
  end if
end do outer
!
!
if(lpatt_full) then                     ! Do full overlay
   do i=-1,1, 2
      do j=-1,1, 2
         do k=-1,1, 2
            pat_vec(1) = 0.5*real(i)
            pat_vec(2) = 0.5*real(j)
            pat_vec(3) = 0.5*real(k)
            pat_vec = matmul(gmat, pat_vec)
            ind_pat(1,:) = min(ind_pat(1,:), nint(pat_vec(:)-0.5))
            ind_pat(2,:) = max(ind_pat(2,:), nint(pat_vec(:)+0.5))
         end do
      end do
   end do
endif
!
m = 0
do i = ind_pat(1,1), ind_pat(2,1)
do j = ind_pat(1,2), ind_pat(2,2)
do k = ind_pat(1,3), ind_pat(2,3)
do l = 1, npat
   pat_vec(1) = pat_pos(1,l) + real(i) !+ 0.5 - 0.5/real(npx(1))
   pat_vec(2) = pat_pos(2,l) + real(j) !+ 0.5 - 0.5/real(npx(2))
   pat_vec(3) = pat_pos(3,l) + real(k) !+ 0.5 - 0.5/real(npx(3))
   pat_vec = matmul(rmat, pat_vec) + 0.5 - 0.5/real(npx)
   if(pat_vec(1)>=0.0 .and. pat_vec(1)<=1.0   .and.   &
      pat_vec(2)>=0.0 .and. pat_vec(2)<=1.0   .and.   &
      pat_vec(3)>=0.0 .and. pat_vec(3)<=1.0           &
     ) then
      m = m + 1
      write(IMRC,'(i5,a,3f11.6,a)') m, ' Pp          PP  1.0000', pat_vec     ,'    1a       1'
      write(IMRC,'(a)') '                            0.000000   0.000000   0.000000  0.00'
   end if
end do
end do
end do
end do
write(IMRC,'(a)') '  0 0 0 0 0 0 0'
write(IMRC,'(a)') 'THERI 1'
do k= 1, m
   write(IMRC,'(i3,a)')  k, '          PP  0.006333'
end do
write(IMRC,'(a)') '  0 0 0'
write(IMRC,'(a)') 'SHAPE'
write(IMRC,'(a)') '  0       0       0       0   0.000000  0   192   192   192   192'
write(IMRC,'(a)') 'BOUND'
write(IMRC,'(a)') '       0        1         0        1         0        1'
write(IMRC,'(a)') '  0   0   0   0  0'
write(IMRC,'(a)') 'SBOND'
write(IMRC,'(a)') '  0 0 0 0'
write(IMRC,'(a)') 'SITET'
k=1
!do k= 1, m
   write(IMRC,'(i3,a)')  k, '          PP  0.0500  50  50 250  50  40  55   0  0'
!end do
write(IMRC,'(a)') '  0 0 0 0 0 0'
write(IMRC,'(a)') 'VECTR'
write(IMRC,'(a)') ' 0 0 0 0 0'
write(IMRC,'(a)') 'VECTT'
write(IMRC,'(a)') ' 0 0 0 0 0'
write(IMRC,'(a)') 'SPLAN'
write(IMRC,'(a)') '  0   0   0   0'
write(IMRC,'(a)') 'LBLAT'
write(IMRC,'(a)') ' -1'
write(IMRC,'(a)') 'LBLSP'
write(IMRC,'(a)') ' -1'
write(IMRC,'(a)') 'DLATM'
write(IMRC,'(a)') ' -1'
write(IMRC,'(a)') 'DLBND'
write(IMRC,'(a)') ' -1'
write(IMRC,'(a)') 'DLPLY'
write(IMRC,'(a)') ' -1'
write(IMRC,'(a)') 'PLN2D'
write(IMRC,'(a)') '  0   0   0   0'
write(IMRC,'(a)') 'ATOMT'
k = 1
!do k= 1, m
   write(IMRC,'(i3,a)')  k, '         Pp  0.0500  50  50 250 183 187 189 204'
!end do
write(IMRC,'(a)') '  0 0 0 0 0 0'
!
deallocate(pat_pos)
!
end subroutine vesta_real
!
!*******************************************************************************
!
END MODULE discus_xplor
