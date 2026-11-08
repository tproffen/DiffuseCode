module lib_detect
!-
!  Procedures to autodetect file formats
!+
!
implicit none
!
integer, parameter :: FILE_HDF5_STRUCTURE   = 1
integer, parameter :: FILE_HDF5_DATA        = 2
integer, parameter :: FILE_HDF5_DISCUS_YELL = 3
integer, parameter :: FILE_HDF5_UNKNOWN     = 4
integer, parameter :: FILE_SUITE_HEADER     = 5  ! DISCUS FILE with detailed header
integer, parameter :: FILE_SPEC             = 6  ! Generic SPEC file
integer, parameter :: FILE_DIFFPY_PDFGETX2  = 7  ! DIFFPY PDFgetx2 file
!
private
public detect_format
!
contains
!
!*******************************************************************************
!
subroutine detect_format(ifil, filename, iformat)
!-
!  Main procedure for autodetect
!+
!
use errlist_mod
use precision_mod
use support_mod,  only: oeffne
!
implicit none
integer         , intent(in)  :: ifil       ! IO number
character(len=*), intent(in)  :: filename
integer         , intent(out) :: iformat    ! Format type number
!
character(len=PREC_STRING) :: line
integer                    :: NMSG
!
NMSG = ubound(ier_msg,1)
iformat = -1   ! Not detected
!
call oeffne(ifil, filename, 'old')
if(ier_num/=0) then
   ier_msg(1) ='Could not open file '
   ier_msg(2) = filename
   return
endif
!
read(ifil,'(a)') line
write(*,'(3a)') '>>',line(1:len_trim(line)),'<<'
!
if(line(2:5)=='HDF') then
   call detect_hdf(ifil, filename, iformat, NMSG, ier_num, ier_msg)
elseif(line(1:8) == '# HEADER') then
   iformat = FILE_SUITE_HEADER
elseif(line(1:9)=='[DEFAULT]') then
   call detect_pdfgetx(ifil, filename, iformat, NMSG, ier_num, ier_msg)
endif
!
write(*,*) ' FILE FORMAT is ', iformat
!
end subroutine detect_format
!
!*******************************************************************************
!
subroutine detect_hdf(ifil, filename, iformat, NMSG, ier_num, ier_msg)
!-
!  Main procedure for autodetect of HDF5 file formats
!+
!
use precision_mod
use support_mod,  only: oeffne
!
!use unified_dump_mod
!
use lib_h5fortran_mod
use lib_unified_chars_mod
use h5fortran, only: hdf5_file
!
implicit none
integer         , intent(in)  :: ifil       ! IO number
character(len=*), intent(in)  :: filename
integer         , intent(out) :: iformat    ! Format type number
integer         , intent(in ) :: NMSG               ! Dimension of ier_msg
integer         , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*), dimension(NMSG), intent(out) :: ier_msg            ! Error messages
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING), dimension(N_META):: crystal_meta
integer                    :: lib_type
integer                    :: i
integer, dimension(-1:4), parameter :: HDF5_FORMAT = &
   (/-1,-1, FILE_HDF5_STRUCTURE, FILE_HDF5_DATA, FILE_HDF5_DISCUS_YELL, FILE_HDF5_UNKNOWN/)
!
type(hdf5_file) :: h5f
!
iformat = -1   ! Not detected
close(unit=ifil)
!
call h5f_open(h5f, filename, 'r', NMSG, ier_num, ier_msg)
if(ier_num/=0) then
   return
endif
!
lib_type = -1
call h5f_check_read_meta(h5f, lib_type, crystal_meta, NMSG, ier_num, ier_msg)
!
iformat = HDF5_FORMAT(lib_type)
!
call h5f%close()
!
end subroutine detect_hdf
!
!*******************************************************************************
!
subroutine detect_pdfgetx(ifil, filename, iformat, NMSG, ier_num, ier_msg)
!-
!  Main procedure for autodetect of DIFFPY.PDFGETX
!+
!
use precision_mod
use support_mod,  only: oeffne
!
!use unified_dump_mod
!
use lib_h5fortran_mod
use lib_unified_chars_mod
use h5fortran, only: hdf5_file
!
implicit none
integer         , intent(in)  :: ifil       ! IO number
character(len=*), intent(in)  :: filename
integer         , intent(out) :: iformat    ! Format type number
integer         , intent(in ) :: NMSG               ! Dimension of ier_msg
integer         , intent(out) :: ier_num            ! an error =0 if all is OK
character(len=*), dimension(NMSG), intent(out) :: ier_msg            ! Error messages
!
character(len=PREC_STRING) :: line
integer :: ios   ! I/O status
!
iformat = -1
!
read(ifil,'(a)', iostat=ios) line
LOOP_MAIN: do while(.not.is_iostat_end(ios))
  if(line(1:7) /= 'version') then
     if(line(1:7) /= 'version = diffpy.pdfgetx') then
        iformat = FILE_DIFFPY_PDFGETX2
        return
     endif
  elseif(line(1:3) == '#S ') then
     iformat = FILE_SPEC
     return
  endif
  read(ifil,'(a)', iostat=ios) line
enddo LOOP_MAIN
  if(is_iostat_end(ios)) then
     ier_num = -31
     ier_msg = 'File starts with [DEFAULT]'
     return
   endif
!

end subroutine detect_pdfgetx
!
!*******************************************************************************
!
end module lib_detect
