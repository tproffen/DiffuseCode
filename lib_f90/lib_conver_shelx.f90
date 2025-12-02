module lib_conv_shelx_mod
!-
!  Programs to convert SHELX fcf files into HKL files
!+
contains
!
!*******************************************************************************
!
subroutine lib_convert_shelx(infile, outfile, hkl_max)
!
! Main interface, receives file
! Calls specific converter routines
!
use blanks_mod
use errlist_mod
use precision_mod
use string_convert_mod
use support_mod
!
implicit none
!
character(len=*), intent(in) :: infile
character(len=*), intent(in) :: outfile   ! output file
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
integer, parameter :: IRD = 55
integer, parameter :: IWR = 56
!
!character(len=PREC_STRING)     ::  infile   !  input file
character(len=PREC_STRING), dimension(:), allocatable :: initial   ! first two lines
integer                        :: ios       ! input error signal
integer                        :: i,j       ! Dummy indices
integer                        :: nlines    ! Number of lines in file
integer                        :: length    ! line length
integer                        :: fcf_type  ! fcf file type
logical                        :: lexist    ! Check if file is present
!
fcf_type = 0
inquire(file=infile, exist=lexist)
if(.not.lexist) then
   ier_num = -1
   ier_num = ER_IO
   ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
   return
endif
!
call oeffne(IRD, infile, 'old')
if(ier_num/=0) then
   ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
   return
endif
!
call oeffne(IWR, outfile, 'unknown')
if(ier_num/=0) then
   ier_msg(1) = outfile((max(1,len_trim(outfile)-39)):len_trim(outfile))
   return
endif
!
nlines = 0
loop_count: do 
   read(IRD,'(a)', iostat=ios)
   if(is_iostat_end(ios)) exit loop_count
   nlines = nlines + 1
enddo loop_count
close(IRD)
!
call oeffne(IRD, infile, 'old')
!
allocate( initial(1:nlines))
initial = ' '
loop_read: do i=1, nlines
   read(IRD,'(a)', iostat=ios) initial(i)
   if(ios/=0) then
      ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
      close(IRD)
      close(IWR)
      return
   endif
   length = len_trim(initial(i))
   call rem_leading_bl(initial(i), length)
   call do_low(initial(i))
enddo loop_read
close(IRD)
!
fcf_type = 2
loop_identify: do i=1, nlines
   j = len_trim(initial(i))
   if(index(initial(i), '_shelx_refln_list_code')>0) then
      read(initial(i)(j:j),*) fcf_type
      exit loop_identify
   elseif(index(initial(i), 'h,k,l, fo-squared, sigma(fo-squared), fc and phi(calc)')>0) then
      fcf_type = 6
      exit loop_identify
   elseif(index(initial(i), 'h,k,l, fc-squared, fo-squared, sigma(Fo-squared) and status flag')>0) then
      fcf_type = 4
      exit loop_identify
   elseif(index(initial(i), 'unique observed reflections after correcting'          )>0) then
      fcf_type = 3
      exit loop_identify
   endif
enddo loop_identify
!
if(fcf_type==6) then
   call lib_convert_shelx_fcf_6(nlines, IWR, initial,hkl_max)
elseif(fcf_type==4) then
   call lib_convert_shelx_fcf_4(nlines, IWR, initial, hkl_max)
elseif(fcf_type==3) then
   call lib_convert_shelx_fcf_3(nlines, IWR, initial, hkl_max)
else
   fcf_type = 2
   call lib_convert_shelx_fcf_2(nlines, IWR, initial, hkl_max)
   if(ier_num/=0) then
      ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
      close(IWR)
      return
   endif
endif
close(IRD)
close(IWR)
!
deallocate(initial)
!
end subroutine lib_convert_shelx
!
!*******************************************************************************
!
subroutine lib_convert_shelx_fcf_6(nlines, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 6 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: nlines
integer, intent(in) :: IWR
character(len=*), dimension(nlines), intent(in) :: initial   ! first two lines
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
!character(len=PREC_STRING)       :: line     !Input lines
integer                          :: ios       ! input error signal
integer                          :: i, j      ! Dummy index
integer, dimension(3)            :: hkl       ! Miller indices
real(kind=PREC_DP), dimension(2) :: fobs      ! F(obs) sigma(Fobs)
!
j=1
loop_header: do, i=1, nlines
   if(index(initial(i), '_refln_phase_calc')>0) then
      j= i
      exit loop_header
   endif
enddo loop_header
!
hkl_max = 0
loop_main: do i=j+1, nlines
   read(initial(i),*, iostat=ios) hkl, fobs
   write(IWR,'(3i4,2f8.2)') hkl, fobs
   hkl_max(1) = max(hkl_max(1), abs(hkl(1)))
   hkl_max(2) = max(hkl_max(2), abs(hkl(2)))
   hkl_max(3) = max(hkl_max(3), abs(hkl(3)))
enddo loop_main 
!
end subroutine lib_convert_shelx_fcf_6
!
!*******************************************************************************
!
subroutine lib_convert_shelx_fcf_4(nlines, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 4 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: nlines
integer, intent(in) :: IWR
character(len=*), dimension(nlines), intent(in) :: initial   ! first two lines
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
integer                          :: ios       ! input error signal
integer                          :: i, j      ! Dummy indices
integer, dimension(3)            :: hkl       ! Miller indices
real(kind=PREC_DP), dimension(3) :: fobs      ! F(obs) sigma(Fobs)
!
j = 1
!
loop_header: do, i=1, nlines
   if(index(initial(i), '_refln_phase_calc')>0) then
      j= i
      exit loop_header
   endif
enddo loop_header
!
loop_main: do i=j+1, nlines
   read(initial(i),*, iostat=ios) hkl, fobs
   write(IWR,'(3i4,2f8.2)') hkl, fobs(2:3)
   hkl_max(1) = max(hkl_max(1), abs(hkl(1)))
   hkl_max(2) = max(hkl_max(2), abs(hkl(2)))
   hkl_max(3) = max(hkl_max(3), abs(hkl(3)))
enddo loop_main 
!
end subroutine lib_convert_shelx_fcf_4
!
!*******************************************************************************
!
subroutine lib_convert_shelx_fcf_3(nlines, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 3 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: nlines
integer, intent(in) :: IWR
character(len=*), dimension(nlines), intent(in) :: initial   ! first two lines
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
integer                          :: ios       ! input error signal
integer                          :: i, j      ! Dummy indices
integer, dimension(3)            :: hkl       ! Miller indices
real(kind=PREC_DP), dimension(2) :: fobs      ! F(obs) sigma(Fobs) A(calc), B(calc)
!
j = 1
!
loop_header: do, i=1, nlines
   if(index(initial(i), '_refln_phase_calc')>0) then
      j= i
      exit loop_header
   endif
enddo loop_header
!
loop_main: do i=j+1, nlines
   read(initial(i),*, iostat=ios) hkl, fobs
   write(IWR,'(3i4,2f8.2)') hkl, fobs(1)**2, fobs(1)*fobs(2)*sqrt(2.0_PREC_DP)
   hkl_max(1) = max(hkl_max(1), abs(hkl(1)))
   hkl_max(2) = max(hkl_max(2), abs(hkl(2)))
   hkl_max(3) = max(hkl_max(3), abs(hkl(3)))
enddo loop_main 
!
end subroutine lib_convert_shelx_fcf_3
!
!*******************************************************************************
!
subroutine lib_convert_shelx_fcf_2(nlines, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 2 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: nlines
integer, intent(in) :: IWR
character(len=*), dimension(nlines), intent(in) :: initial   ! first two lines
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
integer                          :: i         ! dummy index 
integer                          :: ios       ! input error signal
integer, dimension(3)            :: hkl       ! Miller indices
real(kind=PREC_DP), dimension(2) :: fobs      ! F(obs) sigma(Fobs)
integer                          :: phase     ! Reflection phase angle
!
read(initial(1)(1:32),'(3i4,2f8.2,i4)', iostat=ios) hkl, fobs, phase
if(ios/=0) then
   ier_num = -3
   ier_typ = ER_IO
   ier_msg(2) = 'File does not seem to be a LIST 2 type'
   return
endif
write(IWR,'(3i4,2f8.2,i4)') hkl, fobs(1)**2, 2.0_PREC_DP*fobs(2)
read(initial(2)(1:32),'(3i4,2f8.2,i4)', iostat=ios) hkl, fobs, phase
write(IWR,'(3i4,2f8.2,i4)') hkl, fobs(1)**2, 2.0_PREC_DP*fobs(2)
!
loop_main: do i=3, nlines
   read(initial(i),'(3i4,2f8.2,i4)', iostat=ios) hkl, fobs, phase
   if(is_iostat_end(ios)) exit loop_main
   fobs(2) = abs(2.0_PREC_DP*fobs(1)*fobs(2))   ! sigma(fobs) => sigma(I)
   write(IWR,'(3i4,2f8.2,i4)') hkl, fobs
   hkl_max(1) = max(hkl_max(1), abs(hkl(1)))
   hkl_max(2) = max(hkl_max(2), abs(hkl(2)))
   hkl_max(3) = max(hkl_max(3), abs(hkl(3)))
enddo loop_main 
!
end subroutine lib_convert_shelx_fcf_2
!
!*******************************************************************************
!
end module lib_conv_shelx_mod
