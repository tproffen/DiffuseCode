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
use errlist_mod
use precision_mod
use support_mod
!
implicit none
!
character(len=*), intent(in) :: infile
character(len=*), intent(in) :: outfile   ! output file
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
integer, parameter :: IRD = 45
integer, parameter :: IWR = 46
!
!character(len=PREC_STRING)     ::  infile   !  input file
character(len=PREC_STRING), dimension(2) :: initial   ! first two lines
integer                        :: ios       ! input error signal
integer                        :: fcf_type  ! fcf file type
logical                        :: lexist    ! Check if file is present
!
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
read(IRD,'(a)', iostat=ios) initial(1)
if(ios/=0) then
   ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
   close(IRD)
   close(IWR)
   return
endif
!
read(IRD,'(a)', iostat=ios) initial(2)
if(ios/=0) then
   ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
   close(IRD)
   close(IWR)
   return
endif
!
if(initial(2)(1:56)=='# h,k,l, Fo-squared, sigma(Fo-squared), Fc and phi(calc)') then
   fcf_type = 6
   call lib_convert_shelx_fcf_6(IRD, IWR, initial,hkl_max)
elseif(initial(2)(1:66)=='# h,k,l, Fc-squared, Fo-squared, sigma(Fo-squared) and status flag') then
   fcf_type = 4
   call lib_convert_shelx_fcf_4(IRD, IWR, initial, hkl_max)
elseif(initial(2)(1:47)=='# Unique observed reflections after correcting') then
   fcf_type = 3
else
   fcf_type = 2
   call lib_convert_shelx_fcf_2(IRD, IWR, initial, hkl_max)
   if(ier_num/=0) then
      ier_msg(1) = infile((max(1,len_trim(infile)-39)):len_trim(infile))
      close(IRD)
      close(IWR)
      return
   endif
endif
close(IRD)
close(IWR)
!
end subroutine lib_convert_shelx
!
!*******************************************************************************
!
subroutine lib_convert_shelx_fcf_6(IRD, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 6 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: IRD
integer, intent(in) :: IWR
character(len=*), dimension(2), intent(in) :: initial   ! first two lines
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
character(len=PREC_STRING)       :: line     !Input lines
integer                          :: ios       ! input error signal
integer, dimension(3)            :: hkl       ! Miller indices
real(kind=PREC_DP), dimension(2) :: fobs      ! F(obs) sigma(Fobs)
!
loop_header: do
   read(IRD,'(a)') line
   if(line(1:18) == ' _refln_phase_calc') exit loop_header
enddo loop_header
!
hkl_max = 0
loop_main: do
   read(IRD,'(a)', iostat=ios) line
   if(is_iostat_end(ios)) exit loop_main
   read(line,*, iostat=ios) hkl, fobs
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
subroutine lib_convert_shelx_fcf_4(IRD, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 4 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: IRD
integer, intent(in) :: IWR
character(len=*), dimension(2), intent(in) :: initial   ! first two lines
!
character(len=PREC_STRING)       :: line     !Input lines
integer                          :: ios       ! input error signal
integer, dimension(3)            :: hkl       ! Miller indices
real(kind=PREC_DP), dimension(3) :: fobs      ! F(obs) sigma(Fobs)
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
loop_header: do
   read(IRD,'(a)') line
   if(line(1:18) == ' _refln_observed_status') exit loop_header
enddo loop_header
!
loop_main: do
   read(IRD,'(a)', iostat=ios) line
   if(is_iostat_end(ios)) exit loop_main
   read(line,*, iostat=ios) hkl, fobs
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
subroutine lib_convert_shelx_fcf_2(IRD, IWR, initial, hkl_max)
!-
!  Convert a SHELX fcf LIST 2 file
!
use errlist_mod
use precision_mod
!
implicit none
!
integer, intent(in) :: IRD
integer, intent(in) :: IWR
character(len=*), dimension(2), intent(in) :: initial   ! first two lines
integer, dimension(3), intent(out) :: hkl_max   ! maximum abs(h,k,l)
!
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
loop_main: do
   read(IRD,'(3i4,2f8.2,i4)', iostat=ios) hkl, fobs, phase
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
