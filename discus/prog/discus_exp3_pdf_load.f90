module exp3_pdf_load
!-
!   Routines to load / merge data
!+
!
contains
!
!*******************************************************************************
!
subroutine load_merge_hklf4(infile)
!-
!  Load and merge an HKLF4 file
!+
!
use exp3_pdf_data
use errlist_mod
use support_mod
use precision_mod
!
implicit none
!
integer, parameter ::IRD=66
real(kind=PREC_DP), parameter :: NULL = 0.0D0
!
character(len=*), intent(in) :: infile    ! Main input file
!
character(len=32) ::frmt  ! Format string
character(len=PREC_STRING) :: ofile
integer :: i,j,k   ! Dummy loop indices
integer :: ios     ! I/O status
integer :: nref    ! Number of reflections
integer, dimension(3) :: hkl  ! Indices
integer, dimension(3) :: hkl_max  ! Indices
integer, dimension(3,2) :: hkl_min_max  ! Min/Max Indices
integer, dimension(:,:,:), allocatable :: data_mult   ! How many values at this hkl
real(kind=PREC_DP) :: inte, sigma
!
call oeffne(IRD, infile, 'old')
if(ier_num /= 0) then
   ier_msg(1) = ' file :'
   return
endif
nref = 0
frmt = '(3i4, f8.2, f7.2)' 
!
nref = 0
hkl         = 0
hkl_min_max = 0
read(IRD, frmt , iostat=ios) hkl, inte, sigma
loop_first: do
   if(is_iostat_end(ios) ) exit loop_first
   if(hkl(1) == NULL .and. hkl(2)==NULL .and. HKL(3)==NULL) exit LOOP_first
   nref = nref + 1
   hkl_min_max(1,1) = min(hkl(1), hkl_min_max(1,1))
   hkl_min_max(2,1) = min(hkl(2), hkl_min_max(2,1))
   hkl_min_max(3,1) = min(hkl(3), hkl_min_max(3,1))
   hkl_min_max(1,2) = max(hkl(1), hkl_min_max(1,2))
   hkl_min_max(2,2) = max(hkl(2), hkl_min_max(2,2))
   hkl_min_max(3,2) = max(hkl(3), hkl_min_max(3,2))
   read(IRD, frmt , iostat=ios) hkl, inte, sigma
enddo loop_first
!
close(unit=IRD)
!
hkl_max(1) = max(abs(hkl_min_max(1,1)), abs(hkl_min_max(1,2)))
hkl_max(2) = max(abs(hkl_min_max(2,1)), abs(hkl_min_max(2,2)))
hkl_max(3) = max(abs(hkl_min_max(3,1)), abs(hkl_min_max(3,2)))
if(allocated(exp3_data))  deallocate(exp3_data)
if(allocated(exp3_sigma)) deallocate(exp3_sigma)
if(allocated(exp3_mask))  deallocate(exp3_mask)
if(allocated(data_mult))  deallocate(data_mult)
allocate(exp3_data (-hkl_max(1):hkl_max(1), -hkl_max(2):hkl_max(2), -hkl_max(3):hkl_max(3)))
allocate(exp3_sigma(-hkl_max(1):hkl_max(1), -hkl_max(2):hkl_max(2), -hkl_max(3):hkl_max(3)))
allocate(exp3_mask (-hkl_max(1):hkl_max(1), -hkl_max(2):hkl_max(2), -hkl_max(3):hkl_max(3)))
allocate(data_mult (-hkl_max(1):hkl_max(1), -hkl_max(2):hkl_max(2), -hkl_max(3):hkl_max(3)))
exp3_data = 0.0D0
exp3_mask = .false.
data_mult = 0
!
call oeffne(IRD, infile, 'old')
loop_read: do i=1, nref
   read(IRD, frmt , iostat=ios) hkl, inte, sigma
   exp3_data (hkl(1), hkl(2), hkl(3)) = exp3_data (hkl(1), hkl(2), hkl(3)) + inte
   exp3_sigma(hkl(1), hkl(2), hkl(3)) = exp3_sigma(hkl(1), hkl(2), hkl(3)) + sigma
   data_mult (hkl(1), hkl(2), hkl(3)) = data_mult (hkl(1), hkl(2), hkl(3)) + 1
   exp3_mask (hkl(1), hkl(2), hkl(3)) = .true.
! Quick and dirty ... Do Friedels law
   exp3_data (-hkl(1), -hkl(2), -hkl(3)) = exp3_data (-hkl(1), -hkl(2), -hkl(3)) + inte
   exp3_sigma(-hkl(1), -hkl(2), -hkl(3)) = exp3_sigma(-hkl(1), -hkl(2), -hkl(3)) + sigma
   data_mult (-hkl(1), -hkl(2), -hkl(3)) = data_mult (-hkl(1), -hkl(2), -hkl(3)) + 1
   exp3_mask (-hkl(1), -hkl(2), -hkl(3)) = .true.
enddo loop_read
close(unit=IRD)
exp3_data (0,0,0) = 0.0D0
exp3_sigma(0,0,0) = 0.0D0
exp3_mask (0,0,0) = .false.
where(data_mult>1) exp3_data  = exp3_data /data_mult
where(data_mult>1) exp3_sigma = exp3_sigma/data_mult
!
ofile = 'merged.hkl'
call oeffne(IRD, ofile, 'unknown')
do k= -hkl_max(3), hkl_max(3)
   do j= -hkl_max(2), hkl_max(2)
      do i= -hkl_max(1), hkl_max(1)
         if(exp3_mask(i,j,k)) then
            write(IRD, frmt) i,j,k, exp3_data(i,j,k), exp3_sigma(i,j,k)
         endif
      enddo
   enddo
enddo
write(IRD, frmt) 0,0,0,  0.00D0, 0.00D0
close(unit=IRD)
!
end subroutine load_merge_hklf4
!
!*******************************************************************************
!
end module exp3_pdf_load
