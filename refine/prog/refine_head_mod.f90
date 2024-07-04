module refine_head_mod
!
! Variables / routines to create header lines for "refine_best.mac""efine_new.res"
!
use precision_mod
!
private
logical :: refine_head_i = .true.!  Initialize header
logical :: refine_head_a = .true.!  Accumulate command lines into header
integer :: refine_head_n = 0     !  Number of lines in header
character(len=PREC_STRING), dimension(:), allocatable :: refine_header ! Actual header lines
!
public refine_head_a
public accumulate_header
public write_header
public reset_header
!
contains
!
!*******************************************************************************
!
subroutine accumulate_header(line)
!-
!  Copy the line into the header, expand if needed
!+
implicit none
!
character(len=*), intent(in) :: line
character(len=PREC_STRING), dimension(:), allocatable :: tmp ! Actual header lines
!
if(refine_head_i) then                ! Initialize header
   if(allocated(refine_header)) deallocate(refine_header)
   allocate(refine_header(10))
   refine_head_i = .false.
endif
!
if(allocated(refine_header)) then     ! Header exists
   if(ubound(refine_header,1)<refine_head_n+1) then  ! File needs space
      call move_alloc(refine_header, tmp)
      allocate(refine_header(refine_head_n+10))
      refine_header(1:refine_head_n) = tmp(1:refine_head_n)
      refine_header(refine_head_n+1:) = ' '
      deallocate(tmp)
   endif
endif
!
refine_head_n = refine_head_n + 1
refine_header(refine_head_n) = line
!
end subroutine accumulate_header
!
!*******************************************************************************
!
subroutine write_header(iwr)
!-
!  Write the header lines into output file number iwr
!+
!
implicit none
!
integer, intent(in) :: iwr
!
integer :: i   ! Dummy loop index
!
do i=1, refine_head_n
  write(iwr,'(a)') refine_header(i)(1:len_trim(refine_header(i)))
enddo
!
end subroutine write_header
!
!*******************************************************************************
!
subroutine reset_header
!-
!  Deallocate header
!+
!
implicit none
!
if(allocated(refine_header)) then     ! Header exists
   deallocate(refine_header)
endif
refine_head_n = 0
refine_head_i = .true.
!
end subroutine reset_header
!
!*******************************************************************************
!
end module refine_head_mod
