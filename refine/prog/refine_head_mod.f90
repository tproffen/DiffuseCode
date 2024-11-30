module refine_head_mod
!
! Variables / routines to create header lines for "refine_best.mac""efine_new.res"
!
use precision_mod
!
private
logical :: refine_head_i = .true.!  Initialize header
logical :: refine_foot_i = .true.!  Initialize header
logical :: refine_head_a = .true.!  Accumulate command lines into header
logical :: refine_foot_a = .false.!  Accumulate command lines into footer; After newparam
integer :: refine_head_n = 0     !  Number of lines in header
integer :: refine_foot_n = 0     !  Number of lines in header
character(len=PREC_STRING), dimension(:), allocatable :: refine_header ! Actual header lines
character(len=PREC_STRING), dimension(:), allocatable :: refine_footer ! Actual footer lines
!
public refine_head_a
public refine_foot_a
public accumulate_header
public accumulate_footer
public write_footer
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
subroutine accumulate_footer(line)
!-
!  Copy the line into the footer, expand if needed
!+
implicit none
!
character(len=*), intent(in) :: line
character(len=PREC_STRING), dimension(:), allocatable :: tmp ! Actual footer lines
!
if(refine_foot_i) then                ! Initialize footer
   if(allocated(refine_footer)) deallocate(refine_footer)
   allocate(refine_footer(10))
   refine_foot_i = .false.
endif
!
if(allocated(refine_footer)) then     ! Header exists
   if(ubound(refine_footer,1)<refine_foot_n+1) then  ! File needs space
      call move_alloc(refine_footer, tmp)
      allocate(refine_footer(refine_foot_n+10))
      refine_footer(1:refine_foot_n) = tmp(1:refine_foot_n)
      refine_footer(refine_foot_n+1:) = ' '
      deallocate(tmp)
   endif
endif
!
refine_foot_n = refine_foot_n + 1
refine_footer(refine_foot_n) = line
!
end subroutine accumulate_footer
!
!*******************************************************************************
!
subroutine write_header(iwr, refine_mac, refine_plot_mac)
!-
!  Write the header lines into output file number iwr
!+
!
implicit none
!
integer, intent(in) :: iwr
character(len=*), intent(in) :: refine_mac
character(len=*), intent(in) :: refine_plot_mac
!
integer :: i   ! Dummy loop index
!
loop_main:do i=1, refine_head_n
  if(refine_header(i)(2:len_trim(refine_header(i))) == refine_mac(1:len_trim(refine_mac))) cycle loop_main
  if(refine_header(i)(2:len_trim(refine_header(i))) == refine_plot_mac(1:len_trim(refine_plot_mac))) cycle loop_main
  write(iwr,'(a)') refine_header(i)(1:len_trim(refine_header(i)))
enddo loop_main
!
end subroutine write_header
!
!*******************************************************************************
!
subroutine write_footer(iwr, refine_mac, refine_plot_mac)
!-
!  Write the footer lines into output file number iwr
!+
!
implicit none
!
integer, intent(in) :: iwr
character(len=*), intent(in) :: refine_mac
character(len=*), intent(in) :: refine_plot_mac
!
integer :: i   ! Dummy loop index
!
loop_main:do i=1, refine_foot_n
  if(refine_footer(i)(2:len_trim(refine_footer(i))) == refine_mac(1:len_trim(refine_mac))) cycle loop_main
  if(refine_footer(i)(2:len_trim(refine_footer(i))) == refine_plot_mac(1:len_trim(refine_plot_mac))) cycle loop_main
  write(iwr,'(a)') refine_footer(i)(1:len_trim(refine_footer(i)))
enddo loop_main
!
end subroutine write_footer
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
if(allocated(refine_footer)) then     ! Header exists
   deallocate(refine_footer)
endif
refine_foot_n = 0
refine_foot_i = .true.
!
end subroutine reset_header
!
!*******************************************************************************
!
end module refine_head_mod
