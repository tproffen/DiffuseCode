module lib_timer_class
!-
! A generic timer class
!+
!
use precision_mod
!
implicit none
!
type :: timer
   private
   real(kind=prec_DP) :: start_time
   logical            :: lcpu
   logical            :: is_started
!
contains
!
   procedure, public :: reset_timer
   procedure, public :: start_timer
   procedure, public :: stop_timer
   procedure, public :: get_started
end type timer
!
contains
!
!*******************************************************************************
!
! Actual subroutines
!
!*******************************************************************************
!
subroutine reset_timer(this, lcpu)
!-
!  Initiate time to zero, start timer
!
class(timer) :: this
logical, intent(in) :: lcpu
!
!
this%lcpu       = lcpu
this%is_started = .false.
!
end subroutine reset_timer
!
!*******************************************************************************
!
subroutine start_timer(this, lcpu)
!-
!  Initiate time to zero, start timer
!
class(timer) :: this
logical, intent(in) :: lcpu
!
integer, dimension(8) :: wtimes
real(kind=PREC_SP)    :: ctimes
!
if(lcpu) then
   call cpu_time(ctimes)
   this%start_time = ctimes
   this%lcpu       = .true.
else
   call date_and_time(VALUES=wtimes)
   this%start_time = 86400.0_PREC_DP*wtimes(3) + 3600.0_PREC_DP*wtimes(5) &
      + 60.0_PREC_DP*wtimes(6) + real(wtimes(7),kind=PREC_DP) + 0.001_PREC_DP*wtimes(8)
   this%lcpu       = .false.
endif
!
this%is_started = .true.
!
end subroutine start_timer
!
!*******************************************************************************
!
real(kind=PREC_SP) function stop_timer(this)  result(elapsed)
!-
!  Stop timer, evaluate elapsed (CPU/Wall) time
!+
class(timer) :: this
!
integer, dimension(8) :: wtimes
real(kind=PREC_SP)    :: current_time
!
if(this%lcpu) then
   call cpu_time(current_time)
else
   call date_and_time(VALUES=wtimes)
!
   current_time = real( 86400.0_PREC_DP*wtimes(3) + 3600.0_PREC_DP*wtimes(5) &
      + 60.0_PREC_DP*wtimes(6) + real(wtimes(7),kind=PREC_DP) + 0.001_PREC_DP*wtimes(8), &
      kind=PREC_SP)
endif
this%is_started = .false.
elapsed = current_time - real(this%start_time, kind=PREC_SP)
!
end function stop_timer
!
!*******************************************************************************
!
logical function get_started(this) result(is_started)
!-
!  Returns if timer had been started
!+
class(timer) :: this
!
is_started = this%is_started
!
end function get_started
!
!*******************************************************************************
!
end module lib_timer_class
