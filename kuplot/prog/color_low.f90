module kuplot_color_low_mod
!-
!  Dummy routines if compiled wtihout graphics
!+
!
contains
!
!*******************************************************************************
!
subroutine inquire_color_index(i1, i2)
!
integer         , intent(out) :: i1
integer         , intent(out) :: i2
!
call PGQCOL(i1, i2)          ! Inquire color index range
!
end subroutine inquire_color_index
!
!*******************************************************************************
!
end module kuplot_color_low_mod

