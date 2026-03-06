module kuplot_draw_mod
!
! dummy subroutines if compiled without graphics
!+
contains
!
!******7****************************************************************
!
subroutine do_mouse (zeile, lp)
!
character(len=*), intent(in) :: zeile
integer         , intent(in) :: lp
!
end subroutine do_mouse
!

end module kuplot_draw_mod
