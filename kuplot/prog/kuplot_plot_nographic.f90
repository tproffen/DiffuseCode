module kuplot_plot_mod
!-
!  Dummy rou7tines if compiled without graphics option
!
contains
!
!*******************************************************************************
!
subroutine do_plot (lmenu)
!
logical , intent(in) :: lmenu
!
end subroutine do_plot
!
!*******************************************************************************
!
subroutine do_hardcopy (befehl, zeile, lbef, lp)
!
character(len=*), intent(in) :: befehl
character(len=*), intent(in) :: zeile
integer,          intent(in) :: lbef
integer,          intent(in) :: lp

end subroutine do_hardcopy
!
!*******************************************************************************
!
end module kuplot_plot_mod
