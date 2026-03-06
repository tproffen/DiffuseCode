module kuplot_select_device_mod
!
! Dummy routine to select and close PGPLOT
!
contains
!
!*******************************************************************************
!
subroutine kuplot_select_device(dev_id)
!-
! Select and open the device
!+
integer, intent(inout) :: dev_id
!
dev_id = 0
!
end subroutine kuplot_select_device
!
!*******************************************************************************
!
subroutine kuplot_pgend
!
continue
!
end subroutine kuplot_pgend
!
!*******************************************************************************
!
end module kuplot_select_device_mod
