module kuplot_select_device_mod
!
! Select PGPLOT device
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
CALL PGSLCT(dev_id)
CALL PGCLOS
!
end subroutine kuplot_select_device
subroutine kuplot_pgend
!
call PGEND
!
end subroutine kuplot_pgend
!
!*******************************************************************************
!
end module kuplot_select_device_mod
