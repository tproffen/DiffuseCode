subroutine experi_mache_kdo(line, lend, length)
!-
!  Main command handling routine to EXPERI
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(inout) :: line     ! User command line
logical         , intent(inout) :: lend     ! End flag
integer         , intent(inout) :: length   ! command line length
!
lend = .TRUE.
!
end subroutine experi_mache_kdo
