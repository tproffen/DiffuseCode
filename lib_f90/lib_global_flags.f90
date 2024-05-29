module lib_global_flags_mod
!
! A few global flags to communicate accross discus_suite
!
use precision_mod
!
implicit none
integer(kind=PREC_INT_BYTE), parameter    :: LIB_GLOBAL_FLAGS_MAX = 8
integer(kind=PREC_INT_BYTE), dimension(LIB_GLOBAL_FLAGS_MAX) :: lib_global_flags = 0
!
! (1) : KUPLOT 1=> at exit: writes last data set into DISCUS_SUITE_DERIVATIVES/data.iiii.jjjj (REF_KID, REF_INDIV)
! (2) : unused
! (3) : unused
! (4) : unused
! (5) : unused
! (6) : unused
! (7) : unused
! (8) : unused
!
contains
!
!******************************************************************************+
!
subroutine lib_global_flags_reset
!-
!  Everything back to zero
!+
implicit none
lib_global_flags = 0
!
end subroutine lib_global_flags_reset
!
!******************************************************************************+
!
end module lib_global_flags_mod
