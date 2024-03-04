module experi_setup_mod
!
contains
!
!*******************************************************************************
!
subroutine experi_setup
!-
!  Make initial set up for EXPERI
!+
!
use prompt_mod
!
implicit none
!
pname             = 'experi'
pname_cap         = 'EXPERI'
prompt            = pname
prompt_status     = PROMPT_ON
prompt_status_old = PROMPT_ON
!
! call experi_alloc_default
! call experi_initarrays
!
end subroutine experi_setup
!
!*******************************************************************************
!
end module experi_setup_mod
