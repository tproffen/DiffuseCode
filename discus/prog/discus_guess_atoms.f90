module guess_atoms_mod
!-
!  Routines to estimate which chemical element is "hidden" behind user provided atom names like H, H1, H1+, H1N1 etc
!
contains
!
!*******************************************************************************
!
subroutine guess_atom_all
!-
!  Loop over all atom types to perform guess work.
!  Atoms whose chemical element can be guessed are assigned an equivalen 
!  scattering type. 
!
!  Is called by readcell and readstru, prior to fourier and powder run
!+
use crystal_mod
!
use element_data_mod
use errlist_mod
use lib_element_status_mod
use lib_errlist_func
!
implicit none
!
!
!integer, parameter :: IS_FAIL       = -1
!integer, parameter :: IS_NEUTRAL    =  1
!integer, parameter :: IS_ION        =  2
!integer, parameter :: IS_EQUIVALENT =  3
!integer, parameter :: IS_GUESS      =  4
!
character(len=4) :: trial  ! Test atom name
integer          :: is_cond     ! Dummy loop indices
integer          :: i      ! Dummy loop indices
!
loop_main: do i=1, cr_nscat
   call  guess_element(trial, is_cond, cr_at_lis(i), scat_equ=cr_scat_equ(i),   &
                                                     scat_equ_name=cr_at_equ(i))
!write(*,'(2a,l2,2a5,i3)') ' ATOM_NAME ',cr_at_lis(i), cr_scat_equ(i), cr_at_equ(i), trial, is_cond
   if(is_cond==IS_EQUIVALENT) then
      cycle loop_main
   elseif(is_cond==IS_NEUTRAL) then
      cycle loop_main
   elseif(is_cond==IS_ION) then
      cycle loop_main
   elseif(is_cond==IS_GUESS) then
      cr_scat_equ(i) = .true.       ! Automatically assign equivalent (neutral) atom name
      cr_at_equ(i)   = trial
      cycle loop_main
   elseif(is_cond==IS_GUESS_2) then
      cr_scat_equ(i) = .true.       ! Automatically assign equivalent (neutral) atom name
      cr_at_equ(i)   = trial
      cycle loop_main
   elseif(is_cond==IS_GUESS_1) then
      cr_scat_equ(i) = .true.       ! Automatically assign equivalent (neutral) atom name
      cr_at_equ(i)   = trial
      cycle loop_main
   elseif(is_cond==IS_FAIL) then
      ier_num = 9
      ier_typ = ER_APPL
      ier_msg(1) = 'Atom name: ' // cr_at_lis(i)
      call errlist
      cycle loop_main
   endif 
enddo loop_main
if(ier_num>0) then
   call errlist
endif
!
end subroutine guess_atom_all
!
!*******************************************************************************
!
end module guess_atoms_mod
