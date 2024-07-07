module reduce_atoms_mod
!
!  Routines to modify the number of atom/molecule types on the fly
!-
private
!
public reduce_atoms
!
contains
!
!*******************************************************************************
!
subroutine reduce_atoms(MAXMASK, uni_mask)
!-
! Reduce the number of atom types according to the uniqe mask
!+
!
use crystal_mod
use discus_allocate_appl_mod
use atom_line_mod
!
use precision_mod
!
implicit none
!
integer                      , intent(in) :: MAXMASK   ! Array dimension
logical, dimension(0:MAXMASK), intent(in) :: uni_mask  ! The unique mask
!
integer :: lb
integer :: ub
!
character(len=4)  , dimension(:), allocatable :: at_lis
real(kind=PREC_DP), dimension(:), allocatable :: at_dw
real(kind=PREC_DP), dimension(:), allocatable :: at_occ
!
character(len=4)                     :: nw_name
real(kind=PREC_DP)                   :: nw_dw
real(kind=PREC_DP)                   :: nw_occ
!
integer :: i     ! Loop indices
integer :: is_type ! Derived atom type
integer :: nscat ! number atom types as to unique mask
integer :: nscat_start
!
lb = lbound(cr_at_lis,1)
ub = ubound(cr_at_lis,1)+5
allocate(at_lis (lb:ub))
allocate(at_dw  (lb:ub))
allocate(at_occ (lb:ub))
!
if(lb==0) then
   at_lis(lb) = 'VOID'
   at_dw(lb)  = 0.0D0
   at_occ(lb) = 1.0D0
endif
!
nscat = 0
nscat_start = lb
!
do i=1, cr_natoms
   nw_name = cr_at_lis(cr_iscat(1,i))         ! Use present atom name
   nw_dw   = cr_dw(cr_iscat(1,i))             ! use present ADP
   nw_occ  = cr_occ(cr_iscat(1,i))            ! use present ADP
   is_type = atom_get_type(MAXSCAT, nscat_start, nscat, MAXMASK,      &
                       at_lis, at_dw, at_occ,        &
                       nw_name, nw_dw, nw_occ, uni_mask)
   if(is_type==-1) then                     ! New atom type
      nscat = nscat + 1
      if(nscat>ubound(cr_at_lis,1)) then
         call alloc_crystal_scat(nscat)
      endif
      at_lis(nscat) = nw_name
      at_dw (nscat) = nw_dw
      at_occ(nscat) = nw_occ
      is_type = nscat
   endif
   cr_iscat(1,i) = is_type
enddo
!
cr_nscat = nscat
cr_at_lis = at_lis
cr_dw     = at_dw
cr_occ    = at_occ
as_at_lis = at_lis
as_dw     = at_dw
as_occ    = at_occ
!
deallocate(at_lis)
deallocate(at_dw )
deallocate(at_occ)
!
end subroutine reduce_atoms
!
!
!*******************************************************************************
!
end module reduce_atoms_mod
