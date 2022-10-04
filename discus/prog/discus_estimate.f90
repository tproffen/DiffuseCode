module dis_estimate_mod
!
! Procedures to estimate number of unit cells, average atoms per unit cell
!
private
!
public estimate_ncells
public estimate_ncatom
!
contains
!
!*******************************************************************************
!
subroutine estimate_ncells(ncell_out, pdt_dims, pdt_ilow, pdt_ihig, pdt_ncells)
!-
!  Try to estimate the number of unit cells
!+
use crystal_mod
!
use precision_mod
!
implicit none
!
integer, dimension(3), intent(out) :: ncell_out
real(kind=PREC_DP), dimension(3, 2), intent(inout) :: pdt_dims
integer           , dimension(3)   , intent(out)   :: pdt_ilow        ! Unit cell dimensions in periodic
integer           , dimension(3)   , intent(out)   :: pdt_ihig        ! low and high inidce
integer                            , intent(out)   :: pdt_ncells     ! Number of cells in periodic crystal volume
!
integer :: i          ! Dummy counter
!
pdt_dims(:,1) = cr_pos(:,1)           ! Initialize the min/ max dimensions
pdt_dims(:,2) = cr_pos(:,1)
!
do i=1, cr_natoms
   pdt_dims(1,1) = min(pdt_dims(1,1), cr_pos(1,i))
   pdt_dims(1,2) = max(pdt_dims(1,2), cr_pos(1,i))
   pdt_dims(2,1) = min(pdt_dims(2,1), cr_pos(2,i))
   pdt_dims(2,2) = max(pdt_dims(2,2), cr_pos(2,i))
   pdt_dims(3,1) = min(pdt_dims(3,1), cr_pos(3,i))
   pdt_dims(3,2) = max(pdt_dims(3,2), cr_pos(3,i))
enddo
!
!write(*,'(a, f9.5, 2x, f9.5)') ' dims x ', pdt_dims(1,:)
!write(*,'(a, f9.5, 2x, f9.5)') ' dims y ', pdt_dims(2,:)
!write(*,'(a, f9.5, 2x, f9.5)') ' dims z ', pdt_dims(3,:)
!
pdt_ilow = 0
pdt_ihig = 0
pdt_ncells = 1
do i=1, 3
   if(nint(pdt_dims(i,2)-pdt_dims(i,1))>2) then
!  pdt_ilow(i) = nint(pdt_dims(i,1))! + 1
!  pdt_ihig(i) = nint(pdt_dims(i,2))! - 1
      pdt_ilow(i) = nint(pdt_dims(i,1))! + 1
      pdt_ihig(i) = pdt_ilow(i) + int(pdt_dims(i,2)-pdt_dims(i,1))! - 1
   else
      pdt_ilow(i) = nint(pdt_dims(i,1))
      pdt_ihig(i) = pdt_ilow(i) + int(pdt_dims(i,2)-pdt_dims(i,1))
   endif
   pdt_ncells = pdt_ncells * (pdt_ihig(i)-pdt_ilow(i)+1)
   ncell_out(i) = (pdt_ihig(i)-pdt_ilow(i)+1)
!write(*,'( i3,  i3, i3 )' ) i, pdt_ilow(i), pdt_ihig(i)
enddo
!!!pdt_usr_ncell = .false.              ! Unit cells were estimates automatically
!
end subroutine estimate_ncells
!
!*******************************************************************************
!
subroutine estimate_ncatom(aver, sigma, pdt_ilow, pdt_ihig, pdt_ncells)
!-
!  Try to estimate the number of atoms per unit cell
!+
use crystal_mod
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(out) :: aver
real(kind=PREC_DP), intent(out) :: sigma
integer           , dimension(3)   , intent(out)   :: pdt_ilow        ! Unit cell dimensions in periodic
integer           , dimension(3)   , intent(out)   :: pdt_ihig        ! low and high inidce
integer                            , intent(out)   :: pdt_ncells     ! Number of cells in periodic crystal volume
!
integer :: i,j,k, l   ! Dummy counter
integer :: natoms     ! Total atoms in search window
integer :: maxn       ! Max number per unit cell
integer           , dimension(3)                  :: ixyz         ! Atom is in this unit cell
integer           , dimension(0:500)              :: counter
integer           , dimension(:,:,:), allocatable :: pdt_cell     ! Atoms per cell
!
allocate(pdt_cell(pdt_ilow(1):pdt_ihig(1), pdt_ilow(2):pdt_ihig(2),pdt_ilow(3):pdt_ihig(3)))
pdt_cell = 0
do i=1, cr_natoms
   ixyz(1) = int(cr_pos(1,i) + 0.05 - pdt_ilow(1)) + pdt_ilow(1)             ! add 0.5 to avoid loosing atoms at the low edge
   ixyz(2) = int(cr_pos(2,i) + 0.05 - pdt_ilow(2)) + pdt_ilow(2)
   ixyz(3) = int(cr_pos(3,i) + 0.05 - pdt_ilow(3)) + pdt_ilow(3)
   if(ixyz(1)>=pdt_ilow(1) .and. ixyz(1)<=pdt_ihig(1)   .and.          &
      ixyz(2)>=pdt_ilow(2) .and. ixyz(2)<=pdt_ihig(2)   .and.          &
      ixyz(3)>=pdt_ilow(3) .and. ixyz(3)<=pdt_ihig(3)         ) then
      pdt_cell(ixyz(1), ixyz(2), ixyz(3)) = pdt_cell(ixyz(1), ixyz(2), ixyz(3)) + 1
   endif
enddo
!
maxn = maxval(pdt_cell)
natoms = sum(pdt_cell)
pdt_ncells =  (pdt_ihig(1)-pdt_ilow(1)+1)*(pdt_ihig(2)-pdt_ilow(2)+1)*(pdt_ihig(3)-pdt_ilow(3)+1)
!write(*,*) ' Nunit ', minval(pdt_cell), maxval(pdt_cell), sum(pdt_cell),  &
!  (pdt_ihig(1)-pdt_ilow(1)+1)*(pdt_ihig(2)-pdt_ilow(2)+1)*(pdt_ihig(3)-pdt_ilow(3)+1)
counter = 0
do k=pdt_ilow(3), pdt_ihig(3)
   do j=pdt_ilow(2), pdt_ihig(2)
      do i=pdt_ilow(1), pdt_ihig(1)
         l = pdt_cell(i,j,k)
         counter(l) = counter(l) + 1
      enddo
   enddo
enddo
aver  = 0.0
sigma = 0.0
do i=0,maxn
  aver = aver + real(i*counter(i), kind=PREC_DP)
enddo
aver = aver/real(pdt_ncells, kind=PREC_DP)
do i=0,maxn
   sigma = sigma + counter(i)*(aver-i)**2
enddo
sigma = sqrt(sigma)/real(pdt_ncells*(pdt_ncells-1), kind=PREC_DP)
!
deallocate(pdt_cell)
!
end subroutine estimate_ncatom
!
!*******************************************************************************
!
end module dis_estimate_mod
