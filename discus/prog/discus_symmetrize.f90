module symmetrize_mod
!
! Routines to apply the space group symmetry to some data 
!+
!
contains
!
!*******************************************************************************
!
subroutine symmetrize_four_3d(inc, icenter, csf, llaue)
!
! Apply the point/Laue symmetry to 3D complex structure factor; replace in place
!+
!
use wyckoff_mod
!
use matrix_mod
use precision_mod
!
implicit none
!
integer              , dimension(3)                     , intent(in)    :: inc      ! Dimensions
integer              , dimension(3)                     , intent(in)    :: icenter  ! Location of 0,0,0
complex(kind=PREC_DP), dimension(inc(1), inc(2), inc(3)), intent(inout) :: csf      ! Data
logical                                                 , intent(in)    :: llaue    ! Apply centrosymmetry
!
integer :: i, h,k,l   ! Loop indices
real(kind=PREC_DP), dimension(3,3) :: mat    ! Symmetry matrix in direct     space
real(kind=PREC_DP), dimension(3,3) :: rmat   ! Symmetry matrix in reciprocal space
integer           , dimension(3,3) :: irmat  ! Symmetry matrix in reciprocal space
integer           , dimension(3)   :: hkl    ! vector 
integer           , dimension(3)   :: uvw    ! image of vector
integer           , dimension(3,3), parameter :: IMAT= reshape((/-1,0,0, 0,-1,0, 0,0,-1/), shape(IMAT)) ! Unit matrix
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: tcsf   ! Temporary data
integer              , dimension(:,:,:), allocatable :: mask   ! Weight for individual data points
!
!write(*,*) ' IN SYMMETRIZE ', spc_n, llaue
allocate(tcsf(inc(1), inc(2), inc(3)))
allocate(mask(inc(1), inc(2), inc(3)))
tcsf = csf                    ! Apply symmetry "1"
mask = 1
do i=2, spc_n                 ! Loop over all relevant matrices
   mat = spc_mat(1:3, 1:3, i)
   call matinv(mat, rmat)     ! Determine inverse matrix
   irmat = nint(rmat)
!write(*,*) ' Matrix ', irmat(1,:)
!write(*,*) ' Matrix ', irmat(2,:)
!write(*,*) ' Matrix ', irmat(3,:)
!write(*,*) 
   do l = 1, inc(3)           ! Loop over all data points h,k,l
      hkl(3) = l-icenter(3)
      do k = 1, inc(2)
         hkl(2) = k-icenter(2)
         do h = 1, inc(1)
            hkl(1) = h-icenter(1)
            uvw = matmul(irmat, hkl) + icenter    ! Calculate image
            if(uvw(1)>0 .and. uvw(1)<=inc(1)  .and. uvw(2)>0 .and. uvw(2)<=inc(2) .and. &
               uvw(3)>0 .and. uvw(3)<=inc(3) ) then      ! Image is within old data range
               tcsf(uvw(1), uvw(2), uvw(3)) = tcsf(uvw(1), uvw(2), uvw(3)) + csf(h,k,l)
               mask(uvw(1), uvw(2), uvw(3)) = mask(uvw(1), uvw(2), uvw(3)) + 1
            endif
         enddo
      enddo
   enddo
enddo
!
if(llaue) then
! Apply 1bar symmetry
   do l = 1, inc(3)
      hkl(3) = -(l-icenter(3))
      do k = 1, inc(2)
         hkl(2) = -(k-icenter(2))
         do h = 1, inc(1)
            hkl(1) = -(h-icenter(1))
            uvw = matmul(irmat, hkl) + icenter
            if(uvw(1)>0 .and. uvw(1)<=inc(1)  .and. uvw(2)>0 .and. uvw(2)<=inc(2) .and. &
               uvw(3)>0 .and. uvw(3)<=inc(3) ) then
               tcsf(uvw(1), uvw(2), uvw(3)) = tcsf(uvw(1), uvw(2), uvw(3)) + csf(h,k,l)
               mask(uvw(1), uvw(2), uvw(3)) = mask(uvw(1), uvw(2), uvw(3)) + 1
            endif
         enddo
      enddo
   enddo
endif
!
!write(*,*) minval(mask), maxval(mask)
csf = tcsf/ real(mask)
deallocate(tcsf)
!
end subroutine symmetrize_four_3d
!
!*******************************************************************************
!
end module symmetrize_mod
