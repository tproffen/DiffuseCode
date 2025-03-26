module four_finufft
!-
!  FOURIER transformatons with finufft library
!          1/2/3 dimensional non-uniform-FFT
!          At the moment axes must be parallel to a*; b*; c*
!          If less than 1/cr_icc(:) grid points per reciprocal unt cell are
!          used, the coordinates are multiplied with scales(:) 
!          Locally an allocatable array fk is used, data are copied into 
!          the "diffuse" variable tcsf in sequence with top the fastes and
!          abscissa the slowest index.
!
!+
use finufft_mod
!
contains
!
!*******************************************************************************
!
subroutine four_strucf_1d(MAXATOMS, nat, inc, MAXQXY, xpos, occ, iscales, lscales, tcsf)
!-
!  Calculate a 1D-NUFFT via finufft
!+
!
use iso_c_binding
!
use errlist_mod
use precision_mod
use wink_mod
!
implicit none
!
integer                                   , intent(in)  :: MAXATOMS ! Dimension of xpos
integer                                   , intent(in)  :: nat      ! Actual atom number
integer                                   , intent(in)  :: inc      ! Number data points along axes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Actual Reciprocal space dimensions
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer                                   , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
integer                                   , intent(in)  :: lscales  ! Scaling if scale > MAXSCALE
complex(kind=PREC_DP), dimension(MAXQXY(1),1,1), intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: n          ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
!
type(finufft_opts) :: opts
!
iflag = 1
M   = nat
N   = (inc-1) * lscales          + 1
allocate(cj(M))
allocate(fk(N))
!
cj = cmplx(occ, 0.0D0, kind=PREC_DP)
tol = 1d-9
!
call finufft_default_opts(opts)                            ! Use default options 
call finufft1d1(M, xpos(1:nat), cj, iflag, tol, N, fk, opts, ier)   ! Do 1D NUFFT
!
if(ier==0) then
   tcsf(1:MAXQXY(1), 1, 1) = fk( 1: N:lscales)
endif
!
deallocate(cj)
deallocate(fk)
!
end subroutine four_strucf_1d
!
!******************************************************************************
!
subroutine four_strucf_2d(MAXATOMS, nat, inc, MAXQXY, xpos, ypos, occ, iscales, lscales, tcsf)
!-
!  Calculate a 2D-NUFFT via finufft
!+
!
use errlist_mod
use precision_mod
use wink_mod
!
implicit none
!
integer                                   , intent(in)  :: MAXATOMS ! Dimension of xpos
integer                                   , intent(in)  :: nat      ! Actual atom number
integer              , dimension(2)       , intent(in)  :: inc      ! Number data points along axes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Total number in reciprocal space
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: ypos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer              , dimension(2)       , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
integer              , dimension(2)       , intent(in)  :: lscales  ! Scaling if scales > MAXSCALE
complex(kind=PREC_DP), dimension(MAXQXY(1), MAXQXY(2), 1)  , intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: ms, mt     ! Number of uniform points in Fourier space
integer(kind=PREC_INT_LONG) :: h, k, i    ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: hh, kk     ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: nj         ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
complex(kind=PREC_DP), dimension(:,:), allocatable :: tt  ! Complex structure factor
!
type(finufft_opts) :: opts                  ! Use default options
!
iflag = 1
M   = nat
ms  = (inc(1)-1) * lscales(1) + 1
mt  = (inc(2)-1) * lscales(2) + 1
nj = ms*mt
!
allocate(cj(M))
allocate(fk(nj))
!
cj = cmplx(occ, 0.0D0, kind=PREC_DP)
tol = 1d-9
!
call finufft_default_opts(opts)
call finufft2d1(M, xpos(1:nat), ypos(1:nat), cj, iflag, tol, ms, mt, fk, opts, ier)
deallocate(cj)
i = 0
!
! Copy fk onto a 2D grid
!
allocate(tt(ms, mt))
do h=1, ms
   do k=1, mt
      tt(h,k) = fk((k-1)*ms+h)
   enddo
enddo
!
! Cut the 2D grid with steps of lscales(*)
!
hh = 0
do h=1, ms, lscales(1)
   hh = hh + 1
   kk = 0
   do k=1, mt, lscales(2)
      kk = kk + 1
      tcsf(hh,kk,1) = tt(h,k)
   enddo
enddo
!
deallocate(fk)
deallocate(tt)
!
end subroutine four_strucf_2d
!
!******************************************************************************
!
subroutine four_strucf_3d(MAXATOMS, nat, inc, MAXQXY, xpos, ypos, zpos, occ, iscales, lscales, tcsf)
!-
!  Calculate a 3D-NUFFT via finufft
!+
!
use errlist_mod
use precision_mod
use wink_mod
use support_mod
!
implicit none
!
integer                                   , intent(in)  :: MAXATOMS ! Dimension of xpos
integer                                   , intent(in)  :: nat      ! Actual atom number
integer              , dimension(3)       , intent(in)  :: inc      ! Number data points along exes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Total data points reciprocal space
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: ypos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(inout)  :: zpos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer              , dimension(3)       , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
integer              , dimension(3)       , intent(in)  :: lscales  ! Scaling if scale > MAXSCALE
complex(kind=PREC_DP), dimension(MAXQXY(1),MAXQXY(2),MAXQXY(3))  , intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: ms, mt, mu ! Number of uniform points in Fourier space
integer(kind=PREC_INT_LONG) :: h, k, l, i ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: hh,kk,ll   ! Loop indices             in Fourier space
integer(kind=PREC_INT_LONG) :: nj         ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: ttcsf  ! Complex structure factor
!
type(finufft_opts) :: opts                  ! Use default options
!
iflag = 1
M   = nat
ms  = inc(1)
mt  = inc(2)
mu  = inc(3)
ms  = (inc(1)-1) * lscales(1) + 1
mt  = (inc(2)-1) * lscales(2) + 1
mu  = (inc(3)-1) * lscales(3) + 1
nj = ms*mt*mu
allocate(cj(M))
allocate(fk(nj))
!
cj = cmplx(occ, 0.0D0, kind=PREC_DP)
tol = 1.0D-7
!
call finufft_default_opts(opts)
call finufft3d1(M, xpos(1:nat), ypos(1:nat), zpos(1:nat), cj, iflag, tol, ms, mt, mu, fk, opts, ier)
deallocate(cj)
!
! copy fk ont 3D grid
!
allocate(ttcsf(ms,mt,mu))
i = 0
do l=1, mu
   do k=1, mt
      do h=1,ms
         i = i+1
         ttcsf(h ,k ,l ) = fk(i)
      enddo
   enddo
enddo
!
! Cut 3D grid at steps lscales(*)
l = 0
do ll=1,mu,lscales(3)
   l = l + 1
   k = 0
   do kk=1,mt,lscales(2)
      k = k + 1
      h = 0
      do hh=1,ms,lscales(1)
         h = h + 1
         tcsf(h,k,l) = ttcsf(hh,kk,ll)
      enddo
   enddo
enddo
!
deallocate(fk)
deallocate(ttcsf)
!
!
end subroutine four_strucf_3d
!
!******************************************************************************
!
end module four_finufft
