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
subroutine four_strucf_1d(MAXATOMS, nat, icell, inc, MAXQXY, xpos, occ, iscales, scales, tcsf)
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
integer                                   , intent(in)  :: icell    ! Actual unit cell number
integer                                   , intent(in)  :: inc      ! Number data points along axes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Actual Reciprocal space dimensions
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(in)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer                                   , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
real(kind=PREC_DP)                        , intent(in)  :: scales   ! Scaling
complex(kind=PREC_DP), dimension(MAXQXY(1),1,1)  , intent(out) :: tcsf     ! Resulting structure factor
!
integer                     :: ier   = 0  ! Return error flag from finufft
integer                     :: iflag = 1  ! Flag to FINUFFT sign of i
integer(kind=PREC_INT_LONG) :: m          ! Number of non-uniform points == number atoms
integer(kind=PREC_INT_LONG) :: n          ! Number of uniform points in Fourier space
real(kind=PREC_DP)                               :: tol ! Required tolerance
real(kind=PREC_DP)   , dimension(:), allocatable :: xj  ! Coordinates
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
!
type(finufft_opts) :: opts
!
iflag = 1
M   = nat
N   = inc
allocate(xj(M))
allocate(cj(M))
allocate(fk(N))
xj = zpi*xpos(1:nat)/real(icell,kind=PREC_DP)*scales
cj = cmplx(occ, 0.0D0)
tol = 1d-9
!
call finufft_default_opts(opts)                            ! Use default options 
call finufft1d1(M, xj, cj, iflag, tol, N, fk, opts, ier)   ! Do 1D NUFFT
!
if(ier==0) then
   tcsf(1:MAXQXY(1), 1, 1) = fk(1:N:iscales)
endif
!
deallocate(xj)
deallocate(cj)
deallocate(fk)
!
end subroutine four_strucf_1d
!
!******************************************************************************
!
subroutine four_strucf_2d(MAXATOMS, nat, icell, inc, MAXQXY, xpos, ypos, occ, iscales, scales, tcsf)
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
integer              , dimension(2)       , intent(in)  :: icell    ! Actual unit cell numbers
integer              , dimension(2)       , intent(in)  :: inc      ! Number data points along axes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Total number in reciprocal space
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(in)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(in)  :: ypos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer              , dimension(2)       , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
real(kind=PREC_DP)   , dimension(2)       , intent(in)  :: scales   ! Scaling
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
real(kind=PREC_DP)   , dimension(:), allocatable :: xj  ! Coordinates
real(kind=PREC_DP)   , dimension(:), allocatable :: yj  ! Coordinates
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
!
type(finufft_opts) :: opts                  ! Use default options
!
iflag = 1
M   = nat
ms  = inc(1)
mt  = inc(2)
nj = ms*mt
allocate(xj(M))
allocate(yj(M))
allocate(cj(M))
allocate(fk(nj))
!
xj = zpi*xpos(1:nat)/real(icell(1),kind=PREC_DP)*scales(1)
yj = zpi*ypos(1:nat)/real(icell(2),kind=PREC_DP)*scales(2)
cj = cmplx(occ, 0.0D0)
tol = 1d-9
!
call finufft_default_opts(opts)
call finufft2d1(M, xj, yj, cj, iflag, tol, ms, mt, fk, opts, ier)
i = 0
do h=1,ms, iscales(1)
  hh = (h-1)/iscales(1) + 1
  do k=1, mt, iscales(2)
    kk = (k-1)/iscales(2) + 1
!   i = i+1
!   tcsf(i) = fk((k-1)*ms+h)
    tcsf(hh,kk,1) = fk((k-1)*ms+h)
  enddo
enddo
!
deallocate(xj)
deallocate(yj)
deallocate(cj)
deallocate(fk)
!
end subroutine four_strucf_2d
!
!******************************************************************************
!
subroutine four_strucf_3d(MAXATOMS, nat, icell, inc, MAXQXY, xpos, ypos, zpos, occ, iscales, scales, tcsf)
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
integer              , dimension(3)       , intent(in)  :: icell    ! Actual unit cell numbers
integer              , dimension(3)       , intent(in)  :: inc      ! Number data points along exes
integer              , dimension(3)       , intent(in)  :: MAXQXY   ! Total data points reciprocal space
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(in)  :: xpos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(in)  :: ypos     ! Actual atom coordinates
real(kind=PREC_DP)   , dimension(MAXATOMS), intent(in)  :: zpos     ! Actual atom coordinates
real(kind=PREC_DP)                        , intent(in)  :: occ      ! Atom occupancies
integer              , dimension(3)       , intent(in)  :: iscales  ! Scaling if DELTA (hkl) /= 1/cr_icc
real   (kind=PREC_DP), dimension(3)       , intent(in)  :: scales   ! Scaleing
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
real(kind=PREC_DP)   , dimension(:), allocatable :: xj  ! Coordinates
real(kind=PREC_DP)   , dimension(:), allocatable :: yj  ! Coordinates
real(kind=PREC_DP)   , dimension(:), allocatable :: zj  ! Coordinates
complex(kind=PREC_DP), dimension(:), allocatable :: cj  ! Real space amplitudes = 1*occ
complex(kind=PREC_DP), dimension(:), allocatable :: fk  ! Complex structure factor
!
type(finufft_opts) :: opts                  ! Use default options
!
iflag = 1
M   = nat
ms  = inc(1)
mt  = inc(2)
mu  = inc(3)
nj = ms*mt*mu
allocate(xj(M))
allocate(yj(M))
allocate(zj(M))
allocate(cj(M))
allocate(fk(nj))
!
xj = zpi*xpos(1:nat)/real(icell(1),kind=PREC_DP)*scales(1)
yj = zpi*ypos(1:nat)/real(icell(2),kind=PREC_DP)*scales(2)
zj = zpi*zpos(1:nat)/real(icell(3),kind=PREC_DP)*scales(3)
cj = cmplx(occ, 0.0D0)
tol = 1.0D-7
!
call finufft_default_opts(opts)
call finufft3d1(M, xj, yj, zj, cj, iflag, tol, ms, mt, mu, fk, opts, ier)
!
i = 0
do l=1, mu, iscales(3)
   ll = (l-1)/iscales(3) + 1
   do k=1, mt, iscales(2)
      kk = (k-1)/iscales(2) + 1
      do h=1,ms, iscales(1)
         hh = (h-1)/iscales(1) + 1
         i = i+1
!      tcsf(i) = fk((l-1)*ms*mt+(k-1)*ms+h)
!        tcsf(hh,kk,ll) = fk((l-1)*ms*mt+(k-1)*ms+h)
         tcsf(hh,kk,ll) = fk(i)
      enddo
   enddo
enddo
!
deallocate(xj)
deallocate(yj)
deallocate(zj)
deallocate(cj)
deallocate(fk)
!
!
end subroutine four_strucf_3d
!
!*******************************************************************************
!
!
end module four_finufft
