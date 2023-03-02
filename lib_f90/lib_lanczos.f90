module do_lanczos_mod
!-
!  Lanczos filter
!
!  Apply a Lanczos filter to the infield
!  The outfield if "nscale" times larger in all dimensions
!  The +-m data points of the infield are taken into account.
!  rcut allows to adjust the dampening of the infield values away from the central point
!+
!
interface do_lanczos
  module procedure do_lanczos_3d_cart, do_lanczos_2d_cart, do_lanczos_1d_cart,  &
                   do_lanczos_3d_fft,  do_lanczos_2d_fft,  do_lanczos_1d_fft
end interface do_lanczos
!
private
public do_lanczos
!
contains
!
!###############################################################################
!
subroutine do_lanczos_3d_cart(nscale, rcut, m, inc, infield, idims, outfield)
!-
!  3-dimensonal cartesian version
!+
use precision_mod
!
!
implicit none
!
real(kind=PREC_DP)                                       , intent(in )   :: nscale ! Outfield is times larger
real(kind=PREC_DP)                                       , intent(in )   :: rcut
integer                                                  , intent(in )   :: m
integer           , dimension(3)                         , intent(in )   :: inc  !Input dimensions
!real(kind=PREC_DP), dimension(4)                         , intent(in )   :: limits
real(kind=PREC_DP), dimension(inc(1), inc(2), inc(3))    , intent(in )   ::  infield
integer           , dimension(3)                         , intent(in )   :: idims  !Input dimensions
real(kind=PREC_DP), dimension(idims(1), idims(2), idims(3)), intent(out) :: outfield
!
integer :: i1, i2, i3 !,j,k
integer :: hl, hh   ! lower, higher limit
integer :: kl, kh   ! lower, higher limit
integer :: ll, lh   ! lower, higher limit
integer :: h, k, l
integer :: qg1, qg2, qg3
real(kind=PREC_DP), dimension(3)  :: q  ! Q-vector
real(kind=PREC_DP), dimension(3)  :: g  ! G-vector
real(kind=PREC_DP) :: summe
real(kind=PREC_DP) :: factor
real(kind=PREC_DP), dimension(:,:,:), allocatable :: sinc_table
!
outfield = 0.0D0
!
ll = nint(-nscale*m-1.0D0)
lh = nint( nscale*m+1.0D0)
kl = nint(-nscale*m-1.0D0)
kh = nint( nscale*m+1.0D0)
hl = nint(-nscale*m-1.0D0)
hh = nint( nscale*m+1.0D0)
allocate(sinc_table(hl:hh, kl:kh, ll:lh))
q  = 0.0D0
do i3=kl, kh
   g(3) = real(i3,kind=PREC_DP)/nscale
do i2=kl, kh
   g(2) = real(i2,kind=PREC_DP)/nscale
   do i1=hl, hh
      g(1) = real(i1,kind=PREC_DP)/nscale
      sinc_table(i1,i2,i3) = wsinc_3d(q, g, rcut)*lqgm_3d(q,g,m)
   enddo
enddo
enddo
!
do i3 = 1, idims(3)
!write(*,*) ' AT ', i3, idims(3)
   ll = max(1       , nint((i3-1)/nscale-m + 1))
   lh = min( inc (2), nint((i3-1)/nscale+m + 1))
   q(3) = real((i3-1),kind=PREC_DP)/nscale + 1
do i2 = 1, idims(2)
   kl = max(1       , nint((i2-1)/nscale-m + 1))
   kh = min( inc (2), nint((i2-1)/nscale+m + 1))
   q(2) = real((i2-1),kind=PREC_DP)/nscale + 1
   do i1 = 1, idims(1)      ! Loop over all interpolated data points
      hl = max(1       , nint((i1-1)/nscale-m + 1))
      hh = min( inc (1), nint((i1-1)/nscale+m + 1))
      q(1) = real((i1-1),kind=PREC_DP)/nscale + 1
      summe = 0.0D0
      do l=ll, lh
         g(3) = real(l, kind=PREC_DP)
         qg3 = nint((q(3)-g(3))*nscale)
         do k=kl, kh
            g(2) = real(k, kind=PREC_DP)
            qg2 = nint((q(2)-g(2))*nscale)
            do h=hl, hh
               g(1) = real(h, kind=PREC_DP)
               qg1 = nint((q(1)-g(1))*nscale)
               factor = sinc_table(qg1, qg2, qg3)
!              factor = wsinc_3d(q, g, rcut)*lqgm_3D(q,g,m)
               summe  = summe + factor
               outfield(i1,i2, i3) = outfield(i1,i2, i3) + infield(h,k,l)* factor
            enddo
         enddo
      enddo
      if(summe>0) then
         outfield(i1,i2, i3) = outfield(i1,i2, i3)/summe
      endif
   enddo
enddo
enddo
!
deallocate(sinc_table)
!
end subroutine do_lanczos_3d_cart
!
!*******************************************************************************
!
subroutine do_lanczos_2d_cart(nscale, rcut, m, inc, infield, idims, outfield)
!-
!  2-dimensonal cartesian version
!+
use precision_mod
!
!
implicit none
!
real(kind=PREC_DP)                                       , intent(in )   :: nscale ! Outfield is times larger
real(kind=PREC_DP)                                       , intent(in )   :: rcut
integer                                                  , intent(in )   :: m
integer           , dimension(2)                         , intent(in )   :: inc  !Input dimensions
!real(kind=PREC_DP), dimension(4)                         , intent(in )   :: limits
real(kind=PREC_DP), dimension(inc(1), inc(2))            , intent(in )   ::  infield
integer           , dimension(2)                         , intent(in )   :: idims  !Input dimensions
real(kind=PREC_DP), dimension(idims(1), idims(2))        , intent(out) :: outfield
!
integer :: i1, i2 !,j,k
integer :: hl, hh   ! lower, higher limit
integer :: kl, kh   ! lower, higher limit
integer :: h, k
integer :: qg1, qg2
real(kind=PREC_DP), dimension(2)  :: q  ! Q-vector
real(kind=PREC_DP), dimension(2)  :: g  ! G-vector
real(kind=PREC_DP) :: summe
real(kind=PREC_DP) :: factor
real(kind=PREC_DP), dimension(:,:), allocatable :: sinc_table
!
outfield = 0.0D0
!
kl = nint(-nscale*m-1.0D0)
kh = nint( nscale*m+1.0D0)
hl = nint(-nscale*m-1.0D0)
hh = nint( nscale*m+1.0D0)
allocate(sinc_table(hl:hh, kl:kh))
q  = 0.0D0
do i2=kl, kh
   g(2) = real(i2,kind=PREC_DP)/nscale
   do i1=hl, hh
      g(1) = real(i1,kind=PREC_DP)/nscale
      sinc_table(i1,i2) = wsinc_2d(q, g, rcut)*lqgm_2d(q,g,m)
   enddo
enddo
!
!
do i2 = 1, idims(2)
   kl = max(1       , nint((i2-1)/nscale-m + 1   ))
   kh = min( inc (2), nint((i2-1)/nscale+m + 1   ))
   q(2) = real((i2-1),kind=PREC_DP)/nscale + 1
   do i1 = 1, idims(1)      ! Loop over all interpolated data points
      hl = max(1       , nint(i1/nscale-m    ))
      hh = min( inc (1), nint(i1/nscale+m    ))
      q(1) = real((i1-1),kind=PREC_DP)/nscale + 1
      summe = 0.0D0
      do k=kl, kh
         g(2) = real(k, kind=PREC_DP)
         qg2 = nint((q(2)-g(2))*nscale)
         do h=hl, hh
            g(1) = real(h, kind=PREC_DP)
            qg1 = nint((q(1)-g(1))*nscale)
            factor = sinc_table(qg1, qg2)
!           factor = wsinc_2d(q, g, rcut)*lqgm_2d(q,g,m)
            summe  = summe + factor
           outfield(i1,i2) = outfield(i1,i2) + infield(h,k)* factor
         enddo
      enddo
      if(summe>0) then
         outfield(i1,i2) = outfield(i1,i2)/summe
      endif
   enddo
enddo
!
deallocate(sinc_table)
!
end subroutine do_lanczos_2d_cart
!
!*******************************************************************************
!
subroutine do_lanczos_1d_cart(nscale, rcut, m, inc, infield, idims, outfield)
!-
!  2-dimensonal cartesian version
!+
use precision_mod

!
!
implicit none
!
real(kind=PREC_DP)                                       , intent(in )   :: nscale ! Outfield is times larger
real(kind=PREC_DP)                                       , intent(in )   :: rcut
integer                                                  , intent(in )   :: m
integer                                                  , intent(in )   :: inc  !Input dimensions
!real(kind=PREC_DP), dimension(4)                         , intent(in )   :: limits
real(kind=PREC_DP), dimension(inc)                    , intent(in )   ::  infield
integer                                                  , intent(in )   :: idims  !Input dimensions
real(kind=PREC_DP), dimension(idims)                  , intent(out) :: outfield
!
integer :: i1 !,j,k
integer :: hl, hh   ! lower, higher limit
integer :: h
integer :: qg1
real(kind=PREC_DP) :: q  ! Q-vector
real(kind=PREC_DP) :: g  ! G-vector
real(kind=PREC_DP) :: summe
real(kind=PREC_DP) :: factor
real(kind=PREC_DP), dimension(:), allocatable :: sinc_table
!
outfield = 0.0D0
q = 0.0D0
hl = nint(-nscale*m-1.0D0)
hh = nint( nscale*m+1.0D0)
allocate(sinc_table(hl:hh))
do i1=hl, hh
   g = real(i1,kind=PREC_DP)/nscale
   sinc_table(i1) = wsinc_1d(q, g, rcut)*lqgm_1d(q,g,m)
enddo
!
do i1 = 1, idims         ! Loop over all interpolated data points
   hl = max(1       , nint((i1-1)/nscale-m +1   ))
   hh = min( inc    , nint((i1-1)/nscale+m +1   ))
   q    = real((i1-1),kind=PREC_DP)/nscale + 1
   summe = 0.0D0
   do h=hl, hh
      g    = real(h, kind=PREC_DP)
      qg1 = nint((q-g)*nscale)
      factor = sinc_table(qg1)
      summe  = summe + factor
      outfield(i1) = outfield(i1) + infield(h)* factor
   enddo
   if(summe>0) then
      outfield(i1) = outfield(i1)/summe
   endif
enddo
deallocate(sinc_table)
!
end subroutine do_lanczos_1d_cart
!
!*******************************************************************************
!
subroutine do_lanczos_3d_fft(nscale, rcut, m, inc, infield, idims, outfield, flag)
!-
!  2-dimensonal cartesian version
!+
!
use lib_conv_mod
use precision_mod
!
!
implicit none
!
real(kind=PREC_DP)                                       , intent(in )   :: nscale ! Outfield is times larger
real(kind=PREC_DP)                                       , intent(in )   :: rcut
integer                                                  , intent(in )   :: m
integer           , dimension(3)                         , intent(in )   :: inc  !Input dimensions
real(kind=PREC_DP), dimension(inc(1), inc(2), inc(3))    , intent(in )   ::  infield
integer           , dimension(3)                         , intent(in )   :: idims  !Input dimensions
real(kind=PREC_DP), dimension(idims(1), idims(2), idims(3)), intent(out) :: outfield
logical , intent(in) :: flag
!
integer :: i1, i2, i3 !,j,k
integer :: hl, hh   ! lower, higher limit
integer :: kl, kh   ! lower, higher limit
integer :: ll, lh   ! lower, higher limit
real(kind=PREC_DP), dimension(3)  :: q  ! Q-vector
real(kind=PREC_DP), dimension(3)  :: g  ! G-vector
real(kind=PREC_DP) :: factor
real(kind=PREC_DP), dimension(:,:,:), allocatable :: sinc_table
!
allocate(sinc_table (1:idims(1), 1:idims(2), 1:idims(3)))
!
outfield    = 0.0D0
sinc_table  = 0.0D0
!
i1=(idims(1)+1)/2
i2=(idims(2)+1)/2
i3=(idims(3)+1)/2
!
hl = nint(-nscale*m-1.0D0) + i1
hh = nint( nscale*m+1.0D0) + i1
kl = nint(-nscale*m-1.0D0) + i2
kh = nint( nscale*m+1.0D0) + i2
ll = nint(-nscale*m-1.0D0) + i3
lh = nint( nscale*m+1.0D0) + i3
!
q(1) = i1/nscale
q(2) = i2/nscale
q(3) = i3/nscale
!!write(*,*) 
!write(*,*) 
do i3=ll, lh
   g(3) = real(i3,kind=PREC_DP)/nscale
   do i2=kl, kh
      g(2) = real(i2,kind=PREC_DP)/nscale
      do i1=hl, hh
         g(1) = real(i1,kind=PREC_DP)/nscale
         sinc_table(i1,i2,i3) = wsinc_2d(q, g, rcut)*lqgm_2d(q,g,m)
      enddo
   enddo
enddo
!write(*,*) ' SINC TABLE ', sum(sinc_table)
!
outfield(1:idims(1):nint(nscale), 1:idims(2):nint(nscale), 1:idims(3):nint(nscale)) = &
      infield(1:inc(1), 1:inc(2), 1:inc(3))
call do_convolute_3d_real(idims, outfield, sinc_table)
!
factor = sum( infield)/sum(outfield)*nscale
!write(*,*) sum( infield) , sum(outfield), sum(outfield)/nscale, factor
outfield = outfield * factor
!
deallocate(sinc_table)
!
end subroutine do_lanczos_3d_fft
!
!*******************************************************************************
!
subroutine do_lanczos_2d_fft(nscale, rcut, m, inc, infield, idims, outfield, flag)
!-
!  2-dimensonal cartesian FFT version
!+
!
use lib_conv_mod
use precision_mod
!
!
implicit none
!
real(kind=PREC_DP)                                       , intent(in )   :: nscale ! Outfield is times larger
real(kind=PREC_DP)                                       , intent(in )   :: rcut
integer                                                  , intent(in )   :: m
integer           , dimension(2)                         , intent(in )   :: inc  !Input dimensions
real(kind=PREC_DP), dimension(inc(1), inc(2))            , intent(in )   ::  infield
integer           , dimension(2)                         , intent(in )   :: idims  !Input dimensions
real(kind=PREC_DP), dimension(idims(1), idims(2))        , intent(out) :: outfield
logical , intent(in) :: flag
!
integer :: i1, i2 !,j,k
integer :: hl, hh   ! lower, higher limit
integer :: kl, kh   ! lower, higher limit
real(kind=PREC_DP), dimension(2)  :: q  ! Q-vector
real(kind=PREC_DP), dimension(2)  :: g  ! G-vector
real(kind=PREC_DP) :: factor
real(kind=PREC_DP), dimension(:,:), allocatable :: sinc_table
!
allocate(sinc_table (1:idims(1), 1:idims(2)))
!
outfield    = 0.0D0
sinc_table  = 0.0D0
!
i1=(idims(1)+1)/2
i2=(idims(2)+1)/2
!
hl = nint(-nscale*m-1.0D0) + i1
hh = nint( nscale*m+1.0D0) + i1
kl = nint(-nscale*m-1.0D0) + i2
kh = nint( nscale*m+1.0D0) + i2
!
q(1) = i1/nscale
q(2) = i2/nscale
!write(*,*) 
!write(*,*) 
do i2=kl, kh
   g(2) = real(i2,kind=PREC_DP)/nscale
   do i1=hl, hh
      g(1) = real(i1,kind=PREC_DP)/nscale
      sinc_table(i1,i2) = wsinc_2d(q, g, rcut)*lqgm_2d(q,g,m)
   enddo
enddo
!
outfield(1:idims(1):nint(nscale), 1:idims(2):nint(nscale)) = infield(1:inc(1), 1:inc(2))
call do_convolute_2d_real(idims, outfield, sinc_table)
!
factor = sum( infield)/sum(outfield)*nscale
outfield = outfield * factor
!
deallocate(sinc_table)
!
end subroutine do_lanczos_2d_fft
!
!*******************************************************************************
!
subroutine do_lanczos_1d_fft (nscale, rcut, m, inc, infield, idims, outfield, flag)
!-
!  1-dimensonal cartesian FFT version
!+
!
use lib_conv_mod
use lib_write_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP)                                       , intent(in )   :: nscale ! Outfield is times larger
real(kind=PREC_DP)                                       , intent(in )   :: rcut
integer                                                  , intent(in )   :: m
integer                                                  , intent(in )   :: inc  !Input dimensions
!real(kind=PREC_DP), dimension(4)                         , intent(in )   :: limits
real(kind=PREC_DP), dimension(inc)                    , intent(in )   ::  infield
integer                                                  , intent(in )   :: idims  !Input dimensions
real(kind=PREC_DP), dimension(idims)                  , intent(out) :: outfield
logical , intent(in) :: flag
!
integer :: i1 !,j,k
integer :: hl, hh   ! lower, higher limit
!integer :: h
real(kind=PREC_DP) :: q  ! Q-vector
real(kind=PREC_DP) :: g  ! G-vector
!real(kind=PREC_DP) :: summe
real(kind=PREC_DP) :: factor
!real(kind=PREC_DP), dimension(:), allocatable :: infield_aug
real(kind=PREC_DP), dimension(:), allocatable :: sinc_table
!
!allocate(infield_aug(1:idims))
allocate(sinc_table (1:idims))
!
outfield    = 0.0D0
outfield    = 0.0D0
sinc_table  = 0.0D0
!
i1=(idims+1)/2
hl = nint(-nscale*m-1.0D0) + i1
hh = nint( nscale*m+1.0D0) + i1
q = i1/nscale
!
do i1=hl, hh
   g = real(i1,kind=PREC_DP)/nscale
   sinc_table(i1) = wsinc_1d(q, g, rcut)*lqgm_1d(q,g,m)
enddo
!
outfield      (1:idims:nint(nscale)) = infield(1:inc)
!
call do_convolute_1d_real(idims, outfield, sinc_table)
!
factor = sum( infield)/sum(outfield)*nscale
!write(*,*) 'INTEGRAL ', sum(infield), sum(outfield), sum(outfield), factor
outfield = outfield *factor
call tofile(idims, 'convolute.data', outfield, -5.0D0, (10.0D0/(idims-1)) )
!
deallocate(sinc_table)
!deallocate(infield_aug)
!
end subroutine do_lanczos_1d_fft
!
!*******************************************************************************
!
function wsinc_3d(q, g, rcut)  result(val)
!
use lib_functions_mod
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: val
real(kind=PREC_DP), dimension(3), intent(in) :: q
real(kind=PREC_DP), dimension(3), intent(in) :: g
real(kind=PREC_DP) :: rcut
real(kind=PREC_DP), dimension(3) :: x
!
val = 1.0D0
!
x   = PI*(q-g)*rcut
val = sinc(x(1)) * sinc(x(2)) * sinc(x(3))
!
end function wsinc_3d
!
!*******************************************************************************
!
function wsinc_2d(q, g, rcut)  result(val)
!
use lib_functions_mod
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: val
real(kind=PREC_DP), dimension(2), intent(in) :: q
real(kind=PREC_DP), dimension(2), intent(in) :: g
real(kind=PREC_DP) :: rcut
real(kind=PREC_DP), dimension(2) :: x
!
val = 0.0D0
!
x   = PI*(q-g)*rcut
val = sinc(x(1)) * sinc(x(2))
!
end function wsinc_2d
!
!*******************************************************************************
!
function wsinc_1d(q, g, rcut)  result(val)
!
use lib_functions_mod
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: val
real(kind=PREC_DP), intent(in) :: q
real(kind=PREC_DP), intent(in) :: g
real(kind=PREC_DP) :: rcut
real(kind=PREC_DP) :: x
!
val = 0.0D0
!
x   = PI*(q-g)*rcut
val = sinc(x) * sinc(x)
!
end function wsinc_1d
!
!*******************************************************************************
!
function lqgm_3d(q, g, m)  result(val)
!
use lib_functions_mod
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: val
real(kind=PREC_DP), dimension(3), intent(in) :: q
real(kind=PREC_DP), dimension(3), intent(in) :: g
integer            :: m
real(kind=PREC_DP), dimension(3) :: x
!
val = 1.0D0
!
if(all(abs((q-g))<real(m,kind=PREC_DP))) then
   x   = PI*(q-g)/real(m, kind=PREC_DP)
   val = sinc(x(1)) * sinc(x(2)) * sinc(x(3))
else
  val = 0.0D0
endif
!
end function lqgm_3d
!
!*******************************************************************************
!
function lqgm_2d(q2, g, m)  result(val)
!
use lib_functions_mod
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: val
real(kind=PREC_DP), dimension(2), intent(in) :: q2
real(kind=PREC_DP), dimension(2), intent(in) :: g
integer            :: m
real(kind=PREC_DP), dimension(2) :: x
!
val = 0.0D0
!
if(all(abs((q2-g))<real(m,kind=PREC_DP))) then
   x   = PI*(q2-g)/real(m, kind=PREC_DP)
   val = sinc(x(1)) * sinc(x(2))
else
  val = 0.0D0
endif
!
end function lqgm_2d
!
!*******************************************************************************
!
function lqgm_1d(q2, g, m)  result(val)
!
use lib_functions_mod
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP) :: val
real(kind=PREC_DP), intent(in) :: q2
real(kind=PREC_DP), intent(in) :: g
integer            :: m
real(kind=PREC_DP) :: x
!
val = 0.0D0
!
if((abs((q2-g))<real(m,kind=PREC_DP))) then
   x   = PI*(q2-g)/real(m, kind=PREC_DP)
   val = sinc(x)
else
  val = 0.0D0
endif
!
end function lqgm_1d
!
!###############################################################################
!
end module do_lanczos_mod
