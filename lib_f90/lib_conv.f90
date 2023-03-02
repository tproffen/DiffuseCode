module lib_conv_mod
!-
!   Convolute two data sets via FFT
!+
!
contains
!
!###############################################################################
!
subroutine do_convolute_1d_real(idims, infield, secfield)
!-
!  Convolute two real valued 1D data sets via FFT
!+
!
use iso_c_binding
!
use lib_f90_fftw3
use map_1dtofield
use precision_mod
!
implicit none
!
integer                             , intent(in)    :: idims     ! field dimension
real(kind=PREC_DP), dimension(idims), intent(inout) :: infield   ! Will be replaced by result
real(kind=PREC_DP), dimension(idims), intent(in   ) :: secfield  ! The scond field
!
integer, dimension(3) :: num_a
integer, dimension(3) :: dsort
complex(kind=PREC_DP), dimension(:), allocatable :: in_pattern
complex(kind=PREC_DP), dimension(:), allocatable :: out_pattern
complex(kind=PREC_DP), dimension(:), allocatable :: temp
real(kind=PREC_DP) :: weight
!
type(c_ptr) :: plan    ! FFWT3 plan
!
num_a = 1
dsort = 1
num_a(1) = idims
dsort(1) = 1
dsort(2) = 2
dsort(3) = 3
!
allocate(in_pattern(1:idims))
allocate(out_pattern(1:idims))
allocate(temp(1:idims))
!
in_pattern = cmplx(infield, 0.0D0)
plan = fftw_plan_dft_1d(idims, in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call maptofftfd(num_a, dsort,secfield, in_pattern)      ! Prepare second field
call fftw_execute_dft(plan, in_pattern, temp)           ! FFT second field
!
call maptofftfd(num_a, dsort,infield, in_pattern)      ! Prepare input field
call fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT infield
!
out_pattern = temp * out_pattern                        ! Multiply the Fouriers
!
call fftw_execute_dft(plan, out_pattern,  in_pattern)   ! Back FFT the product
!
call mapfftfdtoline(num_a, dsort, infield, in_pattern)  ! Map back to original sequence
!
weight = sqrt(real(num_a(1)**3       ,kind=PREC_DP))    ! Scale for H
infield =  infield/weight
!
call   fftw_destroy_plan(plan)
deallocate(in_pattern)
deallocate(out_pattern)
deallocate(temp)
!
end subroutine do_convolute_1d_real
!
!*******************************************************************************
!
subroutine do_convolute_2d_real(idims, infield, secfield)
!-
!  Convolute two real valued 1D data sets via FFT
!+
!
use iso_c_binding
!
use lib_f90_fftw3
use map_1dtofield
use precision_mod
!
implicit none
!
integer           , dimension(2)                 , intent(in)    :: idims     ! field dimension
real(kind=PREC_DP), dimension(idims(1), idims(2)), intent(inout) :: infield   ! Will be replaced by result
real(kind=PREC_DP), dimension(idims(1), idims(2)), intent(in   ) :: secfield  ! The scond field
!
integer, dimension(3) :: num_a
integer, dimension(3) :: dsort
complex(kind=PREC_DP), dimension(:,:), allocatable :: in_pattern
complex(kind=PREC_DP), dimension(:,:), allocatable :: out_pattern
complex(kind=PREC_DP), dimension(:,:), allocatable :: temp
real(kind=PREC_DP) :: weight
!
type(c_ptr) :: plan    ! FFWT3 plan
!
num_a = 1
dsort = 1
num_a(1) = idims(1)
num_a(2) = idims(2)
dsort(1) = 1
dsort(2) = 2
dsort(3) = 3
!
allocate(in_pattern (1:idims(1), 1:idims(2)))
allocate(out_pattern(1:idims(1), 1:idims(2)))
allocate(temp       (1:idims(1), 1:idims(2)))
!
in_pattern = cmplx(infield, 0.0D0)
plan = fftw_plan_dft_2d(idims(1), idims(2), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call maptofftfd(num_a, dsort,secfield, in_pattern)      ! Prepare second field
!write(*,*) ' CONV2D sec ', sum(secfield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, temp)           ! FFT second field
!
call maptofftfd(num_a, dsort,infield, in_pattern)      ! Prepare input field
!write(*,*) ' CONV2D in  ', sum( infield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT infield
!
out_pattern = temp * out_pattern                        ! Multiply the Fouriers
!
call fftw_execute_dft(plan, out_pattern,  in_pattern)   ! Back FFT the product
call mapfftfdtoline(num_a, dsort, infield, in_pattern)  ! Map back to original sequence
!write(*,*) ' CONV2D res ', sum( infield), sum(in_pattern)
!write(*,*) ' idims ', idims , idims**3
!
weight = sqrt(real(idims(1)**3,kind=PREC_DP))*sqrt(real(idims(2)**3,kind=PREC_DP))
infield = infield/weight
!write(*,*) ' CONV2D fin ', sum( infield),  weight
!
call   fftw_destroy_plan(plan)
deallocate(in_pattern)
deallocate(out_pattern)
deallocate(temp)
!
end subroutine do_convolute_2d_real
!
!*******************************************************************************
!
subroutine do_convolute_3d_real(idims, infield, secfield)
!-
!  Convolute two real valued 3D data sets via FFT
!+
!
use iso_c_binding
!
use lib_f90_fftw3
use map_1dtofield
use precision_mod
!
implicit none
!
integer           , dimension(3)                           , intent(in)    :: idims     ! field dimension
real(kind=PREC_DP), dimension(idims(1), idims(2), idims(3)), intent(inout) :: infield   ! Will be replaced by result
real(kind=PREC_DP), dimension(idims(1), idims(2), idims(3)), intent(in   ) :: secfield  ! The scond field
!
integer, dimension(3) :: num_a
integer, dimension(3) :: dsort
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: in_pattern
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: out_pattern
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: temp
real(kind=PREC_DP) :: weight
!
type(c_ptr) :: plan    ! FFWT3 plan
!
num_a = 1
dsort = 1
num_a(1) = idims(1)
num_a(2) = idims(2)
num_a(3) = idims(3)
dsort(1) = 1
dsort(2) = 2
dsort(3) = 3
!
allocate(in_pattern (1:idims(1), 1:idims(2), 1:idims(3)))
allocate(out_pattern(1:idims(1), 1:idims(2), 1:idims(3)))
allocate(temp       (1:idims(1), 1:idims(2), 1:idims(3)))
!
in_pattern = cmplx(infield, 0.0D0)
plan = fftw_plan_dft_3d(idims(1), idims(2), idims(3), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call maptofftfd(num_a, dsort,secfield, in_pattern)      ! Prepare second field
!write(*,*) ' CONV3D sec ', sum(secfield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, temp)           ! FFT second field
!
call maptofftfd(num_a, dsort,infield, in_pattern)      ! Prepare input field
!write(*,*) ' CONV3D in  ', sum( infield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT infield
!
out_pattern = temp * out_pattern                        ! Multiply the Fouriers
!
call fftw_execute_dft(plan, out_pattern,  in_pattern)   ! Back FFT the product
call mapfftfdtoline(num_a, dsort, infield, in_pattern)  ! Map back to original sequence
!write(*,*) ' CONV3D res ', sum( infield), sum(in_pattern)
!write(*,*) ' idims ', idims , idims**3
!
weight = sqrt(real(idims(1)**3,kind=PREC_DP))*sqrt(real(idims(2)**3,kind=PREC_DP)) * &
         sqrt(real(idims(3)**3,kind=PREC_DP))
infield = infield/weight
!write(*,*) ' CONV3D fin ', sum( infield),  weight
!
call   fftw_destroy_plan(plan)
deallocate(in_pattern)
deallocate(out_pattern)
deallocate(temp)
!
end subroutine do_convolute_3d_real
!
!###############################################################################
!
end module lib_conv_mod
