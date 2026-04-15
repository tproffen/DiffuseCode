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
!use lib_write_mod
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
!real(kind=PREC_DP)   :: xmin
!real(kind=PREC_DP)   :: xstep
!
type(c_ptr) :: plan    ! FFWT3 plan
!
!write(*,*) ' IN do_convolute_1d_real'
!xmin  = 0
!xstep = 1
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
!call tofile(idims, 'in_field', infield,  xmin, xstep)
!call tofile(idims, 'sec_field',  secfield, xmin, xstep)
in_pattern = cmplx(infield, 0.0D0, kind=PREC_DP)
plan = fftw_plan_dft_1d(idims, in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call maptofftfd(num_a, dsort,secfield, in_pattern)      ! Prepare second field
!call tofile(idims, 'sec_field_map',  in_pattern, xmin, xstep)
call fftw_execute_dft(plan, in_pattern, temp)           ! FFT second field
!call tofile(idims, 'sec_field_fft',  temp      , xmin, xstep)
!
call maptofftfd(num_a, dsort,infield, in_pattern)      ! Prepare input field
!call tofile(idims, 'in_field_map', in_pattern,  xmin, xstep)
call fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT infield
!call tofile(idims, 'in_field_fft',out_pattern,  xmin, xstep)
!
out_pattern = temp * out_pattern                        ! Multiply the Fouriers
!call tofile(idims, 'out_field_fft',out_pattern,  xmin, xstep)
!
call fftw_execute_dft(plan, out_pattern,  in_pattern)   ! Back FFT the product
!call tofile(idims, 'out_field_bck', in_pattern,  xmin, xstep)
!
call mapfftfdtoline(num_a, dsort, infield, in_pattern)  ! Map back to original sequence
!call tofile(idims, 'in_field_final', infield,  xmin, xstep)
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
!  Convolute two real valued 2D data sets via FFT
!+
!
use iso_c_binding
!
use lib_f90_fftw3
use map_1dtofield
use precision_mod
!
!use lib_write_mod
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
!real(kind=PREC_DP)   , dimension(3) :: xmin
!real(kind=PREC_DP)   , dimension(3) :: xstep
!xmin = 0.0D0
!xstep = 1.0D0
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
!call tofile(idims(1:2), 'in_field', infield,  xmin, xstep)
!call tofile(idims(1:2), 'sec_field',  secfield, xmin, xstep)
in_pattern = cmplx(infield, 0.0D0, kind=PREC_DP)
plan = fftw_plan_dft_2d(idims(2), idims(1), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call maptofftfd(num_a, dsort,secfield, in_pattern)      ! Prepare second field
!call tofile(idims(1:2), 'sec_field_map', in_pattern, xmin, xstep)
!write(*,*) ' CONV2D sec ', sum(secfield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, temp)           ! FFT second field
!call tofile(idims, 'sec_field_fft', temp    , xmin, xstep)
!
call maptofftfd(num_a, dsort,infield, in_pattern)      ! Prepare input field
!call tofile(idims(1:2), 'in_pattern_map', in_pattern, xmin, xstep)
!write(*,*) ' CONV2D in  ', sum( infield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT infield
!call tofile(idims(1:2), 'in_pattern_fft',out_pattern, xmin, xstep)
!
out_pattern = temp * out_pattern                        ! Multiply the Fouriers
!call tofile(idims(1:2), 'out_pattern_fft',out_pattern, xmin, xstep)
!
call fftw_execute_dft(plan, out_pattern,  in_pattern)   ! Back FFT the product
!call tofile(idims(1:2), 'out_pattern_bck', in_pattern, xmin, xstep)
call mapfftfdtoline(num_a, dsort, infield, in_pattern)  ! Map back to original sequence
!call tofile(idims(1:2), 'in_field_final', infield, xmin, xstep)
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
!use lib_write_mod
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
!complex(kind=PREC_DP), dimension(:,:,:), allocatable :: temp_pattern
complex(kind=PREC_DP), dimension(:,:,:), allocatable :: temp_sec
real(kind=PREC_DP) :: weight
!real(kind=PREC_DP)   , dimension(3) :: xmin
!real(kind=PREC_DP)   , dimension(3) :: xstep
!
type(c_ptr) :: plan    ! FFWT3 plan
!
!xmin   = 0.0D0
!xstep  = 1.0D0
num_a = 1
dsort = 1
num_a(1) = idims(1)
num_a(2) = idims(2)
num_a(3) = idims(3)
dsort(1) = 1
dsort(2) = 2
dsort(3) = 3
!
allocate(in_pattern  (1:idims(1), 1:idims(2), 1:idims(3)))
allocate(out_pattern (1:idims(1), 1:idims(2), 1:idims(3)))
allocate(temp_sec    (1:idims(1), 1:idims(2), 1:idims(3)))
!allocate(temp_pattern(1:idims(1), 1:idims(2), 1:idims(3)))
!
in_pattern = cmplx(infield, 0.0D0, kind=PREC_DP)
!
!call tofile(idims, 'in_field', infield,  xmin, xstep)
!call tofile(idims, 'sec_field',  secfield, xmin, xstep)
!write(*,'(4(2x,f18.2))') real(in_pattern(      1 ,     1 ,       1 )), real(in_pattern(      1 , idims(2),       1 )), &
!                         real(in_pattern(idims(1),     1 ,       1 )), real(in_pattern(idims(1), idims(2),       1 ))
!write(*,'(4(2x,f18.2))') real(in_pattern(      1 ,     1 , idims(3))), real(in_pattern(      1 , idims(2), idims(3))), &
!                         real(in_pattern(idims(1),     1 , idims(3))), real(in_pattern(idims(1), idims(2), idims(3)))
plan = fftw_plan_dft_3d(idims(3), idims(2), idims(1), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call maptofftfd(num_a, dsort,secfield, in_pattern)      ! Prepare second field
!call tofile(idims, 'sec_field_map', in_pattern, xmin, xstep)
!write(*,*) ' CONV3D sec ', sum(secfield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, temp_sec)       ! FFT second field
!call tofile(idims, 'sec_field_fft', temp_sec, xmin, xstep)
!
call maptofftfd(num_a, dsort,infield, in_pattern)      ! Prepare input field
!call tofile(idims, 'in_pattern_map', in_pattern, xmin, xstep)
!!write(*,*) ' CONV3D in  ', sum( infield), sum(in_pattern)
call fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT infield
!call tofile(idims, 'in_pattern_fft',out_pattern, xmin, xstep)
!
out_pattern = temp_sec * out_pattern                    ! Multiply the Fouriers
!call tofile(idims, 'out_pattern_fft',out_pattern, xmin, xstep)
!
call fftw_execute_dft(plan, out_pattern,  in_pattern)   ! Back FFT the product
!call tofile(idims, 'out_pattern_bck', in_pattern, xmin, xstep)
call mapfftfdtoline(num_a, dsort, infield, in_pattern)  ! Map back to original sequence
!call tofile(idims, 'in_field_final', infield, xmin, xstep)
!write(*,*) ' CONV3D res ', sum( infield), sum(in_pattern)
!write(*,*) ' idims ', idims , idims**3
!
weight = sqrt(real(idims(1)**3,kind=PREC_DP))*sqrt(real(idims(2)**3,kind=PREC_DP)) * &
         sqrt(real(idims(3)**3,kind=PREC_DP))
infield = infield/weight
!write(*,'(4(2x,f18.2))') infield(      1 ,     1 ,       1 ), infield(      1 , idims(2),       1 ), &
!                         infield(idims(1),     1 ,       1 ), infield(idims(1), idims(2),       1 )
!write(*,'(4(2x,f18.2))') infield(      1 ,     1 , idims(3)), infield(      1 , idims(2), idims(3)), &
!                         infield(idims(1),     1 , idims(3)), infield(idims(1), idims(2), idims(3))
!write(*,*) ' CONV3D fin ', sum( infield),  weight
!
call   fftw_destroy_plan(plan)
deallocate(in_pattern)
deallocate(out_pattern)
deallocate(temp_sec)
!
end subroutine do_convolute_3d_real
!
!###############################################################################
!
end module lib_conv_mod
