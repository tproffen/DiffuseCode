! benchmark.f90
program benchmark
    use, intrinsic :: omp_lib
    implicit none
    integer, parameter :: dp       = kind(0.0d0)
    integer, parameter :: n        = 2000
    integer            :: nthreads = 4

integer :: MAXTHREADS = 12
    integer       :: i, j, k,l
    real(kind=dp) :: t1, t2, ts
!   real(kind=dp) :: a(n, n), b(n, n), c(n, n)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: a, b, c

!$OMP PARALLEL
!$ IF(OMP_GET_THREAD_NUM()==0) THEN
!$ MAXTHREADS  = OMP_GET_NUM_THREADS()
!$ ENDIF
!$OMP END PARALLEL
!
allocate(a(n,n))
allocate(b(n,n))
allocate(c(n,n))
    ! Set number of threads to use.
!    call omp_set_num_threads(nthreads)

    ! Initialise the PRNG and fill matrices A and B with random numbers.
    call random_seed()
    call random_number(a)
    call random_number(b)

    c = 0.0_dp
    t1 = omp_get_wtime()

    ! Calculate C = A * B sequentially.
    do j = 1, n
        do k = 1, n
            do i = 1, n
                c(i, j) = c(i, j) + a(i, k) * b(k, j)
            end do
        end do
    end do

    t2 = omp_get_wtime()
    print '(a, i3, f7.3, a)', 'single: ', 1, t2 - t1, ' s'
ts = t2 - t1


!!$ MAXTHREAD = OMP_GET_NUM_THREADS()
write(*,*) ' MAXTHREADS ', MAXTHREADS
do l = 1, MAXTHREADS
    c = 0.0_dp
    t1 = omp_get_wtime()
    call omp_set_num_threads(l)
    ! Calculate C = A * B in parallel with OpenMP, using static scheduling.
    !$omp parallel shared(a, b, c) private(i, j, k)
    !$omp do schedule(static)
    do j = 1, n
        do k = 1, n
            do i = 1, n
                c(i, j) = c(i, j) + a(i, k) * b(k, j)
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    t2 = omp_get_wtime()
    print '(a, i3, f7.3, a, f7.3)', 'OpenMP: ',l,  t2 - t1, ' s', ts/(t2-t1)
enddo
end program benchmark
