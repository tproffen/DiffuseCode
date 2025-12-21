MODULE refine_run_mod
!
USE refine_mac_mod
use precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: refine_calc     ! Calculated (ix, iy)
REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: refine_temp     ! Calculated (ix, iy) for derivs
REAL(kind=PREC_DP), DIMENSION(:,:,:, :), ALLOCATABLE :: refine_derivs   ! Derivativs (ix, iy, parameter)
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_run(line, length)
!
! Main Routine to perform the fit
!
USE refine_control_mod
USE refine_data_mod
USE refine_fit_erg
use refine_log_mod
USE refine_params_mod
USE refine_show_mod
!
USE doact_mod
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE prompt_mod
USE take_param_mod
!
USE global_data_mod
use mpi_slave_mod
use diffev_random
!
use run_mpi_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line          ! Input command line
INTEGER         , INTENT(INOUT) :: length        ! length of input command line
!
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara    ! Parameter strings
INTEGER            , DIMENSION(MAXW) :: lpara    ! length of each parameter strign
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte    ! Parameter values
!
INTEGER            , DIMENSION(4   ) :: dimen    ! Dimension of global array
INTEGER                              :: i        ! Dummy loop parameter
INTEGER                              :: ndata    ! number of data points
INTEGER                              :: ianz     ! number of parameters
LOGICAL                              :: lexist   ! File exists yes/no
!CHARACTER(LEN=MAX(PREC_STRING,LEN(line)))                  :: plmac    ! optional plot macro
LOGICAL                              :: ref_do_plot ! Do plot yes/no
LOGICAL                              :: linit       ! Initialize mrq
integer                              :: diffev_l_get_random_state  ! Copy of current random state
logical                              :: l_old_plot_status
!
character(len=PREC_STRING) :: string
integer :: ier_cmd
integer :: exit_msg
character(len=PREC_STRING) :: message
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: OPLOT     = 1
INTEGER, PARAMETER :: OINIT     = 2
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate
!
DATA oname  / 'plot' , 'init' /
DATA loname /  4     ,  4     /
opara  =  (/ '        ', 'no      '  /)
lopara =  (/ 8         ,  8          /)
owerte =  (/  0.000000 ,   0.0000000 /)
!
ref_run_u = line
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ianz<1) THEN                   ! Macro file name is needed
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num/=0) RETURN
!
linit = refine_init .OR. opara(OINIT)=='yes' .or. ubound(refine_beta,1)< refine_par_n
!
IF(.NOT.ALLOCATED(ref_x)) THEN
   ier_num = -1
   ier_typ = ER_APPL
   ier_msg(1) = 'Check if data file was not loaded '
   RETURN
ENDIF
!
refine_mac   = cpara(1)
refine_mac_l = lpara(1)
INQUIRE(FILE=cpara(1)(1:lpara(1)), EXIST=lexist)
IF(.NOT.lexist) THEN              ! Macro does not exist
   ier_num = -2
   ier_typ = ER_APPL
   ier_msg(1) = 'Check if macro name was spelled correctly'
   RETURN
ENDIF
!
if(refine_log) then
   string = 'mkdir -p DISCUS_SUITE_DERIVATIVES'
   call execute_command_line(string, CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg, wait=.true.)
endif
!
ref_do_plot = .FALSE.
refine_plot_mac       = ' '
IF(opara(OPLOT) /= ' ') THEN
  ref_do_plot = .TRUE.
  refine_plot_mac       = opara(OPLOT)
ENDIF
!
ALLOCATE(refine_calc  (ref_dim(1), ref_dim(2), ref_dim(3)))
ALLOCATE(refine_temp  (ref_dim(1), ref_dim(2), ref_dim(3)))
ALLOCATE(refine_derivs(ref_dim(1), ref_dim(2), ref_dim(3), refine_par_n))
IF(linit) THEN
  IF(ALLOCATED(refine_cl)) DEALLOCATE(refine_cl)
  ALLOCATE(refine_cl(refine_par_n, refine_par_n))
  IF(ALLOCATED(refine_alpha)) DEALLOCATE(refine_alpha)
  ALLOCATE(refine_alpha(refine_par_n, refine_par_n))
  IF(ALLOCATED(refine_beta )) DEALLOCATE(refine_beta )
  ALLOCATE(refine_beta (refine_par_n              ))
  refine_cl    = 0.0
  refine_alpha = 0.0
  refine_beta  = 0.0
  lconvergence = .FALSE.
  refine_chisqr  = -1.0000
  refine_rval    = -1.0
ENDIF
!
CALL gl_set_use(.TRUE.)          ! Turn usage of global data on
CALL gl_set_npara(refine_par_n)  ! Define number of parameters
dimen(1) = ref_dim(1)
dimen(2) = ref_dim(2)
dimen(3) = ref_dim(3)
dimen(4) = refine_par_n
CALL gl_alloc(dimen)             ! Turn usage of global data on
call gl_set_data(dimen(1),dimen(2), dimen(3), -2, ref_data( :,:,:))
call gl_set_data(dimen(1),dimen(2), dimen(3), -3, ref_sigma(:,:,:))
call gl_set_x(dimen(1), ref_x)
call gl_set_y(dimen(2), ref_y)
call gl_set_z(dimen(3), ref_z)
!
!
diffev_l_get_random_state = l_get_random_state
l_get_random_state = 1    ! We need the slaves to maintain fixe random seeds
! Initialize MPI interface
!
if(mpi_active) then
   call refine_init_pop
endif
!
l_old_plot_status = l_plot_status   ! save plot status
l_plot_status = .FALSE.             ! Turn plotting off during refinements
!
! Call main refinement routine
!
!write(*,*) 'INIT    ', linit
!write(*,*) ' LAMBDA ', refine_lamda
!write(*,*) 'START MRQ ', linit, REF_MAXPARAM, refine_par_n, REF_MAXPARAM_FIX, refine_fix_n
!write(*,*) 'cycles    ', refine_cycles, ref_kupl
!write(*,*) 'PARAMS    ', lbound(refine_params), ubound(refine_params)
!write(*,*) 'dims      ', ref_dim
!write(*,*) 'DATA      ', lbound(ref_data), ubound(ref_data)
!write(*,*) 'SIGMA     ', lbound(ref_sigma), ubound(ref_sigma)
!write(*,*) 'x         ', lbound(ref_x    ), ubound(ref_x    )
!write(*,*) 'y         ', lbound(ref_y    ), ubound(ref_y    )
!write(*,*) 'conv      ', lconvergence,                                                   &
!                refine_chisqr, refine_conf, refine_lamda, refine_lamda_s,       &
!                refine_lamda_d, refine_lamda_u
!write(*,*) 'RVAL      ', refine_rval,refine_rexp
!write(*,*) 'P         ', lbound(refine_p), ubound(refine_p)
!write(*,*) 'range     ', lbound(refine_range), ubound(refine_range)
!write(*,*) 'shift     ', lbound(refine_shift), ubound(refine_shift)
!write(*,*) 'nderiv    ', lbound(refine_nderiv), ubound(refine_nderiv)
!write(*,*) 'dp        ', lbound(refine_dp), ubound(refine_dp)
!write(*,*) 'cl        ', lbound(refine_cl), ubound(refine_cl)
!write(*,*) 'alpha     ', lbound(refine_alpha), ubound(refine_alpha)
!write(*,*) 'beta      ', lbound(refine_beta ), ubound(refine_beta )
!write(*,*) 'do_plot   ', ref_do_plot, refine_plot_mac(1:len_trim(refine_plot_mac))
!write(*,*) 'Para      ', refine_p(1:refine_par_n)
CALL refine_mrq(linit, REF_MAXPARAM, refine_par_n, refine_cycles, ref_kupl,     &
                refine_params, ref_dim, ref_data, ref_sigma, ref_x, ref_y,      &
                ref_z,                                                          &
                conv_status, conv_dp_sig, conv_dchi2, conv_chi2, conv_conf,     &
                conv_lambda, lconvergence, lconv,                               &
                refine_chisqr, refine_conf, refine_lamda, refine_lamda_s,       &
                refine_lamda_d, refine_lamda_u, refine_rval,                    &
                refine_rexp, refine_p, refine_range, refine_shift, refine_nderiv,&
                refine_dp, refine_cl, refine_alpha, refine_beta, ref_do_plot, refine_plot_mac)
!
l_plot_status = l_old_plot_status  ! Return plot status to old value
!
main:IF(ier_num==0) THEN
   refine_init = .FALSE.
!
! Copy fixed parameters for output
!
   DO i=1, refine_fix_n
      cpara(1) = refine_fixed(i)
      lpara(1) = LEN_TRIM(refine_fixed(i))
      ianz = 1
      CALL ber_params(ianz, cpara, lpara, werte, MAXW)
      IF(ier_num/=0) EXIT main
      refine_f(i) = werte(1)
   ENDDO
   ndata = ref_dim(1)*ref_dim(2)
   CALL show_fit_erg(output_io, REF_MAXPARAM, REF_MAXPARAM_FIX, refine_par_n, refine_fix_n,   &
                     ndata, refine_mac, refine_mac_l, ref_load, ref_kload,         &
                     ref_csigma, ref_ksigma, .FALSE., refine_chisqr, refine_conf,  &
                     refine_lamda, refine_rval, refine_rexp, refine_params,        &
                     refine_p, refine_dp, refine_range, refine_cl, refine_fixed,   &
                        refine_f, lconv,                                           &
                     conv_status, conv_dp_sig, conv_dchi2, conv_chi2, conv_conf, conv_lambda     )
!
ENDIF main
DEALLOCATE(refine_calc)      ! Clean up temporary files
DEALLOCATE(refine_temp)
DEALLOCATE(refine_derivs)
!
lmacro_close = .TRUE.        ! Do close macros in do-loops
!
l_get_random_state = diffev_l_get_random_state   ! Restore diffev random state variable
if(mpi_active) then
   string = 'rm -rf   DISCUS_SUITE_DERIVATIVES'
   if(.not. refine_log) then
      call execute_command_line(string, CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg, wait=.true.)
   endif
endif
!
END SUBROUTINE refine_run
!
!*******************************************************************************
!
SUBROUTINE refine_theory(MAXP, inx, iny, inz, xx, yy, zz, NPARA, p, par_names,  &
                         prange,  p_shift, p_nderiv, &
                         data_dim, & !data_data, data_sigma, data_x, data_y, &
                         kupl_last, &
                         f, df, LDERIV)
!
! Calculates the "theoretical" value through the user macro
! Only if ix and iy == 1 the macro is actually called, otherwise
! the lookup values will be used
!
use refine_params_mod
USE refine_random_mod
USE refine_set_param_mod
!
USE errlist_mod
USE matrix_mod
USE precision_mod
!
use population
!
USE global_data_mod
use mpi_slave_mod
!
use parallel_mod
!OMP !$ use omp_lib
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: inx     ! Point number along x
INTEGER                                              , INTENT(IN)  :: iny     ! Point number along y
INTEGER                                              , INTENT(IN)  :: inz     ! Point number along z
REAL(kind=PREC_DP)                                   , INTENT(IN)  :: xx      ! Point value  along x
REAL(kind=PREC_DP)                                   , INTENT(IN)  :: yy      ! Point value  along y
REAL(kind=PREC_DP)                                   , INTENT(IN)  :: zz      ! Point value  along z
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL(kind=PREC_DP), DIMENSION(MAXP )                 , INTENT(IN)  :: p       ! Parameter values
CHARACTER(LEN=*)  , DIMENSION(MAXP)                  , INTENT(IN)  :: par_names    ! Parameter names
REAL(kind=PREC_DP), DIMENSION(MAXP, 2               ), INTENT(IN)  :: prange      ! Allowed parameter range
REAL(kind=PREC_DP), DIMENSION(MAXP )                 , INTENT(IN)  :: p_shift   ! Parameter shift for deriv
INTEGER           , DIMENSION(MAXP )                 , INTENT(IN)  :: p_nderiv  ! Number of points for derivative
INTEGER           , DIMENSION(3)                     , INTENT(IN)  :: data_dim     ! Data array dimensions
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL(kind=PREC_DP)                                                 , INTENT(OUT) :: f       ! Function value at (ix,iy)
REAL(kind=PREC_DP), DIMENSION(NPARA)                 , INTENT(OUT) :: df      ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!
!
INTEGER              :: k ! Dummy loop variable
!REAL(kind=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: refine_derivs_p   ! Derivativs (ix, iy, iz)
!INTEGER              :: nder          ! Numper of points for derivative
!REAL(KIND=PREC_DP)                 :: delta         ! Shift to calculate derivatives
!REAL(KIND=PREC_DP), DIMENSION(-2:2):: dvec          ! Parameter at P-2delta, P-delta, P, P+delta and P+2delta
!LOGICAL           , DIMENSION(-2:2):: lvec          ! Calculate at P-2delta, P-delta, P, P+delta and P+2delta
!REAL(KIND=PREC_DP), DIMENSION(3)   :: yvec          ! Chisquared at P, P+delta and P-delta
!REAL(KIND=PREC_DP), DIMENSION(3)   :: avec          ! Params for derivative y = a + bx + cx^2
!REAL(KIND=PREC_DP), DIMENSION(3,3) :: xmat          ! Rows are : 1, P, P^2
!REAL(KIND=PREC_DP), DIMENSION(3,3) :: imat          ! Inverse to xmat
!
!REAL(kind=PREC_DP), DIMENSION(:,:,:,:), ALLOCATABLE :: refine_tttt     ! Calculated (ix, iy) for derivs
!
!  OMP variables
integer :: tid
integer :: nthreads
logical :: lserial
!character(len=PREC_STRING) :: string
!integer :: iix
!
tid = 0
nthreads = 1
lserial = .true.
!OMP if(.not.lserial .and. par_omp_use) then
!OMP !$OMP parallel private(tid)
!OMP !$   tid = OMP_GET_THREAD_NUM()
!OMP !$   if(tid == 0) THEN
!OMP !$      if(par_omp_maxthreads == -1) THEN
!OMP !$         nthreads = OMP_GET_NUM_THREADS()
!OMP !$      else
!OMP !$         nthreads = max(1,MIN(par_omp_maxthreads, OMP_GET_NUM_THREADS()))
!OMP !$      endif
!OMP !$   endif
!OMP !$OMP end parallel
!OMP endif
nthreads = 1
!
initial: IF(inx==1 .AND. iny==1 .and. inz==1) THEN   ! Initial point, call user macro
!
   DO k=1, NPARA                      ! Update user defined parameter variables
      CALL refine_set_param(NPARA, par_names(k), k, p(k))
   ENDDO 
!
   CALL refine_save_seeds             ! Save current random number seeds
   CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
                     data_dim, refine_calc)
   CALL gl_set_data(data_dim(1), data_dim(2), data_dim(3), -1,  &
        refine_calc(1:data_dim(1),1:data_dim(2),1:data_dim(3)))
!
   IF(ier_num /= 0) RETURN
!
!  Loop over all parameters to derive derivatives
!
   if_deriv: IF(LDERIV) THEN
      cond_MPI: if(mpi_active) then
         call refine_calc_deriv_mpi(NPARA, MAXP, REF_MAXPARAM_SPC, data_dim,         &
                 p, par_names, refine_params, refine_spc_name, refine_spc_delta, &
                 prange, p_shift, p_nderiv, refine_mac, refine_mac_l, kupl_last, &
                              refine_derivs)
      else  cond_MPI
!
!OMP if(nthreads==1) then                ! Serial loop
         loop_deriv_serial: DO k=1, NPARA
            call refine_calc_deriv(NPARA, MAXP, REF_MAXPARAM_SPC, data_dim, k, &
              p, par_names, refine_params, refine_spc_name, refine_spc_delta, &
              prange, p_shift, p_nderiv, refine_mac, refine_mac_l, kupl_last, &
                           refine_derivs)
         ENDDO loop_deriv_serial
!OMP elseif(nthreads>1) then
!OMP !$OMP PARALLEL
!OMP    if(allocated(refine_derivs_p)) deallocate(refine_derivs_p)
!OMP    allocate(refine_derivs_p(data_dim(1), data_dim(2), data_dim(3)))
!OMP $OMP DO SCHEDULE(STATIC)
!OMP       loop_deriv: DO k=1, NPARA
!OMP          call refine_calc_deriv_p(NPARA, MAXP, REF_MAXPARAM_SPC, data_dim, k, &
!OMP            p, par_names, refine_params, refine_spc_name, refine_spc_delta, &
!OMP            prange, p_shift, p_nderiv, refine_mac, refine_mac_l, kupl_last, &
!OMP                         refine_derivs_p)
!OMP !          refine_temp, refine_derivs)
!OMP !write(* , '(a,6f12.6)') 'SUBROUTINE ', (refine_derivs(1,1,1,k))
!OMP !        ALLOCATE(refine_tttt(1:data_dim(1), 1:data_dim(2), -p_nderiv(k)/2:p_nderiv(k)/2))
!OMP          refine_derivs(:,:,:,k) = refine_derivs(:,:,:,k) + refine_derivs_p(:,:,:)
!OMP       ENDDO loop_deriv
!OMP !$OMP END DO
!OMP !$OMP END PARALLEL
!OMP endif
      endif cond_MPI
!do k=1, npara
!write(string,'(a,i4.4)') 'DISCUS_SUITE_DERIVATIVES/DERIVATIVE', k
!open(55, file=string(1:len_trim(string)), status='unknown')
!do iix=1, data_dim(1)
!write(55,'(2(g20.12e3,2x))') real(iix), refine_derivs(iix, 1  , 1  , k)
!enddo
!close(55)
!enddo
   ENDIF if_deriv
   if(.not.(gl_is_der(1) .and. gl_is_der(2))) then
            CALL refine_restore_seeds
   CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
                     data_dim, refine_calc)
   endif
!write(*,*) ' GOT FINAL   ', lderiv, minval(refine_calc),maxval(refine_calc), p(1:2)
!write(*,*) ' DERIVS Fin  ', 1, minval(refine_derivs(:,:,:,1)), maxval(refine_derivs(:,:,:,1))
!write(*,*) ' DERIVS Fin  ', 2, minval(refine_derivs(:,:,:,2)), maxval(refine_derivs(:,:,:,2))
!write(*,*) ' C11 CMM     ', 0, refine_calc(1,1,1), refine_calc(data_dim(1),data_dim(2),1)
!write(*,*) ' D11 DMM     ', 1, refine_derivs(1,1,1,1), refine_derivs(data_dim(1),data_dim(2),1,1)
!write(*,*) ' D11 DMM     ', 2, refine_derivs(1,1,1,2), refine_derivs(data_dim(1),data_dim(2),1,2)
!write(* , '(6f12.6)') (refine_derivs(1,1,1,k), k=1, NPARA)
!open(33, file='DERIV/deriv.1', status='unknown')
!do l=1, data_dim(1)
!write(33, '(6f12.6)') refine_derivs(l,1,1,k), k=1, NPARA
!enddo
!close(33)
!open(33, file='DERIV/deriv.1', status='unknown')
!write(33, '(2i4)') 21, 21
!write(33, '(4f8.2)') -1.0 , 1.0, -1.0, 1.0
!do i=1, 21
!  write(33, '(5f12.6)') refine_derivs(:,i,1,1)
!enddo
!close(33)
!open(33, file='DERIV/deriv.2', status='unknown')
!write(33, '(2i4)') 21, 21
!write(33, '(4f8.2)') -1.0 , 1.0, -1.0, 1.0
!do i=1, 21
!  write(33, '(5f12.6)') refine_derivs(:,i,1,2)
!enddo
!close( 33)
!write(*,*) ' WROTE DERIVS '
!read(*,*) i
!
ENDIF initial
!
f = refine_calc(inx, iny, inz)           ! Function value to be returned
!iz = 1
IF(LDERIV) THEN                   ! Derivatives are needed
   DO k=1, NPARA
      df(k) = refine_derivs(inx, iny, inz, k)
   ENDDO
ENDIF
!
END SUBROUTINE refine_theory
!
!*******************************************************************************
!
subroutine refine_calc_deriv(NPARA, MAXP, REF_MAXPARAM_SPC, data_dim, k, &
           p, par_names, refine_params, refine_spc_name, refine_spc_delta, &
           prange, p_shift, p_nderiv, refine_mac, refine_mac_l, kupl_last, &
                        refine_derivs)
!          refine_temp, refine_derivs)
!-
! Calculate the k derivatives
!+
!
use refine_random_mod
use refine_set_param_mod
!
use errlist_mod
use global_data_mod
use matrix_mod
use precision_mod
!
implicit none
!
integer , intent(in) :: NPARA    ! Total number of parameter
integer , intent(in) :: MAXP     ! Max  number of parameters
integer , intent(in) :: REF_MAXPARAM_SPC     ! Max  number of special parameters
integer, dimension(3)    , intent(in) :: data_dim         ! Data dimensions
integer , intent(in) :: k        ! Current derivative number
real(kind=PREC_DP)        , dimension(MAXP)            , intent(in) :: p                ! Parameter values
character(len=*)          , dimension(MAXP)            , intent(in)  :: par_names    ! Parameter names
character(len=PREC_STRING), dimension(MAXP)            , intent(in) :: refine_params    ! Parameter names
character(len=PREC_STRING), dimension(REF_MAXPARAM_SPC), intent(in) :: refine_spc_name ! Parameter names, special
real(kind=PREC_DP)        , dimension(REF_MAXPARAM_SPC), intent(in) :: refine_spc_delta ! Special parameter shift
real(kind=PREC_DP)        , dimension(MAXP, 2)         , intent(in) :: prange           ! Parameter range
real(kind=PREC_DP)        , dimension(MAXP)            , intent(in) :: p_shift ! Parameter shift
integer                   , dimension(MAXP)            , intent(in) :: p_nderiv! Parameter uses this many points for derivatives
character(len=*)                                       , intent(inout) :: refine_mac    ! Refine macro name
integer                                                , intent(inout) :: refine_mac_l  ! Macro name length
integer                                                , intent(in) :: kupl_last     
!logical           , dimension(NPARA), intent(in) :: gl_is_der        ! Derivative has been calculated analytically
!real(kind=PREC_DP)        , dimension(data_dim(1), data_dim(2), data_dim(3))  , intent(inout) :: refine_temp           ! caclulated value
real(kind=PREC_DP)        , dimension(data_dim(1), data_dim(2), data_dim(3),NPARA) , intent(out) :: refine_derivs         ! caclulated derivatives
!integer                                , intent(out) :: l_ier_num     
!integer                                , intent(out) :: l_ier_typ     
!character         , dimension(7)       , intent(out) :: l_ier_msg     
!
integer :: i, l    ! Dummy loop indeces
integer :: j2, j3  ! Dummy loop indeces
integer :: iix, iiy, iiz ! Dummy loop indices
integer :: nder    ! Number of derivatives calculated
logical           , dimension(-3:3) :: lvec
real(kind=PREC_DP), dimension(-3:3) :: dvec
real(kind=PREC_DP)                  :: delta
real(kind=PREC_DP), dimension(:,:,:  ), allocatable :: refine_temp           ! Calculated values at derivative l
real(kind=PREC_DP), dimension(:,:,:,:), allocatable :: refine_tttt           ! Calculated values at derivative l
real(kind=PREC_DP), dimension(3,3)  :: xmat   ! matrix with derivatives
real(kind=PREC_DP), dimension(3,3)  :: imat   ! Inverse of xmax
real(kind=PREC_DP), dimension(3)    :: avec   ! Inverse of xmax
real(kind=PREC_DP), dimension(3)    :: yvec   ! Inverse of xmax
!
ALLOCATE(refine_temp(1:data_dim(1), 1:data_dim(2), 1:data_dim(3)))
!
is_deriv: IF(gl_is_der(k)) THEN
   CALL gl_get_data(k,  data_dim(1),   data_dim(2),   data_dim(3), &
        refine_derivs(1:data_dim(1), 1:data_dim(2), 1:data_dim(3), k))
ELSE is_deriv
   ALLOCATE(refine_tttt(1:data_dim(1), 1:data_dim(2), 1:data_dim(3), -2:2))
         nder = 1                     ! First point is at P(k)
         dvec(:) = 0.0
         lvec(:) = .FALSE.
         dvec(0) = p(k)               ! Store parameter value
         if_delta: IF(p(k)/=0.0) THEN
            delta = ABS(p(k)*p_shift(k))  ! A multiplicative variation of the parameter seems best
         ELSE if_delta
            do i=1, REF_MAXPARAM_SPC      ! Check for special names
               if(refine_params  (k)(1:len_trim(refine_spc_name(i)))==          &
                  refine_spc_name(i)(1:len_trim(refine_spc_name(i)))     ) then
                  delta = refine_spc_delta(i)
                  exit if_delta
               endif
            enddo
            delta = 1.0D-4                ! Default delta for unknown variables
         ENDIF if_delta
!
!                                     ! Test at P + DELTA
         IF(prange(k,1)<=prange(k,2)) THEN     ! User provided parameter range
            IF(p(k)==prange(k,1)) THEN         ! At lower limit, use +delta +2delta
               delta = MIN(ABS(delta), 0.5D0*(ABS(prange(k,2) - p(k)))) ! Make sure +2Delta fits
               dvec(1) = p(k) + delta
               dvec(2) = p(k) + 2.0D0*delta
               lvec(1) = .TRUE.
               lvec(2) = .TRUE.
               nder = 3
            ELSEIF(p(k)==prange(k,2)) THEN         ! At upper limit, use -delta -2delta
               delta = MIN(ABS(delta), 0.5D0*(ABS(p(k) - prange(k,1)))) ! Make sure -2Delta fits
               dvec(-2) = p(k) - 2.0D0*delta
               dvec(-1) = p(k) - 1.0D0*delta
               lvec(-2) = .TRUE.
               lvec(-1) = .TRUE.
               nder = 3
            ELSE                               ! within range, use +-delta or +2delta
               IF(p_nderiv(k)==3) THEN         ! Three point derivative
                  dvec(-1) = MIN(prange(k,2),MAX(prange(k,1),p(k)+(delta))) ! +delta
                  dvec( 1) = MIN(prange(k,2),MAX(prange(k,1),p(k)-(delta))) ! -delta
                  lvec(-1) = .TRUE.
                  lvec( 1) = .TRUE.
                  nder = 3
               ELSEIF(p_nderiv(k)==5) THEN     ! Five point derivative
                  delta = MIN(delta, 0.5D0*(ABS(prange(k,2) - p(k))),  &
                                     0.5D0*(ABS(p(k) - prange(k,1))) ) ! Make sure +-2delta fits
                  dvec(-2) = p(k) - 2.0D0*delta
                  dvec(-1) = p(k) - 1.0D0*delta
                  dvec( 1) = p(k) + 1.0D0*delta
                  dvec( 2) = p(k) + 2.0D0*delta
                  lvec(-2) = .TRUE.
                  lvec(-1) = .TRUE.
                  lvec( 1) = .TRUE.
                  lvec( 2) = .TRUE.
                  nder = 5
               ENDIF
            ENDIF
!           p_d     = MIN(prange(k,2),MAX(prange(k,1),p(k)+REAL(delta)))
         ELSE
            IF(p_nderiv(k)==3) THEN         ! Three point derivative
               dvec(-1) = p(k) - 1.0D0*delta
               dvec( 1) = p(k) + 1.0D0*delta
               lvec(-1) = .TRUE.
               lvec( 1) = .TRUE.
               nder = 3
            ELSEIF(p_nderiv(k)==5) THEN     ! Five point derivative
               dvec(-2) = p(k) - 2.0D0*delta
               dvec(-1) = p(k) - 1.0D0*delta
               dvec( 0) = p(k)
               dvec( 1) = p(k) + 1.0D0*delta
               dvec( 2) = p(k) + 2.0D0*delta
               lvec(-2) = .TRUE.
               lvec(-1) = .TRUE.
               lvec( 1) = .TRUE.
               lvec( 2) = .TRUE.
               nder = 5
            ENDIF
!           p_d     = p(k) + delta
         ENDIF
         DO l = -2, 2
            IF(lvec(l)) THEN
            CALL refine_set_param(NPARA, par_names(k), k, (dvec(l)) )  ! Set modified value
            CALL refine_restore_seeds
            CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
                              data_dim, refine_temp)
!write(*,*) ' GOT DERIV   ', k
            IF(ier_num /= 0) THEN
               DEALLOCATE(refine_tttt)
               exit is_deriv
            ENDIF
            do iiz = 1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
                  refine_tttt(iix,iiy,iiz,l) =  refine_temp(iix,iiy,iiz)
               ENDDO
            ENDDO
            enddo
!do iiz=-2, 2
!write(line,'(a,i1.1)') 'CALC/expli.',iiz+3
!open(33,file = line, status='unknown')
!write(33,'(2i5)') data_dim(1:2)
!write(33,'(4f8.2)') -1.0, 1.0, -1.0, 1.0
!do iiy=1, data_dim(2)
!  write(33, '(5f12.6)') refine_tttt(:,iiy,1,iiz)
!enddo
!close(33)
!enddo
            ENDIF
         ENDDO
call refine_rvalue_tttt(data_dim, refine_tttt, dvec)
         CALL refine_set_param(NPARA, par_names(k), k, p(k))  ! Return to original value
!
         IF(nder==5) THEN             ! Got all five  points for derivative
            do iiz=1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
                  refine_derivs(iix, iiy, iiz, k) = (-1.0*refine_tttt(iix,iiy,iiz, 2)   &
                                                     +8.0*refine_tttt(iix,iiy,iiz, 1)   &
                                                     -8.0*refine_tttt(iix,iiy,iiz,-1)   &
                                                     +1.0*refine_tttt(iix,iiy,iiz,-2))/ &
                                                    (12.*delta)
               ENDDO
            ENDDO
            enddo
         ELSEIF(nder==3) THEN             ! Got all three points for derivative
            xmat(:,1) =  1.0
            xmat(1,2) =  dvec(0) !p(k)
            IF(lvec(2)) THEN              ! +delta, + 2delta
              j2 =  1
              j3 =  2
            ELSEIF(lvec(-2)) THEN         ! -delta, - 2delta
              j2 = -1
              j3 = -2
            ELSE
              j2 = -1
              j3 =  1
            ENDIF
              xmat(2,2) =  dvec(j2)       ! p(k) - delta
              xmat(3,2) =  dvec(j3)       ! p(k) + delta
              xmat(2,3) = (dvec(j2))**2   ! (p(k) - delta ) **2
              xmat(3,3) = (dvec(j3))**2   ! (p(k) + delta ) **2
!           xmat(2,2) =  dvec(2) !p(k) + delta
!           xmat(3,2) =  dvec(3) !p(k) - delta
!           xmat(1,3) = (dvec(0))**2 !(p(k)        ) **2
!           xmat(2,3) = (dvec(2))**2 !(p(k) + delta) **2
!           xmat(3,3) = (dvec(3))**2 !(p(k) - delta) **2
            CALL matinv3(xmat, imat)
            IF(ier_num/=0) then
               ier_msg(1) = 'Error determining derivative '
               write(ier_msg(2), '(a,i3)') 'At derivative ',k
               exit is_deriv
            endif
            do iiz=1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
!
!              Derivative is calculated as a fit of a parabola at P, P+delta, P-delta
                  yvec(1) = refine_calc  (iix, iiy, iiz)
                  yvec(2) = refine_tttt  (iix, iiy, iiz, j2)
                  yvec(3) = refine_tttt  (iix, iiy, iiz, j3)
                  avec = MATMUL(imat, yvec)
!
                  refine_derivs(iix, iiy, iiz, k) = avec(2) + 2.*avec(3)*p(k)
               ENDDO
            ENDDO
            enddo
!        ELSEIF(nder==2) THEN
!           IF(dvec(2)==0) THEN          ! P + Delta failed
!              DO iiy=1, data_dim(2)
!                 DO iix=1, data_dim(1)
!                    refine_derivs(iix, iiy, k) = (refine_temp  (iix, iiy)-refine_calc  (iix, iiy))/ &
!                                                 (dvec(3)-dvec(1))
!                 ENDDO
!              ENDDO
!           ELSEIF(dvec(3)==0) THEN      ! P - Delta failed
!              DO iiy=1, data_dim(2)
!                 DO iix=1, data_dim(1)
!                    refine_derivs(iix, iiy, k) = (refine_derivs(iix, iiy,k)-refine_calc  (iix, iiy))/ &
!                                                 (dvec(2)-dvec(1))
!                 ENDDO
!              ENDDO
!           ENDIF
!
         ELSE
            ier_num = -9
            ier_typ = 6
            ier_msg(1) = par_names(k)
            DEALLOCATE(refine_tttt)
            exit is_deriv
         ENDIF
         DEALLOCATE(refine_tttt)
!
         ENDIF is_deriv
deallocate(refine_temp)
!
end subroutine refine_calc_deriv
!
!*******************************************************************************
!
subroutine refine_calc_deriv_p(NPARA, MAXP, REF_MAXPARAM_SPC, data_dim, k, &
           p, par_names, refine_params, refine_spc_name, refine_spc_delta, &
           prange, p_shift, p_nderiv, refine_mac, refine_mac_l, kupl_last, &
                        refine_derivs_p)
!          refine_temp, refine_derivs)
!-
! Calculate the k derivative parallel version
!+
!
use refine_random_mod
use refine_set_param_mod
!
use errlist_mod
use global_data_mod
use matrix_mod
use precision_mod
!
implicit none
!
integer , intent(in) :: NPARA    ! Total number of parameter
integer , intent(in) :: MAXP     ! Max  number of parameters
integer , intent(in) :: REF_MAXPARAM_SPC     ! Max  number of special parameters
integer, dimension(3)    , intent(in) :: data_dim         ! Data dimensions
integer , intent(in) :: k        ! Current derivative number
real(kind=PREC_DP)        , dimension(MAXP)            , intent(in) :: p                ! Parameter values
character(len=*)          , dimension(MAXP)            , intent(in)  :: par_names    ! Parameter names
character(len=PREC_STRING), dimension(MAXP)            , intent(in) :: refine_params    ! Parameter names
character(len=PREC_STRING), dimension(REF_MAXPARAM_SPC), intent(in) :: refine_spc_name ! Parameter names, special
real(kind=PREC_DP)        , dimension(REF_MAXPARAM_SPC), intent(in) :: refine_spc_delta ! Special parameter shift
real(kind=PREC_DP)        , dimension(MAXP, 2)         , intent(in) :: prange           ! Parameter range
real(kind=PREC_DP)        , dimension(MAXP)            , intent(in) :: p_shift ! Parameter shift
integer                   , dimension(MAXP)            , intent(in) :: p_nderiv! Parameter uses this many points for derivatives
character(len=*)                                       , intent(inout) :: refine_mac    ! Refine macro name
integer                                                , intent(inout) :: refine_mac_l  ! Macro name length
integer                                                , intent(in) :: kupl_last     
!logical           , dimension(NPARA), intent(in) :: gl_is_der        ! Derivative has been calculated analytically
!real(kind=PREC_DP)        , dimension(data_dim(1), data_dim(2), data_dim(3))  , intent(inout) :: refine_temp           ! caclulated value
real(kind=PREC_DP)        , dimension(data_dim(1), data_dim(2), data_dim(3)) , intent(out) :: refine_derivs_p         ! caclulated derivatives
!integer                                , intent(out) :: l_ier_num     
!integer                                , intent(out) :: l_ier_typ     
!character         , dimension(7)       , intent(out) :: l_ier_msg     
!
integer :: i, l    ! Dummy loop indeces
integer :: j2, j3  ! Dummy loop indeces
integer :: iix, iiy, iiz ! Dummy loop indices
integer :: nder    ! Number of derivatives calculated
logical           , dimension(-3:3) :: lvec
real(kind=PREC_DP), dimension(-3:3) :: dvec
real(kind=PREC_DP)                  :: delta
real(kind=PREC_DP), dimension(:,:,:  ), allocatable :: refine_temp           ! Calculated values at derivative l
real(kind=PREC_DP), dimension(:,:,:,:), allocatable :: refine_tttt           ! Calculated values at derivative l
real(kind=PREC_DP), dimension(3,3)  :: xmat   ! matrix with derivatives
real(kind=PREC_DP), dimension(3,3)  :: imat   ! Inverse of xmax
real(kind=PREC_DP), dimension(3)    :: avec   ! Inverse of xmax
real(kind=PREC_DP), dimension(3)    :: yvec   ! Inverse of xmax
!
ALLOCATE(refine_temp(1:data_dim(1), 1:data_dim(2), 1:data_dim(3)))
!
is_deriv: IF(gl_is_der(k)) THEN
   CALL gl_get_data(k,  data_dim(1),   data_dim(2),   data_dim(3), &
        refine_derivs_p(1:data_dim(1), 1:data_dim(2), 1:data_dim(3)))
!write(*,*) ' DERIVS Calc ', k, minval(refine_derivs(:,:,:,k)), maxval(refine_derivs(:,:,:,k))
ELSE is_deriv
   ALLOCATE(refine_tttt(1:data_dim(1), 1:data_dim(2), 1:data_dim(3), -2:2))
         nder = 1                     ! First point is at P(k)
         dvec(:) = 0.0
         lvec(:) = .FALSE.
         dvec(0) = p(k)               ! Store parameter value
         if_delta: IF(p(k)/=0.0) THEN
            delta = ABS(p(k)*p_shift(k))  ! A multiplicative variation of the parameter seems best
         ELSE if_delta
            do i=1, REF_MAXPARAM_SPC      ! Check for special names
               if(refine_params  (k)(1:len_trim(refine_spc_name(i)))==          &
                  refine_spc_name(i)(1:len_trim(refine_spc_name(i)))     ) then
                  delta = refine_spc_delta(i)
                  exit if_delta
               endif
            enddo
            delta = 1.0D-4                ! Default delta for unknown variables
         ENDIF if_delta
!
!                                     ! Test at P + DELTA
         IF(prange(k,1)<=prange(k,2)) THEN     ! User provided parameter range
            IF(p(k)==prange(k,1)) THEN         ! At lower limit, use +delta +2delta
               delta = MIN(ABS(delta), 0.5D0*(ABS(prange(k,2) - p(k)))) ! Make sure +2Delta fits
               dvec(1) = p(k) + delta
               dvec(2) = p(k) + 2.0D0*delta
               lvec(1) = .TRUE.
               lvec(2) = .TRUE.
               nder = 3
            ELSEIF(p(k)==prange(k,2)) THEN         ! At upper limit, use -delta -2delta
               delta = MIN(ABS(delta), 0.5D0*(ABS(p(k) - prange(k,1)))) ! Make sure -2Delta fits
               dvec(-2) = p(k) - 2.0D0*delta
               dvec(-1) = p(k) - 1.0D0*delta
               lvec(-2) = .TRUE.
               lvec(-1) = .TRUE.
               nder = 3
            ELSE                               ! within range, use +-delta or +2delta
               IF(p_nderiv(k)==3) THEN         ! Three point derivative
                  dvec(-1) = MIN(prange(k,2),MAX(prange(k,1),p(k)+(delta))) ! +delta
                  dvec( 1) = MIN(prange(k,2),MAX(prange(k,1),p(k)-(delta))) ! -delta
                  lvec(-1) = .TRUE.
                  lvec( 1) = .TRUE.
                  nder = 3
               ELSEIF(p_nderiv(k)==5) THEN     ! Five point derivative
                  delta = MIN(delta, 0.5D0*(ABS(prange(k,2) - p(k))),  &
                                     0.5D0*(ABS(p(k) - prange(k,1))) ) ! Make sure +-2delta fits
                  dvec(-2) = p(k) - 2.0D0*delta
                  dvec(-1) = p(k) - 1.0D0*delta
                  dvec( 1) = p(k) + 1.0D0*delta
                  dvec( 2) = p(k) + 2.0D0*delta
                  lvec(-2) = .TRUE.
                  lvec(-1) = .TRUE.
                  lvec( 1) = .TRUE.
                  lvec( 2) = .TRUE.
                  nder = 5
               ENDIF
            ENDIF
!           p_d     = MIN(prange(k,2),MAX(prange(k,1),p(k)+REAL(delta)))
         ELSE
            IF(p_nderiv(k)==3) THEN         ! Three point derivative
               dvec(-1) = p(k) - 1.0D0*delta
               dvec( 1) = p(k) + 1.0D0*delta
               lvec(-1) = .TRUE.
               lvec( 1) = .TRUE.
               nder = 3
            ELSEIF(p_nderiv(k)==5) THEN     ! Five point derivative
               dvec(-2) = p(k) - 2.0D0*delta
               dvec(-1) = p(k) - 1.0D0*delta
               dvec( 0) = p(k)
               dvec( 1) = p(k) + 1.0D0*delta
               dvec( 2) = p(k) + 2.0D0*delta
               lvec(-2) = .TRUE.
               lvec(-1) = .TRUE.
               lvec( 1) = .TRUE.
               lvec( 2) = .TRUE.
               nder = 5
            ENDIF
!           p_d     = p(k) + delta
         ENDIF
         DO l = -2, 2
            IF(lvec(l)) THEN
            CALL refine_set_param(NPARA, par_names(k), k, (dvec(l)) )  ! Set modified value
            CALL refine_restore_seeds
            CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
                              data_dim, refine_temp)
!write(*,*) ' GOT DERIV   ', k
            IF(ier_num /= 0) THEN
               DEALLOCATE(refine_tttt)
               exit is_deriv
            ENDIF
            do iiz = 1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
                  refine_tttt(iix,iiy,iiz,l) =  refine_temp(iix,iiy,iiz)
               ENDDO
            ENDDO
            enddo
!do iiz=-2, 2
!write(line,'(a,i1.1)') 'CALC/expli.',iiz+3
!open(33,file = line, status='unknown')
!write(33,'(2i5)') data_dim(1:2)
!write(33,'(4f8.2)') -1.0, 1.0, -1.0, 1.0
!do iiy=1, data_dim(2)
!  write(33, '(5f12.6)') refine_tttt(:,iiy,1,iiz)
!enddo
!close(33)
!enddo
            ENDIF
         ENDDO
call refine_rvalue_tttt(data_dim, refine_tttt, dvec)
         CALL refine_set_param(NPARA, par_names(k), k, p(k))  ! Return to original value
!
         IF(nder==5) THEN             ! Got all five  points for derivative
            do iiz=1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
                  refine_derivs_p(iix, iiy, iiz   ) = (-1.0*refine_tttt(iix,iiy,iiz, 2)   &
                                                     +8.0*refine_tttt(iix,iiy,iiz, 1)   &
                                                     -8.0*refine_tttt(iix,iiy,iiz,-1)   &
                                                     +1.0*refine_tttt(iix,iiy,iiz,-2))/ &
                                                    (12.*delta)
               ENDDO
            ENDDO
            enddo
         ELSEIF(nder==3) THEN             ! Got all three points for derivative
            xmat(:,1) =  1.0
            xmat(1,2) =  dvec(0) !p(k)
            IF(lvec(2)) THEN              ! +delta, + 2delta
              j2 =  1
              j3 =  2
            ELSEIF(lvec(-2)) THEN         ! -delta, - 2delta
              j2 = -1
              j3 = -2
            ELSE
              j2 = -1
              j3 =  1
            ENDIF
              xmat(2,2) =  dvec(j2)       ! p(k) - delta
              xmat(3,2) =  dvec(j3)       ! p(k) + delta
              xmat(2,3) = (dvec(j2))**2   ! (p(k) - delta ) **2
              xmat(3,3) = (dvec(j3))**2   ! (p(k) + delta ) **2
!           xmat(2,2) =  dvec(2) !p(k) + delta
!           xmat(3,2) =  dvec(3) !p(k) - delta
!           xmat(1,3) = (dvec(0))**2 !(p(k)        ) **2
!           xmat(2,3) = (dvec(2))**2 !(p(k) + delta) **2
!           xmat(3,3) = (dvec(3))**2 !(p(k) - delta) **2
            CALL matinv3(xmat, imat)
            IF(ier_num/=0) then
               ier_msg(1) = 'Error determining derivative '
               write(ier_msg(2), '(a,i3)') 'At derivative ',k
               exit is_deriv
            endif
            do iiz=1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
!
!              Derivative is calculated as a fit of a parabola at P, P+delta, P-delta
                  yvec(1) = refine_calc  (iix, iiy, iiz)
                  yvec(2) = refine_tttt  (iix, iiy, iiz, j2)
                  yvec(3) = refine_tttt  (iix, iiy, iiz, j3)
                  avec = MATMUL(imat, yvec)
!
                  refine_derivs_p(iix, iiy, iiz   ) = avec(2) + 2.*avec(3)*p(k)
               ENDDO
            ENDDO
            enddo
!        ELSEIF(nder==2) THEN
!           IF(dvec(2)==0) THEN          ! P + Delta failed
!              DO iiy=1, data_dim(2)
!                 DO iix=1, data_dim(1)
!                    refine_derivs(iix, iiy, k) = (refine_temp  (iix, iiy)-refine_calc  (iix, iiy))/ &
!                                                 (dvec(3)-dvec(1))
!                 ENDDO
!              ENDDO
!           ELSEIF(dvec(3)==0) THEN      ! P - Delta failed
!              DO iiy=1, data_dim(2)
!                 DO iix=1, data_dim(1)
!                    refine_derivs(iix, iiy, k) = (refine_derivs(iix, iiy,k)-refine_calc  (iix, iiy))/ &
!                                                 (dvec(2)-dvec(1))
!                 ENDDO
!              ENDDO
!           ENDIF
!
         ELSE
            ier_num = -9
            ier_typ = 6
            ier_msg(1) = par_names(k)
            DEALLOCATE(refine_tttt)
            exit is_deriv
         ENDIF
         DEALLOCATE(refine_tttt)
!
         ENDIF is_deriv
deallocate(refine_temp)
!
end subroutine refine_calc_deriv_p
!
!*******************************************************************************
!
subroutine refine_calc_deriv_mpi(NPARA, MAXP, REF_MAXPARAM_SPC, data_dim,  &
           p, par_names, refine_params, refine_spc_name, refine_spc_delta, &
           prange, p_shift, p_nderiv, refine_mac, refine_mac_l, kupl_last, &
                        refine_derivs)
!          refine_temp, refine_derivs)
!-
! Calculate the k derivatives
!+
!
use refine_log_mod
use refine_random_mod
use refine_set_param_mod
!
use diffev_mpi_mod
use run_mpi_mod    ! contains "run_mpi_senddata"
use population
!
use kuplot_load_mod
!
use errlist_mod
use global_data_mod
use matrix_mod
use precision_mod
use random_state_mod
!
implicit none
!
integer , intent(in) :: NPARA    ! Total number of parameter
integer , intent(in) :: MAXP     ! Max  number of parameters
integer , intent(in) :: REF_MAXPARAM_SPC     ! Max  number of special parameters
integer, dimension(3)    , intent(in) :: data_dim         ! Data dimensions
real(kind=PREC_DP)        , dimension(MAXP)            , intent(in) :: p                ! Parameter values
character(len=*)          , dimension(MAXP)            , intent(in)  :: par_names    ! Parameter names
character(len=PREC_STRING), dimension(MAXP)            , intent(in) :: refine_params    ! Parameter names
character(len=PREC_STRING), dimension(REF_MAXPARAM_SPC), intent(in) :: refine_spc_name ! Parameter names, special
real(kind=PREC_DP)        , dimension(REF_MAXPARAM_SPC), intent(in) :: refine_spc_delta ! Special parameter shift
real(kind=PREC_DP)        , dimension(MAXP, 2)         , intent(in) :: prange           ! Parameter range
real(kind=PREC_DP)        , dimension(MAXP)            , intent(in) :: p_shift ! Parameter shift
integer                   , dimension(MAXP)            , intent(in) :: p_nderiv! Parameter uses this many points for derivatives
character(len=*)                                       , intent(inout) :: refine_mac    ! Refine macro name
integer                                                , intent(inout) :: refine_mac_l  ! Macro name length
integer                                                , intent(in) :: kupl_last     
!logical           , dimension(NPARA), intent(in) :: gl_is_der        ! Derivative has been calculated analytically
!real(kind=PREC_DP)        , dimension(data_dim(1), data_dim(2), data_dim(3))  , intent(inout) :: refine_temp           ! caclulated value
real(kind=PREC_DP)        , dimension(data_dim(1), data_dim(2), data_dim(3),NPARA) , intent(out) :: refine_derivs         ! caclulated derivatives
!integer                                , intent(out) :: l_ier_num     
!integer                                , intent(out) :: l_ier_typ     
!character         , dimension(7)       , intent(out) :: l_ier_msg     
integer :: k        ! Current derivative number
!
integer :: i, j, l    ! Dummy loop indeces
integer :: j2, j3  ! Dummy loop indeces
integer :: iix, iiy, iiz ! Dummy loop indices
integer           , dimension(     NPARA):: nder    ! Number of derivatives calculated
logical           , dimension(-3:3,NPARA) :: lvec
real(kind=PREC_DP), dimension(-3:3,NPARA) :: dvec
real(kind=PREC_DP), dimension(     NPARA) :: delta
real(kind=PREC_DP), dimension(:,:,:  ), allocatable :: refine_temp           ! Calculated values at derivative l
real(kind=PREC_DP), dimension(:,:,:,:), allocatable :: refine_tttt           ! Calculated values at derivative l
real(kind=PREC_DP), dimension(3,3)  :: xmat   ! matrix with derivatives
real(kind=PREC_DP), dimension(3,3)  :: imat   ! Inverse of xmax
real(kind=PREC_DP), dimension(3)    :: avec   ! Inverse of xmax
real(kind=PREC_DP), dimension(3)    :: yvec   ! Inverse of xmax
!
character(len=PREC_STRING)           :: string   ! Dummy character string
integer, dimension(:,:), allocatable :: kid_is   ! Lookup kid is parameter, derivative number
integer, dimension(:,:), allocatable :: para_has ! Lookup para has these kids
!
allocate(kid_is(2, NPARA*4))
kid_is = 0
allocate(para_has(0:5, NPARA*4))
para_has = 0
!
ALLOCATE(refine_temp(1:data_dim(1), 1:data_dim(2), 1:data_dim(3)))
!
! Initial loop, populate all values in pop_t  Trial for DIFFEV
do k=1, NPARA
  pop_t(k, :) = p(k)
enddo
run_mpi_senddata%children = 0
!
loop_npara_init: do k=1, NPARA
   cond_is_deriv: IF(gl_is_der(k)) THEN
      CALL gl_get_data(k,  data_dim(1),   data_dim(2),   data_dim(3), &
           refine_derivs(1:data_dim(1), 1:data_dim(2), 1:data_dim(3), k))
      ELSE cond_is_deriv
         nder(k) = 1                     ! First point is at P(k)
         dvec(:,k) = 0.0
         lvec(:,k) = .FALSE.
         dvec(0,k) = p(k)               ! Store parameter value
         cond_if_delta: IF(p(k)/=0.0) THEN
            delta(k) = ABS(p(k)*p_shift(k))  ! A multiplicative variation of the parameter seems best
         ELSE cond_if_delta
            do i=1, REF_MAXPARAM_SPC      ! Check for special names
               if(refine_params  (k)(1:len_trim(refine_spc_name(i)))==          &
                  refine_spc_name(i)(1:len_trim(refine_spc_name(i)))     ) then
                  delta(k) = refine_spc_delta(i)
                  exit cond_if_delta
               endif
            enddo
            delta(k) = 1.0D-4                ! Default delta for unknown variables
         ENDIF cond_if_delta
!
!                                     ! Test at P + DELTA
      IF(prange(k,1)<=prange(k,2)) THEN     ! User provided parameter range
         IF(p(k)==prange(k,1)) THEN         ! At lower limit, use +delta +2delta
            delta(k) = MIN(ABS(delta(k)), 0.5D0*(ABS(prange(k,2) - p(k)))) ! Make sure +2Delta fits
            dvec(1,k) = p(k) + delta(k)
            dvec(2,k) = p(k) + 2.0D0*delta(k)
            lvec(1,k) = .TRUE.
            lvec(2,k) = .TRUE.
            nder(k) = 3
         ELSEIF(p(k)==prange(k,2)) THEN         ! At upper limit, use -delta -2delta
            delta(k) = MIN(ABS(delta(k)), 0.5D0*(ABS(p(k) - prange(k,1)))) ! Make sure -2Delta fits
            dvec(-2,k) = p(k) - 2.0D0*delta(k)
            dvec(-1,k) = p(k) - 1.0D0*delta(k)
            lvec(-2,k) = .TRUE.
            lvec(-1,k) = .TRUE.
            nder(k) = 3
         ELSE                               ! within range, use +-delta or +2delta
            IF(p_nderiv(k)==3) THEN         ! Three point derivative
               dvec(-1,k) = MIN(prange(k,2),MAX(prange(k,1),p(k)+(delta(k)))) ! +delta
               dvec( 1,k) = MIN(prange(k,2),MAX(prange(k,1),p(k)-(delta(k)))) ! -delta
               lvec(-1,k) = .TRUE.
               lvec( 1,k) = .TRUE.
               nder(k) = 3
            ELSEIF(p_nderiv(k)==5) THEN     ! Five point derivative
               delta(k) = MIN(delta(k), 0.5D0*(ABS(prange(k,2) - p(k))),  &
                                  0.5D0*(ABS(p(k) - prange(k,1))) ) ! Make sure +-2delta fits
               dvec(-2,k) = p(k) - 2.0D0*delta(k)
               dvec(-1,k) = p(k) - 1.0D0*delta(k)
               dvec( 1,k) = p(k) + 1.0D0*delta(k)
               dvec( 2,k) = p(k) + 2.0D0*delta(k)
               lvec(-2,k) = .TRUE.
               lvec(-1,k) = .TRUE.
               lvec( 1,k) = .TRUE.
               lvec( 2,k) = .TRUE.
               nder(k) = 5
            ENDIF
         ENDIF
!           p_d     = MIN(prange(k,2),MAX(prange(k,1),p(k)+REAL(delta)))
      ELSE
         IF(p_nderiv(k)==3) THEN         ! Three point derivative
            dvec(-1,k) = p(k) - 1.0D0*delta(k)
            dvec( 1,k) = p(k) + 1.0D0*delta(k)
            lvec(-1,k) = .TRUE.
            lvec( 1,k) = .TRUE.
            nder(k) = 3
         ELSEIF(p_nderiv(k)==5) THEN     ! Five point derivative
            dvec(-2,k) = p(k) - 2.0D0*delta(k)
            dvec(-1,k) = p(k) - 1.0D0*delta(k)
            dvec( 0,k) = p(k)
            dvec( 1,k) = p(k) + 1.0D0*delta(k)
            dvec( 2,k) = p(k) + 2.0D0*delta(k)
            lvec(-2,k) = .TRUE.
            lvec(-1,k) = .TRUE.
            lvec( 1,k) = .TRUE.
            lvec( 2,k) = .TRUE.
            nder(k) = 5
         ENDIF
!           p_d     = p(k) + delta
      ENDIF
!
!     Populate trial values for DIFFEV
!
      do l=-2,2
         if(lvec(l,k)) then
            run_mpi_senddata%children = run_mpi_senddata%children + 1    ! at p(k) + delta
            pop_t(k, run_mpi_senddata%children) = dvec(l,k)
            kid_is(1,run_mpi_senddata%children) =  k
            kid_is(2,run_mpi_senddata%children) =  l
            para_has(0,k) = para_has(0,k) + 1
            para_has(para_has(0,k),k) = run_mpi_senddata%children
         endif
      enddo
   endif cond_is_deriv
enddo loop_npara_init
pop_n = run_mpi_senddata%children
pop_c = run_mpi_senddata%children
pop_gen = 0
lastgen = -1
!do k=1, run_mpi_senddata%children
   if(refine_log) then
      write(run_mpi_senddata%out,'(a)') 'DISCUS_SUITE_DERIVATIVES/LOGFILE'
      run_mpi_senddata%out_l = 32
   else
      write(run_mpi_senddata%out,'(a)') '/dev/null/'
      run_mpi_senddata%out_l = 10
   endif
!enddo
!
call random_current(run_mpi_senddata%nseeds, run_mpi_senddata%seeds)
run_mpi_senddata%l_get_state=1
call run_mpi_master
ier_num = 0
ier_typ = 0
ALLOCATE(refine_tttt(1:data_dim(1), 1:data_dim(2), 1:data_dim(3), -2:2))
!
! Here we will place call to run_mpi_master
!
if(1==2) then   ! DUMMY OLD STUFF
   DO l = -2, 2
      IF(lvec(l,k)) THEN
         CALL refine_set_param(NPARA, par_names(k), k, (dvec(l,k)) )  ! Set modified value
         CALL refine_restore_seeds
         CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
                              data_dim, refine_temp)
         IF(ier_num /= 0) THEN
            DEALLOCATE(refine_tttt)
!           exit cond_is_deriv_2
            return
         ENDIF
         do iiz = 1, data_dim(3)
            DO iiy=1, data_dim(2)
               DO iix=1, data_dim(1)
                  refine_tttt(iix,iiy,iiz,l) =  refine_temp(iix,iiy,iiz)
               ENDDO
            ENDDO
         enddo
      ENDIF
   ENDDO
   call refine_rvalue_tttt(data_dim, refine_tttt, dvec(:,k))
   CALL refine_set_param(NPARA, par_names(k), k, p(k))  ! Return to original value
endif
!
! Now evaluate all derivatives
!
loop_npara_eval: do k=1, NPARA
   cond_is_deriv_2: if(gl_is_der(k)) then
      continue
   else cond_is_deriv_2
!
! Load derivatives from disk
!
   j2 = 38
   do i=1, para_has(0,k)           ! Loop over all derivativs for this parameter
      j = para_has(i,k)            ! Current REF_KID
      l = kid_is(2,j)
      string =  ' '
      write(string,'(a,i4.4)') 'h5, DISCUS_SUITE_DERIVATIVES/data.', j
      call refine_kupl_last(0)
      call do_load(string, j2, .FALSE.)
      call refine_load_calc(data_dim, refine_temp)
      refine_tttt(:,:,:,l) = refine_temp(:,:,:)
   enddo
!
   IF(nder(k)==5) THEN             ! Got all five  points for derivative
      do iiz=1, data_dim(3)
         DO iiy=1, data_dim(2)
            DO iix=1, data_dim(1)
               refine_derivs(iix, iiy, iiz, k) = (-1.0*refine_tttt(iix,iiy,iiz, 2)   &
                                                  +8.0*refine_tttt(iix,iiy,iiz, 1)   &
                                                  -8.0*refine_tttt(iix,iiy,iiz,-1)   &
                                                  +1.0*refine_tttt(iix,iiy,iiz,-2))/ &
                                                 (12.*delta(k))
            ENDDO
         ENDDO
      enddo
   ELSEIF(nder(k)==3) THEN             ! Got all three points for derivative
      xmat(:,1) =  1.0
      xmat(1,2) =  dvec(0,k) !p(k)
      IF(lvec(2,k)) THEN              ! +delta, + 2delta
         j2 =  1
         j3 =  2
      ELSEIF(lvec(-2,k)) THEN         ! -delta, - 2delta
         j2 = -1
         j3 = -2
      ELSE
         j2 = -1
         j3 =  1
      ENDIF
      xmat(2,2) =  dvec(j2,k)       ! p(k) - delta
      xmat(3,2) =  dvec(j3,k)       ! p(k) + delta
      xmat(2,3) = (dvec(j2,k))**2   ! (p(k) - delta ) **2
      xmat(3,3) = (dvec(j3,k))**2   ! (p(k) + delta ) **2
      CALL matinv3(xmat, imat)
      IF(ier_num/=0) then
         ier_msg(1) = 'Error determining derivative '
         write(ier_msg(2), '(a,i3)') 'At derivative ',k
         exit cond_is_deriv_2
      endif
      do iiz=1, data_dim(3)
         DO iiy=1, data_dim(2)
            DO iix=1, data_dim(1)
!
!              Derivative is calculated as a fit of a parabola at P, P+delta, P-delta
               yvec(1) = refine_calc  (iix, iiy, iiz)
               yvec(2) = refine_tttt  (iix, iiy, iiz, j2)
               yvec(3) = refine_tttt  (iix, iiy, iiz, j3)
               avec = MATMUL(imat, yvec)
!
               refine_derivs(iix, iiy, iiz, k) = avec(2) + 2.*avec(3)*p(k)
            ENDDO
         ENDDO
      enddo
!
   ELSE
      ier_num = -9
      ier_typ = 6
      ier_msg(1) = par_names(k)
      DEALLOCATE(refine_tttt)
      exit cond_is_deriv_2
   ENDIF
!
ENDIF cond_is_deriv_2
enddo loop_npara_eval
!
DEALLOCATE(refine_tttt)
deallocate(refine_temp)
!
end subroutine refine_calc_deriv_mpi
!
!*******************************************************************************
!
SUBROUTINE refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, refine_params, p, &
                        dimen, array)
!
!   Runs the user defined macro to calculate the cost value
!
USE refine_fit_set_sub_mod
!USE refine_setup_mod
!
USE doact_mod
USE do_if_mod
USE do_set_mod
USE errlist_mod
USE lib_errlist_func
USE lib_macro_func
USE precision_mod
USE prompt_mod
USE set_sub_generic_mod
USE sup_mod
!
USE global_data_mod
!
IMPLICIT NONE
!
INTEGER                             , INTENT(IN)    :: MAXP          ! Parameter array size
CHARACTER(LEN=*)                    , INTENT(INOUT) :: refine_mac    ! Macro name
INTEGER                             , INTENT(INOUT) :: refine_mac_l  ! length of macro name
INTEGER                             , INTENT(IN)    :: NPARA         ! Nmber of parameters
INTEGER                             , INTENT(IN)    :: kupl_last     ! Last dat set needed in KUPLOT
CHARACTER(LEN=*)  , DIMENSION(MAXP ), INTENT(IN)    :: refine_params ! parameter names
REAL(kind=PREC_DP), DIMENSION(MAXP ), INTENT(IN)    :: p             ! parameter values
INTEGER           , DIMENSION(3    ), INTENT(IN)    :: dimen         ! Data array dimensions
REAL(kind=PREC_DP), DIMENSION(dimen(1), dimen(2), dimen(3)), INTENT(OUT) :: array  ! The actual data

CHARACTER(LEN=PREC_STRING) :: string           ! dumy string variable
CHARACTER(len=PREC_STRING) :: zeile            ! dummy string
CHARACTER(LEN=4   ) :: befehl           ! Command verb
INTEGER             :: lbef             ! length   of        befehl
INTEGER             :: lp               ! length   of        befehl
INTEGER             :: length           ! Length of a parameter name
INTEGER             :: ier_number, ier_type ! Local error status
LOGICAL             :: lend             ! TRUE if macro is finished
LOGICAL             :: l_prompt_restore
!
INTERFACE
   SUBROUTINE refine_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE refine_mache_kdo
END INTERFACE
!
CALL file_kdo(refine_mac, refine_mac_l)
lmacro_close = .FALSE.       ! Do not close macros in do-loops, return instead
!
CALL refine_kupl_last(kupl_last)   ! Set last data set needed in KUPLOT
!                            ! Turn off output
lturn_off = .true.
!lturn_off = .false.
l_prompt_restore = .FALSE.
IF(output_status == OUTPUT_SCREEN .AND.lturn_off) THEN
   l_prompt_restore = .TRUE.
   string = 'prompt, off, off, save'
   length = 22
   CALL do_set(string, length)
ENDIF
!
CALL refine_fit_set_sub
main: DO
   CALL get_cmd (string, length, befehl, lbef, zeile, lp, prompt)
   IF (ier_num == 0) THEN
      IF (string == ' '.OR.string (1:1)  == '#' .OR. string=='!') CYCLE main
      IF (befehl (1:3)  == 'do '.OR.befehl (1:2)  == 'if') THEN
         CALL do_loop (string, lend, length) !, kuplot_mache_kdo)
         IF(ier_num/=0) EXIT main
      ELSE
         CALL refine_fit_mache_kdo(string, lend, length)
      ENDIF
      IF(lend .OR. ier_num/=0) EXIT main
   ELSE
      EXIT main
   ENDIF
ENDDO main
!
ier_number = 0
ier_type   = 0
IF(ier_num==-9 .AND. ier_typ==1) THEN   ! Error condition in macro
   ier_number = -6
   ier_type   =  ER_APPL
   ier_msg(1) = 'check performance of user macro'
   ier_num = 0
   ier_typ = 0
ELSEIF(ier_num /= 0) THEN         ! Back up error status
   ier_number = ier_num
   ier_type   = ier_typ
   CALL errlist
   ier_number = -6
   ier_type   =  ER_APPL
   ier_msg(1) = 'Check performance of user macro'
   ier_num = 0
   ier_typ = 0
ENDIF
CALL macro_terminate
p_mache_kdo   => refine_mache_kdo
lmacro_close = .TRUE.       ! Do close macros in do-loops
!
!
IF(ier_number == 0 .AND. IER_NUM == 0) THEN
   IF(gl_is_der(0)) THEN       ! Calculated data were set into globalk array
!     dimen_gl(1:2) = dimen(1:2)
!     dimen_gl(3)   = 1
!     dimen_gl(4)   = NPARA
      CALL gl_get_data(0, dimen(1),dimen(2), dimen(3), array(1:dimen(1), 1:dimen(2), 1:dimen(3)))
   ELSE
      CALL refine_load_calc(dimen, array)
   ENDIF
   IF(ier_num/=0) RETURN
ENDIF
!
IF(l_prompt_restore) THEN
   string = 'prompt, on,on'
   length = 13
   CALL do_set(string, length)
ENDIF
!
IF(ier_number /= 0) THEN    ! If necessary restore error status
!  CALL refine_set_sub
   ier_num = ier_number
   ier_typ = ier_type
ENDIF
CALL refine_fit_un_sub
!
END SUBROUTINE refine_macro
!
!*******************************************************************************
!
SUBROUTINE refine_do_plot(plmac)
!
USE do_set_mod
USE set_sub_generic_mod
!
USE class_macro_internal
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*   ), INTENT(IN) :: plmac     ! optional plot macro
!
CHARACTER(LEN=PREC_STRING)                    :: string
CHARACTER(LEN=MAX(PREC_STRING,LEN(plmac)+14)) :: kline
INTEGER                                       :: lcomm, length
!
string = 'prompt, off, off, save'
length = 22
CALL do_set(string, length)
kline = 'kuplot -macro ' // plmac(1:LEN_TRIM(plmac))
lcomm = LEN_TRIM(kline)
CALL p_branch (kline, lcomm, .FALSE., 0     )
string = 'prompt, on,on'
length = 13
CALL do_set(string, length)
!
END SUBROUTINE refine_do_plot
!
!*******************************************************************************
!
SUBROUTINE refine_kupl_last(kupl_last)   ! Set last data set needed in KUPLOT
!
USE kuplot_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: kupl_last     ! Last dat set needed in KUPLOT
!
iz = kupl_last + 1
!
END SUBROUTINE refine_kupl_last
!
!*******************************************************************************
!
SUBROUTINE refine_load_calc(dimen, array)
!
! Transfer calculated data from last KUPLOT data set
! 
USE kuplot_mod
!
USE errlist_mod
use lib_data_struc_h5
use lib_ik_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER           , DIMENSION(3)                 , INTENT(IN)  :: dimen  ! Array dimensions
REAL(kind=PREC_DP), DIMENSION(dimen(1), dimen(2), dimen(3)), INTENT(OUT) :: array  ! The actual data
!
INTEGER :: ndata    ! Data set to be loaded
INTEGER :: iix, iiy, iiz   ! loop variables
!
IF(iz<1) THEN
   ier_num = -3
   ier_typ = ER_APPL
   ier_msg(1) = 'KUPLOT does not have any data sets'
   ier_msg(2) = 'Load calculated data into KUPLOT, '
   ier_msg(3) = 'or remove a ''reset'' from KUPLOT section'
   RETURN
ENDIF
!
ndata = iz-1                           ! Always last data set in KUPLOT
!
!IF(lni(ndata)) THEN                    ! 2D data set
if(ku_ndims(ndata)==3) then            ! 3D data set
   continue
!
!write(*,*) ' DIMEN', dimen
   call data2local(ndata   , ier_num, ier_typ, ik1_node_number, ik1_infile, ik1_data_type,    &
        ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
        ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, ik1_corners, ik1_vectors, ik1_a0, ik1_win,  &
        ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
        ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!write(*,*) ' IK1_DIM ', ik1_dims
!write(*,*) ' ik1_x   ', ik1_x(lbound(ik1_x)), ik1_x(ubound(ik1_x))
!write(*,*) ' ik1_y   ', ik1_y(lbound(ik1_y)), ik1_y(ubound(ik1_y))
!write(*,*) ' ik1_z   ', ik1_z(lbound(ik1_z)), ik1_z(ubound(ik1_z))
   array = ik1_data
   deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   if(allocated(ik1_x    )) deallocate(ik1_x)
   if(allocated(ik1_y    )) deallocate(ik1_y)
   if(allocated(ik1_z    )) deallocate(ik1_z)

elseif(ku_ndims(ndata)==2) then        ! 2D data set
   IF(dimen(1)/=nx(ndata).OR.dimen(2)/=ny(ndata)) THEN
      ier_num = -4
      ier_typ = ER_APPL
      ier_msg(1) = 'Check x/y limits for data calculation '
      RETURN
   ENDIF
!
   iiz = 1
   DO iiy=1,dimen(2)
      DO iix=1,dimen(1)
         array(iix,iiy,iiz)  = z(offz(ndata - 1) + (iix - 1)*ny(ndata) + iiy)
      ENDDO
   ENDDO
elseif(ku_ndims(ndata)==1) then       ! 1D data seelseif(ku_ndims(ndata)==1) then       ! 1D data set
   IF(dimen(1)/=lenc(ndata)) THEN
      ier_num = -4
      ier_typ = ER_APPL
      ier_msg(1) = 'Check x/y limits for data calculation '
      RETURN
   ENDIF
!
   iiz = 1
   DO iix=1,dimen(1)
      array(iix, 1, 1)   = y(offxy(ndata - 1) + iix)
   ENDDO
ENDIF
!
END SUBROUTINE refine_load_calc
!
!*******************************************************************************
!
SUBROUTINE refine_mrq(linit, MAXP, NPARA, ncycle, kupl_last, par_names, data_dim, &
                      data_data, data_sigma, data_x, data_y, data_z, conv_status, &
                      conv_dp_sig, conv_dchi2, conv_chi2, conv_conf, conv_lambda, lconvergence,&
                      lconv, chisq, conf, lamda_fin, lamda_s, lamda_d, lamda_u,  rval,   &
                      rexp, p, prange, p_shift, p_nderiv, dp, cl, alpha, beta,    &
                      ref_do_plot, plmac)
!+                                                                      
!   This routine runs the refinement cycles, interfaces with the 
!   Levenberg-Marquardt routine modified after Numerical Recipes
!-                                                                      
!
USE refine_set_param_mod
!
USE errlist_mod
USE ber_params_mod
USE calc_expr_mod
USE gamma_mod
USE param_mod 
USE precision_mod
USE prompt_mod 
!                                                                       
IMPLICIT NONE
!
LOGICAL                                                , INTENT(IN)    :: linit       ! Initialize MRQ
INTEGER                                                , INTENT(IN)    :: MAXP        ! Parameter array size
INTEGER                                                , INTENT(IN)    :: NPARA       ! Number of parameters
INTEGER                                                , INTENT(IN)    :: ncycle      ! maximum cycle number
INTEGER                                                , INTENT(IN)    :: kupl_last   ! Last KUPLOT DATA that are needed
CHARACTER(LEN=*)  , DIMENSION(MAXP)                    , INTENT(IN)    :: par_names   ! Parameter names
INTEGER           , DIMENSION(3)                       , INTENT(IN)    :: data_dim    ! Data array dimensions
REAL(kind=PREC_DP), DIMENSION(data_dim(1), data_dim(2), data_dim(3)), INTENT(IN)    :: data_data   ! Data array
REAL(kind=PREC_DP), DIMENSION(data_dim(1), data_dim(2), data_dim(3)), INTENT(IN)    :: data_sigma ! Data sigmas
REAL(kind=PREC_DP), DIMENSION(data_dim(1))             , INTENT(IN)    :: data_x      ! Data coordinates x
REAL(kind=PREC_DP), DIMENSION(data_dim(2))             , INTENT(IN)    :: data_y      ! Data coordinates y
REAL(kind=PREC_DP), DIMENSION(data_dim(3))             , INTENT(IN)    :: data_z      ! Data coordinates y
LOGICAL                                                , INTENT(IN)    :: conv_status ! Do convergence test
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: conv_dp_sig ! Max parameter shift
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: conv_dchi2  ! Max Chi^2     shift
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: conv_chi2   ! Min Chi^2     value 
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: conv_conf   ! Min confidence level
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: conv_lambda ! Max lambda    value 
LOGICAL                                                , INTENT(INOUT) :: lconvergence ! Convergence criteria were reached
LOGICAL           , DIMENSION(4)                       , intent(out)   :: lconv       ! Convergence criteria
REAL(kind=PREC_DP)                                     , INTENT(OUT)   :: chisq       ! Chi^2 value
REAL(kind=PREC_DP)                                     , INTENT(OUT)   :: conf        ! Confidence test from Gamma function
REAL(kind=PREC_DP)                                     , INTENT(OUT)   :: lamda_fin   ! Final Marquardt lamda
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: lamda_s    ! Start value lamda
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: lamda_d    ! Multiplier_down
REAL(kind=PREC_DP)                                     , INTENT(IN)    :: lamda_u    ! Multiplier_up
REAL(kind=PREC_DP)                                     , INTENT(OUT)   :: rval        ! Weighted R-value
REAL(kind=PREC_DP)                                     , INTENT(OUT)   :: rexp        ! Expected R-value
!
REAL(kind=PREC_DP), DIMENSION(MAXP)                    , INTENT(INOUT) :: p         ! Parameter array
REAL(kind=PREC_DP), DIMENSION(MAXP,2)                  , INTENT(INOUT) :: prange    ! Parameter range
REAL(kind=PREC_DP), DIMENSION(MAXP )                   , INTENT(IN)    :: p_shift   ! Parameter shift for deriv
INTEGER           , DIMENSION(MAXP )                   , INTENT(IN)    :: p_nderiv  ! Number of derivative points needed
REAL(kind=PREC_DP), DIMENSION(MAXP)                    , INTENT(INOUT) :: dp        ! Parameter sigmas
REAL(kind=PREC_DP), DIMENSION(NPARA, NPARA)            , INTENT(INOUT) :: cl        ! Covariance matrix
REAL(kind=PREC_DP), DIMENSION(NPARA, NPARA)            , INTENT(INOUT) :: alpha
REAL(kind=PREC_DP), DIMENSION(NPARA       )            , INTENT(INOUT) :: beta
LOGICAL                                                , INTENT(IN)    :: ref_do_plot ! Do plot yes/no
CHARACTER(LEN=*   )                                    , INTENT(IN)    :: plmac     ! optional plot macro
!
INTEGER :: icyc
INTEGER :: k
INTEGER :: last_i, prev_i
REAL(kind=PREC_DP), SAVE    :: alamda
!
INTEGER :: MAXW
INTEGER :: ianz
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(par_names))), DIMENSION(:), ALLOCATABLE :: cpara
INTEGER           , DIMENSION(:), ALLOCATABLE :: lpara
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: werte
!
LOGICAL                            :: lsuccess  ! Ture if cycle improved CHI^2
REAL(kind=PREC_DP), DIMENSION(0:3) :: last_chi
REAL(kind=PREC_DP), DIMENSION(0:3) :: last_shift
REAL(kind=PREC_DP), DIMENSION(0:3) :: last_conf 
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: last_p 
INTEGER :: last_ind
!
IF(linit) THEN        ! INITIALIZE MRQ 
   alamda     = -0.01    ! Negative lamda initializes MRQ routine
   rval       = 0.0
   rexp       = 0.0
   chisq      = 0.0
   last_chi(:)   = HUGE(0.0)
   last_shift(:) = -1.0
   last_conf(:)  = HUGE(0.0)
   last_i        = 0
ENDIF
!
! Set initial parameter values
!
MAXW = NPARA
ALLOCATE(cpara(MAXW))
ALLOCATE(lpara(MAXW))
ALLOCATE(werte(MAXW))
DO k = 1, NPARA                  ! Copy parameter names to prepare calculation
   cpara(k) = par_names(k)
   lpara(k) = LEN_TRIM(par_names(k))
ENDDO
ianz = NPARA
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO k = 1, NPARA                  ! Copy values of the parameters
   p(k) = werte(k)
ENDDO
!
ALLOCATE(last_p(0:2,NPARA))
last_p(2,:) = p(1:NPARA)
prev_i = 2
!
WRITE(output_io,'(a,10x,a,7x,a,6x,a)') 'Cyc Chi^2/(N-P)   MAX(dP/sig) Par   Conf',&
                                       'Lambda', 'wRvalue', 'Rexp'
!
lconvergence = .FALSE.
icyc = 0
cycles:DO
   CALL mrqmin(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, p, NPARA, &
               par_names, prange, p_shift, p_nderiv, kupl_last, cl, alpha, beta, chisq, alamda,     &
               lamda_s, lamda_d, lamda_u, lsuccess, dp, rval, rexp)
!
   IF(ier_num/=0) EXIT cycles
   IF(lsuccess) THEN ! .OR. rval == 0.0) THEN
      CALL refine_rvalue(rval, rexp, NPARA)
      IF(ier_num/=0) then
         ier_msg(1)=' Error calculating R-value '
         EXIT cycles
      ENDIF
   ENDIF
   last_shift(:) = -1.0                                             ! Set parameter shift / sigma to negative
   DO k=1, NPARA
      CALL refine_set_param(NPARA, par_names(k), k, p(k))   ! Updata parameters
      IF(ier_num/=0) then
         ier_msg(1)=' Error setting parameters   '
         EXIT cycles
      endif
      last_p(last_i,k) = p(k)
      IF(cl(k,k)/= 0.0) THEN
!        last_shift(last_i) = MAX(last_shift(last_i), ABS((p(k)-last_p(prev_i,k))/SQRT(cl(k,k)))) !/SQRT(cl(k,k))))
         IF(ABS((p(k)-last_p(prev_i,k))/SQRT(cl(k,k))) > last_shift(last_i)) THEN
            last_ind = k
            last_shift(last_i) = ABS((p(k)-last_p(prev_i,k))/SQRT(cl(k,k)))
         ENDIF
      ENDIF
   ENDDO
   conf = gammaq((data_dim(1)*data_dim(2)*data_dim(3)-2)*0.5D0, chisq*0.5D0)
   last_chi( last_i) = chisq/(data_dim(1)*data_dim(2)*data_dim(3)-NPARA)        ! Store Chi^2/(NDATA-NPARA)
   last_conf(last_i) = conf
   WRITE(*,'(i3,2g13.5e2,i4,f9.4,4x,4g13.5e2)') icyc,chisq/(data_dim(1)*data_dim(2)*data_dim(3)-NPARA), &
         last_shift(last_i), last_ind, conf, alamda, rval, rexp
!
!  CALL refine_rvalue(rval, rexp, NPARA)
!  DO k = 1, NPARA
!     dp(k) = SQRT(ABS(cl(k,k)))
!  ENDDO
   CALL refine_best(rval)                                              ! Write best macro
!
   IF(lsuccess) THEN
      IF(ref_do_plot) THEN
         l_plot_status = .TRUE.
         CALL refine_do_plot(plmac)
         l_plot_status = .FALSE.
         IF(ier_num/=0) then
            ier_msg(1)=' Error in interactive plot'
            EXIT cycles
         endif
      ENDIF
   ENDIF
!
   IF(conv_status) THEN 
   lconv(1) = (ABS(last_chi( last_i)-last_chi(prev_i))<conv_dchi2 .AND.   &
                   last_shift(last_i) < conv_dp_sig               .AND.   &
                   last_conf(last_i)  > conv_conf)
   lconv(2) = (last_shift(last_i)>0.0 .AND.                               &
               ABS(last_chi(last_i)-last_chi(prev_i))<conv_dchi2)
   lconv(3) = last_chi(last_i)  < conv_chi2
   lconv(4) = alamda            > conv_lambda
!  lconv(2) =   last_chi( last_i)  < conv_chi2
   IF(lconv(1) .OR. lconv(2) .OR. lconv(3) .or. lconv(4)) THEN
!  IF(ABS(last_chi( last_i)-last_chi(prev_i))<conv_dchi2 .AND.   &
!     last_shift(last_i) < conv_dp_sig                   .AND.   &
!     last_conf(last_i)  > conv_conf                     .OR.    &
!     last_chi( last_i)  < conv_chi2                   ) THEN
      WRITE(output_io, '(/,a)') 'Convergence reached '        
      lconvergence = .TRUE.
      EXIT cycles
   ENDIF
   ENDIF
   IF(alamda > 0.5*HUGE(0.0)) THEN
      lconvergence = .FALSE.
      EXIT cycles
   ENDIF
   prev_i = last_i
   last_i = MOD(last_i+1, 3)
   icyc   = icyc + 1
   IF(icyc==ncycle) EXIT cycles
ENDDO cycles
!
IF(.NOT. lconvergence) THEN
   IF(icyc==ncycle) THEN
      WRITE(output_io,'(/,a)') ' Maximum cycle is reached without convergence'
   ELSEIF(alamda > 0.5*HUGE(0.0)) THEN
      WRITE(output_io,'(/,a)') ' MRQ Parameter reached infinity'
   ENDIF
ENDIF
!
! Make a final call with lamda = 0 to ensure that model is up to date
!
IF(ier_num==0) THEN
   lamda_fin = alamda
   alamda = 0.0
!
   CALL mrqmin(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, p, NPARA, &
               par_names, prange, p_shift, p_nderiv, kupl_last, cl, alpha, beta, chisq, alamda,     &
               lamda_s, lamda_d, lamda_u, lsuccess, dp, rval, rexp)
   CALL refine_rvalue(rval, rexp, NPARA)
!   DO k = 1, NPARA
!      dp(k) = SQRT(ABS(cl(k,k)))
!   ENDDO
   CALL refine_best(rval)                                              ! Write best macro
   alamda = lamda_fin
ENDIF
!
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
DEALLOCATE(last_p)
!
END SUBROUTINE refine_mrq
!
!*******************************************************************************
!
SUBROUTINE mrqmin(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, a, NPARA, &
    par_names, &
    prange, p_shift, p_nderiv, kupl_last, covar, alpha, beta, chisq, alamda, lamda_s, lamda_d, lamda_u, &
    lsuccess, dp, rval, rexp)
!
!  Least squares routine adapted from Numerical Recipes Chapter 14
!  p 526
!
USE refine_constraint_mod
USE errlist_mod
USE gaussj_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER                                               , INTENT(IN)    :: MAXP        ! Array size for parameters
INTEGER           , DIMENSION(3)                      , INTENT(IN)    :: data_dim    ! no of ata points along x, y
REAL(kind=PREC_DP), DIMENSION(data_dim(1),data_dim(2),data_dim(3)), INTENT(IN)    :: data_data   ! Data values
REAL(kind=PREC_DP), DIMENSION(data_dim(1),data_dim(2),data_dim(3)), INTENT(IN)    :: data_sigma ! Data sigmas
REAL(kind=PREC_DP), DIMENSION(data_dim(1)           ) , INTENT(IN)    :: data_x      ! Actual x-coordinates
REAL(kind=PREC_DP), DIMENSION(data_dim(2)           ) , INTENT(IN)    :: data_y      ! Actual y-coordinates
REAL(kind=PREC_DP), DIMENSION(data_dim(3)           ) , INTENT(IN)    :: data_z      ! Actual y-coordinates
INTEGER                                               , INTENT(IN)    :: NPARA       ! number of parameters
REAL(kind=PREC_DP), DIMENSION(MAXP, 2              )  , INTENT(IN)    :: prange      ! Allowed parameter range
REAL(kind=PREC_DP), DIMENSION(MAXP )                  , INTENT(IN)    :: p_shift   ! Parameter shift for deriv
INTEGER           , DIMENSION(MAXP )                  , INTENT(IN)    :: p_nderiv  ! Number of derivative points needed
INTEGER                                               , INTENT(IN)    :: kupl_last   ! Last KUPLOT DATA that are needed
REAL(kind=PREC_DP), DIMENSION(MAXP)                   , INTENT(INOUT) :: a           ! Parameter values
CHARACTER(LEN=*)  , DIMENSION(MAXP)                   , INTENT(IN)    :: par_names   ! Parameter names
REAL(kind=PREC_DP), DIMENSION(NPARA, NPARA)           , INTENT(OUT)   :: covar       ! Covariance matrix
REAL(kind=PREC_DP), DIMENSION(NPARA, NPARA)           , INTENT(INOUT) :: alpha       ! Temp arrays
REAL(kind=PREC_DP), DIMENSION(NPARA     )             , INTENT(INOUT) :: beta        ! Temp arrays
REAL(kind=PREC_DP)                                    , INTENT(INOUT) :: chisq       ! Chi Squared
REAL(kind=PREC_DP)                                    , INTENT(INOUT) :: alamda      ! Levenberg parameter
REAL(kind=PREC_DP)                                    , INTENT(IN)    :: lamda_s    ! Start value for lamda
REAL(kind=PREC_DP)                                    , INTENT(IN)    :: lamda_d    ! Multiplier_down
REAL(kind=PREC_DP)                                    , INTENT(IN)    :: lamda_u    ! Multiplier_up
LOGICAL                                               , INTENT(OUT)   :: lsuccess   ! True if improved
REAL(kind=PREC_DP), DIMENSION(MAXP)                   , INTENT(OUT)   :: dp          ! Parameter sigmas
REAL(kind=PREC_DP)                                    , INTENT(OUT)   :: rval        ! Initial Rvalue
REAL(kind=PREC_DP)                                    , INTENT(OUT)   :: rexp        ! Initial expected Rvalue
!
INTEGER                         :: j     = 0
INTEGER                         :: k     = 0
REAL(kind=PREC_DP), DIMENSION(MAXP)        :: atry != 0.0
REAL(kind=PREC_DP), DIMENSION(NPARA)       :: da != 0.0
REAL(kind=PREC_DP), DIMENSION(NPARA,NPARA) :: cl != 0.0
REAL(kind=PREC_DP)                            :: ochisq = 0.0
logical           , dimension(npara)       :: lderiv_ok
logical                                    :: l_none
LOGICAL, PARAMETER              :: LDERIV = .TRUE.
LOGICAL, PARAMETER              :: NDERIV = .FALSE.
!
!real(kind=PREC_DP) :: e_del           ! Change of Chi^2
!real(kind=PREC_DP) :: e_ran           ! Value for Boltzmann probability
!real(kind=PREC_DP) :: mo_kt           ! Pseudo Temperature
!real(kind=PREC_DP) :: r1              ! Random number
!integer, save :: icyy = 0
!
lderiv_ok = .true.                    ! Assume all derivatives went fine
ochisq = chisq
IF(alamda < 0) THEN                   ! Initialization
   alamda = lamda_s
   CALL mrqcof(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, a, &
               NPARA, par_names, prange, p_shift, p_nderiv, kupl_last, alpha, beta, &
               chisq, LDERIV) !, funcs)
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
   CALL refine_rvalue(rval, rexp, NPARA)          ! Call initial R-value
   ochisq = chisq
   DO j=1, NPARA
      IF(prange(j,1)<=prange(j,2)) THEN
         atry(j) = MIN(prange(j,2),MAX(prange(j,1),a(j)))
      ELSE
         atry(j) = a(j)
      ENDIF
   ENDDO
ENDIF
!
DO j=1, NPARA
   DO k=1, NPARA
      covar(j,k) = alpha(j,k)
   ENDDO
   covar(j,j) = alpha(j,j)*(1.0+alamda)       ! Adjust by Levenberg-Marquard algorithm
   da(j) = beta(j)
ENDDO
l_none = .true.
do j=1, npara
   if(covar(j,j) == 0) then                   ! Failed to get derivative for this param
      covar(j,j) = 1.0D0
      lderiv_ok(j) = .false.
   else
      lderiv_ok(j) = .true.
      l_none       = .false.
   endif
enddo
if(l_none) then
   lsuccess = .FALSE.
   ier_msg(1) = 'All derivatives are zero; check user macro'
   ier_msg(2) = 'try macro with different parameter values' 
   ier_num = -1
   ier_typ =  5
   return
endif
!
CALL gausj(NPARA, covar, da)
IF(ier_num/=0) THEN
   ier_msg(1)  = ' Error in Gauss-Jordan function '
   ier_msg(2)  = ' Are parameters correlated by 100%?'
   RETURN
ENDIF
!
cl(:,:) = covar(:,:)       ! Back up of covariance
DO j=1,npara
  dp(j) = SQRT(ABS(cl(j,j)))
ENDDO
IF(alamda==0) THEN
!
!  Last call to update the structure, no derivative is needed
!
   CALL mrqcof(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, a, &
               NPARA, par_names, prange, p_shift, p_nderiv, kupl_last, covar, da, chisq, NDERIV) !, funcs)
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
   covar(:,:) = cl(:,:)       ! Restore    covariance
   RETURN
ENDIF
!
!open(78,file='shift.dat', status='unknown', position='append')
!if(icyy==0) write(78, '(15a)') 'cyc', (par_names(j),j=1,npara)
!icyy = icyy+1
!write(78,'(i3,14g15.6e3)') icyy,da(1:npara)
!close(78)
if(data_dim(2)==1 .and. data_dim(3)==1) then     ! Only for 1D data
   CALL refine_constrain_fwhm_in(MAXP, NPARA, a, da, data_x(1), data_x(data_dim(1)), ier_num, ier_typ)
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
   CALL refine_constrain_eta_in(MAXP, NPARA, a, da, data_x(1), data_x(data_dim(1)), ier_num, ier_typ)
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
endif
DO j=1,NPARA
   IF(prange(j,1)<=prange(j,2)) THEN
      atry(j) = MIN(prange(j,2),MAX(prange(j,1),a(j)+da(j)))
   ELSE
      atry(     (j)) = a(     (j)) + da(j)
   ENDIF
ENDDO
!
CALL mrqcof(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, atry, &
            NPARA, par_names, prange, p_shift, p_nderiv, kupl_last, covar, da, chisq, LDERIV) !, funcs)
IF(ier_num/=0) THEN
   RETURN
ENDIF
!
!do j=1, npara
!write(*,'(a5,3(f20.8,2x))') 'ond ', a(j), atry(j), -(a(j)-atry(j))
!enddo
!
!write(*,*) ' PAR ', params(1:npara)
!write(*,*) ' dP  ', da    (1:npara)
!write(*,*) 'c, d, u', alamda, lamda_d, lamda_u
IF(chisq < ochisq) THEN            ! Success, accept solution
   alamda = lamda_d*alamda         ! Decrease alamda by facror lamda_u (down)
   ochisq = chisq
   DO j=1,NPARA
      DO k=1,NPARA
         alpha(j,k) = covar(j,k)
      ENDDO
      beta(j) = da(j)
      a(     (j)) = atry(     (j)) ! Adjust parameter values to trial ones
   ENDDO
   lsuccess = .TRUE.
ELSE                               ! Failure, reject, keep old params and 
!  mo_kt = 30.0D0
!  r1 = (data_dim(1)*data_dim(2)*data_dim(3)-npara)
!  e_del  = (chisq - ochisq)/r1    ! Calculate change in Chi^2
!  e_ran = exp ( - e_del / mo_kt)
!  e_ran = e_ran / (1 + e_ran)
!  write(*,'(a, 5f10.2)') ' Failure ', chisq/r1, ochisq/r1, e_del, mo_kt, e_ran
!  CALL RANDOM_NUMBER(r1)
!  if(e_ran > r1) then
!     ochisq = chisq
!     DO j=1,NPARA
!        DO k=1,NPARA
!           alpha(j,k) = covar(j,k)
!        ENDDO
!        beta(j) = da(j)
!        a(     (j)) = atry(     (j)) ! Adjust parameter values to trial ones
!     ENDDO
!     lsuccess = .TRUE.
!     write(*,'(a)') ' ACCEPTED '
!     lsuccess = .TRUE.
!     alamda =  lamda_u*alamda        ! increment alamda by factor lamda_u (up)
!  else
!
   alamda =  lamda_u*alamda        ! increment alamda by factor lamda_u (up)
   chisq = ochisq                  ! Keep old chis squared
   lsuccess = .FALSE.
!  endif
ENDIF
!
!
END SUBROUTINE mrqmin
!
!*******************************************************************************
!
SUBROUTINE mrqcof(MAXP, data_dim, data_data, data_sigma, data_x, data_y, data_z, &
                  params, NPARA, par_names, prange, p_shift, p_nderiv, kupl_last,  &
                  alpha, beta, chisq, LDERIV)
!
! Modified after NumRec 14.4
! Version for 2-d data
!
USE errlist_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER                                               , INTENT(IN)  :: MAXP        ! Maximum parameter number
INTEGER           , DIMENSION(3)                      , INTENT(IN)  :: data_dim    ! Data dimensions
REAL(kind=PREC_DP), DIMENSION(data_dim(1),data_dim(2),data_dim(3)), INTENT(IN)  :: data_data   ! Observables
REAL(kind=PREC_DP), DIMENSION(data_dim(1),data_dim(2),data_dim(3)), INTENT(IN)  :: data_sigma ! Observables
REAL(kind=PREC_DP), DIMENSION(data_dim(1)           ) , INTENT(IN)  :: data_x      ! x-coordinates
REAL(kind=PREC_DP), DIMENSION(data_dim(2)           ) , INTENT(IN)  :: data_y      ! y-coordinates
REAL(kind=PREC_DP), DIMENSION(data_dim(3)           ) , INTENT(IN)  :: data_z      ! y-coordinates
INTEGER                                               , INTENT(IN)  :: NPARA       ! Number of refine parameters
REAL(kind=PREC_DP), DIMENSION(MAXP)                   , INTENT(IN)  :: params      ! Current parameter values
CHARACTER(LEN=*)  , DIMENSION(MAXP)                   , INTENT(IN)  :: par_names   ! Parameter names
REAL(kind=PREC_DP), DIMENSION(MAXP, 2              )  , INTENT(IN)  :: prange      ! Allowed parameter range
REAL(kind=PREC_DP), DIMENSION(MAXP )                  , INTENT(IN)  :: p_shift   ! Parameter shift for deriv
INTEGER           , DIMENSION(MAXP )                  , INTENT(IN)  :: p_nderiv  ! Number of derivative points needed
INTEGER                                               , INTENT(IN)  :: kupl_last   ! Last KUPLOT DATA that are needed
REAL(kind=PREC_DP), DIMENSION(NPARA, NPARA)           , INTENT(OUT) :: alpha
REAL(kind=PREC_DP), DIMENSION(NPARA)                  , INTENT(OUT) :: beta
REAL(kind=PREC_DP)                                    , INTENT(OUT) :: chisq
LOGICAL                                               , INTENT(IN)  :: LDERIV      ! Derivatives are needed?
!
!
INTEGER                  :: j
INTEGER                  :: k
INTEGER                  :: iix, iiy, iiz
!
REAL(kind=PREC_DP)                     :: xx    = 0.0     ! current x-coordinate
REAL(kind=PREC_DP)                     :: yy    = 0.0     ! current y-coordinate
REAL(kind=PREC_DP)                     :: zz    = 0.0     ! current y-coordinate
REAL(kind=PREC_DP)                     :: ymod  = 0.0     ! zobs(calc)
REAL(kind=PREC_DP), DIMENSION(1:NPARA) :: dyda            ! Derivatives
REAL(kind=PREC_DP)                     :: sig2i = 0.0     ! sigma squared
REAL(kind=PREC_DP)                     :: dy    = 0.0     ! Difference zobs(obs) - zobs(calc)
REAL(kind=PREC_DP)                     :: wt    = 0.0     ! Weight 
!
alpha(:,:) = 0.0
beta(:)    = 0.0
chisq      = 0.0
iiz = 1
!
loopix: do iix=1, data_dim(1)
   xx = data_x(iix)
   loopiy: do iiy=1, data_dim(2)
      yy = data_y(iiy)
   loopiz: do iiz=1, data_dim(3)
      zz = data_z(iiz)
!
      CALL refine_theory(MAXP, iix, iiy, iiz, xx, yy, zz, NPARA, params, par_names,   &
                         prange, p_shift, p_nderiv, &
                         data_dim, & ! data_data, data_sigma, data_x, data_y, &
                         kupl_last, &
                         ymod, dyda, LDERIV)
      IF(ier_num/=0) THEN
         ier_msg(1) = ' Error in theory function '
         RETURN
      ENDIF
!
      if(data_sigma(iix,iiy,iiz)>0.0D0) then
         sig2i =1./(data_sigma(iix,iiy,iiz)*data_sigma(iix,iiy,iiz))
      elseif(data_sigma(iix,iiy,iiz)<-1000.0D0) then
         cycle loopiz
      else
         sig2i = 1.0e-9
      endif
!write(*,'(3f6.1, 2f12.4)') xx,yy,zz, ymod, data_sigma(iix,iiy,iiz)
      dy = data_data(iix,iiy,iiz) - ymod
      DO j=1, NPARA
         wt = dyda(j)*sig2i
         DO k=1, j
            alpha(j,k) = alpha(j,k) + wt*dyda(k)
         ENDDO
         beta(j) = beta(j) + dy*wt
      ENDDO
      chisq = chisq + dy*dy*sig2i
!
   ENDDO loopiz
   ENDDO loopiy
ENDDO loopix
!write(*,*) ' BETA ', beta(1:NPARA)
!
DO j=2, NPARA                   ! Fill in upper diagonal
   DO k=1,j-1
      alpha(k,j) = alpha(j,k)
   ENDDO
ENDDO
!write(*,*) ' ALPHA 1: ', alpha(1,:)
!write(*,*) ' ALPHA 2: ', alpha(2,:)
!write(*,*) ' ALPHA 3: ', alpha(3,:)
!
END SUBROUTINE mrqcof
!
!*******************************************************************************
!
SUBROUTINE refine_rvalue(rval, rexp, npar)
!
USE refine_data_mod
!
use lib_data_types_mod
use precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP)   , INTENT(OUT) :: rval
REAL(kind=PREC_DP)   , INTENT(OUT) :: rexp
INTEGER              , INTENT(IN)  :: npar
!
INTEGER :: iix, iiy, iiz
REAL(kind=PREC_DP) :: sumz
REAL(kind=PREC_DP) :: sumn
REAL(kind=PREC_DP) :: wght
!
sumz = 0.0
sumn = 0.0
!write(*,*) ' DATA_TYPE ', ref_type
!write(*,*) ' BOUND DAT ', lbound(ref_data), ubound(ref_data)
!write(*,*) ' BOUND CAL ', lbound(refine_calc), ubound(refine_calc)
!write(*,*) ' X         ', ref_x(1), ref_x(ref_dim(1))
!write(*,*) ' y         ', ref_y(1), ref_y(ref_dim(2))
!write(*,*) ' z         ', ref_z(1), ref_z(ref_dim(3))
!write(*,*) ' VALS        ', minval(ref_data), maxval(ref_data)
!write(*,*) ' CALC        ', minval(refine_calc), maxval(refine_calc)
!write(*,*) ' SIGMAS      ', minval(ref_sigma), maxval(ref_sigma)
!cond_type: if(ref_type==H5_BRAGG_I) then      ! Rvalue for Bragg reflection list
!
!  call refine_rvalue_hkl(sumz, sumn)
!
!else cond_type                     ! All other data sets
!
!
DO iiz = 1, ref_dim(3)
DO iiy = 1, ref_dim(2)
   loop_inner: DO iix = 1, ref_dim(1)
      IF(ref_sigma(iix,iiy,iiz)> 0.0) THEN
         wght = 1./(ref_sigma(iix,iiy,iiz))**2
      elseif(ref_sigma(iix,iiy,iiz)<-1000.0) then
         cycle loop_inner
      else
         wght = 1.0
      ENDIF
      sumz = sumz + wght*(ref_data(iix,iiy, iiz)-refine_calc(iix,iiy, iiz))**2
      sumn = sumn + wght*(ref_data(iix,iiy, iiz)                     )**2
   ENDDO loop_inner
ENDDO
ENDDO
!endif cond_type
IF(sumn/=0.0) THEN
   rval = SQRT(sumz/sumn)
   rexp = SQRT((ref_dim(1)*ref_dim(2)*ref_dim(3)-npar)/sumn)
ELSE
   rval = -1.0
   rexp = -1.0
ENDIF
!
END SUBROUTINE refine_rvalue
!
!*******************************************************************************
!
SUBROUTINE refine_rvalue_tttt(dims, refine_tttt, dvec)
!
USE refine_data_mod
!
use lib_data_types_mod
use precision_mod
!
IMPLICIT NONE
!
integer, dimension(3), intent(in) :: dims
real(kind=prec_DP), dimension(dims(1), dims(2), dims(3),-2:2), intent(in) :: refine_tttt
real(kind=prec_DP), dimension(                          -2:2), intent(in) :: dvec
REAL(kind=PREC_DP)                 :: rval
REAL(kind=PREC_DP)                 :: rexp
!INTEGER              , INTENT(IN)  :: npar
!
integer :: j
INTEGER :: iix, iiy, iiz
REAL(kind=PREC_DP) :: sumz
REAL(kind=PREC_DP) :: sumn
REAL(kind=PREC_DP) :: wght
!
sumz = 0.0
sumn = 0.0
!write(*,*) ' DATA_TYPE ', ref_type
!write(*,*) ' BOUND DAT ', lbound(ref_data), ubound(ref_data)
!write(*,*) ' BOUND CAL ', lbound(refine_calc), ubound(refine_calc)
!write(*,*) ' X         ', ref_x(1), ref_x(ref_dim(1))
!write(*,*) ' y         ', ref_y(1), ref_y(ref_dim(2))
!write(*,*) ' z         ', ref_z(1), ref_z(ref_dim(3))
!write(*,*) ' VALS        ', minval(ref_data), maxval(ref_data)
!write(*,*) ' CALC        ', minval(refine_calc), maxval(refine_calc)
!write(*,*) ' SIGMAS      ', minval(ref_sigma), maxval(ref_sigma)
!cond_type: if(ref_type==H5_BRAGG_I) then      ! Rvalue for Bragg reflection list
!
!  call refine_rvalue_hkl(sumz, sumn)
!
!else cond_type                     ! All other data sets
!
!
do j=-2, 2
DO iiz = 1, ref_dim(3)
DO iiy = 1, ref_dim(2)
   loop_inner: DO iix = 1, ref_dim(1)
      IF(ref_sigma(iix,iiy,iiz)> 0.0) THEN
         wght = 1./(ref_sigma(iix,iiy,iiz))**2
      elseif(ref_sigma(iix,iiy,iiz)<-1000.0) then
         cycle loop_inner
      else
         wght = 1.0
      ENDIF
if(j==0) then
      sumz = sumz + wght*(ref_data(iix,iiy, iiz)-refine_calc(iix,iiy, iiz))**2
else
      sumz = sumz + wght*(ref_data(iix,iiy, iiz)-refine_tttt(iix,iiy, iiz,j))**2
endif
      sumn = sumn + wght*(ref_data(iix,iiy, iiz)                     )**2
   ENDDO loop_inner
ENDDO
ENDDO
!endif cond_type
IF(sumn/=0.0) THEN
   rval = SQRT(sumz/sumn)
!  rexp = SQRT((ref_dim(1)*ref_dim(2)*ref_dim(3)-npar)/sumn)
ELSE
   rval = -1.0
   rexp = -1.0
ENDIF
!write(*,'(a, i3, 2f15.6)') 'RVALUE AT DERIV ', j, rval, dvec(j)
enddo
!
END SUBROUTINE refine_rvalue_tttt
!
!*******************************************************************************
!
subroutine refine_rvalue_hkl(sumz, sumn)
!-
! Calculate R-value for an HKL data set
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP), intent(inout) :: sumz
real(kind=PREC_DP), intent(inout) :: sumn
!
end subroutine refine_rvalue_hkl
!
!*******************************************************************************
!
subroutine refine_init_pop
!-
! Global parameters for DIFFEV interface
!+
!
use refine_mac_mod
use refine_params_mod
!
use diffev_allocate_appl    ! Allocation of DIFFEV population
use run_mpi_mod             ! DIFFEV data structure
use population              ! DIFFEV populationn
!
use blanks_mod
use errlist_mod
use lib_global_flags_mod
use random_state_mod
use support_mod, only:oeffne
!
implicit none
!
character(len=PREC_STRING) :: string
integer, parameter :: IRD=33   ! READ MAC FILE
integer :: i   ! Dummy loop index
integer :: ios ! I/O error status
integer :: ier_cmd
integer :: exit_msg
character(len=PREC_STRING) :: message
!
MAXBACK = 1
!
call alloc_population(refine_par_n * 4, refine_par_n + refine_fix_n)
call alloc_senddata(refine_par_n + refine_fix_n, 1)
!
call oeffne(IRD, refine_mac, 'old')
if(ier_num/=0) then
   return
endif
read(IRD, '(a)', iostat=ios) string
close(IRD)
if(ios/=0) then
   return
endif
ios = len_trim(string)
call rem_dbl_bl(string,ios)
if(string(1:13)=='branch discus') then
   run_mpi_senddata%prog       = 'kuplot'
elseif(string(1:13)=='branch kuplot') then
   run_mpi_senddata%prog       = 'discus'
endif
!
string = 'mkdir -p DISCUS_SUITE_DERIVATIVES'
call execute_command_line(string, CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg, wait=.true.)
!
run_mpi_senddata%generation = 0
run_mpi_senddata%member     = 0
run_mpi_senddata%children   = 0
run_mpi_senddata%parameters = refine_par_n + refine_fix_n
run_mpi_senddata%nindiv     = 1
run_mpi_senddata%kid        = 1
run_mpi_senddata%indiv      = 1
run_mpi_senddata%ierr       = 0
run_mpi_senddata%ierr_typ   = 0
run_mpi_senddata%ierr_msg_l = len(ier_msg)
run_mpi_senddata%ierr_msg_n = ubound(ier_msg,1)
run_mpi_senddata%direc_l    = 1
run_mpi_senddata%prog_l     = 6   ! == len('refine')
run_mpi_senddata%mac_l      = refine_mac_l
run_mpi_senddata%out_l      = 9
run_mpi_senddata%prog_num   = 1
run_mpi_senddata%s_remote   = 0
run_mpi_senddata%port       = 0
call random_current(run_mpi_senddata%nseeds, run_mpi_senddata%seeds)
run_mpi_senddata%n_rvalue_i = 1
run_mpi_senddata%n_rvalue_o = 1
run_mpi_senddata%global_flags = lib_global_flags   ! Use all global flags in current settings
run_mpi_senddata%global_flags(1) = 1 ! Tell KUPLOT to save
run_mpi_senddata%repeat     = .FALSE.
run_mpi_senddata%prog_start = .FALSE.
run_mpi_senddata%l_rvalue   = .FALSE.
run_mpi_senddata%l_get_state= 1 
run_mpi_senddata%l_first_job= .TRUE.
run_mpi_senddata%spacer2    = .FALSE.
run_mpi_senddata%spacer3    = .FALSE.
run_mpi_senddata%ierr_msg   = ' '
run_mpi_senddata%direc      = '.'
!run_mpi_senddata%prog       = 'discus'
run_mpi_senddata%mac        = refine_mac(1:min(len(run_mpi_senddata%mac), len_trim(refine_mac)))
run_mpi_senddata%out        = '/dev/null'
run_mpi_senddata%trial_values = 0.0_PREC_DP
run_mpi_senddata%rvalue       = 0.0_PREC_DP
run_mpi_senddata%trial_names  = ' '
!
run_mpi_senddata%trial_names (1:refine_par_n)  = refine_params(1:refine_par_n)
pop_name                     (             1:refine_par_n               ) = refine_params(1:refine_par_n)
!
if(refine_fix_n>0) then
run_mpi_senddata%trial_names (refine_par_n+1:refine_par_n+refine_fix_n) = refine_fixed(1:refine_fix_n)
run_mpi_senddata%trial_values(refine_par_n+1:refine_par_n+refine_fix_n) = refine_f    (1:refine_fix_n)
!
pop_name                     (refine_par_n+1:refine_par_n+refine_fix_n  ) = refine_fixed (1:refine_fix_n)
do i=1,MAXPOP
pop_t                        (refine_par_n+1:refine_par_n+refine_fix_n,i) = refine_f     (1:refine_fix_n)
enddo
endif
pop_dimx = refine_par_n+refine_fix_n
!
end subroutine refine_init_pop
!
!*******************************************************************************
!
SUBROUTINE refine_best(rval)
!
use refine_control_mod
USE refine_params_mod
USE refine_data_mod
use refine_head_mod
!
use errlist_mod
USE ber_params_mod
use blanks_mod
use build_name_mod
USE get_params_mod
USE precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN) :: rval   ! Current R-value
!
INTEGER, PARAMETER :: IWR=11
INTEGER, PARAMETER :: INP=12
!
CHARACTER(LEN=15), PARAMETER :: ofile='refine_best.mac'
CHARACTER(LEN=15), PARAMETER :: nfile='refine_new.res '
INTEGER             :: i, j
!
INTEGER :: MAXW
INTEGER :: ianz
REAL(kind=PREC_DP)    :: step
REAL(kind=PREC_DP)    :: rval_w
!
character(len=PREC_STRING) ::  string
character(len=PREC_STRING) ::  long_line
integer                    :: length
!
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: cpara
INTEGER            , DIMENSION(:), ALLOCATABLE :: lpara
REAL(KIND=PREC_DP) , DIMENSION(:), ALLOCATABLE :: werte
!
OPEN(UNIT=IWR, FILE=ofile, STATUS='unknown')
!
WRITE(IWR, '(a)') '#@ HEADER'
WRITE(IWR, '(a)') '#@ NAME         refine_best.mac'
WRITE(IWR, '(a)') '#@ '
WRITE(IWR, '(a)') '#@ KEYWORD      refine, best fit'
WRITE(IWR, '(a)') '#@ '
WRITE(IWR, '(a)') '#@ DESCRIPTION  This macro contains the parameters for the current best'
WRITE(IWR, '(a)') '#@ DESCRIPTION  fit status.'
WRITE(IWR, '(a)') '#@ DESCRIPTION  '
WRITE(IWR, '(a)') '#@'
WRITE(IWR, '(a)') '#@ PARAMETER    $0, 0'
WRITE(IWR, '(a)') '#@'
WRITE(IWR, '(a)') '#@ USAGE        @refine_best.mac'
WRITE(IWR, '(a)') '#@'
WRITE(IWR, '(a)') '#@ END'
WRITE(IWR, '(a)') '#'
!
WRITE(IWR,'(a)') 'variable real, Rvalue'
!
WRITE(IWR,'(a)') 'variable integer, F_DATA'
!WRITE(IWR,'(a)') 'variable real, F_XMIN'
!WRITE(IWR,'(a)') 'variable real, F_XMAX'
!WRITE(IWR,'(a)') 'variable real, F_YMIN'
!WRITE(IWR,'(a)') 'variable real, F_YMAX'
!WRITE(IWR,'(a)') 'variable real, F_XSTP'
!WRITE(IWR,'(a)') 'variable real, F_YSTP'
DO i=1, refine_fix_n            ! Make sure each fixed parameter is defined as a variable
   WRITE(IWR,'(a,a)') 'variable real, ',refine_fixed(i)
ENDDO
DO i=1, refine_par_n            ! Make sure each parameter is defined as a variable
   WRITE(IWR,'(a,a)') 'variable real, ',refine_params(i)
ENDDO
! Test for wrong R-values
!
rval_w = rval
IF(ISNAN(rval)) rval_w = -1.0
IF(ABS(rval) >= HUGE(0.0)) rval_w = -1.0
!
WRITE(IWR,'(a,G20.8E3)') 'Rvalue           = ', rval_w
WRITE(IWR,'(a,I15    )') 'F_DATA           = ', ref_kupl
WRITE(IWR,'(a,G20.8E3)') 'F_XMIN           = ', ref_x(1)
WRITE(IWR,'(a,G20.8E3)') 'F_XMAX           = ', ref_x(ref_dim(1))
WRITE(IWR,'(a,G20.8E3)') 'F_YMIN           = ', ref_y(1)
WRITE(IWR,'(a,G20.8E3)') 'F_YMAX           = ', ref_y(ref_dim(2))
step = (ref_x(ref_dim(1))-ref_x(1))/FLOAT(MAX(1,ref_dim(1)-1))
WRITE(IWR,'(a,G20.8E3)') 'F_XSTP           = ', step
step = (ref_y(ref_dim(2))-ref_y(1))/FLOAT(MAX(1,ref_dim(2)-1))
WRITE(IWR,'(a,G20.8E3)') 'F_YSTP           = ', step
!
! Set fixed parameter values
!
MAXW = refine_fix_n
ALLOCATE(cpara(MAXW))
ALLOCATE(lpara(MAXW))
ALLOCATE(werte(MAXW))
DO i = 1, refine_fix_n                  ! Copy parameter names to prepare calculation
   cpara(i) = refine_fixed(i)
   lpara(i) = LEN_TRIM(refine_fixed(i))
ENDDO
ianz = refine_fix_n
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO i=1, refine_fix_n            ! Make sure each parameter is defined as a variable
   WRITE(IWR,'(a,a,G20.8E3)') refine_fixed(i), ' = ', werte(i)
ENDDO
!
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
MAXW = 20
ALLOCATE(cpara(MAXW))
ALLOCATE(lpara(MAXW))
ALLOCATE(werte(MAXW))
!
WRITE(IWR, '(a)') '#'
!
IF(ref_LOAD/= ' ') THEN           ! Data set was loaded 
   i = index(ref_load(1:LEN_TRIM(ref_load)), ',')
   string= ref_load(i+1:LEN_TRIM(ref_load))
   length = len_trim(string)
   call get_params (string, ianz, cpara, lpara, MAXW, length)
   call do_build_name (ianz, cpara, lpara, werte, MAXW, 1)
   write(string,'(a,a,a)') 'load ', ref_load(1:i), cpara(1)(1:lpara(1))
   do i=2, ianz
      string = string(1:len_trim(string)) // ', ' // cpara(i)(1:lpara(i))
   enddo
   WRITE(IWR, '(a)') 'kuplot'
   WRITE(IWR, '(a)') 'rese'
   WRITE(IWR, '(a)') string(1:len_trim(string))
   WRITE(IWR, '(a)') 'exit'
   WRITE(IWR, '(a)') '#'
ENDIF
!
IF(ref_csigma/= ' ') THEN           ! sigma set was loaded 
   i = index(ref_csigma(1:LEN_TRIM(ref_csigma)), ',')
   string= ref_csigma(i+1:LEN_TRIM(ref_csigma))
   length = len_trim(string)
   call get_params (string, ianz, cpara, lpara, MAXW, length)
   call do_build_name (ianz, cpara, lpara, werte, MAXW, 1)
   write(string,'(a,a,a)') 'load ', ref_csigma(1:i), cpara(1)(1:lpara(1))
   do i=2, ianz
      string = string(1:len_trim(string)) // ', ' // cpara(i)(1:lpara(i))
   enddo
   WRITE(IWR, '(a)') 'kuplot'
   WRITE(IWR, '(a)') string(1:len_trim(string))
   WRITE(IWR, '(a)') 'exit'
   WRITE(IWR, '(a)') '#'
ENDIF
WRITE(IWR, '(a)') 'refine'
WRITE(IWR, '(a)') '#'
call write_header(IWR, refine_mac, refine_plot_mac)          ! Write all accumulated header line
WRITE(IWR, '(a)') '#'
!
! Write refined values
!
DO i=1, refine_par_n            ! Write values for all refined parametersWrite values for all refined parameters
   WRITE(IWR,'(a,a,G20.8E3,a, G20.8E3)') refine_params(i), ' = ', refine_p(i) , ' ! +- ', refine_dp(i)
ENDDO
call write_footer(IWR, refine_mac, refine_plot_mac)          ! Write all accumulated footer line
!
WRITE(IWR, '(a)') '#'
WRITE(IWR, '(a,a)') '@',refine_mac(1:LEN_TRIM(refine_mac))
if(refine_plot_mac /= ' ') then   ! User provided a plot macro
   write(IWR, '(a)') 'branch kuplot'
   write(IWR, '(a,a)') '  @',refine_plot_mac(1:LEN_TRIM(refine_plot_mac))
endif
WRITE(IWR, '(a)') '#'
WRITE(IWR, '(a)') 'exit'
!
CLOSE(UNIT=IWR)
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
!
! Write a 'refine_new.res' macro
!

OPEN(UNIT=IWR, FILE=nfile, STATUS='unknown')
OPEN(UNIT=INP, FILE='newpara_new.mac', status='unknown')
write(IWR, '(a)') 'refine'
write(IWR, '(a)') 'reset'
write(IWR, '(a)') '#'
call write_header(IWR, refine_mac, refine_plot_mac)       ! Write accumulated header lines
write(IWR, '(a)') '#'
if(ref_load_u /= ' ') then
   write(IWR, '(2a)') 'data ', ref_load_u(1:len_trim(ref_load_u))
endif
if(ref_csigma_u /= ' ') then
   write(IWR, '(2a)') 'data ', ref_csigma_u(1:len_trim(ref_csigma_u))
endif
write(IWR, '(a)') '#'
write(IWR, '(a)') '@newpara_new.mac'
!
j = 0
do i=1, refine_par_n            ! Write values for all refined parameters
   j = max(j, len_trim(refine_params(i)))
enddo
!
!  Set free parameters
!
do i=1, refine_par_n            ! Write values for all refined parameters
   string = ' '
   if(refine_range(i,1)>refine_range(i,2)) then    ! No range parameter
      string = ' '
   elseif(refine_range(i,1) > -0.5*HUGE(0.0) .and. refine_range(i,2) < 0.5*HUGE(0.0)) then  ! [low, high]
      write(string,'(a,g15.8e3,a,g15.8e3,a)') ', range:[', refine_range(i,1), ',', refine_range(i,2), ']'
   elseif(refine_range(i,1) < -0.5*HUGE(0.0) .and. refine_range(i,2) < 0.5*HUGE(0.0)) then  ! [   , high]
      write(string,'(a,g15.8e3,a)') ', range:[,', refine_range(i,2), ']'
   elseif(refine_range(i,1) > -0.5*HUGE(0.0) .and. refine_range(i,2) > 0.5*HUGE(0.0)) then  ! [low,     ]
      write(string,'(a,g15.8e3,a)') ', range:[', refine_range(i,1), ',]'
   else
      write(string,'(a,g15.8e3,a,g15.8e3,a)') ', range:[', refine_range(i,1), ',', refine_range(i,2), ']'
   endif
   write(long_line, '(a)') refine_params(i)(1:len_trim(refine_params(i)))
   write(long_line(j+1:PREC_STRING), '(a8,G16.8E3,(a,i1),(a,g15.8e3),a,a)') ', value:',refine_p(i), &
   ', points:', refine_nderiv(i), ', shift:',abs(refine_shift(i)), ', status:free', string(1:len_trim(string))
!  write(long_line, '(2a,G20.8E3,(a,i1),(a,g15.8e3),a,a)') refine_params(i)(1:len_trim(refine_params(i))), ' , value:',refine_p(i), &
!  ' , points:', refine_nderiv(i), ' , shift:',abs(refine_shift(i)), ' , status:free', string(1:len_trim(string))
   length = len_trim(long_line)
!  call rem_bl(long_line, length)
   write(INP, '(2a)') 'newpara ', long_line(1:length)
enddo
write(INP, '(a)') '#'
!
! Set fixed parameter values
!
MAXW = refine_fix_n
ALLOCATE(cpara(MAXW))
ALLOCATE(lpara(MAXW))
ALLOCATE(werte(MAXW))
j = 0
DO i = 1, refine_fix_n                  ! Copy parameter names to prepare calculation
   cpara(i) = refine_fixed(i)
   lpara(i) = LEN_TRIM(refine_fixed(i))
   j = max(j, len_trim(refine_fixed(i)))
ENDDO
ianz = refine_fix_n
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO i=1, refine_fix_n            ! Make sure each parameter is defined as a variable
   string = ' '
   if(refine_range_fix(i,1)>refine_range_fix(i,2)) then    ! No range parameter
      string = ' '
   elseif(refine_range_fix(i,1) > -0.5*HUGE(0.0) .and. refine_range_fix(i,2) < 0.5*HUGE(0.0)) then  ! [low, high]
      write(string,'(a,g15.8e3,a,g15.8e3,a)') ', range:[', refine_range_fix(i,1), ',', refine_range_fix(i,2), ']'
   elseif(refine_range_fix(i,1) < -0.5*HUGE(0.0) .and. refine_range_fix(i,2) < 0.5*HUGE(0.0)) then  ! [   , high]
      write(string,'(a,g15.8e3,a)') ', range:[,', refine_range_fix(i,2), ']'
   elseif(refine_range_fix(i,1) > -0.5*HUGE(0.0) .and. refine_range_fix(i,2) > 0.5*HUGE(0.0)) then  ! [low,     ]
      write(string,'(a,g15.8e3,a)') ', range:[', refine_range_fix(i,1), ',]'
   else
      write(string,'(a,g15.8e3,a,g15.8e3,a)') ', range:[', refine_range_fix(i,1), ',', refine_range_fix(i,2), ']'
   endif
   long_line = ' '
   write(long_line, '(a)') refine_fixed(i)(1:len_trim(refine_fixed(i)))
   WRITE(long_line(j+1:PREC_STRING), '(a8,G16.8E3,(a,i1),(a,g15.8e3),a,a)') ', value:', werte(i),  &
   ', points:', refine_nderiv_fix(i), ', shift:',abs(refine_shift_fix(i)), ', status:fixed', string(1:len_trim(string))
!  WRITE(long_line, '(a9,G20.8E3,(a,i1),(a,g15.8e3),a,a)') refine_fixed(i)(1:len_trim(refine_fixed(i))), ', value:', werte(i),  &
!  ' , points:', refine_nderiv_fix(i), ' , shift:',abs(refine_shift_fix(i)), ' , status:fixed', string(1:len_trim(string))
   length = len_trim(long_line)
!  call rem_bl(long_line, length)
   write(INP, '(2a)') 'newpara ', long_line(1:length)
ENDDO
!
write(INP, '(a)') '#'
close(INP)
write(IWR, '(a)') '#'
write(IWR, '(a, i10)') 'set cycle, ', refine_cycles
if(conv_status) then
  string = 'set conver, status:on'
else
  string = 'set conver, status:off'
endif
length = len_trim(string)
!
if(conv_dchi2_u) then          ! User provided dchi:
   write(string(length+1:),'(a,g20.8e3)') ', dchi:', conv_dchi2
   length = len_trim(string)
endif
!
if(conv_dp_sig_u) then          ! User provided pshift:
   write(string(length+1:),'(a,g20.8e3)') ', pshift:', conv_dp_sig
   length = len_trim(string)
endif
!
if(conv_conf_u) then          ! User provided conf:
   write(string(length+1:),'(a,g20.8e3)') ', conf:', conv_conf
   length = len_trim(string)
endif
!
if(conv_chi2_u) then          ! User provided pshift:
   write(string(length+1:),'(a,g20.8e3)') ', chisq:', conv_chi2
   length = len_trim(string)
endif
!
if(conv_lamb_u) then          ! User provided pshift:
   write(string(length+1:),'(a,g20.8e3)') ', lambda:', conv_lambda
   length = len_trim(string)
endif
!
length = len_trim(string)
write(IWR, '(a)') string(1:length)
!
if(refine_lamda_s_u .or. refine_lamda_u_u .or. refine_lamda_d_u) then  ! 'relax'
   string = 'set relax  '
   length = len_trim(string)
!
   if(refine_lamda_s_u) then      ! start:
      write(string(length+1:),'(a,g20.8e3)') ', start:', refine_lamda_s
      length = len_trim(string)
   endif
!
   if(refine_lamda_u_u) then      ! fail:
      write(string(length+1:),'(a,g20.8e3)') ', fail:', refine_lamda_u
      length = len_trim(string)
   endif
!
   if(refine_lamda_d_u) then      ! success:
      write(string(length+1:),'(a,g20.8e3)') ', success:', refine_lamda_d
      length = len_trim(string)
   endif
!
   length = len_trim(string)
   write(IWR, '(a)') string(1:length)
endif
call write_footer(IWR, refine_mac, refine_plot_mac)          ! Write all accumulated footer lines
WRITE(IWR, '(a)') '#'
WRITE(IWR, '(a,a)') '@',refine_mac(1:LEN_TRIM(refine_mac))
if(refine_plot_mac /= ' ') then   ! User provided a plot macro
   write(IWR, '(a)') 'branch kuplot'
   write(IWR, '(a,a)') '  @',refine_plot_mac(1:LEN_TRIM(refine_plot_mac))
endif
write(IWR, '(a)') '#'
write(IWR, '(a,a)') 'run ', ref_run_u(1:len_trim(ref_run_u))
write(IWR, '(a)') '#'
write(IWR, '(a)') 'exit'

close(IWR)
!
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
!
END SUBROUTINE refine_best
!
!*******************************************************************************
!
END MODULE refine_run_mod
