module exp2pdf_run_mod
!-
!  Perform the actual transformation powder pattern to PDF
!
!+
!
contains
!
!*******************************************************************************
!
subroutine exp2pdf_run(string, length)
!-
!  Main transformation routine
!+
!
use exp2pdf_data_mod
!
use errlist_mod
use do_set_mod
use get_params_mod
use take_param_mod
use precision_mod
!
implicit none
!
character(len=*), intent(inout) :: string
integer         , intent(inout) :: length
!
integer, parameter :: MAXW = 1
!
integer :: ianz
character(len=PREC_STRING)                  :: line
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
!real(kind=PREC_DP)        , dimension(MAXW) :: werte
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_INTER   = 1
character(LEN=   4), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'mode'/
data loname /  4    /
opara  =  (/ 'silent' /)   ! Always provide fresh default values
lopara =  (/  6       /)
owerte =  (/  0.0     /)
!
call get_params (string, ianz, cpara, lpara, maxw, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(opara(O_INTER)=='inter') then
   exp_inter = .true.
elseif(opara(O_INTER)=='silent') then
   exp_inter = .false.
   line = 'prompt, off, off, save'
   length = len_trim(line)
   call do_set(line, length)
else
   ier_num = -6
   ier_typ = ER_COMM
   return
endif
!
call exp2pdf_spline           ! Spline experimental data onto equidistant grid
if(ier_num /= 0) return
call exp2pdf_init_aver        ! Calculate average atomic form factors
if(ier_num /= 0) return
call exp2pdf_convert          ! Divide by <f>^2 and convert to PDF
if(ier_num /= 0) return
call exp2pdf_adjust           ! Adjust qmax/qmin limits prepare final F(Q)
if(ier_num /= 0) return
!
call exp2pdf_deallocate
!
line = 'prompt, old'
length = len_trim(line)
call do_set(line, length)
!
end subroutine exp2pdf_run
!
!*******************************************************************************
!
subroutine exp2pdf_spline
!
use exp2pdf_data_mod
!
use errlist_mod
use precision_mod
!
implicit none
!
integer :: d_nequi    ! Data points in data equidistant steps
integer :: b_nequi    ! Data points in background equidistant steps
integer :: id_low
integer :: id_high
integer :: ib_low
integer :: ib_high
!
real(kind=PREC_DP) :: q1   ! Dummy Q values
real(kind=PREC_DP) :: q2
real(kind=PREC_DP) :: q3
real(kind=PREC_DP), dimension(:), allocatable :: d_x   ! Data Q-values
real(kind=PREC_DP), dimension(:), allocatable :: d_y   ! Data Intensities
real(kind=PREC_DP), dimension(:), allocatable :: b_x   ! Back Q-values
real(kind=PREC_DP), dimension(:), allocatable :: b_y   ! Back Intensities
real(kind=PREC_DP), dimension(:), allocatable :: s_x   ! Sigma Data Q-values
real(kind=PREC_DP), dimension(:), allocatable :: s_y   ! Sigma Data Intensities
real(kind=PREC_DP), dimension(:), allocatable :: p_x   ! Sigma Back Q-values
real(kind=PREC_DP), dimension(:), allocatable :: p_y   ! Sigma Back Intensities
!
if(exp_inter) call exp2pdf_plot_init
!
if(allocated(exp_temp_x))  deallocate(exp_temp_x)
if(allocated(exp_temp_y))  deallocate(exp_temp_y)
if(allocated(exp_temp_dy)) deallocate(exp_temp_dy)
!
! 1st task, spline data onto an equidistant grid in Q
!
call exp2pdf_spline_dummy(exp_x(1), exp_x(exp_dim(1)), exp_qstep, exp_dim(1), &
                          exp_x, exp_data(:,1,1) , d_nequi, d_x, d_y)
call exp2pdf_spline_dummy(exp_x(1), exp_x(exp_dim(1)), exp_qstep, exp_dim(1), &
                          exp_x, exp_sigma(:,1,1) , d_nequi, s_x, s_y)
!
! If user provided background use, otherwise take data as is
! The qmin/qmax limits are restricted to the smallest interval
! between maximum Q of qmin in (data, background , user qmin)  
! and     minimum Q of qmax in (data, background , user qmax)  
!
if(allocated(exp_back)) then    ! User provided a background file
!
! Spline Background
!
   call exp2pdf_spline_dummy(exp_xb(1), exp_xb(exp_dim_b(1)), exp_qstep, exp_dim_b(1), &
                             exp_xb, exp_back(:,1,1) , b_nequi, b_x, b_y)
   call exp2pdf_spline_dummy(exp_xb(1), exp_xb(exp_dim_b(1)), exp_qstep, exp_dim_b(1), &
                             exp_xb, exp_sigmab(:,1,1) , b_nequi, p_x, p_y)
!
   q1 = d_x(1)                                          ! Lowest Q in equidistant data
   q2 = (int(exp_qmin_i/exp_qstep) + 1) * exp_qstep     ! User Qmin
   q3 = b_x(1)                                          ! Lowest Q in equidistant background
   exp_qmin = max(q1, q2, q3)                           ! Effective Qmin
   id_low = nint( (exp_qmin-d_x(1))/exp_qstep) + 1      ! First index in equidistant data
   ib_low = nint( (exp_qmin-b_x(1))/exp_qstep) + 1      ! First index in equidistant background
!
   q1 = d_x(d_nequi)                                    ! Highest Q in equidistant data
   q2 = (int(exp_qmax_u/exp_qstep) - 1) * exp_qstep     ! User Qmax
   q3 = b_x(d_nequi)                                    ! Highest Q in equidistant background
   exp_qmax = min(q1, q2, q3)                           ! Effective Qmax
   id_high= nint( (exp_qmax-d_x(1))/exp_qstep) + 1      ! First index in equidistant data
   ib_high= nint( (exp_qmax-b_x(1))/exp_qstep) + 1      ! First index in equidistant background
!
   exp_nstep = id_high -id_low + 1
!
   allocate(exp_temp_x( 1:exp_nstep))
   allocate(exp_temp_y( 1:exp_nstep))
   allocate(exp_temp_dy(1:exp_nstep))
   exp_temp_x = d_x(id_low:id_high)
   exp_temp_y = d_y(id_low:id_high) - b_y(ib_low:ib_high) * exp_bscale
   exp_temp_dy= s_y(id_low:id_high) + p_y(ib_low:ib_high) * exp_bscale
   deallocate(d_x)
   deallocate(d_y)
   deallocate(b_x)
   deallocate(b_y)
   deallocate(s_x)
   deallocate(s_y)
   deallocate(p_x)
   deallocate(p_y)
!
   exp_qmax_f = min(exp_qmax_f, exp_qmax)
!
else
!
! No background 
!
   q1 = d_x(1)                                          ! Lowest Q in equidistant data
   q2 = (int(exp_qmin_i/exp_qstep) + 1) * exp_qstep     ! User Qmin
   exp_qmin = max(q1, q2)                               ! Effective Qmin
   id_low = nint( (exp_qmin-d_x(1))/exp_qstep) + 1      ! First index in equidistant data
!
   q1 = d_x(d_nequi)                                    ! Highest Q in equidistant data
   q2 = (int(exp_qmax_u/exp_qstep) - 1) * exp_qstep     ! User Qmax
   exp_qmax = min(q1, q2)                               ! Effective Qmax
   id_high= nint( (exp_qmax-d_x(1))/exp_qstep) + 1      ! First index in equidistant data
!
   exp_nstep = id_high -id_low + 1
!
   allocate(exp_temp_x( 1:exp_nstep))
   allocate(exp_temp_y( 1:exp_nstep))
   allocate(exp_temp_dy(1:exp_nstep))
   exp_temp_x = d_x(id_low:id_high)
   exp_temp_y = d_y(id_low:id_high)
   exp_temp_dy= s_y(id_low:id_high)
   deallocate(d_x)
   deallocate(d_y)
!
   exp_qmax_f = min(exp_qmax_f, exp_qmax)
!
endif
!
! Replace input data by splined data
!
deallocate(exp_data)
deallocate(exp_x   )
deallocate(exp_sigma)
allocate(exp_data  (exp_nstep, 1, 1))
allocate(exp_sigma (exp_nstep, 1, 1))
allocate(exp_x     (exp_nstep))
exp_x   (:    )  = exp_temp_x
exp_data(:,1,1)  = exp_temp_y
exp_sigma(:,1,1) = exp_temp_dy
!
!
end subroutine exp2pdf_spline
!
!*******************************************************************************
!
subroutine exp2pdf_spline_dummy(x_min, x_max, x_step, npkt, x, y, nequi, xequi, yequi)
!-
!   Spline the data onto an equidistant grid
!+
use spline_mod
use precision_mod
!
implicit none
!
real(kind=PREC_DP)                           , intent(in) :: x_min
real(kind=PREC_DP)                           , intent(in) :: x_max
real(kind=PREC_DP)                           , intent(in) :: x_step
integer                                      , intent(in) :: npkt   ! Original number
real(kind=PREC_DP), dimension(npkt)          , intent(in) :: x      ! Original x
real(kind=PREC_DP), dimension(npkt)          , intent(in) :: y      ! Original y
integer                                      , intent(out) :: nequi  ! Equidistant number
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: xequi  ! Equidistant x
real(kind=PREC_DP), dimension(:), allocatable, intent(out) :: yequi  ! Equidistant y
!
integer :: i
integer :: ier_num
real(kind=PREC_DP), dimension(:), allocatable :: y2a   ! Splined intensities
real(kind=PREC_DP) :: q1
real(kind=PREC_DP) :: q2
real(kind=PREC_DP) :: q
!
q1    = (nint(x_min/x_step) + 1) * x_step
q2    = (nint(x_max/x_step) - 1) * x_step
nequi = nint( (q2-q1)/x_step) + 1
allocate(xequi(nequi))      !Temporary data x
allocate(yequi(nequi))      !Temporary data y
allocate(y2a(1:npkt))
call spline(npkt, x, y, 1.d31, 1.d31, y2a)
!
do i=1, nequi
   q = q1 + x_step*real(i-1, kind=PREC_DP)
   call splint (npkt, x, y, y2a, q, q2, ier_num)
   xequi(i) = q
   yequi(i) = q2
!write(77, '(2(2x,g20.8e3))')  exp_temp_x(i), exp_temp_y(i)
enddo
!
deallocate(y2a)
!
end subroutine exp2pdf_spline_dummy
!
!*******************************************************************************
!
subroutine exp2pdf_init_aver
!
! Calculate average atomic form factors <f>^2
!
use crystal_mod
use chem_mod
use discus_allocate_appl_mod
use diffuse_mod
use exp2pdf_data_mod
use fourier_sup
use output_mod
use pdf_mod
use powder_top_mod , only:powder_run, pow_show, pow_conv_limits
use powder_mod
use powder_write_mod
use qval_mod
use spcgr_apply
use structur , only:do_readfree
!
use discus_show_menu
!
use kuplot_mod
!
USE element_data_mod
use precision_mod
!
implicit none
!
integer, parameter :: MAXW = 7
!
integer                                     :: i       ! Dummy loop index
integer                                     :: ik      ! Dat set number in KUPLOT
integer                                     :: ianz
integer                                     :: outvalue 
logical :: lout = .false.
character(len=PREC_STRING)                  :: string
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
!real(kind=PREC_DP)        , dimension(MAXW) :: werte
real(kind=PREC_DP) :: occ
!
if(.not. exp_comp_current) then     ! We cannot use current structure 
   cpara(1:3) = '1.000'
   lpara(1:3) = 5
   cpara(4:6) = '90.000'
   lpara(4:6) = 6
   cpara(7)   = 'abc'
   lpara(7)   = 3
   ianz = 6
   call do_readfree(ianz, maxw, cpara, lpara, cpara(7), lpara(7))
   CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,      &
        cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
        cr_fmat, cr_cartesian,                                      &
        cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
   CALL get_symmetry_matrices
   occ = sum(exp_atocc(1:exp_natom))
   call alloc_crystal(exp_natom, exp_natom)
   cr_at_lis = ' '
   do i=1, exp_natom
      cr_pos(1,i) = REAL(i, kind=PREC_DP)
      cr_pos(2:3,i) = 0.0D0
      cr_prop(i)    = 1
      cr_mole(i)    = 0
      cr_surf(:,i)  = 0
      cr_magn(:,i)  = 0.0
      cr_iscat(i)   = i
      cr_at_lis(i)(1:2)  = exp_atname(i)(1:2)
      cr_scat_equ (i) = .false.
      cr_dw(i)      = 0.0
      cr_occ(i)     = exp_atocc(i)/occ
!
   enddo
   cr_natoms = exp_natom
   cr_nscat  = exp_natom
   cr_dim = 0.0
   cr_dim(1,1) = cr_pos(1,1)
   cr_dim(1,2) = cr_pos(1, cr_natoms)
   cr_dim(2,1) = 0.0
   cr_dim(2,2) = 1.0
   cr_dim(3,1) = 0.0
   cr_dim(3,2) = 1.0
   rlambda = 0.1
   lambda = ' '
   l_energy = .false.
   chem_period(:) = .FALSE.    ! Turn off periodic boundary
   chem_quick     = .FALSE.    ! Turn off quick search mode
endif
!
if(exp_radiation=='xray') then
   lxray = .true.
   diff_radiation = RAD_XRAY
   diff_table = RAD_WAAS
   rlambda = 0.1
   lambda = ' '
   l_energy = .false.
elseif(exp_radiation=='electron') then
   lxray = .true.
   diff_radiation = RAD_ELEC
   diff_table = RAD_INTER
   renergy  = 200.000
   lambda   = ' '
   l_energy = .true.
elseif(exp_radiation=='neutron') then
   lxray = .false.
   diff_radiation = RAD_NEUT
   diff_table = RAD_INTER
   rlambda = 0.1
   lambda = ' '
   l_energy = .false.
endif
!
pow_axis     = POW_AXIS_Q
ano          = .false.
pdf_clin_a   = 0.0D0
pdf_cquad_a  = 0.0D0
pow_delta    = 0.0D0
pow_ka21     = 0.50
pow_ka21_u   = .FALSE.
pow_period   = 0.00D0
pow_lperiod  = .FALSE.
pow_pref     = .false.
pow_profile  = 0
pow_delta    = 0.0
pow_eta      = 0.5
pow_scale    = 1.0D0
pow_four_type= POW_DEBYE
pow_lp       = POW_LP_NONE
ldbw         = .false.
pow_deltaq   = exp_qstep
pow_qmax     = exp_qmax
pow_qtthmax  = .TRUE.
pow_qmin     = exp_qmin
pow_qtthmin  = .TRUE.
pow_qzero    = 0.0D0
pow_qtthzero = .FALSE.
!
!
call dlink(ano, lambda, rlambda, renergy, l_energy, &
           diff_radiation, diff_table, diff_power)
call pow_conv_limits
if(exp_inter) call pow_show      ! Show powder setting in interactive mode
!
string = ' '
i      = 0
call powder_run(string, i)
!
outfile = 'kuplot.exp2pdf.faver2'
outvalue = val_faver2
cpow_form = 'q  '
out_user_limits      = .false.
CALL powder_out (outvalue, .false.)
!
if(allocated(exp_faver2)) deallocate(exp_faver2)
allocate(exp_faver2(1:exp_nstep))
ik = iz - 1
exp_faver2(1:exp_nstep) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
!
end subroutine exp2pdf_init_aver
!
!*******************************************************************************
!
subroutine exp2pdf_convert
!
use exp2pdf_data_mod
!
use kuplot_mod
use kuplot_fit6
use kuplot_fit6_low_mod
use kuplot_para_mod
use kuplot_extrema_mod
use kuplot_show_mod
!
use do_set_mod
use envir_mod
use errlist_mod
use spline_mod
use precision_mod
use prompt_mod
use set_sub_generic_mod
!
implicit none
!
character(len=PREC_STRING) :: string
integer :: length     ! Dummy index
integer :: ik         ! Data set number in KUPLOT
!
exp_temp_y  = exp_temp_y /exp_faver2            ! Divide by <f>^2
exp_temp_dy = exp_temp_dy/exp_faver2            ! Divide by <f>^2
exp_temp_y  = exp_temp_y *exp_temp_x            ! Multiply by Q 
exp_temp_dy = exp_temp_dy*exp_temp_x            ! Multiply by Q 
!
ik = iz - 1
y( offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) = exp_temp_y   ! Replace <f>^2 in KUPLOT by I/<f>^2
dy(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) = exp_temp_dy  ! Replace <f>^2 in KUPLOT by I/<f>^2
!
! Fit a polynomial through temporary F(Q)
!
ikfit = iz - 1
call fit_init                                        ! Initialize KUPLOT/fit
!
write(string,'(a,i3)') 'poly, ', exp_npoly           ! func : 'poly, 7'
length =  len_trim(string)
call do_fit_fkt(string, length)
!
string = 'para 1, 1, 0.0'
length = len_trim(string)
call do_fit_par(string, length)
!
string = 'one'
length = len_trim(string)
call do_fit_wichtung(string, length)
!
ncycle = 80
!
string = 'status:off'
length = len_trim(string)
call kuplot_set_convergence(string, length)
!
string = ' '
length = 0
call set_skal(string, length)
call get_extrema
!
if(exp_inter) then
   CALL do_fit_info (output_io, .true., .true., .true.)
!else
!   CALL do_fit_info (output_io, .false., .false., .false.)
endif
!
fstart = .false.          ! ! no output
string = 'prompt, off'
length = 11
call do_set(string,length)
string = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE., 1)  ! Initialize into KUPLOT
CALL do_fit
call p_branch_io(string, length, .FALSE.,-1)  ! Return from     KUPLOT
string = 'prompt, on'
length = 10
call do_set(string,length)
!   CALL do_fit_info (output_io, .true., .true., .true.)
!
call get_extrema
!
ik = iz - 1
exp_temp_y(1:exp_nstep) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) !  Get F(Q)-poly
exp_kfq = ik                                          ! Data set F(Q) in KUPLOT
!
!
end subroutine exp2pdf_convert
!
!*******************************************************************************
!
subroutine exp2pdf_adjust
!
! Adjust the lower and upper limits, merge to a common data set
!+
!
use exp2pdf_data_mod
use powder_fft_mod
use discus_output_powder_mod
!
use kuplot_mod
use kuplot_math_mod
use kuplot_extrema_mod
!
use errlist_mod
use envir_mod
use lib_errlist_func
use param_mod
use prompt_mod
use set_sub_generic_mod
!
implicit none
!
integer, parameter :: MAXMAX =  20
!
character(len=PREC_STRING) :: string
integer :: i          ! Dummy index
integer :: ii         ! Dummy index
integer :: j          ! Dummy index
integer :: length     ! Dummy index
integer :: iqmin      ! (rough) index for Qmax
integer :: iqmax      ! (rough) index for Qmax
integer :: iqfirst    !         index for first maximum in F(Q)
integer :: iqerstes   !         index for first maximum in F(Q)
integer :: istep      ! Step direction for "zero" search
integer :: nstep      ! Limit for "zero" search
integer :: ik         ! Data set number in KUPLOT
integer :: npkt_wrt   ! Number of data points in F(Q)
integer   :: npkt_fft    ! number of points in powder pattern for Fast Fourier
logical            :: lsuccess
real(kind=PREC_DP) :: qmax   ! User supplied Qmax
real(kind=PREC_DP) :: qmin   ! qmin into FFT
real(kind=PREC_DP) :: scalex ! Scale factor for Q-axis
real(kind=PREC_DP) :: scalef ! Magic scale factor to approximate absolute scale
real(kind=PREC_DP) :: ydummy ! dummy intensity
real(kind=PREC_DP), dimension(:), allocatable :: x_wrt
real(kind=PREC_DP), dimension(:), allocatable :: y_wrt
real(kind=PREC_DP), dimension(:), allocatable :: xfour
real(kind=PREC_DP), dimension(:), allocatable :: yfour
!
integer :: ikk
real(kind=PREC_DP), dimension(:), allocatable :: wmax
real(kind=PREC_DP), dimension(:), allocatable :: wmin
integer           , dimension(:), allocatable :: ixmax
integer           , dimension(:), allocatable :: ixmin
integer :: ima
!
npkt_fft = 2**18
!
ik = exp_kfq    ! F(Q) in KUPLOT
y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik)) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) ! Copy into previous data set
! Smooth previous dat set
write(string, '(i3, a)') ik-1, ',151'
length = len_trim(string)
call do_glat (string, length, .true.)
!
if(exp_inter) call exp2pdf_plot_fq(1) ! Display F(Q) in interactive mode
!
iqmax = nint((exp_qmax_f-exp_temp_x(1))/exp_qstep) + 1
!
! Get slope in smoothed data
!
if(iqmax>exp_nstep-3) then
   istep = -1
   nstep =  2
else
if((y(offxy(ik-2)+iqmax-1) - y(offxy(ik-2)+iqmax))>=0 .and. exp_temp_y(iqmax)>=0) then ! Previous point higher, above zero
   istep =  1                ! Positive steps
   nstep = exp_nstep-1       ! Upper limit
elseif((y(offxy(ik-2)+iqmax-1) - y(offxy(ik-2)+iqmax))>=0 .and. exp_temp_y(iqmax)< 0.0D0) then   ! Previous point higher, below zero
   istep = -1                ! Negative steps
   nstep =  2                ! Lower limit
elseif((y(offxy(ik-2)+iqmax-1) - y(offxy(ik-2)+iqmax))<=0 .and. exp_temp_y(iqmax)>=0) then       ! Previous point lower, above zero
   istep = -1                ! Negative steps
   nstep =  2                ! Lower limit
elseif((y(offxy(ik-2)+iqmax-1) - y(offxy(ik-2)+iqmax))<=0 .and. exp_temp_y(iqmax)< 0.0D0) then   ! Previous point lower, below zero
   istep =  1                ! Positive steps
   nstep = exp_nstep-1       ! Upper limit
else
   istep = -1                ! Default to negative direction
   nstep =  2
endif
endif
!
!  Get zero point in smoothed data
!
lsuccess = .false.
loop_zero: do i=iqmax, nstep, istep
   if(y(offxy(ik-2)+i)*y(offxy(ik-2)+i+istep)<0.0D0) then
      exp_qmax_f = exp_temp_x(i)
      lsuccess = .true.
      exit loop_zero
   endif
enddo loop_zero
!
if(lsuccess) then
   if(abs(exp_temp_y(i))<=abs(exp_temp_y(i+istep))) then
      iqmax = i
   else
      iqmax = i+istep
   endif
else
   iqmax = nstep
endif
exp_qmax_f = exp_temp_x(iqmax)
!
!  Determine qmin for Fourier
!
!  First search for maxima in broadly smoothed data
!
ifen = 251
allocate(wmax( maxmax))
allocate(ixmax(maxmax))
!
ik = exp_kfq    ! F(Q) in KUPLOT
ikk = ik - 1
y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik)) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) ! Copy into previous data set
iframe = 2
call get_extrema
ex(iwin, iframe, 1) = xmin(ikk)
ex(iwin, iframe, 2) = xmax(ikk)
ey(iwin, iframe, 1) = ymin(ikk)
ey(iwin, iframe, 2) = ymax(ikk)
!
call do_fmax_xy(ikk, wmax, ixmax, maxmax, ima)
call no_error             ! Usually more than the MAXMAX=2 extrema are found ignore this error
iqerstes = 0
loop_qmax: do j=1, ima
  if(exp_temp_y(ixmax(j))>0.0D0 .and. exp_temp_y(ixmax(j))>ymax(ikk)*0.05) then
    iqerstes = ixmax(j)
    exit loop_qmax
  endif
enddo loop_qmax
!
ik = exp_kfq    ! F(Q) in KUPLOT
y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik)) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) ! Copy into previous data set
!
! Smooth previous dat set
write(string, '(i3, a)') ik-1, ', 51'
length = len_trim(string)
call do_glat (string, length, .true.)
!
!
iframe = 1
ikk = ik -1                  ! Search in smoothed data
!ikk = ik -0                  ! Search in smoothed data
y(offxy(ikk-1)+1:offxy(ikk-1)+lenc(ik)) = y(offxy(ikk-1)+1:offxy(ikk-1)+lenc(ik))* (-1.0D0)
!
ifen =  75
allocate(wmin( maxmax))
allocate(ixmin(maxmax))
call do_fmax_xy(ikk, wmin, ixmin, maxmax, ima)
call no_error             ! Usually more than the MAXMAX=2 extrema are found ignore this error
!write(*,*) ' QMIN ', exp_qmin_i, exp_qmin_f, exp_qmin, exp_qmin_il, exp_qmin_fl
if(exp_qmin_fl) then
   iqmin = max(1,nint((exp_qmin_f-exp_temp_x(1))/exp_qstep))
!write(*,*) ' TRUE ', exp_qmin_f, exp_temp_x(1), exp_qstep
else
   ydummy = -1.D9
   iqmin = 1
   loop_qmin: do i=1, ima
      if(x(offxy(ik-2)+ixmin(i))>x(offxy(ik-2)+iqerstes)) exit loop_qmin
      if(y(offxy(ik-2)+ixmin(i))>ydummy .and. y(offxy(ik-2)+ixmin(i))>0.0) then
         iqmin = ixmin(i)
         ydummy = y(offxy(ik-2)+ixmin(i))
!         exit loop_qmin
      endif
   enddo loop_qmin
!write(*,*) ' FALSE ', ixmin(1)
   exp_qmin_f = exp_temp_x(iqmin)
endif
y(offxy(ikk-1)+1:offxy(ikk-1)+lenc(ik)) = y(offxy(ikk-1)+1:offxy(ikk-1)+lenc(ik))*( -1.0D0)
!
if(exp_inter) then
   exp_qmin_f = max(exp_temp_x(1), exp_qmin_i, exp_qmin_f)
   call exp2pdf_plot_qmin        ! Display F(Q) in interactive mode
   iqmin = max(1,nint((exp_qmin_f-exp_temp_x(1))/exp_qstep))
endif
!
!
! Find first maximum above Qmax for electron camera length correction
!
j = 0
if(exp_qfirst_o > 0.0) then    ! Only if user specified a Qfirst
!
! Find first minimum for Qmin extrapolation
   if(exp_inter) call exp2pdf_plot_qscale
   ikk = ik - 1                  ! Search in smoothed data
   call do_fmax_xy(ikk, wmin, ixmin, maxmax, ima)
   call no_error             ! Usually more than the MAXMAX=2 extrema are found ignore this error
   if(exp_qfirst_o< 0.01) then
      iqfirst = nint((exp_qfirst_c-exp_temp_x(1))/exp_qstep) + 1
   else
      iqfirst = nint((exp_qfirst_o-exp_temp_x(1))/exp_qstep) + 1
   endif
   ii      = exp_nstep + 1
   loop_first:do i=1, ima
      if(abs(ixmin(i)-iqfirst) <= ii) then
         j       = ixmin(i)
         ii = abs(ixmin(i)-iqfirst)
      endif
   enddo loop_first
   iqfirst = j
   if(iqfirst> 0) then    ! A maximum was found, correct Q-scale
      scalex = exp_qfirst_c/exp_temp_x(iqfirst)
      exp_qscale = exp_qfirst_c/exp_temp_x(iqfirst)
      exp_qfirst_o  = exp_temp_x(iqfirst)
      x(offxy(ik-1)+1:offxy(ik-1)+lenc(ik  )) = x(offxy(ik-1)+1:offxy(ik-1)+lenc(ik  ))*scalex
      x(offxy(ik-2)+1:offxy(ik-2)+lenc(ik-1)) = x(offxy(ik-2)+1:offxy(ik-2)+lenc(ik-1))*scalex
!     exp_qfirst_o  = exp_qfirst_o * scalex
      exp_temp_x    = exp_temp_x * scalex
      exp_x         = exp_x      * scalex
      exp_qmax_f    = exp_qmax_f * scalex
   endif
endif

deallocate(wmin)
deallocate(ixmin )
deallocate(wmax)
deallocate(ixmax )
!
!
if(allocated(x_wrt)) deallocate(x_wrt)
if(allocated(y_wrt)) deallocate(y_wrt)
!
exp_nstep_f = int(exp_qmax_f/exp_qstep)
allocate(x_wrt(0: iqmax-iqmin+1))
allocate(y_wrt(0: iqmax-iqmin+1))
npkt_wrt = iqmax-iqmin+1
!
!
x_wrt(0:iqmax-iqmin+1) = exp_temp_x(iqmin-1:iqmax)
y_wrt(0:iqmax-iqmin+1) = exp_temp_y(iqmin-1:iqmax)
qmin     = x_wrt(1)
exp_qmin = x_wrt(1)
qmax     = exp_qmax_f
!
! The integral over abs(F(Q)) divided by Qmax and SQRT (maximum value F(Q)) should be about 0.4
! Scale accordingly
scalef = (0.40/(sum(abs(y_wrt))*exp_qstep/sqrt(maxval(y_wrt))/exp_qmax) ) **2
y_wrt = y_wrt * scalef
!
if(exp_inter) call exp2pdf_plot_final ! Display F(Q) in interactive mode
!
!
!open(77, file='TEMP/temp_final.fq', status='unknown')
!do i=0, npkt_wrt
!  write(77, '(2(2x, g30.8e3))') x_wrt(i), y_wrt(i)
!enddo
!close(77)
exp_npdf = int( (exp_rmax-exp_rmin)/exp_rstep) + 1
allocate(xfour(0: exp_npdf     ))
allocate(yfour(0: exp_npdf     ))
!
!write(*,*) exp_qmin, exp_qfirst_o, exp_qmax_f, exp_qmax
!
call fft_fq(npkt_wrt, x_wrt, y_wrt, qmin, qmax, exp_qstep, exp_rmin, exp_rmax, exp_rstep, &
                  npkt_fft, exp_npdf, xfour, yfour)
!
!do i=0, exp_npdf
!  write(77, '(2(2x, g12.8e3))') xfour(i), yfour(i)
!enddo
!close(77)
!
if(exp_outgr_l) then
   call powder_do_write(exp_outgr  , exp_npdf-1, xfour, yfour)
endif
!
if(exp_outfq_l) then
   call powder_do_write(exp_outfq  , npkt_wrt-1, x_wrt, y_wrt)
endif
!
if(exp_outsq_l) then
   do i=1,npkt_wrt
      if(x_wrt(i)>0.0D0) y_wrt(i) = y_wrt(i)/ x_wrt(i)
   enddo
   y_wrt = y_wrt + 1.0D0
   call powder_do_write(exp_outsq  , npkt_wrt-1, x_wrt, y_wrt)
endif
!
if(exp_outiq_l) then
   deallocate(x_wrt)
   deallocate(y_wrt)
   npkt_wrt = ubound(exp_temp_x,1)
   allocate(x_wrt(0:npkt_wrt            ))
   allocate(y_wrt(0:npkt_wrt            ))
   x_wrt(0:npkt_wrt  -1) = exp_x     (1:npkt_wrt  )
   y_wrt(0:npkt_wrt  -1) = exp_data  (1:npkt_wrt,1, 1)
   call powder_do_write(exp_outiq  , npkt_wrt  -1, x_wrt, y_wrt)
endif
!
if(allocated(x_wrt)) deallocate(x_wrt)
if(allocated(x_wrt)) deallocate(x_wrt)
if(allocated(yfour)) deallocate(yfour)
if(allocated(yfour)) deallocate(yfour)
!
end subroutine exp2pdf_adjust
!
!*******************************************************************************
!
subroutine exp2pdf_deallocate
!-
!  After a run deallocate all temporary variables
!+
use exp2pdf_data_mod
!
implicit none
!
if(allocated(exp_data  )) deallocate ( exp_data   ) ! the actual data set
if(allocated(exp_back  )) deallocate ( exp_back   ) ! the actual data set
if(allocated(exp_sigma )) deallocate ( exp_sigma  ) ! sigma at each data point
if(allocated(exp_x     )) deallocate ( exp_x      ) ! x-values of data set
if(allocated(exp_y     )) deallocate ( exp_y      ) ! y-values of data set
if(allocated(exp_xb    )) deallocate ( exp_xb     ) ! x-values of background set
if(allocated(exp_faver2)) deallocate ( exp_faver2 ) ! Calculated <f>**2
if(allocated(exp_temp_x)) deallocate ( exp_temp_x ) ! the actual data set
if(allocated(exp_temp_y)) deallocate ( exp_temp_y ) ! the actual data set
!if(allocated(exp_temp_b)) deallocate ( exp_temp_b ) !
if(allocated(exp_pdf_x )) deallocate ( exp_pdf_x  ) !
if(allocated(exp_pdf_y )) deallocate ( exp_pdf_y  ) !
!
end subroutine exp2pdf_deallocate
!
!*******************************************************************************
!
subroutine  exp2pdf_plot_fq(istep)
!-
!  In interactive mode plot F(Q) al data
!+
use exp2pdf_data_mod
!
use kuplot_mod
use kuplot_extrema_mod
use kuplot_frame_mod
use kuplot_para_mod
use kuplot_plot_mod
!
use errlist_mod
use envir_mod
use prompt_mod
use precision_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
integer, intent(in) :: istep
!
character(len=PREC_STRING) :: string
integer :: length
!integer :: i
integer :: ios
integer, save :: ikk = 0
real(kind=PREC_DP) :: mstep
real(kind=PREC_DP) :: xx
real(kind=PREC_DP), save :: ytic
!
string = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE., 1     )        ! Initialize into KUPLOT 
!
if((exp_temp_x(exp_nstep)-exp_temp_x(1))>40.0D0) then
   mstep = 10.0D0
elseif((exp_temp_x(exp_nstep)-exp_temp_x(1))>20.0D0) then
   mstep = 5.0D0
elseif((exp_temp_x(exp_nstep)-exp_temp_x(1))>10.0) then
   mstep = 2.0D0
elseif((exp_temp_x(exp_nstep)-exp_temp_x(1))>7.0) then
   mstep = 1.0D0
else
   mstep = 0.5D0
endif
if(exp_qfirst_o>0.0) then
  xx = max(exp_temp_x(1), exp_qfirst_o)
else
  xx = max(0.0D0, exp_temp_x(1)) !qmin_f
endif
ikk = exp_kfq
!
!
   iaf(iwin) = 2
   frame(iwin, 1, 1) = 0.00                ! sfra 1, 0.0, 0.45, 1.0, 1.0
   frame(iwin, 1, 2) = 0.45
   frame(iwin, 1, 3) = 1.00
   frame(iwin, 1, 4) = 1.00
   frame(iwin, 2, 1) = 0.00                ! sfra 2, 0.0, 0.00, 1.0, 0.55
   frame(iwin, 2, 2) = 0.00
   frame(iwin, 2, 3) = 1.00
   frame(iwin, 2, 4) = 0.55
!  Allocate vertical line for qscale
   fname(iz) = 'qscale'
   fform(iz) = 'XY'
   lni(iz)   = .false.
   lh5(iz)   = .false.
   lenc(iz)  = 3
   offxy(iz) = offxy(iz-1) + lenc(iz)
   x(offxy(iz - 1) + 1) = xx
   x(offxy(iz - 1) + 2) = xx
   x(offxy(iz - 1) + 3) = xx
   y(offxy(iz - 1) + 1) = ymin(ikk )
   y(offxy(iz - 1) + 3) = 0.0
   y(offxy(iz - 1) + 3) = ymax(ikk )
   offz(iz)             = offz(iz - 1)
   iz = iz + 1
!  Allocate vertical line for qmin
   fname(iz) = 'qmin_u'
   fform(iz) = 'XY'
   lni(iz)   = .false.
   lh5(iz)   = .false.
   lenc(iz)  = 3
   offxy(iz) = offxy(iz-1) + lenc(iz)
   x(offxy(iz - 1) + 1) = max(exp_temp_x(1), exp_qmin_f)
   x(offxy(iz - 1) + 2) = max(exp_temp_x(1), exp_qmin_f)
   x(offxy(iz - 1) + 3) = max(exp_temp_x(1), exp_qmin_f)
   y(offxy(iz - 1) + 1) = ymin(ikk )
   y(offxy(iz - 1) + 3) = 0.0
   y(offxy(iz - 1) + 3) = ymax(ikk )
   offz(iz)             = offz(iz - 1)
   iz = iz + 1
!  Allocate vertical line for qmax_f
   fname(iz) = 'qmax_f'
   fform(iz) = 'XY'
   lni(iz)   = .false.
   lh5(iz)   = .false.
   lenc(iz)  = 3
   offxy(iz) = offxy(iz-1) + lenc(iz)
   x(offxy(iz - 1) + 1) = min(exp_temp_x(exp_nstep), exp_qmax_f)
   x(offxy(iz - 1) + 2) = min(exp_temp_x(exp_nstep), exp_qmax_f)
   x(offxy(iz - 1) + 3) = min(exp_temp_x(exp_nstep), exp_qmax_f)
   y(offxy(iz - 1) + 1) = ymin(ikk )
   y(offxy(iz - 1) + 3) = 0.0
   y(offxy(iz - 1) + 3) = ymax(ikk )
   offz(iz)             = offz(iz - 1)
   iz = iz + 1
   CALL get_extrema
!do i=1, iz-1
!write(*,*) ' DATA SET ', i, xmin(i), xmax(i), ymin(i), ymax(i)
!enddo
   iaf(iwin)        = 2                ! nfra 2
   infra(iwin,1:2, :) = 0              ! Remove all data sets from frame
!  Set range for frame 2
   infra(iwin, 2, 1) = ikk
   iframe = 2
   write(string,'(a,i4,a,i4,a)')  'xmax[', ikk , ']-4.0, xmax[' , ikk , ']+0.1'
   length = len_trim(string)
   call set_skal(string,length)
   iframe = 1
   infra(iwin,1:2, :) = 0              ! Remove all data sets from frame
   infra(iwin,1:2, 1) = iz - 1         ! Set vertival line, F(Q) and F(Q) smoothed
   infra(iwin,1:2, 2) = iz - 2         ! Set vertival line, F(Q) and F(Q) smoothed
   infra(iwin,1:2, 3) = iz - 3         ! Set vertival line, F(Q) and F(Q) smoothed
   infra(iwin,1:2, 4) = iz - 4         ! Set vertival line, F(Q) and F(Q) smoothed
!
!  AFRA 1
   iframe = 1
!  ex(iwin, 1, 1) = -9999.0            ! skal in frame 1
!  ey(iwin, 1, 1) = -9999.0
write(string,'(2(f7.3,a), g12.3e3,a, g12.3e3)') xmin(ikk), ',', xmax(ikk), ',', ymin(ikk), ',', ymax(ikk)
length = len_trim(string)
   call set_skal(string,length)
   ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
   ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
   t(1,iframe,2) = ytic
   t(iwin, 1, 1)     = mstep
!  t(iwin, 1, 2)     = ytic
   titel(iwin, 1, 1) = 'Intermediate F(Q)'
   titel(iwin, 1, 2) = 'vertical lines mark Q\dmin\u and Q\dmax\u'
   achse(iwin, 1:2, 1) = 'Q [\A\u-1\d]'
   achse(iwin, 1:2, 2) = 'F(Q)'
   ifname(iwin, 1)   = .false.                     ! 'fname on'
   ilinecol(iwin, 1:2, iz -1)  = 1   ! lcol qmax_f, red
   ilinecol(iwin, 1:2, iz -2)  = 1   ! lcol qmin_u, red
   ilinecol(iwin, 1:2, iz -3)  = 2   ! lcol qscale, green
   ilinecol(iwin, 1:2, iz -4)  = 3   ! lcol F(Q)  , blue
   ilinecol(iwin,   2, iz -5)  = 1   ! lcol F(Q)  , blue
   ilinetyp(iwin, 1:2, iz-4:iz-1)  = 1   ! ltyp : , 1
   imarktyp(iwin, 1:2, iz-4:iz-1)  = 0   ! mtyp : , 1
   igrid(iwin, 1)    = .false.
!
!  AFRA 2
   iframe = 2
   titel(iwin, iframe, 1) = ' '
   titel(iwin, iframe, 1) = ' '
   t(iwin, iframe, 1)     = 0.5
   ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
   ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
   t(1,iframe,iframe) = ytic
   igrid(iwin, iframe)    = .true.
   ifname(iwin, iframe)   = .false.                    ! 'fname off'
!
   CALL do_plot(.false.)
!
   string = 'kuplot'
   length = 6
   call p_branch_io(string, length, .FALSE.,-1     )        ! Return    from  KUPLOT 
!
!  if(.not. exp_qmax_fl) then        ! User did NOT provide Qmax for Fourier
!
      write(output_io,'(a)') ' Choose approximate Qmax in lower frame'
      write(output_io,'(a)') ' EXP2PDF will set Qmax to nearest Q at which F(Q) = 0.0'
      write(output_io,'(a, f8.3)') ' Current value  : ', exp_qmax_f
      write(output_io, '(a)', advance='no') ' Type new value or Enter to accept current Qmax = '
      string = ' '
      read(*,'(a)', iostat=ios) string
      if(.not.is_iostat_end(ios)) then
         if(string /= ' ') then
            read(string,*) exp_qmax_f
            exp_qmax_f = min(exp_qmax_f, exp_temp_x(exp_nstep), exp_qmax)
         endif
      endif
!
!  else
!     string = 'return'
!     length = 6
!     call do_input(string,length)
!  endif
!
!
end subroutine exp2pdf_plot_fq
!
!*******************************************************************************
!
subroutine  exp2pdf_plot_final
!-
!  In interactive mode plot F(Q) final data
!+
use exp2pdf_data_mod
!
use kuplot_mod
use kuplot_extrema_mod
use kuplot_frame_mod
use kuplot_para_mod
use kuplot_plot_mod
!
use errlist_mod
use envir_mod
use prompt_mod
use precision_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
!
character(len=PREC_STRING) :: string
integer :: length
integer, save :: ikk = 0
real(kind=PREC_DP) :: mstep
real(kind=PREC_DP) :: xx
real(kind=PREC_DP), save :: ytic
!
string = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE., 1     )        ! Initialize into KUPLOT 
!
if((exp_temp_x(exp_nstep)-exp_temp_x(1))>40.0D0) then
   mstep = 10.0D0
elseif((exp_temp_x(exp_nstep)-exp_temp_x(1))>20.0D0) then
   mstep = 5.0D0
elseif((exp_temp_x(exp_nstep)-exp_temp_x(1))>10.0) then
   mstep = 2.0D0
elseif((exp_temp_x(exp_nstep)-exp_temp_x(1))>7.0) then
   mstep = 1.0D0
else
   mstep = 0.5D0
endif
if(exp_qfirst_o>0.0) then
  xx = max(exp_temp_x(1), exp_qfirst_c)
else
!  xx = max(0.0, exp_temp_x(1)) !qmin_f
  xx = 0.0D0
endif
ikk = exp_kfq !max(exp_kload, exp_kback) + 3
ytic   = 10.0**real(int((log(abs((ymax(ikk)-ymin(ikk))*0.5))/log(10.0))))
!
!
iaf(iwin)        = 3                ! nfra 2
frame(iwin, 1, 1) = 0.00                ! sfra 1, 0.0, 0.45, 1.0, 1.0
frame(iwin, 1, 2) = 0.45
frame(iwin, 1, 3) = 1.00
frame(iwin, 1, 4) = 1.00
!
frame(iwin, 2, 1) = 0.00                ! sfra 2, 0.0, 0.00, 1.0, 0.55
frame(iwin, 2, 2) = 0.00
frame(iwin, 2, 3) = 0.58
frame(iwin, 2, 4) = 0.55
!
frame(iwin, 3, 1) = 0.42                ! sfra 2, 0.0, 0.00, 1.0, 0.55
frame(iwin, 3, 2) = 0.00
frame(iwin, 3, 3) = 1.00
frame(iwin, 3, 4) = 0.55
!
infra(iwin,1:3, :) = 0              ! Remove all data sets from frame
!  Set range for frame 2
infra(iwin, 2, 1) = ikk
iframe = 2
write(string,'(f8.3,a, f8.3)' ) max(0.0D0,exp_temp_x(1) - 0.25D0), ',',max(exp_temp_x(1), exp_qmin_f) + 2.5
length = len_trim(string)
call set_skal(string,length)
!
!  Set range for frame 3
infra(iwin, 3, 1) = ikk
iframe = 3
write(string,'(f8.3,a, f8.3)' ) min(exp_qmax_f-2.0D0, xmax(ikk)-4.0), ',', min(exp_qmax_f+2.0D0, xmax(ikk)+0.1)
!write(string,'(a,i4,a,i4,a)')  'xmax[', ikk , ']-4.0, xmax[' , ikk , ']+0.1'
length = len_trim(string)
call set_skal(string,length)
!
if(ey(1,2,2)/ey(1,3,2) > ey(1,2,1)/ey(1,3,1)) then
   ey(1,2,1) = abs(ey(1,2,1))/ey(1,2,1) * abs(ey(1,2,2)) * abs(ey(1,3,1)/ey(1,3,2))
elseif(ey(1,2,2)/ey(1,3,2) < ey(1,2,1)/ey(1,3,1)) then
   ey(1,2,2) = abs(ey(1,2,2))/ey(1,2,2) * abs(ey(1,2,1)) * abs(ey(1,3,2)/ey(1,3,1))
endif
!
iframe = 1
infra(iwin,1:3, :) = 0              ! Remove all data sets from frame
infra(iwin,1:3, 1) = ikk + 3        ! Set vertival line, F(Q) and F(Q) smoothed
infra(iwin,1:3, 2) = ikk + 2        ! Set vertival line, F(Q) and F(Q) smoothed
infra(iwin,1:3, 3) = ikk + 1        ! Set vertival line, F(Q) and F(Q) smoothed
infra(iwin,1:3, 4) = ikk            ! Set vertival line, F(Q) and F(Q) smoothed
!
x(offxy(ikk +0) + 1) = xx
x(offxy(ikk +0) + 2) = xx
x(offxy(ikk +0) + 3) = xx
x(offxy(ikk +1) + 1) = exp_qmin
x(offxy(ikk +1) + 2) = exp_qmin
x(offxy(ikk +1) + 3) = exp_qmin
x(offxy(ikk +2) + 1) = exp_qmax_f
x(offxy(ikk +2) + 2) = exp_qmax_f
x(offxy(ikk +2) + 3) = exp_qmax_f
call get_extrema
!
!  AFRA 1
iframe = 1
write(string,'(2(f7.3,a), g12.3e3,a, g12.3e3)') xmin(ikk), ',', xmax(ikk), ',', ymin(ikk), ',', ymax(ikk)
length = len_trim(string)
call set_skal(string,length)
ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
t(iwin, iframe, 1)     = mstep
t(iwin, iframe, 2)     = ytic
titel(iwin, iframe, 1) = 'Final F(Q)'
if(exp_qfirst_o>0.0D0) then
   titel(iwin, iframe, 2) = 'vertical red lines mark Q\dmin\u and Q\dmax\u; green line marks Q\dscale_o\u'
else
   titel(iwin, iframe, 2) = 'vertical red lines mark Q\dmin\u and Q\dmax\u'
endif
achse(iwin, 1:3, 1) = 'Q [\A\u-1\d]'
achse(iwin, 1:3, 2) = 'F(Q)'
ifname(iwin, iframe)   = .false.                     ! 'fname on'
ilinecol(iwin, 1:3, ikk+3)  = 1   ! lcol qmax_f, red
ilinecol(iwin, 1:3, ikk+2)  = 1   ! lcol qmin_u, red
ilinecol(iwin, 1:3, ikk+1)  = 2   ! lcol qscale, green
ilinecol(iwin, 1:3, ikk  )  = 3   ! lcol F(Q)  , blue
ilinetyp(iwin, 1:3, ikk :ikk+3) = 1   ! ltyp : , 1
imarktyp(iwin, 1:3, ikk :ikk+2) = 0   ! mtyp : , 1
igrid(iwin, iframe)    = .false.
ibox (iwin, iframe)    = 3
!
!  AFRA 2
iframe = 2
titel(iwin, iframe, 1) = ' '
titel(iwin, iframe, 1) = ' '
ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
t(iwin, iframe, 1)     = 0.5
t(iwin, iframe, 2)     = ytic
igrid (iwin, iframe)   = .true.
ifname(iwin, iframe)   = .false.                    ! 'fname off'
ibox  (iwin, iframe)   = 3
!
!  AFRA 3
iframe = 3
titel(iwin, iframe, 1) = ' '
titel(iwin, iframe, 1) = ' '
ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
t(iwin, iframe, 1)     = 1.0
t(iwin, iframe, 2)     = ytic
igrid (iwin, iframe)   = .true.
ifname(iwin, iframe)   = .false.                    ! 'fname off'
ibox  (iwin, iframe)   = -3
!
CALL do_plot(.false.)
!
string = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE.,-1     )        ! Return    from  KUPLOT 
!
!
write(output_io, *)
write(output_io,'(2(a,f7.3),a)') ' Final F(Q) with user defined values for Qmin, Qmax [ ', &
   exp_qmin_f,', ', exp_qmax_f, ' ]'
if(exp_qfirst_o>0.0D0) then
   write(output_io,'(2(a,f7.3),a, f7.3)') ' Q-scale adapted with:   Qscale_obs, Qscale_crystal [ ', &
   exp_qfirst_o,', ', exp_qfirst_c, ' ] ', exp_qscale
endif
write(output_io, *)
!
string = 'return'
length = 6
call do_input(string,length)
!
!
end subroutine exp2pdf_plot_final
!
!*******************************************************************************
!
subroutine exp2pdf_plot_init
!-
!  In interactive mode plot initial data
!+
use exp2pdf_data_mod
!
use kuplot_mod
use kuplot_plot_mod
!
use precision_mod
use prompt_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
!
character(len=PREC_STRING) :: string
integer :: length
real(kind=PREC_DP) :: ytic
!
iaf(iwin)         = 1                  ! 'nfra 1'
iframe = 1
frame(iwin, iframe, 1) = 0.0                !Standard region for a single frame
frame(iwin, iframe, 2) = 0.0
frame(iwin, iframe, 3) = 1.0
frame(iwin, iframe, 4) = 1.0
!
if(exp_kback>0) then                   ! With background
   infra(iwin,iframe, :) = 0                ! Remove all data sets from frame
   infra(iwin,iframe, 1) = exp_kload        ! Set exp_kload and exp_back into frame 1
   infra(iwin,iframe, 2) = exp_kback
   ex(iwin, iframe, 1  ) = xmin(exp_kload)  ! Set skal xmin[exp_kload], xmax[exp_kload], 0.0, ymax[exp_kload]*1.025
   ex(iwin, iframe, 2  ) = xmax(exp_kload)
   ey(iwin, iframe, 1  ) = 0.0
   ey(iwin, iframe, 2  ) = ymax(exp_kload)*1.025
   ilinecol(iwin, iframe, exp_kback)  = 1   ! lcol exp_kback, red
   ilinecol(iwin, iframe, exp_kload)  = 3   ! lcol exp_kload, blue
   ilinetyp(iwin, iframe, exp_kback)  = 1   ! lcol exp_kback, red
   ilinetyp(iwin, iframe, exp_kload)  = 1   ! lcol exp_kload, blue
else
   infra(iwin,iframe, :) = 0                ! Remove all data sets from frame
   infra(iwin,iframe, 1) = exp_kload        ! Set exp_kload and exp_back into frame 1
   ex(iwin, iframe, 1  ) = xmin(exp_kload)  ! Set skal xmin[exp_kload], xmax[exp_kload], 0.0, ymax[exp_kload]*1.025
   ex(iwin, iframe, 2  ) = xmax(exp_kload)
   ey(iwin, iframe, 1  ) = 0.0
   ey(iwin, iframe, 2  ) = ymax(exp_kload)*1.025
   ilinecol(iwin, iframe, exp_kload)  = 3   ! lcol exp_kload, blue
   ilinetyp(iwin, iframe, exp_kload)  = 1   ! lcol exp_kload, blue
endif
!
t(iwin, iframe, 1)     = 5.0
t(iwin, iframe, 2)     = 10.0**real(nint(log(abs(ymax(exp_kload)*1.025)/5.0)/log(10.0)))
ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
t(iwin, iframe, 2)     = ytic
titel(iwin, iframe, 1) = 'Initial intensity'
titel(iwin, iframe, 2) = ' '
achse(iwin, iframe, 1) = 'Q [\A\u-1\d]'
achse(iwin, iframe, 2) = 'Intensity'
ifname(iwin, iframe)   = .true.                     ! 'fname on'
string            = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE., 1)      ! initialize into KUPLOT
CALL do_plot(.false.)
call p_branch_io(string, length, .FALSE.,-1)      ! Retrurn from    KUPLOT
!
write(output_io,'(a)') ' '
write(output_io,'(a)') ' Observed intensity '
if(exp_kback>0) write(output_io,'(a)') '  and unscaled background '
!
string = 'return'
length = 6
call do_input(string,length)
!
end subroutine exp2pdf_plot_init
!
!*******************************************************************************
!
subroutine exp2pdf_plot_qmin
!-
!  Interactive plot for Q-scale
!+
use exp2pdf_data_mod
!
use kuplot_mod
use kuplot_plot_mod
use kuplot_para_mod
!
use envir_mod
use precision_mod
use prompt_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
!
character(len=PREC_STRING) :: string
integer :: length
integer :: ios
integer :: ikk
real(kind=PREC_DP) :: ytic
!
ikk = exp_kfq
!ytic   = (ymax(ikk)-ymin(ikk))*0.25
!
iframe          = 2      ! afra 2
infra(iwin,iframe,:) = 0
infra(iwin,iframe,1) = ikk
!infra(iwin,iframe,3) = ikk-1
write(string,'(f8.3,a, f8.3)' ) max(0.0D0,exp_temp_x(1) - 0.25D0), ',',max(exp_temp_x(1), exp_qmin_f) + 2.5
length = len_trim(string)
call set_skal(string,length)
ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
t(1,iframe,iframe) = ytic
!
infra(iwin,iframe,2) = ikk + 2
ilinecol(iwin, iframe, ikk  ) = 3
ilinecol(iwin, iframe, ikk+2) = 1
x(offxy(ikk +1) + 1) = max(exp_temp_x(1), exp_qmin_f)
x(offxy(ikk +1) + 2) = max(exp_temp_x(1), exp_qmin_f)
x(offxy(ikk +1) + 3) = max(exp_temp_x(1), exp_qmin_f)

string = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE., 1     )     ! Initialize into KUPLOT
call do_plot(.false.)
call p_branch_io(string, length, .FALSE.,-1     )     ! Retrurn from    KUPLOT
!
!if(.not. exp_qfirst_l) then        ! User did NOT provide Qobs for Qscale 
!
   write(output_io,'(a)') ' '
   write(output_io,'(a)') ' Choose approximate Qmin for Fourier in lower frame'
   write(output_io,'(a, f8.3)') ' Current value: ', exp_qmin_f
   write(output_io, '(a)', advance='no') ' Type new value or Enter to accept current Qmin = '
   read(*,'(a)', iostat=ios) string
   if(.not. is_iostat_end(ios)) then
      if(string /= ' ') then
         read(string,*) exp_qmin_f
         exp_qmin_f = max(exp_qmin_f, exp_temp_x(1), exp_qmin)
      endif
   endif
!else
!   string = 'return'
!   length = 6
!   call do_input(string,length)
!endif
!
end subroutine exp2pdf_plot_qmin
!
!*******************************************************************************
!
subroutine exp2pdf_plot_qscale
!-
!  Interactive plot for Q-scale
!+
use exp2pdf_data_mod
!
use kuplot_mod
use kuplot_plot_mod
use kuplot_para_mod
!
use envir_mod
use precision_mod
use prompt_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
!
character(len=PREC_STRING) :: string
integer :: length
integer :: ios
integer :: ikk
real(kind=PREC_DP) :: ytic
!
ikk = exp_kfq
!
iframe          = 2      ! afra 2
infra(iwin,iframe,:) = 0
infra(iwin,iframe,1) = ikk
infra(iwin,iframe,2) = ikk + 1
write(string,'(f8.3,a, f8.3)' ) exp_qfirst_o-1.0D0, ',', exp_qfirst_o+1.0D0
length = len_trim(string)
call set_skal(string,length)
ytic =     10.0**real(int(log(abs(ey(1,iframe,2)-ey(1,iframe,1)    ))/log(10.0)))
ytic = ytic * max(1,int(abs(ey(1,iframe,2)-ey(1,iframe,1))/ytic))*0.25
t(iwin, iframe, 1)     = 0.5
t(iwin, iframe, 2)     = ytic
!
x(offxy(ikk +0) + 1) = max(exp_temp_x(1), exp_qfirst_o)
x(offxy(ikk +0) + 2) = max(exp_temp_x(1), exp_qfirst_o)
x(offxy(ikk +0) + 3) = max(exp_temp_x(1), exp_qfirst_o)
ilinecol(iwin, iframe, ikk  ) = 3
ilinecol(iwin, iframe, ikk+1) = 2
string = 'kuplot'
length = 6
call p_branch_io(string, length, .FALSE., 1     )     ! Initialize into KUPLOT
call do_plot(.false.)
call p_branch_io(string, length, .FALSE.,-1     )     ! Retrurn from    KUPLOT
!
if(.not. exp_qfirst_l) then        ! User did NOT provide Qobs for Qscale 
!
   write(output_io,'(a)') ' '
   write(output_io,'(a)') ' Choose approximate Q for Q-scale in lower frame'
   write(output_io,'(a)') ' EXP2PDF will set Q to highest peak nearby'
   write(output_io,'(a, f8.3)') ' Current value: ', exp_qfirst_o
   write(output_io, '(a)', advance='no') ' Type new value or Enter to accept current Qscale = '
   read(*,'(a)', iostat=ios) string
   if(.not. is_iostat_end(ios)) then
      if(string /= ' ') read(string,*) exp_qfirst_o
   endif
else
   string = 'return'
   length = 6
   call do_input(string,length)
endif
!
end subroutine exp2pdf_plot_qscale
!
!*******************************************************************************
!
subroutine exp2pdf_reset
!
use exp2pdf_data_mod
!
implicit none
!
call exp2pdf_deallocate
!
exp_load         = ' '            ! Load string data
exp_cback        = ' '            ! Load string Background
exp_csigma       = ' '            ! Load string Sigma's
exp_radiation    = 'xray'         ! Radiation 'xray', 'neutron', 'electron'
if(allocated(exp_atname)) deallocate(exp_atname)            ! Atom names
if(allocated(exp_atocc )) deallocate(exp_atocc)             ! Occupancies
exp_natom        = 0              ! Number of atom types
exp_kload        = 0              ! Data set within KUPLOT
exp_kback        = 0              ! Data set within KUPLOT
exp_kfq          = 0              ! F(Q) set within KUPLOT
exp_ksigma       = 0              ! Sigma set within KUPLOT
exp_kupl         = 0              ! Data set within KUPLOT that needs to be kept
exp_dim          = 0              ! Dimensions of data set
exp_dim_b        = 0              ! Dimensions of data set
exp_nstep        = 0              ! Data points on equdistant grid
exp_nstep_f      = 0              ! Data points on equdistant grid
exp_npoly        = 7              ! Order of polynomial through F(Q)
exp_comp_current = .true.         ! Composition is that of current structure
exp_inter        = .true.         ! Show intermittent plots 
exp_bscale        = 1.0D0         ! Background scale
exp_qmin         = 0.0            ! Internal      Qmin
exp_qmax         = 0.0            ! Internal      Qmax
exp_qmin_i       = 0.0D0          ! User supplied Qmin
exp_qmin_f       = 0.0D0          ! User supplied Qmin
exp_qmax_f       = 1.0D9          ! User supplied Qmax for Fourier
exp_qmax_u       = 1.0D9          ! User supplied Qmax for adhoc correction
!
exp_qmin_il      = .false.        ! User supplied Qmin for adhoc correction
exp_qmin_fl      = .false.        ! User supplied Qmin for Fourier
exp_qmax_ul      = .false.        ! User supplied Qmax for adhoc correction
exp_qmax_fl      = .false.        ! User supplied Qmax for Fourier
!
exp_qstep        =   0.001D0      ! Internal Q-step usually 0.001
exp_rmin         =   0.01D0       ! PDF Rmin
!
exp_outgr        = 'discus.grobs' ! Write GROBS
exp_outiq        = 'discus.iqobs' ! Radiation 'xray', 'neutron', 'electron'
exp_outfq        = 'discus.fqobs' ! Radiation 'xray', 'neutron', 'electron'
exp_outsq        = 'discus.sqobs' ! Radiation 'xray', 'neutron', 'electron'
exp_outgr_l      = .true.         ! Radiation 'xray', 'neutron', 'electron'
exp_outiq_l      = .false.        ! Write IQ
exp_outfq_l      = .false.        ! Write FQ
exp_outsq_l      = .false.        ! Write SQ
exp_rmax         = 100.01D0       ! PDF Rmax
exp_rstep        =   0.01D0       ! PDF Rstep
exp_qfirst_o     = 0.00000        ! Q value at first maximum
exp_qfirst_c     = 0.00000        ! Q value at first maximum
exp_qscale       = 1.0D0          ! Q-scale factor
exp_npdf         = 0              ! Number data points in PDF
!
end subroutine exp2pdf_reset
!
!*******************************************************************************
!
end module exp2pdf_run_mod
