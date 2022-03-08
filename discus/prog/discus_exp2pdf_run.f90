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
use spline_mod
use precision_mod
!
implicit none
!
integer :: i          ! Dummy index
integer :: npkt       ! Data points in experimental field
real(kind=PREC_DP) :: q1   ! Dummy Q values
real(kind=PREC_DP) :: q2
real(kind=PREC_DP), dimension(:), allocatable :: y2a   ! Splined intensities
!
! 1st task, spline data onto an equidistant grid in Q
!
q1 = (int(exp_x(1)  /exp_qstep) + 1) * exp_qstep
q2 = (int(exp_qmin_u/exp_qstep) + 1) * exp_qstep
!
exp_qmin = max(q1, q2)
q1 = (int(exp_x(exp_dim(1))/exp_qstep) - 1) * exp_qstep
q2 = (int(exp_qmax_u       /exp_qstep) - 1) * exp_qstep
!
exp_qmax = min(q1, q2)
exp_qmax_f = min(exp_qmax_f, exp_qmax)
!
if(allocated(exp_temp_x)) deallocate(exp_temp_x)
if(allocated(exp_temp_y)) deallocate(exp_temp_y)
if(allocated(exp_temp_b)) deallocate(exp_temp_b)
if(allocated(y2a   )) deallocate(y2a   )
!
npkt = exp_dim(1)
allocate(y2a(1:npkt))

exp_nstep = nint( (exp_qmax-exp_qmin)/exp_qstep) + 1
allocate(exp_temp_x(exp_nstep))
allocate(exp_temp_y(exp_nstep))
allocate(exp_temp_b(exp_nstep))
exp_temp_b = 0.0D0
!
! Spline Data
!
call spline(npkt, exp_x(1:exp_dim(1)), exp_data(1:exp_dim(1), 1, 1), 1.d31, 1.d31, y2a)
!
!open(77,file='temp.equi', status='unknown')
do i=1, exp_nstep
   q1 = exp_qmin + exp_qstep*real(i-1, kind=PREC_DP)
   call splint (npkt, exp_x(1:exp_dim(1)), exp_data(1:exp_dim(1), 1, 1), y2a, q1, q2, ier_num)
   exp_temp_x(i) = q1
   exp_temp_y(i) = q2
!write(77, '(2(2x,g20.8e3))')  exp_temp_x(i), exp_temp_y(i)
enddo
!close (77)
!
! Spline Background
!
if(allocated(exp_back)) then    ! User provided a background file
!
   if(exp_xb(1) > exp_qmin+0.49*exp_qstep .or.   &
      exp_xb(exp_dim_b(1)) < exp_qmax-0.49*exp_qstep ) then
      ier_num = -6
      ier_typ = ER_COMM
      return
   endif
!
   npkt = exp_dim_b(1)
   call spline(npkt, exp_xb(1:exp_dim_b(1)), exp_back(1:exp_dim_b(1), 1, 1), 1.d31, 1.d31, y2a)
!
!open(77,file='temp.equi', status='unknown')
   do i=1, exp_nstep
      q1 = exp_qmin + exp_qstep*real(i-1, kind=PREC_DP)
      call splint (npkt, exp_xb(1:exp_dim_b(1)), exp_back(1:exp_dim_b(1), 1, 1), y2a, q1, q2, ier_num)
!  exp_temp_xb(i) = q1
      exp_temp_b(i) = q2
!write(77, '(2(2x,g20.8e3))'), exp_temp_x(i), exp_temp_y(i)
   enddo
!close (77)
endif
if(exp_inter) call exp2pdf_plot_init                  ! show intermediate plot
!
exp_temp_y = exp_temp_y - exp_bscale*exp_temp_b       ! Subtract background
!
! Replace input data by splined data
!
deallocate(exp_data)
deallocate(exp_x   )
allocate(exp_data  (exp_nstep, 1, 1))
allocate(exp_x     (exp_nstep))
exp_x   (:    ) = exp_temp_x
exp_data(:,1,1) = exp_temp_y
!
deallocate(y2a)
!
end subroutine exp2pdf_spline
!
!*******************************************************************************
!
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
use powder
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
use envir_mod
use errlist_mod
use spline_mod
use precision_mod
use prompt_mod
use set_sub_generic_mod
!
implicit none
!
integer, parameter :: itmp = 77
!
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: tempfile
integer :: length     ! Dummy index
integer :: ik         ! Data set number in KUPLOT
!
exp_temp_y = exp_temp_y/exp_faver2            ! Divide by <f>^2
exp_temp_y = exp_temp_y*exp_temp_x            ! Multiply by Q 
!
ik = iz - 1
y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) = exp_temp_y   ! Replace <f>^2 in KUPLOT by I/<f>^2
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
tempfile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/exp2pdf_fit.mac'
open(unit=itmp, file=tempfile, status='unknown')
write(itmp, '(a)') 'set prompt, off'
write(itmp, '(a)') 'fit n[1] - 2'
write(itmp, '(a)') 'output off'
write(itmp, '(a)') 'run'
write(itmp, '(a)') 'exit'
write(itmp, '(a)') 'exit'
write(itmp, '(a)') 'set prompt, on'
close(unit=itmp)
!
write(string,'(a,a)')  'kuplot -macro ', tempfile(1:len_trim(tempfile))
length = len_trim(string)
call p_branch(string, length, .FALSE.)
call get_extrema
!
ik = iz - 1
exp_temp_y(1:exp_nstep) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik)) !  Get F(Q)-poly
!
! Delete fit macro file
write(string,'(a,a)')  'rm -f ', tempfile(1:len_trim(tempfile))
call execute_command_line(string)
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
integer, parameter :: MAXMAX =  10
!
character(len=PREC_STRING) :: string
integer :: i          ! Dummy index
integer :: ii         ! Dummy index
integer :: j          ! Dummy index
integer :: length     ! Dummy index
integer :: iqmin      ! (rough) index for Qmax
integer :: iqmax      ! (rough) index for Qmax
integer :: iqfirst    !         index for first maximum in F(Q)
integer :: istep      ! Step direction for "zero" search
integer :: nstep      ! Limit for "zero" search
integer :: ik         ! Data set number in KUPLOT
integer :: npkt_wrt   ! Number of data points in F(Q)
integer   :: npkt_fft    ! number of points in powder pattern for Fast Fourier
logical            :: lsuccess
real(kind=PREC_DP) :: qmax   ! User supplied Qmax
real(kind=PREC_DP) :: qmin   ! qmin into FFT
real(kind=PREC_DP) :: scalex ! Scale factr for Q-axis
real(kind=PREC_DP) :: scalef ! Magic scale factor to approximate absolute scale
real(kind=PREC_DP), dimension(:), allocatable :: x_wrt
real(kind=PREC_DP), dimension(:), allocatable :: y_wrt
real(kind=PREC_DP), dimension(:), allocatable :: xfour
real(kind=PREC_DP), dimension(:), allocatable :: yfour
!
integer :: ikk
real(kind=PREC_DP), dimension(:), allocatable :: wmax
integer           , dimension(:), allocatable :: ixm 
integer :: ima
!
npkt_fft = 2**18
!
ik = iz - 1
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
iframe = 1
y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik)) = y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik))* (-1.0D0)
!
ikk = ik -1                  ! Search in smoothed data
ifen = 251
allocate(wmax(maxmax))
allocate(ixm (maxmax))
call do_fmax_xy(ikk, wmax, ixm, maxmax, ima)
call no_error             ! Usually more than the MAXMAX=2 extrema are found ignore this error
iqmin = ixm(1)
!
y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik)) = y(offxy(ik-2)+1:offxy(ik-2)+lenc(ik))*( -1.0D0)
!
! Find first maximum above Qmax for electron camera length correction
!
j = 0
if(exp_qfirst_o > 0.0) then    ! Only if user specified a Qfirst
!
! Find first minimum for Qmin extrapolation
   if(exp_inter) call exp2pdf_plot_qscale
   ikk = ik - 1                  ! Search in smoothed data
   call do_fmax_xy(ikk, wmax, ixm, maxmax, ima)
   call no_error             ! Usually more than the MAXMAX=2 extrema are found ignore this error
   iqfirst = nint((exp_qfirst_c-exp_x(1))/exp_qstep) + 1
   ii      = exp_nstep + 1
   loop_first:do i=1, ima
      if(abs(ixm(i)-iqfirst) <= ii) then
         j       = ixm(1)
         ii = abs(ixm(i)-iqfirst)
      endif
   enddo loop_first
   iqfirst = j
   if(iqfirst> 0) then    ! A maximum was found, correct Q-scale
      scalex = exp_qfirst_c/exp_temp_x(iqfirst)
      exp_qfirst_o  = exp_x(iqfirst)
      x(offxy(ik-1)+1:offxy(ik-1)+lenc(ik  )) = x(offxy(ik-1)+1:offxy(ik-1)+lenc(ik  ))*scalex
      x(offxy(ik-2)+1:offxy(ik-2)+lenc(ik-1)) = x(offxy(ik-2)+1:offxy(ik-2)+lenc(ik-1))*scalex
      exp_qfirst_o  = exp_qfirst_o * scalex
      exp_temp_x    = exp_temp_x * scalex
      exp_x         = exp_x      * scalex
      exp_qmax_f    = exp_qmax_f * scalex
   endif
endif

deallocate(wmax)
deallocate(ixm )
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
if(exp_inter) call exp2pdf_plot_fq(2) ! Display F(Q) in interactive mode
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
if(allocated(exp_temp_b)) deallocate ( exp_temp_b ) !
if(allocated(exp_temp_b)) deallocate ( exp_pdf_x  ) !
if(allocated(exp_temp_b)) deallocate ( exp_pdf_y  ) !
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
use envir_mod
use prompt_mod
use precision_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
integer, intent(in) :: istep
!
integer, parameter :: itmp = 77
character(len=PREC_STRING) :: tempfile
character(len=PREC_STRING) :: string
integer :: length
integer :: ios
integer, save :: ikk = 0
real(kind=PREC_DP) :: mstep
real(kind=PREC_DP) :: xx
real(kind=PREC_DP), save :: ytic
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
  xx = exp_qfirst_o
else
  xx = exp_qmin_u
endif
if(ikk==0) then
   ikk  = iz -1
   ytic   = (ymax(ikk)-ymin(ikk))*0.25
endif
!
if(istep==1) then                   ! F(Q) initial step 
   tempfile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/exp2pdf_qmax.mac'
   open(unit=itmp, file=tempfile, status='unknown')
!
   write(itmp, '(a)') 'set prompt, off'
   write(itmp, '(a)') 'nfra 2'
   write(itmp, '(a)') 'sfra 1, 0.0, 0.45, 1.0, 1.0'
   write(itmp, '(a)') 'sfra 2, 0.0, 0.00, 1.0, 0.55'
   write(itmp, '(a)') 'alloc line, 3'
   write(itmp, '(a)') 'alloc line, 3'
   write(itmp, '(a)') 'alloc line, 3'
   write(itmp, '(a,G12.6e3)') 'x[n[1]-2,1] = ', xx
   write(itmp, '(a,G12.6e3)') 'x[n[1]-2,2] = ', xx
   write(itmp, '(a,G12.6e3)') 'x[n[1]-2,3] = ', xx
   write(itmp, '(a,G12.6e3)') 'y[n[1]-2,1] =    ymin[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'y[n[1]-2,2] = ', 0.0
   write(itmp, '(a,G12.6e3)') 'y[n[1]-2,3] =    ymax[n[1]-3]' !y[n[1]-3,max(1,nint((2.0000000000)-xmin[n[1]-3]/0.001))]'
   write(itmp, '(a,G12.6e3)') 'x[n[1]-1,1] = ', exp_qmin_u
   write(itmp, '(a,G12.6e3)') 'x[n[1]-1,2] = ', exp_qmin_u
   write(itmp, '(a,G12.6e3)') 'x[n[1]-1,3] = ', exp_qmin_u
   write(itmp, '(a,G12.6e3)') 'y[n[1]-1,1] =    ymin[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'y[n[1]-1,2] = ', 0.0
   write(itmp, '(a,G12.6e3)') 'y[n[1]-1,3] =    ymax[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'x[n[1]  ,1] = ', exp_qmax_f
   write(itmp, '(a,G12.6e3)') 'x[n[1]  ,2] = ', exp_qmax_f
   write(itmp, '(a,G12.6e3)') 'x[n[1]  ,3] = ', exp_qmax_f
   write(itmp, '(a,G12.6e3)') 'y[n[1]  ,1] =    ymin[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'y[n[1]  ,2] = ', 0.0
   write(itmp, '(a,G12.6e3)') 'y[n[1]  ,3] =    ymax[n[1]-3]'
   write(itmp, '(a)') 'kfra 1, n[1], n[1]-1, n[1]-2, n[1]-3, n[1]-4'
   write(itmp, '(a)') 'kfra 2, n[1], n[1]-1, n[1]-2, n[1]-3, n[1]-4'
   write(itmp, '(a)') '#####################'
   write(itmp, '(a)') 'afra 1'
   write(itmp, '(a)') 'skal'
   write(itmp, '(a, f4.1, a, G20.2e3)') 'mark ', mstep, ', ',ytic
   write(itmp, '(a)') 'achx Q [\A\u-1\d]'
   write(itmp, '(a)') 'achy F(Q)'
   write(itmp, '(a)') 'fnam off'
   write(itmp, '(a)') 'grid off'
   write(itmp, '(a)') 'tit1 Final F(Q)'
   write(itmp, '(a)') 'tit2 vertical lines mark Q\dmin\u and Q\dmax\u'
   write(itmp, '(a)') 'lcol n[1]-4 , blue'
   write(itmp, '(a)') 'lcol n[1]-3 , blue'
   write(itmp, '(a)') 'lcol n[1]-2 , green'
   write(itmp, '(a)') 'lcol n[1]-1 , red'
   write(itmp, '(a)') 'lcol n[1]   , red'
   write(itmp, '(a)') 'mtyp n[1]-4, 0'
   write(itmp, '(a)') 'mtyp n[1]-3, 0'
   write(itmp, '(a)') 'mtyp n[1]-2, 0'
   write(itmp, '(a)') 'mtyp n[1]-1, 0'
   write(itmp, '(a)') 'mtyp n[1]  , 0'
   write(itmp, '(a)') 'ltyp n[1]-4, 1'
   write(itmp, '(a)') 'ltyp n[1]-3, 1'
   write(itmp, '(a)') 'ltyp n[1]-2, 1'
   write(itmp, '(a)') 'ltyp n[1]-1, 1'
   write(itmp, '(a)') 'ltyp n[1]  , 1'
   write(itmp, '(a)') '#####################'
   write(itmp, '(a)') 'afra 2'
   write(itmp, '(a)') 'kfra 2, n[1]-3 '
   write(itmp, '(a)') 'skal xmax[n[1]-3]-4, xmax[n[1]-3]+0.1'
   write(itmp, '(a)') 'kfra 2, n[1], n[1]-1, n[1]-2, n[1]-3, n[1]-4'
   write(itmp, '(a, G20.2e3)') 'mark 0.5, ', ytic
   write(itmp, '(a)') 'grid on'
   write(itmp, '(a)') 'tit1'
   write(itmp, '(a)') 'tit2'
   write(itmp, '(a)') 'achx Q [\A\u-1\d]'
   write(itmp, '(a)') 'achy'
   write(itmp, '(a)') 'lcol n[1]-4 , red'
   write(itmp, '(a)') 'lcol n[1]-3 , blue'
   write(itmp, '(a)') 'lcol n[1]-1 , red'
   write(itmp, '(a)') 'lcol n[1]   , red'
   write(itmp, '(a)') 'lwid n[1]-4, 0.9'
   write(itmp, '(a)') 'fnam'
   write(itmp, '(a)') 'plot'
   write(itmp, '(a)') '#'
   write(itmp, '(a)') 'exit'
   write(itmp, '(a)') 'set prompt, on'
   close(unit=itmp)
!
   write(string,'(a,a)')  'kuplot -macro ', tempfile(1:len_trim(tempfile))
   length = len_trim(string)
   call p_branch(string, length, .FALSE.)
!
   if(.not. exp_qmax_fl) then        ! User did NOT provide Qmax for Fourier
!
      write(output_io,'(a)') 'Choose approximate Qmax in lower frame'
      write(output_io,'(a)') 'EXP2PDF will set Qmax to nearest Q at which F(Q) = 0.0'
      write(output_io,'(a, f8.3)') 'Current value  : ', exp_qmax_f
      write(output_io, '(a)', advance='no') ' Type new value or Enter to accept current Qmax = '
      string = ' '
      read(*,'(a)', iostat=ios) string
      if(.not.is_iostat_end(ios)) then
         if(string /= ' ') read(string,*) exp_qmax_f
      endif
!
   else
      string = 'return'
      length = 6
      call do_input(string,length)
   endif
!
   write(string,'(a,a)')  'rm -f ', tempfile(1:len_trim(tempfile))
   call execute_command_line(string)
elseif(istep==2) then

   tempfile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/exp2pdf_qmax2.mac'
   open(unit=itmp, file=tempfile, status='unknown')
!
   write(itmp, '(a)') 'set prompt, off'
   write(itmp, '(a)') 'afra 1'
   write(itmp, '(a)') 'skal'
   write(itmp, '(a, f4.1, a, G20.2e3)') 'mark ', mstep, ', ',ytic
   write(itmp, '(a)') 'afra 2'
   write(itmp, '(a)') 'kfra 2, n[1]-3 '
   write(itmp, '(a)') 'skal xmax[n[1]-3]-4, xmax[n[1]-3]'
   write(itmp, '(a)') 'kfra 2, n[1], n[1]-1, n[1]-2, n[1]-3, n[1]-4'
   write(itmp, '(a, G20.2e3)') 'mark 0.5, ', ytic
   write(itmp, '(a,G12.6e3)') 'x[n[1]-2,1] = ', xx
   write(itmp, '(a,G12.6e3)') 'x[n[1]-2,2] = ', xx
   write(itmp, '(a,G12.6e3)') 'x[n[1]-2,3] = ', xx
   write(itmp, '(a,G12.6e3)') 'y[n[1]-2,1] =    ymin[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'y[n[1]-2,2] = ', 0.0
   write(itmp, '(a,G12.6e3)') 'y[n[1]-2,3] =    ymax[n[1]-3]' !y[n[1]-3,max(1,nint((exp_qfirst_o)-xmin[n[1]-3]]/exp_qstep))'
   write(itmp, '(a,G12.6e3)') 'x[n[1]-1,1] = ', exp_qmin
   write(itmp, '(a,G12.6e3)') 'x[n[1]-1,2] = ', exp_qmin
   write(itmp, '(a,G12.6e3)') 'x[n[1]-1,3] = ', exp_qmin
   write(itmp, '(a,G12.6e3)') 'y[n[1]-1,1] =    ymin[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'y[n[1]-1,2] = ', 0.0
   write(itmp, '(a,G12.6e3)') 'y[n[1]-1,3] =    ymax[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'x[n[1]  ,1] = ', exp_qmax_f
   write(itmp, '(a,G12.6e3)') 'x[n[1]  ,2] = ', exp_qmax_f
   write(itmp, '(a,G12.6e3)') 'x[n[1]  ,3] = ', exp_qmax_f
   write(itmp, '(a,G12.6e3)') 'y[n[1]  ,1] =    ymin[n[1]-3]'
   write(itmp, '(a,G12.6e3)') 'y[n[1]  ,2] = ', 0.0
   write(itmp, '(a,G12.6e3)') 'y[n[1]  ,3] =    ymax[n[1]-3]'
   write(itmp, '(a)') 'tit1 Final F(Q)'
   write(itmp, '(a)') 'tit2 vertical lines mark Q\dmin\u and Q\dmax\u'
   write(itmp, '(a)') 'plot'
!  write(itmp, '(a)') 'save ps, fq_final.ps'
   write(itmp, '(a)') 'n[1] = n[1]-3'    ! effectively delete the last three temproary data sets
   write(itmp, '(a)') 'exit'
   write(itmp, '(a)') 'set prompt, on'
   close(unit=itmp)
!
   write(string,'(a,a)')  'kuplot -macro ', tempfile(1:len_trim(tempfile))
   length = len_trim(string)
   call p_branch(string, length, .FALSE.)
!
   string = 'return'
   length = 6
   call do_input(string,length)
!
   write(string,'(a,a)')  'rm -f ', tempfile(1:len_trim(tempfile))
   call execute_command_line(string)
!
   ikk = 0
endif
!
end subroutine exp2pdf_plot_fq
!
!*******************************************************************************
!
subroutine exp2pdf_plot_init
!-
!  In interactive mode plot initial data
!+
use exp2pdf_data_mod
!
use envir_mod
use prompt_mod
use precision_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
!
integer, parameter :: itmp = 77
character(len=PREC_STRING) :: tempfile
character(len=PREC_STRING) :: string
integer :: length
!
tempfile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/exp2pdf_init.mac'
open(unit=itmp, file=tempfile, status='unknown')
write(itmp, '(a)') 'set prompt, off'
write(itmp, '(a)') 'nfra 1'
if(exp_kback>0) then                   ! With background
   write(itmp, '(a)') 'variable integer, exp_kload'
   write(itmp, '(a)') 'variable integer, exp_kback'
   write(itmp, '(a,i3)') 'exp_kload = ',exp_kload
   write(itmp, '(a,i3)') 'exp_kback = ',exp_kback
   write(itmp, '(a)') 'kfra 1,  exp_kload, exp_kback'
   write(itmp, '(a)') 'skal xmin[exp_kload], xmax[exp_kload], 0.0, ymax[exp_kload]*1.025'
   write(itmp, '(a)') 'lcol exp_kback, red '
   write(itmp, '(a)') 'lcol exp_kload, blue '
   write(itmp, '(a)') 'ltyp exp_kload, 1'
   write(itmp, '(a)') 'ltyp exp_kback, 1'
else
!   write(itmp, '(a)') 'kfra 1, n[1]'
   write(itmp, '(a)') 'variable integer, exp_kload'
   write(itmp, '(a,i3)') 'exp_kload = ',exp_kload
   write(itmp, '(a,i3, a, i3)') 'kfra 1, ', exp_kload
   write(itmp, '(a)') 'skal xmin[exp_kload], xmax[exp_kload], 0.0, ymax[exp_kload]*1.025'
   write(itmp, '(a)') 'lcol exp_kload, blue '
   write(itmp, '(a)') 'ltyp exp_kload, 1'
endif
write(itmp, '(a)') 'mark 5.0'
write(itmp, '(a)') 'tit1 Initial intensity'
write(itmp, '(a)') 'tit2                  '
write(itmp, '(a)') 'fnam on               '
write(itmp, '(a)') 'achx Q [\A\u-1\d]'
write(itmp, '(a)') 'achy Intensity'
write(itmp, '(a)') 'plot'
!write(itmp, '(a)') 'save ps, iq_initial.ps'
write(itmp, '(a)') 'exit'
write(itmp, '(a)') 'set prompt, on'
close(unit=itmp)
!
write(string,'(a,a)')  'kuplot -macro ', tempfile(1:len_trim(tempfile))
length = len_trim(string)
call p_branch(string, length, .FALSE.)
!
string = 'return'
length = 6
call do_input(string,length)
!
write(string,'(a,a)')  'rm -f ', tempfile(1:len_trim(tempfile))
call execute_command_line(string)
!
end subroutine exp2pdf_plot_init
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
!
use envir_mod
use precision_mod
use prompt_mod
use set_sub_generic_mod
use do_wait_mod
!
implicit none
!
integer, parameter :: itmp = 77
character(len=PREC_STRING) :: tempfile
character(len=PREC_STRING) :: string
integer :: length
integer :: ios
integer :: ikk
real(kind=PREC_DP) :: ytic
!
ikk  = max(1,iz - 4)
ytic   = (ymax(ikk)-ymin(ikk))*0.25
!
tempfile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/exp2pdf_qscale.mac'
open(unit=itmp, file=tempfile, status='unknown')
write(itmp, '(a)') 'set prompt, off'
write(itmp, '(a)') 'afra 2'
write(itmp, '(a)') 'kfra 2, n[1]-3, n[1]-2'
write(itmp, '(a,G12.6e3,a, G12.6e3)') 'skal ', exp_qfirst_o-1.0D0, ',', exp_qfirst_o+1.0D0
write(itmp, '(a)') 'lcol n[1]-3, blue'
write(itmp, '(a)') 'lcol n[1]-2, green'
write(itmp, '(a)') 'plot '
write(itmp, '(a)') 'set prompt, on'
write(itmp, '(a)') 'exit '
close(unit=itmp)
!
write(string,'(a,a)')  'kuplot -macro ', tempfile(1:len_trim(tempfile))
length = len_trim(string)
call p_branch(string, length, .FALSE.)
!
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
      if(string /= ' ') read(*,*) exp_qfirst_o
   endif
else
   string = 'return'
   length = 6
   call do_input(string,length)
endif
!
!
write(string,'(a,a)')  'rm -f ', tempfile(1:len_trim(tempfile))
call execute_command_line(string)
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
exp_qmin_u       = 0.0D0          ! User supplied Qmin
exp_qmax_u       = 1.0D9          ! User supplied Qmax for adhoc correction
exp_qmax_f       = 1.0D9          ! User supplied Qmax for Fourier
exp_qmax_ul      = .false.        ! User supplied Qmax for adhoc correction
exp_qmax_fl      = .false.        ! User supplied Qmax for Fourier
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
exp_npdf         = 0              ! Number data points in PDF
!
end subroutine exp2pdf_reset
!
!*******************************************************************************
!
end module exp2pdf_run_mod
