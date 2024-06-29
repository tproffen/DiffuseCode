module powder_out_partial_mod
!
contains
!
!*******************************************************************************
!
subroutine powder_out_partial(value, ltemp)
!-                                                                      
!     Convert the powder pattern to a partial PDF
!+                                                                      
use chem_aver_mod
use crystal_mod
use diffuse_mod
use exp2pdf_load_mod
use exp2pdf_run_mod
use exp2pdf_supp_mod
use output_mod
use powder_mod
use powder_write_mod
use save_temp_mod
!
use kuplot_mod
use kuplot_reset_mod
!
use param_mod
use precision_mod
use string_convert_mod
!
integer, intent(in) :: value ! Type of output
logical, intent(in) :: ltemp ! Prepare output in case of periodic boundary conditions
!                            ! True only for periodic boundary and PDF output
!
character(len=PREC_STRING) :: string          ! generic string
character(len=PREC_STRING) :: par_outfile     ! Copy for outputfile from main output menu
character(len=PREC_STRING) :: par_tmpfile     ! Temporary backup of crystal structure
character(len=8), dimension(3)           :: crad
integer                    :: length 
integer                    :: val_temp = 1    ! Value for intensity
integer                    :: ik              ! Kuplot data set number
integer                    :: izz             ! Initial Kuplot data set number
integer,dimension(2) :: par_ipartial
logical                    :: par_user_limits
real(kind=PREC_DP), dimension(3) :: par_rmin
data crad/'xray    ', 'neutron ', 'electron'/
!
! RESET KUPLOT
string = ' '
length = 1
call kuplot_do_reset(string, length)
!
par_tmpfile = 'internal.partial.temp'
call save_temp(par_tmpfile)
!
par_outfile = outfile
if(out_user_limits) then
   par_rmin = out_user_values
else
   par_rmin(1) =   0.01_PREC_DP
   par_rmin(2) = 100.0_PREC_DP
   par_rmin(3) =   0.01_PREC_DP
endif
par_user_limits = out_user_limits
!
out_user_values(1) = pow_qmin_u
out_user_values(2) = pow_qmax_u
out_user_values(3) = pow_deltaq_u
out_user_limits    = .true.
out_user_inc(1)    = int(((pow_qmax_u-pow_qmin_u))/pow_deltaq_u) + 1
!
izz = iz
!
outfile = 'kuplot.inte_temp'    ! temporary file for interim "pseudo intensity"
cpow_form = 'q  '
call powder_out(val_temp, ltemp)
!
!  Transform inte to PDF via exp2pdf
!
!  Reset exp2pdf
!
call exp2pdf_reset
!
!  Load data set 
!
ik = iz - 1
write(string,'(a8,i3)') 'kuplot, ', ik
length = 11
call exp2pdf_load(1, string, length)
!
! Set composition
!
call chem_elem(.false.)     ! Silently get composition

if(pow_ipartial(1)/=pow_ipartial(2)) then
   write(string,'(a17,a2,f5.3,a1, a2, f5.3)') 'composition comp:',&
       cr_at_lis(pow_ipartial(1))(1:2), res_para(pow_ipartial(1)+1), &
   ' ',cr_at_lis(pow_ipartial(2))(1:2), res_para(pow_ipartial(2)+1)
   call do_low(string(19:19))
   call do_low(string(27:27))
   if(string(19:19)=='-' .or. string(19:19)=='-') string(19:19) =  ' ' 
   if(string(27:27)=='-' .or. string(27:27)=='-') string(27:27) =  ' ' 
   length = 32
else
   write(string,'(a17,a2,f5.3)') 'composition comp:',&
      cr_at_lis(pow_ipartial(1))(1:2), 1.0
   call do_low(string(19:19))
   if(string(19:19)=='-' .or. string(19:19)=='-') string(19:19) =  ' ' 
   length = 24
endif
call exp2pdf_composition(string,length)
!
! Set polynomial order
!
string = 'order:7'
length = 7
call exp2pdf_poly(string, length)
!
! set radiation
!
string = 'radiation ' // crad(diff_radiation)
length = len_trim(string)
call exp2pdf_radiation(string, length)
!
! Set limits
write(string,'(a,f8.3,a,f8.3,a)') 'inst:[,',pow_qmax_u-5.0*pow_deltaq_u,'], fourier:[',pow_qmin_u+5.0*pow_deltaq_u,',]'
length = 37
call exp2pdf_qmax(string, length)
!
! Set output file name
!
write(string,'(a,a,3(a,f7.3))') 'gr:', par_outfile(1:len_trim(par_outfile)), &
' , rmin:',par_rmin(1), &
' , rmax:',par_rmin(2), &
' , rstep:',par_rmin(3)
length = 65
!write(*,*) ' OUTF >', string(1:length),'< ', len_trim(string)
call exp2pdf_outfile(string, length)
!
! run command silently
!
par_ipartial = pow_ipartial
pow_ipartial = 0
pow_l_partial = .FALSE.
!
!
string = 'mode:silent'
length = 11
call exp2pdf_run(string,length)
call exp2pdf_init
!
pow_l_partial = .TRUE.
pow_ipartial = par_ipartial
!
call restore_temp(par_tmpfile)
!
! Restore KUPLOT data set numbers
iz = ik + 1
iz = izz
if(iz==1) then
   string = ' '
   length = 1
   call kuplot_do_reset(string, length)
endif
!
end subroutine powder_out_partial
!
end module powder_out_partial_mod
