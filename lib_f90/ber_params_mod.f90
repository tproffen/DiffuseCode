module ber_params_mod
!
private
!
public ber_params
public ber_param
public eva_params
!
contains
!
!*****7***************************************************************  
!
subroutine ber_params (ianz, cpara, lpara, werte, maxpara) 
!-                                                                      
!     Calculated the value of all expressions stored in cpara           
!+                                                                      
use berechne_mod
!use calc_expr_mod
use do_replace_expr_mod
use errlist_mod 
use precision_mod
!
implicit none 
!                                                                       
integer, intent(in) :: ianz
integer, intent(in) :: maxpara 
character(len=*   ), dimension(MAXPARA), intent(in)  :: cpara
integer            , dimension(MAXPARA), intent(in)  :: lpara
real(kind=PREC_DP) , dimension(MAXPARA), intent(out) :: werte
!                                                                       
character(len=max(PREC_STRING,len(cpara))) :: line 
integer             :: ll, i ,j
real(kind=PREC_DP)  :: wert 
!                                                                       
!                                                                       
main: do i = 1, ianz 
   ll   = lpara (i) 
   line = ' ' 
   line = '('//cpara (i) (1:ll) //')' 
   ll   = ll + 2 
   do j=1, ll
      if(iachar(line(j:j))==9) line(j:j) = ' '
   enddo
   call do_replace_expr(line,ll)          ! Initially replace any "EXPR" 
   wert = berechne (line, ll) 
   if (ier_num /= 0) exit main
   werte (i) = wert 
enddo  main
!
end subroutine ber_params                     
!
!*****7***************************************************************  
!
subroutine ber_param  (ipara, cpara, lpara, werte, maxpara) 
!-                                                                      
!     Calculated the value of the expressions stored in cpara ipara
!+                                                                      
use berechne_mod
!use calc_expr_mod
use do_replace_expr_mod
use errlist_mod 
use precision_mod
!
implicit none 
!                                                                       
integer, intent(in) :: ipara
integer, intent(in) :: maxpara 
character(len=*   ), dimension(MAXPARA), intent(in)  :: cpara
integer            , dimension(MAXPARA), intent(in)  :: lpara
real(kind=PREC_DP) , dimension(MAXPARA), intent(out) :: werte
!                                                                       
character(len=max(PREC_STRING,len(CPARA))) :: line 
integer             :: ll, i ,j
real(kind=PREC_DP)  :: wert 
!                                                                       
!                                                                       
i = ipara 
   ll   = lpara (i) 
   line = ' ' 
   line = '('//cpara (i) (1:ll) //')' 
   ll   = ll + 2 
   do j=1, ll
      if(iachar(line(j:j))==9) line(j:j) = ' '
   enddo
   call do_replace_expr(line,ll)          ! Initially replace any "EXPR" 
   wert = berechne (line, ll) 
   if (ier_num /= 0) return
   werte (i) = wert 
!
end subroutine ber_param                     
!
!*****7***************************************************************  
!
subroutine eva_params (ianz, cpara, lpara, werte, names, maxpara, inums, inames)
!-                                                                      
!     Evaluates parameters stored in cpara one by one.
!     Those where the numeric evaluation give no error are 
!     stored in 'werte', the other ones are retained in 'names'
!+                                                                      
use berechne_mod
use do_replace_expr_mod
use errlist_mod 
use lib_errlist_func
use precision_mod
!
implicit none 
!                                                                       
integer, intent(in) :: ianz
integer, intent(in) :: maxpara 
character(len=*   ), dimension(MAXPARA), intent(in)  :: cpara
integer            , dimension(MAXPARA), intent(in)  :: lpara
real(kind=PREC_DP) , dimension(ianz)   , intent(out) :: werte
character(len=*   ), dimension(ianz)   , intent(out) :: names
integer                                , intent(out) :: inums
integer                                , intent(out) :: inames
!
real(kind=PREC_DP), dimension(MAXPARA) :: wwerte   ! Automatic dynamic array
!                                                                       
integer :: i
!
inums =  0  ! no numerical parameters yet
inames = 0  ! no character parameters yet
do i=1, ianz
  call no_error
  call ber_param(i, cpara, lpara, wwerte, MAXPARA)
  if(ier_num == 0) then        ! numerical value
     inums = inums + 1
     werte(inums) = wwerte(i)
  else
     inames = inames + 1
     names(inames) = cpara(i)
  endif
enddo
call no_error
!
end subroutine eva_params
!
!*****7***************************************************************  
!
end module ber_params_mod
