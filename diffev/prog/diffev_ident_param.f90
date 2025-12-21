module ident_param_mod
!
! Procedures to identify a parameter by number or by name
!
private
public ident_param
contains
!
!*******************************************************************************
!
subroutine ident_param(line, POP_DIMX, pop_name, pname, par_number)
!
use ber_params_mod
use errlist_mod
use precision_mod
!
implicit none
!
character(len=*), intent(inout) :: line
integer         , intent(in)    :: POP_DIMX
character(len=*), dimension(POP_DIMX), intent(in   ) :: pop_name
character(len=*), intent(out  ) :: pname
integer         , intent(out  ) :: par_number

integer, parameter :: MAXW = 1
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
character(len=PREC_STRING) :: string
!
integer :: ianz  ! Dummy index
integer :: i     ! Dummy index
!
!  Check if parameter name exists.
!
pname = line
par_number = 0
loop_par: DO i=1,pop_dimx
   IF(pname == pop_name(i)) THEN
      par_number = i               ! Name existed, use this entry
      EXIT loop_par
   ENDIF
ENDDO loop_par
IF(par_number==0) THEN                ! Not found, search for names PARA0000
   loop_def: DO i=1,pop_dimx
      WRITE(string,'(a4,i4.4)') 'PARA',i
      IF(pop_name(i) == string .OR.  &
         pop_name(i) == 'PARA0000') THEN
         par_number = i               ! Name existed, use this entry
         EXIT loop_def
      ENDIF
   ENDDO loop_def
ENDIF
IF(par_number==0) THEN                ! Not found, evaluate  expression
   cpara(1) = line
   lpara(1) = len_trim(line)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, maxw)
   if(ier_num/=0) return
   par_number = nint(werte(1))
   if(par_number<1 .or. par_number > POP_DIMX) then
      ier_num = -38
      ier_typ = ER_APPL
      ier_msg(1) = 'Wrong parameter ist'
      ier_msg(2) = line(1:min(len(ier_msg),len_trim(line)))
      par_number = 0
   endif
endif
!
end subroutine ident_param
!
!
!*******************************************************************************
!
end module ident_param_mod
