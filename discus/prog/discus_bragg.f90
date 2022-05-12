module bragg_mod
!
! Calculations for bragg's equation
!
contains
!
!*******************************************************************************
!
subroutine bragg_eq(string, length)
!-
! Calculate Bragg's equation in various forms:
! lambda = 2 d sin(theta)
! d      = lambda / ( 2 * sin(theta) )
! theta  = asin ( lambda / ( 2 * d ) )
! 
! Calculates the missing term from the other two
!
! lambda | rlambda | energy
! d      | dstar   | Q
! theta  | 2theta
!
use errlist_mod
use get_params_mod
use ber_params_mod
use param_mod
use precision_mod
use prompt_mod
use str_comp_mod
use take_param_mod
use trig_degree_mod
use wink_mod
!
implicit none
!
character(len=*), intent(inout) :: string
integer         , intent(inout) :: length
!
integer, parameter :: MAXW = 8
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
integer                                     :: ianz
integer                                     :: i
integer                                     :: radiation
real(KIND=PREC_DP)  :: theta
real(KIND=PREC_DP)  :: ttheta
real(KIND=PREC_DP)  :: rlambda
real(KIND=PREC_DP)  :: energy 
real(KIND=PREC_DP)  :: dval   
real(KIND=PREC_DP)  :: dstar  
real(KIND=PREC_DP)  :: qstar  
!
character(len=20         ), dimension(3) :: units
!
integer, parameter :: NOPTIONAL = 9
integer, parameter :: O_DSTAR   = 1
integer, parameter :: O_DVAL    = 2
integer, parameter :: O_Q       = 3
integer, parameter :: O_ENERGY  = 4
integer, parameter :: O_THETA   = 5
integer, parameter :: O_2THETA  = 6
integer, parameter :: O_LAMBDA  = 7
integer, parameter :: O_RAD     = 8
integer, parameter :: O_OUT     = 9
character(LEN=   9), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 6 ! Number of values to calculate 
!
data oname  / 'dstar ', 'dval  ',  'q     ',  'energy', 'theta', 'twotheta', 'lambda', 'radiation', 'output   '   /
data loname /  5,        4,         1,         6      ,  5     ,  6        ,  8      ,  9         ,  6            /
opara  =  (/ '1.000000', '1.000000', '6.283185',  '17.39700', '23.80030', '47.60050',  '0.709260', 'xray    ', 'none    ' /) ! Always provide fresh default values
lopara =  (/  8,          8,          8        ,   8        ,  8        ,  8        ,   8         , 8        ,  8         /)
owerte =  (/  1.0,        1.0,        6.283185 ,   17.397   ,  23.8003  ,  47.6005  ,   0.709260  , 1.0E0    ,  0.000     /)
!
units(1) = '     [keV]; xray'
units(2) = '     [meV]; neutron'
units(3) = '     [keV]; electron'
!
call get_params (string, ianz, cpara, lpara, maxw, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
if(str_comp(opara(O_RAD), 'xray', 1, lopara(O_RAD) , 4)) then
   radiation = 1
elseif(str_comp(opara(O_RAD), 'neutron', 1, lopara(O_RAD) , 7)) then
   radiation = 2
elseif(str_comp(opara(O_RAD), 'electron', 1, lopara(O_RAD) , 8)) then
   radiation = 3
else
   ier_num = -183
   ier_typ = ER_APPL
   ier_msg(1) = 'Radiation must be ''xray'', ''neutron'' or ''electron'' '
   return
endif
!
if(lpresent(O_THETA) .and. lpresent(O_2THETA)) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Only one of ''theta'' and ''2theta'' may be present'
   return
elseif(lpresent(O_ENERGY) .and. lpresent(O_LAMBDA)) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Only one of ''energy'' and ''lambda'' may be present'
   return
else
   i = 0
   if(lpresent(O_DVAL)) i = i + 1
   if(lpresent(O_DSTAR)) i = i + 1
   if(lpresent(O_Q)) i = i + 1
   if(i>1) then
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'Only one of ''d'', ''dstar'' and ''q'' may be present'
      return
   endif
endif
if(.not. (lpresent(O_ENERGY) .or. lpresent(O_LAMBDA))) then    ! Lambda is missing
   call get_dvalue(NOPTIONAL, lpresent, O_DVAL, O_DSTAR, O_Q, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, dval)
   if(ier_num/=0) return
   call get_theta (NOPTIONAL, lpresent, O_THETA, O_2THETA, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, theta)
   if(ier_num/=0) return
!
   rlambda = 2.0D0*dval*sind(theta)
!
   dstar   = 1.0D0/dval
   qstar   = zpi  /dval
   ttheta  = 2.0d0 * theta
   lpresent(O_LAMBDA) = .TRUE.
   lpresent(O_ENERGY) = .FALSE.
   write(opara(O_LAMBDA), '(g20.8E3)') rlambda
   call get_lambda(NOPTIONAL, lpresent, O_LAMBDA, O_ENERGY, opara, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, radiation, rlambda, energy)
!  energy  = 0.0
   if(ier_num/=0) return
elseif(.not. (lpresent(O_DVAL) .or. lpresent(O_DSTAR) .or. lpresent(O_Q))) then    !dval is missing
   call get_lambda(NOPTIONAL, lpresent, O_LAMBDA, O_ENERGY, opara, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, radiation, rlambda, energy)
   if(ier_num/=0) return
   call get_theta (NOPTIONAL, lpresent, O_THETA, O_2THETA, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, theta)
   if(ier_num/=0) return
!
   dval    = rlambda / ( 2.0d0 * sind(theta) )
!
   dstar   = 1.0D0/dval
   qstar   = zpi  /dval
   ttheta  = 2.0d0 * theta
elseif(.not. (lpresent(O_THETA) .or. lpresent(O_2THETA))) then    ! theta is missing
   call get_lambda(NOPTIONAL, lpresent, O_LAMBDA, O_ENERGY, opara, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, radiation, rlambda, energy)
   if(ier_num/=0) return
   call get_dvalue(NOPTIONAL, lpresent, O_DVAL, O_DSTAR, O_Q, owerte, ier_num,  &
        ier_typ, ubound(ier_msg,1), ier_msg, ER_FORT, ER_COMM, dval)
   if(ier_num/=0) return
!
   if(rlambda*0.5/dval<=1.0) then
      theta   = asind(rlambda*0.5/dval)
   else
      ier_num = -6
      ier_typ = ER_FORT
      ier_msg(1) = 'Theta would be > 90 deg'
      ier_msg(2) = 'Check d-value and wave length'
   endif
!
   dstar   = 1.0D0/dval
   qstar   = zpi  /dval
   ttheta  = 2.0d0 * theta

else
  ier_num = -6
  ier_typ = ER_COMM
endif
res_para(0) = 7
res_para(1) = rlambda
res_para(2) = dval
res_para(3) = theta
res_para(4) = energy 
res_para(5) = dstar
res_para(6) = qstar
res_para(7) = ttheta
!
if(str_comp(opara(O_OUT), 'screen', 1, lopara(O_OUT) , 6)) then
   write(output_io, *)
   write(output_io, '(a)')                  ' Braggs law '
   write(output_io, '(a,f12.8,a, f9.4, a)') ' Lambda : ',rlambda,' [A]   Energy: ', energy, units(radiation)
   write(output_io, '(3(a,f12.8),a      )') ' d-value: ',dval   ,' [A]   d*    :  ', dstar,  ' [A^-1] Q    : ', qstar, ' [A^-1]'
   write(output_io, '(2(a,f12.8),a      )') ' Theta  : ',theta , ' [deg] 2Theta:  ', ttheta, ' [deg]'
   write(output_io, *)
endif

end subroutine bragg_eq
!
!*******************************************************************************
!
subroutine get_dvalue(NOPTIONAL, lpresent, O_DVAL, O_DSTAR, O_Q, owerte,        &
                      ier_num, ier_typ, msg_dim, ier_msg, ER_FORT, ER_COMM, dval)
!-
!  Determine dvalue from either dval, dstar or q
!+
!
use precision_mod
use wink_mod
!
implicit none
!
integer                                 , intent(in)  :: NOPTIONAL
logical           , dimension(NOPTIONAL), intent(in)  :: lpresent
integer                                 , intent(in)  :: O_DVAL
integer                                 , intent(in)  :: O_DSTAR
integer                                 , intent(in)  :: O_Q
real(kind=PREC_DP), dimension(NOPTIONAL), intent(in)  :: owerte
integer                                 , intent(out) :: ier_num
integer                                 , intent(out) :: ier_typ
integer                                 , intent(in)  :: MSG_DIM
character(len=*), dimension(MSG_DIM)    , intent(out) :: ier_msg
integer                                 , intent(in)  :: ER_FORT
integer                                 , intent(in)  :: ER_COMM
real(kind=PREC_DP)                      , intent(out) :: dval
!
   if(lpresent(O_DVAL)) then
      if(owerte(O_DVAL)==0.0) then
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'd-value must be larger than zero'
         return
      endif
      dval = owerte(O_DVAL)
   elseif(lpresent(O_DSTAR)) then
      if(owerte(O_DSTAR)==0) then
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'd*-value must be larger than zero'
         return
      endif
      dval = 1./owerte(O_DSTAR)
   elseif(lpresent(O_Q)) then
      if(owerte(O_Q)==0.0) then
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'Q-value must be larger than zero'
         return
      endif
      dval = zpi/owerte(O_Q)
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'One of ''dval'', ''dstar'' and ''q'' must be present'
      return
   endif
end subroutine get_dvalue
!
!*******************************************************************************
!
subroutine get_theta (NOPTIONAL, lpresent, O_THETA, O_2THETA, owerte, ier_num,  &
        ier_typ, MSG_DIM, ier_msg, ER_FORT, ER_COMM, theta)
!-
!  Determine theta from either theta, 2theta
!+
!
use precision_mod
use wink_mod
!
implicit none
!
integer                                 , intent(in)  :: NOPTIONAL
logical           , dimension(NOPTIONAL), intent(in)  :: lpresent
integer                                 , intent(in)  :: O_THETA
integer                                 , intent(in)  :: O_2THETA
real(kind=PREC_DP), dimension(NOPTIONAL), intent(in)  :: owerte
integer                                 , intent(out) :: ier_num
integer                                 , intent(out) :: ier_typ
integer                                 , intent(in)  :: MSG_DIM
character(len=*), dimension(MSG_DIM)    , intent(out) :: ier_msg
integer                                 , intent(in)  :: ER_FORT
integer                                 , intent(in)  :: ER_COMM
real(kind=PREC_DP)                      , intent(out) :: theta
!
   if(lpresent(O_THETA)) then
      if(owerte(O_THETA)==0.0) then
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'Theta   must be larger than zero'
         return
      endif
      theta = owerte(O_THETA)
   elseif(lpresent(O_2THETA)) then
      if(owerte(O_2THETA)==0.0) then
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = '2Theta  must be larger than zero'
         return
      endif
      theta = 0.50D0 * owerte(O_2THETA)
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'One of ''theta'' and ''twotheta'' must be present'
      return
   endif
!
end subroutine get_theta
!
!*******************************************************************************
!
subroutine get_lambda(NOPTIONAL, lpresent, O_LAMBDA, O_ENERGY, opara, owerte, ier_num,  &
        ier_typ, MSG_DIM, ier_msg, ER_FORT, ER_COMM, radiation, rlambda, energy)
!-
!  Determine lambda from either lambda (Symbol or wave length) , energy
!+
!
use ber_params_mod
use element_data_mod, only:get_wave
use precision_mod
use string_convert_mod
use wink_mod
!
implicit none
!
integer                                 , intent(in)  :: NOPTIONAL
logical           , dimension(NOPTIONAL), intent(in)  :: lpresent
integer                                 , intent(in)  :: O_LAMBDA
integer                                 , intent(in)  :: O_ENERGY
character(len=*)  , dimension(NOPTIONAL), intent(in)  :: opara
real(kind=PREC_DP), dimension(NOPTIONAL), intent(in)  :: owerte
integer                                 , intent(out) :: ier_num
integer                                 , intent(out) :: ier_typ
integer                                 , intent(in)  :: MSG_DIM
character(len=*), dimension(MSG_DIM)    , intent(out) :: ier_msg
integer                                 , intent(in)  :: ER_FORT
integer                                 , intent(in)  :: ER_COMM
integer                                 , intent(in)  :: radiation
real(kind=PREC_DP)                      , intent(out) :: rlambda
real(kind=PREC_DP)                      , intent(out) :: energy
!
integer, parameter :: MAXW = 1
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
integer                                     :: ianz
real(kind=PREC_DP)        , dimension(MAXW) :: werte   ! Calculated values
!
character(len=4) :: symbol
integer :: lsymbol
logical :: l_energy
!
if(lpresent(O_LAMBDA)) then
   symbol  = opara(O_LAMBDA)(1:len_trim(opara(O_LAMBDA)))
   call do_cap(symbol)
   lsymbol = len_trim(symbol)
   ianz    = 1
   cpara(1) = opara(O_LAMBDA)
   lpara(1) = len_trim(opara(O_LAMBDA))
   CALL ber_params(ianz, cpara, lpara, werte, maxw)
   IF (ier_num == 0) then
      rlambda = werte(1)
      symbol = ' '
      l_energy = .false.
      energy  = 0.0
      call get_wave ( symbol , rlambda, energy, l_energy, &
                radiation,ier_num, ier_typ )
   else
      ier_num = 0
      ier_typ = 0
      l_energy = .false.
      energy  = 0.0
      rlambda = 0.0
      call get_wave ( symbol , rlambda, energy, l_energy, &
                radiation,ier_num, ier_typ )
   endif
elseif(lpresent(O_ENERGY)) then
   ianz    = 1
   cpara(1) = opara(O_ENERGY)
   lpara(1) = len_trim(opara(O_ENERGY))
   CALL ber_params(ianz, cpara, lpara, werte, maxw)
   IF (ier_num /= 0) return
   energy = werte(1)
   symbol = ' '
   rlambda = 0.0
   l_energy = .true.
   call get_wave ( symbol , rlambda, energy, l_energy, &
             radiation,ier_num, ier_typ )
else
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'One of ''lambda'' and ''energy'' must be present'
endif
!
end subroutine get_lambda
!
!*******************************************************************************
!
end module bragg_mod
