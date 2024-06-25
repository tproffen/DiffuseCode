module lib_timer_mod
!-
!  Generic interface to the timer class lib_timer_class
!+
!
use lib_timer_class
!
type(timer), dimension(:), allocatable :: ttt
integer :: ntimer
!
private
public   lib_timer
!
contains
!
!*******************************************************************************
!
subroutine lib_timer(line)
!-
!  Main timer routine
!+
!
use ber_params_mod
use errlist_mod
use get_params_mod
use param_mod
use prompt_mod
use str_comp_mod
use take_param_mod
!
implicit none
!
character(len=*) :: line
!
integer, parameter                          :: MAXW = 4
integer, parameter                          :: MAXWW = 1
!
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
integer                                     :: ianz
character(len=PREC_STRING), dimension(MAXWW) :: ccpara
integer                   , dimension(MAXWW) :: llpara
real(kind=PREC_DP)        , dimension(MAXWW) :: wwerte
integer                                      :: iianz
integer                                     :: i        ! Dummy index
integer                                     :: lp
integer                                     :: itimer
logical                                     :: lstart
logical                                     :: loutput
logical                                     :: lcpu
!real(kind=PREC_DP)        , dimension(MAXW) :: lpara
real(kind=PREC_DP) :: zeit
type(timer), dimension(:), allocatable :: tmp
!
integer, parameter :: NOPTIONAL = 4
integer, parameter :: O_NUMBER  = 1
integer, parameter :: O_MODE    = 2
integer, parameter :: O_OUTPUT  = 3
integer, parameter :: O_TYPE    = 4
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'number', 'mode', 'output',  'type'     /
data loname /  6      ,  4    ,  6      ,   4         /
opara  =  (/ '1.0000', '0.0000', '0.0000', '0.0000' /)   ! Always provide fresh default values
lopara =  (/  6,        6,        6,        6       /)
owerte =  (/  1.0,      0.0,      0.0    ,  0.0     /)
!
lp = len_trim(line)
call get_params (line, ianz, cpara, lpara, maxw, lp)
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
lstart = .true.
lcpu   = .true.
itimer = 1
loutput = .true.
!
if(str_comp(cpara(1), 'reset', 3, lpara(1), 5)) then
   if(allocated(ttt)) deallocate(ttt)
   allocate(ttt(1))
   call ttt(1)%reset_timer(lcpu)
   ntimer = 1
   itimer = 1
   return
endif
!
if(lpresent(O_MODE)) then       ! mode: must always be present
   lstart = str_comp(opara(O_MODE), 'start', 5, lopara(O_MODE), 5)
   if(.not. allocated(ttt)) then
      allocate(ttt(1))
      ntimer = 1
      itimer = 1
   endif
!
!  TYPE
!
   if(lstart .and. lpresent(O_TYPE)) THEN
      lcpu = str_comp(opara(O_TYPE), 'cpu', 3, lopara(O_TYPE), 3)
   endif
!
!  NUMBER
!
   if(lpresent(O_NUMBER)) then
      if(opara(O_NUMBER)=='last') then
         itimer = ntimer
      else
         ccpara(1) = opara(O_NUMBER)
         llpara(1) = lopara(O_NUMBER)
         iianz = 1
         call ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
         if(ier_num/=0) then
            ier_msg(1) = 'Error calculating timer number'
            return
         endif
         itimer = nint(wwerte(1))
         if(itimer>ntimer) then
             call move_alloc(ttt, tmp)
             allocate(ttt(itimer))
             do i=ntimer+1,itimer
                call ttt(i)%reset_timer(lcpu)
             enddo
             ttt(1:ntimer) = tmp(1:ntimer)
             ntimer = itimer
         endif
      endif
   endif
endif
!
if(lstart) then
   call ttt(itimer)%start_timer(lcpu)
else
   if(ttt(itimer)%get_started()) then
      zeit = ttt(itimer)%stop_timer()
      res_para(1) = zeit
      res_para(0) = 1
   else
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'Timer was not started'
   endif
   if(ier_num==0) then
      loutput=.false.
      if(lpresent(O_OUTPUT) .and. str_comp(opara(O_OUTPUT), 'screen', 5, lopara(O_OUTPUT), 5) .or. &
         .not. lpresent(O_OUTPUT)) loutput = .true.
      if(loutput) &
      write(output_io,'(a, i8,2x, f8.3,a)') ' Elapsed  time@timer: ', itimer, zeit,' s'
   endif
endif
!
end subroutine lib_timer
!
!*******************************************************************************
!
end module lib_timer_mod
