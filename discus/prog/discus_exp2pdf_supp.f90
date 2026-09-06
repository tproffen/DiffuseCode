module exp2pdf_supp_mod
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE exp2pdf_composition(line, length)
!
! Set the composition either as current structure or 
! explicitly as list of atoms
!
USE exp2pdf_data_mod
!
use blanks_mod
use charact_mod
USE errlist_mod
use lib_errlist_func, only: no_error
USE ber_params_mod
USE get_params_mod
use element_data_mod, only: guess_element
USE precision_mod
use string_convert_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
character(len=4)           :: atom_full
character(len=4)           :: atom_name
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: wrong_atom
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz, iianz
integer                              :: i, j
logical :: l_number
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_COMP    = 1
character(LEN=   4), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'comp' /
data loname /  4     /
opara  =  (/ 'current' /) ! Always provide fresh default values
lopara =  (/  7        /)
owerte =  (/  0.0      /)
!
string = ' '
string(1:length) = line(1:length)
!
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_COMP)) then
   if(opara(O_COMP)=='current') then
      exp_comp_current = .true.
   else
      exp_comp_current = .false.
      if(opara(O_COMP)(1:1)=='[') opara(O_COMP)(1:1) = ' '
      i = lopara(O_COMP)
      if(opara(O_COMP)(i:i)==']') opara(O_COMP)(i:i) = ' '
      call rem_bl(opara(O_COMP), lopara(O_COMP))
      string = opara(O_COMP)
      length = lopara(O_COMP)
   endif
else                          ! User omitted "comp:"
   length = len_trim(string)
!   call rem_bl(string, length)
   exp_comp_current = .false.
endif
!
if(allocated(exp_atname)) deallocate(exp_atname)
if(allocated(exp_atocc )) deallocate(exp_atocc )
!
i = len_trim(string)
allocate(exp_atname(1:i))
allocate(exp_atocc (1:i))
exp_atname = ' '
exp_atocc  = 0.0D0
!
call get_params(string, ianz, cpara, lpara, MAXW, length)
if(ianz==1) then     ! Just one atom or blanks
   call get_params_blank(string, ianz, cpara, lpara, MAXW, length)
endif
!
exp_natom = 0
iianz     = 1
l_number  = .false.
loop_params: do i=1, ianz      ! Loop over all params
   wrong_atom = cpara(1)
   call rem_bl(cpara(1), lpara(1))
   call do_cap(cpara(1))
!  write(*,*) ' CPARA ', i, cpara(1)(1:10), l_number
   if(l_number) then           ! This might be a number for previous parameter
      call ber_params(iianz, cpara, lpara, werte, maxw)    ! Try to calculate
      if(ier_num==0) then   ! Success to calculate
         exp_atocc(exp_natom) = werte(1)
         l_number  = .false.        ! Next parameter is an atom name
         call del_params(1, ianz, cpara, lpara, MAXW)
         cycle loop_params          ! Work on next parameter
      else
        call no_error
      endif
   endif
!
   atom_name = 'NONE'
   call guess_element(atom_name, atom_full, j, cpara(1))
   if(j==2) atom_name = atom_full    ! This is an ion
!
   cond_is_atom: if(j/=-1) then     ! Found element

         exp_natom = exp_natom+1
         exp_atname(exp_natom) = atom_name
         cpara(1)(1:len_trim(atom_name)) = ' '     ! CLear atom name
         if(cpara(1)/=' ') then
            call ber_params(iianz, cpara, lpara, werte, maxw)
            if(ier_num/=0) then
               ier_msg(1) = 'Error calculating the composition'
               ier_msg(2) = 'Wrong entry: '''//wrong_atom(1:len_trim(wrong_atom))//''''
               return
            endif
            exp_atocc(exp_natom) = werte(1)
            l_number  = .false.        ! Next parameter is an atom name
         else
            exp_atocc(exp_natom) = 1.0D0  ! Default to one
            l_number  = .true.         ! Next param might be number for current atom
         endif
         call del_params(1, ianz, cpara, lpara, MAXW)
         cycle loop_params
   else cond_is_atom
!
      if(l_number) then
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'Error calculating the composition'
         ier_msg(1) = '  or wrong element name'
         ier_msg(2) = 'Wrong entry: '''//wrong_atom(1:len_trim(wrong_atom))//''''
         return
      else
      ier_num = -204
      ier_typ = ER_APPL
      ier_msg(1) = 'Wrong element: '''//wrong_atom(1:len_trim(wrong_atom))//''''
      ier_msg(2) = 'First character of element name must be'
      ier_msg(3) = 'upper case, letter second lower case'
      ier_msg(4) = 'Separate elements by blank or comma'
      return
      endif
   endif cond_is_atom
enddo loop_params
!
!write(*,*) 'FOUND ', exp_natom
!do i=1, exp_natom
!   write(*,*) ' ATOM NAME >', exp_atname(i)(1:len_trim(exp_atname(i))),'<< ', exp_atocc(i)
!enddo
!
!do i=1, exp_natom
!  write(*,*) 'Atom, occupancy ', exp_atname(i), exp_atocc(i)
!enddo
!
end subroutine exp2pdf_composition
!
!*******************************************************************************
!
SUBROUTINE exp2pdf_radiation(line, length)
!-
!  Hande lusers radiation definition
!
!
use exp2pdf_data_mod
!
use errlist_mod
use lib_errlist_func, only: no_error
use element_data_mod, only: get_wave
use get_params_mod
use str_comp_mod
use take_param_mod
!
implicit none
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 3 
!
character(len=max(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
integer                             :: ianz
INTEGER           , DIMENSION(MAXW) :: lpara
!REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
logical :: l_energy   ! User provided energy instead of lambda
integer :: diff_radiation ! = 1, 2, 3
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_LAMBDA  = 1
integer, parameter :: O_ENERGY  = 2
character(len=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(len=PREC_STRING),dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 2 ! Number of values to calculate 
!
data oname  / 'lambda', 'energy'/
data loname /  6      ,  6      /
opara  =  (/ '0.0000', '0.0000' /)   ! Always provide fresh default values
lopara =  (/  6,        6       /)
owerte =  (/  0.1,    123.9842004    /)
!
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num==-1 .and. ier_typ==2) then
   call no_error
endif
!
!
if(str_comp(cpara(1), 'xray',  2, lpara(1), 4) ) then
   exp_radiation = 'xray'
   diff_radiation = 1
elseif(str_comp(cpara(1), 'electron',  2, lpara(1), 8) ) then
   exp_radiation = 'electron'
   diff_radiation = 3
elseif(str_comp(cpara(1), 'neutron',  2, lpara(1), 7) ) then
   exp_radiation = 'neutron'
   diff_radiation = 2
else
   ier_num = -183
   ier_typ = ER_APPL
   ier_msg(1) = 'Must be either ''xray'', ''neutron'' or '
   ier_msg(2) = '               ''electron'' '
   return
endif
!
if(lpresent(O_LAMBDA) .and. .not.lpresent(O_ENERGY)) then
   exp_is_rad = 1
   l_energy = .false.
!  exp_rlambda = owerte(O_LAMBDA)
!  exp_renergy = 12.39842/exp_rlambda
elseif(.not.lpresent(O_LAMBDA) .and. lpresent(O_ENERGY)) then
   exp_is_rad = 2
   l_energy = .true.
!  exp_renergy = owerte(O_ENERGY)
!  exp_rlambda = 12.39842/exp_renergy
elseif(lpresent(O_LAMBDA) .and. lpresent(O_ENERGY)) then
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Only one of ''lambda:''; ''energy'' may be present'
   return
elseif(.not.lpresent(O_LAMBDA) .and. .not.lpresent(O_ENERGY)) then
   exp_is_rad = 1
   l_energy = .false.
   opara(O_LAMBDA) = ' '
else
   exp_is_rad = 0
   l_energy = .false.
endif
!
call get_wave ( opara(O_LAMBDA), exp_rlambda, exp_renergy, l_energy, &
                      diff_radiation,ier_num, ier_typ )
!
end SUBROUTINE exp2pdf_radiation
!
!
!*******************************************************************************
!
SUBROUTINE exp2pdf_outfile(line, length)
!
! Set the outputfile and the PDF range
!
USE exp2pdf_data_mod
!
!use blanks_mod
use build_name_mod
!use charact_mod
USE errlist_mod
!USE ber_params_mod
USE get_params_mod
USE precision_mod
!use string_convert_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
!REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
!
integer, parameter :: NOPTIONAL = 8
integer, parameter :: O_RMIN    = 1
integer, parameter :: O_RMAX    = 2
integer, parameter :: O_RSTEP   = 3
integer, parameter :: O_SCALE   = 4
integer, parameter :: O_GR      = 5
integer, parameter :: O_IQ      = 6
integer, parameter :: O_SQ      = 7
integer, parameter :: O_FQ      = 8
character(LEN=   5), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 4 ! Number of values to calculate 
!
data oname  / 'rmin ', 'rmax ', 'rstep', 'scale', 'gr', 'iq', 'sq', 'fq' /
data loname /  4     ,  4     ,  5     ,  5     ,  2  ,  2  ,  2  ,  2   /
opara  =  (/ '0.0100000000',  '100.01000000', '0.0100000000', '1.0000000000', 'discus.grobs', 'discus.iqobs', 'discus.sqobs', 'discus.fqobs' /) ! Always provide fresh default values
lopara =  (/  12           ,   12           ,  12           ,  12           ,  12           ,  12           ,  12           ,  12            /)
owerte =  (/  0.0100000000 ,   100.01000000 ,  0.0100000000 , 1.00000       ,  0.0000        ,  0.0000       ,  0.0000       ,  0.0000        /)
!
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_GR  )) then
   exp_outgr   = opara(O_GR)
   exp_outgr_l = .true.
   exp_rmin  = nint(owerte(O_RMIN )*10000.0D0)/10000.0D0
   exp_rmax  = nint(owerte(O_RMAX )*10000.0D0)/10000.0D0
   exp_rstep = nint(owerte(O_RSTEP)*10000.0D0)/10000.0D0
   exp_scale = owerte(O_SCALE)
!else
!   call do_build_name(ianz, cpara, lpara, werte, MAXW, 1) 
!   if(ier_num /= 0) return
!   exp_outfile = cpara(1)
endif
!
if(lpresent(O_IQ  )) then
   exp_outiq    = opara(O_IQ)
   exp_outiq_l = .true.
endif
!
if(lpresent(O_SQ  )) then
   exp_outsq    = opara(O_SQ)
   exp_outsq_l = .true.
endif
!
if(lpresent(O_FQ  )) then
   exp_outfq    = opara(O_FQ)
   exp_outfq_l = .true.
endif
!
!
end subroutine exp2pdf_outfile
!
!*******************************************************************************
!
subroutine exp2pdf_qmax(line, length)
!
! Set the outputfile and the PDF range
!
USE exp2pdf_data_mod
!
use blanks_mod
USE errlist_mod
USE get_params_mod
use take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 4
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
integer :: icomma
!
integer, parameter :: NOPTIONAL = 2
integer, parameter :: O_QINST   = 1
integer, parameter :: O_QFOUR   = 2
character(LEN=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 0 ! Number of values to calculate 
!
data oname  / 'inst', 'fourier' /
data loname /  4     , 7        /
opara  =  (/ '[0.0D0,1.0D9]',  '[0.0D0,1.0D9]' /) ! Always provide fresh default values
lopara =  (/  13    ,   13     /)
owerte =  (/  1.0D9  ,  1.0D9  /)
!
call rem_bl(line,length)
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
exp_qmin_il = .false.
exp_qmax_ul = .false.
exp_qmin_fl = .false.
exp_qmax_fl = .false.
!
if(lpresent(O_QINST ) ) then
   if(opara(O_QINST)(1:1)=='[' .and. opara(O_QINST)(lopara(O_QINST):lopara(O_QINST))==']') then
      icomma = index(opara(O_QINST)(1:lopara(O_QINST)),',')
      if(icomma>0) then                    ! User provided a ','
         if(icomma==2) then                ! User provided upper parameter only [,number]
            exp_qmin_il = .false.
            exp_qmax_ul = .true.
            opara(O_QINST) = '[0.0D0' // opara(O_QINST)(2:lopara(O_QINST))
            lopara(O_QINST)= lopara(O_QINST) + 5
         elseif(icomma==lopara(O_QINST)-1) then       ! User provided lower parameter only [number,]
            exp_qmin_il = .true.
            exp_qmax_ul = .false.
            opara(O_QINST) = opara(O_QINST)(1:lopara(O_QINST)-1) // '1.0D9]'
            lopara(O_QINST)= lopara(O_QINST) + 5
         else                              ! User provided both [number,number]
            exp_qmin_il = .true.
            exp_qmax_ul = .true.
         endif
      else                                 ! User provided upper parameter only [number]
         exp_qmin_il = .false.
         exp_qmax_ul = .true.
         opara(O_QINST) = '[0.0D0,' // opara(O_QINST)(2:lopara(O_QINST)) 
         lopara(O_QINST)= lopara(O_QINST) + 6
      endif
   else                                   ! User provided upper parameter only   number
      exp_qmin_il = .false.
      exp_qmax_ul = .true.
      opara(O_QINST) = '[0.0D0,' // opara(O_QINST)(1:lopara(O_QINST)) // ']'
      lopara(O_QINST)= lopara(O_QINST) + 8
   endif
else                                      ! User provided no parameters at all
   exp_qmin_il = .false.
   exp_qmax_ul = .false.
endif
!
call get_optional_multi(MAXW, opara(O_QINST), lopara(O_QINST), werte, ianz)
if(ier_num/=0) return
!write(*,*) 'LL ', werte(1), ' X ', exp_x(1), exp_x(exp_dim(1)), &
!werte(1) < exp_x(1), werte(1) > exp_x(exp_dim(1)), &
!werte(2) < exp_x(1), werte(2) > exp_x(exp_dim(1))
if(werte(1) >= werte(2)) then   ! Illegal limits
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Lower instrumental limit is'
   ier_msg(2) = 'larger or equal to upper limit'
   return
elseif(werte(1) < exp_x(1) .or. werte(1) > exp_x(exp_dim(1)) ) then   ! Illegal limits
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Lower instrumental limit is'
   ier_msg(2) = 'outside data range            '
   return
elseif(werte(2) < exp_x(1) .or. werte(2) > exp_x(exp_dim(1)) ) then   ! Illegal limits
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Upper instrumental limit is'
   ier_msg(2) = 'outside data range            '
   return
endif
!
exp_qmin_i = werte(1)
exp_qmax_u = werte(2)
!
if(lpresent(O_QFOUR ) ) then
   if(opara(O_QFOUR)(1:1)=='[' .and. opara(O_QFOUR)(lopara(O_QFOUR):lopara(O_QFOUR))==']') then
      icomma = index(opara(O_QFOUR)(1:lopara(O_QFOUR)),',')
      if(icomma>0) then                    ! User provided a ','
         if(icomma==2) then                ! User provided upper parameter only [,number]
            exp_qmin_fl = .false.
            exp_qmax_fl = .true.
            opara(O_QFOUR) = '[0.0D0' // opara(O_QFOUR)(2:lopara(O_QFOUR))
            lopara(O_QFOUR)= lopara(O_QFOUR) + 5
         elseif(icomma==lopara(O_QFOUR)-1) then       ! User provided lower parameter only [number,]
            exp_qmin_fl = .true.
            exp_qmax_fl = .false.
            opara(O_QFOUR) = opara(O_QFOUR)(1:lopara(O_QFOUR)-1) // '1.0D9]'
            lopara(O_QFOUR)= lopara(O_QFOUR) + 5
         else                              ! User provided both [number,number]
            exp_qmin_fl = .true.
            exp_qmax_fl = .true.
         endif
      else                                 ! User provided upper parameter only [number]
         exp_qmin_fl = .false.
         exp_qmax_fl = .true.
         opara(O_QFOUR) = '[0.0D0,' // opara(O_QFOUR)(2:lopara(O_QFOUR)) 
         lopara(O_QFOUR)= lopara(O_QFOUR) + 6
      endif
   else                                   ! User provided upper parameter only   number
      exp_qmin_fl = .false.
      exp_qmax_fl = .true.
      opara(O_QFOUR) = '[0.0D0,' // opara(O_QFOUR)(1:lopara(O_QFOUR)) // ']'
      lopara(O_QFOUR)= lopara(O_QFOUR) + 8
   endif
else                                      ! User provided no parameters at all
   exp_qmin_fl = .false.
   exp_qmax_fl = .false.
endif
!
call get_optional_multi(MAXW, opara(O_QFOUR), lopara(O_QFOUR), werte, ianz)
if(ier_num/=0) return
if(werte(1) >= werte(2)) then   ! Illegal limits
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Lower fourier limit is'
   ier_msg(2) = 'larger or equal to upper limit'
   return
elseif(werte(1)<exp_qmin_i) then
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Lower fourier limit is'
   ier_msg(2) = 'less than lower instrumental limit'
   return
elseif(werte(1)>exp_qmax_u) then
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Lower fourier limit is'
   ier_msg(2) = 'larger than upper instrumental limit'
   return
elseif(werte(2)<exp_qmin_i) then
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Upper fourier limit is'
   ier_msg(2) = 'less than lower instrumental limit'
   return
elseif(werte(2)>exp_qmax_u) then
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Upper fourier limit is'
   ier_msg(2) = 'larger than upper instrumental limit'
   return
elseif(werte(1) < exp_x(1) .or. werte(1) > exp_x(exp_dim(1)) ) then   ! Illegal limits
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Lower fourier limit is'
   ier_msg(2) = 'outside data range            '
   return
elseif(werte(2) < exp_x(1) .or. werte(2) > exp_x(exp_dim(1)) ) then   ! Illegal limits
   ier_num = -61
   ier_typ = ER_FORT
   ier_msg(1) = 'Upper fourier limit is'
   ier_msg(2) = 'outside data range            '
   return
endif
!
exp_qmin_f = werte(1)
exp_qmax_f = werte(2)
!
!
end subroutine exp2pdf_qmax
!
!*******************************************************************************
!
subroutine exp2pdf_poly(line, length)
!
! Set the outputfile and the PDF range
!
USE exp2pdf_data_mod
!
USE errlist_mod
use ber_params_mod
USE get_params_mod
use take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 1
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_ORDER   = 1
character(LEN=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'order' /
data loname /  5      /
opara  =  (/ '9.000' /) ! Always provide fresh default values
lopara =  (/  5          /)
owerte =  (/  9.000  /)
!
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_ORDER) ) then
   exp_npoly   = nint(owerte(O_ORDER))
else 
   call ber_params (ianz, cpara, lpara, werte, maxw)
   exp_npoly   = nint(werte(1))
endif
!
end subroutine exp2pdf_poly
!
!*******************************************************************************
!
subroutine exp2pdf_qscale(line, length)
!
! Set the outputfile and the PDF range
!
USE exp2pdf_data_mod
!
USE errlist_mod
use ber_params_mod
USE get_params_mod
use take_param_mod
use wink_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 10
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
!REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
!
integer, parameter :: NOPTIONAL = 3
integer, parameter :: O_QVAL    = 1
integer, parameter :: O_DVAL    = 2
integer, parameter :: O_QOBS    = 3
character(LEN=   8), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 3 ! Number of values to calculate 
!
data oname  / 'qcrystal' , 'dcrystal', 'qobs' /
data loname /  8         ,  8        ,  4  /
opara  =  (/ '0.000', '0.000', '0.000' /) ! Always provide fresh default values
lopara =  (/  5     ,  5     ,  5      /)
owerte =  (/  0.000 ,  0.000 ,  0.00   /)
!
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_QOBS)) then
   exp_qfirst_o   = owerte(O_QOBS)
   exp_qfirst_l   = .true.
   if(lpresent(O_QVAL) ) then
      exp_qfirst_c   = owerte(O_QVAL)
   elseif(lpresent(O_DVAL)) then
      exp_qfirst_c   = 2.0D0*PI/owerte(O_QVAL)
   endif
else
   exp_qfirst_l   = .false.
   if(lpresent(O_QVAL) ) then
      exp_qfirst_o   = owerte(O_QVAL)
      exp_qfirst_c   = owerte(O_QVAL)
   elseif(lpresent(O_DVAL)) then
      exp_qfirst_o   = 2.0D0*PI/owerte(O_QVAL)
      exp_qfirst_c   = 2.0D0*PI/owerte(O_QVAL)
   endif
endif
!
end subroutine exp2pdf_qscale
!
!*******************************************************************************
!
subroutine exp2pdf_init
!-
!  Initializes flags to write I(Q), S(Q), F(Q) and q-scale correction
!+
!
use exp2pdf_data_mod
!
implicit none
!
exp_outiq_l = .false.
exp_outsq_l = .false.
exp_outfq_l = .false.
!
exp_qmin_il  = .false.
exp_qmin_fl  = .false.
exp_qmax_fl  = .false.
exp_qmax_ul  = .false.
!
exp_qfirst_l =.false.
exp_qfirst_o = 0.0D0
exp_qfirst_c = 0.0D0
!
end subroutine exp2pdf_init
!
!*******************************************************************************
!
end module exp2pdf_supp_mod
