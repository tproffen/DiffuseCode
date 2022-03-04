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
USE ber_params_mod
USE get_params_mod
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
character(len=2)           :: atom_name
character(len=PREC_STRING) :: string
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
integer                              :: i, j
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
      call rem_bl(opara(O_COMP), lopara(O_COMP))
      string = opara(O_COMP)
      length = lopara(O_COMP)
   endif
else                          ! User omitted "comp:"
   length = len_trim(string)
   call rem_bl(string, length)
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
i = length
exp_natom = 0
loop_search: do
   if(i==0) exit loop_search
   j = i               ! End of a possible number
   loop_number: do
      if(i==0) exit loop_number
      if((iachar(string(i:i))>= zero .and. iachar(string(i:i))<= nine) .or.  &
          string(i:i)=='.' ) then
         i = i-1
      else
         if(j-i>0) then       ! Found significant number string
            cpara(1) = '(' // string(i+1:j) // ')'
            lpara(1) = j-i + 2
         else
            cpara(1) = '(1.0)'
            lpara(1) = 5
         endif
         exit loop_number
      endif
   enddo loop_number
      ianz = 1
      call ber_params (ianz, cpara, lpara, werte, maxw)
   j = i
   if(iachar(string(i:i))>= a .and. iachar(string(i:i))<= z) then   ! Lower case letter
     i = i-1
   endif
   if(iachar(string(i:i))>= aa .and. iachar(string(i:i))<= zz) then   ! Upper case letter
      atom_name = string(i:j)
   endif
   exp_natom= exp_natom + 1
   call do_cap(atom_name)
   exp_atname(exp_natom) = atom_name
   exp_atocc(exp_natom ) = werte(1)
   i = i-1
enddo loop_search
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
use str_comp_mod
!
implicit none
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
if(str_comp(line, 'xray',  2, length, 4) ) then
   exp_radiation = 'xray'
elseif(str_comp(line, 'electron',  2, length, 8) ) then
   exp_radiation = 'electron'
elseif(str_comp(line, 'neutron',  2, length, 7) ) then
   exp_radiation = 'neutron'
else
   ier_num = -183
   ier_typ = ER_APPL
   ier_msg(1) = 'Must be either ''xray'', ''neutron'' or '
   ier_msg(2) = '               ''electron'' '
endif
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
integer, parameter :: NOPTIONAL = 7
integer, parameter :: O_RMIN    = 1
integer, parameter :: O_RMAX    = 2
integer, parameter :: O_RSTEP   = 3
integer, parameter :: O_GR      = 4
integer, parameter :: O_IQ      = 5
integer, parameter :: O_SQ      = 6
integer, parameter :: O_FQ      = 7
character(LEN=   5), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 3 ! Number of values to calculate 
!
data oname  / 'rmin ', 'rmax ', 'rstep', 'gr', 'iq', 'sq', 'fq' /
data loname /  4     ,  4     ,  5     ,  2  ,  2  ,  2  ,  2   /
opara  =  (/ '0.0100000000',  '100.01000000', '0.0100000000', 'discus.grobs', 'discus.iqobs', 'discus.sqobs', 'discus.fqobs' /) ! Always provide fresh default values
lopara =  (/  12           ,   12           ,  12           ,  12           ,  12           ,  12           ,  12            /)
owerte =  (/  0.0100000000 ,   100.01000000 ,  0.0100000000 , 0.0000        ,  0.0000       ,  0.0000       ,  0.0000        /)
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
USE errlist_mod
USE get_params_mod
use take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 3
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
!REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
!
integer, parameter :: NOPTIONAL = 3
integer, parameter :: O_QMAX_I  = 1
integer, parameter :: O_QMAX_F  = 2
integer, parameter :: O_QMIN    = 3
character(LEN=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 3 ! Number of values to calculate 
!
data oname  / 'inst', 'fourier', 'qmin' /
data loname /  4     , 7       ,  4     /
opara  =  (/ '1.0D9',  '1.0D9', '0.0D0' /) ! Always provide fresh default values
lopara =  (/  5     ,   5     ,  5      /)
owerte =  (/  1.0D9  ,  1.0D9 ,  0.0D0  /)
!
call get_params(line, ianz, cpara, lpara, MAXW, length)
if(ier_num/=0) return
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
if(ier_num/=0) return
!
if(lpresent(O_QMAX_I) ) then
   exp_qmax_u  = owerte(O_QMAX_I)
   exp_qmax_ul = .true. 
endif
!
if(lpresent(O_QMAX_F) ) then
   exp_qmax_f  = owerte(O_QMAX_F)
   exp_qmax_fl = .true. 
endif
!
if(lpresent(O_QMIN  ) ) then
   exp_qmin_u  = owerte(O_QMIN  )
endif
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
opara  =  (/ '7.000' /) ! Always provide fresh default values
lopara =  (/  5          /)
owerte =  (/  7.000  /)
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
end module exp2pdf_supp_mod