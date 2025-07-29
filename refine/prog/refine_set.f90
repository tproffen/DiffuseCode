MODULE refine_set_mod
!
!  Routines to set various user settings
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_set(string, str_length)
!
USE refine_control_mod
use refine_log_mod
!
USE errlist_mod
USE ber_params_mod
use do_set_mod
USE get_params_mod
USE precision_mod
USE str_comp_mod
use prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER         , INTENT(IN) :: str_length
!
INTEGER, PARAMETER :: MAXW = 20
!
character(len=PREC_STRING) :: line
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
integer                              :: length
!
line = string
length = str_length
!
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
!IF(IANZ>1) THEN
!ELSE
!   ier_num = -6
!   ier_typ = ER_FORT
!   RETURN
!ENDIF
!
IF(str_comp (cpara(1), 'cycles', 3, lpara(1), 6) ) THEN
   cpara(1) = '0'
   lpara(1) = 1
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/= 0) RETURN
   refine_cycles = NINT(werte(2))
ELSEIF(str_comp (cpara(1), 'convergence', 3, lpara(1), 11) ) THEN
   CALL refine_set_convergence(line, length)
ELSEIF(str_comp (cpara(1), 'relax', 3, lpara(1), 5) ) THEN
   CALL refine_set_lamda(line, length)
!  cpara(1) = '0'
!  lpara(1) = 1
!  CALL ber_params(ianz, cpara, lpara, werte, MAXW)
!  IF(ier_num/= 0) RETURN
!  refine_lamda_s  = werte(2)
!  IF(ianz>=3) THEN
!     refine_lamda_u  = werte(3)
!     IF(ianz==4) THEN
!        refine_lamda_d  = werte(4)
!     ENDIF
!  ENDIF
elseif(str_comp(cpara(1), 'log', 3, lpara(1), 3)) then
   refine_log = str_comp(cpara(2), 'on', 2, lpara(1), 2)
ELSE
   line = string
   length = str_length
   call do_set(line, length)
!  ier_num = -8
!  ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE refine_set
!
!*******************************************************************************
!
SUBROUTINE refine_set_convergence(line, length)
!
USE refine_control_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE str_comp_mod
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
INTEGER            , DIMENSION(MAXW) :: lpara
!
INTEGER                              :: ianz
!
!
!
!
INTEGER, PARAMETER :: NOPTIONAL = 6
INTEGER, PARAMETER :: O_DCHI    = 1
INTEGER, PARAMETER :: O_PSHIFT  = 2
INTEGER, PARAMETER :: O_CONF    = 3
INTEGER, PARAMETER :: O_CHI     = 4
INTEGER, PARAMETER :: O_LAMBDA  = 5
INTEGER, PARAMETER :: O_STATUS  = 6
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 5 ! Number of values to calculate
logical, save :: INIT = .true.
!
!
DATA oname  / 'dchi  ' , 'pshift'  ,  'conf '   ,  'chisq ' , 'lambda'  , 'status'  /
DATA loname /  4       ,  6        ,   4        ,  5        ,  6        ,  6        /
opara  =  (/ '0.500000', '0.005000',  '0.010000', '0.500000', '1.0E10  ', 'on      '/)   ! Always provide fresh default values
lopara =  (/  8        ,  8        ,   8        ,  8        ,  8        ,  8        /)
owerte =  (/  0.500000 ,  0.005000 ,   0.010000 ,  0.005000 ,  1.0E10   ,  0.000000 /)
if(INIT) then
   conv_dchi2    = owerte(O_DCHI)
   conv_dp_sig   = owerte(O_PSHIFT)
   conv_conf     = owerte(O_CONF)
   conv_chi2     = owerte(O_CHI)
   conv_lambda   = owerte(O_LAMBDA)
   INIT          = .false.
endif
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
!
IF(IANZ>=1) THEN
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   IF(ier_num/=0) RETURN
   if(lpresent(O_STATUS)) conv_status   = str_comp(opara(O_STATUS), 'on', 2, lopara(1), 2)
   if(lpresent(O_DCHI))   then
      conv_dchi2   = owerte(O_DCHI)
      conv_dchi2_u = .TRUE.
   else
      conv_dchi2_u = .FALSE.
   endif
   if(lpresent(O_PSHIFT)) then
      conv_dp_sig   = owerte(O_PSHIFT)
      conv_dp_sig_u = .TRUE.
   else
      conv_dp_sig_u = .FALSE.
   endif
   if(lpresent(O_CONF))   then
      conv_conf   = owerte(O_CONF)
      conv_conf_u = .TRUE.
   else
      conv_conf_u = .FALSE.
   endif
   if(lpresent(O_CHI))    then
      conv_chi2   = owerte(O_CHI)
      conv_chi2_u = .TRUE.
   else
      conv_chi2_u = .FALSE.
   endif
   if(lpresent(O_LAMBDA))    then
      conv_lambda = owerte(O_LAMBDA)
      conv_lamb_u = .TRUE.
   else
      conv_lamb_u = .FALSE.
   endif
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
END SUBROUTINE refine_set_convergence
!
!*******************************************************************************
!
SUBROUTINE refine_set_lamda(line, length)
!
USE refine_control_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
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
INTEGER            , DIMENSION(MAXW) :: lpara
!
INTEGER                              :: ianz
!
!
!
!
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_START   = 1
INTEGER, PARAMETER :: O_FAIL    = 2
INTEGER, PARAMETER :: O_SUCCESS = 3
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 3 ! Number of values to calculate
!
DATA oname  / 'start ' , 'fail  '  ,  'success' /
DATA loname /  5       ,  4        ,   7      /
opara  =  (/ '0.020000', '16.00000',  '0.500000'/)   ! Always provide fresh default values
lopara =  (/  8        ,  8        ,   8        /)
owerte =  (/  0.020000 ,  16.00000 ,   0.500000 /)
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
!
IF(IANZ>=1) THEN
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   IF(ier_num/=0) RETURN
   refine_lamda_s   = owerte(O_START)
   refine_lamda_s_u = lpresent(O_START)
!
   refine_lamda_u   = owerte(O_FAIL)
   refine_lamda_u_u = lpresent(O_FAIL)
!
   refine_lamda_d   = owerte(O_SUCCESS)
   refine_lamda_d_u = lpresent(O_SUCCESS)
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
END SUBROUTINE refine_set_lamda
!
!*******************************************************************************
END MODULE refine_set_mod
