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
SUBROUTINE refine_set(line, length)
!
USE refine_control_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
!
!
LOGICAL, EXTERNAL :: str_comp
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(IANZ>1) THEN
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
IF(str_comp (cpara(1), 'cycles', 3, lpara(1), 6) ) THEN
   cpara(1) = '0'
   lpara(1) = 1
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/= 0) RETURN
   refine_cycles = NINT(werte(2))
ELSEIF(str_comp (cpara(1), 'conver', 3, lpara(1), 6) ) THEN
   CALL refine_set_convergence(line, length)
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
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
!REAL               , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
!
!
LOGICAL, EXTERNAL :: str_comp
!
!
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_DCHI    = 1
INTEGER, PARAMETER :: O_PSHIFT  = 2
INTEGER, PARAMETER :: O_CONF    = 3
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 3 ! Number of values to calculate
!
DATA oname  / 'dchi  ' , 'pshift'  ,  'conf ' /
DATA loname /  4       ,  6        ,   4      /
opara  =  (/ '0.500000', '0.005000',  '0.010000'/)   ! Always provide fresh default values
lopara =  (/  8        ,  8        ,   8        /)
owerte =  (/  0.500000 ,  0.005000 ,   0.010000 /)
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
!
IF(IANZ>=1) THEN
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   IF(ier_num/=0) RETURN
   conv_dchi2    = owerte(O_DCHI)
   conv_dp_sig   = owerte(O_PSHIFT)
   conv_conf     = owerte(O_CONF)
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
END SUBROUTINE refine_set_convergence
!
!*******************************************************************************
END MODULE refine_set_mod
