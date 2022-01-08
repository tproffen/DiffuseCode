MODULE take_param_mod
!
!  Combines all routines that handle the command parameters within the suite
!  At the moment only the code for optional parameters, 
!  get_params  ! gets the parameters from the input line, currently in dkdo.f90
!  ber_params  ! calculates the parameters from the input line, currently in dkdo.f90
!
CONTAINS
!
SUBROUTINE get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL, ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
!
!  Takes any optional parameter out of the list of input parameters.
!  These optional parameters are copied into the array opara and the 
!  original list is cleaned of these optional parameters.
!  An optional parameter takes the form:
!  name:value
!  Here name is a string supplied by the user, and value a string that
!  gives the value. To handle numerical and character values, the 
!  numerical values must always be the first ncalc parameters. The
!  remaining parameters will not be calculated, and the routine returns
!  the strings only.
!
USE ber_params_mod
use charact_mod
USE errlist_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT NONE
!
INTEGER                               , INTENT(INOUT) :: ianz      ! Actual parameter number
INTEGER                               , INTENT(IN)    :: MAXW      ! Max para numbeer
CHARACTER(LEN=*), DIMENSION(1:MAXW)   , INTENT(INOUT) :: cpara     ! input strings
INTEGER,          DIMENSION(1:MAXW)   , INTENT(INOUT) :: lpara     ! input string lengths
INTEGER                               , INTENT(IN)    :: NOPTIONAL ! No opt. params
INTEGER                               , INTENT(IN)    :: ncalc     ! No of params to calculate
CHARACTER(LEN=*), DIMENSION(NOPTIONAL), INTENT(IN)    :: oname     ! Lookup table
INTEGER,          DIMENSION(NOPTIONAL), INTENT(IN)    :: loname    ! lookup table length
CHARACTER(LEN=*), DIMENSION(NOPTIONAL), INTENT(INOUT) :: opara     ! with default values
INTEGER,          DIMENSION(NOPTIONAL), INTENT(INOUT) :: lopara    ! length of results
LOGICAL,          DIMENSION(NOPTIONAL), INTENT(OUT)   :: lpresent  ! Is param present ?
REAL(KIND=PREC_DP),DIMENSION(NOPTIONAL), INTENT(INOUT) :: owerte    ! calc. results with default values
!
INTEGER :: i,j, iopt, istart
INTEGER :: letter
INTEGER :: icolon
INTEGER :: len_look, len_user, l0
LOGICAL :: ascii                   ! Test if we have letters only
!
!
lpresent(:) = .FALSE.
IF(ianz==0) RETURN           ! No parameters at all
!
istart = ianz
search: DO i=istart,1, -1     ! Count backwards
   icolon = INDEX(cpara(i)(1:lpara(i)),':')
   IF(icolon > 0) THEN    ! We might have an optional parameter
      ascii = .TRUE.
      DO j=1, icolon-1    ! Check if all characters are small letters
         letter = IACHAR(cpara(i)(j:j))
         ascii = ascii .AND. (a<=letter .AND. letter<=z)
      ENDDO
      IF(.NOT.ascii) THEN
!CYCLE search  ! String contains non-(a..z) skip this parameter
         ier_num = -12
         ier_typ = ER_COMM
         ier_msg(1) = 'Offending parameter name:'
         ier_msg(2) = cpara(i)(1:MIN(43,lpara(i)))
         RETURN
      ENDIF
      look: DO iopt=1, NOPTIONAL ! Look up optional parameter name
         len_look = loname(iopt)
         len_user = icolon-1
         l0 = MIN(len_look,len_user)
         IF(str_comp(cpara(i)(1:len_user), oname(iopt), l0, len_user, len_look)) THEN  ! Found parameter
            opara(iopt)  = cpara(i)(icolon+1:lpara(i))  ! Copy user provided string
            lopara(iopt) = lpara(i)-icolon              ! record user provided string length
            lpresent(iopt) = .TRUE.
            cpara(i) = ' '      ! clear parameter, this avoids issue for parameter ianz
            lpara(i) = 0
            DO j = i+1, ianz    ! shift subsequent parameters down
               cpara(j-1) = cpara(j)
               lpara(j-1) = lpara(j)
            ENDDO
            ianz = ianz - 1
            CYCLE search
         ENDIF
      ENDDO look
      ier_num = -12
      ier_typ = ER_COMM
      ier_msg(1) = 'Offending parameter name:'
      ier_msg(2) = cpara(i)(1:MIN(43,lpara(i)))
      RETURN
   ENDIF
ENDDO search
!
!  Calculate numerical values for the first ncalc parameters
!
CALL ber_params (ncalc, opara, lopara, owerte, NOPTIONAL)
!
!
END SUBROUTINE get_optional
!
!*******************************************************************************
!
SUBROUTINE get_optional_multi(MAXW, opara, lopara, werte, ianz)
!-
! Calculate multiple values for the optional parameters as for "dim:[3,3]"
!+
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: MAXW     ! Dimension of werte
CHARACTER(LEN=*)                   , INTENT(INOUT) :: opara    ! The string with optional values
INTEGER                            , INTENT(INOUT) :: lopara   ! length of string
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(OUT)   :: werte    ! Numerical values
INTEGER                            , INTENT(OUT)   :: ianz    ! Number of numerical values
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(opara))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER                              :: length
!
IF(opara(1:1)=='[' .AND. opara(lopara:lopara)==']') THEN     ! Matching "[...]"
   length = lopara - 2                                       ! search string is shorter by "[" and "]"
   CALL get_params (opara(2:lopara-1), ianz, cpara, lpara, maxw, length)
   IF (ier_num.ne.0) RETURN
   CALL ber_params (ianz, cpara, lpara, werte, MAXW)
   IF (ier_num.ne.0) RETURN
ELSE
   ier_num = -9
   ier_typ = ER_FORT
   ier_msg(1) = 'Multiple optional values must be '
   ier_msg(2) = 'enclosed by []'
ENDIF
!
END SUBROUTINE get_optional_multi
!
!*******************************************************************************
!
subroutine eva_optional_multi(opara, lopara, MAXWW, ccpara, wwerte, nums, nnames)
!-
!  evaluate an optional parameter that might consist of numerical values and/or character strings
!+
!
use ber_params_mod
use errlist_mod
use precision_mod
!
implicit none
!
character(len=*)  , intent(in) :: opara     ! The optional parameter string
integer           , intent(in) :: lopara    ! its length
integer           , intent(in) :: MAXWW     ! Dimensions of the output arrays
character(len=*)  , dimension(MAXWW), intent(out) ::ccpara     ! The optional parameter string
real(kind=PREC_DP), dimension(MAXWW), intent(out) :: wwerte    ! The optional parameter string
integer                             , intent(out) :: nums
integer                             , intent(out) :: nnames
!
integer, parameter :: MAXW = 10
character(len=PREC_STRING), dimension(MAXW) :: cpara     ! Temporary array for lower routines
integer                   , dimension(MAXW) :: lpara     ! Temporary array for lower routines
!real(kind=PREC_DP)        , dimension(MAXW) :: werte     ! Temporary array for lower routines
integer :: ianz
integer :: i
!
if(opara(1:1)=='[') then     ! We seem to have a [1,2] style
   call sep_optional_multi(MAXW, opara, lopara, cpara, ianz)
   do i=1,ianz
      lpara(i) = len_trim(cpara(i))
   enddo
   if(ier_num/=0) then
      ier_msg(1) = 'Optional parameters erroneous'
      return
   endif
else
   cpara(1) = opara
   lpara(1) = lopara
   ianz = 1
endif
nums = 0
nnames = 0
call eva_params(ianz, cpara, lpara, wwerte, ccpara, MAXW, nums, nnames)
!
end subroutine eva_optional_multi
!
!*******************************************************************************
!
SUBROUTINE sep_optional_multi(MAXW, opara, lopara, cpara, ianz)
!-
! Separates multiple values for the optional parameters as for "dim:[3,3]"
!+
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: MAXW     ! Dimension of cpara
CHARACTER(LEN=*)                   , INTENT(IN   ) :: opara    ! The string with optional values
INTEGER                            , INTENT(IN   ) :: lopara   ! length of string
CHARACTER(LEN=*)  , DIMENSION(MAXW), intent(out)   :: cpara
INTEGER                            , INTENT(OUT)   :: ianz    ! Number of numerical values
!
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER                              :: length
!
IF(opara(1:1)=='[' .AND. opara(lopara:lopara)==']') THEN     ! Matching "[...]"
   length = lopara - 2                                       ! search string is shorter by "[" and "]"
   CALL get_params (opara(2:lopara-1), ianz, cpara, lpara, maxw, length)
   IF (ier_num.ne.0) RETURN
ELSE
   ier_num = -9
   ier_typ = ER_FORT
   ier_msg(1) = 'Multiple optional values must be '
   ier_msg(2) = 'enclosed by []'
ENDIF
!
END SUBROUTINE sep_optional_multi
!
!*******************************************************************************
!
END MODULE take_param_mod
