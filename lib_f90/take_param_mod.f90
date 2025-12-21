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
use build_name_mod
use get_params_mod
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
character(len=PREC_STRING), dimension(:), allocatable :: cfpara    ! optional parameters with format string
integer                   , dimension(:), allocatable :: lfpara    ! optional parameters with format string
real(kind=PREC_DP)        , dimension(:), allocatable ::  fwerte   ! optional parameters with format string
integer                                               :: fanz
!
INTEGER :: i,j, k, iopt, istart
INTEGER :: letter
INTEGER :: icolon
integer :: nper                    ! number of '%' in optional parameter
integer :: maxf                    ! maximum format specifiers found
INTEGER :: len_look, len_user, l0
LOGICAL :: ascii                   ! Test if we have letters only
!
!
lpresent(:) = .FALSE.
IF(ianz==0) RETURN           ! No parameters at all
!
maxf = 0
istart = ianz
search: DO i=istart,1, -1     ! Count backwards
   icolon = INDEX(cpara(i)(1:lpara(i)),':')
   IF(icolon > 0) THEN    ! We might have an optional parameter
      ascii = .TRUE.
      DO j=1, icolon-1    ! Check if all characters are small letters
         letter = IACHAR(cpara(i)(j:j))
         ascii = ascii .AND. ( (a<=letter .AND. letter<=z)  .or.                &
                               (zero<=letter .and. letter<=nine) )
      ENDDO
      IF(.NOT.ascii) THEN
!CYCLE search  ! String contains non-(a..z, 0..9) skip this parameter
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
            nper = 0            ! Assume no '%' in parameter string
!           ! Test if parameter is like "...%..."
            if(opara(iopt)(1:1)=='"' .and. opara(iopt)(lopara(iopt):lopara(iopt))=='"') then
               maxf = 1
               do j=2,lopara(iopt)       ! Search for '%' in parameter string
                  if(opara(iopt)(j:j)=='%') then    ! Found a format specifier
                     maxf = maxf + 1                ! Increment number of format specifiers
                     nper = nper + 1                ! Increment and concatenate
                     opara(iopt) = opara(iopt)(1:lopara(iopt)) //','// cpara(i+nper)(1:lpara(i+nper))
                     lopara(iopt) = len_trim(opara(iopt))
                     cpara(i+nper) = ' '
                     lpara(i+nper) = 0
                  endif
               enddo
            endif
            do k=0, nper           ! Shift format parameters down as well
               DO j = i+1, ianz    ! shift subsequent parameters down
                  cpara(j-1) = cpara(j)
                  lpara(j-1) = lpara(j)
               ENDDO
               ianz = ianz - 1
            enddo
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
!  If parameters ncalc+ start with '"', build name
if(maxf>0) then
   allocate(cfpara(1:maxf))
   allocate(lfpara(1:maxf))
   allocate(fwerte(1:maxf))
!
   loop_form:do i=ncalc+1, NOPTIONAL
      if(opara(i)(1:1) == '"') then      ! Format string
         call get_params(opara(i), fanz, cfpara, lfpara, MAXF, lopara(i))
         if(ier_num /= 0) exit loop_form
         j = 1
         call do_build_name(fanz, cfpara, lfpara, fwerte, MAXF, j)
         if(ier_num /= 0) exit loop_form
         opara(i) = cfpara(j)
         lopara(i) = len_trim(opara(i))
      endif
   enddo loop_form
   deallocate(cfpara)
   deallocate(lfpara)
   deallocate(fwerte)
endif
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
subroutine get_range_multi(lpresent, opara, lopara, iianz, range_low, range_high, LREFINE)
!-
! Set parameter ranges, if "range:" was given
!+
!
use errlist_mod
use ber_params_mod
use get_params_mod
use lib_errlist_func
use precision_mod
!
implicit none
!
logical           , intent(in)  :: lpresent
character(len=*)  , intent(in)  :: opara
integer           , intent(in)  :: lopara
integer           , intent(out) :: iianz
real(kind=PREC_DP), intent(out) :: range_low
real(kind=PREC_DP), intent(out) :: range_high
logical           , intent(in)  :: LREFINE
!
integer, parameter :: MAXF = 2
character(len=PREC_STRING), dimension(2) :: ccpara
character(len=PREC_STRING)               :: string
integer                   , dimension(2) :: llpara
integer :: length
logical, dimension(2) :: lrange
real(kind=PREC_DP)        , dimension(2) :: wwerte

range_low  = +1.0
range_high = -1.0
IF(lpresent) THEN
   IF(opara(1:1) == '[' .AND. opara(lopara:lopara) == ']' .AND. &
      INDEX(opara(2:lopara-1),',')>0) THEN
      string = ' '
      string(1:lopara-2) = opara(2:lopara-1)
      length = lopara-2
      ccpara(:) = ' '
      llpara(:) = 0
      wwerte(:) = 0.0
      CALL get_params (string, iianz, ccpara, llpara, MAXF, length)
      lrange = .TRUE.
      IF(llpara(1)==0 .AND. ier_num==-2 .AND. ier_typ==ER_COMM) THEN
         lrange(1) = .FALSE.
         CALL no_error
         ccpara(1) = '0'
         llpara(1) = 1
         iianz = 2
      ELSEIF(llpara(2)==0 .AND. ier_num==-2 .AND. ier_typ==ER_COMM) THEN
         lrange(2) = .FALSE.
         CALL no_error
         ccpara(2) = '0'
         llpara(2) = 1
         iianz = 2
      ELSEIf(ier_num/=0) THEN
!     IF(ier_num /= 0) THEN
         ier_msg(1) = 'Incorrect ''range:[]'' parameter'
         RETURN
      ENDIF
      CALL ber_params (iianz, ccpara, llpara, wwerte, MAXF)
      IF(ier_num /= 0) THEN
         ier_msg(1) = 'Incorrect ''range:[]'' parameter'
!        ier_msg(2) = 'Variables can only be arrays with'
!        ier_msg(3) = 'one or two dimensions '
         RETURN
      ENDIF
      IF(iianz==2) THEN
         range_low  = wwerte(1)
         range_high = wwerte(2)
         IF(.NOT.lrange(1)) range_low  = -HUGE(0.0)
         IF(.NOT.lrange(2)) range_high =  HUGE(0.0)
         if(LREFINE) then                  ! Limits for REFINE
            IF(range_low>=range_high) THEN
               ier_num = -6
               ier_typ = ER_FORT
               ier_msg(1) = 'Incorrect ''range:[]'' parameter'
               ier_msg(2) = 'Low > = High    '
               RETURN
            ENDIF
         else                              ! Limits for DIFFEV
            IF(range_low>range_high) THEN
               ier_num = -6
               ier_typ = ER_FORT
               ier_msg(1) = 'Incorrect ''range:[]'' parameter'
               ier_msg(2) = 'Low > High    '
               RETURN
            ENDIF
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'Incorrect ''range:[]'' parameter'
         RETURN
      ENDIF
   ELSE
      ier_msg(1) = 'Incorrect ''range:[]'' parameter'
      RETURN
   ENDIF
ENDIF
!
end subroutine get_range_multi
!
!*******************************************************************************
!
END MODULE take_param_mod
