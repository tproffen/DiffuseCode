MODULE refine_add_param_mod
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_add_param(line, length)
!
!  Add a new parameter. values of existing parameters are updated.
!  If an optional parameter is omitted, general defaults are applied.
!  Sepcial parameter names like "P_lat" etc get special defaults.
!
USE refine_allocate_appl
USE refine_params_mod
!
USE errlist_mod
USE define_variable_mod
USE get_params_mod
USE ber_params_mod
USE calc_expr_mod
USE lib_errlist_func
USE precision_mod
USE take_param_mod
!
USE global_data_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE.
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
!REAL(kind=PREC_DP), DIMENSION(MAXW) :: werte   ! Might be used later on
!
INTEGER                              :: n_params    ! temp number of parameters
INTEGER                              :: ianz, iianz
INTEGER                              :: i, j
LOGICAL                              :: is_new = .FALSE.  ! is a new parameter name
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line)))                  :: string
CHARACTER(LEN=  16)                  :: pname
INTEGER                              :: lpname
INTEGER                              :: laenge
INTEGER                              :: indxg
INTEGER                              :: ipar         ! enty number for this parameter
!
LOGICAL                              :: lrefine
!
REAL(kind=PREC_DP)                   :: range_low    ! template for parameter
REAL(kind=PREC_DP)                   :: range_high   ! ranges 
LOGICAL, DIMENSION(2)                :: lrange       ! Range is provided
!
INTEGER, PARAMETER :: MAXF=2
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXF) :: ccpara
INTEGER            , DIMENSION(MAXF) :: llpara
REAL(KIND=PREC_DP) , DIMENSION(MAXF) :: wwerte
REAL(KIND=PREC_DP) :: temp_val     ! temporary value from value: statement
!
INTEGER, PARAMETER :: NOPTIONAL = 5
INTEGER, PARAMETER :: OSHIFT    = 1
INTEGER, PARAMETER :: ONDERIV   = 2
INTEGER, PARAMETER :: OVALUE    = 3
INTEGER, PARAMETER :: OSTATUS   = 4
INTEGER, PARAMETER :: ORANGE    = 5
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate
!
DATA oname  / 'shift'  , 'points' ,  'value ' , 'status'  ,  'range' /
DATA loname /  5       ,  6       ,   5       ,  6        ,   5      /
opara  =  (/ '0.003000', '3.000000', '-1.00000', 'free    ',  '0.000000'/)   ! Always provide fresh default values
lopara =  (/  8        ,  8        ,  8        ,  8        ,   8        /)
owerte =  (/  0.00300  ,  3.0      ,  -1.0     ,  0.0      ,   0.0      /)
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num/=0) RETURN
!
IF(ianz/=1) THEN
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
pname  = ' '
pname  = cpara(1)(1:MIN(lpara(1),LEN(pname)))
lpname = MIN(lpara(1),LEN(pname))
cpara(1) = '0'
lpara(1) = 1
!
string = 'real, '//pname
laenge = 6+lpname
CALL define_variable(string, laenge, IS_DIFFEV)           ! Define as user variable
IF(ier_num/=0) THEN
   ier_msg(1) = 'Could not define variable name'
   RETURN
ENDIF
!
! Check shift for special parameter names
!
if_shift: if(.not.lpresent(OSHIFT)) then                            ! User did not specify shift:
   do i=1,REF_MAXPARAM_SPC                                ! Compare to special parameters
      if(pname             (1:len_trim(refine_spc_name(i))) ==                  &
         refine_spc_name(i)(1:len_trim(refine_spc_name(i)))) then
         owerte(OSHIFT) = refine_spc_shift(i)             ! Take default value for special parameter
         exit if_shift
      endif
   enddo
endif if_shift
!
! Set starting value for parameter, if 'value:' was given
!
temp_val = 0.0D0
IF(lpresent(OVALUE)) THEN
   IF(opara(OVALUE)/='current') THEN                      ! User did not specify "current"
      ccpara(1) = opara(OVALUE)(1:lopara(OVALUE))
      llpara(1) = lopara(OVALUE)
      iianz = 1
      CALL ber_params (iianz, ccpara, llpara, wwerte, MAXF)
      IF(ier_num/=0) THEN
         ier_msg(1) = 'Could not set parameter value'
         RETURN
      ENDIF
      temp_val = wwerte(1)
      WRITE(string,'(a,a,a)') pname(1:lpname), ' = ', opara(OVALUE)(1:lopara(OVALUE))
      indxg = lpname + 2
      CALL do_math (string, indxg, length)
      IF(ier_num/=0) THEN
         ier_msg(1) = 'Could not set parameter value'
         RETURN
      ENDIF
   ENDIF
ENDIF
!
! Check number of derivative points "points:"
!
IF(lpresent(ONDERIV)) THEN                            ! User provided "points:"
   IF(.NOT. (NINT(owerte(ONDERIV))==3 .OR.         &
             NINT(owerte(ONDERIV))==5     ) ) THEN
      ier_num = -8
      ier_typ = ER_APPL
      RETURN
   ENDIF
else
   loop_nderiv: do i=1,REF_MAXPARAM_SPC                                ! Compare to special parameters
      if(pname             (1:len_trim(refine_spc_name(i))) ==                  &
         refine_spc_name(i)(1:len_trim(refine_spc_name(i)))) then
         owerte(ONDERIV) = real(refine_spc_nderiv(i),kind=PREC_DP)     ! Take default value for special parameter
         exit loop_nderiv
      endif
   enddo loop_nderiv
ENDIF
!
! Set parameter ranges, if "range:" was given
!
range_low  = +1.0
range_high = -1.0
IF(lpresent(ORANGE)) THEN
   IF(opara(ORANGE)(1:1) == '[' .AND. opara(ORANGE)(lopara(ORANGE):lopara(ORANGE)) == ']' .AND. &
      INDEX(opara(ORANGE)(2:lopara(ORANGE)-1),',')>0) THEN
      string = ' '
      string(1:lopara(ORANGE)-2) = opara(ORANGE)(2:lopara(ORANGE)-1)
      length = lopara(ORANGE)-2
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
         IF(range_low>=range_high) THEN
            ier_num = -6
            ier_typ = ER_FORT
            ier_msg(1) = 'Incorrect ''range:[]'' parameter'
            RETURN
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
! A parameter is added to the list of refined parameters only if
! the status is set to 'refine' or 'free'. Otherwise it is omited
! from the list of parameters.
! Fixed parameters are added to the list refine_fixed instead.
!
lrefine = .TRUE.
IF(lpresent(OSTATUS)) THEN
  IF(opara(OSTATUS)=='refine' .OR. opara(OSTATUS)=='free') THEN
      lrefine = .TRUE.
   ELSEIF(opara(OSTATUS)=='fix' .OR. opara(OSTATUS)=='fixed') THEN
      lrefine = .FALSE.
   ELSE
      ier_num = -6
      ier_typ = ER_FORT
      ier_msg(1) = 'Unknown keyword for ''status:'' '
      RETURN
   ENDIF
ELSE
   lrefine = .TRUE.
ENDIF
!
ipar = 1
IF(lrefine) THEN
   is_new = .TRUE.
   old: DO i=1, refine_par_n
      IF(pname == refine_params(i)) THEN        ! Found old parameter name
         ipar = i
         is_new = .FALSE.
         EXIT old
      ENDIF
   ENDDO old
!
   IF(is_new) THEN                              ! New parameter, add to list
      IF(refine_par_n==REF_MAXPARAM) THEN
         n_params = REF_MAXPARAM + 10
         CALL alloc_params(n_params)
      ENDIF
      refine_par_n = refine_par_n + 1
      refine_params(refine_par_n)   = pname
      ipar = refine_par_n
   ENDIF
   refine_range(ipar,1) = range_low
   refine_range(ipar,2) = range_high
   refine_shift(ipar)   = owerte(OSHIFT)
   refine_nderiv(ipar)  = owerte(ONDERIV)
   IF(opara(OVALUE)/='current') THEN                      ! User did not specify "current"
!      READ(opara(OVALUE)(1:lopara(OVALUE)), *) refine_p(ipar)
      refine_p(ipar) = temp_val
   ENDIF
!   write(*,'(a,a,f12.6,i3,f12.6)') 'FREE   ',refine_params(ipar), refine_p(ipar), &
!   refine_nderiv(ipar), refine_shift(ipar)
!
   fixed: DO i=1, refine_fix_n                 ! Remove from fixed list
      IF(pname == refine_fixed(i)) THEN        ! Found old parameter name
         IF(i==refine_fix_n) THEN              ! This is the last parameter
            refine_fixed(i) = ' '
            refine_f    (i) = 1.0D0
            refine_shift_fix(i) = owerte(OSHIFT)
            refine_nderiv_fix(i) = owerte(ONDERIV)
            refine_range_fix(ipar,1) = range_low
            refine_range_fix(ipar,2) = range_high
            refine_fix_n = refine_fix_n - 1
         ELSE                                  ! This is not the last parameter
            DO j=i+1, refine_fix_n
               refine_fixed(j-1) = refine_fixed(j)
               refine_f    (j-1) = refine_f    (j)
               refine_nderiv_fix(j-1) = refine_nderiv_fix(j)
               refine_shift_fix(j-1)  = refine_shift_fix(j)
               refine_range_fix(j-1,:)  = refine_range_fix(j,:)
            ENDDO
            refine_fix_n = refine_fix_n - 1
            EXIT fixed
         ENDIF
      ENDIF
   ENDDO fixed
ELSE
   is_new = .TRUE.
   old_f: DO i=1, refine_fix_n
      IF(pname == refine_fixed(i)) THEN        ! Found old parameter name
         ipar = i
         is_new = .FALSE.
         EXIT old_f
      ENDIF
   ENDDO old_f
   IF(opara(OVALUE)/='current') THEN                      ! User did not specify "current"
!     READ(opara(OVALUE)(1:lopara(OVALUE)), *) refine_p(ipar)
      refine_p(ipar) = temp_val
   ENDIF
!
   IF(is_new) THEN                              ! New parameter, add to list
      IF(refine_fix_n==REF_MAXPARAM_FIX) THEN
         n_params = REF_MAXPARAM_FIX + 10
         CALL alloc_params_fix(n_params)
      ENDIF
      refine_fix_n = refine_fix_n + 1
      refine_fixed(refine_fix_n) = pname
      refine_f(refine_fix_n)     = temp_val
      refine_shift_fix(refine_fix_n) = owerte(OSHIFT)
      refine_nderiv_fix(refine_fix_n) = owerte(ONDERIV)
      refine_range_fix(refine_fix_n,1)= range_low
      refine_range_fix(refine_fix_n,2)= range_high
   ENDIF
ENDIF
!
CALL gl_set_pnumber(REF_MAXPARAM, REF_MAXPARAM_FIX, refine_par_n, refine_fix_n, &
     refine_params, refine_fixed)
!!do i=1, max(REF_MAXPARAM,refine_par_n)
!   write(*,'(a,a,f12.6)') 'FREE   ',refine_params(i), refine_p(i), &
!   refine_nderiv(i), refine_shift(i)
!enddo
!do i=1, max(REF_MAXPARAM_FIX,refine_fix_n)
!  write(*,'(a,a,f12.6)') 'FIXED  ',refine_fixed(i), refine_f(i)
!enddo
!
END SUBROUTINE refine_add_param
!
!*******************************************************************************
!
END MODULE refine_add_param_mod
