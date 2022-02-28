MODULE diffev_release_mod
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE do_release(zeile, length)
!
USE compare
USE initialise
USE population
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(INOUT) :: zeile
INTEGER          , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: maxw = 6
!
CHARACTER (LEN=    PREC_STRING            ), DIMENSION(MAXW) :: cpara   = ' '
INTEGER             , DIMENSION(MAXW) :: lpara = 0
REAL(KIND=PREC_DP)  , DIMENSION(MAXW) :: werte = 0.0
!
INTEGER                               :: ianz
INTEGER                               :: k
INTEGER                               :: lb, d_lb, d_ub
LOGICAL                               :: lexist
LOGICAL                               :: l_init_x
LOGICAL                               :: ldismiss
REAL(KIND=PREC_DP)                    :: set_value
REAL(KIND=PREC_DP)                    :: set_xmin
REAL(KIND=PREC_DP)                    :: set_xmax
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 5
INTEGER, PARAMETER :: O_RANGE = 1
INTEGER, PARAMETER :: O_VALUE = 2
INTEGER, PARAMETER :: O_MIN   = 3
INTEGER, PARAMETER :: O_MAX   = 4
INTEGER, PARAMETER :: O_DISMISS=5
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=    PREC_STRING            ), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Length opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Length opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 4 ! Number of values to calculate 
!
DATA oname  / 'range '  , 'value '   , 'min    '  , 'max   '  , 'dismiss'  /
DATA loname /  6        ,  6         ,  6         ,  6        ,  7         /
opara  =  (/ '-9999.00' , '-9999.00' , '-9999.00' , '-9999.00', '0.000000' /)   ! Always provide fresh default values
lopara =  (/  8         ,  8         , 8          , 8         ,  8         /)
owerte =  (/  -9999.    ,  -9999.    , -9999.0    , -9999.0   ,  0.0       /)
!
INQUIRE(FILE='GENERATION', EXIST=lexist)
IF(lexist) THEN                ! A GENERATION FILE EXISTS
   CALL do_read_values(.TRUE.) ! We need to read values as this can be the first command after a continue
ENDIF
IF((pop_gen==0 .AND. .NOT. pop_initialized) .OR. .NOT.lexist) THEN   ! Population was not yet initialized
   IF(pop_trialfile == ' ') pop_trial_file_wrt= .FALSE.
   CALL do_initialise (l_init_x)
   pop_initialized = .TRUE.
ENDIF
IF(ier_num/=0) THEN
   ier_msg(1) = 'Check the GENERATION, the parameter and'
   ier_msg(2) = 'the lastfile for conflicting generation values'
   RETURN
ENDIF
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, length)
IF (ier_num == 0) THEN
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   IF(ier_num==0) THEN
   IF(ianz == 1) THEN
      IF(owerte(O_RANGE) == -9999.00) THEN   ! range: was not provided
         ier_num = -34
         ier_typ = ER_APPL
         ier_msg(1) = 'The range parameter give the width of the'
         ier_msg(2) = 'new initialization window, it is needed'
         RETURN
      ENDIF
      IF(owerte(O_MIN) == -9999.00) THEN ! min was provided
         IF(owerte(O_VALUE) /= -9999.00 ) THEN !value: was provided
            IF(owerte(O_MIN) > owerte(O_VALUE) ) THEN
               ier_num = -30
               ier_typ = ER_APPL
               ier_msg(1) = 'Minimum is larger than value'
               RETURN
            ENDIF
         ENDIF
         IF(owerte(O_MAX) /= -9999.00 ) THEN !max: was provided
            IF(owerte(O_MIN) > owerte(O_MAX) ) THEN
               ier_num = -30
               ier_typ = ER_APPL
               ier_msg(1) = 'Minimum is larger than maximum'
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF(owerte(O_MAX) == -9999.00) THEN ! max was provided
         IF(owerte(O_VALUE) /= -9999.00 ) THEN !value: was provided
            IF(owerte(O_MAX) < owerte(O_VALUE) ) THEN
               ier_num = -30
               ier_typ = ER_APPL
               ier_msg(1) = 'Maximum is smaller than value'
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ldismiss = .FALSE.
      IF(opara(O_DISMISS) == '0.000000') THEN
         d_lb = 0
         d_ub = 0
         ldismiss = .FALSE.
      ELSEIF(opara(O_DISMISS) == 'none') THEN
         d_lb = 0
         d_ub = 0
         ldismiss = .FALSE.
      ELSEIF(opara(O_DISMISS) == 'all') THEN
         d_lb = 1
         d_ub = pop_n
         ldismiss = .TRUE.
      ELSEIF(opara(O_DISMISS) == 'best') THEN
         d_lb = 2
         d_ub = pop_n
         ldismiss = .TRUE.
      ELSE
         ianz = NOPTIONAL
         CALL ber_params (ianz, opara, lopara, owerte, NOPTIONAL)
         IF(ier_num==0) THEN
            d_lb = pop_n + 1 - NINT(owerte(O_DISMISS))
            d_ub = pop_n
            IF ( 0<d_lb .and. d_lb<=d_ub .and. d_ub<=pop_n) THEN
               ldismiss = .TRUE.
            ELSE
               ier_msg(1) = 'Dismiss range must be zero to pop_n[1]'
               ier_num = -6
               ier_typ = ER_COMM
               RETURN
            ENDIF
         ELSE
            ier_msg(1) = 'Error calculating dismiss: value'
            RETURN
         ENDIF
      ENDIF
      ident: DO k=1,pop_dimx
         IF(cpara(1)==pop_name(k)) THEN
            WRITE(cpara(1),'(I4)') k
            EXIT ident
         ENDIF
      ENDDO ident
      ianz = 1
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF (ier_num.eq.0) THEN
         owerte(O_RANGE) = ABS(owerte(O_RANGE))      ! Ensure positive value for the range
         lb = nint(werte(1))
         IF ( 0<lb .and. lb<=pop_dimx) THEN
            set_value = child(lb,pop_best)           ! Default for value
            IF(owerte(O_VALUE) /= -9999.0) THEN      ! 'value' was provided
               set_value = owerte(O_VALUE)
            ENDIF
            set_xmin = set_value - 3.*owerte(O_RANGE)
            IF(owerte(O_MIN  ) /= -9999.0) THEN      ! 'min  ' was provided
               set_xmin = owerte(O_MIN)
            ENDIF
            set_xmax = set_value + 3.*owerte(O_RANGE)
            IF(owerte(O_MAX  ) /= -9999.0) THEN      ! 'min  ' was provided
               set_xmax = owerte(O_MAX)
            ENDIF
            pop_xmin(lb) = set_xmin
            pop_xmax(lb) = set_xmax
            pop_smin(lb) = MAX(pop_xmin(lb), set_value - 1.*REAL(owerte(O_RANGE)))
            pop_smax(lb) = MIN(pop_xmax(lb), set_value + 1.*REAL(owerte(O_RANGE)))
            pop_refine(lb) = .TRUE.
            CALL init_x(lb,lb)
            IF(ldismiss) CALL do_dismiss(d_lb,d_ub)
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = ' The parameter name or number is outside the '
            ier_msg(2) = ' range defined by ''newpara'' or pop_dimx'
            RETURN
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         RETURN
      ENDIF
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = ' The release command takes one parameter name'
      ier_msg(2) = ' and one optional parameter range:'
      ier_msg(3) = ' '
   ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE do_release
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE diffev_release_mod
