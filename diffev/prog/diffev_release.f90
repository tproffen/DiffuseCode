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
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(INOUT) :: zeile
INTEGER          , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: maxw = 2
!
CHARACTER (LEN=1024), DIMENSION(MAXW) :: cpara   = ' '
INTEGER             , DIMENSION(MAXW) :: lpara = 0
REAL                , DIMENSION(MAXW) :: werte = 0.0
!
INTEGER                               :: ianz
INTEGER                               :: k
INTEGER                               :: lb
LOGICAL                               :: lexist
LOGICAL                               :: l_init_x
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
!
DATA oname  / 'range'  /
DATA loname /  5       /
opara  =  (/ '-1.00000' /)   ! Always provide fresh default values
lopara =  (/  8         /)
owerte =  (/  -1.0      /)
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
                     oname, loname, opara, lopara, owerte)
   IF(ianz == 1) THEN
      k = 0
      ident: DO k=1,pop_dimx
         IF(cpara(1)==pop_name(k)) THEN
            WRITE(cpara(1),'(I4)') k
            EXIT ident
         ENDIF
      ENDDO ident
      ianz = 1
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF (ier_num.eq.0) THEN
         lb = nint(werte(1))
         IF ( 0<lb .and. lb<=pop_dimx) THEN
            IF(owerte(1) > 0.0) THEN   ! range: was provided
               pop_xmin(lb) = MIN(pop_xmin(lb), child(lb,pop_best) - 3.*owerte(1))
               pop_xmax(lb) = MAX(pop_xmax(lb), child(lb,pop_best) + 3.*owerte(1))
               pop_smin(lb) = MIN(pop_smin(lb), child(lb,pop_best) - 1.*owerte(1))
               pop_smax(lb) = MAX(pop_smax(lb), child(lb,pop_best) + 1.*owerte(1))
            ENDIF
            pop_refine(lb) = .TRUE.
            CALL init_x(lb,lb)
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
!
END SUBROUTINE do_release
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE diffev_release_mod
