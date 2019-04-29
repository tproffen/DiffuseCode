MODULE refine_fix_mod
! 
!  Routines to set various user settings
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_fix(line, length)
!
USE refine_allocate_appl
USE refine_params_mod
USE refine_set_param_mod
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
REAL               , DIMENSION(MAXW) :: werte
!
INTEGER                              :: i, j       ! Dummy loop index
INTEGER                              :: ianz
INTEGER                              :: n_params   ! Number of parameters for allocate
LOGICAL                              :: lsuccess
LOGICAL                              :: is_new     ! Parameter name is not yet in list of fixed params
REAL                                 :: p_val      ! fixed parameter value
!
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: OVALUE    = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'value' /
DATA loname /  5      /
opara  =  (/ '-1.00000'/)   ! Always provide fresh default values
lopara =  (/  8        /)
owerte =  (/  -1.0     /)
!
lsuccess = .FALSE.
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) THEN
   ier_msg(1) = 'Get params in fix'
   RETURN
ENDIF
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num/=0) THEN
   ier_msg(1) = 'Get params in fix'
RETURN
ENDIF
!
IF(ianz==1) THEN
   params: DO i=1,refine_par_n
      IF(cpara(1)==refine_params(i)) THEN
         IF(lpresent(OVALUE)) THEN        ! Optional parameter was provided
            IF(opara(OVALUE)=='current') THEN
               cpara(1) = refine_params(i)
               lpara(1) = LEN_TRIM(refine_params(i))
               ianz     = 1
            ELSE
               cpara(1) = opara(OVALUE)
               lpara(1) = lopara(OVALUE)
               ianz     = 1
            ENDIF
         ELSE
            cpara(1) = refine_params(i)
            lpara(1) = LEN_TRIM(refine_params(i))
            ianz     = 1
         ENDIF
         CALL ber_params(ianz, cpara, lpara, werte, MAXW)
         p_val = werte(1)
         CALL refine_set_param(refine_par_n, refine_params(i), i, p_val )  ! Set modified value
         lsuccess = .TRUE.
!
   is_new = .TRUE.
   old_f: DO j=1, refine_fix_n
      IF(refine_params(i) == refine_fixed(j)) THEN        ! Found old parameter name
         is_new = .FALSE.
         EXIT old_f
      ENDIF
   ENDDO old_f
!
   IF(is_new) THEN                              ! New parameter, add to list
      IF(refine_fix_n==REF_MAXPARAM_FIX) THEN
         n_params = REF_MAXPARAM_FIX + 10
         CALL ALLOC_params_fix(n_params)
      ENDIF
      refine_fix_n = refine_fix_n + 1
      refine_fixed(refine_fix_n)   = refine_params(i)
   ENDIF
!
         DO j=i+1, refine_par_n
            refine_params(j-1) = refine_params(j)
            refine_p(j-1)      = refine_p(j)
            refine_range(j-1,1)= refine_range(j,1)
            refine_range(j-1,2)= refine_range(j,2)
         ENDDO
         EXIT params
      ENDIF
   ENDDO params   
   refine_par_n = refine_par_n - 1
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
END SUBROUTINE refine_fix
!
!*******************************************************************************
!
END MODULE refine_fix_mod
