!
!*******************************************************************************
!
SUBROUTINE lib_get_var_type(line,length, var_is_type)
!
! Returns the variable type : INTEGER, REAL, CHARACTER, and Scalar versus field
!
USE constants_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)     , INTENT(IN)  :: line
INTEGER              , INTENT(IN)  :: length
INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
INTEGER, PARAMETER :: MAXPAR =  7
CHARACTER(LEN=16), DIMENSION(MAXPAR) :: libf90_names
INTEGER          , DIMENSION(MAXPAR) :: libf90_type
INTEGER          , DIMENSION(MAXPAR) :: libf90_dim
LOGICAL          , DIMENSION(MAXPAR) :: libf90_ro 
!
INTEGER :: i
DATA libf90_names  &
    /'F_DERIV ', 'ref_para', 'F_PARA  ', 'seed    ', 'res     ', &
     'r       ', 'i       '                                      &
    /
DATA libf90_type &
    /  IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL , &
       IS_REAL ,   IS_INTE                                       &
    /
DATA libf90_dim  &
    /  IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC                                        &
    /
DATA libf90_ro  &
    /  .FALSE. ,   .FALSE. ,   .FALSE. ,   .FALSE. ,   .TRUE.  , &
       .FALSE. ,   .FALSE.                                       &
    /
!
main: DO i=1, MAXPAR
   IF(line(1:length) == libf90_names(i)(1:LEN_TRIM(libf90_names(i)))) THEN
      var_is_type(1) = libf90_type(i)
      var_is_type(2) = libf90_dim (i)
      IF(libf90_ro(i)) THEN
         var_is_type(3) = IS_READ
      ELSE
         var_is_type(3) = IS_WRITE
      ENDIF
      RETURN
   ENDIF
ENDDO main
!
! Test user defined variables
!
vari: DO i=1,var_num
   IF(line(1:length) == var_name(i)(1:var_l(i))) THEN
      var_is_type(1) = var_type(i)
      IF(var_entry(i) == 0) THEN
         var_is_type(2) = IS_SCAL
      ELSEIF(var_field(var_entry(i))%var_shape(2)>1) THEN
         var_is_type(2) = IS_ARR 
      ELSEIF(var_field(var_entry(i))%var_shape(2)==1) THEN
         var_is_type(2) = IS_VEC
      ENDIF
      var_is_type(3) = IS_WRITE
      RETURN
   ENDIF
ENDDO vari
!
END SUBROUTINE lib_get_var_type
