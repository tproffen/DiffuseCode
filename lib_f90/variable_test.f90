MODULE variable_test
!
USE variable_mod
!
IMPLICIT NONE
!
CONTAINS
!
   SUBROUTINE variable_exist(c_temp,l_temp, c_type, l_exist, l_type, var_no)
!
   CHARACTER (LEN=*), INTENT(IN)  :: c_temp    ! Name to be tested
   INTEGER          , INTENT(IN)  :: l_temp    ! length of input name
   INTEGER          , INTENT(IN)  :: c_type    ! Type of input variable
   LOGICAL          , INTENT(OUT) :: l_exist   ! True if exists
   LOGICAL          , INTENT(OUT) :: l_type    ! True if correct type
   INTEGER          , INTENT(OUT) :: var_no    ! variable number if exists
!
   INTEGER :: i
!
   l_exist = .false.
   l_type  = .false.
   var_no  = 0
!
   search: DO i=1,var_num
      IF( c_temp(1:l_temp) == var_name(i)) THEN
         l_exist = .true.
         var_no  = i
         IF( c_type == var_type(i)) THEN
            l_type = .true.
         ENDIF
         EXIT search
      ENDIF
   ENDDO search
!
   END SUBROUTINE variable_exist
!
END MODULE variable_test
