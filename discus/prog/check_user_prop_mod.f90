MODULE check_user_prop_mod
!
PRIVATE
PUBLIC check_user_property
!
CONTAINS
!
!*******************************************************************************
!
LOGICAL FUNCTION check_user_property(iatom)
!
USE crystal_mod
USE conn_sup_mod
USE prop_para_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: iatom
!
CHARACTER(LEN=256) :: c_name
INTEGER            :: c_name_l
INTEGER            :: ino
INTEGER            :: i, ll
INTEGER            :: natoms
LOGICAL            :: test
INTEGER, DIMENSION(:), ALLOCATABLE :: c_list
INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs
!
check_user_property = .TRUE.
!
test = .TRUE.
check: DO i=1,prop_user_no
   test = .TRUE.
   IF(prop_user(i)%act/=0) THEN                               ! NOT ignore
      test = test .AND.      cr_iscat(iatom) == prop_user(i)%at_type  & ! same atom type
                        .OR. -1 == prop_user(i)%at_type                 ! Any atom type allowed
      IF(test) THEN
         ino  = prop_user(i)%conn_no
         c_name   = prop_user(i)%conn_name
         c_name_l = LEN_TRIM(c_name)
         ll       = c_name_l
         CALL get_connectivity_identity(cr_iscat(iatom), ino, c_name, c_name_l)
         test = test .AND. prop_user(i)%conn_name(1:ll)==c_name(1:c_name_l)
         CALL get_connectivity_list (iatom, cr_iscat(iatom), ino, c_list, c_offs, natoms )
         test = test .AND. prop_user(i)%n_min<=natoms     &         ! Correct number of neighbors
                     .AND.                     natoms<=prop_user(i)%n_max
         test = test .AND. .NOT. (prop_user(i)%e_min<=natoms       &       ! Correct number of neighbors
                               .AND.               natoms<=prop_user(i)%e_max)
         IF(prop_user(i)%act==-1) test = .NOT.test               ! Absent invert the test
      ENDIF
   ENDIF
   IF(.NOT.test) EXIT check
ENDDO check
check_user_property = test
!
END FUNCTION check_user_property
!
!*******************************************************************************
!
END MODULE check_user_prop_mod
