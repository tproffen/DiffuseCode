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
USE errlist_mod
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
IF(prop_user_no>0) THEN
!
   test = .TRUE.
   check: DO i=1,prop_user_no
      IF(prop_user(i)%act/=0) THEN                               ! NOT ignore
         IF(                    cr_iscat(iatom,1) == prop_user(i)%at_type  & ! same atom type
                           .OR. -1 == prop_user(i)%at_type    ) THEN       ! Any atom type allowed
            ino  = prop_user(i)%conn_no
            c_name   = prop_user(i)%conn_name
            c_name_l = LEN_TRIM(c_name)
            ll       = c_name_l
            CALL get_connectivity_identity(cr_iscat(iatom,1), ino, c_name, c_name_l)
            IF(ier_num == 0) THEN       ! Atom has a requested connectivity
               test = test .AND. prop_user(i)%conn_name(1:ll)==c_name(1:c_name_l)
               CALL get_connectivity_list (iatom, cr_iscat(iatom,1), ino, c_list, c_offs, natoms )
               test = test .AND. prop_user(i)%n_min<=natoms     &         ! Correct number of neighbors
                        .AND.                     natoms<=prop_user(i)%n_max
               test = test .AND. .NOT. (prop_user(i)%e_min<=natoms       &       ! Correct number of neighbors
                                     .AND.               natoms<=prop_user(i)%e_max)
            ELSE                        ! Atom does not have requested connectivity
               test = .FALSE.
               ier_num = 0
               ier_typ = 0
               ier_msg(:) = ' '
            ENDIF
            IF(prop_user(i)%act==-1) test = .NOT.test               ! Absent invert the test
         ENDIF
      ENDIF
   ENDDO check
!
ELSE
   test = .TRUE.
ENDIF
!
check_user_property = test
!
END FUNCTION check_user_property
!
!*******************************************************************************
!
END MODULE check_user_prop_mod
