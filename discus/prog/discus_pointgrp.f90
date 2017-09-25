MODULE point_grp
!
!  Creates the symmetry equivalent faces hkl for an input face
!
PRIVATE
PUBLIC point_init, point_set, point_test
!
INTERFACE point_init
   MODULE PROCEDURE point_init_real, point_init_int
END INTERFACE point_init
INTERFACE point_set
   MODULE PROCEDURE point_set_real, point_set_int
END INTERFACE point_set
INTERFACE point_test
   MODULE PROCEDURE point_test_real, point_test_int
END INTERFACE point_test
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE point_init_real(hkl, point_hkl, point_n)
!
!  Initialize. Delete old arrays make a default size and set input value
!
USE wyckoff_mod
!
IMPLICIT NONE
!
REAL, DIMENSION(4)                , INTENT(IN)    :: hkl
REAL, DIMENSION(:,:), ALLOCATABLE , INTENT(INOUT) :: point_hkl
INTEGER                           , INTENT(OUT)   :: point_n
!
IF(ALLOCATED(point_hkl)) DEALLOCATE(point_hkl)
ALLOCATE(point_hkl(1:3, 1:96))
point_hkl = 0.0
point_n   = 1
point_hkl(:,point_n) = NINT(hkl(1:3))
!
END SUBROUTINE point_init_real
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE point_init_int(hkl, point_hkl, point_n)
!
!  Initialize. Delete old arrays make a default size and set input value
!
USE wyckoff_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(4)                , INTENT(IN)    :: hkl
INTEGER, DIMENSION(:,:), ALLOCATABLE , INTENT(INOUT) :: point_hkl
INTEGER                              , INTENT(OUT)   :: point_n
!
IF(ALLOCATED(point_hkl)) DEALLOCATE(point_hkl)
ALLOCATE(point_hkl(1:3, 1:96))
point_hkl = 0
point_n   = 1
point_hkl(:,point_n) = hkl(1:3)
!
END SUBROUTINE point_init_int
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE point_set_real(hkl, point_hkl, point_n)
!
USE crystal_mod
USE wyckoff_mod
!
IMPLICIT NONE
!
REAL, DIMENSION(4)   , INTENT(IN)  :: hkl
REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: point_hkl
INTEGER              , INTENT(INOUT) :: point_n
!
REAL, PARAMETER :: EPS = 0.000001
!
INTEGER            :: i,j
LOGICAL            :: l_new
REAL, DIMENSION(4) :: hklw
REAL, DIMENSION(:,:), ALLOCATABLE :: temp_hkl
!
matrix_set: DO j=1, spc_n
   hklw = MATMUL(hkl,spc_mat(:,:,j))
   hklw(4) = 0
   l_new = .TRUE.
   search: DO i=1, point_n
      IF(ABS(point_hkl(1,i)-hklw(1)).lt.EPS .AND.   &
         ABS(point_hkl(2,i)-hklw(2)).lt.EPS .AND.   &
         ABS(point_hkl(3,i)-hklw(3)).lt.EPS ) THEN
         l_new = .FALSE.
         EXIT search
      ENDIF
   ENDDO search
   IF(l_new) THEN
      IF(UBOUND(point_hkl,2)==point_n) THEN  ! INCREASE SIZE
         temp_hkl = point_hkl                ! Automatically allocate new array
         DEALLOCATE(point_hkl)
         ALLOCATE  (point_hkl(1:3, 1:point_n + 48))
         point_hkl(1:3,1:UBOUND(temp_hkl,2)) = temp_hkl(:,:)
         DEALLOCATE(temp_hkl)
      ENDIF
      point_n = point_n + 1
      point_hkl(:,point_n) = hklw(1:3)
   ENDIF
ENDDO matrix_set
!
END SUBROUTINE point_set_real
!
!
SUBROUTINE point_set_int(hkl, point_hkl, point_n)
!
USE crystal_mod
USE wyckoff_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(4)               , INTENT(IN)  :: hkl
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: point_hkl
INTEGER                             , INTENT(INOUT) :: point_n
!
REAL, PARAMETER :: EPS = 0.000001
!
INTEGER            :: i,j
LOGICAL            :: l_new
INTEGER, DIMENSION(4) :: hklw
INTEGER, DIMENSION(:,:), ALLOCATABLE :: temp_hkl
!
matrix_set: DO j=1, spc_n
   hklw = INT(MATMUL(hkl,spc_mat(:,:,j)))
   hklw(4) = 0
   l_new = .TRUE.
   search: DO i=1, point_n
      IF(ABS(point_hkl(1,i)-hklw(1)).lt.EPS .AND.   &
         ABS(point_hkl(2,i)-hklw(2)).lt.EPS .AND.   &
         ABS(point_hkl(3,i)-hklw(3)).lt.EPS ) THEN
         l_new = .FALSE.
         EXIT search
      ENDIF
   ENDDO search
   IF(l_new) THEN
      IF(UBOUND(point_hkl,2)==point_n) THEN  ! INCREASE SIZE
         temp_hkl = point_hkl                ! Automatically allocate new array
         DEALLOCATE(point_hkl)
         ALLOCATE  (point_hkl(1:3, 1:point_n + 48))
         point_hkl(1:3,1:UBOUND(temp_hkl,2)) = temp_hkl(:,:)
         DEALLOCATE(temp_hkl)
      ENDIF
      point_n = point_n + 1
      point_hkl(:,point_n) = hklw(1:3)
   ENDIF
ENDDO matrix_set
!
END SUBROUTINE point_set_int
!
LOGICAL FUNCTION point_test_real(hkl, point_hkl, point_n, l_form)
!
USE crystal_mod
USE wyckoff_mod
!
IMPLICIT NONE
!
REAL, DIMENSION(4)   , INTENT(IN)  :: hkl
REAL, DIMENSION(:,:) , INTENT(OUT) :: point_hkl
INTEGER              , INTENT(OUT) :: point_n
LOGICAL              , INTENT(IN)  :: l_form
!
REAL, PARAMETER :: EPS = 0.000001
!
INTEGER :: i
INTEGER :: last
!
last = 1
IF(l_form) last = point_n
!
point_test_real = .FALSE.
main:DO i=1,last
   IF( ABS(hkl(1) - point_hkl(1,i)) < EPS .AND. &
       ABS(hkl(2) - point_hkl(2,i)) < EPS .AND. &
       ABS(hkl(3) - point_hkl(3,i)) < EPS ) THEN
      point_test_real = .TRUE.
      EXIT main
   ENDIF
ENDDO main
!
END FUNCTION point_test_real
!
!
LOGICAL FUNCTION point_test_int(hkl, point_hkl, point_n, l_form)
!
USE crystal_mod
USE wyckoff_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(4)   , INTENT(IN)  :: hkl
INTEGER, DIMENSION(:,:) , INTENT(OUT) :: point_hkl
INTEGER                 , INTENT(OUT) :: point_n
LOGICAL                 , INTENT(IN)  :: l_form
!
REAL, PARAMETER :: EPS = 0.000001
!
INTEGER :: i
INTEGER :: last
!
last = 1
IF(l_form) last = point_n
!
point_test_int = .FALSE.
main:DO i=1,last
   IF( ABS(hkl(1) - point_hkl(1,i)) < EPS .AND. &
       ABS(hkl(2) - point_hkl(2,i)) < EPS .AND. &
       ABS(hkl(3) - point_hkl(3,i)) < EPS ) THEN
      point_test_int = .TRUE.
      EXIT main
   ENDIF
ENDDO main
!
END FUNCTION point_test_int
!
!
END MODULE point_grp
