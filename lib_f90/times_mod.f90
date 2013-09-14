MODULE times_mod
!+
!     This file contains the time related variables
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, DIMENSION(3) :: int_time
   INTEGER, DIMENSION(3) :: int_date
   INTEGER               :: millisec
   INTEGER               :: midnight
!
   CHARACTER (LEN=24)    :: f_modt
   CHARACTER (LEN=24)    :: f_date
!
!
END MODULE times_mod
