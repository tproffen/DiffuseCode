MODULE lib_send_func
!
CONTAINS
!
! These are generic send and get routines that are used by 
! the programs in the interface to python. They allow
! python to send / get parts of the arrays i[] and r[]
! to / from (discus, diffev, kuplot)
!
! Each of '*_py.f90' can have its own specific exchange
! routines as well.
!
! Currently this file is included into '*_py' by an explicit
! INCLUDE 'send_get.f90'
! line. If I decide to drop the module around the 
! '*_py.f90' file, then it might be more elegant to
! add them to the f2py line in CMakeLists.txt. 
! Does not really make any difference
!
SUBROUTINE send_i (iin, lower, upper )
!
! The outer routine sends integer valued numbers for i[:]
!
USE setup_mod
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                         INTENT(IN) :: lower
INTEGER,                         INTENT(IN) :: upper
INTEGER, DIMENSION(lower:upper), INTENT(IN) :: iin
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL setup
ENDIF
!
IF(lower>0 .and. upper<501 .and. lower <= upper) THEN
   inpara(lower:upper) = iin(lower:upper)
ENDIF
!
END SUBROUTINE send_i
!
!
SUBROUTINE send_r (rin, lower, upper )
!
! The outer routine sends real valued numbers for r[:]
!
USE setup_mod
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                      INTENT(IN) :: lower
INTEGER,                      INTENT(IN) :: upper
REAL, DIMENSION(lower:upper), INTENT(IN) :: rin
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL setup
ENDIF
!
IF(lower>0 .and. upper<501 .and. lower <= upper) THEN
   rpara(lower:upper) = rin(lower:upper)
ENDIF
!
END SUBROUTINE send_r
!
SUBROUTINE get_i (iout, lower, upper )
!
! The outer routine gets integer valued numbers from i[:]
!
USE setup_mod
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                         INTENT(IN ) :: lower
INTEGER,                         INTENT(IN ) :: upper
INTEGER, DIMENSION(lower:upper), INTENT(OUT) :: iout
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL setup
ENDIF
!
IF(lower>0 .and. upper<501 .and. lower <= upper) THEN
   iout(lower:upper) = inpara(lower:upper) 
ENDIF
!
END SUBROUTINE get_i
!
!
SUBROUTINE get_r (rout, lower, upper )
!
! The outer routine gets real valued numbers from r[:]
!
USE setup_mod
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                      INTENT(IN ) :: lower
INTEGER,                      INTENT(IN ) :: upper
REAL, DIMENSION(lower:upper), INTENT(OUT) :: rout
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL setup
ENDIF
!
IF(lower>0 .and. upper<501 .and. lower <= upper) THEN
   rout(lower:upper) = rpara(lower:upper) 
ENDIF
!
END SUBROUTINE get_r
END MODULE lib_send_func
