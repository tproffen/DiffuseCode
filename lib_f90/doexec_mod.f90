MODULE doexec_mod
!+
!     Variablves used during do and if execution
!-
   USE doloop_mod
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   CHARACTER (LEN=1024), DIMENSION(0:MAXCOM,0:MAXLEV) :: do_comm  !(0:MAXCOM,0:MAXLEV)

   INTEGER             , DIMENSION(0:MAXCOM,0:MAXLEV) :: do_leng  ! (0:MAXCOM,0:MAXLEV)
   INTEGER                                            :: level
   INTEGER             , DIMENSION(         0:MAXLEV) :: nlevel   ! (0:MAXLEV)
   INTEGER             , DIMENSION(         0:MAXLEV) :: ilevel   ! (0:MAXLEV)
   INTEGER             , DIMENSION(         0:MAXLEV) :: jump     ! (0:MAXLEV)
   LOGICAL             , DIMENSION(         0:MAXLEV) :: ltest    ! (0:MAXLEV)
!
!
END MODULE doexec_mod
