MODULE macro_mod
!+
!     include file for macro and file related variables
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, PARAMETER :: MAC_MAX_LEVEL = 10
   INTEGER, PARAMETER :: MAC_MAX_PARA  = 20
   INTEGER, PARAMETER :: MAC_MAX_IO    = 10
   INTEGER, PARAMETER :: MAC_MAX_FORM  = 50
!
   CHARACTER(LEN=1024), DIMENSION(1:MAC_MAX_LEVEL) :: mac_name      ! (MAC_MAX_LEVEL)
   CHARACTER(LEN=256 ), DIMENSION(0:MAC_MAX_PARA,1:MAC_MAX_LEVEL) :: mac_para ! (0:MAC_MAX_PARA,MAC_MAX_LEVEL)
   INTEGER            , DIMENSION(0:MAC_MAX_PARA,1:MAC_MAX_LEVEL) :: mac_leng ! (0:MAC_MAX_PARA,MAC_MAX_LEVEL)
   INTEGER                                         :: mac_level
   INTEGER            , DIMENSION(1:MAC_MAX_LEVEL) :: mac_line      ! (MAC_MAX_LEVEL)
   INTEGER            , DIMENSION(1:MAC_MAX_LEVEL) :: mac_n_par     ! (MAC_MAX_LEVEL)
   LOGICAL                                         :: lmakro
   LOGICAL            , DIMENSION(1:MAC_MAX_LEVEL) :: lmacro_dbg    ! (MAC_MAX_LEVEL)
!
   CHARACTER(LEN=1024), DIMENSION(1:MAC_MAX_IO   ) :: io_file       ! (MAC_MAX_IO)
   CHARACTER(LEN=20  ), DIMENSION(1:MAC_MAX_FORM ) :: io_out_format ! (MAC_MAX_FORM)
   LOGICAL            , DIMENSION(1:MAC_MAX_IO   ) :: io_open       ! (MAC_MAX_IO)
   INTEGER            , DIMENSION(1:MAC_MAX_IO   ) :: io_unit       ! (MAC_MAX_IO)
   LOGICAL            , DIMENSION(1:MAC_MAX_IO   ) :: io_eof        ! (MAC_MAX_IO)
   INTEGER            , DIMENSION(1:MAC_MAX_IO, 2) :: io_get_sub    ! (MAC_MAX_IO,2)
!
END MODULE macro_mod
