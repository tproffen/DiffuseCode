MODULE macro_mod
!+
!     include file for file related variables, all macro stuff is "macro_internal.f90"
!-
USE precision_mod
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, PARAMETER :: MAC_MAX_LEVEL = 10
   INTEGER, PARAMETER :: MAC_MAX_PARA  = 20
   INTEGER, PARAMETER :: MAC_MAX_IO    = 20
   INTEGER, PARAMETER :: MAC_MAX_FORM  = 50
!
   CHARACTER(LEN=PREC_STRING), DIMENSION(1:MAC_MAX_IO   ) :: io_file       ! (MAC_MAX_IO)
   CHARACTER(LEN=20  ), DIMENSION(1:MAC_MAX_FORM ) :: io_out_format ! (MAC_MAX_FORM)
   LOGICAL            , DIMENSION(1:MAC_MAX_IO   ) :: io_open       ! (MAC_MAX_IO)
   INTEGER            , DIMENSION(1:MAC_MAX_IO   ) :: io_unit       ! (MAC_MAX_IO)
   LOGICAL            , DIMENSION(1:MAC_MAX_IO   ) :: io_eof        ! (MAC_MAX_IO)
   INTEGER            , DIMENSION(1:MAC_MAX_IO, 2) :: io_get_sub    ! (MAC_MAX_IO,2)
!
END MODULE macro_mod
