MODULE errlist_mod
!
!     This file contains the error variables
!
!*****7*****************************************************************
!
!     ER_NONE   ! no error
!     ER_COMM   ! Command language error
!     ER_FORT   ! FORTRAN Interpreter error
!     ER_IO     ! General I/O error
!     ER_MA!    ! Macro handling
!     ER_MATH   ! General mathematical error
!     ER_APPL   ! Application specific error
!     ER_RMC    ! error in RMC module
!     ER_CHEM   ! error in CHEM module
!     ER_FOUR   ! error in FOUR module
!     ER_MMC    ! error in MMC  module
!     ER_PDF    ! error in PDF  module
!
!     ER_S_CONT ! Continuation upon error
!     ER_S_EXIT ! Program termination upon error
!     ER_S_LIVE ! Continuation including Macros upon error
!
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, PARAMETER :: ER_NONE =  0
   INTEGER, PARAMETER :: ER_COMM =  1
   INTEGER, PARAMETER :: ER_FORT =  2
   INTEGER, PARAMETER :: ER_IO   =  3
   INTEGER, PARAMETER :: ER_MAC  =  4
   INTEGER, PARAMETER :: ER_MATH =  5
   INTEGER, PARAMETER :: ER_APPL =  6
   INTEGER, PARAMETER :: ER_RMC  =  7
   INTEGER, PARAMETER :: ER_CHEM =  8
   INTEGER, PARAMETER :: ER_FOUR =  9
   INTEGER, PARAMETER :: ER_MMC  = 10
   INTEGER, PARAMETER :: ER_PDF  = 11
!
   INTEGER, PARAMETER :: ER_S_CONT =  0
   INTEGER, PARAMETER :: ER_S_EXIT =  1
   INTEGER, PARAMETER :: ER_S_LIVE = -1
!
   CHARACTER(LEN= 80), DIMENSION(3) :: ier_msg
   CHARACTER(LEN= 80) :: ier_out
   INTEGER            :: ier_num
   INTEGER            :: ier_typ
   INTEGER            :: ier_sta
!
END MODULE errlist_mod
