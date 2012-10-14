!
!     This file contains COMMON block for the error variables
!
!*****7*****************************************************************
!
!     ER_NONE      ! no error
!     ER_COMM      ! Command language error
!     ER_FORT      ! FORTRAN Interpreter error
!     ER_IO      ! General I/O error
!     ER_MA!     ! Macro handling
!     ER_MATH      ! General mathematical error
!     ER_APPL      ! Application specific error
!     ER_RMC  ! error in RMC module
!     ER_CHEM ! error in CHEM module
!     ER_FOUR ! error in FOUR module
!     ER_PDF  ! error in PDF  module
!
!     ER_S_CONT ! Continuation upon error
!     ER_S_EXIT ! Program termination upon error
!
      INTEGER        ER_NONE
      INTEGER        ER_COMM
      INTEGER        ER_FORT
      INTEGER        ER_IO
      INTEGER        ER_MAC
      INTEGER        ER_MATH
      INTEGER        ER_APPL
      INTEGER        ER_RMC
      INTEGER        ER_MMC
      INTEGER        ER_CHEM
      INTEGER        ER_FOUR
      INTEGER        ER_PDF
!
      PARAMETER     (ER_NONE = 0 )
      PARAMETER     (ER_COMM = 1 )
      PARAMETER     (ER_FORT = 2 )
      PARAMETER     (ER_IO   = 3 )
      PARAMETER     (ER_MAC  = 4 )
      PARAMETER     (ER_MATH = 5 )
      PARAMETER     (ER_APPL = 6 )
      PARAMETER     (ER_RMC  = 7 )
      PARAMETER     (ER_CHEM = 8 )
      PARAMETER     (ER_FOUR = 9 )
      PARAMETER     (ER_MMC  = 10)
      PARAMETER     (ER_PDF  = 11)
!
      INTEGER        ER_S_CONT
      INTEGER        ER_S_EXIT
      INTEGER        ER_S_LIVE
!
      PARAMETER     (ER_S_CONT = 0)
      PARAMETER     (ER_S_EXIT = 1)
      PARAMETER     (ER_S_LIVE =-1)
!
      CHARACTER*80   ier_msg(3)
      CHARACTER*80   ier_out
      INTEGER        ier_num
      INTEGER        ier_typ
      INTEGER        ier_sta
!
      COMMON /errmsg/ ier_msg
      COMMON /errors/ ier_out,ier_num,ier_typ,ier_sta
!
