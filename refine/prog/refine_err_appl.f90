!*****7****************************************************************
!
SUBROUTINE refine_errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      USE errlist_mod
USE lib_errlist_func
      IMPLICIT      none
!
!
      INTEGER       iu,io
      PARAMETER    (IU=-11,IO=0)
!
      CHARACTER(LEN=45) ERROR(IU:IO)
!
!
      DATA ERROR ( IU:  0) /                                       &
     &  'Parameter(s) outside range/constrain        ',            & ! -11 ! refine
     &  'FWHM parameter give negative FWHM           ',            & ! -10 ! refine
     &  'Calculation of derivatives failed           ',            & ! -9  ! refine
     &  'points is limited to 3 or 5                 ',            & ! -8  ! refine
     &  'Sigma in data set are zero                  ',            & ! -7  ! refine
     &  'An error occurred in the user macro         ',            & ! -6  ! refine
     &  'Data set must be loaded prior to sigma      ',            & ! -5  ! refine
     &  'Dimensions of calc. and observed data differ',            & ! -4  ! refine
     &  'Data not present within KUPLOT',                          & ! -3  ! refine
     &  'Macro file does not exist',                               & ! -2  ! refine
     &  'Data dimensions have not yet been defined',               & ! -1  ! refine
     &  ' '                                                        & !  0  ! refine
     &     /
!
      CALL disp_error ('APPL',error,iu,io)
!
      END SUBROUTINE refine_errlist_appl
!*****7****************************************************************
