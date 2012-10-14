!*****7****************************************************************
!
      SUBROUTINE errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      IMPLICIT      none
!
      INCLUDE      'errlist.inc'
!
      INTEGER       iu,io
      PARAMETER    (IU=-22,IO=0)
!
      CHARACTER(LEN=41) ERROR(IU:IO)
!
      DATA ERROR (-22:-21) /                                            &
     &  'Error initializing MPI system',                                & !-22  ! diffev
     &  'No MPI included, run standalon macro via system command'       & !-21  ! diffev
     &     /
!
      DATA ERROR (-20:  0) /                                            &
     &  'Generation number is zero    ',                                & !-20  ! diffev
     &  'Error reading GENERATION file',                                & !-19  ! diffev
     &  'Parameter value outside hard limits',                          & !-18  ! diffev
     &  'Conflicting parameter size in GENERATION file',                & !-17  ! diffev
     &  'Conflicting population size in GENERATION file',               & !-16  ! diffev
     &  'Adjust size of NSTACK in compare.f90',                         & !-15  ! diffev
     &  'Refinement parameter number outside limits',                   & !-14  ! diffev
     &  'Error reading parameters from logfile',                        & !-13  ! diffev
     &  'Error reading parameters for a member',                        & !-12  ! diffev
     &  'Error reading cost function result value',                     & !-11  ! diffev
     &  'Refinement parameter number < zero',                           & !-10  ! diffev
     &  'Member number outside limits',                                 & ! -9  ! diffev
     &  'Could not fulfill constraints',                                & ! -8  ! diffev
     &  'Too many constraints conditions',                              & ! -7  ! diffev
     &  'Multiplier must be >= zero',                                   & ! -6  ! diffev
     &  'The cross over must be in interval [0,1]',                     & ! -5  ! diffev
     &  'Too many refinement parameters',                               & ! -4  ! diffev
     &  'There must be at least 4 members',                             & ! -3  ! diffev
     &  'Too many members in population',                               & ! -2  ! diffev
     &  'The generation number must be >= zero',                        & ! -1  ! diffev
     &  ' '                                                             & !  0  ! diffev
     &     /
!
      CALL disp_error ('APPL',error,iu,io)
!
      END SUBROUTINE errlist_appl
!*****7****************************************************************
