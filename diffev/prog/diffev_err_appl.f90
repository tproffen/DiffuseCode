!*****7****************************************************************
!
SUBROUTINE diffev_errlist_appl
!-
!     Displays error Messages for the error type APPLication
!+
      USE errlist_mod
      IMPLICIT      none
!
!
      INTEGER       iu,io
      PARAMETER    (IU=-37,IO=0)
!
      CHARACTER(LEN=45) ERROR(IU:IO)
!
      DATA ERROR ( IU:-21) /                                            &
     &  'Error reading partial Rvlaues',                                & !-37  ! diffev
     &  'Generation number exceeds current value',                      & !-36  ! diffev
     &  'GENERATION file does not exist',                               & !-35  ! diffev
     &  'THe ""range"" parameter is required',                          & !-34  ! diffev
     &  'Less children*NINDIV than compute nodes*CPUs ',                & !-33  ! diffev
     &  'Less children than compute nodes             ',                & !-32  ! diffev
     &  'Parameter name is too long',                                   & !-31  ! diffev
     &  'Conflicting lower/upper values for the limits',                & !-30  ! diffev
     &  'Parameter name already in use',                                & !-29  ! diffev
     &  'Fixed parameter has non-zero distribution',                    & !-28  ! diffev
     &  'Silent mode within discus_suite only',                         & !-27  ! diffev
     &  'MPI TERMINATED due to slave error ',                           & !-26  ! diffev
     &  'Parameter could not be initialized',                           & !-25  ! diffev
     &  'Diffev not compiled/started with MPI',                         & !-24  ! diffev
     &  'MPI returned slave number, MPI is not active',                 & !-23  ! diffev
     &  'Error initializing MPI system',                                & !-22  ! diffev
     &  'No MPI run standalone macro via system '                       & !-21  ! diffev
     &     /
!
      DATA ERROR (-20:  0) /                                            &
     &  'Generation number is zero    ',                                & !-20  ! diffev
     &  'Error reading GENERATION file',                                & !-19  ! diffev
     &  'Parameter value outside hard limits',                          & !-18  ! diffev
     &  'Conflicting parameter size in GENERATION file',                & !-17  ! diffev
     &  'Conflicting population size in GENERATION',                    & !-16  ! diffev
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
      END SUBROUTINE diffev_errlist_appl
!*****7****************************************************************
