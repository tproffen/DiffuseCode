MODULE precision_command_mod
!
CHARACTER(LEN=2), PARAMETER :: C_WIDTH = '20'  ! '15'
CHARACTER(LEN=2), PARAMETER :: C_EXPO  = '03'  ! '02'
!                     Overall width      WW 20 15
INTEGER, PARAMETER :: PREC_WIDTH = 10*(IACHAR(C_WIDTH(1:1))-IACHAR('0')) &
                                    + (IACHAR(C_WIDTH(2:2))-IACHAR('0'))
INTEGER, PARAMETER :: PREC_EXPO  = 10*(IACHAR(C_EXPO (1:1))-IACHAR('0')) &
                                    + (IACHAR(C_EXPO (2:2))-IACHAR('0'))
!                     Digits in exponent xx  3  2
INTEGER, PARAMETER :: PREC_MANTIS = PREC_WIDTH - PREC_EXPO - 2 ! width left of E+xx    15 11
INTEGER, PARAMETER :: PREC_DIGIT  = PREC_MANTIS - 3            ! Significant digits DD 12  8
!
!  Automatically build format string '(E15.08E02)'             ! Format string EWW.DDExx
CHARACTER(LEN=24), PARAMETER :: PREC_F_REAL = &
   '(E' // C_WIDTH // '.' //                                     &
   ACHAR(PREC_DIGIT/10+IACHAR('0'))      //                      &
   ACHAR(MOD(PREC_DIGIT,10)+IACHAR('0')) // 'E' // C_EXPO // ')'
!
!  Automatically build format string '(I15)'                   ! Format string IWW
CHARACTER(LEN=24), PARAMETER :: PREC_F_INTE = '(I' // C_WIDTH // ')'
!
END MODULE precision_command_mod
