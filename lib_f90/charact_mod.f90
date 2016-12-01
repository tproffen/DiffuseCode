MODULE charact_mod
!+
!     Parameter definitions for constant strings of ASCII equivalents
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   CHARACTER(LEN=1), PARAMETER :: TAB = ACHAR(9)
!
   INTEGER, PARAMETER ::  a      = IACHAR('a')
   INTEGER, PARAMETER ::  z      = IACHAR('z')
   INTEGER, PARAMETER ::  aa     = IACHAR('A')
   INTEGER, PARAMETER ::  zz     = IACHAR('Z')
!
   INTEGER, PARAMETER ::  zero   = IACHAR('0')
   INTEGER, PARAMETER ::  nine   = IACHAR('9')
!
   INTEGER, PARAMETER ::  period = IACHAR('.')
!
   INTEGER, PARAMETER ::  u      = IACHAR('_')
   INTEGER, PARAMETER ::  blank1 = IACHAR(' ')
!
END MODULE charact_mod
