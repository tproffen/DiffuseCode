MODULE learn_mod
!+
!     Learning mode flag
!-
USE precision_mod
!
IMPLICIT NONE
PUBLIC
SAVE
!
CHARACTER(LEN=PREC_STRING) :: fname
LOGICAL ::  llearn = .false.
!
END MODULE learn_mod
