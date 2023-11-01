MODULE constraint
!-
!      Variables needed to define the constraints
!+
!
USE precision_mod
!
SAVE
PUBLIC
!
CHARACTER (LEN=PREC_STRING),DIMENSION(:), ALLOCATABLE ::  constr_line
!
INTEGER                                        ::  MAX_CONSTR    ! Maximum constraint number
INTEGER             ,DIMENSION(:), ALLOCATABLE ::  constr_length
INTEGER                                        ::  constr_number
!
END MODULE constraint
