MODULE constraint
!-
!      Variables needed to define the constraints
!+
!
SAVE
PUBLIC
!
CHARACTER (LEN=1024),DIMENSION(:), ALLOCATABLE ::  constr_line
!
INTEGER                                        ::  MAX_CONSTR    ! Maximum constraint number
INTEGER             ,DIMENSION(:), ALLOCATABLE ::  constr_length
INTEGER                                        ::  constr_number
INTEGER                                        ::  constr_size_of  ! Bytes allocated for constraint
!
!      common /constr/    constr_line,constr_length,constr_number
!
END MODULE constraint
