MODULE recipro_mod
!+
!
!     This file contains variables for inverse Fourier and
!     Patterson input
!-
SAVE
!
INTEGER, PARAMETER  ::  REC_MAX_SYM  =  48
!
INTEGER                             :: rec_n_sym
REAL   , DIMENSION(4,4,REC_MAX_SYM) :: rec_sym
!
END MODULE recipro_mod
