MODULE suite_python_mod
!-
!     Variables for communication with Python
!+
USE precision_mod
!
INTEGER :: xpts = 0
INTEGER :: ypts = 0

REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xdata
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ydata

END MODULE suite_python_mod

