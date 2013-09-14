MODULE output_mod
!-
!     Common Block und Definitionen der Outputvariablen fuers INCLUDE
!+
SAVE
!
CHARACTER(LEN=200)      ::  outfile      = 'fcalc.dat'
INTEGER                 ::  ityp         = 0
INTEGER                 ::  extr_abs     = 1
INTEGER                 ::  extr_ord     = 2
INTEGER                 ::  extr_top     = 3
INTEGER                 ::  rho_extr_abs = 1
INTEGER                 ::  rho_extr_ord = 2
INTEGER                 ::  out_extr_abs = 1
INTEGER                 ::  out_extr_ord = 2
INTEGER                 ::  out_extr_top = 3
INTEGER, DIMENSION(3)   ::  out_inc ! (2)
REAL   , DIMENSION(3,4) ::  out_eck ! (3,4)
REAL   , DIMENSION(3,3) ::  out_vi  ! (3,3)
CHARACTER(LEN=  3)      ::  cpow_form    = 'tth'
!
END MODULE output_mod
