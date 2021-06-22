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
INTEGER, DIMENSION(3)   ::  out_center   = 0  ! Center of map in pixels
INTEGER, DIMENSION(3)   ::  out_pixel    = 0  ! Full size of map in pixels
INTEGER                 ::  out_lrange   = 0   ! 0=No limit; 1=center; 2=quad
INTEGER                 ::  out_lcenter  = 0   ! 0=middle; 1=user
INTEGER                 ::  out_lpixel   = 0   ! 0=all pixels, 1 = user values
CHARACTER(LEN=3)        ::  out_quad = '   '
INTEGER, DIMENSION(3)   ::  out_inc ! (2)
REAL   , DIMENSION(3,4) ::  out_eck ! (3,4)
REAL   , DIMENSION(3,3) ::  out_vi  ! (3,3)
CHARACTER(LEN=  3)      ::  cpow_form    = 'tth'
LOGICAL                 ::  out_user_limits = .false.
REAL   , DIMENSION(3)   ::  out_user_values = (/1.0, 10.0, 0.01/)
integer                 ::  out_mode = 0       ! Output mode if KUPLOT 0;1;2 == new, old, add
!
END MODULE output_mod
