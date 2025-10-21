MODULE gauss_lorentz_pseudo_mod
!
! Functions and data for a Gaussian, Lorentzian and 
! Pseudovoigt function normalized to a width FWHM 1.0
!
USE precision_mod
use wink_mod
!
PRIVATE
PUBLIC  GLP_MAX, GLP_NPT, glp_gauss_tbl, glp_loren_tbl
PUBLIC  glp_setup, glp_gauss, glp_loren, glp_pseud, glp_pseud_indx
public  lognormal
!
INTEGER           , PARAMETER :: GLP_MAX = 40000
REAL(KIND=PREC_DP), PARAMETER :: glp_step = 0.0005d0
INTEGER           , PARAMETER :: GLP_NPT = INT(1./glp_step)
!
REAL(KIND=PREC_DP), DIMENSION(0:GLP_MAX) :: glp_gauss_tbl ! Gauss lookup table
REAL(KIND=PREC_DP), DIMENSION(0:GLP_MAX) :: glp_loren_tbl ! Lorentzian lookup table
LOGICAL            :: glp_init = .TRUE.
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE glp_setup
!
! Set values into Gaussuian and Lorentzian Lookup tables, performed only once
!
! Contants : four_ln2  = 4.*ln(2.0)          = 2.772588722239781237669D
!            sq4ln2_pi = SQRT(4.*ln(2.0)/PI) = 0.9394372786996513337723D0
!            two_pi    = 2.0/PI              = 
!
! The functions are calculated as
! GAUSS = SQRT(4.*ln(2.0)/PI) / FWHM * exp(-4.*ln(2.0)*(x/FWHM)**2)
! LOREN = 2/PI * FWHM / ( FWHM**2 + 4x**2)
! with FWHM = 1.0
!
USE wink_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), PARAMETER :: four_ln2  = 2.772588722239781237669D0
REAL(KIND=PREC_DP), PARAMETER :: sq4ln2_pi = 0.9394372786996513337723D0
REAL(KIND=PREC_DP), PARAMETER :: two_pi    = 2.0D0/PI
!
INTEGER :: i
REAL(KIND=PREC_DP) :: pref_g
REAL(KIND=PREC_DP) :: pref_l
!
pref_g = sq4ln2_pi / 1.0D0
pref_l = two_pi * 1.0D0
!
IF(glp_init) THEN
  DO i=0,GLP_MAX
     glp_gauss_tbl(i) = pref_g * exp(-four_ln2 * (REAL(i,KIND=PREC_DP)*glp_step)**2)
  ENDDO
  DO i=0,GLP_MAX
     glp_loren_tbl(i) = pref_l / ( 1.0D0 + 4.D0*(REAL(i,KIND=PREC_DP)*glp_step)**2)
  ENDDO
ENDIF
!
glp_init = .FALSE.
!
END SUBROUTINE glp_setup
!
!*******************************************************************************
!
REAL(KIND=PREC_DP) FUNCTION glp_gauss(x, fwhm)
!
!  Interface to lookup table for Gaussian from x and FWHM
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), INTENT(IN) :: x     ! Position, must be >= 0.0D0
REAL(KIND=PREC_DP), INTENT(IN) :: fwhm  ! FWHM,     must be >= 0.0D0
!
INTEGER :: i
!
glp_gauss = 0.0D0
!
i =     INT(x/fwhm)*GLP_NPT
glp_gauss = glp_gauss_tbl(i)
!
END FUNCTION glp_gauss
!
!*******************************************************************************
!
REAL(KIND=PREC_DP) FUNCTION glp_loren(x, fwhm)
!
!  Interface to lookup table for Lorentzian from x and FWHM
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), INTENT(IN) :: x     ! Position, must be >= 0.0D0
REAL(KIND=PREC_DP), INTENT(IN) :: fwhm  ! FWHM,     must be >= 0.0D0
!
INTEGER :: i
!
glp_loren = 0.0D0
!
i =     INT(x/fwhm)*GLP_NPT 
glp_loren = glp_loren_tbl(i)
!
END FUNCTION glp_loren
!
!*******************************************************************************
!
REAL(KIND=PREC_DP) FUNCTION glp_pseud(x, eta, fwhm)
!
!  Interface to lookup table for Pseudo-Voigt from x, eta and FWHM
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), INTENT(IN) :: x     ! Position, must be >= 0.0D0
REAL(KIND=PREC_DP), INTENT(IN) :: eta   ! Mixing parameter
REAL(KIND=PREC_DP), INTENT(IN) :: fwhm  ! FWHM,     must be >= 0.0D0
!
INTEGER :: i
!
i = MIN(INT(x/fwhm*GLP_NPT), GLP_MAX)
glp_pseud = eta/fwhm*glp_loren_tbl(i  ) + (1.0D0-eta)/fwhm*glp_gauss_tbl(i  )
!
END FUNCTION glp_pseud
!
!*******************************************************************************
!
REAL(KIND=PREC_DP) FUNCTION glp_pseud_indx(i, eta, fwhm)
!
!  Interface to lookup table for Pseudo-Voigt from index i, eta and FWHM
!  i = MIN(INT(x/fwhm*GLP_NPT), GLP_MAX) calculated externally
!
IMPLICIT NONE
!
INTEGER           , INTENT(IN) :: i     ! Position, must be >= 0.0 is INT(x/fwhm*GLP_NPT)
REAL(KIND=PREC_DP), INTENT(IN) :: eta   ! Mixing parameter
REAL(KIND=PREC_DP), INTENT(IN) :: fwhm  ! FWHM,     must be >= 0.0D0
!
!
glp_pseud_indx =(eta*glp_loren_tbl(i) + (1.0D0-eta)*glp_gauss_tbl(i))/fwhm
!
END FUNCTION glp_pseud_indx
!
!*******************************************************************************
!
real(kind=PREC_DP) function lognormal(x, ln_mean, ln_sigma)
!-
!  Calculate a lognormal function value
!  ln_mean is the ln(mean) of the Gaussian distribution
!  ln_sigma is the ln(sigma) of the Gaussian distribution
!+
!
implicit none
!
real(kind=PREC_DP), intent(in) :: x
real(kind=PREC_DP), intent(in) :: ln_mean
real(kind=PREC_DP), intent(in) :: ln_sigma
!
lognormal = 1./(x*ln_sigma*sq_zpi)*exp(-((log(x)-ln_mean)**2)/(2.0_PREC_DP*ln_sigma**2))
!
end function lognormal
!
!*******************************************************************************
!
END MODULE gauss_lorentz_pseudo_mod
