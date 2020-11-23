MODULE lib_f90_profile
!
!  Generic profile functions
!
CONTAINS
!
!*****7************************************************************************* 
!
REAL(KIND=KIND(0.0D0)) FUNCTION pseudo_voigt(xx, P_eta, P_inte, P_pos,     &
                                                      P_fwhm, P_asym1, P_asym2, &
                                                      axis) RESULT(ww)
!+
!  Calculate a Pseudo-Voigt function at position xx
!  If axis==TRUE on 2Theta axis
!  Else on Q-axis
!-
!
USE gauss_lorentz_pseudo_mod
USE trig_degree_mod
USE wink_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), INTENT(IN)  :: xx
REAL(KIND=PREC_DP), INTENT(IN)  :: P_eta
REAL(KIND=PREC_DP), INTENT(IN)  :: P_inte
REAL(KIND=PREC_DP), INTENT(IN)  :: P_pos
REAL(KIND=PREC_DP), INTENT(IN)  :: P_fwhm
REAL(KIND=PREC_DP), INTENT(IN)  :: P_asym1
REAL(KIND=PREC_DP), INTENT(IN)  :: P_asym2
LOGICAL           , INTENT(IN)  :: axis
!
!REAL(KIND=PREC_DP) :: ww
!
REAL(KIND=PREC_DP), PARAMETER:: vln2 = 4.D0 * LOG (2.D0)
REAL(KIND=PREC_DP), PARAMETER:: gpre = 2.D0 * SQRT (LOG (2.D0) / PI)
!
REAL(KIND=PREC_DP)  :: xw             ! Deviation from position 
REAL(KIND=PREC_DP)  :: zz             ! Value for asymmetry function
REAL(KIND=PREC_DP)  :: fa             ! Value for asymmetry function
REAL(KIND=PREC_DP)  :: fb             ! Value for asymmetry function
REAL(KIND=PREC_DP)  :: asym           ! Value for asymmetry function
!
INTEGER :: i
!
xw = xx - P_pos                 ! Deviation from position
zz = xw / P_fwhm
fa = 2.D0 * zz * EXP( -zz**2)
fb = 2.D0 * (2.D0 * zz**2 - 3.) * fa
asym = 1.0D0
IF(axis) THEN
      asym = asym + (P_asym1*fa + P_asym2*fb) / tand(0.5*P_pos )
ELSE
      asym = asym + (P_asym1*fa + P_asym2*fb) / P_pos
ENDIF
!
!ww = P_inte*(  P_eta       *(2.0/PI*P_fwhm/(P_fwhm*P_fwhm + 4.0D0 * xw * xw))   &
!             +(1.0D0-P_eta)*(gpre / P_fwhm * EXP( -vln2 / P_fwhm**2 * xw**2))   &
!            )                                                                   &
!     * asym
i = MIN(INT(ABS(xw)/P_fwhm*GLP_NPT), GLP_MAX)          ! Index for lookup
ww = P_inte * glp_pseud_indx(i, P_eta, P_fwhm) * asym  ! Look up Pseudo-Voigt
!
END FUNCTION pseudo_voigt
!
!*****7************************************************************************* 
!
END MODULE lib_f90_profile
