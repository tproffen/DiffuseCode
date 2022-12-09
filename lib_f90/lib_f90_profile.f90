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
                                                 axis, lambda) RESULT(ww)
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
REAL(KIND=PREC_DP), INTENT(IN)  :: xx            ! Position at which to calculate Pseduovoigt
REAL(KIND=PREC_DP), INTENT(IN)  :: P_eta         ! Profile Mixing parameter
REAL(KIND=PREC_DP), INTENT(IN)  :: P_inte        ! Integrated intensity
REAL(KIND=PREC_DP), INTENT(IN)  :: P_pos         ! Position
REAL(KIND=PREC_DP), INTENT(IN)  :: P_fwhm        ! FWHM
REAL(KIND=PREC_DP), INTENT(IN)  :: P_asym1       ! Asymmetry parameter 1
REAL(KIND=PREC_DP), INTENT(IN)  :: P_asym2       ! Asymmetry parameter 2
LOGICAL           , INTENT(IN)  :: axis          ! TRUE if TTH Scale
REAL(KIND=PREC_DP), INTENT(IN)  :: lambda        ! Wave length
!
!REAL(KIND=PREC_DP) :: ww
!
!REAL(KIND=PREC_DP), PARAMETER:: vln2 = 4.D0 * LOG (2.D0)
!REAL(KIND=PREC_DP), PARAMETER:: gpre = 2.D0 * SQRT (LOG (2.D0) / PI)
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
      asym = asym + (P_asym1*fa + P_asym2*fb) / tanh(RAD*0.5*P_pos )
ELSE
      asym = asym + (P_asym1*fa + P_asym2*fb) / tanh(asin(P_pos/FPI/lambda))
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
!*******************************************************************************
!
real(kind=PREC_DP) function pearson_vii(xx, P_m, P_inte, P_pos,     &
                                        P_fwhm, P_asym1, P_asym2,   &
                                        axis, lambda) RESULT(ww)
!-
! Caclulate a Pearson Type VII function 
!+
!
use precision_mod
use gamma_mod
use wink_mod
!
implicit none
!
real(KIND=PREC_DP), INTENT(IN)  :: xx            ! Position at which to calculate Pseduovoigt
real(KIND=PREC_DP), INTENT(IN)  :: P_m           ! Profile Mixing parameter
real(KIND=PREC_DP), INTENT(IN)  :: P_inte        ! Integrated intensity
real(KIND=PREC_DP), INTENT(IN)  :: P_pos         ! Position
real(KIND=PREC_DP), INTENT(IN)  :: P_fwhm        ! FWHM
REAL(KIND=PREC_DP), INTENT(IN)  :: P_asym1       ! Asymmetry parameter 1
REAL(KIND=PREC_DP), INTENT(IN)  :: P_asym2       ! Asymmetry parameter 2
LOGICAL           , INTENT(IN)  :: axis          ! TRUE if TTH Scale
REAL(KIND=PREC_DP), INTENT(IN)  :: lambda        ! Wave length
!
real(kind=PREC_DP) :: alpha
!
REAL(KIND=PREC_DP)  :: xw             ! Deviation from position 
REAL(KIND=PREC_DP)  :: zz             ! Value for asymmetry function
REAL(KIND=PREC_DP)  :: fa             ! Value for asymmetry function
REAL(KIND=PREC_DP)  :: fb             ! Value for asymmetry function
REAL(KIND=PREC_DP)  :: asym           ! Value for asymmetry function
!
!lpha = P_fwhm *0.5D0 / sqrt(0.5D0**(-1.0D0/P_m) - 1.0D0)
alpha = P_fwhm *0.5D0 / sqrt(2.0D0**( 1.0D0/P_m) - 1.0D0)
!
xw = xx - P_pos                 ! Deviation from position
zz = xw / P_fwhm
fa = 2.D0 * zz * EXP( -zz**2)
fb = 2.D0 * (2.D0 * zz**2 - 3.) * fa
asym = 1.0D0
if(P_asym1/=0.0D0 .or. P_asym2/=0.0D0) then
IF(axis) THEN
      asym = asym + (P_asym1*fa + P_asym2*fb) / tanh(RAD*0.5*P_pos )
ELSE
      asym = asym + (P_asym1*fa + P_asym2*fb) / tanh(asin(P_pos/FPI/lambda))
ENDIF
endif
!
ww = P_inte * ( (1.0D0+((xx-P_pos)/alpha)**2)**(-P_m) )/ &
              (alpha*func_beta(P_m-0.5D0, 0.5D0)) * asym
!write(*,'(7f10.4)') xx, P_pos, ( (1.0D0+((xx-P_pos)/alpha)**2)**(-P_m) )/ &
!              (alpha*func_beta(P_m-0.5D0, 0.5D0)), (alpha*func_beta(P_m-0.5D0, 0.5D0)), &
!              P_m, P_fwhm, alpha
!
end function pearson_vii
!
!*******************************************************************************
!
real(kind=PREC_DP) function tukey(x, alpha, width)
!-
!  Return a value for the Tukey window function
!
!  t(x,alpha) = 1/2 [ 1 - cos(2PIx/(alpha*width) ]	0             <= x <  alpha*width/2
!             = 1                                       alpha*width/2 <= x <= width/2
!             = 1/2 [ 1 - cos(2PIx/(alpha*width) ]	width/2       <= x <= width
!+
!
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP), intent(in) :: x 
real(kind=PREC_DP), intent(in) :: alpha 
real(kind=PREC_DP), optional, intent(in) :: width 
!
real(kind=PREC_DP) :: xl
real(kind=PREC_DP) :: w
!
if(present(width)) then
  w = width
else
  w = 1.0D0
endif
!
xl = x
if(xl>0.5D0*w) xl = w-xl
!
tukey = 0.0D0
if(xl<0.0D0) then
   tukey = 0.0D0
elseif(xl<0.5D0*w*alpha ) then
   tukey = 0.5D0*(1.0D0 - cos(ZPI*xl/(alpha*w)))
else
   tukey = 1.0D0
endif
!
end function tukey
!
!*****7************************************************************************* 
!
END MODULE lib_f90_profile
