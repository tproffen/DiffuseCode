MODULE trig_degree_mod
!
use precision_mod
!
IMPLICIT NONE
!
INTERFACE sind
   MODULE PROCEDURE sind_real, sind_real8
END INTERFACE sind
INTERFACE cosd
   MODULE PROCEDURE cosd_real, cosd_real8
END INTERFACE cosd
INTERFACE tand
   MODULE PROCEDURE tand_real, tand_real8
END INTERFACE tand
INTERFACE asind
   MODULE PROCEDURE asind_real, asind_real8
END INTERFACE asind
INTERFACE acosd
   MODULE PROCEDURE acosd_real, acosd_real8
END INTERFACE acosd
INTERFACE atand
   MODULE PROCEDURE atand_real, atand_real8
END INTERFACE atand
INTERFACE atan2d
   MODULE PROCEDURE atan2d_real, atan2d_real8
END INTERFACE atan2d
CONTAINS
!                                                                       
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION sind_real (arg) 
!
USE wink_mod

IMPLICIT NONE
REAL(kind=PREC_SP), INTENT(IN) :: arg
!
sind_real = REAL(sin (arg * rad) )
END FUNCTION sind_real
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION cosd_real (arg) 
!
USE wink_mod
REAL(kind=PREC_SP), INTENT(IN) :: arg
!
cosd_real = REAL(cos (arg * rad) )
END FUNCTION cosd_real
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION tand_real (arg) 
!
USE wink_mod
REAL(kind=PREC_SP), INTENT(IN) :: arg
!
tand_real = REAL(tan (arg * rad) )
END FUNCTION tand_real
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION asind_real (arg) 
!
USE wink_mod
REAL(kind=PREC_SP), INTENT(IN) :: arg
!
asind_real = REAL(asin (arg) / rad )
END FUNCTION asind_real                            
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION acosd_real (arg) 
!
USE wink_mod
REAL(kind=PREC_SP), INTENT(IN) :: arg
!
acosd_real = REAL(acos (arg) / rad )
END FUNCTION acosd_real                            
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION atand_real(arg) 
!
USE wink_mod
REAL(kind=PREC_SP), INTENT(IN) :: arg
!
atand_real = REAL(atan (arg) / rad )
END FUNCTION atand_real                            
!*****7**************************************************************** 
REAL(kind=PREC_SP) FUNCTION atan2d_real (arg1, arg2) 
!
USE wink_mod
REAL(kind=PREC_SP), INTENT(IN) :: arg1, arg2
!
atan2d_real = REAL(atan2 (arg1, arg2) / rad )
END FUNCTION atan2d_real
!*****7**************************************************************** 
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION sind_real8 (arg) 
!
USE wink_mod
IMPLICIT NONE
REAL(KIND=PREC_DP), INTENT(IN) :: arg
!
sind_real8 = sin (arg * rad) 
END FUNCTION sind_real8
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION cosd_real8 (arg) 
!
USE wink_mod
REAL(KIND=PREC_DP), INTENT(IN) :: arg
!
cosd_real8 = cos (arg * rad) 
END FUNCTION cosd_real8
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION tand_real8 (arg) 
!
USE wink_mod
REAL(KIND=PREC_DP), INTENT(IN) :: arg
!
tand_real8 = tan (arg * rad) 
END FUNCTION tand_real8
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION asind_real8 (arg) 
!
USE wink_mod
REAL(KIND=PREC_DP), INTENT(IN) :: arg
!
asind_real8 = asin (arg) / rad 
END FUNCTION asind_real8                            
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION acosd_real8 (arg) 
!
USE wink_mod
REAL(KIND=PREC_DP), INTENT(IN) :: arg
!
acosd_real8 = acos (arg) / rad 
END FUNCTION acosd_real8                            
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION atand_real8 (arg) 
!
USE wink_mod
REAL(KIND=PREC_DP), INTENT(IN) :: arg
!
atand_real8 = atan (arg) / rad 
END FUNCTION atand_real8                            
!*****7**************************************************************** 
REAL(KIND=PREC_DP) FUNCTION atan2d_real8 (arg1, arg2) 
!
USE wink_mod
REAL(KIND=PREC_DP), INTENT(IN) :: arg1, arg2
!
atan2d_real8 = atan2 (arg1, arg2) / rad 
END FUNCTION atan2d_real8
!
END MODULE trig_degree_mod
