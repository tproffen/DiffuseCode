MODULE trig_degree_mod
IMPLICIT NONE
!
INTEGER, PARAMETER :: dp = KIND(1.0D0)
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
REAL FUNCTION sind_real (arg) 
!
USE wink_mod
IMPLICIT NONE
REAL, INTENT(IN) :: arg
!
sind_real = REAL(sin (arg * rad) )
END FUNCTION sind_real
!*****7**************************************************************** 
REAL FUNCTION cosd_real (arg) 
!
USE wink_mod
REAL, INTENT(IN) :: arg
!
cosd_real = REAL(cos (arg * rad) )
END FUNCTION cosd_real
!*****7**************************************************************** 
REAL FUNCTION tand_real (arg) 
!
USE wink_mod
REAL, INTENT(IN) :: arg
!
tand_real = REAL(tan (arg * rad) )
END FUNCTION tand_real
!*****7**************************************************************** 
REAL FUNCTION asind_real (arg) 
!
USE wink_mod
REAL, INTENT(IN) :: arg
!
asind_real = REAL(asin (arg) / rad )
END FUNCTION asind_real                            
!*****7**************************************************************** 
REAL FUNCTION acosd_real (arg) 
!
USE wink_mod
REAL, INTENT(IN) :: arg
!
acosd_real = REAL(acos (arg) / rad )
END FUNCTION acosd_real                            
!*****7**************************************************************** 
REAL FUNCTION atand_real(arg) 
!
USE wink_mod
REAL, INTENT(IN) :: arg
!
atand_real = REAL(atan (arg) / rad )
END FUNCTION atand_real                            
!*****7**************************************************************** 
REAL FUNCTION atan2d_real (arg1, arg2) 
!
USE wink_mod
REAL, INTENT(IN) :: arg1, arg2
!
atan2d_real = REAL(atan2 (arg1, arg2) / rad )
END FUNCTION atan2d_real
!*****7**************************************************************** 
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION sind_real8 (arg) 
!
USE wink_mod
IMPLICIT NONE
REAL (KIND=dp), INTENT(IN) :: arg
!
sind_real8 = sin (arg * rad) 
END FUNCTION sind_real8
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION cosd_real8 (arg) 
!
USE wink_mod
REAL (KIND=dp), INTENT(IN) :: arg
!
cosd_real8 = cos (arg * rad) 
END FUNCTION cosd_real8
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION tand_real8 (arg) 
!
USE wink_mod
REAL (KIND=dp), INTENT(IN) :: arg
!
tand_real8 = tan (arg * rad) 
END FUNCTION tand_real8
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION asind_real8 (arg) 
!
USE wink_mod
REAL (KIND=dp), INTENT(IN) :: arg
!
asind_real8 = asin (arg) / rad 
END FUNCTION asind_real8                            
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION acosd_real8 (arg) 
!
USE wink_mod
REAL (KIND=dp), INTENT(IN) :: arg
!
acosd_real8 = acos (arg) / rad 
END FUNCTION acosd_real8                            
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION atand_real8 (arg) 
!
USE wink_mod
REAL (KIND=dp), INTENT(IN) :: arg
!
atand_real8 = atan (arg) / rad 
END FUNCTION atand_real8                            
!*****7**************************************************************** 
REAL (KIND=dp) FUNCTION atan2d_real8 (arg1, arg2) 
!
USE wink_mod
REAL (KIND=dp), INTENT(IN) :: arg1, arg2
!
atan2d_real8 = atan2 (arg1, arg2) / rad 
END FUNCTION atan2d_real8
!
END MODULE trig_degree_mod
