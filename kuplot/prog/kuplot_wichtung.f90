module kuplot_wichtung_mod
!
use precision_mod
!
interface r_wichtung
  module procedure r_wichtung444, r_wichtung484, r_wichtung448, r_wichtung488,    &
                   r_wichtung844, r_wichtung884, r_wichtung848, r_wichtung888
end interface
!
contains
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung444 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_SP), intent(in) :: z
REAL(kind=PREC_SP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_SP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung444                      
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung484 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_SP), intent(in) :: z
REAL(kind=PREC_DP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_SP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung484                      
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung448 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_SP), intent(in) :: z
REAL(kind=PREC_SP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_DP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung448                       
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung488 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_SP), intent(in) :: z
REAL(kind=PREC_DP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_DP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung488                       
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung844 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), intent(in) :: z
REAL(kind=PREC_SP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_SP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung844                       
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung884 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), intent(in) :: z
REAL(kind=PREC_DP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_SP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung884                       
!
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung848 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), intent(in) :: z
REAL(kind=PREC_SP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_DP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung848                       
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) function r_wichtung888 (z, dz, iweight, zc, bck_k)  result( r_wichtung)
!                                                                       
!     calculates the weight                                             
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), intent(in) :: z
REAL(kind=PREC_DP), intent(in) :: dz 
INTEGER           , intent(in) :: iweight 
REAL(kind=PREC_DP), intent(in) :: zc 
REAL(kind=PREC_DP), intent(in) :: bck_k
!                                                                       
INTEGER, parameter :: W_ONE  = 0
INTEGER, parameter :: W_SQUA = 1
INTEGER, parameter :: W_SQRT = 2
INTEGER, parameter :: W_INV  = 3
INTEGER, parameter :: W_LOG  = 4
INTEGER, parameter :: W_ISQ  = 5
INTEGER, parameter :: W_LIN  = 6
INTEGER, parameter :: W_DAT  = 7
INTEGER , parameter :: W_BCK = 8          ! Currently not used
!     PARAMETER (W_ONE = 0) 
!     PARAMETER (W_SQUA = 1) 
!     PARAMETER (W_SQRT = 2) 
!     PARAMETER (W_INV = 3) 
!     PARAMETER (W_LOG = 4) 
!     PARAMETER (W_ISQ = 5) 
!     PARAMETER (W_LIN = 6) 
!     PARAMETER (W_DAT = 7) 
!     PARAMETER (W_BCK = 8) 
!                                                                       
         
r_wichtung = 0.0 
IF (iweight.eq.W_ONE) then 
   r_wichtung = 1.0 
ELSEIF (iweight.eq.W_SQUA) then 
   r_wichtung = z**2 
ELSEIF (iweight.eq.W_SQRT) then 
   IF (z.ge.0) then 
      r_wichtung = sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_INV) then 
   IF (z.ne.0) then 
      r_wichtung = 1.0 / abs (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LOG) then 
   IF (z.gt.0) then 
      r_wichtung = log (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_ISQ) then 
   IF (z.gt.0) then 
      r_wichtung = 1. / sqrt (z) 
   ELSE 
      r_wichtung = 0.0 
   ENDIF 
ELSEIF (iweight.eq.W_LIN) then 
   r_wichtung = z 
ELSEIF (iweight.eq.W_DAT) then 
   IF (dz.ne.0) then 
      r_wichtung = 1.0 / dz**2
   ELSE 
      r_wichtung = 1.0
   ENDIF 
elseif(iweight==W_BCK) then
   if((z-zc)>0.0) then
      r_wichtung = dz*(0.8*bck_k)
   else
      r_wichtung = dz/(0.8*bck_k)
   endif
!           r_wichtung = exp(-abs(bck_k)*(z-zc))
ENDIF 
!                                                                       
END FUNCTION r_wichtung888                       
!
end module kuplot_wichtung_mod
