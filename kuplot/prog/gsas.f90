module kuplot_gsas_mod
!
contains
!
      SUBROUTINE PRPCALC (HTYPE, TTH, PTYPE, PC, DTOF, FUNC, DV, NMAX) 
                                                                        
!PURPOSE: Calculating individual peak profile shape function &          
! derivative                                                            
                                                                        
!Written by: Allen C. Larson      and Robert B. Von Dreele              
!            MS-H805                                                    
!            Los Alamos National Laboratory                             
!            Los Alamos, NM  87545                                      
!                                                                       
!            Copyright, 1984-2000, The Regents of the University        
! of California.                                                        
!                                                                       
! This software was produced under a U.S. Government contract           
! (W-7405-ENG-36)                                                       
! by the Los Alamos National Laboratory, which is operated by the       
! University                                                            
! of California for the U.S. Department of Energy. The U.S.             
! Government is                                                         
! licensed to use, reproduce, and distribute this software.             
! Permission is                                                         
! granted to the public to copy and use this software without           
! charge,                                                               
! provided that this notice and any statement of authorship are         
! reproduced                                                            
! on all copies. Neither the Government nor the University makes        
! any warranty,                                                         
! express or implied, or assumes any liability or responsibility        
! for the use                                                           
! of this software.                                                     
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!PSEUDOCODE:                                                            
                                                                        
!CALLING ARGUMENTS:                                                     
INTEGER, INTENT(IN) :: NMAX
                                                                        
                              !Histogram type                           
      CHARACTER(4) HTYPE 
                              !2-theta in centideg - only for CW        
      REAL TTH 
                              !Profile function type                    
      INTEGER PTYPE 
                              !Coefficients                             
      REAL PC (NMAX) 
                              !=TOF-TOFR                                
      REAL DTOF 
                              !Value of profile function                
      REAL FUNC 
                              !Function derivatives                     
      REAL DV (NMAX) 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
!REAL :: DEXDS, dh1dy, dh2dy, dy1ds, dy2ds, s2sg, sq2pi
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
!CODE:                                                                  
                                                                        
      IF (HTYPE (3:3) .EQ.'T') THEN 
                                                                        
!Back-to-back exponential-gaussian convolution function                 
!R.B. Von Dreele, J. Jorgensen & C. Windsor, J. Appl. Cryst., , , 1982. 
                                                                        
         IF (PTYPE.EQ.1) THEN 
            CALL EXPGAUS1 (DTOF * 1000.0, PC (3), PC (4), PC (5),       &
            FUNC, DV (1), DV (3), DV (4), DV (5) )                      
                                                                        
!W.I.F. David convolution function of Ikeda-Carpenter with pseudo Voight
!Unpublished ( 1986).                                                   
                                                                        
         ELSEIF (PTYPE.EQ.2) THEN 
            CALL IKDAVID (PC (3), PC (4), PC (5), PC (6), PC (7),       &
            DTOF * 1000.0, FUNC, DV (1), DV (3), DV (4), DV (5),        &
            DV (6), DV (7) )                                            
                                                                        
!pseudo Voigt W.I.F.David, J. Appl. Cryst. 19, 63-64                    
                                                                        
         ELSEIF (PTYPE.EQ.3) THEN 
            CALL EPSVOIGT (DTOF * 1000.0, PC (3), PC (4), PC (5),       &
            PC (6), FUNC, DV (1), DV (3), DV (4), DV (5), DV (6) )      
                                                                        
         ELSEIF (PTYPE.EQ.4) THEN 
            CALL EPSVOIGT (DTOF * 1000.0, PC (3), PC (4), PC (5),       &
            PC (6), FUNC, DV (1), DV (3), DV (4), DV (5), DV (6) )      
                                                                        
         ELSE 
            PRINT 1, PTYPE 
    1 FORMAT    (' **** Unknown peak shape function',I3,                &
     &      ' - Programmer error ****')                                 
            STOP 
         ENDIF 
         DV (1) = DV (1) * 1000.0 
                                                                        
      ELSEIF (HTYPE (2:3) .EQ.'NC') THEN 
                                                                        
!Simpson's rule integration of Gaussian function                        
!C.J. Howard (1982). J. Appl. Cryst.,15 ,615-620.                       
! Cooper & Sayer (1975). J. Appl. Cryst., 8, 615-618.                   
! Thomas (1977). J. Appl. Cryst., 10, 12-13.                            
! All calcs in centidegrees                                             
                                                                        
         IF (PTYPE.EQ.1) THEN 
            CALL SIMGAUS1 (DTOF * 100.0, PC (3), PC (4), FUNC, DV (1),  &
            DV (3), DV (4) )                                            
                                                                        
                                                                        
!Simpson's rule integration of pseudovoigt function                     
!pseudo Voigt W.I.F.David, J. Appl. Cryst. 19, 63-64 &                  
!asymmetry correction C.J. Howard, J. Appl. Cryst. (1982), 15, 615-620. 
                                                                        
         ELSEIF (PTYPE.EQ.2) THEN 
            CALL PSCVGT1 (DTOF * 100.0, PC (3), PC (4), PC (5), FUNC,   &
            DV (1), DV (3), DV (4), DV (5) )                            
                                                                        
!pseudo Voigt W.I.F.David, J. Appl. Cryst. 19, 63-64 and                
! P. Thompson, D.E. Cox & J.B. Hastings (1987) J. Appl. Cryst.,20,79-83.
!asymmetry correction L.W. Finger, D.E. Cox & A.P. Jephcoat (1994)      
! J. Appl. Cryst.,27,892-900.                                           
                                                                        
         ELSEIF (PTYPE.EQ.3) THEN 
            CALL PSVFCJ (DTOF * 100.0, TTH, PC (3), PC (4), PC (5),     &
            PC (6), FUNC, DV (1), DV (3), DV (4), DV (5), DV (6) )      
                                                                        
         ELSEIF (PTYPE.EQ.4) THEN 
            CALL PSVFCJ (DTOF * 100.0, TTH, PC (3), PC (4), PC (5),     &
            PC (6), FUNC, DV (1), DV (3), DV (4), DV (5), DV (6) )      
                                                                        
         ELSE 
            PRINT 1, PTYPE 
            STOP 
         ENDIF 
         DV (1) = DV (1) * 100.0 
                                                                        
      ELSEIF (HTYPE (2:3) .EQ.'XC') THEN 
                                                                        
!Simpson's rule integration of Gaussian function                        
!C.J. Howard (1982). J. Appl. Cryst.,15 ,615-620.                       
! Cooper & Sayer (1975). J. Appl. Cryst., 8, 615-618.                   
! Thomas (1977). J. Appl. Cryst., 10, 12-13.                            
! All calcs in centidegrees                                             
                                                                        
         IF (PTYPE.EQ.1) THEN 
            CALL SIMGAUS1 (DTOF * 100.0, PC (3), PC (4), FUNC, DV (1),  &
            DV (3), DV (4) )                                            
                                                                        
!Simpson's rule integration of pseudovoigt function                     
!pseudo Voigt W.I.F.David, J. Appl. Cryst. 19, 63-64 &                  
!asymmetry correction C.J. Howard, J. Appl. Cryst. (1982), 15, 615-620. 
                                                                        
         ELSEIF (PTYPE.EQ.2) THEN 
            CALL PSCVGT1 (DTOF * 100.0, PC (3), PC (4), PC (5), FUNC,   &
            DV (1), DV (3), DV (4), DV (5) )                            
                                                                        
!pseudo Voigt W.I.F.David, J. Appl. Cryst. 19, 63-64 and                
! P. Thompson, D.E. Cox & J.B. Hastings (1987) J. Appl. Cryst.,20,79-83.
!asymmetry correction L.W. Finger, D.E. Cox & A.P. Jephcoat (1994)      
! J. Appl. Cryst.,27,892-900.                                           
                                                                        
         ELSEIF (PTYPE.EQ.3) THEN 
            CALL PSVFCJ (DTOF * 100.0, TTH, PC (3), PC (4), PC (5),     &
            PC (6), FUNC, DV (1), DV (3), DV (4), DV (5), DV (6) )      
                                                                        
         ELSEIF (PTYPE.EQ.4) THEN 
            CALL PSVFCJ (DTOF * 100.0, TTH, PC (3), PC (4), PC (5),     &
            PC (6), FUNC, DV (1), DV (3), DV (4), DV (5), DV (6) )      
                                                                        
         ELSE 
            PRINT 1, PTYPE 
            STOP 
         ENDIF 
         DV (1) = DV (1) * 100.0 
                                                                        
      ELSEIF (HTYPE (2:3) .EQ.'XE') THEN 
                                                                        
         IF (PTYPE.EQ.1) THEN 
            CALL EDSGAUS1 (DTOF, PC (3), FUNC, DV (1), DV (3) ) 
                                                                        
         ELSE 
            PRINT 1, PTYPE 
            STOP 
         ENDIF 
                                                                        
      ENDIF 
      DV (2) = FUNC 
                                                                        
      RETURN 
      END SUBROUTINE PRPCALC                        
                                                                        
!-----------------------------------------------------------------      
                                                                        
      SUBROUTINE EXPGAUS1 (DT, ALP, BET, SIG, PRFUNC, DPRDT, ALPART,    &
      BEPART, SGPART)                                                   
                                                                        
!PURPOSE:                                                               
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta T for profile point      
      REAL DT 
                                        !Profile shape coefficients     
      REAL ALP, BET, SIG 
                                        !Profile function value         
      REAL PRFUNC 
                                        !Position derivative            
      REAL DPRDT 
                                        !Profile shape derivatives      
      REAL ALPART, BEPART, SGPART 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
      REAL TJ1, TJ2, NORM, SQSG, Y1, Y2, H1, H2, EX, EX1 
                              !Intermediate values                      
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
                                        !erfc(x)exp(u^2)                
!     REAL HFUNC 
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA SQ2PI / 2.506628275 /, SQ2 / 1.414213562 / 
                                                                        
!CODE:                                                                  
                                                                        
      TJ1 = ALP * SIG 
      TJ2 = BET * SIG 
      NORM = 0.5 * ALP * BET / (ALP + BET) 
      SQSG = SQRT (2.0 * SIG) 
      S2SG = 2.0 * SQSG * SIG 
      Y1 = (TJ1 + DT) / SQSG 
      DY1DS = (TJ1 - DT) / S2SG 
      Y2 = (TJ2 - DT) / SQSG 
      DY2DS = (TJ2 + DT) / S2SG 
      EX1 = - 0.5 * (DT**2) / SIG 
      DEXDS = - EX1 / SIG 
      H1 = HFUNC (EX1, Y1, 0, DH1DY) 
      H2 = HFUNC (EX1, Y2, 0, DH2DY) 
                                                                        
      PRFUNC = NORM * (H1 + H2) 
                                                                        
                                                  !Compute derivatives  
      EX = 0.0 
      IF (EX1.GE. - 50.0) EX = EXP (EX1) 
                                                                        
      ALPART = 2.0 * NORM * PRFUNC / (ALP * ALP) 
      ALPART = ALPART + 0.5 * NORM * SQSG * DH1DY 
      BEPART = 2.0 * NORM * PRFUNC / (BET * BET) 
      BEPART = BEPART + 0.5 * NORM * SQSG * DH2DY 
      SGPART = PRFUNC * DEXDS + NORM * (DH1DY * DY1DS + DH2DY * DY2DS) 
                                               !df/dt                   
      DPRDT = - NORM * (ALP * H1 - BET * H2) 
                                                                        
      RETURN 
      END SUBROUTINE EXPGAUS1                       
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE IKDAVID (ALP, BET, RAT, SIG, GAM, DELTOF, PRFUNC,      &
      DPRDT, ALPT, BEPT, RATPT, SIGPT, GAMPT)                           
                                                                        
!PURPOSE: W.I.F.DAVID convolution function of                           
                            !Ikeda-Carpenter with pseudo Voight         
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Peak coeffs.                   
      REAL ALP, BET, SIG, GAM, RAT 
                                        !Delta TTOF                     
      REAL DELTOF 
                                        !Value of function at DELTOF    
      REAL PRFUNC 
                                        !partial df(t)/dt               
      REAL DPRDT 
                                                !Derivatives            
      REAL ALPT, BEPT, RATPT, SIGPT, GAMPT 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
      REAL COFT (6), COFN (3) 
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
                                         !SQRT(8LN2)                    
      DATA STOFW / 2.35482005 / 
      DATA COFT / 1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0 / 
      DATA COFN / 1.36603, - 0.47719, 0.11116 / 
                                                                        
!CODE:                                                                  
                                                                        
      SIG2 = SQRT (SIG) 
      FWHG = STOFW * SIG2 
      PGL = FWHG**5 
      SUMHM = PGL 
      DSDL = 0.0 
      DSDG = 0.0 
      DO ITRM = 1, 5 
      PGL = PGL / FWHG 
      DSDL = DSDL + REAL(ITRM) * COFT (ITRM + 1) * PGL 
      DSDG = DSDG + REAL(6 - ITRM) * COFT (ITRM) * PGL 
      PGL = PGL * GAM 
      SUMHM = SUMHM + COFT (ITRM + 1) * PGL 
      ENDDO 
      FWHM = EXP (0.2 * LOG (SUMHM) ) 
      FRAC = GAM / FWHM 
      DEDF = 0.0 
      PF = 1.0 
      ETA = 0.0 
      DO ITRM = 1, 3 
      DEDF = DEDF + REAL(ITRM) * COFN (ITRM) * PF 
      PF = PF * FRAC 
      ETA = ETA + COFN (ITRM) * PF 
      ENDDO 
      SIGP = (FWHM / STOFW) **2 
      CALL IKGAUSS (DELTOF, ALP, BET, SIGP, RAT, PRFUNCG, DPRDTG, ALPTG,&
      BEPTG, SIGPTG, RATPTG)                                            
      CALL IKLRNTZ (DELTOF, ALP, BET, FWHM, RAT, PRFUNCL, DPRDTL, ALPTL,&
      BEPTL, GAMPTL, RATPTL)                                            
                                                                        
      SIGPTG = 2.0 * SIGPTG * FWHM / STOFW**2 
      ALPT = (1.0 - ETA) * ALPTG + ETA * ALPTL 
      BEPT = (1.0 - ETA) * BEPTG + ETA * BEPTL 
      RATPT = (1.0 - ETA) * RATPTG + ETA * RATPTL 
                                                                        
      DFWDG = 0.2 * DSDG * FWHM / SUMHM 
      DFRDG = - FRAC * DFWDG / FWHM 
      SIGPT = DEDF * DFRDG * (PRFUNCL - PRFUNCG) + DFWDG * (ETA *       &
      GAMPTL + (1.0 - ETA) * SIGPTG)                                    
      SIGPT = 0.5 * SIGPT * STOFW / SIG2 
                                                                        
      DFWDL = 0.2 * DSDL * FWHM / SUMHM 
      DFRDL = (1.0 - FRAC * DFWDL) / FWHM 
      GAMPT = DEDF * DFRDL * (PRFUNCL - PRFUNCG) + DFWDL * (ETA *       &
      GAMPTL + (1.0 - ETA) * SIGPTG)                                    
                                                                        
      DPRDT = (1.0 - ETA) * DPRDTG + ETA * DPRDTL 
      DPRDT = - DPRDT 
      PRFUNC = (1.0 - ETA) * PRFUNCG + ETA * PRFUNCL 
                                                                        
      RETURN 
      END SUBROUTINE IKDAVID                        
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE IKGAUSS (DT, ALP, BET, SIG, RAT, FUNC, DFDT, DFDA,     &
      DFDB, DFDS, DFDR)                                                 
                                                                        
!PURPOSE:Compute func. & derivs. for Von Dreele convolution             
!      of modified                                                      
!      Ikeda-Carpenter function with Gaussian & WIFD function           
!      derivatives                                                      
!      All checked for correctness - 12/5/1986                          
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta TOF in micro secs.       
      REAL DT 
                                        !Fast rise coefficient          
      REAL ALP 
                                        !Slow rise coefficient          
      REAL BET 
                                        !Gaussian variance              
      REAL SIG 
                                        !Slow fraction                  
      REAL RAT 
                                        !Value of function at DT        
      REAL FUNC 
      REAL DFDT, DFDA, DFDB, DFDS, DFDR 
                              !Derivatives of variables                 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
REAL :: A1, A2, A3, A4, AMD, APD, C0, C1, C2, C3, C4, C1R
REAL :: D1DQ, D2DV, D3DU, D4DP, DEDS, DEDT, DFDB1, DFDB2
REAL :: DK, DKM, DKP, DKSQ, DPDS, DQDS, DUDS, DVDS
REAL :: EX, RSQ2PI, SQ2, SQPI, SQSG
REAL :: TS, X, XYZ, Y, Z, ZP, ZQ, ZU, ZV
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
                            !Computes exp(-t*t/2s)exp(y*y)erfc(y)       
!     REAL HFUNC 
                            !and derivs                                 
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
                                        !Sqrt(pi)                       
      DATA SQPI / 1.772453851 / 
                                         !1/sqrt(2pi)                   
      DATA RSQ2PI / 0.39894228 / 
                                       !Sqrt(2)                         
      DATA SQ2 / 1.414213562 / 
                                     !Alpha interval                    
      DATA DK / 0.05 / 
                                   !DK**2                               
      DATA DKSQ / 0.0025 / 
                                         !1-DK,1+DK                     
      DATA DKM, DKP / 0.95, 1.05 / 
                                                                        
!CODE:                                                                  
                                                                        
      DFDT = 0.0 
      DFDA = 0.0 
      DFDB = 0.0 
      DFDS = 0.0 
      DFDR = 0.0 
      FUNC = 0.0 
                                                                        
      AMD = ALP * DKM 
      APD = ALP * DKP 
      X = AMD-BET 
      Y = ALP - BET 
      Z = APD-BET 
      XYZ = X * Y * Z 
      SQSG = 1.0 / SQRT (2.0 * SIG) 
      C0 = ALP * (1.0 - DKSQ) / (4.0 * DKSQ) 
      C1R = 2.0 * BET * DKSQ * ALP**2 / XYZ 
      C1 = C1R * RAT 
      C2 = RAT * APD / Z - 1.0 
      C3 = RAT * AMD / X - 1.0 
      C4 = 2.0 * (RAT * ALP / Y - 1.0) 
      TS = DT * SQSG 
      EX = - TS * TS 
      ZU = AMD * SIG * SQSG - TS 
      ZV = APD * SIG * SQSG - TS 
      ZP = ALP * SIG * SQSG - TS 
      ZQ = BET * SIG * SQSG - TS 
      IF (ZQ.GE.6.375) THEN 
         A1 = HFUNC (EX, ZQ, 1, D1DQ) 
      ELSE 
         A1 = HFUNC (EX, ZQ, 0, D1DQ) 
      ENDIF 
      IF (AMAX1 (ZU, ZV, ZP) .GE.6.375) THEN 
         A2 = HFUNC (EX, ZV, 1, D2DV) 
         A3 = HFUNC (EX, ZU, 1, D3DU) 
         A4 = HFUNC (EX, ZP, 1, D4DP) 
      ELSE 
         A2 = HFUNC (EX, ZV, 0, D2DV) 
         A3 = HFUNC (EX, ZU, 0, D3DU) 
         A4 = HFUNC (EX, ZP, 0, D4DP) 
      ENDIF 
      FUNC = C0 * (C1 * A1 + C4 * A4 - C2 * A2 - C3 * A3) 
                                                                        
      DEDT = - DT / SIG 
      DEDS = 0.5 * DEDT**2 
      DUDS = AMD-DEDT 
      DVDS = APD-DEDT 
      DPDS = ALP - DEDT 
      DQDS = BET - DEDT 
                                                                        
      DFDR = C0 * (C1R * A1 + 2.0 * ALP * A4 / Y - APD * A2 / Z - AMD * &
      A3 / X)                                                           
      DFDT = FUNC * DEDT - C0 * SQSG * (C1 * D1DQ + C4 * D4DP - C2 *    &
      D2DV - C3 * D3DU)                                                 
                                                                        
      DFDS = FUNC * DEDS + C0 * SQSG * (C1 * D1DQ * DQDS + C4 * D4DP *  &
      DPDS - C2 * D2DV * DVDS - C3 * D3DU * DUDS) / 2.0                 
                                                                        
      DFDB1 = C1 * (A1 * (1.0 / BET + 1.0 / X + 1.0 / Y + 1.0 / Z)      &
      + SIG * SQSG * D1DQ)                                              
      DFDB2 = APD * A2 / Z**2 + AMD * A3 / X**2 - 2.0 * ALP * A4 / Y**2 
      DFDB = C0 * (DFDB1 - RAT * DFDB2) 
                                                                        
      DFDA = FUNC / ALP + C0 * (C1 * A1 * (2.0 / ALP - DKM / X - 1.0 /  &
      Y - DKP / Z) + RAT * BET * ( - 2.0 * A4 / Y**2 + DKP * A2 / Z**2 +&
      DKM * A3 / X**2) + SIG * SQSG * (C4 * D4DP - C2 * D2DV * DKP - C3 &
      * D3DU * DKM) )                                                   
                                                                        
      RETURN 
      END SUBROUTINE IKGAUSS                        
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE IKLRNTZ (DT, ALP, BET, GAM2, RAT, FUNC, DFDT, DFDA,    &
      DFDB, DFDG, DFDR)                                                 
                                                                        
!PURPOSE: Calculate function and derivatives of                         
! Ikeda-Carpenter-Lorentzian function                                   
!Adapted from similar code written by W.I.F. David 26-Jan-1985          
!Derivatives added by R.B. Von Dreele 1-Aug-1986                        
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta TOF in musec             
      REAL DT 
                                        !Fast decay constant            
      REAL ALP 
                                        !Slow decay constant            
      REAL BET 
                                        !Lorentzian fwhm                
      REAL GAM2 
                                        !Slow fraction                  
      REAL RAT 
                                        !Value of function              
      REAL FUNC 
      REAL DFDT, DFDA, DFDB 
                                        !Derivatives                    
      REAL DFDG, DFDR 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                     
REAL :: AMD, APD
REAL :: C0, C1, C1R, CX, CY, CZ, DDT, DGT, DK, DKM, DKP                                  
REAL :: DKSQ
REAL :: FIX, FIY, FIZ, FRX, FRY, FRZ, GAM, PI, SI, SR
REAL :: X, XY, XYZ, XZ, Y, YZ, Z
REAL :: ZFIX, ZFIY, ZFIZ, ZFRX, ZFRY, ZFRZ, ZSR, ZSI
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA PI / 3.141592654 / 
                                     !Alpha interval                    
      DATA DK / 0.05 / 
                                   !DK**2                               
      DATA DKSQ / 0.0025 / 
                                         !1-DK,1+DK                     
      DATA DKM, DKP / 0.95, 1.05 / 
                                                                        
!CODE:                                                                  
                                                                        
      AMD = ALP * DKM 
      APD = ALP * DKP 
      X = AMD-BET 
      Y = ALP - BET 
      Z = APD-BET 
      XY = X * Y 
      XZ = X * Z 
      YZ = Y * Z 
      XYZ = X * Y * Z 
      GAM = GAM2 / 2.0 
      C0 = ALP * (1.0 - DKSQ) / (2.0 * PI * DKSQ) 
      C1R = 2.0 * BET * DKSQ * ALP**2 
      C1 = C1R * RAT 
      CX = YZ * (X - RAT * AMD) 
      CY = 2.0 * XZ * (Y - RAT * ALP) 
      CZ = XY * (Z - RAT * APD) 
      FRX = - DT * AMD 
      FIX = GAM * AMD 
      FRY = - DT * ALP 
      FIY = GAM * ALP 
      FRZ = - DT * APD 
      FIZ = GAM * APD 
      SR = - DT * BET 
      SI = GAM * BET 
      DDT = DT**2 + GAM**2 
      IF (DDT.GT.0.0) THEN 
         DGT = GAM / DDT 
         DDT = DT / DDT 
      ELSE 
         DGT = 0.0 
         DDT = 0.0 
      ENDIF 
      CALL EXPINT (FRX, FIX, ZFRX, ZFIX) 
      CALL EXPINT (FRY, FIY, ZFRY, ZFIY) 
      CALL EXPINT (FRZ, FIZ, ZFRZ, ZFIZ) 
      CALL EXPINT (SR, SI, ZSR, ZSI) 
                                                                        
                                                           !OK          
      FUNC = C0 * ( - C1 * ZSI - CX * ZFIX - CZ * ZFIZ + CY * ZFIY)     &
      / XYZ                                                             
                                                                        
      DFDR = C0 * ( - C1R * ZSI / XYZ + (AMD * ZFIX / X + APD * ZFIZ /  &
      Z - 2.0 * ALP * ZFIY / Y) )                                       
                                                 !OK                    
                                                                        
      DFDT = C0 * (C1 * (BET * ZSI + DGT) + CX * (AMD * ZFIX + DGT)     &
      + CZ * (APD * ZFIZ + DGT) - CY * (ALP * ZFIY + DGT) ) / XYZ       
                                                       !OK              
                                                                        
      DFDA = (FUNC / ALP) - C0 * (C1 * ZSI * (2.0 / ALP - DKM / X - DKP &
      / Z - 1.0 / Y) / XYZ + RAT * BET * (ZFIX * DKM / X**2 + ZFIZ *    &
      DKP / Z**2 - 2.0 * ZFIY / Y**2) + (CX * DKM * ( - DT * ZFIX + GAM &
      * ZFRX) + CZ * DKP * ( - DT * ZFIZ + GAM * ZFRZ) - CY * ( - DT *  &
      ZFIY + GAM * ZFRY) ) / XYZ)                                       
                                                                        
      DFDB = C1 * (ZSI * (1.0 / BET + 1.0 / X + 1.0 / Y + 1.0 / Z - DT) &
      + GAM * ZSR) / XYZ                                                
                                                    !OK                 
      DFDB = - C0 * (DFDB - RAT * (AMD * ZFIX / X**2 + APD * ZFIZ / Z** &
      2 - 2.0 * ALP * ZFIY / Y**2) )                                    
                                                                        
      DFDG = 0.5 * C0 * (C1 * ( - BET * ZSR + DDT) + CX * ( - AMD *     &
      ZFRX + DDT) + CZ * ( - APD * ZFRZ + DDT) - CY * ( - ALP * ZFRY +  &
      DDT) ) / XYZ                                                      
                                                        !OK             
      RETURN 
      END SUBROUTINE IKLRNTZ                        
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE EXPINT (AX, AY, ANSX, ANSY) 
                                                                        
!PURPOSE: Compute the real & imaginary components and derivatives of    
!the function exp(z)E(z) where E(z) is the exponential                  
! integral function.                                                    
!Originally written by       W.I.F.David      1-SEP-84                  
!Minor changes by R.B. Von Dreele 1-Aug-1986                            
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                !Real & imaginary parts of argument z   
      REAL AX, AY 
                                !Real & imaginary parts of result       
      REAL ANSX, ANSY 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
      REAL ZX, ZY 
      REAL CLNX, CLNY 
      REAL TEMX, TEMY 
      REAL SUBX, SUBY 
      REAL EULERX 
!     REAL ZANSX, ZANSY 
!     REAL DEDZX, DEDZY 
      REAL RATX, RATY 
      REAL ADDX, ADDY 
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA EULERX / 0.5772156649 / 
                                                                        
!CODE:                                                                  
                                                                        
                                           !Z=CMPLX(AX,AY)              
      ZX = AX 
      ZY = AY 
      AYCH = - 1. 
      IF (AX.LT.0.) AYCH = 0.125 * AX + 4. 
      A = SQRT (AX * AX + AY * AY) 
      IF (A.LT.4.0.OR.ABS (AY) .LT.AYCH) THEN 
         ITER = 10 + 4 * NINT (A) 
                                               !CLN=CLOG(Z)             
         CLNX = 0.5 * LOG (ZX**2 + ZY**2) 
         CLNY = ATAN2 (ZY, ZX) 
         TEMX = REAL(ITER) / REAL( (ITER + 1) * (ITER + 1) ) 
         TEMY = ZY * TEMX 
         TEMX = ZX * TEMX 
                                           !SUM=TEM                     
         SUMX = TEMX 
         SUMY = TEMY 
         DO J = ITER, 2, - 1 
         T = 1. / REAL(J) 
                                           !SUB=1.0-SUM                 
         SUBX = 1. - SUMX 
         SUBY = - SUMY 
         T = (1. + T) * (1 - T * T) 
         TEMX = TEMX * T 
         TEMY = TEMY * T 
                                                !SUM=TEM*SUB            
         SUMX = TEMX * SUBX - TEMY * SUBY 
         SUMY = TEMX * SUBY + TEMY * SUBX 
         ENDDO 
         TEMX = 1.0 - SUMX 
         TEMY = - SUMY 
         SUMX = TEMX * ZX - TEMY * ZY 
         SUMY = TEMX * ZY + TEMY * ZX 
         SUMX = - EULERX - CLNX + SUMX 
         SUMY = - CLNY + SUMY 
         EX = 0.0 
         IF (ZX.GT. - 75.0) EX = EXP (ZX) 
         EY = EX * SIN (ZY) 
         EX = EX * COS (ZY) 
                                          !ZANS=SUM*EXP(Z)              
         ANSX = SUMX * EX - SUMY * EY 
         ANSY = SUMX * EY + SUMY * EX 
      ELSE 
         IF (AX.LT.0.) A = SQRT ( (AX + 29.) * (AX + 29.) / 9.0 + AY *  &
         AY)                                                            
         ITER = 4 + NINT (128. / A) 
         TEMX = REAL(ITER) 
         ADDX = TEMX 
         ADDY = 0.0 
         DO I = 1, ITER - 1 
         TEMX = TEMX - 1. 
                                            !SUM = Z+ADD                
         SUMX = ZX + ADDX 
         SUMY = ZY + ADDY 
         X2PY2 = SUMX**2 + SUMY**2 
                                      !RAT = TEM/SUM                    
         RATX = TEMX * SUMX / X2PY2 
         RATY = - TEMX * SUMY / X2PY2 
         RATX = 1.0 + RATX 
         X2PY2 = RATX**2 + RATY**2 
                                      !ADD = TEM/(1.+RAT))              
         ADDX = TEMX * RATX / X2PY2 
         ADDY = - TEMX * RATY / X2PY2 
         ENDDO 
                                            !Z=A+ADD                    
         ZX = ZX + ADDX 
         ZY = ZY + ADDY 
         X2PY2 = ZX**2 + ZY**2 
                                              !ZANS= 1./Z               
         ANSX = ZX / X2PY2 
         ANSY = - ZY / X2PY2 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE EXPINT                         
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE EPSVOIGT (DT, ALP, BET, SIG, GAM, FUNC, DFDX, DFDA,    &
      DFDB, DFDS, DFDG)                                                 
                                                                        
!PURPOSE: Compute function & derivatives exponential X pseudovoigt      
!P.Tompson, D.E. Cox & J.B. Hastings (1987). J. Appl. Cryst.,20,79-83   
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta-TOF from center          
      REAL DT 
                                        !Exponential rise               
      REAL ALP 
                                        !Exponential decay              
      REAL BET 
                                        !Gaussian variance              
      REAL SIG 
                                        !Lorentzian FWHM                
      REAL GAM 
                                        !Value of pseudo-Voigt at DX    
      REAL FUNC 
                                        !dF/dta                         
      REAL DFDX 
                                        !dF/da                          
      REAL DFDA 
                                        !dF/db                          
      REAL DFDB 
                                        !dF/ds                          
      REAL DFDS 
                                        !dF/dg                          
      REAL DFDG 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
                                        !Linear combination coeffs      
!     REAL COFG (6), COFL (6) 
!     REAL ACOFG (7), ACOFL (7) 
                                        !Gaussian Normalization constant
!     REAL GNORM 
      REAL COFT (6), COFN (3) 
                                        !Normalization constant         
      REAL NORM 
                                        !Exp-integral arguements        
      REAL RXA, IXA, RXB, IXB 
                                        !Exp-integral results           
      REAL RFA, IFA, RFB, IFB 
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
                                        !2*SQRT(2LN2)                   
      DATA STOFW / 2.35482005 / 
                                         !SQRT(2PI)                     
      DATA SQ2PI / 2.506628275 / 
                                      !PI                               
      DATA PI / 3.141592654 / 
      DATA COFT / 1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0 / 
      DATA COFN / 1.36603, - 0.47719, 0.11116 / 
                                                                        
!CODE:                                                                  
                                                                        
      SQSGT = MAX (SQRT (SIG), 0.01) 
      GAM = MAX (GAM, 0.1) 
      FWHG = STOFW * SQSGT 
      PGL = FWHG**5 
      SUMHM = PGL 
      DSDL = 0.0 
      DSDG = 0.0 
      DO ITRM = 1, 5 
      PGL = PGL / FWHG 
      DSDL = DSDL + REAL(ITRM) * COFT (ITRM + 1) * PGL 
      DSDG = DSDG + REAL(6 - ITRM) * COFT (ITRM) * PGL 
      PGL = PGL * GAM 
      SUMHM = SUMHM + COFT (ITRM + 1) * PGL 
      ENDDO 
      FWHM = EXP (0.2 * LOG (SUMHM) ) 
      FRAC = GAM / FWHM 
      DEDF = 0.0 
      PF = 1.0 
      ETA = 0.0 
      DO ITRM = 1, 3 
      DEDF = DEDF + REAL(ITRM) * COFN (ITRM) * PF 
      PF = PF * FRAC 
      ETA = ETA + COFN (ITRM) * PF 
      ENDDO 
                                                                        
      NORM = 0.5 * ALP * BET / (ALP + BET) 
      SIGP = (FWHM / STOFW) **2 
      TJ1 = ALP * SIGP 
      TJ2 = BET * SIGP 
      SQSG = SQRT (2.0 * SIGP) 
      S2SG = 2.0 * SQSG * SIGP 
      Y1 = (TJ1 + DT) / SQSG 
      DY1DS = (TJ1 - DT) / S2SG 
      Y2 = (TJ2 - DT) / SQSG 
      DY2DS = (TJ2 + DT) / S2SG 
      EX1 = - 0.5 * (DT**2) / SIGP 
      DEXDS = - EX1 / SIGP 
      H1 = HFUNC (EX1, Y1, 0, DH1DY) 
      H2 = HFUNC (EX1, Y2, 0, DH2DY) 
                                                                        
      TG = NORM * (H1 + H2) 
                                                                        
      DTGDA = 2.0 * NORM * TG / (ALP * ALP) 
                                                         !OK            
      DTGDA = DTGDA + 0.5 * NORM * SQSG * DH1DY 
      DTGDB = 2.0 * NORM * TG / (BET * BET) 
                                                         !OK            
      DTGDB = DTGDB + 0.5 * NORM * SQSG * DH2DY 
      DTGDFW = TG * DEXDS + NORM * (DH1DY * DY1DS + DH2DY * DY2DS) 
                                                         !OK            
      DTGDFW = 2.0 * DTGDFW * FWHM / STOFW**2 
                                                         !OK            
      DTGDT = - NORM * (ALP * H1 - BET * H2) 
                                                                        
      RXA = - ALP * DT 
      RXB = - BET * DT 
      IXA = ALP * FWHM / 2.0 
      IXB = BET * FWHM / 2.0 
      CALL EXPINT (RXA, IXA, RFA, IFA) 
      CALL EXPINT (RXB, IXB, RFB, IFB) 
      TL = - 2.0 * NORM * (IFA + IFB) / PI 
      DIVSOR = DT**2 + FWHM**2 / 4.0 
                                                          !OK           
      DTLDT = - 2.0 * NORM * (ALP * IFA + BET * IFB + FWHM / DIVSOR)    &
      / PI                                                              
                                                                        
      DTLDFW = NORM * ( - ALP * RFA - BET * RFB - 2.0 * DT / DIVSOR)    &
      / PI                                                              
                                                                        
      DTLDA = 2.0 * NORM * TL / ALP**2 
                                                        !OK             
      DTLDA = DTLDA + 2.0 * NORM * (DT * IFA - FWHM * RFA / 2.0)        &
      / PI                                                              
      DTLDB = 2.0 * NORM * TL / BET**2 
                                                        !OK             
      DTLDB = DTLDB + 2.0 * NORM * (DT * IFB - FWHM * RFB / 2.0)        &
      / PI                                                              
                                                                        
                                                        !OK             
      FUNC = ETA * TL + (1.0 - ETA) * TG 
                                                                        
      DFDX = ETA * DTLDT + (1.0 - ETA) * DTGDT 
                                                                        
      DFWDG = 0.2 * DSDG * FWHM / SUMHM 
      DFRDG = - FRAC * DFWDG / FWHM 
                                                                        
      DFDS = DEDF * DFRDG * (TL - TG) 
      DFDS = DFDS + (ETA * DTLDFW + (1.0 - ETA) * DTGDFW) * DFWDG 
                                                         !OK?           
      DFDS = 0.5 * DFDS * STOFW / SQSGT 
                                                                        
      DFWDL = 0.2 * DSDL * FWHM / SUMHM 
      DFRDL = (1.0 - FRAC * DFWDL) / FWHM 
                                                                        
      DFDG = DEDF * DFRDL * (TL - TG) 
                                                        !WRONG          
      DFDG = DFDG + (ETA * DTLDFW + (1.0 - ETA) * DTGDFW) * DFWDL 
                                                                        
                                                        !OK             
      DFDA = ETA * DTLDA + (1.0 - ETA) * DTGDA 
                                                                        
                                                        !OK             
      DFDB = ETA * DTLDB + (1.0 - ETA) * DTGDB 
                                                                        
      RETURN 
      END SUBROUTINE EPSVOIGT                       
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE SIMGAUS1 (DTT, PCOT, SIG, PRFUNC, DPRDT, PCOTPRT,      &
      SIGPRT)                                                           
                                                                        
!PURPOSE:                                                               
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta 2-theta                  
      REAL DTT 
                                        !Asymm coeff.                   
      REAL PCOT 
                                        !Gaussian coeff.                
      REAL SIG 
                                        !Value of function at DTT       
      REAL PRFUNC 
                                        !Position derivative            
      REAL DPRDT 
                                        !Asymm. derivative              
      REAL PCOTPRT 
                                        !Gauss. derivative              
      REAL SIGPRT 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
                                        !Normalization                  
      REAL NORM 
                                        !Simpson's rule coefficients    
      REAL SIMC (7, 3) 
                                        !Simpson's rule multipliers     
      REAL SIMK (7, 3) 
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA SIMC / 0.0, 0.25, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0625, 0.25,&
      0.5625, 1.0, 0.0, 0.0, 0.0, 0.0277778, 0.1111111, 0.25, 0.4444444,&
      0.6944444, 1.0 /                                                  
      DATA SIMK / 1.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 2.0, 4.0,&
      1.0, 0.0, 0.0, 1.0, 4.0, 2.0, 4.0, 2.0, 4.0, 1.0 /                
      DATA SQ2PI / 2.506628275 / 
                                                                        
!CODE:                                                                  
                                                                        
      SQSG = SQRT (SIG) 
      IF (ABS (PCOT / SQSG) .GT.3.5) THEN 
         NT = 7 
         KT = 3 
         NORM = 1.0 / (18.0 * SQ2PI * SQSG) 
      ELSEIF (ABS (PCOT / SQSG) .GT.2.5) THEN 
         NT = 5 
         KT = 2 
         NORM = 1.0 / (12.0 * SQ2PI * SQSG) 
      ELSE 
         NT = 3 
         KT = 1 
         NORM = 1.0 / (6.0 * SQ2PI * SQSG) 
      ENDIF 
      PRFUNC = 0.0 
      SIGPRT = 0.0 
      DPRDT = 0.0 
      PCOTPRT = 0.0 
      DO IT = 1, NT 
      DT = DTT + SIMC (IT, KT) * PCOT 
      T = 0.0 
      EX = - 0.5 * DT**2 / SIG 
      IF (EX.GE. - 50.0) T = EXP (EX) 
      PRFUNC = PRFUNC + SIMK (IT, KT) * T 
      SIGPRT = SIGPRT + SIMK (IT, KT) * T * EX 
      DPRDT = DPRDT + SIMK (IT, KT) * T * DT 
      IF (IT.NE.1) PCOTPRT = PCOTPRT + SIMK (IT, KT) * SIMC (IT, KT)    &
      * T * DT                                                          
      ENDDO 
                                                                        
      PRFUNC = NORM * PRFUNC 
                                                                        
      SIGPRT = - PRFUNC / (2.0 * SIG) - NORM * SIGPRT / SIG 
      DPRDT = NORM * DPRDT / SIG 
      PCOTPRT = - NORM * PCOTPRT / SIG 
      RETURN 
      END SUBROUTINE SIMGAUS1                       
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE PSCVGT1 (DTT, PCOT, SIG, GAM, PRFUNC, DPRDT, PCOTPRT,  &
      SIGPRT, GAMPRT)                                                   
                                                                        
!PURPOSE:                                                               
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta 2-theta                  
      REAL DTT 
                                        !Asymmetry coeff.               
      REAL PCOT 
                                        !Gauss & Lorentz. coeffs.       
      REAL SIG, GAM 
                                        !Value of function at DTT       
      REAL PRFUNC 
                                        !Position derivative            
      REAL DPRDT 
                                        !Asymm. derivative              
      REAL PCOTPRT 
                                        !SIG & GAM derivatives          
      REAL SIGPRT, GAMPRT 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
                                        !Normalization constant         
      REAL NORM 
                                        !Full hw of max.                
      REAL FWHM 
                                        !No. Simpson's rule terms       
      INTEGER KT, NT 
                                        !Float of NT                    
      REAL TNT 
                                        !Simpson's rule coeff.          
      REAL SIMC, SIMK 
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA STOFW / 2.35482005 / 
                                                                        
!CODE:                                                                  
                                                                        
      FWHM = STOFW * SQRT (SIG) + GAM 
      KT = MAX (1, IFIX (2.5 * ABS (PCOT / FWHM) ) ) 
                                                                        
      NT = 2 * KT + 1 
      NORM = 1.0 / (6.0 * REAL(KT) ) 
      PRFUNC = 0.0 
      SIGPRT = 0.0 
      GAMPRT = 0.0 
      DPRDT = 0.0 
      PCOTPRT = 0.0 
      TNT = (NT - 1) **2 
      DO IT = 1, NT 
      IF (MOD (IT, 2) .EQ.0) THEN 
         SIMK = 4.0 
      ELSEIF (IT.EQ.1.OR.IT.EQ.NT) THEN 
         SIMK = 1.0 
      ELSE 
         SIMK = 2.0 
      ENDIF 
      SIMC = (IT - 1) **2 / TNT 
      DT = DTT + PCOT * SIMC 
      CALL PSVOIGT (DT, SIG, GAM, FUNC, DFDT, DFDS, DFDG) 
      PRFUNC = PRFUNC + SIMK * FUNC 
      SIGPRT = SIGPRT + SIMK * DFDS 
      GAMPRT = GAMPRT + SIMK * DFDG 
      DFDT = SIMK * DFDT 
      DPRDT = DPRDT + DFDT 
      PCOTPRT = PCOTPRT + SIMC * DFDT 
      ENDDO 
                                                                        
      PRFUNC = NORM * PRFUNC 
      SIGPRT = NORM * SIGPRT 
      GAMPRT = NORM * GAMPRT 
      DPRDT = - NORM * DPRDT 
      PCOTPRT = NORM * PCOTPRT 
      RETURN 
      END SUBROUTINE PSCVGT1                        
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE LORENTZ (DT, GAM, FUNC, DLDT, DLDG) 
                                                                        
!PURPOSE:Calculate Lorentzian function & derivatives                    
                                                                        
!PSEUDOCODE:                                                            
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
                                                                        
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta                          
      REAL DT 
                                        !Coefficient                    
      REAL GAM 
                                        !Function                       
      REAL FUNC 
                                        !df/dt                          
      REAL DLDT 
                                        !df/dg                          
      REAL DLDG 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
REAL :: PI, FNORM, DENOM
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA PI / 3.1415292654 / 
!CODE:                                                                  
                                                                        
      FNORM = 2.0 * GAM / PI 
      DENOM = GAM**2 + 4.0 * DT**2 
      IF (ABS (GAM) .GT.5.0E-6) THEN 
         FUNC = FNORM / DENOM 
         DLDT = - 8.0 * DT * FUNC / DENOM 
         DLDG = (8.0 * DT**2 - 2.0 * GAM**2) / (PI * DENOM**2) 
      ELSE 
         FUNC = 0.0 
         DLDT = 0.0 
         IF (ABS (DT) .LT.1.0E-5) THEN 
            DLDG = 0.0 
         ELSE 
            DLDG = 2.0 / (PI * DT**2) 
         ENDIF 
      ENDIF 
      RETURN 
      END SUBROUTINE LORENTZ                        
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE PSVFCJ (DTT, TTHETA, SL, HL, SIG, GAM, PRFUNC, DPRDT,  &
      SLPART, HLPART, SIGPART, GAMPART)                                 
                                                                        
!PURPOSE: Compute function & derivatives for Pseudovoigt profile        
!   [W.I.F.David (1986), J. Appl. Cryst. 19, 63-64 &                    
!    P. Thompson, D.E. Cox & J.B. Hastings (1987)                       
!                                         J. Appl. Cryst.,20,79-83.]    
! Finger-Cox-Jephcoat (FCJ94) asymmetry correction                      
!   [L.W. Finger, D.E. Cox & A.P. Jephcoat (1994)                       
!                                         J. Appl. Cryst.,27,892-900.]  
! coded 11/95 by B. H. Toby (NIST). revised version                     
! parameterized as asym1=S/L asym2=H/L                                  
                                                                        
      USE trig_degree_mod
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
      REAL DTT 
      REAL TTHETA 
      REAL SL, HL 
      REAL SIG, GAM 
      REAL PRFUNC 
      REAL DPRDT 
      REAL SLPART, HLPART 
      REAL SIGPART, GAMPART 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
                                        ! pseudo-voight intensity       
      REAL R 
                                        ! deriv R w/r theta             
      REAL DRDT 
                                        ! deriv R w/r sig               
      REAL DRDS 
                                        ! deriv R w/r gam               
      REAL DRDG 
      REAL F 
      REAL G 
      REAL DFDA 
      REAL DGDA 
      REAL DGDB 
      REAL DYDA 
      REAL DYDB 
                                        ! cos(2theta)**2                
      REAL COS2THETA2 
                                        ! cos(2theta)                   
      REAL COS2THETA 
                                        ! sin(2THETA)                   
      REAL SIN2THETA 
                                        ! sin(Delta)                    
      REAL SINDELTA 
                                        ! cos(Delta)                    
      REAL COSDELTA 
                                        ! 1/cos(Delta)                  
      REAL RCOSDELTA 
                                        ! tan(Delta)                    
      REAL TANDELTA 
                                        ! cos(Delta)**2                 
      REAL COSDELTA2 
                                        ! asym1 [coff(7)]               
      REAL A 
                                        ! asym2 [coff(8)]               
      REAL B 
                                        ! (A+B)                         
      REAL APB 
                                        ! (A-B)                         
      REAL AMB 
                                        ! (A+B)**2                      
      REAL APB2 
                                                                        
! Intermediate variables                                                
                                                                        
                                        !      (sum of w G)**2          
      REAL SUMWG2 
                                        !      sum of w G               
      REAL SUMWG 
                                        !      sum of w G               
      REAL SUMWRG 
                                        !      sum of w dGdA            
      REAL SUMWDGDA 
                                        !      sum of w R dGdA          
      REAL SUMWRDGDA 
                                        !      sum of w dGdB            
      REAL SUMWDGDB 
                                        !      sum of w R dGdB          
      REAL SUMWRDGDB 
                                        !      sum of w G dRd(2theta)   
      REAL SUMWGDRD2T 
                                        !      sum of w G dRdp(n)       
      REAL SUMWGDRDSIG 
                                        !      sum of w G dRdp(n)       
      REAL SUMWGDRDGAM 
      REAL SUMWGDRDA 
      REAL SUMWGDRDB 
                                        ! 2phi minimum                  
      REAL EMIN 
                                        ! 2phi of inflection point      
      REAL EINFL 
                                        ! Derivative of Emin wrt A      
      REAL DEMINDA 
                                        ! Angle of integration for      
      REAL DELTA 
!                                         convolution                   
      REAL DDELTADA 
                                        ! intermediates                 
      REAL TMP, TMP1, TMP2 
                                        ! Miscellaneous loop variables  
      INTEGER K, IT 
!                                                                       
!       Local Variables for Gaussian Integration                        
!                                                                       
                                        !number of terms in Gaussian q  
      INTEGER NGT 
!                                         uadrature                     
                                        ! number of pre-computed        
      INTEGER NUMTAB 
!                                         Gaussian tables               
      PARAMETER (numtab = 14) 
                                        ! number of terms in each       
      INTEGER NTERMS (NUMTAB) 
!                                         table - must be even          
                                        ! location of 1st term:         
      INTEGER FSTTERM (NUMTAB) 
!                                         N.B. N/2 terms                
                                        ! true if table has             
      LOGICAL CALCLFG (NUMTAB) 
!                                         previously been calculated    
                                        ! number of selected array      
      INTEGER ARRAYNUM 
                                        ! size of complete array        
      INTEGER ARRAYSIZE 
      PARAMETER (arraysize = 1883) 
                                        !Gaussian abscissas             
      REAL XP (ARRAYSIZE) 
                                        !Gaussian weights               
      REAL WP (ARRAYSIZE) 
                                        !temporary Gaussian abscissas   
      REAL XPT (1000) 
                                        !temporary Gaussian weights     
      REAL WPT (1000) 
!     REAL STOFW 
!     PARAMETER (STOFW = 2.35482005) 
      REAL TODEG 
      PARAMETER (todeg = 57.2957795) 
                                       !Values to be saved across calls 
      SAVE calclfg, xp, wp 
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
                                      ! number of terms in each table   
      DATA NTERMS / 6, 10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 600,&
      800, 1000 /                                                       
!                                         - must be even                
!       note that nterms determines both arraysize and fstterm          
                                       ! loc. of 1st term:              
      DATA FSTTERM / 0, 3, 8, 18, 38, 68, 108, 158, 233, 333, 483, 683, &
      983, 1383 /                                                       
!                                         N.B. N/2 terms are saved      
                                        ! true if table entry           
      DATA calclfg / numtab * .false. / 
!                                         has been calculated           
                                                                        
!CODE:                                                                  
                                                                        
!                                                                       
! f(2theta) intermediates                                               
!                                                                       
      SIN2THETA = SIND (TTHETA / 100.0) 
      cos2THETA = COSD (TTHETA / 100.0) 
      cos2THETA2 = cos2THETA * cos2THETA 
!                                                                       
! Asymmetry terms                                                       
!                                                                       
                        ! A = S/L in FCJ                                
      A = SL 
                        ! B = H/L in FCJ                                
      B = HL 
      ApB = A + B 
      AmB = A - B 
      ApB2 = ApB * ApB 
!                                                                       
! handle the case where there is asymmetry                              
!                                                                       
      IF (A.ne.0.0.or.B.ne.0.0) then 
         Einfl = Acosd (SQRT (1.0 + AmB**2) * cos2THETA) 
! 2phi(infl) FCJ eq 5 (degrees)                                         
         tmp2 = 1.0 + ApB2 
         tmp = SQRT (tmp2) * cos2THETA 
!                                                                       
! Treat case where A or B is zero - Set Einfl = 2theta                  
!                                                                       
         IF (A.eq.0.0.or.B.eq.0.0) Einfl = Acosd (cos2THETA) 
         IF (abs (tmp) .le.1.0) then 
                                 ! 2phi(min) FCJ eq 4 (degrees)         
            Emin = Acosd (tmp) 
            tmp1 = tmp2 * (1.0 - tmp2 * cos2THETA2) 
         ELSE 
            tmp1 = 0.0 
            IF (tmp.gt.0.0) then 
               Emin = 0.0 
            ELSE 
               Emin = 180.0 
            ENDIF 
         ENDIF 
         IF (tmp1.gt.0.and.abs (tmp) .le.1.0) then 
                                              ! N. B. symm w/r A,B      
            dEmindA = - ApB * cos2THETA / SQRT (tmp1) 
         ELSE 
            dEmindA = 0.0 
         ENDIF 
!                                                                       
! Compute Gaussian Quadrature interval                                  
! Note that the number of points must be large enough so that the       
! interval between 2phi(min) and 2theta must be divided into steps no   
! larger than 0.005 degrees.  LWF  8/10/95                              
!                                                                       
! Determine which Gauss-Legendre Table to use                           
!                                                                       
         arraynum = 1 
!                                                                       
! Calculate desired number of intervals                                 
!                                                                       
         tmp = abs (TTHETA - emin * 100.) 
         IF (GAM.LE.0.0) THEN 
            K = INT (TMP * 8.0) 
         ELSE 
            k = max (int (tmp * 8.0), int (0.2 * tmp / GAM) ) 
         ENDIF 
!                                                                       
! Find the next largest set                                             
!                                                                       
         DO while (arraynum.lt.numtab.and.k.gt.nterms (arraynum) ) 
         arraynum = arraynum + 1 
         ENDDO 
         ngt = nterms (arraynum) 
!                                                                       
! calculate the terms, if they have not been used before                
!                                                                       
         IF (.not.calclfg (arraynum) ) then 
            calclfg (arraynum) = .true. 
            CALL gauleg ( - 1., 1., xpt, wpt, ngt) 
            it = fstterm (arraynum) - ngt / 2 
!                                                                       
! copy the terms from our working array to the stored array             
!                                                                       
            DO k = ngt / 2 + 1, ngt 
            xp (k + it) = xpt (k) 
            wp (k + it) = wpt (k) 
            ENDDO 
         ENDIF 
         sumWG2 = 0. 
         sumWG = 0. 
         sumWRG = 0. 
         sumWdGdA = 0. 
         sumWRdGdA = 0. 
         sumWdGdB = 0. 
         sumWRdGdB = 0. 
         sumWGdRd2t = 0. 
         sumWGdRdsig = 0. 
         sumWGdRdgam = 0. 
         sumWGdRdA = 0. 
         sumWGdRdB = 0. 
!                                                                       
! Compute Convolution integral for 2phi(min) <= delta <= 2theta         
!                                                                       
         it = fstterm (arraynum) - ngt / 2 
         DO k = ngt / 2 + 1, ngt 
         delta = emin + (TTHETA / 100. - emin) * xp (k + it) 
! Delta in degrees                                                      
                                             ! N. B. symm w/r A,B       
         dDELTAdA = (1. - xp (k + it) ) * dEmindA 
         sinDELTA = sind (Delta) 
         cosDELTA = cosd (Delta) 
         IF (abs (cosDELTA) .lt.1.0e-15) cosDELTA = 1.0e-15 
         RcosDELTA = 1. / cosDELTA 
         tanDELTA = tand (Delta) 
                                        ! cos(Delta)**2                 
         cosDELTA2 = cosDELTA * cosDELTA 
         tmp = cosDELTA2 - cos2THETA2 
         IF (tmp.gt.0) then 
            tmp1 = SQRT (tmp) 
            F = abs (cos2THETA) / tmp1 
            dFdA = cosDELTA * cos2THETA * sinDELTA * dDELTAdA / (tmp1 * &
            tmp1 * tmp1)                                                
         ELSE 
            F = 0.0 
            dFdA = 0.0 
         ENDIF 
!                                                                       
! Calculate G(Delta,2theta) [G = W /(h cos(delta) ]                     
! [ FCJ eq. 7(a) and 7(b) ]                                             
!                                                                       
         IF (abs (delta - emin) .gt.abs (einfl - emin) ) then 
            IF (A.ge.B) then 
!                                                                       
! N.B. this is the only place where d()/dA <> d()/dB                    
!                                                                       
               G = 2.0 * B * F * RcosDELTA 
               dGdA = 2.0 * B * RcosDELTA * (dFdA + F * tanDELTA *      &
               dDELTAdA)                                                
               dGdB = dGdA + 2.0 * F * RcosDELTA 
            ELSE 
               G = 2.0 * A * F * RcosDELTA 
               dGdB = 2.0 * A * RcosDELTA * (dFdA + F * tanDELTA *      &
               dDELTAdA)                                                
               dGdA = dGdB + 2.0 * F * RcosDELTA 
            ENDIF 
                                ! delta .le. einfl .or. min(A,B) .eq. 0 
         ELSE 
            G = ( - 1.0 + ApB * F) * RcosDELTA 
            dGdA = RcosDELTA * (F - tanDELTA * dDELTAdA + ApB * F *     &
            tanDELTA * dDELTAdA + ApB * dFdA)                           
            dGdB = dGdA 
         ENDIF 
         CALL PSVOIGT (DTT + TTHETA - delta * 100., SIG, GAM, R, dRdT,  &
         dRdS, dRdG)                                                    
         sumWG = sumWG + wp (k + it) * G 
         sumWRG = sumWRG + wp (k + it) * R * G 
         sumWdGdA = sumWdGdA + wp (k + it) * dGdA 
         sumWRdGdA = sumWRdGdA + wp (k + it) * R * dGdA 
         sumWdGdB = sumWdGdB + wp (k + it) * dGdB 
         sumWRdGdB = sumWRdGdB + wp (k + it) * R * dGdB 
         sumWGdRd2t = sumWGdRd2t + wp (k + it) * G * dRdT 
!                     N.B. 1/centidegrees                               
         sumWGdRdsig = sumWGdRdsig + wp (k + it) * G * dRdS 
         sumWGdRdgam = sumWGdRdgam + wp (k + it) * G * dRdG 
         sumWGdRdA = sumWGdRdA + wp (k + it) * G * dRdT * dDELTAdA 
!                     N. B. symm w/r A,B                                
         ENDDO 
         sumWG2 = sumWG * sumWG 
         PRFUNC = sumWRG / sumWG 
         dydA = ( - (sumWRG * sumWdGdA) + sumWG * (sumWRdGdA - 100.0 *  &
         todeg * sumWGdRdA) ) / sumWG2                                  
         dydB = ( - (sumWRG * sumWdGdB) + sumWG * (sumWRdGdB - 100.0 *  &
         todeg * sumWGdRdA) ) / sumWG2                                  
         sigpart = sumWGdRdsig / sumWG 
         gampart = sumWGdRdgam / sumWG 
         DPRDT = - SumWGdRd2T / sumWG 
      ELSE 
!                                                                       
! no asymmetry -- nice and simple!                                      
!                                                                       
         CALL PSVOIGT (DTT, SIG, GAM, R, dRdT, dRdS, dRdG) 
         PRFUNC = R 
         dydA = 0.002 * sign (1.0, TTHETA - DTT) 
         dydB = dydA 
         sigpart = dRdS 
         gampart = dRdG 
         DPRDT = - dRdT 
      ENDIF 
      SLPART = DYDA 
      HLPART = DYDB 
                                                                        
      RETURN 
      END SUBROUTINE PSVFCJ                         
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE PSVOIGT (DX, SIG, GAM, FUNC, DFDX, DFDS, DFDG) 
                                                                        
!PURPOSE: Compute function & derivatives pseudovoigt                    
!pseudo Voigt P.Tompson, D.E. Cox & J.B. Hastings (1987).               
!J. Appl. Cryst.,20,79-83                                               
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta-x from center            
      REAL DX 
                                        !Gaussian variance              
      REAL SIG 
                                        !Lorentzian FWHM                
      REAL GAM 
                                        !Value of pseudo-Voigt at DX    
      REAL FUNC 
                                        !dF/dx                          
      REAL DFDX 
                                        !dF/ds                          
      REAL DFDS 
                                        !dF/dg                          
      REAL DFDG 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES
INTEGER :: ITRM, ITYPE                                                       
REAL DEDF, DFRDG, DFWDG, DFWDL, DSDG, DSDL
REAL :: DTLDFW, DTLDT, ETA, EX, FRAC, FWHG, FWHM, PF, PGL
REAL :: SIGP, SQSG, SQ2PI, STOFW, SUMHM, TG, TL, TS
                                                                        
                                        !Linear combination coeffs      
                                        !Gaussian Normalization constant
      REAL COFT (6), COFN (3) 
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
                                        !2*SQRT(2LN2)                   
      DATA STOFW / 2.35482005 / 
                                         !SQRT(2PI)                     
      DATA SQ2PI / 2.506628275 / 
      DATA COFT / 1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0 / 
      DATA COFN / 1.36603, - 0.47719, 0.11116 / 
      DATA ITYPE / 1 / 
                                                                        
!CODE:                                                                  
                                                                        
      SQSG = MAX (SQRT (SIG), 0.001) 
      FWHG = STOFW * SQSG 
      PGL = FWHG**5 
      SUMHM = PGL 
      DSDL = 0.0 
      DSDG = 0.0 
      DO ITRM = 1, 5 
      PGL = PGL / FWHG 
      DSDL = DSDL + REAL(ITRM) * COFT (ITRM + 1) * PGL 
      DSDG = DSDG + REAL(6 - ITRM) * COFT (ITRM) * PGL 
      PGL = PGL * GAM 
      SUMHM = SUMHM + COFT (ITRM + 1) * PGL 
      ENDDO 
      FWHM = EXP (0.2 * LOG (SUMHM) ) 
      FRAC = GAM / FWHM 
      DEDF = 0.0 
      PF = 1.0 
      ETA = 0.0 
      DO ITRM = 1, 3 
      DEDF = DEDF + REAL(ITRM) * COFN (ITRM) * PF 
      PF = PF * FRAC 
      ETA = ETA + COFN (ITRM) * PF 
      ENDDO 
      CALL LORENTZ (DX, FWHM, TL, DTLDT, DTLDFW) 
      SIGP = (FWHM / STOFW) **2 
      EX = MAX ( - 20.0, - 0.5 * DX**2 / SIGP) 
      TG = STOFW * EXP (EX) / (SQ2PI * FWHM) 
      FUNC = ETA * TL + (1.0 - ETA) * TG 
                                                                        
      TS = - 2.0 * (1.0 - ETA) * TG * (EX + 0.5) / FWHM 
      DFDX = ETA * DTLDT - (1.0 - ETA) * TG * DX / SIGP 
                                                                        
      DFWDG = 0.2 * DSDG * FWHM / SUMHM 
      DFRDG = - FRAC * DFWDG / FWHM 
      DFDS = DEDF * DFRDG * (TL - TG) + (ETA * DTLDFW + TS) * DFWDG 
      DFDS = 0.5 * DFDS * STOFW / SQSG 
                                                                        
      DFWDL = 0.2 * DSDL * FWHM / SUMHM 
      DFRDL = (1.0 - FRAC * DFWDL) / FWHM 
      DFDG = DEDF * DFRDL * (TL - TG) + (ETA * DTLDFW + TS) * DFWDL 
                                                                        
      RETURN 
      END SUBROUTINE PSVOIGT                        
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE gauleg (x1, x2, x, w, n) 
! Routine from Numerical Recipes (Press, Flannery,                      
!  Teukolsky and Vetterling,                                            
!    1986, Cambridge University Press, ISBN 0 521 30811 9)              
!                                                                       
! Given the lower and upper limits of integration (X1, X2)              
!  and the number                                                       
!  of intervals (N), this routine returns arrays X and W of length N,   
!  containing the abscissas and weights of the Gauss-Legendre N-point   
!  quadrature formula.                                                  
!                                                                       
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
      REAL X1, X2, X (N), W (N) 
      PARAMETER (eps = 3.e-7) 
!                                                                       
      pp = 1.0
      m = (n + 1) / 2 
      xm = 0.5 * (x2 + x1) 
      xl = 0.5 * (x2 - x1) 
      DO i = 1, m 
      z = cos (3.141592654 * (i - .25) / (n + .5) ) 
      z1 = 0.0 
      DO while (abs (z - z1) .gt.eps) 
      p1 = 1.0 
      p2 = 0.0 
      DO j = 1, n 
      p3 = p2 
      p2 = p1 
      p1 = ( (2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j 
      ENDDO 
      pp = n * (z * p1 - p2) / (z * z - 1.0) 
      z1 = z 
      z = z1 - p1 / pp 
      ENDDO 
      x (i) = xm - xl * z 
      x (n + 1 - i) = xm + xl * z 
      w (i) = 2.0 * xl / ( (1.0 - z * z) * pp * pp) 
      w (n + 1 - i) = w (i) 
      ENDDO 
      RETURN 
      END SUBROUTINE gauleg                         
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE EDSGAUS1 (DT, SGSQ, PRFUNC, DPRDT, SIGPRT) 
                                                                        
!PURPOSE:                                                               
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Delta KEV                      
      REAL DT 
                                        !Gauss. coeff.                  
      REAL SGSQ 
                                        !Value of function              
      REAL PRFUNC 
                                        !Position derivative            
      REAL DPRDT 
                                        !Gauss. derivative              
      REAL SIGPRT 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
REAL :: EX, GNORM, SIG, SQ2PI, T
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA SQ2PI / 2.506628275 / 
                                                                        
!CODE:                                                                  
      SIG = SQRT (SGSQ) 
      GNORM = 1.0 / (SQ2PI * SIG) 
      EX = - 0.5 * DT**2 / SGSQ 
      T = 0.0 
      IF (EX.GE. - 50.0) T = EXP (EX) 
      PRFUNC = GNORM * T 
      DPRDT = DT * PRFUNC / SGSQ 
                                                                        
      SIGPRT = - PRFUNC * (1 + 2.0 * EX) / (2.0 * SGSQ) 
                                                                        
      RETURN 
      END SUBROUTINE EDSGAUS1                       
                                                                        
!-----------------------------------------------------------------------
                                                                        
      real FUNCTION HFUNC (X, Y, IFLAG, DHDY) 
                                                                        
!PURPOSE: Compute exp(x)exp(y**2)erfc(y) and partial derivatives        
                                                                        
!PSEUDOCODE:                                                            
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !Arguments                      
      REAL X, Y 
                                        !Control for calculation method 
      INTEGER IFLAG 
                                        !Partial derivative wrt y       
      REAL DHDY 
                                        !Value of function              
!     REAL HFUNC 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
REAL :: T, T1, CONST, CONS2
                                                                        
!SUBROUTINES CALLED:                                                    
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
                                        !Complementary error function   
!     REAL GERFC 
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
                                                  !2/SQRT(PI),4/PI      
      DATA CONST / 1.128379167 /, CONS2 / 1.273239545 / 
                                                                        
!CODE:                                                                  
                                                                        
      T = Y**2 + X 
      HFUNC = 0.0 
      DHDY = 0.0 
                                               !Compute function        
      IF (Y.GE.6.000.OR.IFLAG.EQ.1) THEN 
!                                   - protect against underflows        
         T1 = SQRT (Y**2 + CONS2) 
                                               !Upper inequality        
         IF (X.GT. - 75.0) HFUNC = EXP (X) 
         HFUNC = HFUNC * CONST / (Y + T1) 
         DHDY = - HFUNC / T1 
      ELSEIF (Y.LE. - 6.375) THEN 
                                               !Upper inequality        
         IF (T.GT. - 75.0) HFUNC = 2.0 * EXP (T) 
         DHDY = 2.0 * Y * HFUNC 
      ELSE 
         IF (T.GT. - 75.0) HFUNC = EXP (T) 
         HFUNC = HFUNC * GERFC (Y) 
         IF (X.GT. - 75.0) DHDY = EXP (X) 
         DHDY = 2.0 * Y * HFUNC - CONST * DHDY 
      ENDIF 
      RETURN 
      END FUNCTION HFUNC                            
                                                                        
!-----------------------------------------------------------------------
                                                                        
      real FUNCTION GERFC (Y) 
                                                                        
!PURPOSE: This routine returns the Complementary ERROR Function         
!         of a single precision variable (Y)                            
                                                                        
!CALLING ARGUMENTS:                                                     
                                                                        
                                        !                               
      REAL Y 
                                         !Result                        
!     REAL GERFC 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
      DIMENSION P (5), Q (4), P1 (9), Q1 (8), P2 (6), Q2 (5) 
      REAL P, Q, P1, Q1, P2, Q2, XMIN, XLARGE, SQRPI, X, RES, XSQ, XNUM,&
      XDEN, XI                                                          
      INTEGER ISW, I 
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
      DATA P (1) / 113.86415 /, P (2) / 377.48524 /, P (3) / 3209.3776 /&
      , P (4) / .18577771 /, P (5) / 3.1611237 /                        
      DATA Q (1) / 244.02464 /, Q (2) / 1282.6165 /, Q (3) / 2844.2368 /&
      , Q (4) / 23.601291 /                                             
      DATA P1 (1) / 8.8831498 /, P1 (2) / 66.119190 /, P1 (3) /         &
      298.63514 /, P1 (4) / 881.95222 /, P1 (5) / 1712.0476 /, P1 (6)   &
      / 2051.0784 /, P1 (7) / 1230.3394 /, P1 (8) / 2.1531154E-8 /,     &
      P1 (9) / .56418850 /                                              
      DATA Q1 (1) / 117.69395 /, Q1 (2) / 537.18110 /, Q1 (3) /         &
      1621.3896 /, Q1 (4) / 3290.7992 /, Q1 (5) / 4362.6191 /, Q1 (6)   &
      / 3439.3677 /, Q1 (7) / 1230.3394 /, Q1 (8) / 15.744926 /         
      DATA P2 (1) / - 3.6034490E-01 /, P2 (2) / - 1.2578173E-01 /,      &
      P2 (3) / - 1.6083785E-02 /, P2 (4) / - 6.5874916E-04 /, P2 (5)    &
      / - 1.6315387E-02 /, P2 (6) / - 3.0532663E-01 /                   
      DATA Q2 (1) / 1.8729528 /, Q2 (2) / 5.2790510E-01 /, Q2 (3)       &
      / 6.0518341E-02 /, Q2 (4) / 2.3352050E-03 /, Q2 (5) / 2.5685202 / 
      DATA XMIN / 1.0E-10 /, XLARGE / 6.375 / 
      DATA SQRPI / .56418958 / 
                                                                        
!CODE:                                                                  
                                                                        
      IF (Y.GT.0.0) THEN 
         X = Y 
         ISW = 1 
      ELSE 
         X = - Y 
         ISW = - 1 
      ENDIF 
      XSQ = X * X 
      IF (X.LT.0.477) THEN 
         IF (X.LT.XMIN) THEN 
            RES = X * P (3) / Q (3) 
         ELSE 
            XNUM = P (4) * XSQ + P (5) 
            XDEN = XSQ + Q (4) 
            DO I = 1, 3 
            XNUM = XNUM * XSQ + P (I) 
            XDEN = XDEN * XSQ + Q (I) 
            ENDDO 
            RES = X * XNUM / XDEN 
         ENDIF 
         IF (ISW.EQ. - 1) RES = - RES 
         RES = 1.0 - RES 
      ELSE 
         IF (X.LT.XLARGE) THEN 
            IF (X.LE.4.0) THEN 
               XNUM = P1 (8) * X + P1 (9) 
               XDEN = X + Q1 (8) 
               DO I = 1, 7 
               XNUM = XNUM * X + P1 (I) 
               XDEN = XDEN * X + Q1 (I) 
               ENDDO 
               RES = XNUM / XDEN 
            ELSE 
               XI = 1.0 / XSQ 
               XNUM = P2 (5) * XI + P2 (6) 
               XDEN = XI + Q2 (5) 
               DO I = 1, 4 
               XNUM = XNUM * XI + P2 (I) 
               XDEN = XDEN * XI + Q2 (I) 
               ENDDO 
               RES = (SQRPI + XI * XNUM / XDEN) / X 
            ENDIF 
            RES = RES * EXP ( - XSQ) 
         ELSE 
            RES = 0.0 
         ENDIF 
         IF (ISW.EQ. - 1) RES = 2.0 - RES 
      ENDIF 
      GERFC = RES 
      RETURN 
      END FUNCTION GERFC                            
                                                                        
!-----------------------------------------------------------------------
                                                                        
      SUBROUTINE CNVPTP1 (PTYPE, PCOF, NCOF, PPCOF, DSP, STHETA, MAXFIT) 
                                                                        
!PURPOSE: Convert TOF profile coeff. to single peak coeff.              
                                                                        
IMPLICIT INTEGER (i-n)
IMPLICIT REAL (a - h, o - z) 
!CALLING ARGUMENTS:                                                     
                                                                        
INTEGER, INTENT(IN) :: MAXFIT
                                        !Function type                  
      INTEGER PTYPE 
                                        !Input profile coeff            
      REAL PCOF (MAXFIT) 
                                        !No. single peak coeff.         
      INTEGER NCOF 
                                        !Output single peak coeff.      
      REAL PPCOF (MAXFIT) 
                                        !D-spacing                      
      REAL DSP 
                                        !Sin(theta) for bank            
      REAL STHETA 
                                                                        
!INCLUDE STATEMENTS:                                                    
                                                                        
!LOCAL VARIABLES:                                                       
                                                                        
                                        !Powers of d-spacing            
      REAL DSP2, DSP4 
                                        !Wavelength                     
      REAL ALAM 
                                                                        
!FUNCTION DEFINITIONS:                                                  
                                                                        
!DATA STATEMENTS:                                                       
                                                                        
!CODE:                                                                  
                                                                        
      DSP2 = DSP * DSP 
      DSP4 = DSP2 * DSP2 
      IF (PTYPE.EQ.1) THEN 
         NCOF = 5 
                                                       !ALP             
         PPCOF (3) = PCOF (1) + PCOF (2) / DSP 
                                                       !BET             
         PPCOF (4) = PCOF (3) + PCOF (4) / DSP4 
                                                       !SIG             
         PPCOF (5) = PCOF (5) + PCOF (6) * DSP2 + PCOF (7) * DSP4 
      ELSEIF (PTYPE.EQ.2) THEN 
         NCOF = 7 
         ALAM = ABS (2.0 * DSP * STHETA) 
                                                       !ALP             
         PPCOF (3) = 1.0 / (PCOF (1) + PCOF (2) * ALAM) 
                                                       !BET             
         PPCOF (4) = 1.0 / PCOF (3) 
                                                       !SWITCH -> RAT   
         PPCOF (5) = EXP ( - 81.799 / (PCOF (4) * ALAM**2) ) 
                                                       !SIG             
         PPCOF (6) = PCOF (5) + PCOF (6) * DSP2 + PCOF (7) * DSP4 
                                                       !GAM             
         PPCOF (7) = PCOF (8) + PCOF (9) * DSP + PCOF (10) * DSP2 
      ELSEIF (PTYPE.EQ.3) THEN 
         NCOF = 6 
                                                       !ALP             
         PPCOF (3) = PCOF (1) 
                                                       !BET             
         PPCOF (4) = PCOF (2) + PCOF (3) / DSP4 
                                                       !SIG             
         PPCOF (5) = PCOF (4) + PCOF (5) * DSP2 + PCOF (6) * DSP4 
                                                       !GAM             
         PPCOF (6) = PCOF (7) + PCOF (8) * DSP + PCOF (9) * DSP2 
      ELSEIF (PTYPE.EQ.4) THEN 
         NCOF = 6 
                                                       !ALP             
         PPCOF (3) = PCOF (1) 
                                                       !BET             
         PPCOF (4) = PCOF (2) + PCOF (3) / DSP4 
                                                       !SIG             
         PPCOF (5) = PCOF (4) * DSP2 + PCOF (5) * DSP4 
                                                       !GAM             
         PPCOF (6) = PCOF (6) * DSP2 
      ENDIF 
      RETURN 
      END SUBROUTINE CNVPTP1                        
!
end module kuplot_gsas_mod
