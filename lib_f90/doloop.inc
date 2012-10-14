!+
!     This file contains subroutines for: INCLUDE BLOCK fuers dkdo.f
!-
      INTEGER        MAXLEV,MAXCOM
      PARAMETER     (MAXLEV=20)
        PARAMETER     (MAXCOM=2000)
!
      INTEGER        nloop(0:MAXLEV),iloop(0:MAXLEV),do_kpara(1)
!
      LOGICAL        ldostart(0:MAXLEV)
!
      REAL           ghigh(0:MAXLEV),glow(0:MAXLEV)
      REAL           ginc(0:MAXLEV)
!
      COMMON /doloop/    iloop,nloop,glow,ghigh,ginc,do_kpara,ldostart
