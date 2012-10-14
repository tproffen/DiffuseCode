!+
!     Include file for free variables
!
!     Warning!      MAXPAR_RES  in "param.inc" must always be of 
!                                 identical size to
!                 CHEM_MAX_NEIG in DISCUS/"config.inc"
!-
      INTEGER     MAXPAR
      INTEGER     MAXPAR_RES
!
      PARAMETER   (MAXPAR     = 500)
      PARAMETER   (MAXPAR_RES = 3000 )
!
      INTEGER     inpara(0:MAXPAR)
      REAL        rpara(0:MAXPAR)
      REAL        res_para(0:MAXPAR_RES)
!
      COMMON  /params/ inpara,rpara,res_para

