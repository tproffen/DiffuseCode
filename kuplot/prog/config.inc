!#######################################################################
!     KUPLOT configuration file
!#######################################################################
!
!     MAXARRAY   : Defines the total number of points KUPL can read
!     MAXKURVTOT : Defines the maximum number of different data sets
!     MAXPARA    : Defines the maximum number of fitting PARAMETERs
!     MAXHL      : Defines the maximum number of contour line sets
!     MAXZ       : Defines maximum bitmap size bitmap(maxz,maxz)
!     MAXFRAME   : Defines the maximum number of frames
!     MAXWIN     : Defines the maximum number of windows 
!     MAXSP      : Defines the number of points used for spline/steps
!     MAXAN      : Defines the maximum number of annotations
!     MAXBOND    : Defines the maximum number of diff. bonds
!     MAXCOL     : Defines the maximum number of bitmap colours
!     MAXB       : Defines the maximum number of menu buttons
!     MAXGSAS    : Define maximum number of points in GSAS files
!
!#######################################################################
!
      INTEGER MAXPARA,MAXKURVTOT,MAXARRAY,MAXHL
      INTEGER MAXFRAME,MAXSP,MAXAN,MAXCOL,MAXWIN
      INTEGER MAXBOND,MAXZ,MAXB,MAXGSAS
!
      PARAMETER      (MAXARRAY   = 1500000)
      PARAMETER      (MAXKURVTOT = 100)
      PARAMETER      (MAXPARA    = 150)
      PARAMETER      (MAXHL      = 5)
      PARAMETER      (MAXZ       = 501)
      PARAMETER      (MAXFRAME   = 16)
      PARAMETER      (MAXWIN     = 3)
      PARAMETER      (MAXSP      = 500)
      PARAMETER      (MAXAN      = 500)
      PARAMETER      (MAXBOND    = 50)
      PARAMETER      (MAXCOL     = 255)
      PARAMETER      (MAXB       = 12) 
      PARAMETER      (MAXGSAS    = 25000) 
!
!#######################################################################
