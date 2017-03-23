MODULE kuplot_config
!
IMPLICIT NONE
PUBLIC
SAVE
!
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
      INTEGER, PARAMETER :: MAXPARA    =      150
      INTEGER, PARAMETER :: MAXKURVTOT =      200
      INTEGER, PARAMETER :: MAXARRAY   = 13000000
      INTEGER, PARAMETER :: MAXHL      =        5

      INTEGER, PARAMETER :: MAXFRAME   =       16
      INTEGER, PARAMETER :: MAXSP      =      500
      INTEGER, PARAMETER :: MAXAN      =      500
      INTEGER, PARAMETER :: MAXCOL     =      255
      INTEGER, PARAMETER :: MAXWIN     =        3

      INTEGER, PARAMETER :: MAXBOND    =       50
      INTEGER, PARAMETER :: MAXZ       =     1501
      INTEGER, PARAMETER :: MAXB       =       12
      INTEGER, PARAMETER :: MAXGSAS    =    25000
!
!
!#######################################################################
END MODULE kuplot_config
