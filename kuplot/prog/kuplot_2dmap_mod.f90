MODULE kuplot_2dm_mod
!
USE precision_mod
!
PUBLIC
!
INTEGER, PARAMETER :: K2DM_XY  = 1
INTEGER, PARAMETER :: K2DM_CSV = 2
INTEGER, PARAMETER :: K2DM_YF_LOOP = 1
INTEGER, PARAMETER :: K2DM_YF_INC  = 2
INTEGER, PARAMETER :: K2DM_YF_FUNC = 3
INTEGER, PARAMETER :: K2DM_ERROR   = -1
INTEGER, PARAMETER :: K2DM_BLANK   =  0
INTEGER, PARAMETER :: K2DM_IGNORE  =  1
!
CHARACTER(LEN=4)    :: k2dm_ctype  = ' '
CHARACTER(LEN=PREC_STRING) :: k2dm_line   = ' '
CHARACTER(LEN=PREC_STRING) :: k2dm_line_b = ' '
CHARACTER(LEN=PREC_STRING) :: k2dm_line_yf = ' '
INTEGER          :: k2dm_type      = 0
INTEGER          :: k2dm_miss      = -1
LOGICAL          :: k2dm_miss_set  = .FALSE.
INTEGER, DIMENSION(0:2)          :: k2dm_start     = 0
INTEGER, DIMENSION(0:2)          :: k2dm_end       = 0
INTEGER, DIMENSION(0:2)          :: k2dm_step      = 0
LOGICAL          :: k2dm_lxmin     = .TRUE.
LOGICAL          :: k2dm_lxmax     = .TRUE.
REAL             :: k2dm_scale     = 1.0
REAL             :: k2dm_xmin      = 0.0
REAL             :: k2dm_xmax      = 0.0
INTEGER          :: k2dm_type_yf   = K2DM_YF_LOOP
!
END MODULE kuplot_2dm_mod
