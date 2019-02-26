MODULE kuplot_3dm_mod
!
PUBLIC
!
INTEGER            :: k3dm_ik       = 0
LOGICAL            :: k3dm_lflat    = .TRUE.
REAL               :: k3dm_phi      = 0.0
REAL               :: k3dm_rho      = 0.0
!
INTEGER, PARAMETER :: K3DM_XY  = 1
INTEGER, PARAMETER :: K3DM_CSV = 2
INTEGER, PARAMETER :: K3DM_YF_LOOP = 1
INTEGER, PARAMETER :: K3DM_YF_INC  = 2
INTEGER, PARAMETER :: K3DM_YF_FUNC = 3
INTEGER, PARAMETER :: K3DM_ERROR   = -1
INTEGER, PARAMETER :: K3DM_BLANK   =  0
INTEGER, PARAMETER :: K3DM_IGNORE  =  1
!
CHARACTER(LEN=4)    :: k3dm_ctype  = ' '
CHARACTER(LEN=1024) :: k3dm_line   = ' '
CHARACTER(LEN=1024) :: k3dm_line_b = ' '
CHARACTER(LEN=1024) :: k3dm_line_yf = ' '

INTEGER          :: k3dm_type      = 0
INTEGER          :: k3dm_miss      = -1
INTEGER          :: k3dm_start     = 0
INTEGER          :: k3dm_end       = 0
INTEGER          :: k3dm_step      = 0
LOGICAL          :: k3dm_lxmin     = .TRUE.
LOGICAL          :: k3dm_lxmax     = .TRUE.
REAL             :: k3dm_scale     = 1.0
REAL             :: k3dm_xmin      = 0.0
REAL             :: k3dm_xmax      = 0.0
INTEGER          :: k3dm_type_yf   = K3DM_YF_LOOP
!
END MODULE kuplot_3dm_mod
