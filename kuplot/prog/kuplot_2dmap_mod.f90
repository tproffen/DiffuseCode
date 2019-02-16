MODULE kuplot_2dm_mod
!
PUBLIC
!
INTEGER, PARAMETER :: K2DM_XY  = 1
INTEGER, PARAMETER :: K2DM_CSV = 2
!
CHARACTER(LEN=4)    :: k2dm_ctype
CHARACTER(LEN=1024) :: k2dm_line
INTEGER          :: k2dm_type
INTEGER          :: k2dm_loop1
INTEGER          :: k2dm_nanz
CHARACTER(LEN=1024), DIMENSION(:),  ALLOCATABLE :: k2dm_cpara
INTEGER            , DIMENSION(:),  ALLOCATABLE :: k2dm_lpara
INTEGER          :: k2dm_start
INTEGER          :: k2dm_end
INTEGER          :: k2dm_step
LOGICAL          :: k2dm_lxmin
LOGICAL          :: k2dm_lxmax
REAL             :: k2dm_xmin
REAL             :: k2dm_xmax
!
END MODULE kuplot_2dm_mod
