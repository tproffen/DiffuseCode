MODULE inter_readstru
!
IMPLICIT NONE
!
INTEGER rd_NMAX 
INTEGER rd_MAXSCAT 
!
CHARACTER (LEN=1024)                  :: rd_strucfile 
CHARACTER (LEN=  80)                  :: rd_cr_name 
CHARACTER (LEN=  16)                  :: rd_cr_spcgr 
REAL                , DIMENSION(3)    :: rd_cr_a0
REAL                , DIMENSION(3)    :: rd_cr_win
INTEGER                               :: rd_cr_natoms
INTEGER                               :: rd_cr_nscat 
REAL                , DIMENSION(:),   ALLOCATABLE :: rd_cr_dw     ! (0:MAXSCAT) 
CHARACTER (LEN=   4), DIMENSION(:),   ALLOCATABLE :: rd_cr_at_lis ! (0:MAXSCAT) 
REAL                , DIMENSION(:,:), ALLOCATABLE :: rd_cr_pos    ! (1:3,1:NMAX)
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_cr_iscat  ! (1:NMAX)
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_cr_prop   ! (1:NMAX)
REAL                , DIMENSION(3,2)  :: rd_cr_dim
INTEGER                               :: rd_as_natoms 
CHARACTER (LEN=   4), DIMENSION(:),   ALLOCATABLE :: rd_as_at_lis ! (0:MAXSCAT) 
REAL                , DIMENSION(:),   ALLOCATABLE :: rd_as_dw     ! (0:MAXSCAT) 
REAL                , DIMENSION(:,:), ALLOCATABLE :: rd_as_pos    ! (3, MAXSCAT) 
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_as_iscat  ! (1:MAXSCAT)
INTEGER             , DIMENSION(:),   ALLOCATABLE :: rd_as_prop   ! (1:MAXSCAT)
INTEGER             , DIMENSION(3)    :: rd_sav_ncell ! (3) 
LOGICAL                               :: rd_sav_r_ncell 
INTEGER                               :: rd_sav_ncatoms 
INTEGER                               :: rd_spcgr_ianz 
INTEGER                               :: rd_spcgr_para 
!                                                                       
!
END MODULE inter_readstru
