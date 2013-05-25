MODULE domain_mod
!
!
!     Common block for the cluster distributions
!
INTEGER             ::  clu_increment     = 200
INTEGER             ::  CLU_MAX_TYPE      =   1
!
INTEGER, PARAMETER  ::  CLU_IN_CLUSTER    = 0
INTEGER, PARAMETER  ::  CLU_IN_PSEUDO     = 1
!
INTEGER, PARAMETER  ::  CLU_CHAR_CUBE     = -1
INTEGER, PARAMETER  ::  CLU_CHAR_CYLINDER = -2
INTEGER, PARAMETER  ::  CLU_CHAR_SPHERE   = -3
INTEGER, PARAMETER  ::  CLU_CHAR_FUZZY    = -4
!
CHARACTER(LEN=200)                    ::  clu_infile
CHARACTER(LEN=200), DIMENSION(:),     ALLOCATABLE  ::  clu_content   ! (CLU_MAX_TYPE)
CHARACTER(LEN=  4), DIMENSION(:),     ALLOCATABLE  ::  clu_name      ! (CLU_MAX_TYPE)
!
INTEGER,            DIMENSION(:)    , ALLOCATABLE  ::  clu_character ! (CLU_MAX_TYPE)
REAL   ,            DIMENSION(:)    , ALLOCATABLE  ::  clu_fuzzy  ! (CLU_MAX_TYPE)
REAL   ,            DIMENSION(:,:,:), ALLOCATABLE  ::  clu_orient ! (CLU_MAX_TYPE,3,4)
REAL   ,            DIMENSION(:,:,:), ALLOCATABLE  ::  clu_shape  ! (CLU_MAX_TYPE,3,4)
INTEGER                               ::  clu_index
INTEGER                               ::  clu_mode
INTEGER                               ::  clu_number
!
LOGICAL                               ::  clu_surface
LOGICAL                               ::  clu_infile_internal = .false. ! Is infile an internal file ?
INTEGER                               ::  clu_iatom = 0
!
INTEGER                               ::  clu_size_of! Bytes allocated for DOMAIN
!
!     COMMON /cluster/ clu_infile,clu_content,clu_name,clu_mode,        &
!    &                 clu_index,clu_number,clu_CHARACTER,              &
!    &                 clu_fuzzy,clu_orient,clu_shape
!
END MODULE domain_mod
