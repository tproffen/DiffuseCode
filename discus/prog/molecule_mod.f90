MODULE molecule_mod
!-
!     Variables needed to describe molecules. 
!+
use precision_mod
!
INTEGER             ::  MOLE_MAX_MOLE  =       1
INTEGER             ::  MOLE_MAX_TYPE  =       1
INTEGER             ::  MOLE_MAX_GENE  =       1
INTEGER             ::  MOLE_MAX_SYMM  =       1
INTEGER             ::  MOLE_MAX_ATOM  =       1
!
INTEGER, PARAMETER  ::  MOLE_ATOM          =   0
INTEGER, PARAMETER  ::  MOLE_CUBE          =   1
INTEGER, PARAMETER  ::  MOLE_CYLINDER      =   2
INTEGER, PARAMETER  ::  MOLE_SPHERE        =   3
INTEGER, PARAMETER  ::  MOLE_EDGE          =   4
INTEGER, PARAMETER  ::  MOLE_DOM_CUBE      =  -1
INTEGER, PARAMETER  ::  MOLE_DOM_CYLINDER  =  -2
INTEGER, PARAMETER  ::  MOLE_DOM_SPHERE    =  -3
INTEGER, PARAMETER  ::  MOLE_DOM_FUZZY     =  -4
!
LOGICAL                                 ::  mole_l_on     = .false.
LOGICAL                                 ::  mole_l_first  = .false.
!
INTEGER                                 ::  mole_num_mole = 0
INTEGER                                 ::  mole_num_curr = 0
INTEGER                                 ::  mole_num_act  = 0
INTEGER                                 ::  mole_num_type = 0
INTEGER                                 ::  mole_num_unit = 0
INTEGER                                 ::  mole_gene_n   = 0
INTEGER                                 ::  mole_symm_n   = 0
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_gene_power           ! (1:MOLE_MAX_GENE)
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_symm_power           ! (1:MOLE_MAX_SYMM)
!
INTEGER                                 ::  mole_num_atom = 0
INTEGER                                 ::  mole_num_acur = 0
!
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_len                  ! (  0:MOLE_MAX_MOLE)
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_off                  ! (  0:MOLE_MAX_MOLE)
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_type                 ! (  0:MOLE_MAX_MOLE)
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_char                 ! (  0:MOLE_MAX_MOLE)
CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE :: mole_file            ! (  0:MOLE_MAX_MOLE)
INTEGER, DIMENSION(:), ALLOCATABLE      ::  mole_cont                 ! (  0:MOLE_MAX_ATOM)
!
REAL(kind=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  ::  mole_gene                 ! (4,4,1:MOLE_MAX_GENE)
!        RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,MOLE_MAX_GENE)/),SHAPE(mole_gene ))
REAL(kind=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  ::  mole_symm                 ! (4,4,1:MOLE_MAX_SYMM)
!        RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,MOLE_MAX_SYMM)/),SHAPE(mole_symm ))
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE      ::  mole_dens                 ! (  0:MOLE_MAX_MOLE)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE      ::  mole_biso                 ! (  0:MOLE_MAX_TYPE)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE      ::  mole_clin                 ! (  0:MOLE_MAX_TYPE)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE      ::  mole_cqua                 ! (  0:MOLE_MAX_TYPE)
REAL(kind=PREC_DP)   , DIMENSION(:), ALLOCATABLE      ::  mole_fuzzy                ! (  0:MOLE_MAX_MOLE)
!
INTEGER                                 :: mol_size_of = 0
!
END MODULE molecule_mod
