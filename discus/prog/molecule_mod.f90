MODULE molecule_mod
!-
!     Variables needed to describe molecules. 
!+
INTEGER, PRIVATE    ::  ik
INTEGER, PRIVATE    ::  il
INTEGER, PARAMETER  ::  MOLE_MAX_MOLE  =   64000
INTEGER, PARAMETER  ::  MOLE_MAX_TYPE  =     500
INTEGER, PARAMETER  ::  MOLE_MAX_GENE  =      10
INTEGER, PARAMETER  ::  MOLE_MAX_SYMM  =      48
INTEGER, PARAMETER  ::  MOLE_MAX_ATOM  =  512000
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
INTEGER, DIMENSION(MOLE_MAX_GENE)       ::  mole_gene_power = 1
INTEGER, DIMENSION(MOLE_MAX_SYMM)       ::  mole_symm_power = 1
!
INTEGER                                 ::  mole_num_atom = 0
INTEGER                                 ::  mole_num_acur = 0
!
INTEGER  , DIMENSION(  0:MOLE_MAX_MOLE) ::  mole_len      = 0
INTEGER  , DIMENSION(  0:MOLE_MAX_MOLE) ::  mole_off      = 0
INTEGER*1, DIMENSION(  0:MOLE_MAX_MOLE) ::  mole_type     = 0
INTEGER*1, DIMENSION(  0:MOLE_MAX_MOLE) ::  mole_char     = MOLE_ATOM
CHARACTER(LEN=200), DIMENSION(0:MOLE_MAX_MOLE)  ::  mole_file = ' '
INTEGER, DIMENSION(    0:MOLE_MAX_ATOM) ::  mole_cont     = 0
!
REAL   , DIMENSION(4,4,0:MOLE_MAX_GENE) ::  mole_gene     = &
         RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,MOLE_MAX_GENE)/),SHAPE(mole_gene ))
REAL   , DIMENSION(4,4,0:MOLE_MAX_SYMM) ::  mole_symm     = &
         RESHAPE((/(1.,(0.,0.,0.,0.,1.,ik=1,3),il=0,MOLE_MAX_SYMM)/),SHAPE(mole_symm ))
REAL   , DIMENSION(    0:MOLE_MAX_MOLE) ::  mole_dens     = 0.0
REAL   , DIMENSION(    0:MOLE_MAX_MOLE) ::  mole_fuzzy    = 0.0
!
!     COMMON /molecule/ mole_num_mole,mole_num_curr,mole_num_act,       &
!    &                  mole_num_type,mole_num_unit,mole_l_on,          &
!    &                  mole_l_first,mole_gene_n,mole_symm_n,           &
!    &                  mole_gene_power,mole_symm_power,                &
!    &                  mole_num_atom,mole_num_acur,                    &
!    &                  mole_gene,mole_symm,                            &
!    &                  mole_len,mole_off,mole_cont,                    &
!    &                  mole_dens,mole_fuzzy,mole_file,                 &
!    &                  mole_type,mole_char
!
END MODULE molecule_mod
