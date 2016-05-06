MODULE discus_save_mod
!-
!     Variables needed to determine which values are written as 
!     keywords to the structure file, and also to signal, which
!     values were read from the structure file
!+
USE discus_config_mod
!
SAVE
!
INTEGER               ::  SAV_T_MAXSCAT,SAV_MAXSCAT  =  1
LOGICAL, DIMENSION(:), ALLOCATABLE  ::  sav_latom ! (0:MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE  ::  sav_t_latom ! (0:MAXSCAT)
!
LOGICAL               ::  sav_t_sel_atom, sav_sel_atom = .true.
!
CHARACTER(LEN=200)    ::  sav_t_file    , sav_file     = 'crystal.stru'
!
LOGICAL               ::  sav_t_keyword , sav_keyword  = .true.
!
LOGICAL               ::  sav_t_w_scat  , sav_w_scat   = .false.
LOGICAL               ::  sav_t_w_adp   , sav_w_adp    = .false.
LOGICAL               ::  sav_t_r_ncell , sav_r_ncell  = .false.
LOGICAL               ::  sav_t_w_ncell , sav_w_ncell  = .false.
LOGICAL               ::  sav_t_w_gene  , sav_w_gene   = .true.
LOGICAL               ::  sav_t_w_symm  , sav_w_symm   = .true.
LOGICAL               ::  sav_t_w_mole  , sav_w_mole   = .true.
LOGICAL               ::  sav_t_w_obje  , sav_w_obje   = .true.
LOGICAL               ::  sav_t_w_doma  , sav_w_doma   = .true.
LOGICAL               ::  sav_t_w_prop  , sav_w_prop   = .true.
!
INTEGER               ::  sav_t_start   , sav_start    =  1
INTEGER               ::  sav_t_end     , sav_end      = -1
INTEGER, DIMENSION(3) ::  sav_t_ncell   , sav_ncell    =  1
INTEGER, DIMENSION(0:1):: sav_t_sel_prop, sav_sel_prop =  (/0,0/)
INTEGER               ::  sav_t_ncatoms , sav_ncatoms  =  1
INTEGER               ::  sav_t_size_of , sav_size_of  =  0 ! Bytes allocated for save
!
END MODULE discus_save_mod
