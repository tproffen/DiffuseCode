MODULE save_mod
!-
!     Variables needed to determine which values are written as 
!     keywords to the structure file, and also to signal, which
!     values were read from the structure file
!+
USE config_mod
!
SAVE
!

INTEGER             ::  SAV_MAXSCAT  =  1
LOGICAL, DIMENSION(:), ALLOCATABLE  ::  sav_latom ! (0:MAXSCAT)
!
LOGICAL               ::  sav_sel_atom = .true.
!
CHARACTER(LEN=200)    ::  sav_file     = 'crystal.stru'
!
LOGICAL               ::  sav_keyword  = .true.
!
LOGICAL               ::  sav_w_scat   = .false.
LOGICAL               ::  sav_w_adp    = .false.
LOGICAL               ::  sav_r_ncell  = .false.
LOGICAL               ::  sav_w_ncell  = .false.
LOGICAL               ::  sav_w_gene   = .true.
LOGICAL               ::  sav_w_symm   = .true.
LOGICAL               ::  sav_w_mole   = .true.
LOGICAL               ::  sav_w_obje   = .true.
LOGICAL               ::  sav_w_doma   = .true.
LOGICAL               ::  sav_w_prop   = .true.
!
INTEGER               ::  sav_start    =  1
INTEGER               ::  sav_end      = -1
INTEGER, DIMENSION(3) ::  sav_ncell(3) =  1
INTEGER               ::  sav_ncatoms  =  1
INTEGER               ::  sav_size_of  =  0 ! Bytes allocated for save
!
END MODULE save_mod
