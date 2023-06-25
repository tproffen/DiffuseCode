MODULE crystal_mod
!
! This is the module for the crystal data, atom positions etc.
!
!-
!     Definition of all variables for the crystal
!+
!
!
USE discus_config_mod
use precision_mod
!
INTEGER, PARAMETER               ::  cr_triclinic   = 1
INTEGER, PARAMETER               ::  cr_monoclinicB = 2
INTEGER, PARAMETER               ::  cr_monoclinicC = 3
INTEGER, PARAMETER               ::  cr_ortho       = 4
INTEGER, PARAMETER               ::  cr_tetragonal  = 5
INTEGER, PARAMETER               ::  cr_rhombohed   = 6
INTEGER, PARAMETER               ::  cr_trigonal    = 7
INTEGER, PARAMETER               ::  cr_hexagonal   = 8
INTEGER, PARAMETER               ::  cr_cubic       = 9
!
INTEGER                          ::  cr_spcgrno = 1  ! Space group number
INTEGER                          ::  cr_syst    = 1  ! Crystal system number
INTEGER                          ::  cr_ncatoms = 1  ! NUmber of atoms in Unit cell
INTEGER                          ::  cr_ncreal  = 1  ! NUmber of non-void atoms in Unit cell
INTEGER, DIMENSION(3)            ::  cr_icc     = 1  ! Number of unit cells
!
INTEGER                          ::  spcgr_ianz  = 0
INTEGER                          ::  spcgr_para  = 0
INTEGER, DIMENSION(0:1)          ::  cr_sel_prop = 0
!
REAL(kind=PREC_DP)                          ::  cr_v        = 1.0D0
REAL(kind=PREC_DP), DIMENSION(3,3)          ::  cr_gten     = &
   RESHAPE((/1.D0,0.D0,0.D0, 0.D0,1.D0,0.D0, 0.D0,0.D0,1.D0/),SHAPE(cr_gten))
REAL(kind=PREC_DP), DIMENSION(3,3,3)        ::  cr_eps
REAL(kind=PREC_DP), DIMENSION(3,3)          ::  cr_gmat     = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),SHAPE(cr_gmat))
REAL(kind=PREC_DP), DIMENSION(3)            ::  cr_ar       = 1.0D0
REAL(kind=PREC_DP), DIMENSION(3)            ::  cr_wrez     = 90.00D0
REAL(kind=PREC_DP)                          ::  cr_vr       = 1.0D0
REAL(kind=PREC_DP), DIMENSION(3,3)          ::  cr_rten     = &
   RESHAPE((/1.D0,0.D0,0.D0, 0.D0,1.D0,0.D0, 0.D0,0.D0,1.D0/),SHAPE(cr_rten))
REAL(kind=PREC_DP), DIMENSION(3,3,3)        ::  cr_reps
REAL(kind=PREC_DP), DIMENSION(3,3)          ::  cr_fmat     = &
   RESHAPE((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(cr_fmat))
REAL(kind=PREC_DP), DIMENSION(3,3)          ::  cr_dmat     = &
   RESHAPE((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(cr_fmat))
REAL(kind=PREC_DP), DIMENSION(3,3)          ::  cr_dimat     = &
   RESHAPE((/1.0D0,0.0D0,0.0D0, 0.0D0,1.0D0,0.0D0, 0.0D0,0.0D0,1.0D0/),SHAPE(cr_fmat))
REAL(kind=PREC_DP), DIMENSION(4,4)          ::  cr_tran_g   = 0.0D0
REAL(kind=PREC_DP), DIMENSION(4,4)          ::  cr_tran_gi  = 0.0D0
REAL(kind=PREC_DP), DIMENSION(4,4)          ::  cr_tran_f   = 0.0D0
REAL(kind=PREC_DP), DIMENSION(4,4)          ::  cr_tran_fi  = 0.0D0
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE ::  cr_scat   ! (11,0:MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE ::  cr_delfr  ! (  0:MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE ::  cr_delfi  ! (  0:MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE ::  cr_delfr_u! (  0:MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(  :), ALLOCATABLE ::  cr_delfi_u! (  0:MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(3,2)          ::  cr_dim0     = &
   RESHAPE((/0.0D0,0.0D0,0.0D0, 1.0D0,1.0D0,1.0D0/),SHAPE(cr_dim0))
!
LOGICAL                               ::  cr_acentric = .true.
LOGICAL, DIMENSION(  :), ALLOCATABLE  ::  cr_scat_int  ! (  0:MAXSCAT)
LOGICAL, DIMENSION(  :), ALLOCATABLE  ::  cr_scat_equ  ! (  0:MAXSCAT)
LOGICAL, DIMENSION(  :), ALLOCATABLE  ::  cr_delf_int  ! (  0:MAXSCAT)
LOGICAL                               ::  cr_newtype = .true.
!
!
CHARACTER (LEN=80)                               ::  cr_name  = 'crystal'
CHARACTER (LEN=16)                               ::  cr_spcgr = 'P1'
CHARACTER (LEN= 3)                               ::  cr_set   = 'abc'
INTEGER                                          ::  cr_iset  =  1
CHARACTER (LEN=16)                               ::  cr_spcgr_set = 'P1'
CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_lis  ! (  0:MAXSCAT)
CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_equ  ! (  0:MAXSCAT)
CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  as_at_lis  ! (  0:MAXSCAT)
!
INTEGER                               ::  cr_nscat  = 0
INTEGER                               ::  as_natoms = 0
INTEGER, DIMENSION(  :), ALLOCATABLE  ::  as_iscat  ! (  1:MAXSCAT)
INTEGER, DIMENSION(  :), ALLOCATABLE  ::  as_mole   ! (  1:MAXSCAT)
INTEGER, DIMENSION(  :), ALLOCATABLE  ::  as_prop   ! (  1:MAXSCAT)
INTEGER, DIMENSION(  :), ALLOCATABLE  ::  cr_niscat ! (  1:MAXSCAT)  Number of atoms ot type iscat
!
LOGICAL                               :: cr_cartesian = .false.
LOGICAL                               :: cr_magnetic  = .FALSE.   ! Cyrystal is magnetic YES / NO
!
REAL(kind=PREC_DP)   , DIMENSION(3)                 ::  cr_a0  = 1.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)                 ::  cr_win = 90.0D0
REAL(kind=PREC_DP)   , DIMENSION(3,2)               ::  cr_dim = RESHAPE((/0D0,0D0,0D0, 1D0,1D0,1D0/),SHAPE(cr_dim))
REAL(kind=PREC_DP)   , DIMENSION(  :), ALLOCATABLE  ::  cr_dw  ! (  0:MAXSCAT)
REAL(kind=PREC_DP)   , DIMENSION(  :), ALLOCATABLE  ::  cr_occ ! (  0:MAXSCAT)
REAL(kind=PREC_DP)   , DIMENSION(  :), ALLOCATABLE  ::  as_occ ! (  0:MAXSCAT)
REAL(kind=PREC_DP)   , DIMENSION(:,:), ALLOCATABLE  ::  as_pos ! (3,1:MAXSCAT)
REAL(kind=PREC_DP)   , DIMENSION(  :), ALLOCATABLE  ::  as_dw  ! (  0:MAXSCAT)
integer              , dimension(  :), allocatable  ::  cr_is_sym ! Site was created by Sym.Elem no. i 
REAL(kind=PREC_DP)   , DIMENSION(:,:), ALLOCATABLE  ::  cr_anis! U(ij) ADP anis  ! 0: MAXSCAT
REAL(kind=PREC_DP)   , DIMENSION(:,:), ALLOCATABLE  ::  cr_anis_full! U(ij) ADP anis  ! 0: number of atoms per unit cell
REAL(kind=PREC_DP)   , DIMENSION(:,:,:), ALLOCATABLE  ::  cr_prin ! Eigenvectors    ! 0: number of atoms per unit cell
REAL(kind=PREC_DP)   , DIMENSION(:,  :), ALLOCATABLE  ::  cr_u2   ! <u^2> Anisotrop ! 0: number of atoms per unit cell
!
INTEGER, DIMENSION(  :), ALLOCATABLE  ::  cr_amount ! (  0:MAXSCAT)
REAL(kind=PREC_DP)                    ::  cr_u2aver
!
REAL(kind=PREC_DP)                    :: cr_mass   = 0.0D0 ! Crystal mass in u
REAL(kind=PREC_DP)                    :: cr_nreal  = 0.0D0 ! Real amount of atoms including occupancy
!
INTEGER                              ::  cr_natoms       = 0
INTEGER                              ::  cr_n_REAL_atoms = 0
INTEGER, DIMENSION(  :), ALLOCATABLE ::  cr_iscat  ! (  1:NMAX)  !Atom type 0 to cr_nscat
INTEGER, DIMENSION(  :), ALLOCATABLE ::  cr_prop   ! (  1:NMAX)  !Property flag
INTEGER, DIMENSION(  :), ALLOCATABLE ::  cr_mole   ! (  1:NMAX)  !Atom is in this molecule
INTEGER, DIMENSION(:,:), ALLOCATABLE ::  cr_surf   ! (  1:NMAX)  !Atom is on this surface 
REAL(kind=PREC_DP)   , DIMENSION(:,:), ALLOCATABLE ::  cr_magn   ! (  1:NMAX)  !Magnetic moment 
!
REAL(kind=PREC_DP), DIMENSION(:,:), ALLOCATABLE ::  cr_pos    ! (3,1:NMAX)  !Atom coordinates
!
INTEGER                               ::  cry_size_of = 0 ! Bytes allocated fo
!
END MODULE crystal_mod
