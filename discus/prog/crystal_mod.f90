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
   REAL                             ::  cr_v        = 1.
   REAL   , DIMENSION(3,3)          ::  cr_gten     = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),SHAPE(cr_gten))
   REAL   , DIMENSION(3,3,3)        ::  cr_eps
   REAL   , DIMENSION(3,3)          ::  cr_gmat     = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),SHAPE(cr_gmat))
   REAL   , DIMENSION(3)            ::  cr_ar       = 1.0
   REAL   , DIMENSION(3)            ::  cr_wrez     = 90.00
   REAL                             ::  cr_vr       = 1.
   REAL   , DIMENSION(3,3)          ::  cr_rten     = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),SHAPE(cr_rten))
   REAL   , DIMENSION(3,3,3)        ::  cr_reps
   REAL   , DIMENSION(3,3)          ::  cr_fmat     = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),SHAPE(cr_fmat))
   REAL   , DIMENSION(4,4)          ::  cr_tran_g   = 0.0
   REAL   , DIMENSION(4,4)          ::  cr_tran_gi  = 0.0
   REAL   , DIMENSION(4,4)          ::  cr_tran_f   = 0.0
   REAL   , DIMENSION(4,4)          ::  cr_tran_fi  = 0.0
   REAL   , DIMENSION(:,:), ALLOCATABLE ::  cr_scat   ! (11,0:MAXSCAT)
   REAL   , DIMENSION(  :), ALLOCATABLE ::  cr_delfr  ! (  0:MAXSCAT)
   REAL   , DIMENSION(  :), ALLOCATABLE ::  cr_delfi  ! (  0:MAXSCAT)
   REAL   , DIMENSION(3,2)          ::  cr_dim0     =RESHAPE((/0,0,0, 1,1,1/),SHAPE(cr_dim0))
!
   LOGICAL                          ::  cr_acentric = .true.
   LOGICAL, DIMENSION(  :), ALLOCATABLE  ::  cr_scat_int  ! (  0:MAXSCAT)
   LOGICAL, DIMENSION(  :), ALLOCATABLE  ::  cr_scat_equ  ! (  0:MAXSCAT)
   LOGICAL, DIMENSION(  :), ALLOCATABLE  ::  cr_delf_int  ! (  0:MAXSCAT)
   LOGICAL                          ::  cr_newtype = .true.
!
!
   CHARACTER (LEN=80)                          ::  cr_name  = 'crystal'
   CHARACTER (LEN=16)                          ::  cr_spcgr = 'P1'
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_lis  ! (  0:MAXSCAT)
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_equ  ! (  0:MAXSCAT)
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  as_at_lis  ! (  0:MAXSCAT)
!
   INTEGER                          ::  cr_nscat  = 0
   INTEGER                          ::  as_natoms = 0
   INTEGER, DIMENSION(  :), ALLOCATABLE  ::  as_iscat  ! (  1:MAXSCAT)
   INTEGER, DIMENSION(  :), ALLOCATABLE  ::  as_mole   ! (  1:MAXSCAT)
   INTEGER, DIMENSION(  :), ALLOCATABLE  ::  as_prop   ! (  1:MAXSCAT)
!
   LOGICAL                               :: cr_cartesian = .false.
!
   REAL   , DIMENSION(3)                 ::  cr_a0  = 1.0
   REAL   , DIMENSION(3)                 ::  cr_win = 90.0
   REAL   , DIMENSION(3,2)               ::  cr_dim = RESHAPE((/0,0,0, 1,1,1/),SHAPE(cr_dim))
   REAL   , DIMENSION(  :), ALLOCATABLE  ::  cr_dw  ! (  0:MAXSCAT)
   REAL   , DIMENSION(:,:), ALLOCATABLE  ::  as_pos ! (3,1:MAXSCAT)
   REAL   , DIMENSION(  :), ALLOCATABLE  ::  as_dw  ! (  0:MAXSCAT)
!
   INTEGER, DIMENSION(  :), ALLOCATABLE  ::  cr_amount ! (  0:MAXSCAT)
   REAL                                  ::  cr_u2aver
!
   INTEGER                              ::  cr_natoms       = 0
   INTEGER                              ::  cr_n_REAL_atoms = 0
   INTEGER, DIMENSION(  :), ALLOCATABLE ::  cr_iscat  ! (  1:NMAX)  !Atom type 0 to cr_nscat
   INTEGER, DIMENSION(  :), ALLOCATABLE ::  cr_prop   ! (  1:NMAX)  !Property flag
   INTEGER, DIMENSION(  :), ALLOCATABLE ::  cr_mole   ! (  1:NMAX)  !Atom is in this molecule
   INTEGER, DIMENSION(:,:), ALLOCATABLE ::  cr_surf   ! (  1:NMAX)  !Atom is on this surface 
!
   REAL   , DIMENSION(:,:), ALLOCATABLE ::  cr_pos    ! (3,1:NMAX)  !Atom coordinates
!
INTEGER                               ::  cry_size_of = 0 ! Bytes allocated fo
!
END MODULE crystal_mod
