MODULE cryst_class
!
!  Defines a type "cl_cryst".
!  This very incomplete example just has the natoms i.e. 
!  the number of atoms and an array of type "cl_atom"
!
!  Just a few procedures are bound to this class:
!      alloc_atoms      allocate the atom array
!      get_natoms       Returns the number of atoms in this crystal
!      set_cryst_atom   Sets all values for an atom
!      get_cryst_atom   Gets all values for an atom
!
!  A working class will need the full suite of methods....
!
USE atom_class          ! this way we know all interfaces to the "cl_atom" class
!
IMPLICIT none 
!
TYPE :: cl_cryst        ! Define a type "cl_cryst"
   PRIVATE              ! hide local variables from outside; access through procedures only
!
   CHARACTER (LEN=200)                              ::  cr_file     = 'crystal.stru' ! Crystal file name
   CHARACTER (LEN=80)                               ::  cr_name     = 'crystal' ! Crystal name
   CHARACTER (LEN=16)                               ::  cr_spcgr    = 'P1'      ! Space group symbol
   INTEGER                                          ::  cr_spcgrno  = 1 ! Space group number
   INTEGER                                          ::  cr_syst     = 1 ! Crystal system
   LOGICAL                                          ::  cr_acentric = .true.
   INTEGER                                          ::  cr_ncatoms  = 1 ! Atoms in unit cell
   INTEGER           , DIMENSION(3)                 ::  cr_icc      = 1 ! Number of unit cells
!
   INTEGER                                          ::  spcgr_ianz  = 0 ! Number of parameters on spcg line
   INTEGER                                          ::  spcgr_para  = 0 ! 
   INTEGER           , DIMENSION(0:1)               ::  cr_sel_prop = 0 ! Global selection properties
!
   REAL              , DIMENSION(3)                 ::  cr_a0   =  1.0  ! unit cell lengths
   REAL              , DIMENSION(3)                 ::  cr_win  = 90.0  ! unit cell angles
   REAL              , DIMENSION(3,2)               ::  cr_dim  = RESHAPE((/0,0,0, 1,1,1/),(/3,2/))
   REAL              , DIMENSION(3,2)               ::  cr_dim0 = RESHAPE((/0,0,0, 1,1,1/),(/3,2/))
   REAL                                             ::  cr_v    = 1.    ! Unit cell volume
   REAL              , DIMENSION(3,3)               ::  cr_gten = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(3,3,3)             ::  cr_eps          ! Epsilon tensor
   REAL              , DIMENSION(3,3)               ::  cr_gmat = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(3)                 ::  cr_ar   = 1.0   ! reciprocal unit cell lengths
   REAL              , DIMENSION(3)                 ::  cr_wrez = 90.00 ! reciprocal unit cell angles
   REAL                                             ::  cr_vr   = 1.    ! ceiprocal unit cell volume
   REAL              , DIMENSION(3,3)               ::  cr_rten = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(3,3,3)             ::  cr_reps         ! reciprocel epsilon tensor
   REAL              , DIMENSION(3,3)               ::  cr_fmat = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
!
   REAL              , DIMENSION(:,:), ALLOCATABLE  ::  cr_scat   ! Scattering curves (9,0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_delfr  ! anomalous scattering real part (  0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_delfi  ! anomalous scattering imag part (  0:MAXSCAT)
!
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_scat_int  ! (  0:MAXSCAT)
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_scat_equ  ! (  0:MAXSCAT)
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_delf_int  ! (  0:MAXSCAT)
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_equ    ! (  0:MAXSCAT)
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_sav_atom  ! (  0:MAXSCAT)
   LOGICAL                                          ::  cr_newtype = .true.
   LOGICAL                                          ::  cr_cartesian = .false.
   LOGICAL                                          ::  cr_sav_scat  = .false.
   LOGICAL                                          ::  cr_sav_adp   = .false.
   LOGICAL                                          ::  cr_sav_gene  = .true.
   LOGICAL                                          ::  cr_sav_symm  = .true.
   LOGICAL                                          ::  cr_sav_ncell = .false.
   LOGICAL                                          ::  cr_sav_obje  = .true.
   LOGICAL                                          ::  cr_sav_doma  = .true.
   LOGICAL                                          ::  cr_sav_mole  = .true.
!
   INTEGER                                          ::  cr_natoms  !the number of atoms
   INTEGER                                          ::  cr_nscat   !the number of atom types
   INTEGER                                          ::  cr_n_REAL_atoms = 0

   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_lis  ! (  0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_dw      ! (  0:MAXSCAT)
   TYPE(cl_atom)     , DIMENSION(:)  , POINTER      ::  atoms            !the actual atoms
   LOGICAL                                          ::  latom = .false.  !have atoms been allocated?
!
CONTAINS
!
   PROCEDURE, PUBLIC, PASS :: alloc_arrays         ! allocate a set of atoms
   PROCEDURE, PUBLIC, PASS :: get_natoms           ! Return number of atoms in this crystal
   PROCEDURE, PUBLIC, PASS :: get_nscat            ! Return number of atom types in this crystal
   PROCEDURE, PUBLIC, PASS :: set_cryst_atom       ! Set iscat, posit, property
   PROCEDURE, PUBLIC, PASS :: get_cryst_atom       ! Get iscat, posit, property
   PROCEDURE, PUBLIC, PASS :: read_crystal         ! Read from disk file
   PROCEDURE, PUBLIC, PASS :: set_crystal_from_standard ! Copy from DISCUS standard into a crystal type
   PROCEDURE, PUBLIC, PASS :: set_crystal_save_flags    ! Set the "save" flags in the crystal type
   PROCEDURE, PUBLIC, PASS :: get_standard_from_crystal ! Copy from a crystal type into DISCUS standard
   PROCEDURE, PUBLIC, PASS :: is_alloc_atom        ! Has atom array been allocated ?
   PROCEDURE, PUBLIC       :: finalize_atoms
!   FINAL :: finalize_atoms                        ! NOT YET IN GFORTRAN!
END TYPE cl_cryst
!
CONTAINS
!
!******************************************************************************
   SUBROUTINE alloc_arrays   ( this, natoms, nscat)
!
!  Allocate the arrays for "this" crystal 
!  Initialize all variables
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
   INTEGER,              INTENT(IN) :: natoms      ! Number of atoms for this crystal
   INTEGER,              INTENT(IN) :: nscat       ! Number of atom types for this crystal
!
   INTEGER  :: istatus
!
   IF(this%latom) THEN                             ! Crystal was allocated 
      DEALLOCATE ( this%cr_scat    , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfr   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfi   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_equ, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delf_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_at_equ  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_sav_atom, STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_at_lis  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_dw      , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%atoms      , STAT=istatus ) ! Always deallocate
   ENDIF
   ALLOCATE ( this%cr_scat    (9,0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_delfr   (  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_delfi   (  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_scat_int(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_scat_equ(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_delf_int(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_at_equ  (  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_sav_atom(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_at_lis  (  0:nscat), STAT=istatus ) ! Allocate atom names
   ALLOCATE ( this%cr_dw      (  0:nscat), STAT=istatus ) ! Allocate Debye Waller terms
   ALLOCATE ( this%atoms      (natoms   ), STAT=istatus ) ! Allocate list of atoms
!
   this%cr_natoms = natoms                            ! Store number of atoms
   this%cr_nscat  = nscat                             ! Store number of atom types
   this%latom  = .true.                            ! Flag that crystal is allocated
!
   END SUBROUTINE alloc_arrays
!******************************************************************************
   INTEGER FUNCTION get_natoms (this )
!
!  Return the number of atoms in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_natoms = this%cr_natoms
!
   END FUNCTION get_natoms
!******************************************************************************
   INTEGER FUNCTION get_nscat (this )
!
!  Return the number of atoms in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_nscat = this%cr_nscat
!
   END FUNCTION get_nscat
!******************************************************************************
   SUBROUTINE set_cryst_atom ( this, inum, itype, posit, iprop)
!
!  Fully define atom No. <inum> for "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this
   INTEGER,              INTENT(IN) :: inum
   INTEGER,              INTENT(IN) :: itype
   REAL   , DIMENSION(3),INTENT(IN) :: posit
   INTEGER,              INTENT(IN) :: iprop
!
   CALL this%atoms(inum)%set_atom ( itype, posit, iprop )
!
   END SUBROUTINE set_cryst_atom 
!******************************************************************************
   SUBROUTINE get_cryst_atom ( this, inum, itype, posit, iprop)
!
!  Return everything about atom <inum> in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                  :: this
   INTEGER,              INTENT(OUT) :: inum
   INTEGER,              INTENT(OUT) :: itype
   REAL   , DIMENSION(3),INTENT(OUT) :: posit
   INTEGER,              INTENT(OUT) :: iprop
!
   CALL this%atoms(inum)%get_atom ( itype, posit, iprop )
!
   END SUBROUTINE get_cryst_atom 
!******************************************************************************
   SUBROUTINE read_crystal ( this, infile)
!
!  Read a crystal structure from file
!  This procedure interfaces to the old "reastru" in "structur.f90"
!
   USE inter_readstru
   USE structur, ONLY: readstru
!
   IMPLICIT none
!
   CLASS (cl_cryst)                  :: this   ! The current crystal
   CHARACTER (LEN=*   ), INTENT(IN)  :: infile ! Disk file
   INTEGER                           :: inum   ! dummy index
   REAL   , DIMENSION(3)             :: posit  ! dummy position vector
   INTEGER                           :: istat  ! status variable
!
   rd_strucfile = infile
   rd_NMAX      = this%cr_natoms
   rd_MAXSCAT   = this%cr_nscat
!
!  Allocate temporary arrays
!
write(*,*) ALLOCATED(rd_cr_dw    )
write(*,*) ALLOCATED(rd_cr_at_lis)
   ALLOCATE ( rd_cr_dw    (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_at_lis(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_pos   (1:3,1:rd_NMAX)   , STAT = istat )
   ALLOCATE ( rd_cr_iscat (1:rd_NMAX)       , STAT = istat )
   ALLOCATE ( rd_cr_prop  (1:rd_NMAX)       , STAT = istat )
   ALLOCATE ( rd_as_at_lis(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_as_dw    (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_as_pos   (1:3,1:rd_MAXSCAT), STAT = istat )
   ALLOCATE ( rd_as_iscat (1:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_as_prop  (1:rd_MAXSCAT)    , STAT = istat )
write(*,*) 'STATUS ',istat
!
   rd_cr_dw   = 0.0
   rd_cr_at_lis = ' '
   rd_cr_pos    = 0.0
   rd_cr_iscat  = 0
   rd_cr_prop   = 0
   rd_as_at_lis = ' '
   rd_as_dw     = 0.0
   rd_as_pos    = 0.0
   rd_as_iscat  = 0
   rd_as_prop   = 0
!
!WRITE(*,*) ' 0   ', rd_cr_iscat(1)
!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1)
!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1),rd_cr_pos(2,1)
!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1),rd_cr_pos(2,1), rd_cr_pos(3,1)
!WRITE(*,*) ' 0   ', rd_cr_iscat(1), rd_cr_pos(1,1),rd_cr_pos(2,1), rd_cr_pos(3,1), rd_cr_prop(1)
!write(*,*) ' IN read_crystal ', rd_NMAX, rd_MAXSCAT
   CALL readstru (rd_NMAX, rd_MAXSCAT, rd_strucfile, rd_cr_name,        &
               rd_cr_spcgr, rd_cr_a0, rd_cr_win, rd_cr_natoms, rd_cr_nscat, rd_cr_dw,     &
               rd_cr_at_lis, rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_dim, rd_as_natoms, &
               rd_as_at_lis, rd_as_dw, rd_as_pos, rd_as_iscat, rd_as_prop, rd_sav_ncell,  &
               rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para)
write(*,*) ' BACK '
WRITE(*,*) ' nmax', rd_NMAX, rd_MAXSCAT
write(*,*) ' rd_cr_name  ', rd_cr_name
write(*,*) ' rd_cr_spcgr ', rd_cr_spcgr
write(*,*) ' rd_cr_a0    ', rd_cr_a0   
write(*,*) ' rd_cr_win   ', rd_cr_win  

write(*,*) ' spcgr_ianz   ',rd_spcgr_ianz
write(*,*) ' spcgr_para   ',rd_spcgr_para

write(*,*) ' cr_at_lis ', rd_cr_at_lis
WRITE(*,*) ' 1   ', rd_cr_iscat(1)
WRITE(*,*) ' 1   ', rd_cr_pos(1,1)
WRITE(*,*) ' 1   ', rd_cr_pos(2,1)
WRITE(*,*) ' 1   ', rd_cr_pos(3,1)
WRITE(*,*) ' 1   ', rd_cr_prop(1)
write(*,*)
   DO inum=1, rd_NMAX
     posit = rd_cr_pos(:,inum)
write(*,*) rd_cr_iscat(inum), posit, rd_cr_prop(inum)
     CALL this%atoms(inum)%set_atom ( rd_cr_iscat(inum), posit, rd_cr_prop(inum) )
   ENDDO
!
   DEALLOCATE ( rd_cr_dw    , STAT = istat )
   DEALLOCATE ( rd_cr_at_lis, STAT = istat )
   DEALLOCATE ( rd_cr_pos   , STAT = istat )
   DEALLOCATE ( rd_cr_iscat , STAT = istat )
   DEALLOCATE ( rd_cr_prop  , STAT = istat )
   DEALLOCATE ( rd_as_at_lis, STAT = istat )
   DEALLOCATE ( rd_as_dw    , STAT = istat )
   DEALLOCATE ( rd_as_pos   , STAT = istat )
   DEALLOCATE ( rd_as_iscat , STAT = istat )
   DEALLOCATE ( rd_as_prop  , STAT = istat )
!
   END SUBROUTINE read_crystal
!******************************************************************************
   SUBROUTINE set_crystal_from_standard   ( this, strucfile)
!
!  Set all values for the crystal to those of the main DISCUS crystal
!
   USE crystal_mod
   USE save_mod
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
   CHARACTER (LEN=*), INTENT(IN)    :: strucfile
!
   INTEGER               :: inum
   INTEGER               :: itype
   REAL   , DIMENSION(3) :: posit
   INTEGER               :: iprop
   
!
   this%cr_file         = strucfile
   this%cr_name         = cr_name
   this%cr_spcgr        = cr_spcgr
   this%cr_spcgrno      = cr_spcgrno
   this%cr_syst         = cr_syst
   this%cr_acentric     = cr_acentric
   this%cr_ncatoms      = cr_ncatoms
   this%cr_icc          = cr_icc
!
   this%spcgr_ianz      = spcgr_ianz
   this%spcgr_para      = spcgr_para
   this%cr_sel_prop     = cr_sel_prop
!
   this%cr_a0           = cr_a0
   this%cr_win          = cr_win
   this%cr_dim          = cr_dim
   this%cr_dim0         = cr_dim0
   this%cr_v            = cr_v
   this%cr_gten         = cr_gten
   this%cr_eps          = cr_eps
   this%cr_gmat         = cr_gmat
   this%cr_ar           = cr_ar
   this%cr_wrez         = cr_wrez
   this%cr_vr           = cr_vr
   this%cr_rten         = cr_rten
   this%cr_reps         = cr_reps
   this%cr_fmat         = cr_fmat
!
   this%cr_scat         = cr_scat
   this%cr_delfr        = cr_delfr
   this%cr_delfi        = cr_delfi
!
   this%cr_scat_int     = cr_scat_int
   this%cr_scat_equ     = cr_scat_equ
   this%cr_delf_int     = cr_delf_int
   this%cr_at_equ       = cr_at_equ
   this%cr_sav_atom     = sav_latom
   this%cr_newtype      = cr_newtype
   this%cr_cartesian    = cr_cartesian
!
   this%cr_natoms       = cr_natoms
   this%cr_nscat        = cr_nscat
   this%cr_n_REAL_atoms = cr_n_REAL_atoms

   this%cr_at_lis       = cr_at_lis
   this%cr_dw           = cr_dw
!
   DO inum=1,cr_natoms
      itype = cr_iscat(inum)
      posit = cr_pos(:,inum)
      iprop = cr_prop (inum)
      CALL this%atoms(inum)%set_atom ( itype, posit, iprop )
   ENDDO
!
   END SUBROUTINE set_crystal_from_standard
!******************************************************************************
   SUBROUTINE set_crystal_save_flags ( this, sav_scat, sav_adp, sav_gene, sav_symm, &
                                       sav_ncell, sav_obje, sav_doma, sav_mole)
!
   USE crystal_mod
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   LOGICAL, INTENT(IN)              ::  sav_scat
   LOGICAL, INTENT(IN)              ::  sav_adp
   LOGICAL, INTENT(IN)              ::  sav_gene
   LOGICAL, INTENT(IN)              ::  sav_symm
   LOGICAL, INTENT(IN)              ::  sav_ncell
   LOGICAL, INTENT(IN)              ::  sav_obje
   LOGICAL, INTENT(IN)              ::  sav_doma
   LOGICAL, INTENT(IN)              ::  sav_mole

   this%cr_sav_scat  = sav_scat
   this%cr_sav_adp   = sav_adp
   this%cr_sav_gene  = sav_gene
   this%cr_sav_symm  = sav_symm
   this%cr_sav_ncell = sav_ncell
   this%cr_sav_obje  = sav_obje
   this%cr_sav_doma  = sav_doma
   this%cr_sav_mole  = sav_mole
!
   END SUBROUTINE set_crystal_save_flags
!******************************************************************************
   SUBROUTINE get_standard_from_crystal   ( this)
!
!  Get all values from "this" crystal and place into the main DISCUS crystal!
!  The standard crystal is completely overwritten!
!
   USE crystal_mod
   USE save_mod    ! should be a read_mod !?!?
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER               :: inum
   INTEGER               :: itype
   REAL   , DIMENSION(3) :: posit
   INTEGER               :: iprop
   
!
   cr_name         = this%cr_name
   cr_spcgr        = this%cr_spcgr
   cr_spcgrno      = this%cr_spcgrno
   cr_syst         = this%cr_syst
   cr_acentric     = this%cr_acentric
   cr_ncatoms      = this%cr_ncatoms
   cr_icc          = this%cr_icc
!
   spcgr_ianz      = this%spcgr_ianz
   spcgr_para      = this%spcgr_para
   cr_sel_prop     = this%cr_sel_prop
!
   cr_a0           = this%cr_a0
   cr_win          = this%cr_win
   cr_dim          = this%cr_dim
   cr_dim0         = this%cr_dim0
   cr_v            = this%cr_v
   cr_gten         = this%cr_gten
   cr_eps          = this%cr_eps
   cr_gmat         = this%cr_gmat
   cr_ar           = this%cr_ar
   cr_wrez         = this%cr_wrez
   cr_vr           = this%cr_vr
   cr_rten         = this%cr_rten
   cr_reps         = this%cr_reps
   cr_fmat         = this%cr_fmat
!
   cr_scat         = this%cr_scat
   cr_delfr        = this%cr_delfr
   cr_delfi        = this%cr_delfi
!
   cr_scat_int     = this%cr_scat_int
   cr_scat_equ     = this%cr_scat_equ
   cr_delf_int     = this%cr_delf_int
   cr_at_equ       = this%cr_at_equ
   sav_latom       = this%cr_sav_atom
   cr_newtype      = this%cr_newtype
   cr_cartesian    = this%cr_cartesian
!
   cr_natoms       = this%cr_natoms
   cr_nscat        = this%cr_nscat
   cr_n_REAL_atoms = this%cr_n_REAL_atoms

   cr_at_lis       = this%cr_at_lis
   cr_dw           = this%cr_dw
!
   DO inum=1,cr_natoms
      CALL this%atoms(inum)%get_atom ( itype, posit, iprop )
      cr_iscat(inum) = itype
      cr_pos(:,inum) = posit
      cr_prop (inum) = iprop
   ENDDO
!
   END SUBROUTINE get_standard_from_crystal
!******************************************************************************
   LOGICAL FUNCTION is_alloc_atom (this)
!
   IMPLICIT NONE
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   is_alloc_atom = this%latom
!
   END FUNCTION is_alloc_atom 
!******************************************************************************
   SUBROUTINE finalize_atoms   ( this)
!
!  Finalize the class cl_cryst
!  Deallocate the array of atoms for "this" crystal 
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER  :: istatus
!
   IF(this%latom) THEN
      DEALLOCATE ( this%cr_scat    , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfr   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfi   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_equ, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delf_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_at_equ  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_sav_atom, STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_at_lis  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_dw      , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%atoms      , STAT=istatus )
   ENDIF
   this%cr_natoms = 0
   this%latom  = .false.
!
   END SUBROUTINE finalize_atoms
!******************************************************************************
END MODULE cryst_class
