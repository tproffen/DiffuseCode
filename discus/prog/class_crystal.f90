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
   INTEGER                              :: natoms  !the number of atoms
   INTEGER                              :: nscat   !the number of atom types
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cl_cr_at_lis  ! (  0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cl_cr_dw      ! (  0:MAXSCAT)
   TYPE(cl_atom), DIMENSION(:), POINTER :: atoms   !the actual atoms
   LOGICAL                              :: latom = .false.!have atoms been allocated?
CONTAINS
!
   PROCEDURE, PUBLIC, PASS :: alloc_atoms          ! allocate a set of atoms
   PROCEDURE, PUBLIC, PASS :: get_natoms           ! Return number of atoms in this crystal
   PROCEDURE, PUBLIC, PASS :: set_cryst_atom       ! Set iscat, posit, property
   PROCEDURE, PUBLIC, PASS :: get_cryst_atom       ! Get iscat, posit, property
   PROCEDURE, PUBLIC, PASS :: read_crystal         ! Read from disk file
   PROCEDURE, PUBLIC, PASS :: is_alloc_atom        ! Has atom array been allocated ?
   PROCEDURE, PUBLIC       :: finalize_atoms
!   FINAL :: finalize_atoms                        ! NOT YET IN GFORTRAN!
END TYPE cl_cryst
!
CONTAINS
!
!******************************************************************************
   SUBROUTINE alloc_atoms   ( this, natoms, nscat)
!
!  Allocate the array of atoms for "this" crystal 
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
      DEALLOCATE ( this%atoms,     STAT=istatus)   ! Always deallocate
      DEALLOCATE ( this%cl_cr_at_lis, STAT=istatus)   ! Always deallocate
      DEALLOCATE ( this%cl_cr_dw,     STAT=istatus)   ! Always deallocate
      DEALLOCATE ( this%atoms, STAT=istatus)       ! Always deallocate
   ENDIF
   ALLOCATE ( this%atoms(natoms),     STAT=istatus ) ! Create new crystal
   ALLOCATE ( this%cl_cr_at_lis(natoms), STAT=istatus ) ! Create new crystal
   ALLOCATE ( this%cl_cr_dw(natoms),     STAT=istatus ) ! Create new crystal
   this%natoms = natoms                            ! Store number of atoms
   this%nscat  = nscat                             ! Store number of atom types
   this%latom  = .true.                            ! Flag that crystal is allocated
!
   END SUBROUTINE alloc_atoms
!******************************************************************************
   INTEGER FUNCTION get_natoms (this )
!
!  Return the number of atoms in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_natoms = this%natoms
!
   END FUNCTION get_natoms
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
!  Return everything about atom >inum> in "this" crystal
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
   rd_NMAX      = this%natoms
   rd_MAXSCAT   = this%nscat
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
      DEALLOCATE ( this%atoms, STAT=istatus)
   ENDIF
   this%natoms = 0
   this%latom  = .false.
!
   END SUBROUTINE finalize_atoms
!******************************************************************************
END MODULE cryst_class
