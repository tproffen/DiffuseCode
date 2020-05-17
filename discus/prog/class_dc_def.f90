MODULE dc_def_class
!
!  Defines a class dc_def_class
USE prompt_mod
USE errlist_mod
USE precision_mod
!
!  This class contains the definitions for a molecule decoration.
!
IMPLICIT NONE
!
INTEGER, PARAMETER        :: DC_MAXMODE  = 6  ! We have six decoration modes
!
TYPE dc_con
   INTEGER, DIMENSION(:), ALLOCATABLE :: dc_con_surf   ! Surface atom type that is connected to molecule
   INTEGER                :: dc_con_mole   ! Molecule atom no that is connected to surface
   REAL                   :: dc_con_dist   ! Distance to surface atoms
   TYPE (dc_con), POINTER :: next          ! next connection entry
END TYPE dc_con
!
!  This is the basic definition TYPE:
TYPE dc_def
   INTEGER                  :: dc_def_index    ! Current definition number
   INTEGER                  :: dc_def_lname    ! Length of the definition name
   CHARACTER(LEN=PREC_STRING)      :: dc_def_name     ! the definition name
   INTEGER                  :: dc_def_lfile    ! Length of the molecule file name
   CHARACTER(LEN=PREC_STRING)      :: dc_def_file     ! the molecule file name
   INTEGER                  :: dc_def_molnum   = 0       ! This molecule in dc_molecules belongs to me
   INTEGER                  :: dc_def_type     = 1       ! Connection type (NORMAL, BRIDGE, ...)
   CHARACTER(LEN=8)         :: dc_def_ntype    = 'NORMAL'! Connection type name (NORMAL, BRIDGE, ...)
   INTEGER                  :: dc_def_ltype    = 6       ! Connection type name length (NORMAL, BRIDGE, ...)
   INTEGER,DIMENSION(1:2)   :: dc_def_axis     ! Molecule axis defined by two atoms
   INTEGER,DIMENSION(1:20)  :: dc_def_surfnew  ! Molecule surface atoms
   REAL                     :: dc_def_dens     = 0.100   ! Molecule density per square Angstroem
   LOGICAL                  :: dc_def_restrict = .FALSE. ! Restriction yes / no
   LOGICAL                  :: dc_def_l_form   = .FALSE. ! Single hkl or form
   INTEGER                  :: dc_def_n_hkl    = 0 ! Number of hkls that make special form
!  INTEGER,DIMENSION(1:3,2) :: dc_def_hkl      = 0 ! Surface restrictions
   INTEGER,ALLOCATABLE,DIMENSION(:,:) :: dc_def_hkl! Surface restrictions
   INTEGER                  :: dc_def_ncon     = 0 ! number of connections defined
   TYPE (dc_con), POINTER   :: dc_def_con      ! Chain of Connections associated to this definition
   TYPE (dc_def), POINTER   :: next            ! next definition entry
!
END TYPE dc_def
!
!
CHARACTER (LEN=8), DIMENSION(0:DC_MAXMODE) :: type_name
INTEGER          , DIMENSION(0:DC_MAXMODE) :: type_length
DATA type_name /'none def','normal  ','bridge  ','double','multiple','acceptor', 'donor'/
DATA type_length /8,6,6,6,8,8,5/
PUBLIC dc_find_def
PUBLIC dc_set_file
!
CONTAINS
!
!  Methods to work on the molecule decoration definitions
!
   RECURSIVE SUBROUTINE dc_find_def(this, search, temp_lname, temp_name, temp_id, lnew, success)
!
   TYPE (dc_def), POINTER :: this
   TYPE (dc_def), POINTER :: search
   INTEGER             , INTENT(IN)    :: temp_lname
   CHARACTER (LEN=*)   , INTENT(IN)    :: temp_name
   INTEGER             , INTENT(INOUT) :: temp_id
   LOGICAL             , INTENT(IN)    :: lnew
   INTEGER             , INTENT(INOUT) :: success
!
   IF (ASSOCIATED(this)) THEN
      IF ( this%dc_def_name(1:this%dc_def_lname) == temp_name(1:temp_lname) ) THEN
!        Found old entry
         search => this
         temp_id = this%dc_def_index
         ier_num = 0
         IF(lnew) THEN
            success = -1   ! Found old although told to make a new
         ELSE
            success =  0   ! Found old and told not to make a new
         ENDIF
      ELSE
         IF(ASSOCIATED(this%next)) THEN
!           Found a different node, keep searching
            temp_id = this%dc_def_index
            CALL dc_find_def(this%next,search,temp_lname, temp_name, temp_id, lnew, success)
         ELSE
            IF ( lnew ) THEN
!              Node not found, user wants a new node
               ALLOCATE(this%next   )
               this%next%dc_def_index = this%dc_def_index + 1
               this%next%dc_def_name  = temp_name
               this%next%dc_def_lname = temp_lname
               this%next%dc_def_type  = 0
               this%next%dc_def_lfile  = 0
               this%next%dc_def_file   = ' '
               this%next%dc_def_axis(:)= -1
               this%next%dc_def_surfnew(:)= 0
               NULLIFY(this%next%next)
               NULLIFY(this%next%dc_def_con)
               this%next%dc_def_ncon = 0
!              CALL dc_set_con_default(this%next%dc_def_con)
               search => this%next
               temp_id = this%dc_def_index
               success =  0   ! Did not find old and told to make a new
            ELSE
!              Node not found, user does NOT want a new node
               success =  -1  ! Did not find old and told not to make a new
               ier_num = -1116
            ENDIF
         ENDIF
      ENDIF
   ELSE
      IF ( lnew ) THEN
!        Node not found, user wants a new node
         ALLOCATE(this)
         this%dc_def_index  = temp_id + 1
         this%dc_def_name   = temp_name
         this%dc_def_lname  = temp_lname
         this%dc_def_type   = 0
         this%dc_def_lfile  = 0
         this%dc_def_file   = ' '
         this%dc_def_axis(:)= -1
         this%dc_def_surfnew(:)= 0
         NULLIFY(this%next)
         NULLIFY(this%dc_def_con)
         this%dc_def_ncon = 0
!        CALL dc_set_con_default(this%dc_def_con)
         search => this
         temp_id = this%dc_def_index
         success =  0   ! Did not find old and told to make a new
      ELSE
!        Node not found, user does NOT want a new node
         success =  -1  ! Did not find old and told not to make a new
         ier_num = -1116
      ENDIF
   ENDIF
   END SUBROUTINE dc_find_def
!
!*******************************************************************************
!
   SUBROUTINE dc_set_file(this, dc_temp_lfile, dc_temp_file)
!
   TYPE (dc_def), POINTER :: this
   INTEGER             , INTENT(IN) :: dc_temp_lfile
   CHARACTER (LEN=*   ), INTENT(IN) :: dc_temp_file
!
   IF(ASSOCIATED(this)) THEN
      this%dc_def_file  = dc_temp_file
      this%dc_def_lfile = dc_temp_lfile
   ELSE
      write(*,*) ' SET FILE, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_file
!
!*******************************************************************************
!
   SUBROUTINE dc_set_axis(this, dc_temp_axis)
!
   TYPE (dc_def), POINTER :: this
   INTEGER, DIMENSION(1:2), INTENT(IN) :: dc_temp_axis
!
   IF(ASSOCIATED(this)) THEN
      this%dc_def_axis(:)  = dc_temp_axis(:)
   ELSE
      write(*,*) ' SET AXIS, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_axis
!
!*******************************************************************************
!
   SUBROUTINE dc_set_surfnew(this, dc_temp_surfnew)
!
   TYPE (dc_def), POINTER :: this
   INTEGER, DIMENSION(1:20), INTENT(IN) :: dc_temp_surfnew
!
   IF(ASSOCIATED(this)) THEN
      this%dc_def_surfnew(:)  = dc_temp_surfnew(:)
   ELSE
      write(*,*) ' SET SURFNEW, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_surfnew
!
!*******************************************************************************
!
   SUBROUTINE dc_set_molnum(this, dc_temp_molnum)
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(IN)    :: dc_temp_molnum
!
   IF(ASSOCIATED(this)) THEN
      this%dc_def_molnum  = dc_temp_molnum
   ELSE
      write(*,*) ' SET molnum, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_molnum
!
!*******************************************************************************
!
   SUBROUTINE dc_set_dens(this, dc_temp_dens)
!
   TYPE (dc_def), POINTER :: this
   REAL         , INTENT(IN) :: dc_temp_dens
!
   IF(ASSOCIATED(this)) THEN
      this%dc_def_dens  = dc_temp_dens
   ELSE
      write(*,*) ' SET DENS, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_dens
!
!*******************************************************************************
!
   SUBROUTINE dc_set_type(this, dc_temp_type)
!
   TYPE (dc_def), POINTER :: this
   INTEGER             , INTENT(IN) :: dc_temp_type
!
   IF(ASSOCIATED(this)) THEN
      this%dc_def_type  = dc_temp_type
      this%dc_def_ntype = type_name  (dc_temp_type)
      this%dc_def_ltype = type_length(dc_temp_type)
   ELSE
      write(*,*) ' SET TYPE, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_type
!
!*******************************************************************************
!
   SUBROUTINE dc_set_hkl(this, dc_temp_restrict, dc_temp_n_hkl, dc_temp_hkl, l_form)
!
   TYPE (dc_def), POINTER :: this
   LOGICAL                   , INTENT(IN) :: dc_temp_restrict
   INTEGER                   , INTENT(IN) :: dc_temp_n_hkl
   INTEGER, DIMENSION(:,:)   , INTENT(IN) :: dc_temp_hkl
   LOGICAL                   , INTENT(IN) :: l_form
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: temp_hkl
   INTEGER :: nold
!
   nold = MIN(0,this%dc_def_n_hkl)
   IF(ASSOCIATED(this)) THEN
      IF(ALLOCATED(this%dc_def_hkl)) THEN
         nold = MIN(this%dc_def_n_hkl, UBOUND(this%dc_def_hkl,2))
         IF(UBOUND(this%dc_def_hkl,2)<dc_temp_n_hkl) THEN
            ALLOCATE(temp_hkl(1:3,1:dc_temp_n_hkl + 5))
            temp_hkl = 0
            temp_hkl(:,1:this%dc_def_n_hkl) = this%dc_def_hkl(:,1:this%dc_def_n_hkl)
            DEALLOCATE(this%dc_def_hkl)
            ALLOCATE  (this%dc_def_hkl(1:3,1:dc_temp_n_hkl + 5))
            this%dc_def_hkl(:,:) = temp_hkl(:,:)
            DEALLOCATE(temp_hkl)
         ENDIF
      ELSE
         ALLOCATE(this%dc_def_hkl(1:3,1:dc_temp_n_hkl + 5))
         this%dc_def_hkl = 0
      ENDIF
      this%dc_def_restrict   = dc_temp_restrict
      this%dc_def_l_form     = l_form
      this%dc_def_n_hkl      = dc_temp_n_hkl
      this%dc_def_hkl(:,nold+1:dc_temp_n_hkl) = dc_temp_hkl(:,nold+1:dc_temp_n_hkl)
   ELSE
      write(*,*) ' SET HKL, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_hkl
!
!*******************************************************************************
!
   SUBROUTINE dc_set_connection(this)
!
   TYPE (dc_def), POINTER :: this
!
   IF(ASSOCIATED(this)) THEN
!     this%dc_def_file  = dc_temp_file
!     this%dc_def_lfile = dc_temp_lfile
   ELSE
      write(*,*) ' SET FILE, NODE NOT associated!'
   ENDIF
!
   END SUBROUTINE dc_set_connection
! 
!*******************************************************************************
!
   RECURSIVE SUBROUTINE dc_show_def (this, ier_num )
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(OUT)   :: ier_num
   INTEGER :: i
!
   IF ( ASSOCIATED(this)) THEN
      WRITE(output_io,*)
      WRITE(output_io,1100) this%dc_def_index,this%dc_def_name(1:this%dc_def_lname)
      WRITE(output_io,1200) this%dc_def_file(1:this%dc_def_lfile)
      IF(this%dc_def_axis(1) == -1) THEN
         WRITE(output_io,1300) this%dc_def_type, this%dc_def_ntype(1:this%dc_def_ltype)
      ELSE
         IF(this%dc_def_axis(2) == -1) THEN
            WRITE(output_io,1310) this%dc_def_type,  &
               this%dc_def_ntype(1:this%dc_def_ltype),  &
               this%dc_def_axis(1)
         ELSE
            WRITE(output_io,1320) this%dc_def_type,  &
               this%dc_def_ntype(1:this%dc_def_ltype),  &
               this%dc_def_axis(1), this%dc_def_axis(2)
         ENDIF
         IF(.NOT.this%dc_def_restrict) THEN
            WRITE(output_io, 1330)
         ELSE
            IF(this%dc_def_l_form) THEN
               WRITE(output_io, 1340)
            ELSE
               WRITE(output_io, 1345)
            ENDIF
            DO i=1, this%dc_def_n_hkl
               WRITE(output_io, 1350) this%dc_def_hkl(:,i)
            ENDDO
         ENDIF
      ENDIF
      WRITE(output_io, 1500 ) this%dc_def_dens
      WRITE(output_io, 1400 ) this%dc_def_ncon
      CALL dc_show_con(this%dc_def_con)
      IF ( ASSOCIATED(this%next)) THEN
         CALL dc_show_def(this%next, ier_num)
      ENDIF
   ENDIF
!
1100 FORMAT(' Definition        :',i4,' ',4a)
1200 FORMAT('   Molecule file   :     ',a)
1300 FORMAT('   Connection type :     ',i4,' ',a8)
1310 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand axis ',i4,' to last')
1320 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand axis ',i4,' to ',i4 )
1330 FORMAT('   Connection open to all surfaces')
1340 FORMAT('   Connection restricted to following forms')
1345 FORMAT('   Connection restricted to following surfaces')
1350 FORMAT('              H K L:     ',2(i4,2x),i4) 
1500 FORMAT('   Surface density :     ',f6.3,' Ligands/A^2')
1400 FORMAT('   No. of Connect. :     ',i4)
   END SUBROUTINE dc_show_def
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE dc_reset_def (this)
!
!  Remove the current definitions
!
   TYPE (dc_def), POINTER :: this
   TYPE (dc_con), POINTER :: connection
!
   reset_loop: DO 
      IF ( ASSOCIATED(this)) THEN
         IF(ASSOCIATED(this%next)) CALL dc_reset_def( this%next)
         connection => this%dc_def_con
         CALL dc_reset_con(connection)
         IF(ALLOCATED(this%dc_def_hkl)) DEALLOCATE(this%dc_def_hkl)
         DEALLOCATE(this)
      ELSE
         EXIT reset_loop
      ENDIF
   ENDDO reset_loop
   END SUBROUTINE dc_reset_def
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE dc_reset_con (this)
!
!  Remove the current definitions
!
   TYPE (dc_con), POINTER :: this
!
   reset_loop: DO 
      IF ( ASSOCIATED(this)) THEN
         IF(ASSOCIATED(this%next)) CALL dc_reset_con( this%next)
         IF(ALLOCATED(this%dc_con_surf)) DEALLOCATE(this%dc_con_surf)
         DEALLOCATE(this)
      ELSE
         EXIT reset_loop
      ENDIF
   ENDDO reset_loop
   END SUBROUTINE dc_reset_con
!
!*******************************************************************************
!
   SUBROUTINE dc_set_con_default(this)
!
   IMPLICIT NONE
!
   TYPE (dc_con), POINTER :: this
!
   ALLOCATE(this%dc_con_surf(0:1))
   this%dc_con_surf(0)  =  0
   this%dc_con_surf(1)  = -1
   this%dc_con_mole  =  0
   this%dc_con_dist  = -1.0
!  this%dc_con_restrict = .FALSE.
!  this%dc_con_l_form   = .FALSE.
!  this%dc_con_n_hkl    = 1
!  this%dc_con_hkl      = 0
   NULLIFY(this%next)
!
   END SUBROUTINE dc_set_con_default
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE dc_set_con(this, surf, mole, dist)
!
   IMPLICIT NONE
!
   TYPE (dc_con), POINTER :: this
   INTEGER, DIMENSION(0:4), INTENT(IN) :: surf
   INTEGER, INTENT(IN) :: mole
   REAL   , INTENT(IN) :: dist
!
   REAL, PARAMETER     :: EPS = 1.E-6
!
   INTEGER             :: i
!
   IF(ASSOCIATED(this)) THEN
         IF(this%dc_con_surf(1)==surf(1) .and. this%dc_con_mole==mole .and. &
            ABS(this%dc_con_dist-dist) < EPS) THEN
            ier_num = -1117
         ELSE
            CALL dc_set_con(this%next, surf, mole, dist)
         ENDIF
   ELSE
         ALLOCATE(this)
         ALLOCATE(this%dc_con_surf(0:surf(0)))
         DO i= 1,surf(0)
            this%dc_con_surf(i)  = surf(i)
         ENDDO
         this%dc_con_surf(0)  = surf(0)
         this%dc_con_mole  = mole
         this%dc_con_dist  = dist
         NULLIFY(this%next)
         ier_num = 0
   ENDIF
!
   END SUBROUTINE dc_set_con
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE dc_show_con(this)
!
   USE crystal_mod
   USE atom_name
!   USE atom_env_mod
   IMPLICIT NONE
!
   TYPE (dc_con), POINTER :: this
!
   CHARACTER(LEN=*   ) :: line
   INTEGER             :: i
!
   IF(ASSOCIATED(this)) THEN
      line = ' '
      DO i=1,this%dc_con_surf(0)
         line((i-1)*9+1:(i-1)*9+9) = at_name (this%dc_con_surf(i)) //' '
      ENDDO
      WRITE(output_io, 1000) line(1:(this%dc_con_surf(0)*10)), this%dc_con_mole, this%dc_con_dist
      IF(ASSOCIATED(this%next)) THEN
         CALL dc_show_con(this%next)
      ENDIF
   ENDIF
!
1000 FORMAT('   Type , lig, dist: ',a,1x,i4,1x,f8.5)
   END SUBROUTINE dc_show_con
!
!*******************************************************************************
!
   SUBROUTINE dc_get_con(this, surf, mole, dist)
!
   IMPLICIT NONE
!
   TYPE (dc_con), POINTER :: this
   INTEGER, DIMENSION(0:4), INTENT(OUT) :: surf
   INTEGER, INTENT(OUT)   :: mole
   REAL   , INTENT(OUT)   :: dist
!
   INTEGER                :: i
!
   IF(ASSOCIATED(this)) THEN
      DO i= 0,this%dc_con_surf(0)
         surf(i) = this%dc_con_surf(i)
      ENDDO
      surf(0) = this%dc_con_surf(0)
      mole = this%dc_con_mole
      dist = this%dc_con_dist
   ENDIF
!
   END SUBROUTINE dc_get_con
!
!*******************************************************************************
!
   SUBROUTINE dc_get_id(this, dc_temp_id)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(OUT)   :: dc_temp_id
!
   dc_temp_id = this%dc_def_index
!
   END SUBROUTINE dc_get_id
!
!*******************************************************************************
!
   SUBROUTINE dc_get_dens(this, dc_temp_dens)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   REAL   , INTENT(OUT)   :: dc_temp_dens
!
   dc_temp_dens = this%dc_def_dens 
!
   END SUBROUTINE dc_get_dens
!
!*******************************************************************************
!
   SUBROUTINE dc_get_molnum(this, dc_temp_molnum)
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(OUT)   :: dc_temp_molnum
!
   dc_temp_molnum = this%dc_def_molnum
!
   END SUBROUTINE dc_get_molnum
!
!*******************************************************************************
!
   SUBROUTINE dc_get_type(this, dc_temp_type)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(OUT)   :: dc_temp_type
!
   dc_temp_type = this%dc_def_type 
!
   END SUBROUTINE dc_get_type
!
!*******************************************************************************
!
   SUBROUTINE dc_get_axis(this, axis)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   INTEGER, DIMENSION(1:2), INTENT(OUT) :: axis
!
   axis(:) = this%dc_def_axis(:) 
!
   END SUBROUTINE dc_get_axis
!
!*******************************************************************************
!
   SUBROUTINE dc_get_surfnew(this, surfnew)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   INTEGER, DIMENSION(1:20), INTENT(OUT) :: surfnew
!
   surfnew(:) = this%dc_def_surfnew(:) 
!
   END SUBROUTINE dc_get_surfnew
!
!*******************************************************************************
!
   SUBROUTINE dc_get_mole_name(this, mole_name, length)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   CHARACTER (LEN=*   ), INTENT(OUT)   :: mole_name
   INTEGER, INTENT(OUT)   :: length
!
   mole_name = this%dc_def_file 
   length    = this%dc_def_lfile
!
   END SUBROUTINE dc_get_mole_name
!
!*******************************************************************************
!
   SUBROUTINE dc_get_ncon(this, ncon)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(OUT)   :: ncon
!
   ncon = this%dc_def_ncon
!
   END SUBROUTINE dc_get_ncon
!
!*******************************************************************************
!
   SUBROUTINE dc_inc_ncon(this,ncon)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(out)   :: ncon
!
   this%dc_def_ncon = this%dc_def_ncon + 1 
   ncon = this%dc_def_ncon
!
   END SUBROUTINE dc_inc_ncon
!
END MODULE dc_def_class
