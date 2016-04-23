MODULE dc_def_class
!
!  Defines a class dc_def_class
USE prompt_mod
USE errlist_mod
!
!  This class contains the definitions for a molecule decoration.
!
IMPLICIT NONE
!
!PUBLIC dc_def
!
TYPE dc_con
!  INTEGER                :: dc_con_lname  ! Length of the definition name
!  CHARACTER(LEN=1024)    :: dc_con_name   ! the definition name
   INTEGER, DIMENSION(:), ALLOCATABLE :: dc_con_surf   ! Surface atom type that is connected to molecule
   INTEGER                :: dc_con_mole   ! Molecule atom no that is connected to surface
!   INTEGER                :: dc_con_numb   ! Number of surface atoms that are connected to molecule
   REAL                   :: dc_con_dist   ! Distance to surface atoms
   TYPE (dc_con), POINTER :: next          ! next connection entry
END TYPE dc_con
!
!  This is the basic definition TYPE:
TYPE dc_def
   INTEGER                :: dc_def_index    ! Current definition number
   INTEGER                :: dc_def_lname    ! Length of the definition name
   CHARACTER(LEN=1024)    :: dc_def_name     ! the definition name
   INTEGER                :: dc_def_lfile    ! Length of the molecule file name
   CHARACTER(LEN=1024)    :: dc_def_file     ! the molecule file name
   INTEGER                :: dc_def_type     ! Connection type (NORMAL, BRIDGE, ...)
   CHARACTER(LEN=8)       :: dc_def_ntype    ! Connection type name (NORMAL, BRIDGE, ...)
   INTEGER                :: dc_def_ltype    ! Connection type name length (NORMAL, BRIDGE, ...)
   INTEGER,DIMENSION(1:2) :: dc_def_axis     ! Molecule axis defined by two atoms
   INTEGER                :: dc_def_ncon     ! number of connections defined
   TYPE (dc_con), POINTER :: dc_def_con      ! Chain of Connections associated to this definition
   TYPE (dc_def), POINTER :: next            ! next definition entry
!
!   CONTAINS
!   PROCEDURE, PUBLIC, PASS :: dc_find_def 
!   PROCEDURE, PUBLIC, PASS :: dc_set_file
END TYPE dc_def
!
!TYPE (dc_def), POINTER    :: dc_def_head => NULL()
!TYPE (dc_def), POINTER    :: dc_def_tail => NULL()
!TYPE (dc_def), POINTER    :: dc_def_temp => NULL()
!
CHARACTER (LEN=8), DIMENSION(0:4) :: type_name
INTEGER          , DIMENSION(0:4) :: type_length
DATA type_name /'none def','normal  ','bridge  ','double','multiple'/
DATA type_length /8,6,6,6,8/
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
   CHARACTER (LEN=1024), INTENT(IN)    :: temp_name
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
               this%next%dc_def_axis(:)=-1.0
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
         this%dc_def_axis(:)=-1.0
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
!
!   RECURSIVE SUBROUTINE dc_find_con(this, search, temp_lname, temp_name, temp_id, lnew, ier_num)
!
!   TYPE (dc_con), POINTER :: this
!   TYPE (dc_con), POINTER :: search
!   INTEGER             , INTENT(IN)  :: temp_lname
!   CHARACTER (LEN=1024), INTENT(IN)  :: temp_name
!   INTEGER             , INTENT(INOUT) :: temp_id
!   LOGICAL             , INTENT(IN)  :: lnew
!   INTEGER             , INTENT(OUT) :: ier_num
!
!write(*,*) ' TEMP_ID ', temp_id
!   IF (ASSOCIATED(this)) THEN
!      IF ( this%dc_con_name(1:this%dc_con_lname) == temp_name(1:temp_lname) ) THEN
!         search => this
!         ier_num = 0
!         RETURN
!      ELSE
!         IF(ASSOCIATED(this%next)) THEN
!            CALL dc_find_con(this%next,search,temp_lname, temp_name, temp_id, lnew, ier_num)
!         ELSE
!write(*,*) ' Next does not exist, make new ',lnew
!            IF ( lnew ) THEN
!               ALLOCATE(this%next   )
!               this%next%dc_con_name  = temp_name
!               this%next%dc_con_lname = temp_lname
!               NULLIFY(this%next%next)
!               search => this%next
!            ELSE
!               ier_num = -116
!            ENDIF
!         ENDIF
!      ENDIF
!   ELSE
!write(*,*) ' NODE does not exist, make new ',lnew
!      IF ( lnew ) THEN
!         ALLOCATE(this)
!         this%dc_con_name  = temp_name
!         this%dc_con_lname = temp_lname
!         NULLIFY(this%next)
!         search => this
!      ELSE
!         ier_num = -116
!      ENDIF
!   ENDIF
!   END SUBROUTINE dc_find_con
!
   SUBROUTINE dc_set_file(this, dc_temp_lfile, dc_temp_file)
!
   TYPE (dc_def), POINTER :: this
   INTEGER             , INTENT(IN) :: dc_temp_lfile
   CHARACTER (LEN=1024), INTENT(IN) :: dc_temp_file
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
   RECURSIVE SUBROUTINE dc_show_def (this, ier_num )
!
   TYPE (dc_def), POINTER :: this
   INTEGER, INTENT(OUT)   :: ier_num
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
      ENDIF
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
1400 FORMAT('   No. of Connect. :     ',i4)
   END SUBROUTINE dc_show_def
!
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
         DEALLOCATE(this)
      ELSE
         EXIT reset_loop
      ENDIF
   ENDDO reset_loop
   END SUBROUTINE dc_reset_def
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
         DEALLOCATE(this)
      ELSE
         EXIT reset_loop
      ENDIF
   ENDDO reset_loop
   END SUBROUTINE dc_reset_con
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
   NULLIFY(this%next)
!
   END SUBROUTINE dc_set_con_default
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
   RECURSIVE SUBROUTINE dc_show_con(this)
!
   USE crystal_mod
   USE atom_name
!   USE atom_env_mod
   IMPLICIT NONE
!
   TYPE (dc_con), POINTER :: this
!
   CHARACTER(9) at_name_d
   CHARACTER(LEN=1024) :: line
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
   SUBROUTINE dc_get_con(this, surf, mole, dist)
!
   IMPLICIT NONE
!
   TYPE (dc_con), POINTER :: this
   INTEGER, DIMENSION(0:4), INTENT(OUT) :: surf
   INTEGER, INTENT(OUT)   :: mole
   REAL   , INTENT(OUT)   :: dist
!
   REAL, PARAMETER        :: EPS = 1.E-6
!
   INTEGER                :: i
!
   IF(ASSOCIATED(this)) THEN
      DO i= 0,this%dc_con_surf(0)
         surf(i) = this%dc_con_surf(i)
      ENDDO
      mole = this%dc_con_mole
      dist = this%dc_con_dist
   ENDIF
!
   END SUBROUTINE dc_get_con
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
   SUBROUTINE dc_get_mole_name(this, mole_name, length)
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER :: this
   CHARACTER (LEN=1024), INTENT(OUT)   :: mole_name
   INTEGER, INTENT(OUT)   :: length
!
   mole_name = this%dc_def_file 
   length    = this%dc_def_lfile
!
   END SUBROUTINE dc_get_mole_name
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
