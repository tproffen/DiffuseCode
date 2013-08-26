MODULE dc_def_class
!
!  Defines a class dc_def_class
!
!  This class contains the definitions for a molecule decoration.
!
IMPLICIT NONE
!
!PUBLIC dc_def
!
TYPE dc_con
   INTEGER                :: dc_con_lname  ! Length of the definition name
   CHARACTER(LEN=1024)    :: dc_con_name   ! the definition name
   INTEGER                :: dc_con_mole   ! Molecule atom no that is connected to surface
   REAL                   :: dc_con_dist   ! Distance to surface atoms
   INTEGER                :: dc_con_surf   ! Surface atom no that is connected to molecule
   INTEGER                :: dc_con_numb   ! Number of surface atoms that are connected to molecule
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
PUBLIC dc_find_def
PUBLIC dc_set_file
!
CONTAINS
!
!  Methods to work on the molecule decoration definitions
!
   RECURSIVE SUBROUTINE dc_find_def(this, search, temp_lname, temp_name, temp_id, lnew, ier_num)
!
   TYPE (dc_def), POINTER :: this
   TYPE (dc_def), POINTER :: search
   INTEGER             , INTENT(IN)  :: temp_lname
   CHARACTER (LEN=1024), INTENT(IN)  :: temp_name
   INTEGER             , INTENT(INOUT) :: temp_id
   LOGICAL             , INTENT(IN)  :: lnew
   INTEGER             , INTENT(OUT) :: ier_num
!
write(*,*) ' TEMP_ID ', temp_id
   IF (ASSOCIATED(this)) THEN
      IF ( this%dc_def_name(1:this%dc_def_lname) == temp_name(1:temp_lname) ) THEN
         search => this
         ier_num = 0
         RETURN
      ELSE
         IF(ASSOCIATED(this%next)) THEN
            temp_id = this%dc_def_index
            CALL dc_find_def(this%next,search,temp_lname, temp_name, temp_id, lnew, ier_num)
         ELSE
write(*,*) ' Next does not exist, make new ',lnew
            IF ( lnew ) THEN
               ALLOCATE(this%next   )
               this%next%dc_def_index = this%dc_def_index + 1
               this%next%dc_def_name  = temp_name
               this%next%dc_def_lname = temp_lname
               NULLIFY(this%next%next)
               search => this%next
            ELSE
               ier_num = -116
            ENDIF
         ENDIF
      ENDIF
   ELSE
write(*,*) ' NODE does not exist, make new ',lnew
      IF ( lnew ) THEN
         ALLOCATE(this)
         this%dc_def_index = temp_id + 1
         this%dc_def_name  = temp_name
         this%dc_def_lname = temp_lname
         NULLIFY(this%next)
         search => this
      ELSE
         ier_num = -116
      ENDIF
   ENDIF
   END SUBROUTINE dc_find_def
!
!
   RECURSIVE SUBROUTINE dc_find_con(this, search, temp_lname, temp_name, temp_id, lnew, ier_num)
!
   TYPE (dc_con), POINTER :: this
   TYPE (dc_con), POINTER :: search
   INTEGER             , INTENT(IN)  :: temp_lname
   CHARACTER (LEN=1024), INTENT(IN)  :: temp_name
   INTEGER             , INTENT(INOUT) :: temp_id
   LOGICAL             , INTENT(IN)  :: lnew
   INTEGER             , INTENT(OUT) :: ier_num
!
write(*,*) ' TEMP_ID ', temp_id
   IF (ASSOCIATED(this)) THEN
      IF ( this%dc_con_name(1:this%dc_con_lname) == temp_name(1:temp_lname) ) THEN
         search => this
         ier_num = 0
         RETURN
      ELSE
         IF(ASSOCIATED(this%next)) THEN
            CALL dc_find_con(this%next,search,temp_lname, temp_name, temp_id, lnew, ier_num)
         ELSE
write(*,*) ' Next does not exist, make new ',lnew
            IF ( lnew ) THEN
               ALLOCATE(this%next   )
               this%next%dc_con_name  = temp_name
               this%next%dc_con_lname = temp_lname
               NULLIFY(this%next%next)
               search => this%next
            ELSE
               ier_num = -116
            ENDIF
         ENDIF
      ENDIF
   ELSE
write(*,*) ' NODE does not exist, make new ',lnew
      IF ( lnew ) THEN
         ALLOCATE(this)
         this%dc_con_name  = temp_name
         this%dc_con_lname = temp_lname
         NULLIFY(this%next)
         search => this
      ELSE
         ier_num = -116
      ENDIF
   ENDIF
   END SUBROUTINE dc_find_con
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
!
   SUBROUTINE dc_set_connection(this, MAXP, ianz, cpara, lpara, ier_num)
!
   TYPE (dc_def), POINTER :: this
   INTEGER             , INTENT(IN) :: MAXP 
   INTEGER             , INTENT(IN) :: ianz 
   INTEGER             , DIMENSION(MAXP), INTENT(IN) :: lpara
   CHARACTER (LEN=1024), DIMENSION(MAXP), INTENT(IN) :: cpara
   INTEGER             , INTENT(OUT) :: ier_num 
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
      WRITE(*,1100) this%dc_def_index,this%dc_def_name(1:this%dc_def_lname)
      WRITE(*,1200) this%dc_def_file(1:this%dc_def_lfile)
      IF ( ASSOCIATED(this%next)) THEN
         CALL dc_show_def(this%next, ier_num)
      ENDIF
   ENDIF
!
1100 FORMAT(' Definition     :',i4,' ',4a)
1200 FORMAT('   Molecule file:     ',a)
   END SUBROUTINE dc_show_def
!
!
   SUBROUTINE dc_reset_def (this,ier_num)
!
!  Remove the current definitions
!
   TYPE (dc_def), POINTER :: this
   TYPE (dc_def), POINTER :: search
   INTEGER, INTENT(OUT)   :: ier_num
!
   reset_loop: DO 
      IF ( ASSOCIATED(this)) THEN
         SEARCH => this%next
         DEALLOCATE(this)
      ELSE
         EXIT reset_loop
      ENDIF
   ENDDO reset_loop
   END SUBROUTINE dc_reset_def
END MODULE dc_def_class
