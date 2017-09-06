MODULE class_macro_internal
!
USE macro_class
!
IMPLICIT NONE
!
PRIVATE
PUBLIC  :: macro_internal
PUBLIC  :: macro_root, macro_temp, macro_read_temp
PUBLIC  :: macro_add_node, macro_find_node, macro_write_node
PUBLIC  :: mac_tree_root, mac_tree_temp, mac_tree_active, mac_tree_tail
PUBLIC  :: macro_level
PUBLIC  :: lmakro, lmakro_error
PUBLIC  :: sprompt
!

TYPE :: macro_internal
   INTEGER                           :: number
   CHARACTER (LEN=1024)              :: macrofile   ! "file" name for this internal storage
   TYPE (cl_macro)        , POINTER  :: macros      ! The actual data structure for the macros 
   TYPE (macro_internal), POINTER    :: before      ! Pointer to prev stored macros 
   TYPE (macro_internal), POINTER    :: after       ! Pointer to next stored macros 
END TYPE macro_internal
!
TYPE(macro_internal), POINTER        :: macro_root  ! The main root of the internal storage tree
TYPE(macro_internal), POINTER        :: macro_temp  ! a temporary pointer to build the tree
TYPE(macro_internal), POINTER        ::  macro_read_temp  ! a temporary pointer to migrate the tree
!
TYPE :: macro_tree                                      ! Indicates which macro we work on
   INTEGER                             :: nparams=0     ! number of parameters
   CHARACTER (LEN= 256),DIMENSION(0:20):: params =' '   ! parameters
   INTEGER             ,DIMENSION(0:20):: lparams=0     ! parameter lengths
   INTEGER                           :: current  =0     ! Current line
   INTEGER                           :: level    =0     ! Current macro depth
   TYPE (macro_tree), POINTER        :: parent          ! Parent of current macro
   TYPE (macro_tree), POINTER        :: kid             ! kid to the current macro
   TYPE (macro_internal), POINTER    :: active          ! The associated storage node
END TYPE macro_tree
!
TYPE(macro_tree), POINTER            :: mac_tree_root   ! Main root for execution tree
TYPE(macro_tree), POINTER            :: mac_tree_temp   ! temporary pointer to build the tree
TYPE(macro_tree), POINTER            :: mac_tree_active ! temporary pointer that migrates the tree
TYPE(macro_tree), POINTER            :: mac_tree_tail   ! temporary pointer that migrates the tree
!
INTEGER                              :: macro_level = 0 ! Current macro depth
LOGICAL                              :: lmakro = .false.! Macro status is on/off
LOGICAL                              :: lmakro_error = .false.! Macro termination with status is off
CHARACTER(LEN=40)                    :: sprompt         ! Prompt at macro start
!
CONTAINS
!*******************************************************************************
   RECURSIVE SUBROUTINE macro_add_node ( ptr, new_node )
!
   IMPLICIT NONE
!
   TYPE(macro_internal), POINTER :: ptr
   TYPE(macro_internal), POINTER :: new_node
   TYPE(macro_internal), POINTER :: temp
!
   IF ( .not. ASSOCIATED(ptr))  THEN                     ! Pointer does not exist
      ptr  => new_node                                   ! Add here
   ELSEIF ( LLT(new_node%macrofile, ptr%macrofile)) THEN ! new "macrofile" is < old
      IF ( ASSOCIATED(ptr%before) ) THEN                 ! before node exists, 
         CALL macro_add_node(ptr%before, new_node)       !   recursively add new node
      ELSE                                               ! before node does not exist
         ptr%before  => new_node                         ! Add new node here
      ENDIF                                              ! 
   ELSEIF ( new_node%macrofile  ==   ptr%macrofile) THEN ! New "macrofile" is = old
      CALL ptr%macros%finalize_macro                     ! Clean up old macro
      new_node%before => ptr%before                      ! Retain old before
      new_node%after  => ptr%after                       !        and after 
      temp => ptr                                        ! Retain address to node
      ptr => new_node                                    ! Place new node instead of old
      DEALLOCATE(temp)                                   ! delete old node
!     new_node => ptr                                    ! new node points to old
   ELSEIF ( LGT(new_node%macrofile, ptr%macrofile)) THEN ! New "macrofile" is > old
      IF ( ASSOCIATED(ptr%after) ) THEN                  ! after node exists
         CALL macro_add_node(ptr%after, new_node)        !    recursively add new node
      ELSE                                               ! afet node does not exist
         ptr%after  => new_node                          ! Add new node here
      ENDIF
   ENDIF
   END SUBROUTINE macro_add_node
!*******************************************************************************
   RECURSIVE SUBROUTINE macro_find_node ( ptr, macrofile, search, ier_typ )
!
   IMPLICIT NONE
!
   TYPE(macro_internal), POINTER :: ptr     ! Pointer to current position in tree
   CHARACTER(LEN=*)             , INTENT(IN) :: macrofile
   TYPE(macro_internal), POINTER :: search  ! The structure file to be found
   INTEGER, INTENT(INOUT)        :: ier_typ
!

   IF ( LLT(macrofile, ptr%macrofile )) THEN
      IF ( ASSOCIATED(ptr%before) ) THEN
         CALL macro_find_node ( ptr%before, macrofile, search, ier_typ )
      ELSE
         ier_typ = -113
      ENDIF
   ELSEIF ( macrofile == ptr%macrofile ) THEN
      search  => ptr
      ier_typ = 0
   ELSE
      IF ( ASSOCIATED(ptr%after) ) THEN
         CALL macro_find_node ( ptr%after, macrofile, search, ier_typ )
      ELSE
         ier_typ = -113
      ENDIF
   ENDIF
   END SUBROUTINE macro_find_node
!*******************************************************************************
   RECURSIVE SUBROUTINE macro_write_node ( ptr )
!
   USE prompt_mod
   IMPLICIT NONE
!
   TYPE(macro_internal), POINTER :: ptr     ! Pointer to current position in tree
!
   IF ( ASSOCIATED(ptr) ) THEN
      IF ( ASSOCIATED(ptr%before)) THEN
         CALL macro_write_node ( ptr%before )
      ENDIF
      WRITE(output_io,1000) ptr%number, ptr%macrofile(1:70)
      IF ( ASSOCIATED(ptr%after)) THEN
         CALL macro_write_node ( ptr%after )
      ENDIF
   ELSE
      WRITE(output_io,*) 'Pointer is not associated '
      WRITE(output_io,*) ' Error in write_node'
   ENDIF
1000 FORMAT(i4,1x, a70)
   END SUBROUTINE macro_write_node
END MODULE class_macro_internal
