MODULE class_internal
!
USE cryst_class
!
IMPLICIT NONE
!
PRIVATE
PUBLIC  :: internal_storage
PUBLIC  :: store_root, store_temp, read_temp
PUBLIC  :: store_add_node, store_find_node  !, store_write_node
!
TYPE :: internal_storage
!   INTEGER                           :: number      ! numer the files, not needed at the moment 
   CHARACTER (LEN=200)               :: strucfile   ! "file" name for this internal storage
   TYPE (cl_cryst)        , POINTER  :: crystal     ! The actual data structure for the crystal
   TYPE (internal_storage), POINTER  :: before      ! Pointer to prev stored crystal
   TYPE (internal_storage), POINTER  :: after       ! Pointer to next stored crystal
END TYPE internal_storage
!
TYPE(internal_storage), POINTER      :: store_root  ! The main root of the internal storage tree
TYPE(internal_storage), POINTER      :: store_temp  ! a temporary pointer to migrate the tree
TYPE(internal_storage), POINTER      ::  read_temp  ! a temporary pointer to migrate the tree
!
CONTAINS
!*******************************************************************************
   RECURSIVE SUBROUTINE store_add_node ( ptr, new_node )
!
   IMPLICIT NONE
!
   TYPE(internal_storage), POINTER :: ptr
   TYPE(internal_storage), POINTER :: new_node
   TYPE(internal_storage), POINTER :: temp    
!
   IF ( .not. ASSOCIATED(ptr))  THEN                     ! Pointer does not exist
      ptr  => new_node                                   ! Add here
   ELSEIF ( LLT(new_node%strucfile, ptr%strucfile)) THEN ! new "strucfile" is < old
      IF ( ASSOCIATED(ptr%before) ) THEN                 ! before node exists, 
         CALL store_add_node(ptr%before, new_node)       !   recursively add new node
      ELSE                                               ! before node does not exist
         ptr%before  => new_node                         ! Add new node here
      ENDIF                                              ! 
   ELSEIF ( new_node%strucfile  ==   ptr%strucfile) THEN ! New "strucfile" is = old
!
      new_node%before => ptr%before                      ! Retain old before
      new_node%after  => ptr%after                       !        and after
      CALL ptr%crystal%finalize_atoms                    ! Clean up old arrays
      DEALLOCATE(ptr%crystal)                            ! Remove old crystal
      temp => ptr                                        ! Retain address to node
      ptr  => new_node                                   ! Place new node instead of old
      DEALLOCATE(temp)                                   ! delete old node
!
   ELSEIF ( LGT(new_node%strucfile, ptr%strucfile)) THEN ! New "strucfile" is > old
      IF ( ASSOCIATED(ptr%after) ) THEN                  ! after node exists
         CALL store_add_node(ptr%after, new_node)        !    recursively add new node
      ELSE                                               ! afet node does not exist
         ptr%after  => new_node                          ! Add new node here
      ENDIF
   ENDIF
   END SUBROUTINE store_add_node
!*******************************************************************************
   RECURSIVE SUBROUTINE store_find_node ( ptr, search, ier_typ )
!
   IMPLICIT NONE
!
   TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
   TYPE(internal_storage), POINTER :: search  ! The structure file to be found
   INTEGER, INTENT(INOUT)          :: ier_typ
!
   IF ( ASSOCIATED(ptr)) THEN
      IF ( LLT(search%strucfile, ptr%strucfile )) THEN
         IF ( ASSOCIATED(ptr%before) ) THEN
            CALL store_find_node ( ptr%before, search, ier_typ )
         ELSE
            ier_typ = -113
         ENDIF
      ELSEIF ( search%strucfile == ptr%strucfile ) THEN
         search  = ptr
         ier_typ = 0
      ELSE
         IF ( ASSOCIATED(ptr%after) ) THEN
            CALL store_find_node ( ptr%after, search, ier_typ )
         ELSE
            ier_typ = -113
         ENDIF
      ENDIF
   ELSE
      ier_typ = -113
   ENDIF
   END SUBROUTINE store_find_node
!*******************************************************************************
!   RECURSIVE SUBROUTINE store_write_node ( ptr )
!!
!   USE prompt_mod
!   IMPLICIT NONE
!!
!   TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
!!
!   IF ( ASSOCIATED(ptr) ) THEN
!!      IF ( ASSOCIATED(ptr%before)) THEN
!         CALL store_write_node ( ptr%before )
!      ENDIF
!      WRITE(output_io,1000) ptr%number, ptr%strucfile
!      IF ( ASSOCIATED(ptr%after)) THEN
!         CALL store_write_node ( ptr%after )
!      ENDIF
!!   ELSE
!      WRITE(output_io,*) 'Pointer is not associated '
!      WRITE(output_io,*) ' Error in write_node'
!   ENDIF
!1000 FORMAT(i4,1x, a40)
!   END SUBROUTINE store_write_node
!*******************************************************************************
END MODULE class_internal
