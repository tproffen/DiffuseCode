MODULE class_internal
!
USE cryst_class
!
IMPLICIT NONE
!
PRIVATE
PUBLIC  :: internal_storage
PUBLIC  :: store_root, store_temp, read_temp, read_from, read_parent
PUBLIC  :: store_add_node, store_find_node, store_list_node, store_remove_all
PUBLIC  :: store_remove_single
!
TYPE :: internal_storage
!   INTEGER                           :: number      ! numer the files, not needed at the moment 
   CHARACTER (LEN=200)               :: strucfile   ! "file" name for this internal storage
   TYPE (cl_cryst)        , POINTER  :: crystal     ! The actual data structure for the crystal
   TYPE (internal_storage), POINTER  :: before      ! Pointer to prev stored crystal
   TYPE (internal_storage), POINTER  :: after       ! Pointer to next stored crystal
END TYPE internal_storage
!
TYPE(internal_storage), POINTER      :: store_root=>NULL()  ! The main root of the internal storage tree
TYPE(internal_storage), POINTER      :: store_temp=>NULL()  ! a temporary pointer to migrate the tree
TYPE(internal_storage), POINTER      ::  read_temp=>NULL()  ! a temporary pointer to migrate the tree
TYPE(internal_storage), POINTER      ::  read_from=>NULL()  ! a temporary pointer to migrate the tree
TYPE(internal_storage), POINTER      ::read_parent=>NULL()  ! a temporary pointer to migrate the tree
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
      NULLIFY(ptr%before)
      NULLIFY(ptr%after)
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
      CALL finalize_atoms(ptr%crystal)                   ! Clean up old arrays
      DEALLOCATE(ptr%crystal)                            ! Remove old crystal
      temp => ptr                                        ! Retain address to node
      ptr  => new_node                                   ! Place new node instead of old
!     DEALLOCATE(temp)                                   ! delete old node
!
   ELSEIF ( LGT(new_node%strucfile, ptr%strucfile)) THEN ! New "strucfile" is > old
      IF ( ASSOCIATED(ptr%after) ) THEN                  ! after node exists
         CALL store_add_node(ptr%after, new_node)        !    recursively add new node
      ELSE                                               ! afet node does not exist
         ptr%after  => new_node                          ! Add new node here
      ENDIF
   ENDIF
   END SUBROUTINE store_add_node
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE store_find_node ( ptr, from, search, parent, ier_num )
!
   IMPLICIT NONE
!
   TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
   TYPE(internal_storage), POINTER :: from    ! Pointer to current parent
   TYPE(internal_storage), POINTER :: search  ! The structure file to be found
   TYPE(internal_storage), POINTER :: parent  ! The structure's parent
   INTEGER, INTENT(INOUT)          :: ier_num
!
   IF ( ASSOCIATED(ptr)) THEN
      IF ( LLT(search%strucfile, ptr%strucfile )) THEN
         IF ( ASSOCIATED(ptr%before) ) THEN
            CALL store_find_node ( ptr%before, ptr, search, parent, ier_typ )
         ELSE
            ier_num = -113
         ENDIF
      ELSEIF ( search%strucfile == ptr%strucfile ) THEN
         search  = ptr
         parent => from
         ier_typ = 0
      ELSE
         IF ( ASSOCIATED(ptr%after) ) THEN
            CALL store_find_node ( ptr%after, ptr, search, parent, ier_typ )
         ELSE
            ier_num = -113
         ENDIF
      ENDIF
   ELSE
      ier_num = -113
   ENDIF
   END SUBROUTINE store_find_node
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE store_find_before ( ptr, search, ier_num )
!
   IMPLICIT NONE
!
   TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
   TYPE(internal_storage), POINTER :: search  ! The structure file to be found
   INTEGER, INTENT(INOUT)          :: ier_num
!
   IF ( ASSOCIATED(ptr)) THEN
      IF ( ASSOCIATED(ptr%before) ) THEN
         CALL store_find_before ( ptr%before, search, ier_num )
      ELSE
         search  => ptr
         ier_num = 0
      ENDIF
   ELSE
      ier_num = -113
   ENDIF
   END SUBROUTINE store_find_before
!
!*******************************************************************************
!
   RECURSIVE SUBROUTINE store_find_after ( ptr, search, ier_num )
!
   IMPLICIT NONE
!
   TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
   TYPE(internal_storage), POINTER :: search  ! The structure file to be found
   INTEGER, INTENT(INOUT)          :: ier_num
!
   IF ( ASSOCIATED(ptr)) THEN
      IF ( ASSOCIATED(ptr%after) ) THEN
         CALL store_find_after ( ptr%after, search, ier_num )
      ELSE
         search  =>  ptr
         ier_num = 0
      ENDIF
   ELSE
      ier_num = -113
   ENDIF
   END SUBROUTINE store_find_after
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
SUBROUTINE store_list_node( flag_all, node_name, display )
!
USE prompt_mod
USE param_mod
USE errlist_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: flag_all
CHARACTER(LEN=*), INTENT(IN) :: node_name
CHARACTER(LEN=*), INTENT(IN) :: display
!
TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
TYPE(internal_storage), POINTER :: from    ! Pointer to current parent
TYPE(internal_storage), POINTER :: search  ! The structure file to be found
TYPE(internal_storage), POINTER :: parent  ! The structure's parent
!
IF(flag_all) THEN
   IF(.NOT.ASSOCIATED(store_root)) THEN
      WRITE(output_io, '(a)') 'No internal files exist'
      res_para(0:1) = 0
   ELSE
      ptr => store_root
      NULLIFY(from)
      CALL store_list_all( ptr, from, display)
      res_para(0:1) = 0
   ENDIF
ELSE
  ALLOCATE(search)
   search%strucfile = node_name
   NULLIFY(from)
   NULLIFY(parent)
   CALL store_find_node(store_root, from, search, parent, ier_num ) ! Find the proper node
   IF(ier_num == 0) THEN
      CALL store_show_single(search, parent, display)
   ELSE
      ier_typ = ER_APPL
   ENDIF
  DEALLOCATE(search)
ENDIF
!
END SUBROUTINE store_list_node
!
!*******************************************************************************
!
RECURSIVE SUBROUTINE store_list_all( ptr, daddy, display)
!
TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
TYPE(internal_storage), POINTER :: daddy   ! Pointer to current position in tree
CHARACTER(LEN=*), INTENT(IN) :: display
!
IF(ASSOCIATED(ptr) ) THEN
   IF(ASSOCIATED(ptr%before)) THEN
      CALL store_list_all( ptr%before, ptr, display)
   ENDIF
   CALL store_show_single(ptr, daddy, display)
   IF(ASSOCIATED(ptr%after)) THEN
      CALL store_list_all( ptr%after, ptr, display)
   ENDIF
ENDIF
END SUBROUTINE store_list_all
!
!*******************************************************************************
!
SUBROUTINE store_show_single( ptr, daddy, display)
!
USE prompt_mod
USE param_mod
!
TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
TYPE(internal_storage), POINTER :: daddy   ! Pointer to current position in tree
CHARACTER(LEN=*), INTENT(IN) :: display
!
INTEGER :: natoms
INTEGER :: nscat 
INTEGER :: n_mole
INTEGER :: n_type
INTEGER :: n_atom
!
IF(ptr%strucfile /= ' ') THEN
   WRITE(output_io, '(a,a)') ' Internal file : ',ptr%strucfile(1:LEN_TRIM(ptr%strucfile))
   res_para(0:1) = 0
   IF(display=='full' .OR. display=='crystal') THEN
      natoms = ptr%crystal%get_natoms()
      nscat  = ptr%crystal%get_nscat ()
      n_mole = ptr%crystal%get_n_mole()
      n_type = ptr%crystal%get_n_type()
      WRITE(output_io, '(a,a)') '     Content   : ', ' Atoms   Types   Molecules MolTypes'
      WRITE(output_io, '(a,4I7)')   '               : ', natoms, nscat, n_mole, n_type
      res_para(1) = FLOAT(natoms)
      res_para(2) = FLOAT(nscat)
      res_para(3) = FLOAT(n_mole)
      res_para(4) = FLOAT(n_type)
      res_para(0) = 4.0
   ENDIF
   IF(display=='full' .OR. display=='connect') THEN
   IF(associated(daddy)) THEN
      WRITE(output_io,'(a,a)')  '     Parent is : ', daddy%strucfile(1:LEN_TRIM(daddy%strucfile))
   ENDIF
   IF(ASSOCIATED(ptr%before)) THEN
      WRITE(output_io,'(a,a)') '    before kid : ', ptr%before%strucfile(1:LEN_TRIM(ptr%before%strucfile))
   ENDIF
   IF(ASSOCIATED(ptr%after)) THEN
      WRITE(output_io,'(a,a)') '    after  kid : ', ptr%after%strucfile(1:LEN_TRIM(ptr%after%strucfile))
   ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE store_show_single
!
!*******************************************************************************
!
RECURSIVE SUBROUTINE store_remove_all( ptr)
!
TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
!
IF(ASSOCIATED(ptr) ) THEN
   IF(ASSOCIATED(ptr%before)) THEN
      CALL store_remove_all( ptr%before)
   ENDIF
   IF(ASSOCIATED(ptr%after)) THEN
      CALL store_remove_all( ptr%after)
   ENDIF
!
   CALL finalize_atoms(ptr%crystal)
   DEALLOCATE(ptr%crystal)
   ptr%strucfile = ' '
   NULLIFY(ptr) !DEALLOCATE(ptr)
ENDIF
END SUBROUTINE store_remove_all
!
!*******************************************************************************
!
SUBROUTINE store_remove_single( strucfile, ier_num)
!
CHARACTER (LEN=*), INTENT(IN)   :: strucfile   ! "file" name for this internal storage
TYPE(internal_storage), POINTER :: ptr     ! Pointer to current position in tree
TYPE(internal_storage), POINTER :: from    ! Pointer to current parent
!TYPE(internal_storage), POINTER :: temp    ! Pointer to current position in tree
TYPE(internal_storage), POINTER :: search  ! Pointer to current position in tree
TYPE(internal_storage), POINTER :: parent  ! The structure's parent
INTEGER,         INTENT(OUT)    :: ier_num ! Pointer to current position in tree
!
ALLOCATE(search)
search%strucfile = strucfile
NULLIFY(from)
NULLIFY(parent)
CALL store_find_node ( store_root, from, search, parent, ier_num )
IF(ier_num == 0) THEN
!
!  No children at all
!
   IF(.NOT. ASSOCIATED(search%before) .AND. .NOT.ASSOCIATED(search%after)) THEN
      CALL finalize_atoms(search%crystal)
      DEALLOCATE(search%crystal)
      IF(ASSOCIATED(parent)) THEN     ! We are not at root
!
         IF(associated(parent%before)) then
            IF(parent%before%strucfile == search%strucfile) THEN
               NULLIFY(parent%before)
            ENDIF
         ENDIF
         IF(associated(parent%after)) then
            IF(parent%after%strucfile == search%strucfile) THEN
               NULLIFY(parent%after)
            ENDIF
         ENDIF
      ELSE
         NULLIFY(store_root)    ! Node had no parents, nullify root
      ENDIF
      DEALLOCATE(search)
!
!  Exactly one BEFORE child
!
   ELSEIF(      ASSOCIATED(search%before) .AND. .NOT.ASSOCIATED(search%after)) THEN
      CALL finalize_atoms(search%crystal)
      DEALLOCATE(search%crystal)
      IF(ASSOCIATED(parent)) THEN     ! We are not at root
         IF(associated(parent%before)) then
            IF(parent%before%strucfile == search%strucfile) THEN
               parent%before => search%before
            ENDIF
         ENDIF
         IF(associated(parent%after)) then
            IF(parent%after%strucfile == search%strucfile) THEN
               parent%after => search%before
            ENDIF
         ENDIF
      ELSE
         store_root=> search%before
      ENDIF
!
!  Exactly one AFTER  child
!
   ELSEIF( .NOT.ASSOCIATED(search%before) .AND.      ASSOCIATED(search%after)) THEN
      CALL finalize_atoms(search%crystal)
      DEALLOCATE(search%crystal)
      IF(ASSOCIATED(parent)) THEN     ! We are not at root
         IF(associated(parent%before)) then
            IF(parent%before%strucfile == search%strucfile) THEN
               parent%before => search%after
            ENDIF
         ENDIF
         IF(associated(parent%after)) then
            IF(parent%after%strucfile == search%strucfile) THEN
               parent%after => search%after
            ENDIF
         ENDIF
      ELSE
         store_root=> search%before
      ENDIF
!
!  Exactly both children BEFORE and AFTER
!
   ELSE
      CALL finalize_atoms(search%crystal)
      DEALLOCATE(search%crystal)
      CALL store_find_after ( search%before, ptr, ier_typ ) ! Find last node in left branch
      ptr%after => search%after       ! connect rightmost node in left branch to right child
      IF(ASSOCIATED(parent)) THEN     ! We are not at root, connect parent to proper child
         IF(associated(parent%before)) then
            IF(parent%before%strucfile == search%strucfile) THEN
               parent%before => search%before
            ENDIF
         ENDIF
         IF(associated(parent%after)) then
            IF(parent%after%strucfile == search%strucfile) THEN
               parent%after => search%before
            ENDIF
         ENDIF
      ELSE                            ! root was removed, point to left child
         store_root=> search%before
      ENDIF
   ENDIF
ENDIF
NULLIFY(search)
!
END SUBROUTINE store_remove_single
!
!*******************************************************************************
!
END MODULE class_internal
