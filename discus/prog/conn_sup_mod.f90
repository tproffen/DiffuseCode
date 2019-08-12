MODULE conn_sup_mod
!
USE conn_def_mod
USE conn_type_mod
!
CONTAINS
!
   SUBROUTINE get_connectivity_list (jatom, is1, ino, c_list, c_offs, natoms )
!-                                                                      
!     Get the list of neighbors for central atom jatom of type is1
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      IMPLICIT none 
!
      INTEGER, INTENT(IN)  :: jatom   ! central atom number
      INTEGER, INTENT(IN)  :: is1     ! central atom type
      INTEGER, INTENT(INOUT)  :: ino     ! Connectivity def. no.
      CHARACTER(LEN=256)   :: c_name  ! Connectivity name
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: c_list    ! Size of array c_list 
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: c_offs ! Offsets from periodic boundary
      INTEGER, INTENT(OUT) :: natoms  ! number of atoms in connectivity list
!
      INTEGER    :: i,k
!
      natoms = 0
!
!
      IF ( ALLOCATED(at_conn) ) THEN
         i = jatom 
        IF ( ASSOCIATED(at_conn(i)%liste) ) THEN      ! If neighborhood was created
          hood_temp => at_conn(i)%liste               ! point to the first neighborhood
!         hood_head => at_conn(i)%liste
          DO WHILE ( ASSOCIATED(hood_temp) )          ! While there are further neighborhood
             temp => hood_temp%nachbar                    ! point to the first neighbor within current neigborhood
             IF ( hood_temp%neigh_type == ino  .OR. &
                  hood_temp%conn_name  == c_name   ) THEN   ! This is the right neighborhood
                ino    = hood_temp%neigh_type         ! Return actual number and name
                c_name = hood_temp%conn_name
                natoms = hood_temp%natoms
                IF(ALLOCATED(c_list)) THEN
!                  IF(UBOUND(c_list).lt.natoms) THEN
                      DEALLOCATE(c_list)
                      ALLOCATE(c_list(0:natoms))
!                  ENDIF
                ELSE
                   ALLOCATE(c_list(0:natoms))
                ENDIF
                IF(ALLOCATED(c_offs)) THEN
!                  IF(UBOUND(c_offs).lt.natoms) THEN
                      DEALLOCATE(c_offs)
                      ALLOCATE(c_offs(1:3,0:natoms))
!                  ENDIF
                ELSE
                   ALLOCATE(c_offs(1:3,0:natoms))
                ENDIF
                c_list(:)   = 0   ! clear connectivity list
                c_offs(:,:) = 0
                IF(natoms == 0) RETURN                      ! Empty list
                k = 0
                DO WHILE ( ASSOCIATED(temp) )               ! While there are further neighbors
                    k         = k+ 1
                    c_list(k) = temp%atom_number
                    c_offs(1,k) = temp%offset(1)
                    c_offs(2,k) = temp%offset(2)
                    c_offs(3,k) = temp%offset(3)
                  temp => temp%next                         ! Point to next neighbor
                END DO
                RETURN                                      ! End of connectivity list
             ENDIF
             hood_temp => hood_temp%next_neighborhood     ! Point to next neighborhood
          END DO
        END IF
      ENDIF
!
   END SUBROUTINE get_connectivity_list
!
!
   SUBROUTINE get_connectivity_identity (is1, work_id, work_name, work_name_l)
!-                                                                      
!     Get the identity of a connectivity from central atom and number or name
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      IMPLICIT none 
!
!
      INTEGER,            INTENT(IN)      :: is1        ! central atom type
      INTEGER,            INTENT(INOUT)   :: work_id    ! Connectivity def. no.
      CHARACTER(LEN=256), INTENT(INOUT)   :: work_name  ! Connectivity name
      INTEGER           , INTENT(INOUT)   :: work_name_l! Connectivity name length
!
      IF ( ALLOCATED(def_main) ) THEN
         IF(is1>=0 .AND. is1<=UBOUND(def_main,1)) THEN
         is_there: IF ( ASSOCIATED(def_main(is1)%def_liste) ) THEN  ! A list of definitions exists
            def_head => def_main(is1)%def_liste
            def_temp => def_main(is1)%def_liste
            search: DO                                           ! search for working definition
               IF ( .NOT. ASSOCIATED(def_temp)) THEN             ! target is not associated ERROR
                  ier_num = -109
                  ier_typ = ER_APPL
                  RETURN
               ENDIF
               IF ( work_id   == def_temp%valid_id  .OR. &
                    work_name == def_temp%def_name       ) THEN  ! Found working definition
                  work_id     = def_temp%valid_id                ! Make sure ID matches
                  work_name   = def_temp%def_name                ! Make sure name matches
                  work_name_l = def_temp%def_name_l              ! Make sure name matches
                  EXIT search
               ENDIF
               def_head => def_temp
               def_temp => def_temp%def_next
            ENDDO search
         ELSE
            ier_num = -109
            ier_typ = ER_APPL
         ENDIF is_there
         ELSE
            ier_num = -109
            ier_typ = ER_APPL
         ENDIF
      ELSE
         ier_num = -110
         ier_typ = ER_APPL
      ENDIF
!
   END SUBROUTINE get_connectivity_identity
!
!
   INTEGER FUNCTION get_connectivity_numbers (is1)
!-                                                                      
!     Return the number of connectivities defined for an atom type
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod
      IMPLICIT none 
!
      INTEGER,            INTENT(IN)      :: is1        ! central atom type
!
      INTEGER :: numbers
!
      numbers = 0
      IF ( ALLOCATED(def_main) ) THEN
         IF(is1>=0 .AND. is1<=UBOUND(def_main,1)) THEN
         is_there: IF ( ASSOCIATED(def_main(is1)%def_liste) ) THEN  ! A list of definitions exists
            numbers = def_main(is1)%def_number
         ENDIF is_there
      ELSE
         ier_num = -109
         ier_typ = ER_APPL
      ENDIF
      ELSE
         ier_num = -110
         ier_typ = ER_APPL
      ENDIF
!
      get_connectivity_numbers = numbers   ! assign return value
!
   END FUNCTION get_connectivity_numbers
!
END MODULE conn_sup_mod
