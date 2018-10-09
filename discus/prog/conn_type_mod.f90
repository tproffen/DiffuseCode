MODULE conn_type_mod
!
! TYPE NEIGHBORS  is a linear chain to be filled with the actual neighboring atoms
!                 currently only the neighbor number is stored no further info, as
!                 I do not want to update this list all the time.
TYPE :: NEIGHBORS
   INTEGER                      :: atom_number
   INTEGER, DIMENSION(1:3)      :: offset
   TYPE (NEIGHBORS), POINTER    :: next
END TYPE
!
! TYPE NEIGHBORHOOD is a linear chain of the possible neighborhoods.
!                   Each node branches off to one side with the actual neighboring atoms
TYPE :: NEIGHBORHOOD
   INTEGER                      :: central_number     ! absolute number of central atom
   INTEGER                      :: central_type       ! central atom is of this type
   INTEGER                      :: neigh_type         ! this neighbors belongs to this definition
   CHARACTER (LEN=256)          :: conn_name          ! Connectivity name
   INTEGER                      :: conn_name_l        ! Connectivity name length
   INTEGER                      :: mmc_sel            ! This connectivity may be used by mmc to select
   INTEGER                      :: mmc_ene            ! This connectivity may be used by mmc energy
   REAL                         :: distance_min       ! minimum distance to neighbors
   REAL                         :: distance_max       ! maximum distance to neighbors
   INTEGER                      :: natoms             ! number of neighbors
   TYPE (NEIGHBORS), POINTER    :: nachbar            ! The actual list of neighboring atoms
   TYPE (NEIGHBORHOOD), POINTER :: next_neighborhood  ! A next neighborhood
END TYPE
!
! TYPE MAIN_LIST is a structure that contains info on the central atom, as
!                well as a pointer to the NEIGHBORHOOD
TYPE :: MAIN_LIST
   INTEGER                      :: number
   TYPE (NEIGHBORHOOD), POINTER :: liste
END TYPE
!
! In order to have FAST access to the neighborhood of any atom, an
! allocatable array is defined. Entry at_conn(i) gives access to the 
! neighborhood of atom i
TYPE (main_list), DIMENSION(:), ALLOCATABLE :: at_conn
!
! (temporary) pointers of TYPE NEIGHBORS. This allows to move along the 
! neighbors in an individual neighborhood.
TYPE (NEIGHBORS), POINTER       :: head, tail, temp
!
! (temporary) pointers of TYPE NEIGHBORHOOD. This allows to move along the 
! neighborhoods of an individual atom.
TYPE (NEIGHBORHOOD), POINTER       :: hood_head
TYPE (NEIGHBORHOOD), POINTER       :: hood_temp
TYPE (NEIGHBORHOOD), POINTER       :: hood_prev
TYPE (NEIGHBORHOOD), POINTER       :: hood_central
TYPE (NEIGHBORHOOD), POINTER       :: hood_second
!
INTEGER, PARAMETER                 :: STATUS_ON  =  1
INTEGER, PARAMETER                 :: STATUS_OFF = -1
INTEGER, PARAMETER                 :: STATUS_IGN =  0
LOGICAL                            :: conn_status = .FALSE.
INTEGER                            :: conn_mmc_sel = STATUS_IGN ! For later use 
INTEGER                            :: conn_mmc_ene = STATUS_IGN ! For later use 
!
END MODULE conn_type_mod
