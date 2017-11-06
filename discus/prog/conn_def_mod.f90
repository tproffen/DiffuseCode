MODULE conn_def_mod
!
IMPLICIT NONE
!
! Connectivity settings:
!
INTEGER             :: conn_nmax  =  1
!
INTEGER, PARAMETER  :: code_res   = -2                ! Reset=remove the whole list
INTEGER, PARAMETER  :: code_del   = -1                ! Remove one existing definition
INTEGER, PARAMETER  :: code_add   =  1                ! Add a new definition
INTEGER, PARAMETER  :: code_set   =  2                ! Set one existing definition
!
INTEGER, PARAMETER  :: MOLE_SCOPE_IGN  =  0
INTEGER, PARAMETER  :: MOLE_SCOPE_WITH =  1
INTEGER, PARAMETER  :: MOLE_SCOPE_OUT  = -1
!
!!!!!!!!!
TYPE :: CONN_DEFS
   INTEGER                            :: valid_id     ! ID = no of current definition
   CHARACTER (LEN=256)                :: def_name     ! Name of current definition
   INTEGER                            :: def_name_l   ! Length of connectivity name
   INTEGER                            :: mmc_sel      ! This definition may be used by mmc to select
   INTEGER                            :: mmc_ene      ! This definition may be used by mmc energy
   INTEGER                            :: valid_no     ! Number of valid atom types
   INTEGER                            :: intend_no    ! Intended number of neighbors
   INTEGER, DIMENSION(:), ALLOCATABLE :: valid_types  ! Valid atom types in current definition
   INTEGER                            :: mole_scope   ! Scope is limited to a molecule
   REAL                               :: def_rmin     ! minimum distance for current definition
   REAL                               :: def_rmax     ! maximum distance for current definition
   TYPE(CONN_DEFS), POINTER           :: def_next     ! next definition
END TYPE
!
TYPE :: CONN_MAIN
   INTEGER                            :: def_id       ! ID = no of Scattering type
   INTEGER                            :: def_number   ! no of valid definitions for this type
   TYPE(CONN_DEFS), POINTER           :: def_liste    ! The list of settings
END TYPE
!
TYPE(CONN_MAIN), DIMENSION(:), ALLOCATABLE :: def_main  ! Array of size MAXSCAT
!
TYPE(CONN_DEFS), POINTER              :: def_head, def_tail, def_temp
TYPE(CONN_MAIN), POINTER              :: main_head, main_tail, main_temp
!
!
END MODULE conn_def_mod
