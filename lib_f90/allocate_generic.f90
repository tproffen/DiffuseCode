MODULE allocate_generic
!-
!     Contains data and routines to allocate generic arrays
!     Calling interface:
!     Each alloc_* requires the following arguments:
!
!     array               !  The actual array to be allocated 
!     lower bound         !  The new lower boundary
!     upper bound         !  The new upper boundary
!     old lower bound     !  The old lower boundary
!     old upper bound     !  The old upper boundary
!     allocation status   !  will be 0 at success, 1 if new array is too small for old values
!     default value       !  Initializing value for the array
!     size of array       !  Tha allocated size in bytes
!
!     The array and the default value have to be of identical type.
!     The boundary block has to be repeated for each dimension 
!+
PRIVATE
PUBLIC  :: alloc_arr
!
integer, parameter:: dp=kind(0.d0)  ! double precision
!
  INTERFACE alloc_arr                        ! Define a generic interface to allocate n-dim arrays
    MODULE PROCEDURE alloc_1D_char, alloc_1D_int, alloc_1D_log, alloc_1D_real,  alloc_1D_cmplx,  &
                                                                alloc_1D_real8, alloc_1D_cmplx8, &
                     alloc_2D_char, alloc_2D_int, alloc_2D_log, alloc_2D_real,  alloc_2D_cmplx,  &
                                                                                alloc_2D_cmplx8, &
                     alloc_3D_char, alloc_3D_int, alloc_3D_log, alloc_3D_real,  alloc_3D_cmplx,  &
                                                                alloc_3D_real8, alloc_3D_cmplx8, &
                                    alloc_4D_int,               alloc_4D_real,                   &
                                                                alloc_5D_real
  END INTERFACE alloc_arr
!
  CONTAINS
!
!
!
   SUBROUTINE alloc_1D_char ( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, CHARACTER version
!
      IMPLICIT NONE
!
      CHARACTER (LEN=*   ), INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      CHARACTER (LEN=LEN(array))         , DIMENSION (:),ALLOCATABLE :: temp
      CHARACTER (LEN=*), INTENT(IN) :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary

      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_char
!
   SUBROUTINE alloc_1D_int ( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, INTEGER version
!
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      INTEGER,                DIMENSION (:),ALLOCATABLE :: temp
      INTEGER, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_int
!
   SUBROUTINE alloc_1D_log ( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, LOGICAL version
!
      IMPLICIT NONE
!
      LOGICAL, INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      LOGICAL,                DIMENSION (:),ALLOCATABLE :: temp
      LOGICAL, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array, stat = all_status)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_log
!
!
   SUBROUTINE alloc_1D_real ( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, REAL version
!
      IMPLICIT NONE
!
      REAL   , INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      REAL   ,                DIMENSION (:),ALLOCATABLE :: temp
      REAL   , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_real
!
!
   SUBROUTINE alloc_1D_real8( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, REAL*8 version
!
      IMPLICIT NONE
!
      REAL(dp) , INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      REAL(dp) ,                DIMENSION (:),ALLOCATABLE :: temp
      REAL(dp) , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_real8
!
   SUBROUTINE alloc_1D_cmplx ( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, COMPLEX version
!
      IMPLICIT NONE
!
      COMPLEX, INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      COMPLEX,                DIMENSION (:),ALLOCATABLE :: temp
      COMPLEX, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_cmplx
!
   SUBROUTINE alloc_1D_cmplx8 ( array, lb, ub, all_status, def_value, size_of)
!
!     Subroutine to allocate a 1-D array, COMPLEX version
!
      IMPLICIT NONE
!
      COMPLEX (KIND=KIND(0.0D0)), INTENT(INOUT), DIMENSION (:),ALLOCATABLE :: array
      COMPLEX (KIND=KIND(0.0D0)),                DIMENSION (:),ALLOCATABLE :: temp
      COMPLEX (KIND=KIND(0.0D0)), INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb
      INTEGER, INTENT(IN)    :: ub
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:1):: o_lower
      INTEGER, DIMENSION(1:1):: o_upper
      INTEGER                :: i,length
      INTEGER                :: tlb       ! temporary lower boundary
      INTEGER                :: tub       ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub < lb ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb      = MAX ( o_lower(1),lb)                    ! temporary boundaries to save old data
      tub      = MIN ( o_upper(1),ub)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb==o_lower(1) .and. ub==o_upper(1)) THEN         ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb <= tub
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb:tub), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb:tub) = array(tlb:tub)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               WRITE( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb:ub), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb:tub) = temp(tlb:tub)
               DEALlOCATE ( temp )
            ELSE
               READ ( UNIT=IT, REC=1) (array(i),i=tlb,tub)
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_1D_cmplx8
!
   SUBROUTINE alloc_2D_char ( array, lb1, ub1, lb2, ub2, all_status, def_value, size_of)
!
!     Subroutine to allocate a 2-D array, CHARACTER version
!
      IMPLICIT NONE
!
      CHARACTER (LEN=*         ), INTENT(INOUT), DIMENSION (:,:),ALLOCATABLE :: array
      CHARACTER (LEN=LEN(array)),                DIMENSION (:,:),ALLOCATABLE :: temp
      CHARACTER (LEN=*), INTENT(IN) :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2
      INTEGER, INTENT(IN)    :: ub1,  ub2
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:2):: o_lower
      INTEGER, DIMENSION(1:2):: o_upper
      INTEGER                :: i,j,length
      INTEGER                :: tlb1,tlb2 ! temporary lower boundary
      INTEGER                :: tub1,tub2 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 ) THEN                         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2) = array(tlb1:tub1, tlb2:tub2)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO j = tlb2,tub2
                  WRITE( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2) = temp(tlb1:tub1,tlb2:tub2)
               DEALlOCATE ( temp )
            ELSE
               DO j = tlb2,tub2
                  READ ( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_2D_char
!
   SUBROUTINE alloc_2D_int ( array, lb1, ub1, lb2, ub2, all_status, def_value, size_of)
!
!     Subroutine to allocate a 2-D array, INTEGER version
!
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT), DIMENSION (:,:),ALLOCATABLE :: array
      INTEGER,                DIMENSION (:,:),ALLOCATABLE :: temp
      INTEGER, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2
      INTEGER, INTENT(IN)    :: ub1,  ub2
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:2):: o_lower
      INTEGER, DIMENSION(1:2):: o_upper
      INTEGER                :: i,j,length
      INTEGER                :: tlb1,tlb2 ! temporary lower boundary
      INTEGER                :: tub1,tub2 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2) THEN         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
!
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2) = array(tlb1:tub1, tlb2:tub2)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO j = tlb2,tub2
                  WRITE( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array , stat = all_status )                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2) = temp(tlb1:tub1,tlb2:tub2)
               DEALlOCATE ( temp )
            ELSE
               DO j = tlb2,tub2
                  READ ( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_2D_int
!
   SUBROUTINE alloc_2D_log ( array, lb1, ub1, lb2, ub2, all_status, def_value, size_of)
!
!     Subroutine to allocate a 2-D array, LOGICAL version
!
      IMPLICIT NONE
!
      LOGICAL, INTENT(INOUT), DIMENSION (:,:),ALLOCATABLE :: array
      LOGICAL,                DIMENSION (:,:),ALLOCATABLE :: temp
      LOGICAL, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2
      INTEGER, INTENT(IN)    :: ub1,  ub2
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:2):: o_lower
      INTEGER, DIMENSION(1:2):: o_upper
      INTEGER                :: i,j,length
      INTEGER                :: tlb1,tlb2 ! temporary lower boundary
      INTEGER                :: tub1,tub2 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2) THEN         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2) = array(tlb1:tub1, tlb2:tub2)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO j = tlb2,tub2
                  WRITE( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2) = temp(tlb1:tub1,tlb2:tub2)
               DEALlOCATE ( temp )
            ELSE
               DO j = tlb2,tub2
                  READ ( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_2D_log
!
   SUBROUTINE alloc_2D_real ( array, lb1, ub1, lb2, ub2, all_status, def_value, size_of)
!
!     Subroutine to allocate a 2-D array, REAL version
!
      IMPLICIT NONE
!
      REAL   , INTENT(INOUT), DIMENSION (:,:),ALLOCATABLE :: array
      REAL   ,                DIMENSION (:,:),ALLOCATABLE :: temp
      REAL   , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2
      INTEGER, INTENT(IN)    :: ub1,  ub2
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:2):: o_lower
      INTEGER, DIMENSION(1:2):: o_upper
      INTEGER                :: i,j,length
      INTEGER                :: tlb1,tlb2 ! temporary lower boundary
      INTEGER                :: tub1,tub2 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2) THEN         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2) = array(tlb1:tub1, tlb2:tub2)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO j = tlb2,tub2
                  WRITE( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2) = temp(tlb1:tub1,tlb2:tub2)
               DEALlOCATE ( temp )
            ELSE
               DO j = tlb2,tub2
                  READ ( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_2D_real
!
   SUBROUTINE alloc_2D_cmplx ( array, lb1, ub1, lb2, ub2, all_status, def_value, size_of)
!
!     Subroutine to allocate a 2-D array, COMPLEX version
!
      IMPLICIT NONE
!
      COMPLEX, INTENT(INOUT), DIMENSION (:,:),ALLOCATABLE :: array
      COMPLEX,                DIMENSION (:,:),ALLOCATABLE :: temp
      COMPLEX, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2
      INTEGER, INTENT(IN)    :: ub1,  ub2
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:2):: o_lower
      INTEGER, DIMENSION(1:2):: o_upper
      INTEGER                :: i,j,length
      INTEGER                :: tlb1,tlb2 ! temporary lower boundary
      INTEGER                :: tub1,tub2 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2) THEN         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2) = array(tlb1:tub1, tlb2:tub2)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO j = tlb2,tub2
                  WRITE( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2) = temp(tlb1:tub1,tlb2:tub2)
               DEALlOCATE ( temp )
            ELSE
               DO j = tlb2,tub2
                  READ ( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_2D_cmplx
!
   SUBROUTINE alloc_2D_cmplx8 ( array, lb1, ub1, lb2, ub2, all_status, def_value, size_of)
!
!     Subroutine to allocate a 2-D array, COMPLEX version
!
      IMPLICIT NONE
!
      COMPLEX (KIND=KIND(0.0D0)), INTENT(INOUT), DIMENSION (:,:),ALLOCATABLE :: array
      COMPLEX (KIND=KIND(0.0D0)),                DIMENSION (:,:),ALLOCATABLE :: temp
      COMPLEX (KIND=KIND(0.0D0)), INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2
      INTEGER, INTENT(IN)    :: ub1,  ub2
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:2):: o_lower
      INTEGER, DIMENSION(1:2):: o_upper
      INTEGER                :: i,j,length
      INTEGER                :: tlb1,tlb2 ! temporary lower boundary
      INTEGER                :: tub1,tub2 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2) THEN         ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2) = array(tlb1:tub1, tlb2:tub2)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO j = tlb2,tub2
                  WRITE( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2) = temp(tlb1:tub1,tlb2:tub2)
               DEALlOCATE ( temp )
            ELSE
               DO j = tlb2,tub2
                  READ ( UNIT=IT, REC=j) (array(i,j),i=tlb1,tub1)
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_2D_cmplx8
!
!
   SUBROUTINE alloc_3D_char ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, CHARACTER version
!
      IMPLICIT NONE
!
      CHARACTER (LEN=*         ), INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      CHARACTER (LEN=LEN(array)),                DIMENSION (:,:,:),ALLOCATABLE :: temp
      CHARACTER (LEN=*), INTENT(IN) :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_char
!
!
   SUBROUTINE alloc_3D_int  ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, INTEGER version
!
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      INTEGER,                DIMENSION (:,:,:),ALLOCATABLE :: temp
      INTEGER, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_int
!
!
   SUBROUTINE alloc_3D_log  ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, LOGICAL version
!
      IMPLICIT NONE
!
      LOGICAL, INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      LOGICAL,                DIMENSION (:,:,:),ALLOCATABLE :: temp
      LOGICAL, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_log
!
!
   SUBROUTINE alloc_3D_real ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, REAL version
!
      IMPLICIT NONE
!
      REAL   , INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      REAL   ,                DIMENSION (:,:,:),ALLOCATABLE :: temp
      REAL   , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_real
!
!
   SUBROUTINE alloc_3D_cmplx ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, COMPLEX version
!
      IMPLICIT NONE
!
      COMPLEX, INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      COMPLEX,                DIMENSION (:,:,:),ALLOCATABLE :: temp
      COMPLEX, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_cmplx
!
!
   SUBROUTINE alloc_3D_real8 ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, REAL version
!
      IMPLICIT NONE
!
      REAL(dp) , INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      REAL(dp) ,                DIMENSION (:,:,:),ALLOCATABLE :: temp
      REAL(dp) , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_real8
!
!
   SUBROUTINE alloc_3D_cmplx8( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, all_status, def_value, size_of)
!
!     Subroutine to allocate a 3-D array, COMPLEX version
!
      IMPLICIT NONE
!
      COMPLEX (KIND=KIND(0.0D0)), INTENT(INOUT), DIMENSION (:,:,:),ALLOCATABLE :: array
      COMPLEX (KIND=KIND(0.0D0)),                DIMENSION (:,:,:),ALLOCATABLE :: temp
      COMPLEX (KIND=KIND(0.0D0)), INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:3):: o_lower
      INTEGER, DIMENSION(1:3):: o_upper
      INTEGER                :: i,j,k,length
      INTEGER                :: tlb1,tlb2,tlb3 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp(tlb1:tub1, tlb2:tub2, tlb3:tub3) = array(tlb1:tub1, tlb2:tub2, tlb3:tub3)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3) = temp(tlb1:tub1,tlb2:tub2,tlb3:tub3)
               DEALlOCATE ( temp )
            ELSE
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k),i=tlb1,tub1)
                  ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_3D_cmplx8
!
!
   SUBROUTINE alloc_4D_int  ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, lb4, ub4, &
                              all_status, def_value, size_of)
!
!     Subroutine to allocate a 4-D array, REAL version
!
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT), DIMENSION (:,:,:,:),ALLOCATABLE :: array
      INTEGER,                DIMENSION (:,:,:,:),ALLOCATABLE :: temp
      INTEGER, INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3,  lb4
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3,  ub4
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:4):: o_lower
      INTEGER, DIMENSION(1:4):: o_upper
      INTEGER                :: i,j,k,l,length
      INTEGER                :: tlb1,tlb2,tlb3,tlb4 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3,tub4 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!

      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 .or. ub4 < lb4 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      tlb4     = MAX ( o_lower(4),lb4)
      tub4     = MIN ( o_upper(4),ub4)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)   .and. &
              lb4==o_lower(4) .and. ub4==o_upper(4)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3  .and. &
                     o_upper(4) >= o_lower(4) .and. tlb4 <= tub4
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp (tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4) = &
               array(tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO l = tlb4,tub4
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k,l),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3,tlb4:tub4) = &
               temp (tlb1:tub1,tlb2:tub2,tlb3:tub3,tlb4:tub4)
               DEALlOCATE ( temp )
            ELSE
               DO l = tlb4,tub4
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k,l),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_4D_int
!
!
   SUBROUTINE alloc_4D_real ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, lb4, ub4, &
                              all_status, def_value, size_of)
!
!     Subroutine to allocate a 4-D array, REAL version
!
      IMPLICIT NONE
!
      REAL   , INTENT(INOUT), DIMENSION (:,:,:,:),ALLOCATABLE :: array
      REAL   ,                DIMENSION (:,:,:,:),ALLOCATABLE :: temp
      REAL   , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3,  lb4
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3,  ub4
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:4):: o_lower
      INTEGER, DIMENSION(1:4):: o_upper
      INTEGER                :: i,j,k,l,length
      INTEGER                :: tlb1,tlb2,tlb3,tlb4 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3,tub4 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!

      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 .or. ub4 < lb4 ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      tlb4     = MAX ( o_lower(4),lb4)
      tub4     = MIN ( o_upper(4),ub4)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)   .and. &
              lb4==o_lower(4) .and. ub4==o_upper(4)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3  .and. &
                     o_upper(4) >= o_lower(4) .and. tlb4 <= tub4
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp (tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4) = &
               array(tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO l = tlb4,tub4
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k,l),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3, lb4:ub4), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3,tlb4:tub4) = &
               temp (tlb1:tub1,tlb2:tub2,tlb3:tub3,tlb4:tub4)
               DEALlOCATE ( temp )
            ELSE
               DO l = tlb4,tub4
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k,l),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_4D_real
!
!
   SUBROUTINE alloc_5D_real ( array, lb1, ub1, lb2, ub2, &
                                     lb3, ub3, lb4, ub4, &
                                     lb5, ub5, &
                              all_status, def_value, size_of)
!
!     Subroutine to allocate a 5-D array, REAL version
!
      IMPLICIT NONE
!
      REAL   , INTENT(INOUT), DIMENSION (:,:,:,:,:),ALLOCATABLE :: array
      REAL   ,                DIMENSION (:,:,:,:,:),ALLOCATABLE :: temp
      REAL   , INTENT(IN)    :: def_value
      INTEGER, INTENT(IN)    :: lb1,  lb2,  lb3,  lb4,  lb5
      INTEGER, INTENT(IN)    :: ub1,  ub2,  ub3,  ub4,  ub5
      INTEGER, INTENT(INOUT) :: all_status
      INTEGER, INTENT(OUT)   :: size_of
!
      INTEGER, PARAMETER     :: IT = 88
      INTEGER, DIMENSION(1:5):: o_lower
      INTEGER, DIMENSION(1:5):: o_upper
      INTEGER                :: i,j,k,l,m,length
      INTEGER                :: tlb1,tlb2,tlb3,tlb4,tlb5 ! temporary lower boundary
      INTEGER                :: tub1,tub2,tub3,tub4,tub5 ! temporary upper boundary
      LOGICAL                :: ltemp
      LOGICAL                :: lrestore
!
      all_status = 0
      size_of    = 0
      o_lower    = 0
      o_upper    = 0
!
      IF ( ub1 < lb1 .or. ub2 < lb2 .or. ub3 < lb3 .or.  &
           ub4 < lb4 .or. ub5 < lb5                    ) THEN  ! boundaries are wrong; return
        all_status = -1
        RETURN
      ENDIF
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         o_lower = LBOUND(array)
         o_upper = UBOUND(array)
      ENDIF
!
      tlb1     = MAX ( o_lower(1),lb1)                    ! temporary boundaries to save old data
      tub1     = MIN ( o_upper(1),ub1)
      tlb2     = MAX ( o_lower(2),lb2)
      tub2     = MIN ( o_upper(2),ub2)
      tlb3     = MAX ( o_lower(3),lb3)
      tub3     = MIN ( o_upper(3),ub3)
      tlb4     = MAX ( o_lower(4),lb4)
      tub4     = MIN ( o_upper(4),ub4)
      tlb5     = MAX ( o_lower(5),lb4)
      tub5     = MIN ( o_upper(5),ub4)
      ltemp    = .true.
!
      lrestore = .false.
!
      IF ( allocated(array) ) THEN                ! The array is allocated
         IF ( lb1==o_lower(1) .and. ub1==o_upper(1)   .and. &
              lb2==o_lower(2) .and. ub2==o_upper(2)   .and. &
              lb3==o_lower(3) .and. ub3==o_upper(3)   .and. &
              lb4==o_lower(4) .and. ub4==o_upper(4)   .and. &
              lb5==o_lower(5) .and. ub5==o_upper(5)         ) THEN ! Boundaries are same
            RETURN
         ENDIF
         lrestore =  o_upper(1) >= o_lower(1) .and. tlb1 <= tub1  .and. &
                     o_upper(2) >= o_lower(2) .and. tlb2 <= tub2  .and. &
                     o_upper(3) >= o_lower(3) .and. tlb3 <= tub3  .and. &
                     o_upper(4) >= o_lower(4) .and. tlb4 <= tub4  .and. &
                     o_upper(5) >= o_lower(5) .and. tlb5 <= tub5
         IF ( lrestore ) THEN                     ! There are old data to be saved
            ALLOCATE ( temp(tlb1:tub1, tlb2:tub2, tlb3:tub3,  &
                            tlb4:tub4, tlb5:tub5), stat = all_status )
            IF ( all_status == 0 ) THEN           ! Success, use temporary array
               temp (tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4, tlb5:tub5) = &
               array(tlb1:tub1, tlb2:tub2, tlb3:tub3, tlb4:tub4, tlb5:tub5)
               ltemp = .true.
            ELSE                                  ! Could not allocate temp, try to write to disk
               INQUIRE ( iolength=length) array(:,1,1,1,1)
               OPEN( UNIT=IT, STATUS='scratch', ACCESS='direct', RECL=length, FORM='unformatted')
               DO m = tlb5,tub5
               DO l = tlb4,tub4
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     WRITE( UNIT=IT, REC=j) (array(i,j,k,l,m),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
               ltemp = .false.
            END IF
         END IF
         DEALLOCATE ( array)                      ! Deallocate the old array
      END IF
      ALLOCATE ( array(lb1:ub1, lb2:ub2, lb3:ub3,  &
                       lb4:ub4, lb5:ub5), stat = all_status) ! Allocate with new boundaries
      IF ( all_status == 0 ) THEN                 ! Success
         array = def_value
         IF ( lrestore ) THEN                     ! There are old data to be saved
            IF ( ltemp ) THEN
               array(tlb1:tub1,tlb2:tub2,tlb3:tub3,tlb4:tub4,tlb5:tub5) = &
               temp (tlb1:tub1,tlb2:tub2,tlb3:tub3,tlb4:tub4,tlb5:tub5)
               DEALlOCATE ( temp )
            ELSE
               DO m = tlb5,tub5
               DO l = tlb4,tub4
               DO k = tlb3,tub3
                  DO j = tlb2,tub2
                     READ ( UNIT=IT, REC=j) (array(i,j,k,l,m),i=tlb1,tub1)
                  ENDDO
               ENDDO
               ENDDO
               ENDDO
               CLOSE( UNIT=IT)
            ENDIF
         ENDIF
      ENDIF
!
      size_of = 1 ! SIZEOF(array)
!
      RETURN
!
   END SUBROUTINE alloc_5D_real
!
!
END MODULE allocate_generic
