MODULE nexus_discus
!

PRIVATE
PUBLIC  nexus_write
!
CONTAINS
!
! Dummy routines, as no NeXus is available
!
   SUBROUTINE nexus_write ( value, laver)
!
   IMPLICIT NONE
   include'errlist.inc'
   INTEGER, INTENT(IN) :: value
   LOGICAL, INTENT(IN) :: laver
!
   ier_num    = -117
   ier_typ    = ER_APPL
   ier_msg(1) = 'To write a NeXus file you need to'
   ier_msg(2) = 'install the NeXus library        '
   ier_msg(3) = 'and the NeXus version of DISCUS'
!
!
   END SUBROUTINE nexus_write
END MODULE nexus_discus
