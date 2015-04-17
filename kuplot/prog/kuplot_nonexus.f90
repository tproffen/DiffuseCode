MODULE nexus_kuplot
!*****7**************************************************************** 
!     Dummy routines in case of no NeXus support compiled in ..         
!*****7**************************************************************** 
CONTAINS
      SUBROUTINE do_nxinit 
!                                                                       
      RETURN 
      END SUBROUTINE do_nxinit                      
!*****7**************************************************************** 
      SUBROUTINE do_nxopen (line, ll) 
!                                                                       
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) line 
      INTEGER ll 
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxopen                      
!*****7**************************************************************** 
      SUBROUTINE do_nxdir 
!                                                                       
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxdir                       
!*****7**************************************************************** 
      SUBROUTINE do_nxclose 
!                                                                       
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxclose                     
!*****7**************************************************************** 
      SUBROUTINE do_nxload (line, ll) 
!                                                                       
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) line 
      INTEGER ll 
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxload                      
END MODULE nexus_kuplot
