!*****7**************************************************************** 
!     Dummy routines in case of no NeXus support compiled in ..         
!*****7**************************************************************** 
      SUBROUTINE do_nxinit 
!                                                                       
      RETURN 
      END SUBROUTINE do_nxinit                      
!*****7**************************************************************** 
      SUBROUTINE do_nxopen (line, ll) 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
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
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxdir                       
!*****7**************************************************************** 
      SUBROUTINE do_nxclose 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxclose                     
!*****7**************************************************************** 
      SUBROUTINE do_nxload (line, ll) 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      CHARACTER ( * ) line 
      INTEGER ll 
!                                                                       
      ier_num = - 53 
      ier_typ = ER_APPL 
!                                                                       
      END SUBROUTINE do_nxload                      
