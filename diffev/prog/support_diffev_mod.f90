!
!
MODULE support_diffev_mod
!
PRIVATE
PUBLIC  :: make_file
!
! Contains general support routines for diffev
!
CONTAINS
!*****7**************************************************************** 
   SUBROUTINE make_file (infile, linfile, lext, n) 
!
!    builds a filename of type "base.0000"
!
   USE build_name_mod
!
   IMPLICIT none 
!                                                                       
   INTEGER, PARAMETER             :: maxw = 2
!                                                                       
   CHARACTER (LEN= * ),INTENT(INOUT) :: infile   ! Filename to be build
   INTEGER            ,INTENT(INOUT) :: linfile  ! file name length
   INTEGER            ,INTENT(IN   ) :: lext     ! extension length
   INTEGER            ,INTENT(IN   ) :: n        ! extension number
!
   CHARACTER (LEN=5)              :: string 
   CHARACTER (LEN=1024)           :: cpara (maxw) 
   INTEGER                        :: lpara (maxw) 
   INTEGER                        :: fpara 
   INTEGER                        :: ianz, ldot 
   REAL                           :: werte (maxw) 
!                                                                       
   ldot = index (infile, '.') 
   IF (ldot.eq.0) then 
      IF (lext.eq.0) then 
         cpara (1) = '"'//infile (1:linfile) //'.%d"' 
         lpara (1) = linfile+5 
      ELSE 
         WRITE (string, 1000) lext 
         cpara (1) = '"'//infile (1:linfile) //'.%'//string//'D"' 
         lpara (1) = linfile+10 
      ENDIF 
   ELSE 
      IF (lext.eq.0) then 
         cpara (1) = '"'//infile (1:ldot) //'%d"' 
         lpara (1) = ldot + 4 
      ELSE 
         WRITE (string, 1000) lext 
         cpara (1) = '"'//infile (1:ldot) //'%'//string//'D"' 
         lpara (1) = ldot + 9 
      ENDIF 
   ENDIF 
   WRITE (cpara (2), 1000) n 
   lpara (2) = 5 
   ianz = 2 
   fpara = 1 
   CALL do_build_name (ianz, cpara, lpara, werte, maxw, fpara) 
   infile = cpara (1) 
   linfile = lpara (1) 
!                                                                       
    1000 FORMAT    (i5) 
   END SUBROUTINE make_file                      
!
END MODULE support_diffev_mod
