!
! Version for stand alone DISCUS
! The file is always written to hard disk
!
SUBROUTINE powder_do_write (outfile, npkt_wrt, xwrt, ywrt)
!
USE errlist_mod
use precision_mod
IMPLICIT NONE
!
CHARACTER (LEN=*)                , INTENT(IN) :: outfile
INTEGER                          , INTENT(IN) :: npkt_wrt
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_wrt  ) , INTENT(IN) :: xwrt
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_wrt  ) , INTENT(IN) :: ywrt
!
INTEGER, PARAMETER                            :: iff = 2
INTEGER   :: ii
!
CALL oeffne (iff, outfile, 'unknown') 
IF(ier_num == 0) THEN
   DO ii = 0,npkt_wrt
      WRITE( iff, *) xwrt(ii),ywrt(ii)
   ENDDO
   CLOSE(iff)
ENDIF
!
END SUBROUTINE powder_do_write
