!
! Version for stand alone DISCUS
! The file is always written to hard disk
!
SUBROUTINE powder_do_write (outfile, npkt_wrt, POW_MAXPKT, xwrt, ywrt)
!
USE errlist_mod
IMPLICIT NONE
!
CHARACTER (LEN=*)                , INTENT(IN) :: outfile
INTEGER                          , INTENT(IN) :: npkt_wrt
INTEGER                          , INTENT(IN) :: POW_MAXPKT
REAL   , DIMENSION(1:POW_MAXPKT) , INTENT(IN) :: xwrt
REAL   , DIMENSION(1:POW_MAXPKT) , INTENT(IN) :: ywrt
!
INTEGER, PARAMETER                            :: iff = 2
INTEGER   :: ii
!
CALL oeffne (iff, outfile, 'unknown') 
IF(ier_num == 0) THEN
   DO ii = 1,npkt_wrt
      WRITE( iff, *) xwrt(ii),ywrt(ii)
   ENDDO
   CLOSE(iff)
ENDIF
!
END SUBROUTINE powder_do_write
