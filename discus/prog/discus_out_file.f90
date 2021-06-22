!
!  Version for stand alone DISCUS
!  The file is always written to hard disk
!  Separate routines for 1D and 2D files
!
SUBROUTINE output_save_file_1d( outfile, npkt1, xwrt, ywrt, mode )
!
USE errlist_mod
IMPLICIT NONE
!
CHARACTER (LEN=*),           INTENT(IN) :: outfile
INTEGER,                     INTENT(IN) :: npkt1
REAL   , DIMENSION(1:npkt1), INTENT(IN) :: xwrt
REAL   , DIMENSION(1:npkt1), INTENT(IN) :: ywrt
integer                    , intent(in) :: mode
!
INTEGER, PARAMETER :: IFF = 2
INTEGER            :: i
!
CALL oeffne(IFF, outfile, 'unknown')
IF(ier_num == 0) THEN
   DO i=1, npkt1
      WRITE(iff, 1000) xwrt(i), ywrt(i)
   ENDDO
ENDIF
CLOSE(IFF)
1000 FORMAT(2(1x,E12.5))
!
END SUBROUTINE output_save_file_1d
!
!  2D output file
!
SUBROUTINE output_save_file_2d( outfile, ranges, npkt1, npkt2, zwrt, &
           header_lines, nheader)
!
USE errlist_mod
IMPLICIT NONE
!
CHARACTER (LEN=*),           INTENT(IN) :: outfile
REAL   , DIMENSION(1:4),     INTENT(IN) :: ranges
INTEGER,                     INTENT(IN) :: npkt1
INTEGER,                     INTENT(IN) :: npkt2
REAL   , DIMENSION(1:npkt1, 1:npkt2), INTENT(IN) :: zwrt
CHARACTER (LEN=160), DIMENSION(:), INTENT(IN) :: header_lines
INTEGER,                           INTENT(IN) :: nheader! number of lines in header

!
INTEGER, PARAMETER :: IFF = 2
INTEGER            :: i, j
!
CALL oeffne(IFF, outfile, 'unknown')
IF(ier_num == 0) THEN
   WRITE(iff, 1000) npkt1, npkt2
   WRITE(iff, 1100) ranges
   DO j=1, npkt2
      WRITE(iff, 2000) (zwrt(i,j), i=1,npkt1)
   ENDDO
   WRITE(iff, *)
ENDIF
CLOSE(IFF)
1000 FORMAT(2(1x,I8   ))
1100 FORMAT(4(1x,E12.5))
2000 FORMAT(5(1x,E12.5))
!
END SUBROUTINE output_save_file_2d
