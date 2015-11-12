!
! Version for DISCUS SUITE
! The file is may be written directly into kuplot if the
! filename starts with 'kuplot'
!
!  Separate routines for 1D and 2D files
!
SUBROUTINE output_save_file_1d( outfile, npkt1, xwrt, ywrt )
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*),           INTENT(IN) :: outfile
INTEGER,                     INTENT(IN) :: npkt1
REAL   , DIMENSION(1:npkt1), INTENT(IN) :: xwrt
REAL   , DIMENSION(1:npkt1), INTENT(IN) :: ywrt
!
INTEGER, PARAMETER :: IFF = 2
INTEGER            :: i
!
INTEGER   :: lname
LOGICAL   :: lkuplot
INTEGER   :: nr
INTEGER   :: maxpp
INTEGER   :: ik
!
INTEGER   :: len_str
!
lname   = len_str(outfile)
lkuplot = .false.
!
IF(lname >= 6) THEN   ! File name long enough for 'kuplot' string?
   IF(outfile(1:6) == 'kuplot') THEN
     lkuplot = .TRUE.
   ENDIF
ENDIF
IF(lkuplot) THEN      ! 'write' into kuplot array
   maxpp = maxarray - offxy (iz - 1)   ! available points in kuplot
   IF(npkt1 > maxpp) THEN  ! Too many points abort
      ier_num = -29
      ier_typ =ER_IO
      RETURN
   ENDIF
   DO i=1, npkt1
         WRITE(iff, 1000) xwrt(i), ywrt(i)
      x (offxy (iz - 1) + nr) = xwrt(i)
      y (offxy (iz - 1) + nr) = ywrt(i)
      dx(offxy (iz - 1) + nr) = 0.00
      dy(offxy (iz - 1) + nr) = 1.00
      nr = nr + 1
   ENDDO
   len   (iz) = nr - 1                      ! set length
   offxy (iz) = offxy (iz - 1) + len (iz) ! set offset
   offz  (iz) = offz (iz - 1)
   xmin  (iz) = MINVAL(xwrt)
   xmax  (iz) = MAXVAL(xwrt)
   ymin  (iz) = MINVAL(ywrt)
   ymax  (iz) = MAXVAL(ywrt)
   iz = iz + 1                            ! increment number of data sets
   ik = iz - 1
   fname(ik) = outfile                    ! store filename
ELSE                                      ! normal write to disk
   CALL oeffne(IFF, outfile, 'unknown')
   IF(ier_num == 0) THEN
      DO i=1, npkt1
         WRITE(iff, 1000) xwrt(i), ywrt(i)
      ENDDO
   ENDIF
   CLOSE(IFF)
ENDIF
1000 FORMAT(2(1x,E12.5))
!
END SUBROUTINE output_save_file_1d
!
!  2D output file
!
SUBROUTINE output_save_file_2d( outfile, ranges, npkt1, npkt2, zwrt)
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*),           INTENT(IN) :: outfile
REAL   , DIMENSION(1:4),     INTENT(IN) :: ranges
INTEGER,                     INTENT(IN) :: npkt1
INTEGER,                     INTENT(IN) :: npkt2
REAL   , DIMENSION(1:npkt1, 1:npkt2), INTENT(IN) :: zwrt
!
INTEGER, PARAMETER :: IFF = 2
INTEGER            :: i, j
!
INTEGER   :: lname
LOGICAL   :: lkuplot
INTEGER   :: nr
INTEGER   :: maxpkt
INTEGER   :: maxzz
INTEGER   :: ik 
REAL      :: dxx, dyy
!
INTEGER   :: len_str
!
lname   = len_str(outfile)
lkuplot = .false.
!
IF(lname >= 6) THEN   ! File name long enough for 'kuplot' string?
   IF(outfile(1:6) == 'kuplot') THEN
     lkuplot = .TRUE.
   ENDIF
ENDIF
IF(lkuplot) THEN      ! 'write' into kuplot array
   maxpkt = maxarray - offxy(iz-1)
   maxzz  = maxarray - offz (iz-1)
   IF(MAX(npkt1,npkt2) > maxzz  .AND.  &
      MAX(npkt1,npkt2) > maxpkt) THEN  ! Too many points abort
      ier_num = -29
      ier_typ =ER_IO
      RETURN
   ENDIF
   nx(iz)   = npkt1
   ny(iz)   = npkt2
   xmin(iz) = ranges(1)
   xmax(iz) = ranges(2)
   ymin(iz) = ranges(3)
   ymax(iz) = ranges(4)
   DO j=1,npkt2
      DO i=1,npkt1
         z(offz(iz-1) + (i-1)*ny(iz)+j) = zwrt(i,j)
      ENDDO
   ENDDO
!                                                                       
!------ set values for X and Y                                          
!                                                                       
   dxx = (xmax (iz) - xmin (iz) ) / float (nx (iz) - 1) 
   dyy = (ymax (iz) - ymin (iz) ) / float (ny (iz) - 1) 
   DO i = 1, nx (iz) 
      x(offxy(iz - 1) + i) = xmin(iz) + (i-1)*dxx 
   ENDDO 
   DO i = 1, ny (iz) 
      y(offxy(iz - 1) + i) = ymin(iz) + (i-1)*dyy 
   ENDDO 
   lni(iz)   = .true.
   len(iz)   = MAX(nx(iz),ny(iz))
   offxy(iz) = offxy(iz-1) + len(iz)
   offz (iz) = offz (iz-1) + nx(iz)*ny(iz)
   iz        = iz + 1
   zmax (iz-1) = MAXVAL(zwrt)
   zmin (iz-1) = MINVAL(zwrt)
   fname(iz-1) = outfile                    ! store filename
ELSE
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
ENDIF
1000 FORMAT(2(1x,I8   ))
1100 FORMAT(4(1x,E12.5))
2000 FORMAT(5(1x,E12.5))
!
END SUBROUTINE output_save_file_2d
