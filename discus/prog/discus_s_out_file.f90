module discus_output_save_mod
!
contains
!
! Version for DISCUS SUITE
! The file is may be written directly into kuplot if the
! filename starts with 'kuplot'
!
!  Separate routines for 1D and 2D files
!
SUBROUTINE output_save_file_1d( outfile, npkt1, xwrt, ywrt, mode )
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
USE lib_length
use precision_mod
USE support_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*) ,                     INTENT(IN) :: outfile
INTEGER           ,                     INTENT(IN) :: npkt1
REAL(kind=PREC_DP), DIMENSION(1:npkt1), INTENT(IN) :: xwrt
REAL(kind=PREC_DP), DIMENSION(1:npkt1), INTENT(IN) :: ywrt
integer                               , intent(in) :: mode
!
integer, parameter :: NEW = 0
integer, parameter :: OLD = 1
integer, parameter :: ADD = 2
INTEGER, PARAMETER :: IFF = 2
INTEGER            :: i, j                 ! Dummy indices
!
INTEGER   :: lname
LOGICAL   :: lkuplot
INTEGER   :: nr
INTEGER   :: maxpp
INTEGER   :: ik
logical   :: lold_exist
!
!
lname   = len_str(outfile)
lkuplot = .false.
!mode  = ADD
ik = 1
!
IF(lname >= 6) THEN   ! File name long enough for 'kuplot' string?
   IF(outfile(1:6) == 'kuplot') THEN
     lkuplot = .TRUE.
   ENDIF
ENDIF
IF(lkuplot) THEN      ! 'write' into kuplot array
   lold_exist = .false.
   loop_fname: do j=1, iz-1            ! Loop to search for old file
      if(fname(j) == outfile) then
         lold_exist = .true.
         ik = j
         exit loop_fname
      endif
   enddo loop_fname
   if((mode  ==   NEW) .or. (.not. lold_exist)) then                  ! Create a new file
      maxpp = maxarray - offxy (iz - 1)   ! available points in kuplot
      IF(npkt1 > maxpp) THEN  ! Too many points abort
         ier_num = -29
         ier_typ =ER_IO
         RETURN
      ENDIF
      nr = 1
      DO i=1, npkt1
         x (offxy (iz - 1) + nr) = xwrt(i)
         y (offxy (iz - 1) + nr) = ywrt(i)
         dx(offxy (iz - 1) + nr) = 0.00
         dy(offxy (iz - 1) + nr) = 1.00
         nr = nr + 1
      ENDDO
      lenc  (iz) = nr - 1                      ! set length
      offxy (iz) = offxy (iz - 1) + lenc(iz) ! set offset
      offz  (iz) = offz (iz - 1)
      xmin  (iz) = MINVAL(xwrt)
      xmax  (iz) = MAXVAL(xwrt)
      ymin  (iz) = MINVAL(ywrt)
      ymax  (iz) = MAXVAL(ywrt)
      iz = iz + 1                            ! increment number of data sets
      ik = iz - 1
      fname(ik) = outfile                    ! store filename
   else                                   ! overwrite existing file
      if(lenc(ik) == npkt1) then
         if(mode==OLD) then               ! Replace old data set
            nr = 1
            do i=1, npkt1
               x (offxy(ik-1) + nr) = xwrt(i)
               y (offxy(ik-1) + nr) = ywrt(i)
               dx(offxy(ik-1) + nr) = 0.0
               dy(offxy(ik-1) + nr) = 1.0
               nr = nr + 1
            enddo
         elseif(mode==ADD) then           ! Add into old data set
            nr = 1
            do i=1, npkt1
               y (offxy(ik-1) + nr) = y (offxy(ik-1) + nr) + ywrt(i)
               nr = nr + 1
            enddo
         endif
      else
         ier_num = -176
         ier_typ = ER_APPL
      endif
   endif
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
!*******************************************************************************
!
SUBROUTINE output_save_file_2d( outfile, ranges, npkt1, npkt2, zwrt,&
           header_lines, nheader, mode)
!-
!  Write 2D output file
!+
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
USE lib_length
USE support_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*),           INTENT(IN) :: outfile
REAL(kind=PREC_DP)   , DIMENSION(1:4),     INTENT(IN) :: ranges
INTEGER,                     INTENT(IN) :: npkt1
INTEGER,                     INTENT(IN) :: npkt2
REAL(kind=PREC_DP)   , DIMENSION(1:npkt1, 1:npkt2), INTENT(IN) :: zwrt
INTEGER,                                 INTENT(IN) :: nheader! number of lines in header
CHARACTER (LEN=160), DIMENSION(nheader), INTENT(IN) :: header_lines
integer , intent(in) :: mode
!
integer, parameter :: NEW = 0
integer, parameter :: OLD = 1
integer, parameter :: ADD = 2
INTEGER, PARAMETER :: IFF = 2
INTEGER            :: i, j
integer            :: ik       ! Data set number for old
!
INTEGER   :: lname
LOGICAL   :: lkuplot
logical   :: lold_exist
INTEGER   :: maxpkt
INTEGER   :: maxzz
REAL(kind=PREC_DP)      :: dxx, dyy
!
!
lname   = len_str(outfile)
lkuplot = .false.
ik      = 1
!
IF(lname >= 6) THEN   ! File name long enough for 'kuplot' string?
   IF(outfile(1:6) == 'kuplot') THEN
     lkuplot = .TRUE.
   ENDIF
ENDIF
IF(lkuplot) THEN      ! 'write' into kuplot array
   lold_exist = .false.
   loop_fname: do j=1, iz-1            ! Loop to search for old file
      if(fname(j) == outfile) then
         lold_exist = .true.
         ik = j
         exit loop_fname
      endif
   enddo loop_fname
   if((mode  ==   NEW) .or. (.not. lold_exist)) then                  ! Create a new file
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
      dxx = (xmax (iz) - xmin (iz) ) / REAL(nx (iz) - 1) 
      dyy = (ymax (iz) - ymin (iz) ) / REAL(ny (iz) - 1) 
      DO i = 1, nx (iz) 
         x(offxy(iz - 1) + i) = xmin(iz) + (i-1)*dxx 
      ENDDO 
      DO i = 1, ny (iz) 
         y(offxy(iz - 1) + i) = ymin(iz) + (i-1)*dyy 
      ENDDO 
      lni(iz)   = .true.
      lenc(iz)   = MAX(nx(iz),ny(iz))
      offxy(iz) = offxy(iz-1) + lenc(iz)
      offz (iz) = offz (iz-1) + nx(iz)*ny(iz)
      iz        = iz + 1
      zmax (iz-1) = MAXVAL(zwrt)
      zmin (iz-1) = MINVAL(zwrt)
      fname(iz-1) = outfile                    ! store filename
   else                                   ! overwrite existing file
      if(nx(ik) == npkt1 .and. ny(ik) == npkt2) then
         if(mode==OLD) then               ! Replace old data set
            do j=1,npkt2
               do i=1,npkt1
                  z(offz(ik-1) + (i-1)*ny(ik)+j) = zwrt(i,j)
               enddo
            enddo
         elseif(mode==ADD) then           ! Add into old data set
            do j=1,npkt2
               do i=1,npkt1
                  z(offz(ik-1) + (i-1)*ny(ik)+j) = z(offz(ik-1) + (i-1)*ny(ik)+j) + zwrt(i,j)
               enddo
            enddo
         endif
      else
         ier_num = -176
         ier_typ = ER_APPL
      endif
   endif
ELSE
!
   CALL oeffne(IFF, outfile, 'unknown')
   IF(ier_num == 0) THEN
!!!!!!!!!!!!!!!!!!!      CALL write_discus_nipl_header(iff)
      DO i=1, nheader
         WRITE(iff, '(a)') header_lines(i)(1:LEN_TRIM(header_lines(i)))
      ENDDO
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
!
end module discus_output_save_mod
