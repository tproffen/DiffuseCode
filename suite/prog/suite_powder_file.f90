!
! OBSOLETE
! content has been moved to discus_s_powder_file.f90 
! Reference has been removed from suite CMakeLists.txt
!*******************************************************************************
!
! Version for stand alone DISCUS
! The file is may be written directly into kuplot if the
! filename starts with 'kuplot'
!
SUBROUTINE powder_do_write (outfile, npkt_wrt, xwrt, ywrt)
!
USE kuplot_config
USE kuplot_mod
use kuplot_extrema_mod
USE errlist_mod
USE lib_length
use precision_mod
USE support_mod
IMPLICIT NONE
!
CHARACTER (LEN=*)                , INTENT(IN) :: outfile
INTEGER                          , INTENT(IN) :: npkt_wrt
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_wrt  ) , INTENT(IN) :: xwrt
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_wrt  ) , INTENT(IN) :: ywrt
!
INTEGER, PARAMETER                            :: iff = 2
INTEGER :: ii
INTEGER :: lname
LOGICAL :: lkuplot
!
INTEGER :: nr
INTEGER :: maxpp
INTEGER :: ik
!
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
   nr = 1
   maxpp = maxarray - offxy (iz - 1)   ! available points in kuplot
   IF(npkt_wrt > maxpp) THEN           ! Too many points abort
      ier_num = -29
      ier_typ =ER_IO
      RETURN
   ENDIF
   DO ii = 0,npkt_wrt
      x (offxy (iz - 1) + nr)  = xwrt(ii)
      y (offxy (iz - 1) + nr)  = ywrt(ii)
      dx (offxy (iz - 1) + nr) = 0.0
      dy (offxy (iz - 1) + nr) = 1.0
      nr = nr + 1
   ENDDO
   lenc(iz) = nr - 1                      ! set length
   offxy (iz) = offxy (iz - 1) + lenc(iz) ! set offset
   offz (iz) = offz (iz - 1)
   iz = iz + 1                            ! increment number of data sets
   ik = iz - 1
   CALL get_extrema_xy (x, ik, lenc(ik), xmin, xmax)
   fname(ik) = outfile                    ! store filename
ELSE
   CALL oeffne (iff, outfile, 'unknown') 
   IF(ier_num == 0) THEN
      DO ii = 0,npkt_wrt
         WRITE( iff, *) xwrt(ii),ywrt(ii)
      ENDDO
      CLOSE(iff)
   ENDIF
ENDIF
!
END SUBROUTINE powder_do_write
