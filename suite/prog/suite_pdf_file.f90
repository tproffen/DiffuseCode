!
! Version for DISCUS_SUITE
! The file is may be written directly into kuplot if the
! filename starts with 'kuplot'
!
SUBROUTINE pdf_save_file(cdummy, pdf_rfmin, pdf_rfmax, pdf_deltar, pdf_us_int,&
                         pdf_calc_l, pdf_calc_u, pdf_skal,pdf_calc)
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: cdummy
REAL,              INTENT(IN) :: pdf_rfmin
REAL,              INTENT(IN) :: pdf_rfmax
REAL,              INTENT(IN) :: pdf_deltar
INTEGER,           INTENT(IN) :: pdf_us_int
INTEGER,           INTENT(IN) :: pdf_calc_l
INTEGER,           INTENT(IN) :: pdf_calc_u
REAL,              INTENT(IN) :: pdf_skal
REAL(KIND(0.d0)), DIMENSION(pdf_calc_l:pdf_calc_u), INTENT(IN) :: pdf_calc
!
INTEGER   :: nmi
INTEGER   :: nma
INTEGER   :: nmd
INTEGER   :: i
INTEGER   :: lname
LOGICAL   :: lkuplot
REAL      :: r
!
INTEGER   :: nr
INTEGER   :: maxpp
INTEGER   :: ik
!
INTEGER :: len_str
!
lname   = len_str(cdummy)
lkuplot = .false.
!
IF(lname >= 6) THEN   ! File name long enough for 'kuplot' string?
   IF(cdummy(1:6) == 'kuplot') THEN
     lkuplot = .TRUE.
   ENDIF
ENDIF
IF(lkuplot) THEN      ! 'write' into kuplot array
   nr = 1
   maxpp = maxarray - offxy (iz - 1)   ! available points in kuplot
   nmi = nint (pdf_rfmin / pdf_deltar) 
   nma = nint (pdf_rfmax / pdf_deltar) 
   nmd = pdf_us_int   ! step width = (delta r user)/(deltar internal)
   IF(INT((nma-nmi+1)/nmd) > maxpp) THEN  ! Too many points abort
      ier_num = -29
      ier_typ =ER_IO
      RETURN
   ENDIF
   DO i = nmi, nma, nmd                   ! Write all data points
      r = float (i) * pdf_deltar 
      x (offxy (iz - 1) + nr) = r
      y (offxy (iz - 1) + nr) = pdf_calc(i)
      dx (offxy (iz - 1) + nr) = 0.0
      dy (offxy (iz - 1) + nr) = 1.0
      nr = nr + 1
   ENDDO
   len (iz) = nr - 1                      ! set length
   offxy (iz) = offxy (iz - 1) + len (iz) ! set offset
   offz (iz) = offz (iz - 1)
   iz = iz + 1                            ! increment number of data sets
   ik = iz - 1
   CALL get_extrema_xy (x, ik, len (ik), xmin, xmax)
   fname(ik) = cdummy                     ! store filename
ELSE                                      ! Normal write to disk
   CALL oeffne (57, cdummy, 'unknown') 
   IF (ier_num.eq.0) then 
      nmi = nint (pdf_rfmin / pdf_deltar) 
      nma = nint (pdf_rfmax / pdf_deltar) 
      nmd = pdf_us_int   ! step width = (delta r user)/(deltar internal)
      DO i = nmi, nma, nmd 
         r = float (i) * pdf_deltar 
         WRITE (57, 5000) r, pdf_skal * pdf_calc (i), 0.0, 1.0 
      ENDDO 
      CLOSE (57) 
   ENDIF 
ENDIF 
!
5000 FORMAT (F9.4,3X,F21.10,5X,2(F6.2,1X))
!
END SUBROUTINE pdf_save_file
