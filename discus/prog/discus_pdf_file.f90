!
! Version for stan alone DISCUS
! The file is always written to hard disk
!
SUBROUTINE pdf_save_file(cdummy, pdf_rfmin, pdf_rfmax, pdf_deltar, pdf_us_int,&
                         pdf_calc_l, pdf_calc_u, pdf_skal,pdf_calc)
!
USE discus_config_mod
USE errlist_mod
USE precision_mod
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
REAL(PREC_DP), DIMENSION(pdf_calc_l:pdf_calc_u), INTENT(IN) :: pdf_calc
!
INTEGER   :: nmi
INTEGER   :: nma
INTEGER   :: nmd
INTEGER   :: i
REAL      :: r
!
CALL oeffne (57, cdummy, 'unknown') 
IF (ier_num.eq.0) then 
   nmi = nint (pdf_rfmin / pdf_deltar) 
   nma = nint (pdf_rfmax / pdf_deltar) 
   nmd = pdf_us_int   ! step width = (delta r user)/(deltar internal)
   DO i = nmi, nma, nmd 
      r = REAL(i) * pdf_deltar 
      WRITE (57, 5000) r, pdf_skal * pdf_calc (i), 0.0, 1.0 
   ENDDO 
   CLOSE (57) 
ENDIF 
!
5000 FORMAT (F9.4,3X,F21.10,5X,2(F6.2,1X))
!
END SUBROUTINE pdf_save_file
