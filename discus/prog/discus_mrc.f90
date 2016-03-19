MODULE discus_mrc
!
USE errlist_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC mrc_write
!
!
CONTAINS
!
!
   SUBROUTINE mrc_write ( value, laver)
!
   USE crystal_mod
   USE diffuse_mod
   USE output_mod
   USE qval_mod
!
   IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: value
   LOGICAL, INTENT(IN) :: laver
!
   INTEGER, PARAMETER :: IMRC = 97
!
   INTEGER, PARAMETER :: mode = 2
!
   INTEGER, PARAMETER :: nxstart = 0
   INTEGER, PARAMETER :: nystart = 0
   INTEGER, PARAMETER :: nzstart = 0
   INTEGER            :: i, j,l, irec
   REAL               :: sqq
!
   CALL no_error
!
   OPEN(UNIT=IMRC,FILE=outfile,STATUS='unknown', &
        FORM='unformatted', ACCESS='direct',RECL=4)
!
!  Write header is identical to all 
!
   WRITE(IMRC,REC= 1) inc(1)     ! Dimension column
   WRITE(IMRC,REC= 2) inc(2)     ! Dimension row
   WRITE(IMRC,REC= 3) inc(3)     ! Dimension slices
   WRITE(IMRC,REC= 4) mode       ! Data type = 2 32-bit real
   WRITE(IMRC,REC= 5) nxstart    ! Number of first column in map = 0
   WRITE(IMRC,REC= 6) nystart    ! Number of first row    in map = 0
   WRITE(IMRC,REC= 7) nzstart    ! Number of first slice  in map = 0
   WRITE(IMRC,REC= 8) inc(1)-1   ! Number of intervals along x
   WRITE(IMRC,REC= 9) inc(2)-1   ! Number of intervals along y
   WRITE(IMRC,REC=10) inc(3)-1   ! Number of intervals along z
   WRITE(IMRC,REC=11) cr_ar(1)   ! Reciprocal lattice parameter a
   WRITE(IMRC,REC=12) cr_ar(2)   ! Reciprocal lattice parameter b
   WRITE(IMRC,REC=13) cr_ar(3)   ! Reciprocal lattice parameter c
   WRITE(IMRC,REC=14) cr_wrez(1) ! Reciprocal lattice parameter alpha
   WRITE(IMRC,REC=15) cr_wrez(2) ! Reciprocal lattice parameter beta
   WRITE(IMRC,REC=16) cr_wrez(3) ! Reciprocal lattice parameter gamma
   WRITE(IMRC,REC=17) extr_abs   ! Reciprocal axis along columns 
   WRITE(IMRC,REC=18) extr_ord   ! Reciprocal axis along rows 
   WRITE(IMRC,REC=19) extr_top   ! Reciprocal axis along sections 
   WRITE(IMRC,REC=20) 0.0        ! Minimum density
   WRITE(IMRC,REC=21) 0.0        ! Maximum density
   WRITE(IMRC,REC=22) 0.0        ! Average density
   WRITE(IMRC,REC=23) 0          ! 0 = image
   WRITE(IMRC,REC=24) 0          ! No extended header
   DO I = 25,49
      WRITE(IMRC,REC= i) 0          ! 0
   ENDDO
   WRITE(IMRC,REC=50) 0.0        ! Origin
   WRITE(IMRC,REC=51) 0.0        ! Origin
   WRITE(IMRC,REC=52) 0.0        ! Origin
   DO I = 53,256
      WRITE(IMRC,REC= i) 0          ! 0
   ENDDO
!
   irec = 1024
   DO l = 1, out_inc(3)
      DO j = 1, out_inc(2)
         DO i = 1, out_inc(1)
            irec = irec + 1
            sqq       =  qval ( (i - 1) * out_inc(3)*out_inc (2) +        &
                                (j - 1) * out_inc(3)             + l,     &
                                value,  i, j, laver)
            WRITE(IMRC,rec=irec) sqq
         ENDDO
      ENDDO
   ENDDO
!
   CLOSE(IMRC)
!
   END SUBROUTINE mrc_write
!
END MODULE discus_mrc
