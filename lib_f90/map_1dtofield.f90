MODULE map_1dtofield
!
!  Contains routines to map a 1D line onto an N-dimensional field and reverse
!
INTERFACE maptofftfd
  MODULE PROCEDURE maptofftfd_1D_R_C_DP   ! 2D real    => complex Double Precision
  MODULE PROCEDURE maptofftfd_1D_C_C_DP   ! 2D complex => complex Double Precision
  MODULE PROCEDURE maptofftfd_2D_R_C_DP   ! 2D real    => complex Double Precision
  MODULE PROCEDURE maptofftfd_2D_C_C_DP   ! 2D complex => complex Double Precision
  MODULE PROCEDURE maptofftfd_3D_R_C_DP   ! 3D real    => complex Double Precision
  MODULE PROCEDURE maptofftfd_3D_C_C_DP   ! 3D complex => complex Double Precision
END INTERFACE maptofftfd
!
INTERFACE mapfftfdtoline
  MODULE PROCEDURE fftfdtoline_1D_C_R_DP  ! 2D complex => 1D REAL    Double Precision
  MODULE PROCEDURE fftfdtoline_1D_C_C_DP  ! 2D complex => 1D COMPLEX Double Precision
  MODULE PROCEDURE fftfdtoline_2D_C_R_DP  ! 2D complex => 1D REAL    Double Precision
  MODULE PROCEDURE fftfdtoline_2D_C_C_DP  ! 2D complex => 1D COMPLEX Double Precision
  MODULE PROCEDURE fftfdtoline_3D_C_R_DP  ! 3D complex => 1D REAL    Double Precision
  MODULE PROCEDURE fftfdtoline_3D_C_C_DP  ! 3D complex => 1D COMPLEX Double Precision
END INTERFACE mapfftfdtoline
!
!*******************************************************************************
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE maptofftfd_1D_R_C_DP(dimen, dsort, linear, fftfd)
!-
!   Map the 1D real linear array onto a 1D complex array. 
!   Shift the center to point (1,1,1)
!+
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dsort
REAL   (KIND=PREC_DP), DIMENSION(1:dimen(1))         , INTENT(IN) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)))  , INTENT(OUT) :: fftfd
!
INTEGER :: loop, ii
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
DO loop = 1, dimen(1)
  ii = mod((mod(loop-1, dimen(dsort(1)))                 + 1) + INT(dimen(dsort(1))/2) - l(dsort(1)), dimen(dsort(1))) + 1
  fftfd(ii   ) = CMPLX(linear(loop), 0.0D0)
ENDDO
!
END SUBROUTINE maptofftfd_1D_R_C_DP
!
!*******************************************************************************
!
SUBROUTINE maptofftfd_1D_C_C_DP(dimen, dsort, linear, fftfd)
!-
!   Map the 1D complex linear array onto a 1D complex array. 
!   Shift the center to point (1,1,1)
!+
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dsort
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(1))         , INTENT(IN) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)))  , INTENT(OUT) :: fftfd
!
INTEGER :: loop, ii
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
DO loop = 1, dimen(1)
  ii = mod((mod(loop-1, dimen(dsort(1)))                 + 1) + INT(dimen(dsort(1))/2) - l(dsort(1)), dimen(dsort(1))) + 1
  fftfd(ii   ) = linear(loop)
ENDDO
!
END SUBROUTINE maptofftfd_1D_C_C_DP
!
!*******************************************************************************
!
SUBROUTINE maptofftfd_2D_R_C_DP(dimen, dsort, linear, fftfd)
!-
!   Map the 1D linear array onto a 2D complex array. 
!   Shift the center to point (1,1,1)
!+
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dsort
REAL   (KIND=KIND(0.0D0)), DIMENSION(1:dimen(1)*dimen(2)                        &
                                               )     , INTENT(IN) :: linear
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(1:dimen(dsort(1)),                         &
                                     1:dimen(dsort(2))), INTENT(OUT) :: fftfd
!
INTEGER :: loop, ii,jj, i,j
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
DO loop = 1, dimen(1)*dimen(2)
   i = MOD((loop-1)/(dimen(2)*dimen(3)),dimen(1)) + 1            ! Indices: i == H
   j = MOD((loop-1)/(         dimen(3)),dimen(2)) + 1            ! Indices: j == K
   ientry(dsort(1)) = i                 ! i ==> goes into :dsort(1)
   ientry(dsort(2)) = j                 ! j ==> goes into :dsort(2)
   ii = MOD(ientry(1) + INT(dimen(dsort(1))/2) - l(dsort(1)), dimen(dsort(1))) + 1
   jj = MOD(ientry(2) + INT(dimen(dsort(2))/2) - l(dsort(2)), dimen(dsort(2))) + 1
   fftfd(ii,jj) = CMPLX(linear(loop), 0.0D0)
ENDDO
!
END SUBROUTINE maptofftfd_2D_R_C_DP
!
!*******************************************************************************
!
SUBROUTINE maptofftfd_2D_C_C_DP(dimen, dsort, linear, fftfd)
!-
!   Map the 1D complex array onto a 2D complex array. 
!   Shift the center to point (1,1,1)
!+
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dsort
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(1:dimen(1)*dimen(2)                        &
                                               )     , INTENT(IN) :: linear
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(1:dimen(dsort(1)),                         &
                                     1:dimen(dsort(2))), INTENT(OUT) :: fftfd
!
INTEGER :: loop, ii,jj, i,j
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
DO loop = 1, dimen(1)*dimen(2)
   i = MOD((loop-1)/(dimen(2)*dimen(3)),dimen(1)) + 1            ! Indices: i == H
   j = MOD((loop-1)/(         dimen(3)),dimen(2)) + 1            ! Indices: j == K
   ientry(dsort(1)) = i                 ! i ==> goes into :dsort(1)
   ientry(dsort(2)) = j                 ! j ==> goes into :dsort(2)
   ii = MOD(ientry(1) + INT(dimen(dsort(1))/2) - l(dsort(1)), dimen(dsort(1))) + 1
   jj = MOD(ientry(2) + INT(dimen(dsort(2))/2) - l(dsort(2)), dimen(dsort(2))) + 1
   fftfd(ii,jj) =         linear(loop)
ENDDO
!
END SUBROUTINE maptofftfd_2D_C_C_DP
!
!*******************************************************************************
!
!
!*******************************************************************************
!
SUBROUTINE maptofftfd_3D_R_C_DP(dimen, dsort, linear, fftfd)
!-
!   Map the 1D linear array onto a 3D complex array. 
!   Shift the center to point (1,1,1)
!+
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
REAL   (KIND=KIND(0.0D0)), DIMENSION(1:dimen(1)*dimen(2)*                       &
                                       dimen(3))     , INTENT(IN)  :: linear
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(1:dimen(dsort(1)),                         &
                                     1:dimen(dsort(2)),                         &
                                     1:dimen(dsort(3))), INTENT(OUT) :: fftfd
!
INTEGER :: loop, i,j,k, ii,jj, kk
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
l(dsort(3)) = MOD(dimen(dsort(3))-1,2)
DO loop = 1, dimen(1)*dimen(2)*dimen(3)
   i = MOD((loop-1)/(dimen(2)*dimen(3)),dimen(1)) + 1            ! Indices: i == H
   j = MOD((loop-1)/(         dimen(3)),dimen(2)) + 1            ! Indices: j == K
   k = MOD((loop-1)                    ,dimen(3)) + 1            ! Indices: k == L
   ientry(dsort(1)) = i                 ! i ==> goes into :dsort(1)
   ientry(dsort(2)) = j                 ! j ==> goes into :dsort(2)
   ientry(dsort(3)) = k                 ! k ==> goes into :dsort(3)
   ii = MOD(ientry(1) + INT(dimen(dsort(1))/2) - l(dsort(1)), dimen(dsort(1))) + 1
   jj = MOD(ientry(2) + INT(dimen(dsort(2))/2) - l(dsort(2)), dimen(dsort(2))) + 1
   kk = MOD(ientry(3) + INT(dimen(dsort(3))/2) - l(dsort(3)), dimen(dsort(3))) + 1
   fftfd(ii,jj, kk) = CMPLX(linear(loop), 0.0D0)
ENDDO
!
END SUBROUTINE maptofftfd_3D_R_C_DP
!
!*******************************************************************************
!
SUBROUTINE maptofftfd_3D_C_C_DP(dimen, dsort, linear, fftfd)
!-
!   Map the 1D linear array onto a 3D complex array. 
!   Shift the center to point (1,1,1)
!+
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3)                                , INTENT(IN) :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(1:dimen(1)*dimen(2)*                       &
                                       dimen(3))     , INTENT(IN)  :: linear
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(1:dimen(dsort(1)),                         &
                                     1:dimen(dsort(2)),                         &
                                     1:dimen(dsort(3))), INTENT(OUT) :: fftfd
!
INTEGER :: loop, i,j,k, ii,jj, kk
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
l(dsort(3)) = MOD(dimen(dsort(3))-1,2)
DO loop = 1, dimen(1)*dimen(2)*dimen(3)
   i = MOD((loop-1)/(dimen(2)*dimen(3)),dimen(1)) + 1            ! Indices: i == H
   j = MOD((loop-1)/(         dimen(3)),dimen(2)) + 1            ! Indices: j == K
   k = MOD((loop-1)                    ,dimen(3)) + 1            ! Indices: k == L
   ientry(dsort(1)) = i                 ! i ==> goes into :dsort(1)
   ientry(dsort(2)) = j                 ! j ==> goes into :dsort(2)
   ientry(dsort(3)) = k                 ! k ==> goes into :dsort(3)
   ii = MOD(ientry(1) + INT(dimen(dsort(1))/2) - l(dsort(1)), dimen(dsort(1))) + 1
   jj = MOD(ientry(2) + INT(dimen(dsort(2))/2) - l(dsort(2)), dimen(dsort(2))) + 1
   kk = MOD(ientry(3) + INT(dimen(dsort(3))/2) - l(dsort(3)), dimen(dsort(3))) + 1
   fftfd(ii,jj, kk) = linear(loop)
ENDDO
!
END SUBROUTINE maptofftfd_3D_C_C_DP
!
!*******************************************************************************
!
!  REVERSE ROUTINES
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_1D_C_R_DP(dimen, dsort, linear, fftfd)
!
!
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
REAL   (KIND=PREC_DP), DIMENSION(1:dimen(1)         ), INTENT(OUT) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, ii
!
DO ii = 1, dimen(dsort(1))
   loop =   + mod(ii-1+INT(dimen(dsort(1))/2) , dimen(dsort(1))) + 1
   linear(loop ) = REAL(fftfd(ii), KIND=PREC_DP)
ENDDO
!
END SUBROUTINE fftfdtoline_1D_C_R_DP
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_1D_C_C_DP(dimen, dsort, linear, fftfd)
!
!
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(1)         ), INTENT(OUT) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, ii
!
DO ii = 1, dimen(dsort(1))
   loop =   + mod(ii-1+INT(dimen(dsort(1))/2) , dimen(dsort(1))) + 1
   linear(loop ) =      fftfd(ii)
ENDDO
!
END SUBROUTINE fftfdtoline_1D_C_C_DP
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_2D_C_R_SP(dimen, dsort, linear, fftfd)
!
!
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
REAL   (KIND=PREC_DP), DIMENSION(1:dimen(1)*dimen(2)), INTENT(OUT) :: linear
COMPLEX(KIND=PREC_SP), DIMENSION(1:dimen(dsort(1)),                             &
                                 1:dimen(dsort(2)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, i,j, ii,jj
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
DO jj=1, dimen(dsort(2))
  DO ii = 1, dimen(dsort(1))
     ientry(1) = MOD(ii + INT(dimen(dsort(1))/2)             , dimen(dsort(1))) + 1
     ientry(2) = MOD(jj + INT(dimen(dsort(2))/2)             , dimen(dsort(2))) + 1
     i = ientry(dsort(1))                  ! i ==> goes into :dsort(1)
     j = ientry(dsort(2))                  ! j ==> goes into :dsort(2)
     loop = (i-1)*dimen(2) + j
     linear(loop ) = REAL(fftfd(ii,jj), KIND=PREC_DP)
  ENDDO
ENDDO
!
END SUBROUTINE fftfdtoline_2D_C_R_SP
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_2D_C_R_DP(dimen, dsort, linear, fftfd)
!+
! Backward mapping 2D Complex to 1D Real;   Double Precision
!-
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
REAL   (KIND=PREC_DP), DIMENSION(1:dimen(1)*dimen(2)), INTENT(OUT) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)),                             &
                                 1:dimen(dsort(2)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, i,j, ii,jj
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
DO jj=1, dimen(dsort(2))
  DO ii = 1, dimen(dsort(1))
     ientry(1) = MOD(ii + INT(dimen(dsort(1))/2)- l(dsort(1)) - 1, dimen(dsort(1))) + 1
     ientry(2) = MOD(jj + INT(dimen(dsort(2))/2) -l(dsort(2)) - 1, dimen(dsort(2))) + 1
     i = ientry(dsort(1))                  ! i ==> goes into :dsort(1)
     j = ientry(dsort(2))                  ! j ==> goes into :dsort(2)
     loop = (i-1)*dimen(2) + j
     linear(loop ) = REAL(fftfd(ii,jj), KIND=PREC_DP)
  ENDDO
ENDDO
!
END SUBROUTINE fftfdtoline_2D_C_R_DP
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_2D_C_C_DP(dimen, dsort, linear, fftfd)
!+
! Backward mapping 2D Complex to 1D Complex;  Double Precision
!-
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(1)*dimen(2)), INTENT(OUT) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)),                             &
                                 1:dimen(dsort(2)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, i,j, ii,jj
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
DO jj=1, dimen(dsort(2))
  DO ii = 1, dimen(dsort(1))
     ientry(1) = MOD(ii + INT(dimen(dsort(1))/2)- l(dsort(1)) - 1, dimen(dsort(1))) + 1
     ientry(2) = MOD(jj + INT(dimen(dsort(2))/2) -l(dsort(2)) - 1, dimen(dsort(2))) + 1
     i = ientry(dsort(1))                  ! i ==> goes into :dsort(1)
     j = ientry(dsort(2))                  ! j ==> goes into :dsort(2)
     loop = (i-1)*dimen(2) + j
     linear(loop ) = fftfd(ii,jj)
  ENDDO
ENDDO
!
END SUBROUTINE fftfdtoline_2D_C_C_DP
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_3D_C_R_DP(dimen, dsort, linear, fftfd)
!
!
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
REAL   (KIND=PREC_DP), DIMENSION(1:dimen(1)*dimen(2)*dimen(3))     , INTENT(OUT) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)),                             &
                                 1:dimen(dsort(2)),                              &
                                 1:dimen(dsort(3)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, i,j,k, ii,jj, kk
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
l(dsort(3)) = MOD(dimen(dsort(3))-1,2)
!
DO kk=1, dimen(dsort(3))
   DO jj=1, dimen(dsort(2))
      DO ii = 1, dimen(dsort(1))
         ientry(1) = MOD(ii + INT(dimen(dsort(1))/2)- l(dsort(1)) - 1, dimen(dsort(1))) + 1
         ientry(2) = MOD(jj + INT(dimen(dsort(2))/2) -l(dsort(2)) - 1, dimen(dsort(2))) + 1
         ientry(3) = MOD(kk + INT(dimen(dsort(3))/2) -l(dsort(3)) - 1, dimen(dsort(3))) + 1
         i = ientry(dsort(1))                  ! i ==> goes into :dsort(1)
         j = ientry(dsort(2))                  ! j ==> goes into :dsort(2)
         k = ientry(dsort(3))                  ! k ==> goes into :dsort(3)
         loop = (i-1)*dimen(2)*dimen(3) + (j-1)*dimen(3) + k
         linear(loop ) = REAL(fftfd(ii,jj, kk), KIND=PREC_DP)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE fftfdtoline_3D_C_R_DP
!
!*******************************************************************************
!
SUBROUTINE fftfdtoline_3D_C_C_DP(dimen, dsort, linear, fftfd)
!
!
USE precision_mod
!
IMPLICIT NONE
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dimen
INTEGER, DIMENSION(3)                                , INTENT(IN)  :: dsort
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(1)*dimen(2)*dimen(3))     , INTENT(OUT) :: linear
COMPLEX(KIND=PREC_DP), DIMENSION(1:dimen(dsort(1)),                             &
                                 1:dimen(dsort(2)),                              &
                                 1:dimen(dsort(3)))  , INTENT(IN) :: fftfd
!
INTEGER :: loop, i,j,k, ii,jj, kk
INTEGER, DIMENSION(3) :: l              ! Additional shift of 1 for even dimensions
INTEGER, DIMENSION(3) :: ientry         ! Target index for i,j
!
l(dsort(1)) = MOD(dimen(dsort(1))-1,2)
l(dsort(2)) = MOD(dimen(dsort(2))-1,2)
l(dsort(3)) = MOD(dimen(dsort(3))-1,2)
!
DO kk=1, dimen(dsort(3))
   DO jj=1, dimen(dsort(2))
      DO ii = 1, dimen(dsort(1))
         ientry(1) = MOD(ii + INT(dimen(dsort(1))/2)- l(dsort(1)) - 1, dimen(dsort(1))) + 1
         ientry(2) = MOD(jj + INT(dimen(dsort(2))/2) -l(dsort(2)) - 1, dimen(dsort(2))) + 1
         ientry(3) = MOD(kk + INT(dimen(dsort(3))/2) -l(dsort(3)) - 1, dimen(dsort(3))) + 1
         i = ientry(dsort(1))                  ! i ==> goes into :dsort(1)
         j = ientry(dsort(2))                  ! j ==> goes into :dsort(2)
         k = ientry(dsort(3))                  ! k ==> goes into :dsort(3)
         loop = (i-1)*dimen(2)*dimen(3) + (j-1)*dimen(3) + k
         linear(loop ) = fftfd(ii,jj, kk)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE fftfdtoline_3D_C_C_DP
!
!*******************************************************************************
!
END MODULE map_1dtofield 
