MODULE transfrm_mod
!+
!
!     variables needed for the unit cell transformations
!-
!
INTEGER, PRIVATE    :: ik
SAVE
!
INTEGER, PARAMETER  ::       TRAN_INP_F   =  0
INTEGER, PARAMETER  ::       TRAN_INP_FI  =  1
INTEGER, PARAMETER  ::       TRAN_INP_G   =  2
INTEGER, PARAMETER  ::       TRAN_INP_GI  =  3
!
INTEGER             ::  TRAN_MAXSCAT = 1
!
INTEGER             ::  tran_start   =  1
INTEGER             ::  tran_end     = -1
INTEGER             ::  tran_inp     = TRAN_INP_G
LOGICAL, DIMENSION(:), ALLOCATABLE  ::  tran_latom  ! (0:TRAN_MAXSCAT)
!
LOGICAL             ::  tran_oold      = .true.
LOGICAL             ::  tran_sel_atom  = .true.
REAL                ::  tran_orig(3)   = 0.0
REAL                ::  tran_det       = 1.0
REAL                ::  tran_g   (4,4) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_g ))
REAL                ::  tran_gi  (4,4) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_gi))
REAL                ::  tran_f   (4,4) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_f ))
REAL                ::  tran_fi  (4,4) = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(tran_fi))
REAL                ::  tran_deltahkl  = 0.001
!
INTEGER             ::  tran_size_of   = 0 ! Bytes alloacted for transformation
!
END MODULE transfrm_mod
