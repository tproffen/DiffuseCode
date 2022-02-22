MODULE kuplot_toglobal
!-
!   Transfer a data set to/from the globale storage
!+
USE kuplot_mod
USE global_data_mod
!
PRIVATE
PUBLIC kuplot_to_global
!
!*******************************************************************************
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE kuplot_to_global(line) 
!
USE kuplot_config
USE kuplot_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*) , INTENT(INOUT) :: line         ! Command line
!
INTEGER, PARAMETER :: MAXW = 2
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara    ! Parameter strings
INTEGER            , DIMENSION(MAXW) :: lpara    ! length of each parameter strign
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte    ! Parameter values
!

INTEGER :: ik      ! Data set to be transfered to ig
INTEGER :: ig      ! Number in global data set
INTEGER :: length  ! Number in global data set
INTEGER :: ianz    ! Number in global data set
INTEGER :: iianz   ! Number in global data set
INTEGER :: nnpara   ! Number refined parameters used by refine
INTEGER :: nnfix    ! Number fixed   parameters used by refine
INTEGER, DIMENSION(3) :: dimen
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: ext_data
!
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: O_REFINE  = 1
INTEGER, PARAMETER :: O_KUPL    = 2
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate
!
DATA oname  / 'refine' , 'kuplot' /
DATA loname /  6       ,  6     /
opara  =  (/ 'cost    ', 'last    '  /)
lopara =  (/ 4         ,  4          /)
owerte =  (/ -1.000000 ,  -1.0000000 /)
!
CALL gl_get_npara(NNPARA, NNFIX)      ! Check number of refined parameters
IF(NNPARA==0) RETURN                  ! if zero, silently leave
!
length = LEN_TRIM(line)
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num/=0) RETURN
!
ik = iz-1                               ! default to last data set
IF(lpresent(O_KUPL)) THEN
   IF(opara(O_KUPL)=='last') THEN
      ik = iz-1
   ELSE
      cpara(1) = opara(O_KUPL)
      lpara(1) = lopara(O_KUPL)
      iianz = 1
      CALL ber_params(iianz, cpara, lpara, werte, MAXW)
      ik = NINT(werte(1))
   ENDIF
ENDIF
!
ig = 0                                  ! default to last data set
IF(lpresent(O_REFINE)) THEN
   IF(opara(O_REFINE)=='cost') THEN
      ig = 0
   elseif(opara(O_REFINE)=='obs' .or. opara(O_REFINE)=='exp') THEN
      ig = -2
   elseif(opara(O_REFINE)=='sigma') then
      ig = -3
   elseif(opara(O_REFINE)=='opti') then
      ig = -1
   ELSE
      cpara(1) = opara(O_REFINE)
      lpara(1) = lopara(O_REFINE)
      iianz = 1
      CALL ber_params(iianz, cpara, lpara, werte, MAXW)
      ig = NINT(werte(1))
   ENDIF
ENDIF
!
IF(ig<-3          .OR. ig> nnpara) THEN
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Global data set number outside limits'
   RETURN
ELSE
   IF(ig<0) RETURN          ! Fixed data set silently ignore
ENDIF
!
IF(lni(ik)) THEN            ! Data set is 2D (NIPL)
ELSE
   dimen(1) = lenc(ik)
   ALLOCATE(ext_data(1:lenc(ik)))
   ext_data(1:lenc(ik)) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
   CALL gl_set_data(dimen(1), NNPARA, ig, ext_data)
   DEALLOCATE(ext_data)
ENDIF
END SUBROUTINE kuplot_to_global
!
!*******************************************************************************
!
END MODULE kuplot_toglobal
