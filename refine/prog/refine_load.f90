MODULE refine_load_mod
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_load(line, length)
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
INTEGER                              :: ndata
INTEGER                              :: nsigma
!
LOGICAL, EXTERNAL :: str_comp
!
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_SIGMA   = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate
!
DATA oname  / 'sigma ' /
DATA loname /  5       /
opara  =  (/ '0.000000'/)   ! Always provide fresh default values
lopara =  (/  8        /)
owerte =  (/  0.000000 /)
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(IANZ>1) THEN
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num/=0) RETURN
nsigma = NINT(owerte(O_SIGMA))
!
IF(str_comp (cpara(1), 'kuplot', 3, lpara(1), 6) ) THEN
   cpara(1) = '0'
   lpara(2) = 1
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/= 0) RETURN
   ndata = NINT(werte(2))
   CALL refine_load_kuplot(ndata, nsigma)
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   RETURN
ENDIF
!
END SUBROUTINE refine_load
!
!*******************************************************************************
!
SUBROUTINE refine_load_kuplot(ndata, nsigma)
!
! Transfers data set no. ndata from KUPLOT into REFINE memory
!
USE refine_data_mod
!
USE kuplot_mod
!
USE errlist_mod
USE define_variable_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ndata           ! no of data set to be transfered
INTEGER, INTENT(IN) :: nsigma          ! no of data set with sigma's to be transfered
!
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE. ! Prevents user from deleting variables
INTEGER :: ix, iy                      ! Dummy loop variables
REAL    :: step
!
IF(ndata<1 .OR. ndata>(iz - 1) ) THEN
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Data set number outside KUPLOT range'
   RETURN
ENDIF
!
IF(ALLOCATED(ref_data))   DEALLOCATE(ref_data)
IF(ALLOCATED(ref_weight)) DEALLOCATE(ref_weight)
IF(ALLOCATED(ref_x     )) DEALLOCATE(ref_x     )
IF(ALLOCATED(ref_y     )) DEALLOCATE(ref_y     )
!
IF(lni(ndata)) THEN                    ! 2D data set
   ref_dim(1) = nx(ndata )
   ref_dim(2) = ny(ndata )
   ALLOCATE(ref_data  (ref_dim(1),ref_dim(2)))
   ALLOCATE(ref_weight(ref_dim(1),ref_dim(2)))
   ALLOCATE(ref_x     (ref_dim(1)))
   ALLOCATE(ref_y     (ref_dim(2)))
   DO iy=1,ref_dim(2)
      DO ix=1,ref_dim(1)
         ref_data(ix,iy)  = z (offz(ndata - 1) + (ix - 1)*ny(ndata) + iy)
         ref_weight(ix,iy) = 1.0000    ! dz(offxy(iz - 1) + ix) TEMPORARY unit weights
      ENDDO
      ref_y(iy)      = y(offxy(ndata - 1) + iy)
   ENDDO
   DO ix=1,ref_dim(1)
      ref_x(ix)      = x(offxy(ndata - 1) + ix)
   ENDDO
   CALL def_set_variable('real', 'F_XMIN', ref_x(1),          IS_DIFFEV)
   CALL def_set_variable('real', 'F_XMAX', ref_x(ref_dim(1)), IS_DIFFEV)
   CALL def_set_variable('real', 'F_YMIN', ref_y(1),          IS_DIFFEV)
   CALL def_set_variable('real', 'F_YMAX', ref_y(ref_dim(2)), IS_DIFFEV)
   step = (ref_x(ref_dim(1))-ref_x(1))/FLOAT(ref_dim(1)-1)
   CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
   step = (ref_y(ref_dim(2))-ref_y(1))/FLOAT(ref_dim(2)-1)
   CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
   IF(nsigma>0) THEN                    ! User provided a data set with Sigmas
      IF(nsigma>(iz - 1) ) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA Data set number outside KUPLOT range'
         RETURN
      ENDIF
!
      IF(ref_dim(1) /= nx(nsigma) .OR. ref_dim(2) /= ny(nsigma)) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA Data set differs in size'
         RETURN
      ENDIF
      DO iy=1,ref_dim(2)
         DO ix=1,ref_dim(1)
            ref_weight(ix,iy) = z (offz(nsigma - 1) + (ix - 1)*ny(nsigma) + iy)
         ENDDO
      ENDDO
   ENDIF
ELSE
   ref_dim(1) = len(ndata )
   ref_dim(2) = 1
   ALLOCATE(ref_data  (ref_dim(1),ref_dim(2)))
   ALLOCATE(ref_weight(ref_dim(1),ref_dim(2)))
   ALLOCATE(ref_x     (ref_dim(1)))
   ALLOCATE(ref_y     (1         ))
   DO ix=1,ref_dim(1)
      ref_data(ix,1)   = y (offxy(ndata - 1) + ix)
      ref_weight(ix,1) = dy(offxy(ndata - 1) + ix)
      ref_x(ix)        = x (offxy(ndata - 1) + ix)
   ENDDO
   ref_y(1) = 1.0
   CALL def_set_variable('real', 'F_XMIN', ref_x(1),          IS_DIFFEV)
   CALL def_set_variable('real', 'F_XMAX', ref_x(ref_dim(1)), IS_DIFFEV)
   CALL def_set_variable('real', 'F_YMIN', ref_y(1),          IS_DIFFEV)
   CALL def_set_variable('real', 'F_YMAX', ref_y(ref_dim(2)), IS_DIFFEV)
   step = (ref_x(ref_dim(1))-ref_x(1))/FLOAT(ref_dim(1)-1)
   CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
   step = 1.0
   CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
ENDIF
!
!
!
END SUBROUTINE refine_load_kuplot
!
!*******************************************************************************
!
END MODULE refine_load_mod
