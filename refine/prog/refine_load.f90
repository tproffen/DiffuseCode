MODULE refine_load_mod
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_load(LDATA, line, length)
!
! Loads the Data set and/or the Sigmas
! Either explicitly or with reference to a KUPLOT data set
!
USE refine_control_mod
USE refine_data_mod
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!
LOGICAL         , INTENT(IN)    :: LDATA   ! Datai==TRUE ot SIGMA == FALSE
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
INTEGER                              :: ndata 
!
LOGICAL, EXTERNAL :: str_comp
!
IF(line(1:6) == 'kuplot') THEN
   CALL get_params(line, ianz, cpara, lpara, MAXW, length)
   IF(ier_num/= 0) RETURN
   cpara(1) = '0'
   lpara(1) = 1
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/= 0) RETURN
   ndata = NINT(werte(2))
   IF(LDATA) THEN             ! Data loaded from KUPLOT
      ref_load = ' '
      ref_kload = ndata
   ELSE
      ref_csigma = ' '        ! Sigma loaded from KUPLOT
      ref_ksigma = ndata
   ENDIF
ELSE                               ! Presume a "data xy, filename "
   IF(LDATA) THEN
      ref_load = line
      ref_kload = 0
   ELSE
      ref_csigma = line
      ref_ksigma = 0
   ENDIF
   CALL do_load(line, length,.TRUE.)
   IF(ier_num/= 0) RETURN
   ndata = -1                 ! Will be updated to correct value in refine_load_kuplot
ENDIF
CALL refine_load_kuplot(LDATA, ndata)
ref_kupl = MAX(ref_kupl, ndata)   ! This is the last KUPLOT data set that needs to be kept
refine_init = .TRUE.              ! Force initialization, as we have a new data set
!
END SUBROUTINE refine_load
!
!*******************************************************************************
!
SUBROUTINE refine_load_kuplot(LDATA, ndata) !, is_data, is_sigma)
!
! Transfers data set no. ndata from KUPLOT into REFINE memory
! If LDATA== TRUE, it is the data, else it is the sigma
! If ndata==-1, the last KUPLOT data set is taken
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
LOGICAL, INTENT(IN) :: LDATA   ! Datai==TRUE ot SIGMA == FALSE
INTEGER, INTENT(INOUT) :: ndata   ! no of data set to be transfered from KUPLOT
!
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE. ! Prevents user from deleting variables
INTEGER :: ix, iy                      ! Dummy loop variables
REAL    :: step
!
IF(ndata==-1) ndata = iz-1            ! -1 signals last data set
!
IF(ndata<1 .OR. ndata>(iz - 1) ) THEN
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Data set number outside KUPLOT range'
   RETURN
ENDIF
!
IF(LDATA) THEN                         ! This is the data set
   IF(ALLOCATED(ref_data))   DEALLOCATE(ref_data)
   IF(ALLOCATED(ref_sigma )) DEALLOCATE(ref_sigma )
   IF(ALLOCATED(ref_x     )) DEALLOCATE(ref_x     )
   IF(ALLOCATED(ref_y     )) DEALLOCATE(ref_y     )
ENDIF
!
IF(lni(ndata)) THEN                    ! 2D data set
   IF(LDATA) THEN                      ! This is the data set
      ref_dim(1) = nx(ndata )
      ref_dim(2) = ny(ndata )
      ALLOCATE(ref_data  (ref_dim(1),ref_dim(2)))
      ALLOCATE(ref_sigma (ref_dim(1),ref_dim(2)))
      ALLOCATE(ref_x     (ref_dim(1)))
      ALLOCATE(ref_y     (ref_dim(2)))
!
      DO iy=1,ref_dim(2)
         DO ix=1,ref_dim(1)
            ref_data(ix,iy)  = z (offz(ndata - 1) + (ix - 1)*ny(ndata) + iy)
            ref_sigma (ix,iy) = 1.0000    ! dz(offxy(iz - 1) + ix) TEMPORARY unit weights
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
   ELSE
      IF(.NOT.ALLOCATED(ref_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         RETURN
      ENDIF
!
      IF(ref_dim(1) /= nx(ndata) .OR. ref_dim(2) /= ny(ndata)) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA data set differs in size'
         RETURN
      ENDIF
      DO iy=1,ref_dim(2)
         DO ix=1,ref_dim(1)
            ref_sigma (ix,iy) = z (offz(ndata - 1) + (ix - 1)*ny(ndata) + iy)
         ENDDO
      ENDDO
   ENDIF
ELSE                                     ! !D data set
   IF(LDATA) THEN                      ! This is the data set
      ref_dim(1) = len(ndata )
      ref_dim(2) = 1
      ALLOCATE(ref_data  (ref_dim(1),ref_dim(2)))
      ALLOCATE(ref_sigma (ref_dim(1),ref_dim(2)))
      ALLOCATE(ref_x     (ref_dim(1)))
      ALLOCATE(ref_y     (1         ))
      DO ix=1,ref_dim(1)
         ref_data(ix,1)   = y(offxy(ndata - 1) + ix)
         ref_sigma (ix,1) = ABS(dy(offxy(ndata - 1) + ix))
         ref_x(ix)        = x(offxy(ndata - 1) + ix)
      ENDDO
      IF(MINVAL(ref_sigma (:,1))==0.0) THEN
         ier_num = -7
         ier_typ = ER_APPL
         ier_msg(1) = ' Check data and define non-zeo sigma'
         RETURN
      ENDIF
      ref_y(1) = 1.0
      CALL def_set_variable('real', 'F_XMIN', ref_x(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_XMAX', ref_x(ref_dim(1)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMIN', ref_y(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMAX', ref_y(ref_dim(2)), IS_DIFFEV)
      step = (ref_x(ref_dim(1))-ref_x(1))/FLOAT(ref_dim(1)-1)
      CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
      step = 1.0
      CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
   ELSE
      IF(.NOT.ALLOCATED(ref_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         RETURN
      ENDIF
!
      IF(ref_dim(1) /= nx(ndata)) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA data set differs in size'
         RETURN
      ENDIF
!
      DO ix=1,ref_dim(1)
         ref_sigma (ix,1) = y (offxy(ndata - 1) + ix)
      ENDDO
   ENDIF
ENDIF
!
! Scratch all KUPLOT data beyond current data set
!
IF(LDATA) THEN
   CALL def_set_variable('integer', 'F_DATA', FLOAT(ndata),  IS_DIFFEV)
   CALL def_set_variable('integer', 'F_SIGMA', FLOAT(ndata), IS_DIFFEV)
ELSE
   CALL def_set_variable('integer', 'F_SIGMA', FLOAT(ndata), IS_DIFFEV)
ENDIF
iz = ndata + 1
!
!
END SUBROUTINE refine_load_kuplot
!
!*******************************************************************************
!
END MODULE refine_load_mod
