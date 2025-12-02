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
use kuplot_mod
use kuplot_load_mod
!
USE refine_control_mod
USE refine_data_mod
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
use lib_data_struc_h5
!
IMPLICIT NONE
!
LOGICAL         , INTENT(IN)    :: LDATA   ! Data==TRUE or SIGMA == FALSE
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
!
character(len=PREC_STRING) :: string
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
INTEGER                              :: ndata 
!
string = line
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
      ref_load_u = string
   ELSE
      ref_csigma = ' '        ! Sigma loaded from KUPLOT
      ref_ksigma = ndata
      ref_csigma_u = string
   ENDIF
ELSE                               ! Presume a "data xy, filename "
   call dgl5_reset            ! Queitly reset global data storage
   iz = 1                     ! Quiety clear all KUPLOT arrays
   lni = .false.
   lh5 = .false.
   ikfirst = .true.
   ku_ndims = 1
   IF(LDATA) THEN
      ref_load = line
      ref_kload = 0
      ref_load_u = string
   ELSE
      ref_csigma = line
      ref_ksigma = 0
      ref_csigma_u = string
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
use lib_data_struc_h5
use lib_data_types_mod
use lib_ik_mod
use precision_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN)    :: LDATA   ! Datai==TRUE ot SIGMA == FALSE
INTEGER, INTENT(INOUT) :: ndata   ! no of data set to be transfered from KUPLOT
!
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE. ! Prevents user from deleting variables
INTEGER :: iix, iiy, iiz                 ! Dummy loop variables
integer :: i1, i2, j1, j2, k1, k2        ! Dummy loop variables
REAL(kind=PREC_DP)    :: step
!
iiz = 1
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
   IF(ALLOCATED(ref_z     )) DEALLOCATE(ref_z     )
ENDIF
!
!IF(lni(ndata)) THEN                    ! 2D data set
IF(ku_ndims(ndata)==3) THEN             ! 3D data set
!
      call data2local(ndata   , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
           ik1_data_type,    &
           ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
           ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, ik1_corners,   &
           ik1_vectors, ik1_a0, ik1_win,  &
           ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
           ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
   cond_data3: if(LDATA) then                      ! This is the data set
!
!write(*,*) ' LOADED DATA ', ik1_data_type==H5_BRAGG_I
      if(ik1_data_type==H5_BRAGG_I) then
!
         ref_dim(1) = 2*(max(abs(nint(ik1_x(1))),abs(nint(ik1_x(ik1_dims(1)))))) + 1
         ref_dim(2) = 2*(max(abs(nint(ik1_y(1))),abs(nint(ik1_y(ik1_dims(2)))))) + 1
         ref_dim(3) = 2*(max(abs(nint(ik1_z(1))),abs(nint(ik1_z(ik1_dims(3)))))) + 1
         ref_dim    = ik1_dims
         ALLOCATE(ref_data  (ref_dim(1),ref_dim(2),ref_dim(3)))
         ALLOCATE(ref_sigma (ref_dim(1),ref_dim(2),ref_dim(3)))
         ALLOCATE(ref_x     (ref_dim(1)))
         ALLOCATE(ref_y     (ref_dim(2)))
         ALLOCATE(ref_z     (ref_dim(3)))
         ref_data  =  0.0D0
         ref_sigma = -10000.0D0     ! Flags a missing data point
         i1 =  (ref_dim(1)-1)/2 + nint(ik1_x(1)) + 1
         j1 =  (ref_dim(2)-1)/2 + nint(ik1_y(1)) + 1
         k1 =  (ref_dim(3)-1)/2 + nint(ik1_z(1)) + 1
         i2 = i1 + ik1_dims(1) - 1
         j2 = j1 + ik1_dims(2) - 1
         k2 = k1 + ik1_dims(3) - 1
!        ref_data (i1:i2, j1:j2, k1:k2) = ik1_data
!        ref_sigma(i1:i2, j1:j2, k1:k2) = ik1_sigma
         ref_data                       = ik1_data
         ref_sigma                      = ik1_sigma
         where(ref_sigma<=0.0D0 .and. ref_sigma>-1000.0D0) 
            ref_sigma = 9D19
         end where
!        do iix=1,ref_dim(1)
!           ref_x(iix) = -(ref_dim(1)-1)/2 + iix - 1
!        enddo
!        do iix=1,ref_dim(2)
!           ref_y(iix) = -(ref_dim(2)-1)/2 + iix - 1
!        enddo
!        do iix=1,ref_dim(3)
!           ref_z(iix) = -(ref_dim(3)-1)/2 + iix - 1
!        enddo
         ref_x = ik1_x
         ref_y = ik1_y
         ref_z = ik1_z
!
      else
!
         ref_dim    = ik1_dims
         ALLOCATE(ref_data  (ref_dim(1),ref_dim(2),ref_dim(3)))
         ALLOCATE(ref_sigma (ref_dim(1),ref_dim(2),ref_dim(3)))
         ALLOCATE(ref_x     (ref_dim(1)))
         ALLOCATE(ref_y     (ref_dim(2)))
         ALLOCATE(ref_z     (ref_dim(3)))
         ref_data = ik1_data
         if(ik1_has_dval) then
            ref_sigma = ik1_sigma
         else
            ref_sigma = 1.0D0
         endif
         ref_x = ik1_x
         ref_y = ik1_y
         ref_z = ik1_z
      endif
      ref_type = ik1_data_type
      CALL def_set_variable('real', 'F_XMIN', ref_x(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_XMAX', ref_x(ref_dim(1)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMIN', ref_y(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMAX', ref_y(ref_dim(2)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_ZMIN', ref_z(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_ZMAX', ref_z(ref_dim(3)), IS_DIFFEV)
      step = (ref_x(ref_dim(1))-ref_x(1))/FLOAT(ref_dim(1)-1)
      CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
      step = (ref_y(ref_dim(2))-ref_y(1))/FLOAT(ref_dim(2)-1)
      CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
      step = (ref_z(ref_dim(3))-ref_z(1))/FLOAT(ref_dim(3)-1)
      CALL def_set_variable('real', 'F_ZSTP', step             , IS_DIFFEV)

   else cond_data3                     ! This is sigma
      IF(.NOT.ALLOCATED(ref_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         exit cond_data3
      ENDIF
!
      IF(ref_dim(1) /= ik1_dims(1) .OR. ref_dim(2) /= ik1_dims(2) .or. &
         ref_dim(3) /= ik1_dims(3) ) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA and DATA set differ in size'
         exit cond_data3
      ENDIF
      ref_sigma = ik1_data    ! copy these data as sigma
!
      deallocate(ik1_data)
      if(allocated(ik1_sigma)) deallocate(ik1_sigma)
      if(allocated(ik1_x))     deallocate(ik1_x)
      if(allocated(ik1_y))     deallocate(ik1_y)
      if(allocated(ik1_z))     deallocate(ik1_z)
      if(ier_num/=0) return
!
   endif cond_data3
!write(*,*) ' REF_DIM ', ref_dim
!write(*,*) ' bound x ', lbound(ref_x), ubound(ref_x)
!write(*,*) ' bound y ', lbound(ref_y), ubound(ref_y)
!write(*,*) ' bound z ', lbound(ref_z), ubound(ref_z)
!write(*,*) ' X       ', ref_x(lbound(ref_x)), ref_x(ubound(ref_x))
!write(*,*) ' Y       ', ref_y(lbound(ref_y)), ref_y(ubound(ref_y))
!write(*,*) ' Z       ', ref_z(lbound(ref_z)), ref_z(ubound(ref_z))
!
elseif(ku_ndims(ndata)==2) THEN         ! 2D data set
   IF(LDATA) THEN                      ! This is the data set
      ref_dim(1) = nx(ndata )
      ref_dim(2) = ny(ndata )
      ref_dim(3) = 1
      ALLOCATE(ref_data  (ref_dim(1),ref_dim(2),ref_dim(3)))
      ALLOCATE(ref_sigma (ref_dim(1),ref_dim(2),ref_dim(3)))
      ALLOCATE(ref_x     (ref_dim(1)))
      ALLOCATE(ref_y     (ref_dim(2)))
      ALLOCATE(ref_z     (ref_dim(3)))
!
      DO iiy=1,ref_dim(2)
         DO iix=1,ref_dim(1)
            ref_data(iix,iiy, 1)  = z (offz(ndata - 1) + (iix - 1)*ny(ndata) + iiy)
            ref_sigma (iix,iiy, 1) = 1.0000D0  ! dz(offxy(iz - 1) + iix) TEMPORARY unit weights
         ENDDO
         ref_y(iiy)      = y(offxy(ndata - 1) + iiy)
      ENDDO
      DO iix=1,ref_dim(1)
         ref_x(iix)      = x(offxy(ndata - 1) + iix)
      ENDDO
      ref_z = 1.0
      ref_type = H5_2D_GEN
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
         ier_msg(1) = 'SIGMA and DATA set differ in size'
         RETURN
      ENDIF
      DO iiy=1,ref_dim(2)
         DO iix=1,ref_dim(1)
            ref_sigma (iix,iiy,iiz) = z (offz(ndata - 1) + (iix - 1)*ny(ndata) + iiy)
         ENDDO
      ENDDO
   ENDIF
ELSEif(ku_ndims(ndata)==1) then          ! 1D data set
   IF(LDATA) THEN                      ! This is the data set
      ref_dim(1) = lenc(ndata )
      ref_dim(2) = 1
      ref_dim(3) = 1
      ref_type = H5_1D_GEN
      ALLOCATE(ref_data  (ref_dim(1),ref_dim(2), ref_dim(3)))
      ALLOCATE(ref_sigma (ref_dim(1),ref_dim(2), ref_dim(3)))
      ALLOCATE(ref_x     (ref_dim(1)))
      ALLOCATE(ref_y     (1         ))
      ALLOCATE(ref_z     (1         ))
      DO iix=1,ref_dim(1)
         ref_data(iix,1, 1)   = y(offxy(ndata - 1) + iix)
         ref_sigma (iix,1, 1) = ABS(dy(offxy(ndata - 1) + iix))
         ref_x(iix)        = x(offxy(ndata - 1) + iix)
      ENDDO
      IF(MINVAL(ref_sigma (:,1,1))==0.0) THEN
         ier_num = -7
         ier_typ = ER_APPL
         ier_msg(1) = ' Check data and define non-zero sigma'
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
      DO iix=1,ref_dim(1)
         ref_sigma (iix,1,1) = y (offxy(ndata - 1) + iix)
      ENDDO
   ENDIF
ENDIF
!
! Scratch all KUPLOT data beyond current data set
!
IF(LDATA) THEN
   CALL def_set_variable('integer', 'F_DATA', real(ndata,kind=PREC_DP),  IS_DIFFEV)
   CALL def_set_variable('integer', 'F_SIGMA', real(ndata, kind=PREC_DP), IS_DIFFEV)
ELSE
   CALL def_set_variable('integer', 'F_SIGMA', real(ndata, kind=PREC_DP), IS_DIFFEV)
ENDIF
iz = ndata + 1
!
!
END SUBROUTINE refine_load_kuplot
!
!*******************************************************************************
!
END MODULE refine_load_mod
