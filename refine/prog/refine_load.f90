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
use refine_allocate_appl, only:alloc_weights
!
use kuplot_mod
use kuplot_load_mod
!
use lib_data_struc_h5
USE refine_control_mod
USE refine_data_mod
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
use lib_data_struc_h5
use define_variable_mod
!
IMPLICIT NONE
!
LOGICAL         , INTENT(IN)    :: LDATA   ! Data==TRUE or SIGMA == FALSE
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 20
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE. ! Prevents user from deleting variables
!
character(len=PREC_STRING) :: string
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
!
INTEGER                              :: ianz
integer :: i
INTEGER                              :: idata   ! Current data set
INTEGER                              :: ndata 
REAL(kind=PREC_DP)    :: step
!
integer, parameter :: NOPTIONAL = 3
integer, parameter :: O_DATASET = 1
integer, parameter :: O_WEIGHT  = 1
integer, parameter :: O_MAXDATA = 3
character(len=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(len=max(PREC_STRING,len(line))), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent  !opt. para present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 3 ! Number of values to calculate 
!
DATA oname  / 'dataset', 'weight ', 'maxdata' /
DATA loname /  7       ,  6       ,  7        /
opara  =  (/ '1.0000  ', '1.0000  ', '10.000  ' /)   ! Always provide fresh default values
lopara =  (/  6        ,  6        , 6          /)
owerte =  (/  1.0D0    ,  1.0D0    , 1.0D0      /)
!
!write(*,*) ' REF_NDATA ', ref_ndata
string = line
length = len_trim(string)
!write(*,*) 'LINE  ', line(1:length)
call get_params(line, ianz, cpara, lpara, MAXW, length)
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL, ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
idata = nint(owerte(O_DATASET))
!write(*,*) 'LINE  ', line(1:length)
!write(*,*) 'OPARA ', opara(O_DATASET)(1:lopara(O_DATASET)), ' <>', lopara(O_DATASET)
!write(*,*) ' IER ', ier_num, ier_typ
call purge_optional(line, NOPTIONAL, oname, loname, opara, lopara, lpresent)
!write(*,*) 'LINE  ', line(1:len_trim(line)), len_trim(line), length
!write(*,*) ' IER ', ier_num, ier_typ
if(lpresent(O_DATASET)) then
   if(idata==1) then     ! First of multiple data sets
      if(lpresent(O_MAXDATA)) then               ! User provided maxdata:
         ref_maxdata = nint(owerte(O_MAXDATA))
      else
         ref_maxdata = 10                        ! Use default 10 dData sets
      endif
      if(allocated(ref_data_name))   deallocate(ref_data_name)
      if(allocated(ref_data_number)) deallocate(ref_data_number)
      if(allocated(ref_data_ptr   )) deallocate(ref_data_ptr   )
      if(allocated(ref_calc_ptr   )) deallocate(ref_data_ptr   )
      allocate(ref_data_name(ref_maxdata))
      allocate(ref_data_number(ref_maxdata))
      allocate(ref_data_ptr   (ref_maxdata))
      allocate(ref_calc_ptr   (ref_maxdata))
      ref_data_name = ' '
      ref_data_number = 0
      do i=1, ref_maxdata
         nullify(ref_data_ptr(i)%data_ptr)
      enddo
!QMULTwrite(*,'(2(a, i3))') ' DATA SET ', idata, ' max: ', ref_maxdata
   else
!QMULTwrite(*,'(2(a, i3))') ' DATA SET ', idata, ' max: ', ref_maxdata
   endif
else
!write(*,*) ' NO DATA SET ', ref_ndata
   if(ref_ndata==0) then            ! First data set
      if(lpresent(O_MAXDATA)) then               ! User provided maxdata:
         ref_maxdata = nint(owerte(O_MAXDATA))
      else
         ref_maxdata = 1                         ! Use default 1 data sets single refinement
      endif
      if(allocated(ref_data_name)) deallocate(ref_data_name)
      if(allocated(ref_data_number)) deallocate(ref_data_number)
      if(allocated(ref_data_ptr   )) deallocate(ref_data_ptr   )
      if(allocated(ref_calc_ptr   )) deallocate(ref_data_ptr   )
      allocate(ref_data_name(ref_maxdata))
      allocate(ref_data_number(ref_maxdata))
      allocate(ref_data_ptr   (ref_maxdata))
      allocate(ref_calc_ptr   (ref_maxdata))
      ref_data_name = ' '
      ref_data_number = 0
      do i=1, ref_maxdata
         nullify(ref_data_ptr(i)%data_ptr)
      enddo
   endif
endif
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
      ref_ndata = ref_ndata + 1
   ELSE
      ref_csigma = ' '        ! Sigma loaded from KUPLOT
      ref_ksigma = ndata
      ref_csigma_u = string
   ENDIF
ELSE                               ! Presume a "data xy, filename "
   if(ref_ndata == 0) call dgl5_reset   ! Quietly reset global data storage
   iz = 1                               ! Quiety clear all KUPLOT arrays
   lni = .false.
   lh5 = .false.
   ikfirst = .true.
   ku_ndims = 1
   IF(LDATA) THEN
      ref_load = line
      ref_kload = 0
      ref_load_u = string
      ref_ndata = ref_ndata + 1
   ELSE
      ref_csigma = line
      ref_ksigma = 0
      ref_csigma_u = string
   ENDIF
!write(*,*) ' GOING INTO KUPLOT_LOAD '
   CALL do_load(line, length,.TRUE.)
   IF(ier_num/= 0) RETURN
   ndata = -1                 ! Will be updated to correct value in refine_load_kuplot
ENDIF
!write(*,*) ' Do   LOAD '
!
!********************
!
!CALL refine_load_kuplot(LDATA, ndata, idata)
!
if(ldata) then
   ref_data_number(idata) = dgl5_get_number()          ! Get current node number
!write(*,*) ' CALL dgl5_get_pointer ', ref_data_number(idata), ier_num, ier_typ
   ref_data_ptr(idata)%data_ptr => dgl5_get_pointer()  ! Associate pointer to current node
!write(*,*) ' Associated ptr '
   call dgl5_set_temp(.false.)                   ! Turn temporary status off
   if(.not.ref_data_ptr(idata)%data_ptr%has_dval .or. &
      .not.allocated(ref_data_ptr(idata)%data_ptr%sigma)) then     ! Currently no sigma
      if(allocated(ref_data_ptr(idata)%data_ptr%sigma)) then       ! Sigma exists, check boundary
         if(ubound(   ref_data_ptr(idata)%data_ptr%sigma  ,1) /=     &
            ubound(   ref_data_ptr(idata)%data_ptr%datamap,1)   .or. &
            ubound(   ref_data_ptr(idata)%data_ptr%sigma  ,2) /=     &
            ubound(   ref_data_ptr(idata)%data_ptr%datamap,2)   .or. &
            ubound(   ref_data_ptr(idata)%data_ptr%sigma  ,3) /=     &
            ubound(   ref_data_ptr(idata)%data_ptr%datamap,3) )    then
            ier_num = -14
            ier_num = ER_APPL
            return
         endif
      else
         allocate(ref_data_ptr(idata)%data_ptr%sigma(ref_data_ptr(idata)%data_ptr%dims(1),  &
                                                     ref_data_ptr(idata)%data_ptr%dims(2),  &
                                                     ref_data_ptr(idata)%data_ptr%dims(3)    ))
         ref_data_ptr(idata)%data_ptr%sigma = 1.0_PREC_DP
         ref_data_ptr(idata)%data_ptr%has_dval = .true.
      endif
   endif
endif
!QMULTwrite(*,*) ' Dataset ', idata , 'is at node ', ref_data_number(idata), ref_data_ptr(idata)%data_ptr%infile
!
IF(ndata==-1) ndata = iz-1            ! -1 signals last data set
!
!write(*,*) ' START OF REPLACEMENT '
!write(*,*) ' LDATA, ndata, idata  ', LDATA, ndata, idata
!write(*,*) ' ku_ndims             ', ku_ndims(1)
!write(*,*) ' REF_DATA '
!write(*,*) ' data_ptr', allocated(ref_data_ptr), ubound(ref_data_ptr)
!write(*,*) ' data_ptr', associated(ref_data_ptr(1)%data_ptr)
!write(*,*) ' data_nr ',           (ref_data_ptr(1)%data_ptr%data_num)
!write(*,*) ' infile  ',           (ref_data_ptr(1)%data_ptr%infile(1:30)  )
!write(*,*) ' dims    ',           (ref_data_ptr(1)%data_ptr%dims)
!write(*,*) ' data_x  ', allocated (ref_data_ptr(1)%data_ptr%x)
!write(*,*) ' x(1)    ',           (ref_data_ptr(1)%data_ptr%x(1))
!write(*,*) ' x(n)    ',           (ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1)))
!read(*,*) i
!
!  Set the variables for data boundaries
!
if(idata==1 .and. ldata )  then
   if(ref_maxdata>1) then                              ! First data set  of multiple data sets
   !QMULTwrite(*,*) 'Define multiple variables '
      call def_variable('real', 'M_XMIN', IS_DIFFEV, ref_maxdata)
      call def_variable('real', 'M_XMAX', IS_DIFFEV, ref_maxdata)
      call def_variable('real', 'M_XSTP', IS_DIFFEV, ref_maxdata)
      if(ku_ndims(ndata)==2) then          ! 2D data set
         call def_variable('real', 'M_YMIN', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_YMAX', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_YSTP', IS_DIFFEV, ref_maxdata)
      elseif(ku_ndims(ndata)==3) then          ! 3D data set
         call def_variable('real', 'M_YMIN', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_YMAX', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_YSTP', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_ZMIN', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_ZMAX', IS_DIFFEV, ref_maxdata)
         call def_variable('real', 'M_ZSTP', IS_DIFFEV, ref_maxdata)
      endif
   endif
!  elseif(ref_maxdata==1) then                         ! Single data set refinement
      call def_variable('real', 'F_XMIN', IS_DIFFEV)
      call def_variable('real', 'F_XMAX', IS_DIFFEV)
      call def_variable('real', 'F_XSTP', IS_DIFFEV)
      if(ku_ndims(ndata)==2) then          ! 2D data set
         call def_variable('real', 'F_YMIN', IS_DIFFEV)
         call def_variable('real', 'F_YMAX', IS_DIFFEV)
         call def_variable('real', 'F_YSTP', IS_DIFFEV)
      elseif(ku_ndims(ndata)==3) then          ! 3D data set
         call def_variable('real', 'F_YMIN', IS_DIFFEV)
         call def_variable('real', 'F_YMAX', IS_DIFFEV)
         call def_variable('real', 'F_YSTP', IS_DIFFEV)
         call def_variable('real', 'F_ZMIN', IS_DIFFEV)
         call def_variable('real', 'F_ZMAX', IS_DIFFEV)
         call def_variable('real', 'F_ZSTP', IS_DIFFEV)
      endif
endif
if(ldata .and. ref_maxdata>1) then                     ! Multiple data sets, set values
   call set_variable('M_XMIN', ref_data_ptr(idata)%data_ptr%x(1), idata)
   call set_variable('M_XMAX', ref_data_ptr(idata)%data_ptr%x(ref_data_ptr(idata)%data_ptr%dims(1) ), idata)
   step = (ref_data_ptr(idata)%data_ptr%x(ref_data_ptr(idata)%data_ptr%dims(1))-&
           ref_data_ptr(idata)%data_ptr%x(1))/real(max(1,ref_data_ptr(idata)%data_ptr%dims(1)-1), kind=PREC_DP)
   call set_variable('M_XSTP', step, idata)
   if(ku_ndims(ndata)==2) then          ! 2D data set
      call set_variable('M_YMIN', ref_data_ptr(idata)%data_ptr%y(1), idata)
      call set_variable('M_YMAX', ref_data_ptr(idata)%data_ptr%y(ref_data_ptr(idata)%data_ptr%dims(2) ), idata)
      step = (ref_data_ptr(idata)%data_ptr%y(ref_data_ptr(idata)%data_ptr%dims(2))-&
              ref_data_ptr(idata)%data_ptr%y(1))/real(max(1,ref_data_ptr(idata)%data_ptr%dims(2)-1), kind=PREC_DP)
      call set_variable('M_YSTP', step, idata)
   elseif(ku_ndims(ndata)==3) then          ! 3D data set
      call set_variable('M_YMIN', ref_data_ptr(idata)%data_ptr%y(1), idata)
      call set_variable('M_YMAX', ref_data_ptr(idata)%data_ptr%y(ref_data_ptr(idata)%data_ptr%dims(2) ), idata)
      step = (ref_data_ptr(idata)%data_ptr%y(ref_data_ptr(idata)%data_ptr%dims(2))-&
              ref_data_ptr(idata)%data_ptr%y(1))/real(max(1,ref_data_ptr(idata)%data_ptr%dims(2)-1), kind=PREC_DP)
      call set_variable('M_YSTP', step, idata)
      call set_variable('M_ZMIN', ref_data_ptr(idata)%data_ptr%z(1), idata)
      call set_variable('M_ZMAX', ref_data_ptr(idata)%data_ptr%z(ref_data_ptr(idata)%data_ptr%dims(3) ), idata)
      step = (ref_data_ptr(idata)%data_ptr%z(ref_data_ptr(idata)%data_ptr%dims(3))-&
              ref_data_ptr(idata)%data_ptr%z(1))/real(max(1,ref_data_ptr(idata)%data_ptr%dims(3)-1), kind=PREC_DP)
      call set_variable('M_ZSTP', step, idata)
   endif
endif
!elseif(ldata .and. ref_maxdata==1) then
      CALL set_variable('F_XMIN', ref_data_ptr(1)%data_ptr%x(1))
      CALL set_variable('F_XMAX', ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1)))
      step = (ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1))- &
              ref_data_ptr(1)%data_ptr%x(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(1)-1), kind=PREC_DP)
      CALL set_variable('F_XSTP', step)
   if(ku_ndims(ndata)==2) then          ! 2D data set
      CALL set_variable('F_YMIN', ref_data_ptr(1)%data_ptr%y(1))
      CALL set_variable('F_YMAX', ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2)))
      step = (ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2))- &
              ref_data_ptr(1)%data_ptr%y(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(2)-1), kind=PREC_DP)
      CALL set_variable('F_YSTP', step)
   elseif(ku_ndims(ndata)==3) then          ! 3D data set
      CALL set_variable('F_YMIN', ref_data_ptr(1)%data_ptr%y(1))
      CALL set_variable('F_YMAX', ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2)))
      CALL set_variable('F_ZMIN', ref_data_ptr(1)%data_ptr%z(1))
      CALL set_variable('F_ZMAX', ref_data_ptr(1)%data_ptr%z(ref_data_ptr(1)%data_ptr%dims(3)))
      step = (ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2))- &
              ref_data_ptr(1)%data_ptr%y(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(2)-1), kind=PREC_DP)
      CALL set_variable('F_YSTP', step)
      step = (ref_data_ptr(1)%data_ptr%z(ref_data_ptr(1)%data_ptr%dims(3))- &
              ref_data_ptr(1)%data_ptr%z(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(3)-1), kind=PREC_DP)
      CALL set_variable('F_ZSTP', step)
   endif
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
!********************
!
!write(*,*) ' DONE LOAD '
ref_kupl = MAX(ref_kupl, ndata)   ! This is the last KUPLOT data set that needs to be kept
refine_init = .TRUE.              ! Force initialization, as we have a new data set
if(.not. allocated(ref_weight)) then
   call alloc_weights(ref_ndata)
else
   if(idata>ubound(ref_weight,1)) then
      call alloc_weights(max(idata, ref_ndata))
   endif
endif
ref_weight(idata) = owerte(O_WEIGHT)
!write(*,*) ' INIT LOADED DATA '
!write(*,*) ' CHECK ', ref_data_ptr(1    )%data_ptr%dims
!write(*,*) ' DATA  ', ref_data_ptr(1    )%data_ptr%datamap(1  ,1  ,1  ), &
!                      ref_data_ptr(1    )%data_ptr%datamap(51 ,1  ,1  )
!write(*,*) ' SIGMA ', ref_data_ptr(1    )%data_ptr%sigma(1  ,1  ,1  ), &
!                      ref_data_ptr(1    )%data_ptr%sigma(51 ,1  ,1  )
!read(*,*) i
!
END SUBROUTINE refine_load
!
!*******************************************************************************
!
SUBROUTINE refine_load_kuplot(LDATA, ndata, idata) !, is_data, is_sigma)
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
INTEGER, INTENT(INOUT) :: idata   ! Current Refine data set
!
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE. ! Prevents user from deleting variables
INTEGER :: iix, iiy, iiz                 ! Dummy loop variables
integer :: i1, i2, j1, j2, k1, k2        ! Dummy loop variables
integer :: i, j                          ! Dummy loop variables
REAL(kind=PREC_DP)    :: step
!
!
!write(*,*) ' IN REF_KUPLOT'
!read(*,*) i
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
!   IF(ALLOCATED(ref_x     )) DEALLOCATE(ref_x     )
!   IF(ALLOCATED(ref_y     )) DEALLOCATE(ref_y     )
!   IF(ALLOCATED(ref_z     )) DEALLOCATE(ref_z     )
ENDIF
!
!write(*,*)  ' LINE 237 ', ku_ndims(ndata), ndata
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
         ref_dim(1,1) = 2*(max(abs(nint(ik1_x(1))),abs(nint(ik1_x(ik1_dims(1)))))) + 1
         ref_dim(2,1) = 2*(max(abs(nint(ik1_y(1))),abs(nint(ik1_y(ik1_dims(2)))))) + 1
         ref_dim(3,1) = 2*(max(abs(nint(ik1_z(1))),abs(nint(ik1_z(ik1_dims(3)))))) + 1
         ref_dim(:,1)    = ik1_dims
         ALLOCATE(ref_data  (ref_dim(1,1),ref_dim(2,1),ref_dim(3,1),1))
         ALLOCATE(ref_sigma (ref_dim(1,1),ref_dim(2,1),ref_dim(3,1),1))
!         ALLOCATE(ref_x     (ref_dim(1,1),1))
!         ALLOCATE(ref_y     (ref_dim(2,1),1))
!         ALLOCATE(ref_z     (ref_dim(3,1),1))
         ref_data  =  0.0D0
         ref_sigma = -10000.0D0     ! Flags a missing data point
         i1 =  (ref_dim(1,1)-1)/2 + nint(ik1_x(1)) + 1
         j1 =  (ref_dim(2,1)-1)/2 + nint(ik1_y(1)) + 1
         k1 =  (ref_dim(3,1)-1)/2 + nint(ik1_z(1)) + 1
         i2 = i1 + ik1_dims(1) - 1
         j2 = j1 + ik1_dims(2) - 1
         k2 = k1 + ik1_dims(3) - 1
!        ref_data (i1:i2, j1:j2, k1:k2) = ik1_data
!        ref_sigma(i1:i2, j1:j2, k1:k2) = ik1_sigma
         ref_data (:,:,:,1)             = ik1_data
         ref_sigma(:,:,:,1)             = ik1_sigma
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
!         ref_x(:,1) = ik1_x
!         ref_y(:,1) = ik1_y
!         ref_z(:,1) = ik1_z
!
      else
!
         ref_dim(:,1)    = ik1_dims
         ALLOCATE(ref_data  (ref_dim(1,1),ref_dim(2,1),ref_dim(3,1),1))
         ALLOCATE(ref_sigma (ref_dim(1,1),ref_dim(2,1),ref_dim(3,1),1))
!         ALLOCATE(ref_x     (ref_dim(1,1),1))
!         ALLOCATE(ref_y     (ref_dim(2,1),1))
!         ALLOCATE(ref_z     (ref_dim(3,1),1))
         ref_data(:,:,:,1) = ik1_data
         if(ik1_has_dval) then
            ref_sigma(:,:,:,1) = ik1_sigma
         else
            ref_sigma(:,:,:,1) = 1.0D0
         endif
!         ref_x(:,1) = ik1_x
!         ref_y(:,1) = ik1_y
!         ref_z(:,1) = ik1_z
      endif
      ref_type = ik1_data_type
      CALL def_set_variable('real', 'F_XMIN', ref_data_ptr(1)%data_ptr%x(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_XMAX', ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMIN', ref_data_ptr(1)%data_ptr%y(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMAX', ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_ZMIN', ref_data_ptr(1)%data_ptr%z(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_ZMAX', ref_data_ptr(1)%data_ptr%z(ref_data_ptr(1)%data_ptr%dims(3)), IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1))- &
              ref_data_ptr(1)%data_ptr%x(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(1)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2))- &
              ref_data_ptr(1)%data_ptr%y(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(2)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%z(ref_data_ptr(1)%data_ptr%dims(3))- &
              ref_data_ptr(1)%data_ptr%z(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(3)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_ZSTP', step             , IS_DIFFEV)

   else cond_data3                     ! This is sigma
      IF(.NOT.ALLOCATED(ref_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         exit cond_data3
      ENDIF
!
      IF(ref_dim(1,1) /= ik1_dims(1) .OR. ref_dim(2,1) /= ik1_dims(2) .or. &
         ref_dim(3,1) /= ik1_dims(3) ) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA and DATA set differ in size'
         exit cond_data3
      ENDIF
      ref_sigma(:,:,:,1) = ik1_data    ! copy these data as sigma
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
      ref_dim(1,1) = nx(ndata )
      ref_dim(2,1) = ny(ndata )
      ref_dim(3,1) = 1
      ALLOCATE(ref_data  (ref_dim(1,1),ref_dim(2,1),ref_dim(3,1),1))
      ALLOCATE(ref_sigma (ref_dim(1,1),ref_dim(2,1),ref_dim(3,1),1))
!      ALLOCATE(ref_x     (ref_dim(1,1),1))
!      ALLOCATE(ref_y     (ref_dim(2,1),1))
!      ALLOCATE(ref_z     (ref_dim(3,1),1))
!
      DO iiy=1,ref_dim(2,1)
         DO iix=1,ref_dim(1,1)
            ref_data(iix,iiy, 1,1)  = z (offz(ndata - 1) + (iix - 1)*ny(ndata) + iiy)
            ref_sigma (iix,iiy, 1,1) = 1.0000D0  ! dz(offxy(iz - 1) + iix) TEMPORARY unit weights
         ENDDO
!        ref_y(iiy,1)      = y(offxy(ndata - 1) + iiy)
      ENDDO
!     DO iix=1,ref_dim(1,1)
!        ref_x(iix,1)      = x(offxy(ndata - 1) + iix)
!     ENDDO
!     ref_z = 1.0
      ref_type = H5_2D_GEN

      CALL def_set_variable('real', 'F_XMIN', ref_data_ptr(1)%data_ptr%x(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_XMAX', ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMIN', ref_data_ptr(1)%data_ptr%y(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMAX', ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2)), IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1))- &
              ref_data_ptr(1)%data_ptr%x(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(1)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2))- &
              ref_data_ptr(1)%data_ptr%y(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(2)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
   ELSE
      IF(.NOT.ALLOCATED(ref_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         RETURN
      ENDIF
!
      IF(ref_dim(1,1) /= nx(ndata) .OR. ref_dim(2,1) /= ny(ndata)) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA and DATA set differ in size'
         RETURN
      ENDIF
      DO iiy=1,ref_dim(2,1)
         DO iix=1,ref_dim(1,1)
            ref_sigma (iix,iiy,iiz,1) = z (offz(ndata - 1) + (iix - 1)*ny(ndata) + iiy)
         ENDDO
      ENDDO
   ENDIF
ELSEif(ku_ndims(ndata)==1) then          ! 1D data set
!write(*,*) 'REFINE_LOAD_KUPLOT ku_ndims ', ku_ndims(ndata)
   IF(LDATA) THEN                      ! This is the data set
      ref_dim(1,1) = lenc(ndata )
      ref_dim(2,1) = 1
      ref_dim(3,1) = 1
      ref_type = H5_1D_GEN
      ALLOCATE(ref_data  (ref_dim(1,1),ref_dim(2,1), ref_dim(3,1),1))
      ALLOCATE(ref_sigma (ref_dim(1,1),ref_dim(2,1), ref_dim(3,1),1))
!      ALLOCATE(ref_x     (ref_dim(1,1),1))
!      ALLOCATE(ref_y     (1           ,1))
!      ALLOCATE(ref_z     (1           ,1))
      DO iix=1,ref_dim(1,1)
         ref_data(iix,1, 1,1)   = y(offxy(ndata - 1) + iix)
         ref_sigma (iix,1, 1,1) = ABS(dy(offxy(ndata - 1) + iix))
!         ref_x(iix,1)        = x(offxy(ndata - 1) + iix)
      ENDDO
      IF(MINVAL(ref_sigma (:,1,1,1))==0.0) THEN
         ier_num = -7
         ier_typ = ER_APPL
         ier_msg(1) = ' Check data and define non-zero sigma'
         RETURN
      ENDIF
!write(*,*) ' SIGMA ', ubound(ref_sigma)
!write(*,*) ' data_ptr', allocated(ref_data_ptr)
!write(*,*) ' data_ptr', associated(ref_data_ptr(1)%data_ptr)
!write(*,*) ' data_nr ',           (ref_data_ptr(1)%data_ptr%data_num)
!write(*,*) ' data_x  ', allocated (ref_data_ptr(1)%data_ptr%x)

      CALL def_set_variable('real', 'F_XMIN', ref_data_ptr(1)%data_ptr%x(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_XMAX', ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1)), IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMIN', ref_data_ptr(1)%data_ptr%y(1),          IS_DIFFEV)
      CALL def_set_variable('real', 'F_YMAX', ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2)), IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%x(ref_data_ptr(1)%data_ptr%dims(1))- &
              ref_data_ptr(1)%data_ptr%x(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(1)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
      step = (ref_data_ptr(1)%data_ptr%y(ref_data_ptr(1)%data_ptr%dims(2))- &
              ref_data_ptr(1)%data_ptr%y(1))/ &
              real(max(1,ref_data_ptr(1)%data_ptr%dims(2)-1), kind=PREC_DP)
      CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
   ELSE
      IF(.NOT.ALLOCATED(ref_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         RETURN
      ENDIF
!
      IF(ref_dim(1,1) /= nx(ndata)) THEN
         ier_num = -6
         ier_typ = ER_FORT
         ier_msg(1) = 'SIGMA data set differs in size'
         RETURN
      ENDIF
!
      DO iix=1,ref_dim(1,1)
         ref_sigma (iix,1,1,1) = y (offxy(ndata - 1) + iix)
      ENDDO
   ENDIF
ENDIF
!
cond_h5: if(.not.lh5(ndata)) then    ! Data set is not in data_struc   ! WORK)
!  Add to data_struc
if(ldata)  ref_data_name(idata) = ik1_infile
!
if(ku_ndims(ndata)==3) then
   if(ldata) then                                   ! Data set
!   write(*,*) ' Node Number ', ik1_node_number
      call dgl5_set_temp(.false.)                   ! Turn temporary status off
      return
   else
      continue                                      ! WORK add sigma to original data set
      return
   endif
else
!
call dgl5_new_node
ik1_node_number = dgl5_get_number()
!write(*,*) ' added node ', ik1_node_number
!
! Copy data into data_struc
!
if(ku_ndims(ndata)==2) then          ! 2D data set
   ik1_infile = fname(ndata)
   ref_dim(1,1) = nx(ndata)
   ref_dim(2,1) = ny(ndata)
   ref_dim(3,1) = 1
   ik1_dims(1) = nx(ndata)
   ik1_dims(2) = ny(ndata)
   ik1_dims(3) = 1
   if(allocated(ik1_x)) deallocate(ik1_x)
   if(allocated(ik1_y)) deallocate(ik1_y)
   if(allocated(ik1_z)) deallocate(ik1_z)
   if(allocated(ik1_dx)) deallocate(ik1_dx)
   if(allocated(ik1_dy)) deallocate(ik1_dy)
   if(allocated(ik1_dz)) deallocate(ik1_dz)
   if(allocated(ik1_data))  deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   allocate(ik1_x(1:ref_dim(1,1)))
   allocate(ik1_y(1:ref_dim(2,1)))
   allocate(ik1_z(1:1))
   allocate(ik1_dx(1:1))
   allocate(ik1_dy(1:1))
   allocate(ik1_dz(1:1))
   allocate(ik1_data(1:ref_dim(1,1), 1:ref_dim(2,1), 1:ref_dim(3,1)))
   allocate(ik1_sigma(1,1,1))
   ik1_data_type = H5_2D_GEN
   ik1_nlayer    = 1
   ik1_is_direct = .true.
   ik1_is_grid   = .true.
   ik1_has_dxyz  = .false.
   ik1_has_dval  = .false.
   ik1_calc_coor = .false.
   ik1_use_coor  = 1
   ik1_corners   = 0.0_PREC_DP
   ik1_vectors   = 0.0_PREC_DP
   ik1_x   (1:ref_dim(1,1))     = x(offxy(ndata-1)+1:offxy(ndata-1)+lenc(ndata))
   ik1_y   (1:ref_dim(1,1))     = y(offxy(ndata-1)+1:offxy(ndata-1)+lenc(ndata))
   do i=1, nx(ndata)
      do j=1,ny(ndata)
         ik1_data(i,j,1) = z(offz(ndata-1) + (i - 1) * ny (ndata) + j)
      enddo
   enddo
   ik1_sigma     = 0.0_PREC_DP
   ik1_llims     = 0.0_PREC_DP
   ik1_steps     = 0.0_PREC_DP
   ik1_steps_full= 0.0_PREC_DP
!
   call dgl5_set_node(ik1_infile, ik1_data_type, ik1_nlayer, ik1_is_direct, ku_ndims(ndata), ref_dim ,         &
                   ik1_is_grid, ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, &
                   ik1_corners, ik1_vectors,&
                   ik1_a0(1:3), ik1_win(1:3), ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy,  &
                   ik1_dz,     ik1_data               ,ik1_sigma, ik1_llims,      &
                   ik1_steps, ik1_steps_full)
   call dgl5_set_temp(.false.)                    ! Turn temporary status off
elseif(ku_ndims(ndata)==1) then          ! 1D data set
!
!  Define multiple coordinate variable M_XMIN etc
!
   ik1_infile = fname(ndata)
   ref_dim(1,1) = lenc(ndata)
   ref_dim(2,1) = 1
   ref_dim(3,1) = 1
   ik1_dims(1) = lenc(ndata)
   ik1_dims(2) = 1
   ik1_dims(3) = 1
   if(allocated(ik1_x)) deallocate(ik1_x)
   if(allocated(ik1_y)) deallocate(ik1_y)
   if(allocated(ik1_z)) deallocate(ik1_z)
   if(allocated(ik1_dx)) deallocate(ik1_dx)
   if(allocated(ik1_dy)) deallocate(ik1_dy)
   if(allocated(ik1_dz)) deallocate(ik1_dz)
   if(allocated(ik1_data))  deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   allocate(ik1_x(1:ref_dim(1,1)))
   allocate(ik1_y(1:1))
   allocate(ik1_z(1:1))
   allocate(ik1_dx(1:1))
   allocate(ik1_dy(1:1))
   allocate(ik1_dz(1:1))
   allocate(ik1_data(1:ref_dim(1,1), 1:ref_dim(2,1), 1:ref_dim(3,1)))
   allocate(ik1_sigma(1:ref_dim(1,1), 1:ref_dim(2,1), 1:ref_dim(3,1)))
   ik1_data_type = H5_1D_GEN
   ik1_nlayer    = 1
   ik1_is_direct = .true.
   ik1_is_grid   = .true.
   ik1_has_dxyz  = .false.
   ik1_has_dval  = .true.
   ik1_calc_coor = .false.
   ik1_use_coor  = 1
   ik1_corners   = 0.0_PREC_DP
   ik1_vectors   = 0.0_PREC_DP
   ik1_x   (1:ref_dim(1,1))     = x(offxy(ndata-1)+1:offxy(ndata-1)+lenc(ndata))
   ik1_y         = 0.0_PREC_DP
   ik1_z         = 0.0_PREC_DP
   ik1_data(1:ref_dim(1,1),1,1) = y(offxy(ndata-1)+1:offxy(ndata-1)+lenc(ndata))
   ik1_sigma(1:ref_dim(1,1),1,1) = dy(offxy(ndata-1)+1:offxy(ndata-1)+lenc(ndata))
!  ik1_sigma     = 0.0_PREC_DP
   ik1_llims     = 0.0_PREC_DP
   ik1_steps     = 0.0_PREC_DP
   ik1_steps_full= 0.0_PREC_DP
!
   call dgl5_set_node(ik1_infile, ik1_data_type, ik1_nlayer, ik1_is_direct, ku_ndims(ndata), ref_dim ,         &
                   ik1_is_grid, ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, &
                   ik1_corners, ik1_vectors,&
                   ik1_a0(1:3), ik1_win(1:3), ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy,  &
                   ik1_dz,     ik1_data               ,ik1_sigma, ik1_llims,      &
                   ik1_steps, ik1_steps_full)
   call dgl5_set_temp(.false.)                    ! Turn temporary status off
endif
endif
endif cond_h5
!
      call dgl5_set_temp(.false.)                   ! Turn temporary status off
!
if(ldata) then 
   ref_data_number(idata) = dgl5_get_number()
!QMULTwrite(*,*) ' CALL dgl5_get_pointer '
   ref_data_ptr(idata)%data_ptr => dgl5_get_pointer()
endif
!QMULTwrite(*,*) ' Dataset ', idata , 'is at node ', ref_data_number(idata), ref_data_ptr(idata)%data_ptr%infile
!
if(idata==1 .and. ldata .and. ref_maxdata>1) then   ! First data set
   !QMULTwrite(*,*) 'Define multiple variables '
   call def_variable('real', 'M_XMIN', .true., ref_maxdata)
   call def_variable('real', 'M_XMAX', .true., ref_maxdata)
   call def_variable('real', 'M_XSTP', .true., ref_maxdata)
   if(ku_ndims(ndata)==2) then          ! 2D data set
      call def_variable('real', 'M_YMIN', .true., ref_maxdata)
      call def_variable('real', 'M_YMAX', .true., ref_maxdata)
      call def_variable('real', 'M_YSTP', .true., ref_maxdata)
   elseif(ku_ndims(ndata)==3) then          ! 3D data set
      call def_variable('real', 'M_YMIN', .true., ref_maxdata)
      call def_variable('real', 'M_YMAX', .true., ref_maxdata)
      call def_variable('real', 'M_YSTP', .true., ref_maxdata)
      call def_variable('real', 'M_ZMIN', .true., ref_maxdata)
      call def_variable('real', 'M_ZMAX', .true., ref_maxdata)
      call def_variable('real', 'M_ZSTP', .true., ref_maxdata)
   endif
endif
if(ldata .and. ref_maxdata>1) then
   call set_variable('M_XMIN', ik1_x(1), idata)
   call set_variable('M_XMAX', ik1_x(ik1_dims(1) ), idata)
   step = (ik1_x(ik1_dims(1))-ik1_x(1))/FLOAT(ik1_dims(1)-1)
   call set_variable('M_XSTP', step, idata)
   if(ku_ndims(ndata)==2) then          ! 2D data set
      call set_variable('M_YMIN', ik1_y(1), idata)
      call set_variable('M_YMAX', ik1_y(ik1_dims(2) ), idata)
      step = (ik1_y(ik1_dims(2))-ik1_y(1))/FLOAT(ik1_dims(2)-1)
      call set_variable('M_YSTP', step, idata)
   elseif(ku_ndims(ndata)==3) then          ! 3D data set
      call set_variable('M_YMIN', ik1_y(1), idata)
      call set_variable('M_YMAX', ik1_y(ik1_dims(2) ), idata)
      step = (ik1_y(ik1_dims(2))-ik1_y(1))/FLOAT(ik1_dims(2)-1)
      call set_variable('M_YSTP', step, idata)
      call set_variable('M_ZMIN', ik1_z(1), idata)
      call set_variable('M_ZMAX', ik1_z(ik1_dims(3) ), idata)
      step = (ik1_z(ik1_dims(3))-ik1_z(1))/FLOAT(ik1_dims(3)-1)
      call set_variable('M_ZSTP', step, idata)
   endif
endif
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
END SUBROUTINE refine_load_kuplot
!
!*******************************************************************************
!
END MODULE refine_load_mod
