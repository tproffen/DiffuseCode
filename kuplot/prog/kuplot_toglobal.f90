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
use refine_params_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
use lib_ik_mod
use lib_data_struc_h5
use lib_data_types_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*) , INTENT(INOUT) :: line         ! Command line
!
INTEGER, PARAMETER :: MAXW = 3
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara    ! Parameter strings
INTEGER            , DIMENSION(MAXW) :: lpara    ! length of each parameter strign
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte    ! Parameter values
!
integer :: i,j      ! Dummy loop indices
INTEGER :: ik      ! Data set to be transfered to ig
INTEGER :: ig      ! Number in global data set
integer :: idata   ! Data set number
INTEGER :: length  ! Number in global data set
INTEGER :: ianz    ! Number in global data set
INTEGER :: iianz   ! Number in global data set
INTEGER :: nnpara   ! Number refined parameters used by refine
INTEGER :: nnfix    ! Number fixed   parameters used by refine
INTEGER, DIMENSION(3) :: dimen
INTEGER, DIMENSION(4) :: dimen_alloc
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: ext_data
!
integer :: node_number      ! Node in data_struc
!
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_DATASET = 1
INTEGER, PARAMETER :: O_REFINE  = 2
INTEGER, PARAMETER :: O_KUPL    = 3
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate
!
DATA oname  / 'dataset', 'refine'  , 'kuplot' /
DATA loname /  7       ,  6        ,  6       /
opara  =  (/ '1       ', 'cost    ', 'last    '  /)
lopara =  (/ 1         , 4         ,  4          /)
owerte =  (/ 1.000000  , -1.000000 ,  -1.0000000 /)
!
!write(*,*) 'IN TO REFINE'
CALL gl_get_npara(NNPARA, NNFIX)      ! Check number of refined parameters
!IF(NNPARA==0) RETURN                  ! if zero, silently leave
!
length = LEN_TRIM(line)
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!write(*,*) 'IN TO REFINE', ier_num, ier_typ
IF(ier_num/=0) RETURN
!
idata = nint(owerte(O_DATASET))
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
if_parname: IF(lpresent(O_REFINE)) THEN
   IF(opara(O_REFINE)=='cost') THEN     ! Main cost function
      ig = 0
   elseif(opara(O_REFINE)=='obs' .or. opara(O_REFINE)=='exp') THEN   ! Experimental data
      ig = -2
   elseif(opara(O_REFINE)=='sigma') then  ! Uncertainties
      ig = -3
   elseif(opara(O_REFINE)=='opti') then   ! Derivative for a parameter
      ig = -1
   ELSE
      do ig=1, refine_par_n
         if(opara(O_REFINE)==refine_params(ig)) exit if_parname
      enddo
      cpara(1) = opara(O_REFINE)
      lpara(1) = lopara(O_REFINE)
      iianz = 1
      CALL ber_params(iianz, cpara, lpara, werte, MAXW)
      IF(ier_num/=0) RETURN
      ig = NINT(werte(1))
   ENDIF
ENDIF if_parname
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
!IF(lni(ik)) THEN            ! Data set is 2D (NIPL)
if(NNPARA>0) then
if(ku_ndims(ik)==3) then         ! Data set is 3D
   call data2local(ik      , ier_num, ier_typ, ik1_node_number, ik1_infile,     &
     ik1_data_type, &
     ik1_nlayer, ik1_is_direct, ik1_ndims, ik1_dims, ik1_is_grid,            &
     ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, ik1_corners,   &
     ik1_vectors, ik1_a0, ik1_win,  &
     ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy, ik1_dz, ik1_data, ik1_sigma,       &
     ik1_llims, ik1_steps,  ik1_steps_full, ik1_minmaxval, ik1_minmaxcoor)
!
   dimen = ik1_dims
   dimen_alloc(1:3) = dimen
   dimen_alloc(4)   = refine_par_n
   if(.not. gl_is_allocated()) call gl_alloc(dimen_alloc)
!
   CALL gl_set_data(dimen(1), dimen(2), dimen(3), ig, ik1_data )
   deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   if(allocated(ik1_x    )) deallocate(ik1_x    )
   if(allocated(ik1_y    )) deallocate(ik1_y    )
   if(allocated(ik1_z    )) deallocate(ik1_z    )
elseif(ku_ndims(ik)==2) then     ! Data set is 2D (NIPL)
   dimen(1) = nx(ik)
   dimen(2) = ny(ik)
   dimen(3) = 1
   allocate(ext_data(1:nx(ik), 1:ny(ik), 1))
   do i=1, nx(ik)
      do j=1,ny(ik)
         ext_data(i,j,1) = z(offz(ik-1) + (i - 1) * ny (ik) + j)
      enddo
   enddo
   dimen_alloc(1:3) = dimen
   dimen_alloc(4)   = refine_par_n
   if(.not. gl_is_allocated()) call gl_alloc(dimen_alloc)
!
   CALL gl_set_data(dimen(1), dimen(2), ig, ext_data(:,:,1))
   DEALLOCATE(ext_data)
elseif(ku_ndims(ik)==1) then      ! 1D data set
   dimen(1) = lenc(ik)
   dimen(2) = 1
   dimen(3) = 1
   ALLOCATE(ext_data(1:lenc(ik), 1, 1))
   ext_data(1:lenc(ik), 1, 1) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
   dimen_alloc(1:3) = dimen
   dimen_alloc(4)   = refine_par_n
   if(.not. gl_is_allocated()) call gl_alloc(dimen_alloc)
!
   CALL gl_set_data(dimen(1), ig, ext_data(:,1,1))
   DEALLOCATE(ext_data)
ENDIF
endif
!
!  Add to data_struc
if(ku_ndims(ik)==3) then
   return
endif
!
call dgl5_new_node
!write(ik1_infile,'(a,a,i4.4)') opara(O_REFINE)(1:len_trim(opara(O_REFINE))), '.',nint(owerte(O_DATASET))
node_number = dgl5_get_number()
!write(*,*) ' added node ', node_number
if(ku_ndims(ik)==2) then      ! 2D data set  (NIPL)
   dimen(1) = nx(ik)
   dimen(2) = ny(ik)
   dimen(3) = 1
   if(allocated(ik1_x)) deallocate(ik1_x)
   if(allocated(ik1_y)) deallocate(ik1_y)
   if(allocated(ik1_z)) deallocate(ik1_z)
   if(allocated(ik1_dx)) deallocate(ik1_dx)
   if(allocated(ik1_dy)) deallocate(ik1_dy)
   if(allocated(ik1_dz)) deallocate(ik1_dz)
   if(allocated(ik1_data))  deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   allocate(ik1_x(1:dimen(1)))
   allocate(ik1_y(1:dimen(2)))
   allocate(ik1_z(1:1))
   allocate(ik1_dx(1:1))
   allocate(ik1_dy(1:1))
   allocate(ik1_dz(1:1))
   allocate(ik1_data(1:dimen(1), 1:dimen(2), 1:dimen(3)))
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
   ik1_x   (1:dimen(1))     = x(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
   ik1_y   (1:dimen(1))     = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
   do i=1, nx(ik)
      do j=1,ny(ik)
         ik1_data(i,j,1) = z(offz(ik-1) + (i - 1) * ny (ik) + j)
      enddo
   enddo
   ik1_sigma     = 0.0_PREC_DP
   ik1_llims     = 0.0_PREC_DP
   ik1_steps     = 0.0_PREC_DP
   ik1_steps_full= 0.0_PREC_DP
!
   call dgl5_set_node(ik1_infile, ik1_data_type, ik1_nlayer, ik1_is_direct, ku_ndims(ik), dimen ,         &
                   ik1_is_grid, ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, &
                   ik1_corners, ik1_vectors,&
                   ik1_a0(1:3), ik1_win(1:3), ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy,  &
                   ik1_dz,     ik1_data               ,ik1_sigma, ik1_llims,      &
                   ik1_steps, ik1_steps_full)
elseif(ku_ndims(ik)==1) then     ! Data set is 1D (NIPL)
   dimen(1) = lenc(ik)
   dimen(2) = 1
   dimen(3) = 1
   if(allocated(ik1_x)) deallocate(ik1_x)
   if(allocated(ik1_y)) deallocate(ik1_y)
   if(allocated(ik1_z)) deallocate(ik1_z)
   if(allocated(ik1_dx)) deallocate(ik1_dx)
   if(allocated(ik1_dy)) deallocate(ik1_dy)
   if(allocated(ik1_dz)) deallocate(ik1_dz)
   if(allocated(ik1_data))  deallocate(ik1_data)
   if(allocated(ik1_sigma)) deallocate(ik1_sigma)
   allocate(ik1_x(1:dimen(1)))
   allocate(ik1_y(1:1))
   allocate(ik1_z(1:1))
   allocate(ik1_dx(1:1))
   allocate(ik1_dy(1:1))
   allocate(ik1_dz(1:1))
   allocate(ik1_data(1:dimen(1), 1:dimen(2), 1:dimen(3)))
   allocate(ik1_sigma(1,1,1))
   ik1_data_type = H5_1D_GEN
   ik1_nlayer    = 1
   ik1_is_direct = .true.
   ik1_is_grid   = .true.
   ik1_has_dxyz  = .false.
   ik1_has_dval  = .false.
   ik1_calc_coor = .false.
   ik1_use_coor  = 1
   ik1_corners   = 0.0_PREC_DP
   ik1_vectors   = 0.0_PREC_DP
   ik1_x   (1:dimen(1))     = x(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
   ik1_y         = 0.0_PREC_DP
   ik1_z         = 0.0_PREC_DP
   ik1_data(1:dimen(1),1,1) = y(offxy(ik-1)+1:offxy(ik-1)+lenc(ik))
   ik1_sigma     = 0.0_PREC_DP
   ik1_llims     = 0.0_PREC_DP
   ik1_steps     = 0.0_PREC_DP
   ik1_steps_full= 0.0_PREC_DP
!
   call dgl5_set_node(ik1_infile, ik1_data_type, ik1_nlayer, ik1_is_direct, ku_ndims(ik), dimen ,         &
                   ik1_is_grid, ik1_has_dxyz, ik1_has_dval, ik1_calc_coor, ik1_use_coor, &
                   ik1_corners, ik1_vectors,&
                   ik1_a0(1:3), ik1_win(1:3), ik1_x, ik1_y, ik1_z, ik1_dx, ik1_dy,  &
                   ik1_dz,     ik1_data               ,ik1_sigma, ik1_llims,      &
                   ik1_steps, ik1_steps_full)
endif
END SUBROUTINE kuplot_to_global
!
!*******************************************************************************
!
END MODULE kuplot_toglobal
