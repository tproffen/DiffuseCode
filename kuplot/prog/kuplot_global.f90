module kuplot_global
!+
!  Interface KUPLOT <==> GLOBAL data structure
!+
!
contains
!
!*******************************************************************************
!
subroutine data2kuplot(ik, data_name, lout)
!-
! interface to the dimension specific routines
!+
use kuplot_mod
!
use errlist_mod
!
implicit none
!
integer         , intent(in) :: ik         ! Intended data number in KUPLOT
character(len=*), intent(in) :: data_name  ! 'file name'
logical         , intent(in) :: lout       ! Do show_data
!
if(ku_ndims(ik)==1) then
   call data2line(ik, data_name, lout   )
elseif(ku_ndims(ik)==2) then
   call data2nipl(ik, data_name, lout   )
elseif(ku_ndims(ik)==3) then
   call data3nipl(ik, data_name, lout   )
endif
!
end subroutine data2kuplot
!
!*******************************************************************************
!
subroutine data2line(ik, data_name, lout)
!-
!   Transfer the essential parts of a 1D data from global into the (old) KUPLOT 
!   arrays
!+
!
use kuplot_mod
use kuplot_show_mod
!
use errlist_mod
use lib_errlist_func
use lib_data_struc_h5
use precision_mod
use prompt_mod
!
implicit none
!
integer, intent(in) :: ik    ! Intended data number in KUPLOT
character(len=*), intent(in) :: data_name
logical         , intent(in) :: lout       ! Do show_data
!
integer :: i
integer :: length
!
character(len=PREC_STRING)                         :: infile
integer                                            :: node_number  ! Node in global data
integer                                            :: nlayer       ! Current layer (3-D only
logical                                            :: is_direct   ! Data are on direct / reciprocal scale
integer                                            :: ndims        ! Number of dimensions
integer, dimension(3)                              :: dims         ! Dimensions global array
logical                                            :: is_grid     ! Data are on direct / reciprocal scale
logical                                            :: has_dxyz    ! Data are on direct / reciprocal scale
logical                                            :: has_dval    ! Data are on direct / reciprocal scale
real(kind=PREC_DP), dimension(3,4)                 :: corners      ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: vectors      ! Steps        global array
real(kind=PREC_DP), dimension(3)                   :: a0              ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: win             ! Lower limits global array
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_x      ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_y        ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_z        ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dx       ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dy       ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dz       ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: odata        ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: sigma       ! Global data array for real
real(kind=PREC_DP), dimension(3)                   :: llims        ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: steps        ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: steps_full   ! Steps        global array
real(kind=PREC_DP), dimension(2)                   :: minmaxval    ! Lower limits global array
real(kind=PREC_DP), dimension(3,2)                 :: minmaxcoor   ! Lower limits global array
!
call no_error
!
node_number = dgl5_get_h5_is_ku(ik)
call data2local(ik, ier_num, ier_typ, node_number, infile, nlayer, is_direct,   &
     ndims, dims, is_grid, has_dxyz, has_dval, corners, vectors, a0, win,       &
     g_x, g_y, g_z, g_dx, g_dy, g_dz, odata, sigma, llims, steps,               &
     steps_full, minmaxval, minmaxcoor)
!
length =     dims(1)
if(ik<iz) then     ! Old data set to be overwritten, check dimensions
   if(nx(ik)>=dims(1)) then
      ier_num = -73
      ier_typ =  ER_APPL
   endif
endif
!
if(ier_num == 0) then
   offxy(ik  ) = offxy(ik-1) + length
!
   ku_ndims(ik  ) = 1
!
   do i =1, dims(1)
      x(offxy(ik-1) +i) = llims(1)   + (i-1)*steps(1)
      y(offxy(ik-1) +i) = odata(i, 1, 1)
   enddo
!
   fname(ik) = data_name
   fform(ik) = 'XY'
   nx  (ik)   = dims(1)
   xmin(ik  ) = x(offxy(ik-1) + 1)
   xmax(ik  ) = x(offxy(ik-1) + 1) + (dims(1) -1) * steps(1)
   ymin(ik  ) = minval(odata)
   ymax(ik  ) = maxval(odata)

   lni (ik)   = .false.
   lh5 (ik)   = .true.
   lenc(ik  ) = nx(ik)
   if(ik==iz) iz = iz + 1
   if(lout) then
      CALL show_data(ik)
   endif
endif
!
if(allocated(g_x)) deallocate(g_x)
if(allocated(g_y)) deallocate(g_y)
if(allocated(g_z)) deallocate(g_z)
if(allocated(g_dx)) deallocate(g_dx)
if(allocated(g_dy)) deallocate(g_dy)
if(allocated(g_dz)) deallocate(g_dz)
if(allocated(odata)) deallocate(odata)
if(allocated(sigma)) deallocate(sigma)
!
end subroutine data2line
!
!*******************************************************************************
!
subroutine data2nipl(ik, data_name, lout)
!-
!   Transfer the essential parts of a 2D data from global into the (old) KUPLOT 
!   arrays
!+
!
use kuplot_mod
use kuplot_show_mod
!
use errlist_mod
use lib_errlist_func
use lib_data_struc_h5
use precision_mod
use prompt_mod
!
implicit none
!
integer, intent(in) :: ik    ! Intended data number in KUPLOT
character(len=*), intent(in) :: data_name
logical         , intent(in) :: lout       ! Do show_data
!
integer :: i, j, ii
integer :: length
!
character(len=PREC_STRING)                         :: infile
integer                                            :: node_number  ! Node in global data
integer                                            :: nlayer       ! Current layer (3-D only
logical                                            :: is_direct   ! Data are on direct / reciprocal scale
integer                                            :: ndims        ! Number of dimensions
integer, dimension(3)                              :: dims         ! Dimensions global array
logical                                            :: is_grid     ! Data are on direct / reciprocal scale
logical                                            :: has_dxyz    ! Data are on direct / reciprocal scale
logical                                            :: has_dval    ! Data are on direct / reciprocal scale
real(kind=PREC_DP), dimension(3,4)                 :: corners      ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: vectors      ! Steps        global array
real(kind=PREC_DP), dimension(3)                   :: a0              ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: win             ! Lower limits global array
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_x      ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_y        ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_z        ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dx       ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dy       ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dz       ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: odata        ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: sigma       ! Global data array for real
real(kind=PREC_DP), dimension(3)                   :: llims        ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: steps        ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: steps_full   ! Steps        global array
real(kind=PREC_DP), dimension(2)                   :: minmaxval    ! Lower limits global array
real(kind=PREC_DP), dimension(3,2)                 :: minmaxcoor   ! Lower limits global array
!
call no_error
node_number = dgl5_get_h5_is_ku(ik)
!
call data2local(ik, ier_num, ier_typ, node_number, infile, nlayer, is_direct,   &
     ndims, dims, is_grid, has_dxyz, has_dval, corners, vectors, a0, win,       &
     g_x, g_y, g_z, g_dx, g_dy, g_dz, odata, sigma, llims, steps,               &
     steps_full, minmaxval, minmaxcoor)
!
length = max(dims(1), dims(2))
if(ik<iz .and. iz>1) then     ! Old data set to be overwritten, check dimensions
   if(nx(ik)>dims(1) .or. ny(ik)>dims(2)) then
      ier_num = -73
      ier_typ =  ER_APPL
   endif
endif
!
if(ier_num == 0) then
   offxy(ik  ) = offxy(ik-1) + length
   offz (ik  ) = offz (ik-1) + dims(1)*dims(2)
!
   ku_ndims(ik  ) = 2
!
   do i =1, dims(1)
      x(offxy(ik-1) +i) = llims(1)   + (i-1)*steps(1)
   enddo
!
   do i =1, dims(2)
      y(offxy(ik-1) +i) = llims(2)   + (i-1)*steps(2)
   enddo
!
!  Reminder, data are stored in z along y-collumns!
   ii = 0
   do i=1, dims(1)
      do j=1, dims(2)
         ii = ii + 1
         z(offz(ik-1)+ii) = (odata(i, j, 1)) !
      enddo
   enddo
!
   fname(ik) = data_name
   fform(ik) = 'NI'
   nx  (ik)   = dims(1)
   ny  (ik)   = dims(2)
   xmin(ik  ) = x(offxy(ik-1) + 1)
   ymin(ik  ) = y(offxy(ik-1) + 1)
   xmax(ik  ) = x(offxy(ik-1) + 1) + (dims(1) -1) * steps(1)
   ymax(ik  ) = y(offxy(ik-1) + 1) + (dims(2) -1) * steps(2)
   lni (ik)   = .TRUE.
   lh5 (ik)   = .true.
   lenc(ik  ) = max(nx  (ik), ny(ik))
   if(ik==iz) iz = iz + 1
   if(lout) then
      CALL show_data(ik)
   endif
endif
!
if(allocated(g_x)) deallocate(g_x)
if(allocated(g_y)) deallocate(g_y)
if(allocated(g_z)) deallocate(g_z)
if(allocated(g_dx)) deallocate(g_dx)
if(allocated(g_dy)) deallocate(g_dy)
if(allocated(g_dz)) deallocate(g_dz)
if(allocated(odata)) deallocate(odata)
if(allocated(sigma)) deallocate(sigma)
!
end subroutine data2nipl
!
!*******************************************************************************
!
subroutine data3nipl(ik, data_name, lout)
!-
!   Transfer the essential parts of a 3D data from global into the (old) KUPLOT 
!   arrays
!+
!
use kuplot_mod
use kuplot_show_mod
!
use errlist_mod
use lib_errlist_func
use lib_data_struc_h5
use precision_mod
use prompt_mod
!
implicit none
!
integer         , intent(in) :: ik         ! Intended data number in KUPLOT
character(len=*), intent(in) :: data_name
logical         , intent(in) :: lout       ! Do show_data
!
integer :: i, j, ii
integer :: length
!
character(len=PREC_STRING)                         :: infile
integer                                            :: node_number  ! Node in global data
integer                                            :: nlayer       ! Current layer (3-D only
logical                                            :: is_direct   ! Data are on direct / reciprocal scale
integer                                            :: ndims        ! Number of dimensions
integer, dimension(3)                              :: dims         ! Dimensions global array
logical                                            :: is_grid     ! Data are on direct / reciprocal scale
logical                                            :: has_dxyz    ! Data are on direct / reciprocal scale
logical                                            :: has_dval    ! Data are on direct / reciprocal scale
real(kind=PREC_DP), dimension(3,4)                 :: corners      ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: vectors      ! Steps        global array
real(kind=PREC_DP), dimension(3)                   :: a0              ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: win             ! Lower limits global array
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_x      ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_y        ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_z        ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dx       ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dy       ! Global data array for real
real(kind=PREC_DP), dimension(:)    , allocatable  :: g_dz       ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: odata        ! Global data array for real
real(kind=PREC_DP), dimension(:,:,:), allocatable  :: sigma       ! Global data array for real
real(kind=PREC_DP), dimension(3)                   :: llims        ! Lower limits global array
real(kind=PREC_DP), dimension(3)                   :: steps        ! Steps        global array
real(kind=PREC_DP), dimension(3,3)                 :: steps_full   ! Steps        global array
real(kind=PREC_DP), dimension(2)                   :: minmaxval    ! Lower limits global array
real(kind=PREC_DP), dimension(3,2)                 :: minmaxcoor   ! Lower limits global array
!
call no_error
!
node_number = dgl5_get_h5_is_ku(ik)
call data2local(ik, ier_num, ier_typ, node_number, infile, nlayer, is_direct,   &
     ndims, dims, is_grid, has_dxyz, has_dval, corners, vectors, a0, win,       &
     g_x, g_y, g_z, g_dx, g_dy, g_dz, odata, sigma, llims, steps,               &
     steps_full, minmaxval, minmaxcoor)
!
length = max(dims(1), dims(2))
if(ik<iz) then     ! Old data set to be overwritten, check dimensions
   if(nx(ik)>dims(1) .or. ny(ik)>dims(2)) then
      ier_num = -73
      ier_typ =  ER_APPL
   endif
endif
!
if(ier_num == 0) then
   offxy(ik  ) = offxy(ik-1) + length
   offz (ik  ) = offz (ik-1) + dims(1)*dims(2)
!
   ku_ndims(ik  ) = 3
!
   do i =1, dims(1)
      x(offxy(ik-1) +i) = llims(1)   + (i-1)*steps(1)
   enddo
!
   do i =1, dims(2)
      y(offxy(ik-1) +i) = llims(2)   + (i-1)*steps(2)
   enddo
!
!  Reminder, data are stored in z along y-collumns!
   ii = 0
   do i=1, dims(1)
      do j=1, dims(2)
         ii = ii + 1
         z(offz(ik-1)+ii) = (odata(i, j, nlayer)) !
      enddo
   enddo
!
   fname(ik) = data_name
   fform(ik) = 'NI'
   nx  (ik)   = dims(1)
   ny  (ik)   = dims(2)
   xmin(ik  ) = x(offxy(ik-1) + 1)
   ymin(ik  ) = y(offxy(ik-1) + 1)
   xmax(ik  ) = x(offxy(ik-1) + 1) + (dims(1) -1) * steps(1)
   ymax(ik  ) = y(offxy(ik-1) + 1) + (dims(2) -1) * steps(2)
   lni (ik)   = .TRUE.
   lh5 (ik)   = .true.
   lenc(ik  ) = max(nx  (ik), ny(ik))
   if(ik==iz) iz = iz + 1
!
   if(lout) then
      CALL show_data(ik)
      write(output_io,'(a,i7)') '   At  layer:', nlayer
   endif
endif
!
if(allocated(g_x)) deallocate(g_x)
if(allocated(g_y)) deallocate(g_y)
if(allocated(g_z)) deallocate(g_z)
if(allocated(g_dx)) deallocate(g_dx)
if(allocated(g_dy)) deallocate(g_dy)
if(allocated(g_dz)) deallocate(g_dz)
if(allocated(odata)) deallocate(odata)
if(allocated(sigma)) deallocate(sigma)
!
end subroutine data3nipl
!
!*******************************************************************************
!
end module kuplot_global
