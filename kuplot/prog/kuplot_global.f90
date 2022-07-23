module kuplot_global
!+
!  Interface KUPLOT <==> GLOBAL data structure
!+
!
contains
!
!*******************************************************************************
!
subroutine data2local(ik, ier_num, ier_typ, node_number, nlayer, ndims, dims, odata)
!-
! Transfer a data set "ik" into the N-Dim local array odata
!+
!
use lib_data_struc_h5
use precision_mod
!
implicit none
!
integer, intent(in ) :: ik
integer, intent(out) :: ier_num
integer, intent(out) :: ier_typ
integer, intent(out) :: node_number
integer, intent(out) :: nlayer
integer,               intent(out) :: ndims
integer, dimension(3), intent(out) :: dims
real(kind=PREC_DP), dimension(:,:,:), allocatable, intent(out) :: odata

call hdf5_set_pointer(ik, ier_num, ier_typ, node_number)
nlayer = hdf5_get_layer()
if(ier_num/=0) return
!
ndims = hdf5_get_ndims()
call hdf5_get_dims(node_number, dims)
allocate(odata(dims(1), dims(2), dims(3)))
call hdf5_get_map(dims, odata)
!
end subroutine data2local
!
!*******************************************************************************
!
end module kuplot_global
