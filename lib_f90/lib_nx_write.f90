module nx_write_mod
!-
!  Subroutines to write into a NEXUS file
!  DiffuseDevelopers Software convention
!+
private
!public 
public nx_write_scattering  ! Subroutine; write the provided scattering data into Nexus file
!
contains
!
!*******************************************************************************
!
subroutine nx_write_scattering(value, laver, outfile, out_inc, out_eck, out_vi, &
           extr_abs, extr_ord, extr_top,                                        &
           cr_a0, cr_win, qvalues, VAL_PDF, VAL_3DPDF, valmax,                  &
           ier_num, ier_typ, ER_IO, ER_APPL)
!-
!  Main writing routine for Nexus Format DiffuseDevelopers Scattering format
!+
!
!use errlist_mod
use lib_errlist_func
use precision_mod
!
use lib_forpython_mod
use forpy_mod
!
implicit none
!
integer                           , intent(in)  :: value     ! Output value type (inte, pdf etc)
logical                           , intent(in)  :: laver     ! 
character(len=*)                  , intent(in)  :: outfile   ! Output file name
integer           , dimension(3)  , intent(in)  :: out_inc   ! Dimensions of qvalue
real(kind=PREC_DP), dimension(3,4), intent(in)  :: out_eck   ! (:,1) LowerLeft
real(kind=PREC_DP), dimension(3,3), intent(in)  :: out_vi    ! (:,1) Abscissa
integer                           , intent(in)  :: extr_abs  ! Abscissa is 1=h 2=k 3=l
integer                           , intent(in)  :: extr_ord  ! Ordinate is "
integer                           , intent(in)  :: extr_top  ! Top axis is "
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_a0     ! Lattice params (a, b, c)
real(kind=PREC_DP), dimension(3)  , intent(in)  :: cr_win    ! Lattice params (alpha, beta, gamma)
real(kind=PREC_DP), dimension(out_inc(1), out_inc(2), out_inc(3)), intent(in) :: qvalues
integer                           , intent(in)  :: VAL_PDF   ! Indicator for a 3D-PDF
integer                           , intent(in)  :: VAL_3DPDF ! Indicator for a 3D-Delta-PDF
real(KinD=PREC_DP)                , intent(IN)  :: valmax    ! Intended maximum output value
integer                           , intent(out) :: ier_num   ! Error number
integer                           , intent(out) :: ier_typ   ! Error type
integer                           , intent(in)  :: ER_IO     ! Error type input/output
integer                           , intent(in)  :: ER_APPL   ! Error type Application
!
type(tuple)     :: p_args          ! Tuple of arguments for write_diffuse_scattering
type(object )   :: p_outfile       ! Output filename in python interface
type(object )   :: p_program       ! Output program  in python interface
type(object )   :: p_author        ! Output author   in python interface
type(object )   :: p_qvalues_dim   ! Data dimensionality in python interface
type(object )   :: p_radiation     ! radiation type 
type(object )   :: p_space         ! reciprocal, direct space
type(object )   :: p_content       ! intensity, 3d-delta-PDF, etc
type(object )   :: p_axes          ! axes labels ["h", "k","l"] etc
type(ndarray)   :: p_qvalues       ! Intensities     in python interface
type(ndarray)   :: p_lower_limits  ! Left_lower_bottom corner in python  interface
type(ndarray)   :: p_step_vectors  ! Increment vi along abscissa python  interface
type(ndarray)   :: p_h_indices     ! H indices       in python interface
type(ndarray)   :: p_k_indices     ! K indices       in python interface
type(ndarray)   :: p_l_indices     ! L indices       in python interface
type(ndarray)   :: p_unit_cell_length  ! Unit cell length in python interface
type(ndarray)   :: p_unit_cell_angle   ! Unit cell angles in python interface
!
type(module_py) :: interface_module   ! python script name
type(list)      :: paths_to_module    ! python script path
type(object)    :: return_value       ! forpy return value
!character(len=:), allocatable :: return_string         ! Result from forpy
!
character(len=PREC_STRING) :: radiation
character(len=PREC_STRING) :: space
character(len=PREC_STRING) :: content
character(len=PREC_STRING) :: axes
integer :: i, j, k   ! Dummy indices
integer :: qvalues_dim  ! Data dimensionality
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: values        ! Values in C sequence
real(kind=PREC_DP), dimension(out_inc(extr_abs)) :: h_indices      ! Actual H indices along a*
real(kind=PREC_DP), dimension(out_inc(extr_ord)) :: k_indices      ! Actual K indices along b*
real(kind=PREC_DP), dimension(out_inc(extr_top)) :: l_indices      ! Actual L indices along c*
!
write(*,*) ' WRITE NX '
qvalues_dim = 3     ! For right now fixed 3D
radiation = 'xray'
space     = 'reciprocal'
content   = 'intensity'
axes      = '["h", "k", "l"]'
write(*,*) ' transpose axes ', out_inc
!
! Transpose input values 
!
allocate(values(out_inc(3), out_inc(2), out_inc(1)))
write(*,*) ' transpose axes ', ubound(qvalues), ubound(values)
DO i = 1, out_inc(1)
   DO j = 1, out_inc(2)
      DO k = 1, out_inc(3)
         values(k,j,i) = qvalues(i,j,k)
      ENDDO
   ENDDO
ENDDO
write(*,*) ' Wrote original   array 321 ', qvalues(3,2,1)
write(*,*)
write(*,*) ' Wrote original   array llb ', qvalues(1         ,1         ,1         )
write(*,*) ' Wrote original   array rlb ', qvalues(out_inc(1),1         ,1         )
write(*,*) ' Wrote original   array lub ', qvalues(1         ,out_inc(2),1         )
write(*,*) ' Wrote original   array rub ', qvalues(out_inc(1),out_inc(2),1         )
!
write(*,*) ' Wrote original   array llt ', qvalues(1         ,1         ,out_inc(3))
write(*,*) ' Wrote original   array rlt ', qvalues(out_inc(1),1         ,out_inc(3))
write(*,*) ' Wrote original   array lut ', qvalues(1         ,out_inc(2),out_inc(3))
write(*,*) ' Wrote original   array rut ', qvalues(out_inc(1),out_inc(2),out_inc(3))
!
! Write full indices H, K, L along axes
!
write(*,*) ' extr_abs ', extr_abs, extr_ord, extr_top
do i=1, out_inc(extr_abs)
   h_indices(i) = out_eck(1,1) + (i-1)*out_vi(extr_abs,1)
enddo
do i=1, out_inc(extr_ord)
   k_indices(i) = out_eck(2,1) + (i-1)*out_vi(extr_ord,2)
enddo
do i=1, out_inc(extr_top)
   l_indices(i) = out_eck(3,1) + (i-1)*out_vi(extr_top,3)
enddo
!
call no_error
call forpy_start()
!
!
ier_num = cast(          p_outfile,     outfile(1:len_trim(outfile)) )
ier_num = cast(          p_program,     'DISCUS 6.17.02')
ier_num = cast(          p_author,      'rbn')
ier_num = cast(          p_qvalues_dim, qvalues_dim)
ier_num = ndarray_create(p_qvalues,     qvalues)
ier_num = cast(          p_radiation,   radiation(1:len_trim(radiation)) )
ier_num = cast(          p_space,       space(1:len_trim(space)) )
ier_num = cast(          p_content,     content(1:len_trim(content)) )
ier_num = ndarray_create(p_lower_limits,out_eck(:,1))
ier_num = ndarray_create(p_step_vectors,(out_vi))
ier_num = cast          (p_axes        ,axes(1:len_trim(axes)))
ier_num = ndarray_create(p_h_indices,   h_indices)
ier_num = ndarray_create(p_k_indices,   k_indices)
ier_num = ndarray_create(p_l_indices,   l_indices)
ier_num = ndarray_create(p_unit_cell_length, cr_a0)
ier_num = ndarray_create(p_unit_cell_angle,  cr_win)
!
!  Collect the 16 arguments into a tuple: p_args
!
ier_num = tuple_create(p_args, 16)
ier_num = p_args%setitem( 0, p_outfile)
ier_num = p_args%setitem( 1, p_program)
ier_num = p_args%setitem( 2, p_author)
ier_num = p_args%setitem( 3, p_qvalues_dim)
ier_num = p_args%setitem( 4, p_qvalues)
ier_num = p_args%setitem( 5, p_radiation)
ier_num = p_args%setitem( 6, p_space)
ier_num = p_args%setitem( 7, p_content)
ier_num = p_args%setitem( 8, p_lower_limits)
ier_num = p_args%setitem( 9, p_step_vectors)
ier_num = p_args%setitem(10, p_axes)
ier_num = p_args%setitem(11, p_h_indices)
ier_num = p_args%setitem(12, p_k_indices)
ier_num = p_args%setitem(13, p_l_indices)
ier_num = p_args%setitem(14, p_unit_cell_length)
ier_num = p_args%setitem(15, p_unit_cell_angle)
!ier_num = print_py(p_args)
!
! Append current directory to paths
!
ier_num = get_sys_path(paths_to_module)
ier_num = paths_to_module%append('.')
ier_num = import_py(interface_module, 'write_diffuse_scattering')
write(*,*) ' ier_num C ', ier_num
ier_num = call_py(return_value, interface_module, 'write_diffuse_scattering', p_args)
write(*,*) ' ier_num D ', ier_num
!
call p_outfile%destroy
call p_program%destroy
call p_author%destroy
call p_qvalues_dim%destroy
call p_qvalues%destroy
call p_radiation%destroy
call p_space%destroy
call p_content%destroy
call p_lower_limits%destroy
call p_step_vectors%destroy
call p_axes%destroy
call p_h_indices%destroy
call p_k_indices%destroy
call p_l_indices%destroy
call p_unit_cell_length%destroy
call p_unit_cell_angle%destroy
call p_args%destroy
call interface_module%destroy
!
end subroutine nx_write_scattering
!
!*******************************************************************************
end module nx_write_mod
