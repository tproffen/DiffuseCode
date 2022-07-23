MODULE lib_data_struc_h5
!-
!  Contains routines to build the lib_f90 data structure, based on HDF5
!+
USE lib_hdf5_params_mod
!
!USE kuplot_config
!
use hdf5_def_mod
USE precision_mod
!
IMPLICIT NONE
!
PRIVATE
!PUBLIC hdf5_read
PUBLIC hdf5_new_node
PUBLIC hdf5_set_node
!PUBLIC hdf5_place_kuplot
PUBLIC hdf5_set_pointer
PUBLIC hdf5_find_node
PUBLIC hdf5_copy_node
public hdf5_get_h5_is_ku
public hdf5_get_ku_is_h5
PUBLIC hdf5_get_layer
PUBLIC hdf5_get_height
PUBLIC hdf5_get_direct
PUBLIC hdf5_get_dims
PUBLIC hdf5_get_llims
PUBLIC hdf5_get_steps
PUBLIC hdf5_get_map
PUBLIC hdf5_get_tmap
public hdf5_get_data
public hdf5_get_number
public hdf5_get_ndims
PUBLIC hdf5_get_infile
public hdf5_set_h5_is_ku
public hdf5_set_ku_is_h5
PUBLIC hdf5_set_layer
PUBLIC hdf5_set_map
PUBLIC hdf5_reset
!
integer, parameter ::  maxkurvtot = 200      ! Relic from KUPLOT
!
!HARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: h5_datasets       ! Names of the data set in file
!CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
INTEGER, DIMENSION(MAXKURVTOT)                        :: h5_ku_is_h5 = 0   ! Pointer from kuplot number to h5 number
INTEGER, DIMENSION(MAXKURVTOT)                        :: h5_h5_is_ku = 0   ! Pointer from h5 number to kuplot number
INTEGER                                               :: h5_number   = 0   ! Currently loaded h5 data sets
!
TYPE :: h5_data_struc
   INTEGER                                               :: h5_data_num       ! Current data set number
   CHARACTER(LEN=PREC_STRING)                            :: h5_infile         ! input file
   INTEGER                                               :: h5_layer=1        ! Current layer in data set
   LOGICAL                                               :: h5_direct         ! Direct space == TRUE
   INTEGER                                               :: ndims             ! Number of dimensions
   INTEGER                                               :: one_ndims         ! Number of dimensions
   INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: h5_dims           ! Actual dimensions
   INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: maxdims           ! Maximum dimensions
   INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_dims          ! Actual dimensions
   INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_maxdims       ! Maximum dimensions
   REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_data           ! Actual diffraction data
   REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_llims          ! Lower limits
   REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_steps          ! steps in H, K, L
   REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_steps_full     ! steps in H, K, L
   TYPE(h5_data_struc), POINTER                          :: after
END TYPE h5_data_struc
!
TYPE(h5_data_struc), POINTER                          :: h5_root => NULL()
TYPE(h5_data_struc), POINTER                          :: h5_temp => NULL()
TYPE(h5_data_struc), POINTER                          :: h5_find => NULL()
!
!INTEGER                                               :: H5_MAX_DATASETS   ! Current MAX data sets
!INTEGER                                               :: h5_n_datasets     ! Current actual data sets
!INTEGER                                               :: h5_layer=1        ! Current layer in data set
!LOGICAL                                               :: h5_direct         ! Direct space == TRUE
!INTEGER                                               :: ndims             ! Number of dimensions
!INTEGER                                               :: one_ndims         ! Number of dimensions
!INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: h5_dims           ! Actual dimensions
!INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: maxdims           ! Maximum dimensions
!INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_dims          ! Actual dimensions
!INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3)               :: one_maxdims       ! Maximum dimensions
!REAL(KIND=PREC_SP)   , DIMENSION(:,:,:), ALLOCATABLE  :: h5_data           ! Actual diffraction data
!REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_llims          ! Lower limits
!REAL(KIND=PREC_DP)   , DIMENSION(3)                   :: h5_steps          ! steps in H, K, L
!REAL(KIND=PREC_DP)   , DIMENSION(3,3)                 :: h5_steps_full     ! steps in H, K, L
!                                                                          ! 1.st dim: [hkl], 2nd:[abs, ord, top]
!
contains
!
!*******************************************************************************
!
SUBROUTINE hdf5_new_node
!-
! Create a new node
!+
!
IF(ASSOCIATED(h5_root)) THEN                                ! A root node exists
   h5_temp => h5_root                                       ! Point to current node
   find_node: DO WHILE(ASSOCIATED(h5_temp%after))           ! Does next node exist?
      h5_temp => h5_temp%after                              ! Next node exists, point to this next node
   ENDDO find_node
   ALLOCATE(h5_temp%after)                                  ! Create next node
   h5_temp => h5_root%after                                 ! Point to Current working node
ELSE
   ALLOCATE(h5_root)                                        ! Create first node
   NULLIFY(h5_root%after)
   h5_temp => h5_root                                       ! Point to Current working node
ENDIF
! Work on current node
NULLIFY(h5_temp%after)
h5_temp%h5_data_num = h5_number + 1                         ! Increment the data number
h5_number = h5_number + 1                                   ! Increment the global data number
!
!write(*,*) ' ROOT NW', ASSOCIATED(h5_root), ASSOCIATED(h5_temp), h5_temp%h5_data_num, h5_root%h5_data_num, h5_number
!
END SUBROUTINE hdf5_new_node
!
!*******************************************************************************
!
SUBROUTINE hdf5_copy_node(old, new)
!-
!  Copies old node to new node
!+
INTEGER, INTENT(IN)  :: old
INTEGER, INTENT(OUT) :: new
INTEGER :: ier_num
INTEGER :: ier_typ
!
CALL hdf5_find_node(old, ier_num, ier_typ)
CALL hdf5_new_node
h5_temp%h5_infile   = h5_find%h5_infile         ! input file
h5_temp%h5_layer    = h5_find%h5_layer          ! Current layer in data set
h5_temp%h5_direct   = h5_find%h5_direct         ! Direct space == TRUE
h5_temp%ndims       = h5_find%ndims          ! Number of dimensions
h5_temp%one_ndims   = h5_find%one_ndims      ! Number of dimensions
h5_temp%h5_dims     = h5_find%h5_dims           ! Actual dimensions
h5_temp%maxdims     = h5_find%maxdims        ! Maximum dimensions
h5_temp%one_dims    = h5_find%one_dims       ! Actual dimensions
h5_temp%one_maxdims = h5_find%one_maxdims    ! Maximum dimensions
h5_temp%h5_data     = h5_find%h5_data           ! Actual diffraction data
h5_temp%h5_llims    = h5_find%h5_llims          ! Lower limits
h5_temp%h5_steps    = h5_find%h5_steps          ! steps in H, K, L
!
new = h5_temp%h5_data_num
!
END SUBROUTINE hdf5_copy_node
!
!*******************************************************************************
!
SUBROUTINE hdf5_set_node(l_infile, l_layer, l_direct, l_ndims, l_one_ndims, l_dims, &
                   l_maxdims, l_one_dims, l_one_maxdims, l_data, l_llims, l_steps,  &
                   l_steps_full)
!-
!  Place the temporary values into the current hdf5 node
!+
CHARACTER(LEN=*)                   , INTENT(IN)       :: l_infile         ! Input file
INTEGER                            , INTENT(IN)       :: l_layer          ! Current layer in data set
LOGICAL                            , INTENT(IN)       :: l_direct         ! Direct space == TRUE
INTEGER                            , INTENT(IN)       :: l_ndims          ! Number of dimensions
INTEGER                            , INTENT(IN)       :: l_one_ndims      ! Number of dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3), INTENT(IN)   :: l_dims           ! Actual dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3), INTENT(IN)   :: l_maxdims        ! Maximum dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3), INTENT(IN)   :: l_one_dims       ! Actual dimensions
INTEGER(KIND=LIB_HSIZE_T), DIMENSION(3), INTENT(IN)   :: l_one_maxdims    ! Maximum dimensions
REAL(KIND=PREC_SP)   , DIMENSION(l_dims(1), l_dims(2), l_dims(3)), INTENT(IN):: l_data           ! Actual diffraction data
REAL(KIND=PREC_DP)   , DIMENSION(3), INTENT(IN)       :: l_llims          ! Lower limits
REAL(KIND=PREC_DP)   , DIMENSION(3), INTENT(IN)       :: l_steps          ! steps in H, K, L
REAL(KIND=PREC_DP)   , DIMENSION(3,3),INTENT(IN)       :: l_steps_full     ! steps in H, K, L
!
h5_temp%h5_infile   = l_infile         ! input file
h5_temp%h5_layer    = l_layer          ! Current layer in data set
h5_temp%h5_direct   = l_direct         ! Direct space == TRUE
h5_temp%ndims       = l_ndims          ! Number of dimensions
h5_temp%one_ndims   = l_one_ndims      ! Number of dimensions
h5_temp%h5_dims     = l_dims           ! Actual dimensions
h5_temp%maxdims     = l_maxdims        ! Maximum dimensions
h5_temp%one_dims    = l_one_dims       ! Actual dimensions
h5_temp%one_maxdims = l_one_maxdims    ! Maximum dimensions
!ALLOCATE(h5_temp%h5_data(h5_temp%h5_dims(1), h5_temp%h5_dims(2), h5_temp%h5_dims(3)))
h5_temp%h5_data     = l_data           ! Actual diffraction data
h5_temp%h5_llims    = l_llims          ! Lower limits
h5_temp%h5_steps    = l_steps          ! steps in H, K, L
h5_temp%h5_steps_full    = l_steps_full          ! steps in H, K, L
!
END SUBROUTINE hdf5_set_node
!
!*******************************************************************************
!
SUBROUTINE hdf5_set_pointer(izz, ier_num, ier_typ, node_number)
!-
!  Find the node associated to kuplot data set number izz
!+
INTEGER, INTENT(IN ) :: izz
INTEGER, INTENT(OUT) :: ier_num
INTEGER, INTENT(OUT) :: ier_typ
INTEGER, INTENT(OUT) :: node_number
!
!write(*,*) ' ROOT ST', ASSOCIATED(h5_root), ASSOCIATED(h5_temp), node_number, h5_number
IF(.NOT. ASSOCIATED(h5_root)) THEN
   ier_num = -74         ! Root node does not exist !
   ier_typ =   6         ! ER_APPL
   RETURN
ENDIF
!
h5_temp => h5_root
!write(*,*) ' DATA NUMB ', h5_root%h5_data_num, h5_temp%h5_data_num 
find_node: DO            ! Search for node
!write(*,*) ' SEARCHING ', izz, h5_ku_is_h5(izz), h5_temp%h5_data_num, h5_ku_is_h5(izz) == h5_temp%h5_data_num
   IF(h5_ku_is_h5(izz) == h5_temp%h5_data_num) THEN
      node_number = h5_temp%h5_data_num
      EXIT find_node
   ELSE
      IF(ASSOCIATED(h5_temp%after)) THEN    ! A next node exists 
         h5_temp => h5_temp%after
      ELSE
         ier_num = -74         ! Root node does not exist !
         ier_typ =   6         ! ER_APPL
         RETURN
      ENDIF
   ENDIF
ENDDO find_node
!
END SUBROUTINE hdf5_set_pointer
!
!*******************************************************************************
!
SUBROUTINE hdf5_find_node(node_number, ier_num, ier_typ)
!-
!  Find the node with node_number
!+
INTEGER, INTENT(IN)  :: node_number
INTEGER, INTENT(OUT) :: ier_num
INTEGER, INTENT(OUT) :: ier_typ
!
h5_find => h5_root
find_node: DO            ! Search for node
   IF(node_number == h5_find%h5_data_num) THEN
      EXIT find_node
   ELSE
      IF(ASSOCIATED(h5_find%after)) THEN    ! A next node exists 
         h5_find => h5_find%after
      ELSE
         ier_num = -74         ! Root node does not exist !
         ier_typ =   6         ! ER_APPL
         RETURN
      ENDIF
   ENDIF
ENDDO find_node
!
END SUBROUTINE hdf5_find_node
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_layer()
!
IMPLICIT NONE
!
hdf5_get_layer = h5_temp%h5_layer
!
END FUNCTION hdf5_get_layer
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_number()
!
IMPLICIT NONE
!
hdf5_get_number = h5_number
!
END FUNCTION hdf5_get_number
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_ndims()
!
IMPLICIT NONE
!
hdf5_get_ndims = h5_temp%ndims
!
END FUNCTION hdf5_get_ndims
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_ku_is_h5(izz)
!
IMPLICIT NONE
!
integer,intent(in) :: izz
!
hdf5_get_ku_is_h5 = h5_ku_is_h5(izz)
!
END FUNCTION hdf5_get_ku_is_h5
!
!*******************************************************************************
!
INTEGER FUNCTION hdf5_get_h5_is_ku(inumber)
!
IMPLICIT NONE
!
integer,intent(in) :: inumber
!
hdf5_get_h5_is_ku = h5_h5_is_ku(inumber)
!
END FUNCTION hdf5_get_h5_is_ku
!
!*******************************************************************************
!
subroutine hdf5_set_layer(h5_layer)
!
IMPLICIT NONE
!
integer, intent(in) :: h5_layer
!
h5_temp%h5_layer = h5_layer
!
END subroutine hdf5_set_layer
!
!*******************************************************************************
!
LOGICAL FUNCTION hdf5_get_direct()
!
IMPLICIT NONE
!
hdf5_get_direct = h5_temp%h5_direct
!
END FUNCTION hdf5_get_direct
!
!*******************************************************************************
!
REAL FUNCTION hdf5_get_height()
!
IMPLICIT NONE
!
hdf5_get_height = h5_temp%h5_llims(3) + (h5_temp%h5_layer-1)*h5_temp%h5_steps(3)
!
END FUNCTION hdf5_get_height
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_dims(idata, dims)
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN)  :: idata
INTEGER, DIMENSION(3), INTENT(OUT) :: dims
!
dims = h5_temp%h5_dims
!write(*,*) 'GOT DIMS ', ubound(h5_temp%h5_data)
!
END SUBROUTINE hdf5_get_dims
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_llims(idata, llims)
!
use precision_mod
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN)  :: idata
real(kind=PREC_DP), DIMENSION(3), INTENT(OUT) :: llims
!
llims = h5_temp%h5_llims
!
END SUBROUTINE hdf5_get_llims
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_steps(idata, steps)
!
use hdf5_def_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN)  :: idata
real(kind=PREC_DP), DIMENSION(3,3), INTENT(OUT) :: steps
!
!write(*,*) yd_present(YD_step_sizes_abs:YD_step_sizes_TOP), &
!       ALL(yd_present(YD_step_sizes_abs:YD_step_sizes_TOP))
if(ALL(yd_present(YD_step_sizes_abs:YD_step_sizes_TOP))) then
   steps = h5_temp%h5_steps_full
else
   steps = 0.0D0
   steps(1,1) = h5_temp%h5_steps(1)
   steps(2,2) = h5_temp%h5_steps(2)
   steps(3,3) = h5_temp%h5_steps(3)
endif
!write(*,*) ' GETS ', steps(:,1)
!write(*,*) ' GETS ', steps(:,2)
!write(*,*) ' GETS ', steps(:,3)
!
END SUBROUTINE hdf5_get_steps
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_map(dims, odata)
!
IMPLICIT NONE
!
INTEGER,            DIMENSION(3),                         INTENT(IN)  :: dims
REAL(KIND=PREC_DP), DIMENSION(dims(1), dims(2), dims(3)), INTENT(OUT) :: odata
!
INTEGER :: i,j,k
!
DO i=1, dims(1)
   DO j=1, dims(2)
      DO k=1, dims(3)
         odata(i,j,k) = h5_temp%h5_data(i,j,k)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE hdf5_get_map
!
!*******************************************************************************
!
SUBROUTINE hdf5_get_tmap(dims, odata)
!-
! Get the matrix in transposed form, which is the regular Fortran style
!+
!
IMPLICIT NONE
!
INTEGER,            DIMENSION(3),                         INTENT(IN)  :: dims
REAL(KIND=PREC_DP), DIMENSION(dims(1), dims(2), dims(3)), INTENT(OUT) :: odata
!
INTEGER :: i,j,k
!
DO i=1, dims(1)
   DO j=1, dims(2)
      DO k=1, dims(3)
         odata(i,j,k) = h5_temp%h5_data(k,j,i)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE hdf5_get_tmap
!
!*******************************************************************************
!
real(kind=PREC_DP) function hdf5_get_data(i,j,k)
!
!  Get a single data point from the data structure
!
integer, intent(in) :: i
integer, intent(in) :: j
integer, intent(in) :: k
!
hdf5_get_data = h5_temp%h5_data(k,j,i)
!
end function hdf5_get_data
!
!*******************************************************************************
!
subroutine hdf5_get_infile(idata, infile)
!
! Get the file name for the current data set
!
integer         , intent(in)  :: idata
character(len=*), intent(out) :: infile
!
infile = h5_temp%h5_infile
!
end subroutine hdf5_get_infile
!
!*******************************************************************************
!
subroutine hdf5_set_ku_is_h5(izz,ku_is_h5)
!
IMPLICIT NONE
!
integer, intent(in) :: izz
integer, intent(in) :: ku_is_h5
!
h5_ku_is_h5(izz) = ku_is_h5
!
END subroutine hdf5_set_ku_is_h5
!
!*******************************************************************************
!
subroutine hdf5_set_h5_is_ku(inumber, h5_is_ku)
!
IMPLICIT NONE
!
integer, intent(in) :: inumber
integer, intent(in) :: h5_is_ku
!
h5_h5_is_ku(inumber) = h5_is_ku
!
END subroutine hdf5_set_h5_is_ku
!
!*******************************************************************************
!
SUBROUTINE hdf5_set_map(dims, odata)
!
IMPLICIT NONE
!
INTEGER,            DIMENSION(3),                         INTENT(IN) :: dims
REAL(KIND=PREC_DP), DIMENSION(dims(1), dims(2), dims(3)), INTENT(IN) :: odata
!
INTEGER :: i,j,k
!
DO i=1, dims(1)
   DO j=1, dims(2)
      DO k=1, dims(3)
         h5_temp%h5_data(i,j,k)= odata(i,j,k) 
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE hdf5_set_map
!
!*******************************************************************************
!
SUBROUTINE hdf5_reset
!
TYPE(h5_data_struc), POINTER :: h5_current => NULL()
!
IF(ASSOCIATED(h5_root)) THEN       ! A storage does exist
   h5_temp => h5_root
   IF(ALLOCATED(h5_temp%h5_data)) DEALLOCATE(h5_temp%h5_data)
   find_node: DO 
      IF(ASSOCIATED(h5_temp%after)) THEN   ! A next node exists
         h5_current => h5_temp             ! Point to current
         h5_temp    => h5_temp%after       ! Point to next node
         DEALLOCATE(h5_current)            ! Clean up current node
      ELSE
         h5_current => h5_temp             ! Point to current
         DEALLOCATE(h5_current)            ! Clean up current node
         EXIT find_node                    ! We are done
      ENDIF
   ENDDO find_node
ENDIF
NULLIFY(h5_temp)
NULLIFY(h5_root)
h5_number   = 0
h5_h5_is_ku = 0
h5_ku_is_h5 = 0
!
END SUBROUTINE hdf5_reset
!
!*******************************************************************************
!
end MODULE lib_data_struc_h5
