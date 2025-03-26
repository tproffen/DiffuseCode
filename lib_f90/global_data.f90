MODULE global_data_mod
!-
!   Contains a global data set pool, needed by REFINE and KUPLOT
!+
!
USE precision_mod
!
IMPLICIT NONE
!
PRIVATE
PUBLIC   gl_alloc           ! Allocate the global array to (dim1,dim2, dim3, NPARA)
PUBLIC   gl_set_pnumber     ! Populate the lookup for function 'pnumber()'
PUBLIC   gl_get_pnumber     ! Evaluate function p_number
PUBLIC   gl_set_use         ! Turn us of global data on
PUBLIC   gl_is_der          ! Test if derivative was set
PUBLIC   gl_set_npara       ! Set number of refine parameters
PUBLIC   gl_get_npara       ! Get number of refine parameters
public   gl_get_maxdims     ! Get maximum dimensions in gl_data
public   gl_set_x           ! Set x-coordinates
public   gl_set_y           ! Set x-coordinates
public   gl_set_z           ! Set x-coordinates
public   gl_get_x           ! Set x-coordinates
public   gl_get_y           ! Set x-coordinates
public   gl_get_z           ! Set x-coordinates
!
INTERFACE gl_set_data
   MODULE PROCEDURE gl_set_data_1d, gl_set_data_2d, gl_set_data_3d
END INTERFACE
!
INTERFACE gl_get_data
   MODULE PROCEDURE gl_get_data_1d, gl_get_data_2d, gl_get_data_3d
END INTERFACE
!
PUBLIC gl_set_data
PUBLIC gl_get_data
!
CHARACTER(LEN=1024), DIMENSION(:),       ALLOCATABLE :: gl_params
REAL(KIND=PREC_DP),  DIMENSION(:,:,:,:), ALLOCATABLE :: gl_data ! data sets
REAL(KIND=PREC_DP),  DIMENSION(:)      , ALLOCATABLE :: gl_x    ! data sets
REAL(KIND=PREC_DP),  DIMENSION(:)      , ALLOCATABLE :: gl_y    ! data sets
REAL(KIND=PREC_DP),  DIMENSION(:)      , ALLOCATABLE :: gl_z    ! data sets
INTEGER           ,  DIMENSION(4)                    :: gl_dims =(/0,0,0,0/)  ! Dimensions
LOGICAL           ,  DIMENSION(:),       ALLOCATABLE :: gl_lderiv
LOGICAL                                              :: gl_use=.FALSE.        ! Use if allowed by refine
INTEGER                                              :: gl_npara=0            ! Number or refined parameters
INTEGER                                              :: gl_nfix =0            ! Number or fixed   parameters
!
!*******************************************************************************
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE gl_alloc(in_dims)
!-
!   Allocate the global data sets
!   Data dimensions are in gl_data(1,2,3, )
!   gl_data(:,:,:,-3)    ! Experimental sigma
!   gl_data(:,:,:,-2)    ! Experimental data
!   gl_data(:,:,:,-1)    ! Data from refinement current optimum
!   gl_data(:,:,:, 0)    ! Data from Cost function current values, may depend on derivative
!   gl_data(:,:,:, 1:)   ! Data from Cost function derivatives 1, 2, 3, ...
!   
!+
IMPLICIT NONE
!
INTEGER, DIMENSION(4), INTENT(IN) :: in_dims
!
IF(ALLOCATED(gl_data))   DEALLOCATE(gl_data)
IF(ALLOCATED(gl_lderiv)) DEALLOCATE(gl_lderiv)
IF(ALLOCATED(gl_x     )) DEALLOCATE(gl_x)
IF(ALLOCATED(gl_y     )) DEALLOCATE(gl_y)
IF(ALLOCATED(gl_z     )) DEALLOCATE(gl_z)
!
ALLOCATE(gl_data(in_dims(1), in_dims(2), in_dims(3), -3:in_dims(4)))
ALLOCATE(gl_lderiv(-1:in_dims(4)))
allocate(gl_x(in_dims(1)))
allocate(gl_y(in_dims(2)))
allocate(gl_z(in_dims(3)))
!
gl_dims = in_dims
gl_data = 0.0D0
gl_lderiv = .FALSE.
gl_x = 0.0D0
gl_y = 0.0D0
gl_z = 0.0D0
!
END SUBROUTINE gl_alloc
!
!*******************************************************************************
!
SUBROUTINE gl_set_pnumber(MAXPAR, MAXFIX, NPAR, NFIX, refine_params, refine_fixed)
!-
!   Populate the global values for pnumber() function, called by refine
!+
IMPLICIT NONE
!
INTEGER , INTENT(IN) :: MAXPAR
INTEGER , INTENT(IN) :: MAXFIX
INTEGER , INTENT(IN) :: NPAR
INTEGER , INTENT(IN) :: NFIX
CHARACTER(LEN=*), DIMENSION(MAXPAR), INTENT(IN) :: refine_params
CHARACTER(LEN=*), DIMENSION(MAXFIX), INTENT(IN) :: refine_fixed
!
INTEGER :: i
!
IF(ALLOCATED(gl_params)) DEALLOCATE(gl_params)
!
ALLOCATE(gl_params(-ABS(MAXFIX):ABS(MAXPAR)))
gl_params    = ' '
DO i=1, NFIX
   gl_params(-i) = refine_fixed(i)
ENDDO
DO i=1, NPAR
   gl_params( i) = refine_params(i)
ENDDO
gl_npara = NPAR
gl_nfix  = NFIX
!
END SUBROUTINE gl_set_pnumber
!
!*******************************************************************************
!
SUBROUTINE gl_set_use(luse)
!-
!   Sets the use of global data set on/off
!-
LOGICAL, INTENT(IN) :: luse
!
gl_use = luse
!
END SUBROUTINE gl_set_use
!
!*******************************************************************************
!
SUBROUTINE gl_set_npara(npara)
!-
!   Sets the use of global data set on/off
!-
INTEGER, INTENT(IN) :: npara
!
gl_npara = npara
!
END SUBROUTINE gl_set_npara
!
!*******************************************************************************
!
LOGICAL FUNCTION gl_is_der(i)
!-
!   Tests if a derivative number i was set
!-
!
use errlist_mod
INTEGER, INTENT(IN ) :: i
!
gl_is_der = .FALSE.
if(allocated(gl_lderiv)) then
   if(i>lbound(gl_lderiv,1) .and. i<=ubound(gl_lderiv,1)) then
      gl_is_der = gl_lderiv(i)
   else
      ier_num = -6
      ier_typ = ER_COMM
      write(ier_msg(1),'(a,i4)') 'Derivative not stored in global ', i
   endif
else
   ier_num = -6
   ier_typ = ER_COMM
   ier_msg(1) = 'Global data not allocated'
endif
!
END FUNCTION gl_is_der
!
!*******************************************************************************
!
SUBROUTINE gl_get_npara(npara, nfix)
!-
!   Sets the use of global data set on/off
!-
INTEGER, INTENT(OUT) :: npara
INTEGER, INTENT(OUT) :: nfix
!
npara = gl_npara
nfix  = gl_nfix
!
END SUBROUTINE gl_get_npara
!
!*******************************************************************************
!
INTEGER FUNCTION gl_get_pnumber(string)  RESULT(iww)
!-
!   Get the entry number in function pnumber('par_name')
!+
USE errlist_mod
USE build_name_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
!
CHARACTER(LEN=*), INTENT(INOUT) :: string
!
INTEGER, PARAMETER :: MAXW=2
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))), DIMENSION(MAXW) :: cpara
INTEGER                                    , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP),                          DIMENSION(MAXW) :: werte
LOGICAL  :: lstring   ! pnumber argument is a string
LOGICAL  :: lfound    ! pnumber entry was found
INTEGER  :: i
INTEGER  :: length
INTEGER  :: ianz
!
iww = 0
!
IF(.NOT. ALLOCATED(gl_params)) THEN
   ier_num = -52
   ier_typ = ER_FORT
   ier_msg(1) = 'Check if any ''newparam'' command was used in'
   ier_msg(2) = 'the REFINE section'
   RETURN
ENDIF
!
lstring = .FALSE.
IF(string(1:1)=='"') THEN
   length = LEN_TRIM(string)
   CALL get_params(string, ianz, cpara, lpara, MAXW, length)
   CALL do_build_name (ianz, cpara, lpara, werte, MAXW, 1)
   lstring = .TRUE.
ELSEIF(string(1:1) == '''') THEN
   length = LEN_TRIM(string)-2
   CALL get_params(string(2:LEN_TRIM(string)-1), ianz, cpara, lpara, MAXW, length)
!         IF(ianz==1) THEN
!         ENDIF
   lstring = .TRUE.
ELSE
   CALL get_params(string, ianz, cpara, lpara, MAXW, length)
   CALL ber_param(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/=0) RETURN
   lstring = .FALSE.
ENDIF
IF(lstring) THEN      ! parameter is a string, look for refined parameter name
   lfound = .FALSE.
   searchp: DO i=LBOUND(gl_params,1), UBOUND(gl_params,1)
      IF(cpara(1) ==     gl_params(i)) THEN
         iww = i
         lfound=.TRUE.
         EXIT searchp
      ENDIF
   ENDDO searchp
   IF(.NOT.lfound) THEN
      ier_num = -51
      ier_typ = ER_FORT
      ier_msg(1) = 'Check ''newparam'' commands in REFINE'
      RETURN
   ENDIF
ELSE
   iww = NINT(werte(1))
   IF(iww>=LBOUND(gl_params,1) .AND. iww<=UBOUND(gl_params,1)) THEN
      ier_num = 0
      ier_typ = ER_NONE
   ELSE
      ier_num = -50
      ier_typ = ER_FORT
      ier_msg(1) = 'Check ''newparam'' commands in REFINE'
      ier_msg(2) = 'Check range of refined/fixed parameters'
   ENDIF
ENDIF
!
END FUNCTION gl_get_pnumber
!
!*******************************************************************************
!
subroutine gl_get_maxdims(dims)
!
! Read the maximum dimensions of the global data array
!
integer, dimension(4), intent(out) :: dims
!
dims(1) = ubound(gl_data, 1)
dims(2) = ubound(gl_data, 2)
dims(3) = ubound(gl_data, 3)
dims(4) = ubound(gl_data, 4)
!
end subroutine gl_get_maxdims
!
!*******************************************************************************
!
SUBROUTINE gl_set_data_1d(idim1, ipara, ext_data)
!
! popluate global data, 1D 
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN) :: idim1
!INTEGER                            , INTENT(IN) :: NPARA
INTEGER                            , INTENT(IN) :: ipara
REAL(KIND=PREC_DP), DIMENSION(idim1), INTENT(IN) :: ext_data
!
!
gl_data(1:idim1,1, 1, ipara) = ext_data(1:idim1)
if(ipara>=0) gl_lderiv(ipara) = .TRUE.
!
END SUBROUTINE gl_set_data_1d
!
!*******************************************************************************
!
SUBROUTINE gl_set_data_2d(idim1, idim2, ipara, ext_data)
!
! popluate global data, 2D 
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN) :: idim1
INTEGER                            , INTENT(IN) :: idim2
!INTEGER                            , INTENT(IN) :: NPARA
INTEGER                            , INTENT(IN) :: ipara
REAL(KIND=PREC_DP), DIMENSION(idim1, idim2), INTENT(IN) :: ext_data
!
!
gl_data(1:idim1,1:idim2, 1, ipara) = ext_data(1:idim1, 1:idim2)
if(ipara>=0) gl_lderiv(ipara) = .TRUE.
!
END SUBROUTINE gl_set_data_2d
!
!*******************************************************************************
!
SUBROUTINE gl_set_data_3d(idim1, idim2, idim3, ipara, ext_data)
!
! popluate global data, 3D 
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN) :: idim1
INTEGER                            , INTENT(IN) :: idim2
INTEGER                            , INTENT(IN) :: idim3
!INTEGER                            , INTENT(IN) :: NPARA
INTEGER                            , INTENT(IN) :: ipara
REAL(KIND=PREC_DP), DIMENSION(idim1, idim2, idim3), INTENT(IN) :: ext_data
!
!
gl_data(1:idim1,1:idim2, 1:idim3, ipara) = ext_data(1:idim1, 1:idim2, 1:idim3)
if(ipara>=0) gl_lderiv(ipara) = .TRUE.
!
END SUBROUTINE gl_set_data_3d
!
!*******************************************************************************
!
subroutine gl_set_x(idimen, ext_x)
!-
!  Set the x-values of the global data array
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: idimen
real(kind=PREC_DP), dimension(1:idimen), intent(in) :: ext_x
!
gl_x(1:idimen) = ext_x(1:idimen)
!
end subroutine gl_set_x
!
!*******************************************************************************
!
subroutine gl_set_y(idimen, ext_y)
!-
!  Set the x-values of the global data array
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: idimen
real(kind=PREC_DP), dimension(1:idimen), intent(in) :: ext_y
!
gl_y(1:idimen) = ext_y(1:idimen)
!
end subroutine gl_set_y
!
!*******************************************************************************
!
subroutine gl_set_z(idimen, ext_z)
!-
!  Set the x-values of the global data array
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: idimen
real(kind=PREC_DP), dimension(1:idimen), intent(in) :: ext_z
!
gl_z(1:idimen) = ext_z(1:idimen)
!
end subroutine gl_set_z
!
!*******************************************************************************
!
SUBROUTINE gl_get_data_1d(ipara, idim1, ext_data)
!
! retrieve global data, 1D 
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                             , INTENT(IN)  :: ipara
INTEGER                             , INTENT(IN)  :: idim1
REAL(KIND=PREC_DP), DIMENSION(idim1), INTENT(OUT) :: ext_data
!
ext_data(1:idim1) = gl_data(1:idim1,1, 1, ipara)
!
END SUBROUTINE gl_get_data_1d
!
!*******************************************************************************
!
SUBROUTINE gl_get_data_2d(ipara, idim1, idim2, ext_data)
!
! retrieve global data, 2D 
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                             , INTENT(IN)  :: ipara
INTEGER                             , INTENT(IN)  :: idim1
INTEGER                             , INTENT(IN)  :: idim2
REAL(KIND=PREC_DP), DIMENSION(idim1,idim2), INTENT(OUT) :: ext_data
!
ext_data(1:idim1, 1:idim2) = gl_data(1:idim1,1:idim2, 1, ipara)
!
END SUBROUTINE gl_get_data_2d
!
!*******************************************************************************
!
SUBROUTINE gl_get_data_3d(ipara, idim1, idim2, idim3, ext_data)
!
! retrieve global data, 3D 
USE errlist_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)  :: ipara
INTEGER                            , INTENT(IN)  :: idim1
INTEGER                            , INTENT(IN)  :: idim2
INTEGER                            , INTENT(IN)  :: idim3
REAL(KIND=PREC_DP), DIMENSION(idim1, idim2, idim3), INTENT(OUT) :: ext_data
!
ext_data(1:idim1,1:idim2,1:idim3) = gl_data(1:idim1,1:idim2, 1:idim3, ipara)
!
END SUBROUTINE gl_get_data_3d
!
!*******************************************************************************
!
subroutine gl_get_x(idimen, ext_x)
!-
!  Set the x-values of the global data array
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: idimen
real(kind=PREC_DP), dimension(1:idimen), intent(out) :: ext_x
!
ext_x(1:idimen) = gl_x(1:idimen)
!
end subroutine gl_get_x
!
!*******************************************************************************
!
subroutine gl_get_y(idimen, ext_y)
!-
!  Set the x-values of the global data array
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: idimen
real(kind=PREC_DP), dimension(1:idimen), intent(out) :: ext_y
!
ext_y(1:idimen) = gl_y(1:idimen)
!
end subroutine gl_get_y
!
!*******************************************************************************
!
subroutine gl_get_z(idimen, ext_z)
!-
!  Set the z-values of the global data array
!+
!
use precision_mod
!
implicit none
!
integer, intent(in) :: idimen
real(kind=PREC_DP), dimension(1:idimen), intent(out) :: ext_z
!
ext_z(1:idimen) = gl_z(1:idimen)
!
end subroutine gl_get_z
!
!*******************************************************************************
!
END MODULE global_data_mod
