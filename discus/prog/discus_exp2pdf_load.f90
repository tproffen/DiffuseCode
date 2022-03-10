MODULE exp2pdf_load_mod
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE exp2pdf_load(itype, line, length)
!
! Loads the Data set and/or the Sigmas
! Either explicitly or with reference to a KUPLOT data set
!
use kuplot_mod
use kuplot_load_mod
!
USE exp2pdf_data_mod
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!
integer         , INTENT(IN)    :: itype   ! Data==1 Background = 1 or SIGMA == 3
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
integer                              :: lstr
integer                              :: i
INTEGER                              :: ianz
INTEGER                              :: ndata 
!
integer, parameter :: NOPTIONAL = 1
integer, parameter :: O_BACK    = 1
character(LEN=   5), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'scale'   /
data loname /  5        /
opara  =  (/ '1.0000' /)   ! Always provide fresh default values
lopara =  (/  6       /)
owerte =  (/  1.0     /)
!
if(itype==2) then     ! For backgrounds only, check for "scale:"
   string = line
   lstr   = length
   CALL get_params(string, ianz, cpara, lpara, MAXW, lstr)
   if(ier_num /= 0) return
   do i=1, ianz
      if(cpara(i)(1:5)/='scale') then
         cpara(i) = ' '
         lpara(i) = 0
      endif
   enddo
   call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   if(lpresent(O_BACK)) then
      line = ' '
      CALL get_params(string, ianz, cpara, lpara, MAXW, lstr)
      do i=1, ianz
         if(cpara(i)(1:5)/='scale' .and. cpara(i)/=' ') then
            line = line(1:len_trim(line)) // cpara(i)(1:lpara(i))
            if(i<ianz) line = line(1:len_trim(line)) // ','
         endif
      enddo
      length = len_trim(line)
      if(line(length:length)==',') then
         line(length:length) = ' '
         length = length -1
      endif
   else 
      line = string
   endif
endif
!
if(ier_num /= 0) return
!
IF(line(1:6) == 'kuplot') THEN
   CALL get_params(line, ianz, cpara, lpara, MAXW, length)
   IF(ier_num/= 0) RETURN
   cpara(1) = '0'
   lpara(1) = 1
   CALL ber_params(ianz, cpara, lpara, werte, MAXW)
   IF(ier_num/= 0) RETURN
   ndata = NINT(werte(2))
   IF(itype==1) THEN             ! Data loaded from KUPLOT
      exp_load = ' '
      exp_kload = ndata
   elseIF(itype==2) THEN             ! Data loaded from KUPLOT
      exp_cback = ' '
      exp_kback = ndata
   ELSEif(itype==3) then
      exp_csigma = ' '        ! Sigma loaded from KUPLOT
      exp_ksigma = ndata
   ENDIF
ELSE                               ! Presume a "data xy, filename "
   IF(itype==1) THEN
      exp_load = line
      exp_kload = iz - 1
   elseif(itype==2) THEN
      exp_cback = line
      exp_kback = iz - 1
   elseif(itype==3) then
      exp_csigma = line
      exp_ksigma = iz - 1
   ENDIF
   CALL do_load(line, length,.TRUE.)
   IF(ier_num/= 0) RETURN
   ndata = -1                 ! Will be updated to correct value in refine_load_kuplot
ENDIF
CALL exp2pdf_load_kuplot(itype, ndata)
!
exp_kupl = MAX(exp_kupl, ndata)   ! This is the last KUPLOT data set that needs to be kept
!
if(itype == 2) exp_bscale = owerte(O_BACK)
!
END SUBROUTINE exp2pdf_load
!
!*******************************************************************************
!
SUBROUTINE exp2pdf_load_kuplot(itype, ndata) !, is_data, is_sigma)
!
! Transfers data set no. ndata from KUPLOT into REFINE memory
! If itype== 1     it is data
! If itype== 2     it is background
! If itype== 3     it is sigma
! If ndata==-1, the last KUPLOT data set is taken
!
use exp2pdf_data_mod
!
USE kuplot_mod
!
USE errlist_mod
USE define_variable_mod
use precision_mod
!
IMPLICIT NONE
!
integer, INTENT(IN)    :: itype   ! Data==1 Background = 1 or SIGMA == 3
INTEGER, INTENT(INOUT) :: ndata   ! no of data set to be transfered from KUPLOT
!
!LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE. ! Prevents user from deleting variables
INTEGER :: ix! , iy                      ! Dummy loop variables
!REAL(kind=PREC_DP)    :: step
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
IF(itype==1) THEN                         ! This is the data set
   IF(ALLOCATED(exp_data))   DEALLOCATE(exp_data)
   IF(ALLOCATED(exp_sigma )) DEALLOCATE(exp_sigma )
   IF(ALLOCATED(exp_x     )) DEALLOCATE(exp_x     )
   IF(ALLOCATED(exp_y     )) DEALLOCATE(exp_y     )
ENDIF
!
! 3D PDF to follow ... ... 
!IF(lni(ndata)) THEN                    ! 2D data set
!   IF(LDATA) THEN                      ! This is the data set
!      exp_dim(1) = nx(ndata )
!      exp_dim(2) = ny(ndata )
!      exp_dim(3) = ny(ndata )
!      ALLOCATE(exp_data  (exp_dim(1),exp_dim(2)))
!      ALLOCATE(exp_sigma (exp_dim(1),exp_dim(2)))
!      ALLOCATE(exp_x     (exp_dim(1)))
!      ALLOCATE(exp_y     (exp_dim(2)))
!!
!      DO iy=1,exp_dim(2)
!         DO ix=1,exp_dim(1)
!            exp_data(ix,iy)  = z (offz(ndata - 1) + (ix - 1)*ny(ndata) + iy)
!            exp_sigma (ix,iy) = 1.0000    ! dz(offxy(iz - 1) + ix) TEMPORARY unit weights
!         ENDDO
!         exp_y(iy)      = y(offxy(ndata - 1) + iy)
!      ENDDO
!      DO ix=1,exp_dim(1)
!         exp_x(ix)      = x(offxy(ndata - 1) + ix)
!      ENDDO
!      CALL def_set_variable('real', 'F_XMIN', exp_x(1),          IS_DIFFEV)
!      CALL def_set_variable('real', 'F_XMAX', exp_x(exp_dim(1)), IS_DIFFEV)
!      CALL def_set_variable('real', 'F_YMIN', exp_y(1),          IS_DIFFEV)
!      CALL def_set_variable('real', 'F_YMAX', exp_y(exp_dim(2)), IS_DIFFEV)
!      step = (exp_x(exp_dim(1))-exp_x(1))/FLOAT(exp_dim(1)-1)
!      CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
!      step = (exp_y(exp_dim(2))-exp_y(1))/FLOAT(exp_dim(2)-1)
!      CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
!   ELSE
!      IF(.NOT.ALLOCATED(exp_sigma )) THEN 
!         ier_num = -5
!         ier_typ = ER_APPL
!         RETURN
!      ENDIF
!!
!      IF(exp_dim(1) /= nx(ndata) .OR. exp_dim(2) /= ny(ndata)) THEN
!         ier_num = -6
!         ier_typ = ER_FORT
!         ier_msg(1) = 'SIGMA data set differs in size'
!         RETURN
!      ENDIF
!      DO iy=1,exp_dim(2)
!         DO ix=1,exp_dim(1)
!            exp_sigma (ix,iy) = z (offz(ndata - 1) + (ix - 1)*ny(ndata) + iy)
!         ENDDO
!      ENDDO
!   ENDIF
!ELSE                                     ! 1D data set
   IF(itype==1) THEN                      ! This is the data set
      exp_dim(1) = lenc(ndata )
      exp_dim(2) = 1
      exp_dim(3) = 1
      ALLOCATE(exp_data  (exp_dim(1),exp_dim(2), exp_dim(3)))
      ALLOCATE(exp_sigma (exp_dim(1),exp_dim(2), exp_dim(3)))
      ALLOCATE(exp_x     (exp_dim(1)))
      ALLOCATE(exp_y     (1         ))
      DO ix=1,exp_dim(1)
         exp_data(ix,1,1)   = y(offxy(ndata - 1) + ix)
         exp_sigma (ix,1,1) = ABS(dy(offxy(ndata - 1) + ix))
         exp_x(ix)          = x(offxy(ndata - 1) + ix)
      ENDDO
      exp_kload = ndata 
      exp_y(1) = 1.0
!     CALL def_set_variable('real', 'F_XMIN', exp_x(1),          IS_DIFFEV)
!     CALL def_set_variable('real', 'F_XMAX', exp_x(exp_dim(1)), IS_DIFFEV)
!     CALL def_set_variable('real', 'F_YMIN', exp_y(1),          IS_DIFFEV)
!     CALL def_set_variable('real', 'F_YMAX', exp_y(exp_dim(2)), IS_DIFFEV)
!     step = (exp_x(exp_dim(1))-exp_x(1))/FLOAT(exp_dim(1)-1)
!     CALL def_set_variable('real', 'F_XSTP', step             , IS_DIFFEV)
!     step = 1.0
!     CALL def_set_variable('real', 'F_YSTP', step             , IS_DIFFEV)
   ELSE                      ! Background or sigma
      IF(.NOT.ALLOCATED(exp_sigma )) THEN 
         ier_num = -5
         ier_typ = ER_APPL
         RETURN
      ENDIF
!
      if(itype==2) then
         exp_dim_b(1) = lenc(ndata )
         exp_dim_b(2) = 1
         exp_dim_b(3) = 1
         IF(ALLOCATED(exp_back))   DEALLOCATE(exp_back)
         IF(ALLOCATED(exp_sigmab)) DEALLOCATE(exp_sigmab)
         IF(ALLOCATED(exp_xb  ))   DEALLOCATE(exp_xb  )
         ALLOCATE(exp_back  (exp_dim_b(1),exp_dim_b(2), exp_dim_b(3)))
         ALLOCATE(exp_sigmab(exp_dim_b(1),exp_dim_b(2), exp_dim_b(3)))
         ALLOCATE(exp_xb    (exp_dim_b(1)))
         DO ix=1,exp_dim_b(1)
            exp_back(ix,1,1) = y(offxy(ndata - 1) + ix)
            exp_sigmab(ix,1,1) = ABS(dy(offxy(ndata - 1) + ix))
            exp_xb( ix)      = x(offxy(ndata - 1) + ix)
         ENDDO
         exp_kback = ndata
      elseif(itype==3) then
!
         IF(exp_dim(1) /= nx(ndata)) THEN
            ier_num = -6
            ier_typ = ER_FORT
            ier_msg(1) = 'SIGMA/Background data set differs in size'
            RETURN
         ENDIF
         DO ix=1,exp_dim(1)
            exp_sigma (ix,1,1) = y(offxy(ndata - 1) + ix)
         ENDDO
         IF(MINVAL(exp_sigma (:,1, 1))==0.0) THEN
            ier_num = -7
            ier_typ = ER_APPL
            ier_msg(1) = ' Check data and define non-zero sigma'
            RETURN
         ENDIF
      endif
   ENDIF
!ENDIF
!
! Scratch all KUPLOT data beyond current data set
!
!IF(LDATA) THEN
!   CALL def_set_variable('integer', 'F_DATA', real(ndata,kind=PREC_DP),  IS_DIFFEV)
!   CALL def_set_variable('integer', 'F_SIGMA', real(ndata, kind=PREC_DP), IS_DIFFEV)
!ELSE
!   CALL def_set_variable('integer', 'F_SIGMA', real(ndata, kind=PREC_DP), IS_DIFFEV)
!ENDIF
!iz = ndata + 1
!
!
END SUBROUTINE exp2pdf_load_kuplot
!
!*******************************************************************************
!
end MODULE exp2pdf_load_mod
