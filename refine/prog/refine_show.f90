MODULE refine_show_mod
!
IMPLICIT NONE
!
PUBLIC show_fit_erg
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_do_show(line, length)
!
USE refine_fit_erg
USE refine_data_mod
USE refine_mac_mod
USE refine_params_mod
!
USE errlist_mod
USE get_params_mod
USE ber_params_mod
USE build_name_mod
USE precision_mod
USE prompt_mod
USE string_convert_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXWW = 2       ! Max parameters on line
CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: cpara
INTEGER            , DIMENSION(:), ALLOCATABLE :: lpara
REAL(KIND=PREC_DP) , DIMENSION(:), ALLOCATABLE :: werte
INTEGER :: k
INTEGER :: ianz
INTEGER :: MAXW
INTEGER :: iounit
!
LOGICAL :: lcovar
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: OCOVAR    = 1
INTEGER, PARAMETER :: OFILE     = 2
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'covar'  , 'output'/
DATA loname /  5       ,  6      /
opara  =  (/ 'no    ', 'screen' /)   ! Always provide fresh default values
lopara =  (/  6      ,  6         /)
owerte =  (/   0.0   ,   0.0      /)
!
iounit = output_io
!
! Get and interpret (optional) parameters
!
ALLOCATE(cpara(1:MAXWW))
ALLOCATE(lpara(1:MAXWW))
ALLOCATE(werte(1:MAXWW))
!
CALL get_params(line, ianz, cpara, lpara, MAXWW, length)
CALL get_optional(ianz, MAXWW, cpara, lpara, NOPTIONAL, ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
CALL do_low(opara(OCOVAR))
lcovar = opara(OCOVAR)=='yes'     ! Covariance yes/no
!
IF(opara(OFILE)/='screen') THEN   ! Output to screen or to file
   IF(opara(OFILE)(1:1)=='[' .AND. opara(OFILE)(lopara(OFILE):lopara(OFILE))==']') THEN
      k = lopara(OFILE)-2
      CALL get_params(opara(OFILE)(2:lopara(OFILE)-1), ianz, cpara, lpara, MAXWW, k)
      IF(ier_num==0) THEN
         CALL do_build_name(ianz, cpara, lpara, werte, MAXWW, 1)
      ENDIF
   ELSE
      cpara(1) = opara(OFILE)
   ENDIF
   IF(ier_num==0) THEN            ! Open output unit
      iounit = 77
      CALL oeffne(iounit, cpara(1), 'unknown')
   ENDIF
ELSE
   iounit = output_io             ! Default is output into "output_io"
ENDIF
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
IF(ier_num/=0) RETURN
!
!  Caclulate current parameter values
!
ALLOCATE(cpara(1:refine_par_n))
ALLOCATE(lpara(1:refine_par_n))
ALLOCATE(werte(1:refine_par_n))
DO k=1, refine_par_n
   cpara(k) = refine_params(k)
   lpara(k) = LEN_TRIM(refine_params(k))
ENDDO
ianz = refine_par_n
MAXW = refine_par_n
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO k=1, refine_par_n
   refine_p(k) = werte(k)
ENDDO
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
!
!  Caclulate current fixed parameter values
!
ALLOCATE(cpara(1:refine_fix_n))
ALLOCATE(lpara(1:refine_fix_n))
ALLOCATE(werte(1:refine_fix_n))
DO k=1, refine_fix_n
   cpara(k) = refine_fixed(k)
   lpara(k) = LEN_TRIM(refine_fixed(k))
ENDDO
ianz = refine_fix_n
MAXW = refine_fix_n
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO k=1, refine_fix_n
   refine_f(k) = werte(k)
ENDDO
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
!
IF(.NOT.ALLOCATED(refine_CL)) THEN
   ALLOCATE(refine_cl(refine_par_n, refine_par_n))
ENDIF
!
CALL show_fit_erg(iounit, REF_MAXPARAM, REF_MAXPARAM_FIX, refine_par_n,   &
           refine_fix_n, &
           ref_dim(1)*ref_dim(2),      &
           refine_mac, refine_mac_l,   &
           ref_load, ref_kload, ref_sigma, ref_ksigma, lcovar,            &
           refine_chisqr, refine_conf, &
           refine_lamda, refine_rval, refine_rexp,                        &
           refine_params, refine_p, refine_dp, refine_range, refine_cl,   &
           refine_fixed, refine_f)
!
IF(iounit/=output_io) CLOSE(iounit)
!
END SUBROUTINE refine_do_show
!
!*******************************************************************************
!
SUBROUTINE show_fit_erg(iounit, MAXP, MAXF, npara, nfixed, ndata, mac, mac_l,   &
           load, kload, sigma, ksigma, lcovar, chisq, conf, lamda,  &
           r4, re, &
           params, pp, dpp, prange, cl, fixed, pf)
!+                                                                      
!     Display fit results
!-                                                                      
!
!                                                                       
IMPLICIT none 
!
INTEGER, INTENT(IN) :: iounit         ! Output screen / file
INTEGER, INTENT(IN) :: MAXP           ! Parameter array size
INTEGER, INTENT(IN) :: MAXF           ! Fixed parameter array size
INTEGER, INTENT(IN) :: npara          ! Number of refined parameters
INTEGER, INTENT(IN) :: nfixed         ! Number of fixed parameters
INTEGER, INTENT(IN) :: ndata          ! Number of data points = ref_dim(1)*ref_dim(2)
CHARACTER(LEN=*), INTENT(IN) :: mac   ! Refinement macro name
INTEGER, INTENT(IN)          :: mac_l ! Length of macro name
CHARACTER(LEN=*), INTENT(IN) :: load  ! Data set loaded as:
INTEGER, INTENT(IN)          :: kload ! Datas set in KUPLOT number
CHARACTER(LEN=*), INTENT(IN) :: sigma ! Data set loaded as:
INTEGER, INTENT(IN)          :: ksigma! Datas set in KUPLOT number
LOGICAL, INTENT(IN)          :: lcovar! Display full covariance matrix
REAL, INTENT(IN) :: chisq             ! Chi^2
REAL, INTENT(IN) :: conf              ! Confidence level
REAL, INTENT(IN) :: lamda             ! Levenberg-Marquardt Lamda
REAL, INTENT(IN) :: r4                ! Weighted R-value
REAL, INTENT(IN) :: re                ! expected R-value
CHARACTER(LEN=*), DIMENSION(MAXP), INTENT(IN) :: params    ! Parameter names
REAL            , DIMENSION(MAXP), INTENT(IN) :: pp        ! Parameter values
REAL            , DIMENSION(MAXP), INTENT(IN) :: dpp       ! Parameter sigmas
REAL            , DIMENSION(MAXP,2), INTENT(IN) :: prange  ! Parameter range
REAL            , DIMENSION(NPARA,NPARA), INTENT(IN) :: cl ! Covariance matrix
CHARACTER(LEN=*), DIMENSION(MAXF), INTENT(IN) :: fixed     ! Fixed parameter names
REAL            , DIMENSION(MAXF), INTENT(IN) :: pf        ! Fixed parameter values
!                                                                       
INTEGER :: i, j 
LOGICAL :: kor 
!
IF(load==' ') THEN
   IF(kload>0) THEN
      WRITE(iounit, 3000) kload
   ELSE
      WRITE(iounit, 3010)
   ENDIF
ELSE
   WRITE(iounit, 3020) load(1:LEN_TRIM(load))
ENDIF
!
IF(sigma==' ') THEN
   IF(ksigma>0) THEN
      WRITE(iounit, 3100) ksigma
   ELSE
      WRITE(iounit, 3110)
   ENDIF
ELSE
   WRITE(iounit, 3120) sigma(1:LEN_TRIM(sigma))
ENDIF
!
WRITE(iounit, 2000) mac(1:mac_l)
!
IF(chisq>=0.0) THEN
   WRITE(iounit, 1040) chisq, chisq/ndata, conf, chisq/(ndata-npara),  &
                          ndata, npara, lamda, r4, re 
   WRITE(iounit, 1050) 
   kor = .false. 
   DO i = 2, npara 
      DO j = 1, i - 1 
         IF(ABS(cl(i, j)) >=  0.8) THEN 
            WRITE(iounit, 1070) i, j, cl(i, j) 
            kor = .TRUE. 
         ENDIF 
      ENDDO 
   ENDDO 
   IF (.NOT.kor) WRITE(iounit, 1060) 
   IF(lcovar) THEN
      WRITE(iounit,'(a)') ' Correlation matrix' 
      DO i=1, npara
         WRITE(iounit,5000) params(i),(cl(i,j),j=1,i)
      ENDDO
   ENDIF
ENDIF
WRITE(iounit, * ) ' ' 
!
WRITE(iounit,'(a)') ' Refined parameters' 
DO i=1, npara
   IF(prange(i,1)<=prange(i,2)) THEN
      WRITE(iounit,1110) params(i), pp(i), dpp(i), prange(i,1), prange(i,2)
   ELSE
      WRITE(iounit,1100) params(i), pp(i), dpp(i)
   ENDIF
ENDDO
IF(nfixed>0) THEN
   WRITE(iounit,'(a)') ' Fixed   parameters' 
   DO i=1, nfixed
      WRITE(iounit,1200) fixed(i), pf(i)
   ENDDO
ENDIF
WRITE(iounit, * ) ' ' 
!
3000 FORMAT(' Data in KUPLOT no. : ',i3)
3010 FORMAT(' Data not defined   '   )
3020 FORMAT(' Data loaded as     : ',a )
3100 FORMAT(' Sigma in KUPLOT no.: ',i3)
3110 FORMAT(' Sigma not defined  '   )
3120 FORMAT(' Sigma loaded as    : ',a )
2000 FORMAT(' Refinement macro   : ',a)
 1040 FORMAT (/,                                                        &
              ' Information about the fit : ',/,                        &
              3x,'Chi^2      : ',g12.6, 5x,' Chi^2/N : ',g12.6 /,       &
              3x,'Conf. level: ',g12.6, 1x,' Chi^2/(N-P) : ',g12.6 /,   &
              3x,'No.Data    : ',i12  , 4x,' No.Params: ',i12/,          &
              3x,'MRQ final  : ',g12.6,/,                               &
              3x,'wR value   : ',g12.6, 5x,' R exp   : ',g12.6/)        
 1050 FORMAT (' Correlations larger than 0.8 :') 
 1060 FORMAT (3x,'** none **') 
 1070 FORMAT (3x,'Between p(',i2,') - p(',i2,') : ',f6.3) 
 5000 FORMAT (3x,a16,' : ',10g11.3e2)
 1100 FORMAT (3x,a16,' : ',g15.6E3, ' +- ',g15.6,4x)
 1110 FORMAT (3x,a16,' : ',g15.6E3, ' +- ',g15.6,4x, '[',g15.6E3,',',g15.6E3,']')
 1200 FORMAT (3x,a16,' : ',g15.6E3)
!                                                                       
END SUBROUTINE show_fit_erg                   
!
!*******************************************************************************
!
END MODULE refine_show_mod
