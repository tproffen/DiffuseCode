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
USE refine_params_mod
USE ber_params_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: cpara
INTEGER            , DIMENSION(:), ALLOCATABLE :: lpara
REAL               , DIMENSION(:), ALLOCATABLE :: werte
INTEGER :: k
INTEGER :: ianz
INTEGER :: MAXW
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
CALL show_fit_erg(REF_MAXPARAM, REF_MAXPARAM_FIX, refine_par_n, refine_fix_n, &
           ref_dim(1)*ref_dim(2),      &
           refine_chisqr, refine_conf, &
           refine_lamda, refine_rval, refine_rexp,                        &
           refine_params, refine_p, refine_dp, refine_range, refine_cl,   &
           refine_fixed, refine_f)
!
END SUBROUTINE refine_do_show
!
!*******************************************************************************
!
SUBROUTINE show_fit_erg(MAXP, MAXF, npara, nfixed, ndata, chisq, conf, lamda,  &
           r4, re, &
           params, pp, dpp, prange, cl, fixed, pf)
!+                                                                      
!     Display fit results
!-                                                                      
!
USE prompt_mod
!                                                                       
IMPLICIT none 
!
INTEGER, INTENT(IN) :: MAXP
INTEGER, INTENT(IN) :: MAXF
INTEGER, INTENT(IN) :: npara
INTEGER, INTENT(IN) :: nfixed
INTEGER, INTENT(IN) :: ndata
REAL, INTENT(IN) :: chisq
REAL, INTENT(IN) :: conf
REAL, INTENT(IN) :: lamda
REAL, INTENT(IN) :: r4
REAL, INTENT(IN) :: re
CHARACTER(LEN=*), DIMENSION(MAXP), INTENT(IN) :: params
REAL            , DIMENSION(MAXP), INTENT(IN) :: pp
REAL            , DIMENSION(MAXP), INTENT(IN) :: dpp
REAL            , DIMENSION(MAXP,2), INTENT(IN) :: prange
REAL            , DIMENSION(NPARA,NPARA), INTENT(IN) :: cl
CHARACTER(LEN=*), DIMENSION(MAXF), INTENT(IN) :: fixed
REAL            , DIMENSION(MAXF), INTENT(IN) :: pf
!                                                                       
INTEGER :: i, j 
LOGICAL :: kor 
!
IF(chisq>=0.0) THEN
   WRITE(output_io, 1040) chisq, chisq/ndata, conf, chisq/(ndata-npara),  &
                          ndata, npara, lamda, r4, re 
   WRITE(output_io, 1050) 
   kor = .false. 
   DO i = 2, npara 
      DO j = 1, i - 1 
         IF(ABS(cl(i, j)) >=  0.8) THEN 
            WRITE(output_io, 1070) i, j, cl(i, j) 
            kor = .TRUE. 
         ENDIF 
      ENDDO 
   ENDDO 
   IF (.NOT.kor) WRITE(output_io, 1060) 
ENDIF
WRITE(output_io, * ) ' ' 
!
DO i=1, npara
   IF(prange(i,1)<=prange(i,2)) THEN
      WRITE(output_io,1110) params(i), pp(i), dpp(i), prange(i,1), prange(i,2)
   ELSE
      WRITE(output_io,1100) params(i), pp(i), dpp(i)
   ENDIF
ENDDO
DO i=1, nfixed
   WRITE(output_io,1200) fixed(i), pf(i)
ENDDO
WRITE(output_io, * ) ' ' 
!                                                                       
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
 1100 FORMAT (3x,a16,' : ',g15.6E3, ' +- ',g15.6,4x)
 1110 FORMAT (3x,a16,' : ',g15.6E3, ' +- ',g15.6,4x, '[',g15.6E3,',',g15.6E3,']')
 1200 FORMAT (3x,a16,' : ',g15.6E3)
!                                                                       
END SUBROUTINE show_fit_erg                   
!
!*******************************************************************************
!
END MODULE refine_show_mod
