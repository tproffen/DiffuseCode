MODULE diffev_show_mod
!
CONTAINS
!*****7*****************************************************************
!                                                                       
!     This routine shows all settings of DIFFEV parameters. Some        
!     parameter settings can also be displayed by entering the          
!     corresponding command without parameters.                         
!                                                                       
!*****7*****************************************************************
SUBROUTINE diffev_do_show (line, lp) 
!                                                                       
!     Main show menu                                                    
!                                                                       
USE diffev_allocate_appl
!
USE errlist_mod 
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER   :: maxw = 2
!                                                                       
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
INTEGER             , INTENT(INOUT) :: lp 
!
CHARACTER (LEN=1024), DIMENSION(maxw)  :: cpara (maxw) 
INTEGER             , DIMENSION(maxw)  :: lpara (maxw)
INTEGER                                :: ianz
LOGICAL                                :: str_comp 
!
!                                                                       
CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz>=  1) then 
   IF (str_comp (cpara (1) , 'config', 1, lpara (1) , 6) )  THEN
      CALL diffev_show_config
   ELSEIF (str_comp (cpara (1) , 'param', 3, lpara (1) , 5) )  THEN
      CALL diffev_show_param(ianz, cpara, lpara, MAXW)
   ELSEIF (str_comp (cpara (1) , 'population', 3, lpara (1) , 10) )  THEN
      CALL diffev_show_population(ianz, cpara, lpara, MAXW)
   ELSE
!                                                                       
!     -- try generic show commands                                      
!                                                                       
      CALL do_show_generic (cpara, lpara, maxw) 
   ENDIF
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE diffev_do_show
!
!*****7*****************************************************************3*******
!
SUBROUTINE diffev_show_param(ianz, cpara, lpara, MAXW)
!
USE population
!
USE errlist_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: ianz
INTEGER                            , INTENT(IN)    :: MAXW
CHARACTER(LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
INTEGER,          DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
!
INTEGER            , PARAMETER          :: MAXWW=1
CHARACTER(LEN=1024), DIMENSION(1:MAXWW) :: ccpara
INTEGER,             DIMENSION(1:MAXWW) :: llpara
REAL   ,             DIMENSION(1:MAXWW) :: wwerte
INTEGER :: i, j
INTEGER :: iianz
LOGICAL :: str_comp
!
!WRITE(output_io,1000)
!
IF(ianz==1 .or. str_comp (cpara(2), 'all', 3, lpara(2), 3) ) THEN
   DO i=1, pop_dimx
      CALL diffev_show_para_x(i)
   ENDDO
ELSE
   params:DO j=2,ianz
      DO i=1,pop_dimx
         IF(cpara(j) == pop_name(i)) THEN
            CALL diffev_show_para_x(i)
            CYCLE params
         ENDIF
      ENDDO
!     string was not found in parameter name list try to calculate
      ccpara(1) = cpara(j)
      llpara(1) = lpara(j)
      iianz     = 1
      CALL ber_params (iianz, ccpara, llpara, wwerte, MAXWW)
      IF(ier_num==0) THEN
         i = NINT(wwerte(1))
         IF(0<i .AND. i<=pop_dimx) THEN
            CALL diffev_show_para_x(i)
            CYCLE params
         ELSE
            ier_num = -14
            ier_typ = ER_APPL
            RETURN
         ENDIF
      ELSE
         RETURN
      ENDIF
   ENDDO params
ENDIF
!
!1000 FORMAT(' Parameters')
!
END SUBROUTINE diffev_show_param
SUBROUTINE diffev_show_para_x(i)
!
USE population
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=7), DIMENSION(1:2) :: ctype
CHARACTER(LEN=7), DIMENSION(1:2) :: cref
CHARACTER(LEN=7)                 :: string
INTEGER, INTENT(in) :: i
!
DATA ctype /'Reel   ','Integer' /
DATA cref  /'Refined','Fixed  ' /
!
IF(pop_refine(i)) THEN
   string = cref(1)
ELSE
   string = cref(2)
ENDIF
WRITE(output_io,2000) i, pop_name(i)(1:LEN_TRIM(pop_name(i)))
WRITE(output_io,2100) pop_xmin(i), pop_xmax(i), pop_smin(i), pop_smax(i), &
                      MINVAL(pop_t(i,:)), MAXVAL(pop_t(i,:)),             &
                      ctype(pop_type(i)), string
2000 FORMAT(/'Parameter nr., Name', I4, 1x, A)
2100 FORMAT( 'Limits hard, start ', 3(G18.10E3,2x),G18.10E3 / &
             'Current min/max    ', 2(G18.10E3,2x),A7, 2x, A7)
END SUBROUTINE diffev_show_para_x
!
!*****7*****************************************************************3*******
!
SUBROUTINE diffev_show_population(ianz, cpara, lpara, MAXW)
!
USE population
USE diff_evol
!
USE errlist_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: ianz
INTEGER                            , INTENT(IN)    :: MAXW
CHARACTER(LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
INTEGER,          DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
!
CHARACTER(LEN=20), DIMENSION(0:1) :: cdonor
CHARACTER(LEN=35), DIMENSION(0:2) :: csel
INTEGER :: i
DATA cdonor /'Add to random member', 'Add to best member  '/
DATA csel   /'Select parent/child                ',  &
             'Select best of all parents+children',  &
             'Select best child                  '/
!
WRITE(output_io, 1000) pop_gen, pop_n, pop_c, pop_dimx
IF(pop_trial_file_wrt) THEN
   WRITE(output_io, 1100) pop_trialfile (1:LEN_TRIM(pop_trialfile))
ELSE
   WRITE(output_io, 1100) 'Silent mode, not written'
ENDIF
IF(pop_result_file_rd) THEN
   WRITE(output_io, 1200) trial_results (1:LEN_TRIM(trial_results))
ELSE
   WRITE(output_io, 1100) 'Silent mode, not read'
ENDIF
WRITE(output_io, 1300) parent_results(1:LEN_TRIM(parent_results))
WRITE(output_io, 1400) parent_summary(1:LEN_TRIM(parent_summary))
WRITE(output_io, 1500) parent_current(1:LEN_TRIM(parent_current))
DO i= 1, pop_back_number
   WRITE (output_io, 1520) pop_back_fil(i)(1:LEN_TRIM(pop_back_fil(i))), &
                     pop_back_ext(i)(1:LEN_TRIM(pop_back_ext(i))), &
                     pop_back_trg(i)(1:LEN_TRIM(pop_back_trg(i)))
ENDDO
WRITE(output_io, 2000) cdonor(diff_donor_mode), csel(diff_sel_mode)
WRITE(output_io, 2100) diff_cr, diff_f, diff_local, diff_k
!
1000 FORMAT ('generation members children parameters', 4(i8,2x))
1100 FORMAT ('trial file  : ',a)
1200 FORMAT ('result file : ',a)
1300 FORMAT ('log file    : ',a)
1400 FORMAT ('summary file: ',a)
1500 FORMAT ('last    file: ',a)
1520 FORMAT ('backup source/ext/target: ',a,1x,a,1x,a)
2000 FORMAT ('Donor, selection mode   : ',a,'; ',a)
2100 FORMAT (4x,'Cross_over',10x,'Factor',14x,'Local probability ',2x,&
             'Position'/,4(2x, E18.10))
!
END SUBROUTINE diffev_show_population
!
END MODULE diffev_show_mod
