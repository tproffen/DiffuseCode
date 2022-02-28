MODULE kuplot_theory_macro_mod
!!
!!*******************************************************************************
!!
CONTAINS
!!
!!***********************************************************************
!!
!OLD SUBROUTINE theory_macro_n(MAXP, ix, iy, xx, yy, NPARA, params, par_names,          &
!OLD                           prange, l_do_deriv, data_dim, &
!OLD                           data_data, data_sigma, data_x, data_y, &
!OLD                           data_calc, kupl_last,      &
!OLD                           ymod, dyda, LDERIV)
!OLD !
!OLD use precision_mod
!OLD USE kuplot_config
!OLD !
!OLD ! Calculate a theory function using a user supplied macro
!OLD !
!OLD IMPLICIT NONE
!OLD !
!OLD INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
!OLD INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
!OLD INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
!OLD REAL(kind=PREC_DP)                                   , INTENT(IN)  :: xx      ! Point value  along x
!OLD REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
!OLD INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
!OLD REAL(kind=PREC_DP),  DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
!OLD CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
!OLD REAL(kind=PREC_DP),  DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
!OLD LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
!OLD INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
!OLD REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
!OLD REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_sigma  ! Data sigmas
!OLD REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
!OLD REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
!OLD REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
!OLD INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
!OLD REAL(kind=PREC_DP)                                   , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
!OLD REAL(kind=PREC_DP),  DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
!OLD LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!OLD !
!OLD !REAL                     :: xx
!OLD REAL(kind=PREC_DP)                     :: f
!OLD REAL(kind=PREC_DP), DIMENSION(MAXPARA) :: df
!OLD !INTEGER                  :: ix
!OLD !
!OLD CALL theory_macro(xx, f, df, ix)
!OLD data_calc(ix, iy) = f
!OLD !
!OLD END SUBROUTINE theory_macro_n
!
!*******************************************************************************
!
SUBROUTINE theory_macro(xx, f, df, i)
!
! Theory function that executes a macro
!
USE kuplot_config
USE kuplot_mod
use kuplot_gsas_mod
USE fit_params_mod
USE fit_set_sub_mod
!
USE errlist_mod
USE doact_mod
USE do_if_mod
USE do_set_mod
USE lib_macro_func
USE param_mod
USE precision_mod
USE prompt_mod
USE sup_mod
USE set_sub_generic_mod
USE variable_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP)      , INTENT(IN)  :: xx
REAL(kind=PREC_DP)      , INTENT(OUT) :: f
REAL(kind=PREC_DP), DIMENSION(MAXPARA), INTENT(OUT) :: df
INTEGER                 , INTENT(IN)  :: i
!
CHARACTER(LEN=PREC_STRING) :: macro_name
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: zeile
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=4   ) :: befehl
INTEGER             :: length
INTEGER             :: lp
INTEGER             :: lbef
!INTEGER             :: ianz
INTEGER             :: lstring
INTEGER             :: j, ix, iy
LOGICAL             :: lend
LOGICAL             :: l_prompt_restore
!
!
INTERFACE
   SUBROUTINE kuplot_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE kuplot_mache_kdo
END INTERFACE
!
iwert = i
!
l_prompt_restore = .FALSE.
IF(output_status == OUTPUT_SCREEN .AND.lturn_off) THEN
   l_prompt_restore = .TRUE.
   string = 'prompt, off, off, save'
   lstring = 22
   CALL do_set(string, lstring)
ENDIF
df (:) = 0.0
!                                                                       
!--Determine x and y for 2D and 3D data
!                                                                       
IF (lni (ikfit) ) THEN
   ix = (ABS (i) - 1) / ny (ikfit) + 1
   iy = ABS (i) - (ix - 1) * ny (ikfit)
   rpara (0) = x (offxy (ikfit - 1) + ix)
   rpara (1) = y (offxy (ikfit - 1) + iy)
ELSE
   rpara (0) = xx
ENDIF
!
macro_name = fit_func
length = LEN_TRIM(macro_name)
!
CALL file_kdo (macro_name, length)
IF(ier_num/=0) RETURN
!
lmacro_close = .FALSE.        ! Do not close macros in do-loops, return instead
!
CALL fit_set_sub
!
lend = .FALSE.
main: DO
   CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt)
   IF (ier_num == 0) THEN
      IF (line == ' '.OR.line (1:1)  == '#' .OR. LINE=='!') CYCLE main
      IF (befehl (1:3)  == 'do '.OR.befehl (1:2)  == 'if') THEN
         CALL do_loop (line, lend, length) !, kuplot_mache_kdo)
         IF(ier_num/=0) EXIT main
      ELSE
         CALL fit_mache_kdo(line, lend, length)
      ENDIF
      IF(lend .OR. ier_num/=0) EXIT main
   ENDIF
ENDDO main
!
CALL macro_terminate
!
IF(l_prompt_restore) THEN
   string = 'prompt, on,on'
   lstring = 13
   CALL do_set(string, lstring)
ENDIF
!    Transfer function and derivatives
f     = var_val(VAR_FIT_VALUE)
DO j=1,npara
   IF (pinc (j) /=  0) THEN
      df(j) = kupl_deriv(j)
   ELSE
      df(j) = 0
   ENDIF
ENDDO
!
p_mache_kdo   => kuplot_mache_kdo
lmacro_close = .TRUE.        ! Do close macros in do-loops
!
END SUBROUTINE theory_macro
!
!*******************************************************************************
!
END MODULE kuplot_theory_macro_mod
