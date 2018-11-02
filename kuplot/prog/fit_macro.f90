SUBROUTINE theory_macro(xx, f, df, i)
!
! Theory function that executes a macro
!
USE kuplot_config
USE kuplot_mod
USE fit_params_mod
USE fit_set_sub_mod
!
USE errlist_mod
USE doact_mod
USE do_if_mod
USE do_set_mod
USE param_mod
USE prompt_mod
USE sup_mod
USE set_sub_generic_mod
USE variable_mod
!
IMPLICIT NONE
!
REAL                    , INTENT(IN)  :: xx
REAL                    , INTENT(OUT) :: f
REAL, DIMENSION(MAXPARA), INTENT(OUT) :: df
INTEGER                 , INTENT(IN)  :: i
!
!NTEGER, PARAMETER :: MAXP = 2
!HARACTER(LEN=1024), DIMENSION(MAXP) :: cpara
!NTEGER            , DIMENSION(MAXP) :: lpara
!EAL               , DIMENSION(MAXP) :: werte
CHARACTER(LEN=1024) :: macro_name
CHARACTER(LEN=1024) :: line
CHARACTER(LEN=1024) :: zeile
CHARACTER(LEN=1024) :: string
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
