SUBROUTINE refine_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
!
USE refine_params_mod
!
USE blanks_mod
USE errlist_mod 
USE param_mod 
USE lib_errlist_func
USE lib_upd_mod
USE lib_length
USE precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN= * )  , INTENT(INOUT) :: string 
INTEGER              , INTENT(IN   ) :: ikl
INTEGER              , INTENT(IN   ) :: iklz
INTEGER              , INTENT(INOUT) :: ll
INTEGER              , INTENT(IN   ) :: maxw
INTEGER              , INTENT(IN   ) :: ianz
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN   ) :: ww
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) :: zeile
INTEGER :: laenge
INTEGER :: lcomm
INTEGER :: ltyp
INTEGER :: kpara
INTEGER :: kpara2
!
!
CALL lib_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
IF(ier_num == 0) RETURN
CALL no_error
!
laenge = ll
ltyp = 1
zeile = ' '
kpara = NINT(ww(1))
kpara2 = 0
IF (maxw.ge.2) THEN
   kpara2 = NINT(ww(2))
ENDIF
!                                                                       
lcomm = length_com (string, ikl)
!                                                                       
IF (lcomm.eq.3) THEN
!
   IF(ikl > lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
!
   IF(string(ikl-3:ikl - 1) .eq.'par') THEN
      IF(ianz==1) THEN
         IF(1<=kpara .AND.kpara<=refine_par_n ) THEN
            WRITE(zeile(ikl-3:ikl+PREC_WIDTH-2), PREC_F_REAL) refine_p(kpara)
            zeile(ikl+PREC_MANTIS-lcomm:ikl+PREC_MANTIS-lcomm) = 'e'
         ELSEIF(1<=ABS(kpara) .AND.ABS(kpara)<=refine_par_n ) THEN
            WRITE(zeile(ikl-3:ikl+PREC_WIDTH-2), PREC_F_REAL) refine_f(ABS(kpara))
            zeile(ikl+PREC_MANTIS-lcomm:ikl+PREC_MANTIS-lcomm) = 'e'
         ELSE
            ier_num = -8
            ier_typ = ER_FORT
         ENDIF
      ELSE
         ier_num = -13
         ier_typ = ER_FORT
      ENDIF
   ELSEIF(string(ikl-3:ikl - 1) .eq.'sig') THEN
      IF(ianz==1) THEN
         IF(1<=kpara .AND.kpara<=refine_par_n ) THEN
            WRITE(zeile(ikl-3:ikl+PREC_WIDTH-2), PREC_F_REAL) refine_dp(kpara)
            zeile(ikl+PREC_MANTIS-lcomm:ikl+PREC_MANTIS-lcomm) = 'e'
         ELSE
            ier_num = -8
            ier_typ = ER_FORT
         ENDIF
      ELSE
         ier_num = -13
         ier_typ = ER_FORT
      ENDIF
   ELSE
      ier_num = -2
      ier_typ = ER_FORT
      RETURN
   ENDIF
!
ELSEIF (lcomm.eq.5) THEN
!
   IF(ikl > lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
!
   IF(string(ikl-5:ikl - 1) .eq.'covar') THEN
      IF(ALLOCATED(refine_cl)) THEN
         IF(ianz==2) THEN
            IF(1<=kpara .AND.kpara<=refine_par_n.AND. &
               1<=kpara2.AND.kpara2<=refine_par_n  ) THEN
               WRITE(zeile(ikl-5:ikl+PREC_WIDTH-2), PREC_F_REAL) refine_cl(kpara,kpara2)
               zeile(ikl+PREC_MANTIS-lcomm:ikl+PREC_MANTIS-lcomm) = 'e'
            ELSE
               ier_num = -8
               ier_typ = ER_FORT
            ENDIF
         ELSE
            ier_num = -13
            ier_typ = ER_FORT
         ENDIF
      ELSE
         ier_num = -12
         ier_typ = ER_APPL
      ENDIF
   ELSE
      ier_num = -2
      ier_typ = ER_FORT
      RETURN
   ENDIF
ELSE
   ier_num = - 2
   ier_typ = ER_FORT
ENDIF
!
IF (ier_num.eq.0) THEN
   ll = laenge+PREC_WIDTH - ltyp - (iklz - ikl + 1)
   IF(iklz + 1.le.laenge) zeile(ikl + PREC_WIDTH-1:ll) = string(iklz + 1:laenge)
   string = zeile
ELSE
   ll = MIN(40, laenge)
   WRITE (ier_msg (1), '(a)') string (1:ll)
ENDIF
ll  = LEN_TRIM(string)
!
END SUBROUTINE refine_ersetz_para                    
!*****7*****************************************************************
SUBROUTINE refine_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
USE errlist_mod 
USE param_mod 
USE lib_upd_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=* ), INTENT(IN   )    :: ctype 
INTEGER           , INTENT(IN   )    :: maxw
INTEGER           , INTENT(IN   )    :: ianz 
INTEGER           , INTENT(IN   )    :: ww (maxw)
REAL(KIND=PREC_DP), INTENT(IN   )    :: wert 
CHARACTER (LEN=* ), INTENT(IN   )    :: cstring
INTEGER, DIMENSION(2), INTENT(IN)    :: substr ! Indices of substring
!
CALL lib_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr)
!
END SUBROUTINE refine_upd_para                       
!*****7***************************************************************  
SUBROUTINE refine_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!-                                                                      
!     These are special intrinsic function for the DISCSU_SUITE. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!     Currently empty, needed for formal reasons.
!+                                                                      
!
USE refine_params_mod
!
USE ber_params_mod
USE build_name_mod
USE calc_expr_mod
USE errlist_mod 
USE ersetz_mod
USE get_params_mod
USE lib_length
USE param_mod 
USE precision_mod
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER   :: MAXW = 9
!                                                                       
CHARACTER (LEN= * ), INTENT(INOUT) :: string
CHARACTER (LEN= * ), INTENT(INOUT) :: line 
INTEGER            , INTENT(IN   ) :: ikl
INTEGER            , INTENT(IN   ) :: iklz
INTEGER            , INTENT(INOUT) :: laenge
INTEGER            , INTENT(INOUT) :: lp
REAL(KIND=PREC_DP) , INTENT(INOUT) :: ww
!
INTEGER              :: i, lcomm
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
!                                                                       
!                                                                       
lcomm = length_com (string(1:lp), ikl) 
ier_num = - 1 
ier_typ = ER_FORT 
DO i = 1, maxw 
   werte (i) = 0.0 
ENDDO 
!                                                                 
IF (lcomm.eq.0) THEN 
   CALL ersetz2 (string, ikl, iklz, ww, 0, laenge) 
ELSE 
   ier_num = - 3 
   ier_typ = ER_FORT 
ENDIF 
!                                                                 
!IF (ier_num.ne.0) THEN 
!   WRITE ( *, * ) string (1:len_trim(string))
!   WRITE ( *, * ) line   (1:len_trim(line))
!ENDIF 
!                                                                       
END SUBROUTINE refine_calc_intr_spec                 
!
!*****7**************************************************************** 
!
SUBROUTINE refine_calc_intr_log_spec(string, length)
!
IMPLICIT NONE
CHARACTER(LEN=*) , INTENT(INOUT) :: string
INTEGER          , INTENT(INOUT) :: length
!
END SUBROUTINE refine_calc_intr_log_spec
!*****7**************************************************************** 
SUBROUTINE refine_validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, DIFFEV specific part 
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@krist.uni-erlangen.de)    
!+                                                                      
!
USE reserved_mod
USE errlist_mod 
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN= * ), INTENT(IN   ) :: zeile 
INTEGER            , INTENT(IN   ) :: lp 
!                                                                       
INTEGER                                  :: i , length
!                                                                       
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
!                                                                       
DO i = 1, refine_reserved_n 
!  IF (index (refine_reserved (i), zeile (1:lp) ) .ne.0) then 
   length = MAX(LEN_TRIM(refine_reserved(i)), LEN_TRIM(zeile(1:lp)))
   length = MIN(length, LEN(refine_reserved), LEN(zeile))
   IF(refine_reserved (i)(1:length)== zeile(1:length) ) THEN           
      ier_num = - 25 
      ier_typ = ER_FORT 
   ENDIF 
ENDDO 
!                                                                       
END SUBROUTINE refine_validate_var_spec              
!
!*******************************************************************************
!
SUBROUTINE refine_get_var_type(line,length, var_is_type)
!
! Returns the variable type : INTEGER, REAL, CHARACTER, and Scalar versus field
! Currently the refine does not offer local variables.
!
USE constants_mod
USE lib_get_var
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)     , INTENT(IN)  :: line
INTEGER              , INTENT(IN)  :: length
INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
!INTEGER, PARAMETER :: IS_SCAL = 0
!INTEGER, PARAMETER :: IS_VEC  = 1
!INTEGER, PARAMETER :: IS_ARR  = 2
!INTEGER, PARAMETER :: MAXPAR = 24
!CHARACTER(LEN=16), DIMENSION(MAXPAR) :: refine_names
!INTEGER          , DIMENSION(MAXPAR) :: refine_type
!INTEGER          , DIMENSION(MAXPAR) :: refine_dim
!LOGICAL          , DIMENSION(MAXPAR) :: refine_ro 
!INTEGER :: i
!
!DATA refine_names  &
!    /'pdf_scal', 'pdf_dens', 'mol_type', 'mol_dens', 'mol_cont', &
!     'mol_biso', 'mol_len ', 'in_mole ', 'at_type ', 'at_name ', &
!     'sym_n   ', 'rvol    ', 'menv    ', 'cdim    ', 'vol     ', &
!     'occ     ', 'lat     ', 'env     ', 'z       ', 'y       ', &
!     'x       ', 'n       ', 'm       ', 'b       '              &
!    /
!DATA refine_type &
!    /  IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_INTE , &
!       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_CHAR ,   IS_CHAR , &
!       IS_INTE ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL , &
!       IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL , &
!       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_REAL               &
!    /
!DATA refine_dim  &
!    /  IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  ,   IS_VEC  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC                &
!    /
!DATA refine_ro  &
!    /  .FALSE. ,   .FALSE. ,   .TRUE.  ,   .FALSE. ,   .TRUE.  , &
!       .FALSE. ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  , &
!       .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  , &
!       .FALSE. ,   .FALSE. ,   .TRUE.  ,   .FALSE. ,   .FALSE. , &
!       .FALSE. ,   .TRUE.  ,   .TRUE.  ,   .FALSE.               &
!    /
!
var_is_type(:) = IS_UNKNOWN
!
!main: DO i=1, MAXPAR
!   IF(line(1:length) == refine_names(i)(1:LEN_TRIM(refine_names(i)))) THEN
!      var_is_type(1) = refine_type(i)
!      var_is_type(2) = refine_dim (i)
!      IF(refine_ro(i)) THEN
!         var_is_type(3) = IS_READ
!      ELSE
!         var_is_type(3) = IS_WRITE
!      ENDIF
!      RETURN
!   ENDIF
!ENDDO main
!
CALL lib_get_var_type(line, length, var_is_type)
!
!
END SUBROUTINE refine_get_var_typE
