SUBROUTINE suite_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
USE blanks_mod
USE errlist_mod 
USE param_mod 
USE lib_upd_mod
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
REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
CALL lib_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
END SUBROUTINE suite_ersetz_para                    
!*****7*****************************************************************
SUBROUTINE suite_upd_para (ctype, ww, maxw, wert, ianz) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
USE errlist_mod 
USE param_mod 
USE lib_upd_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=* ), INTENT(IN   )    :: ctype 
INTEGER           , INTENT(IN   )    :: maxw
INTEGER           , INTENT(IN   )    :: ianz 
INTEGER           , INTENT(IN   )    :: ww (maxw)
REAL              , INTENT(IN   )    :: wert 
!
CALL lib_upd_para (ctype, ww, maxw, wert, ianz)
!
END SUBROUTINE suite_upd_para                       
!*****7***************************************************************  
SUBROUTINE suite_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!-                                                                      
!     These are special intrinsic function for the DISCSU_SUITE. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!     Currently empty, needed for formal reasons.
!+                                                                      
!
USE calc_expr_mod
USE errlist_mod 
USE ersetz_mod
USE param_mod 
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER   :: maxw = 9
!                                                                       
CHARACTER (LEN= * ), INTENT(INOUT) :: string
CHARACTER (LEN= * ), INTENT(INOUT) :: line 
INTEGER            , INTENT(IN   ) :: ikl
INTEGER            , INTENT(IN   ) :: iklz
INTEGER            , INTENT(INOUT) :: laenge
INTEGER            , INTENT(INOUT) :: lp
REAL               , INTENT(INOUT) :: ww
!
INTEGER              :: i, lcomm
REAL                 :: werte (maxw)
!                                                                       
INTEGER              :: length_com 
!                                                                       
lcomm = length_com (string(1:lp), ikl) 
ier_num = - 1 
ier_typ = ER_FORT 
DO i = 1, maxw 
   werte (i) = 0.0 
ENDDO 
!                                                                 
IF (lcomm.eq.0) then 
   CALL ersetz2 (string, ikl, iklz, ww, 0, laenge) 
ELSE 
   ier_num = - 3 
   ier_typ = ER_FORT 
ENDIF 
!                                                                 
IF (ier_num.ne.0) then 
   WRITE ( *, * ) string 
   WRITE ( *, * ) line 
ENDIF 
!                                                                       
END SUBROUTINE suite_calc_intr_spec                 
!*****7**************************************************************** 
SUBROUTINE suite_validate_var_spec (zeile, lp) 
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
DO i = 1, suite_reserved_n 
!  IF (index (suite_reserved (i), zeile (1:lp) ) .ne.0) then 
   length = MAX(LEN_TRIM(suite_reserved(i)), LEN_TRIM(zeile(1:lp)))
   IF(suite_reserved (i)(1:length)== zeile(1:length) ) THEN           
      ier_num = - 25 
      ier_typ = ER_FORT 
   ENDIF 
ENDDO 
!                                                                       
END SUBROUTINE suite_validate_var_spec              
!
!*******************************************************************************
!
SUBROUTINE suite_get_var_type(line,length, var_is_type)
!
! Returns the variable type : INTEGER, REAL, CHARACTER, and Scalar versus field
! Currently the suite does not offer local variables.
!
USE constants_mod
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
!CHARACTER(LEN=16), DIMENSION(MAXPAR) :: suite_names
!INTEGER          , DIMENSION(MAXPAR) :: suite_type
!INTEGER          , DIMENSION(MAXPAR) :: suite_dim
!LOGICAL          , DIMENSION(MAXPAR) :: suite_ro 
!INTEGER :: i
!
!DATA suite_names  &
!    /'pdf_scal', 'pdf_dens', 'mol_type', 'mol_dens', 'mol_cont', &
!     'mol_biso', 'mol_len ', 'in_mole ', 'at_type ', 'at_name ', &
!     'sym_n   ', 'rvol    ', 'menv    ', 'cdim    ', 'vol     ', &
!     'occ     ', 'lat     ', 'env     ', 'z       ', 'y       ', &
!     'x       ', 'n       ', 'm       ', 'b       '              &
!    /
!DATA suite_type &
!    /  IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_INTE , &
!       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_CHAR ,   IS_CHAR , &
!       IS_INTE ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL , &
!       IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL , &
!       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_REAL               &
!    /
!DATA suite_dim  &
!    /  IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  ,   IS_VEC  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
!       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC                &
!    /
!DATA suite_ro  &
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
!   IF(line(1:length) == suite_names(i)(1:LEN_TRIM(suite_names(i)))) THEN
!      var_is_type(1) = suite_type(i)
!      var_is_type(2) = suite_dim (i)
!      IF(suite_ro(i)) THEN
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
END SUBROUTINE suite_get_var_typE
