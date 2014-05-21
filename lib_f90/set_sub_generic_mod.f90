MODULE set_sub_generic_mod
!
!  sets generic interfaces for the functions/subroutines
!  in DISCUS/DIFFEV/KUPLOT/MIXSCAT that are referenced 
!  within lib_f90 under identical name 
!  A procedure pointer is defined for eah to allow flexible 
!  switching between the routines
!
!  SUBROUTINE mache_kdo
!  SUBROUTINE ersetz_para
!  SUBROUTINE upd_para 
!  SUBROUTINE calc_intr_spec
!  SUBROUTINE validate_var_spec
!
!
INTERFACE
   SUBROUTINE mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE errlist_appl
   END SUBROUTINE errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL, DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE upd_para (ctype, ww, maxw, wert, ianz)
!
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL              , INTENT(IN   )    :: wert
!
   END SUBROUTINE upd_para 
END INTERFACE
!
INTERFACE
   SUBROUTINE calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL               , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE calc_intr_spec
END INTERFACE
!
INTERFACE
   SUBROUTINE validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE validate_var_spec
END INTERFACE
!
PROCEDURE(mache_kdo     )   , POINTER :: p_mache_kdo      => NULL()
PROCEDURE(errlist_appl  )   , POINTER :: p_errlist_appl   => NULL()
PROCEDURE(ersetz_para   )   , POINTER :: p_ersetz_para    => NULL()
PROCEDURE(upd_para      )   , POINTER :: p_upd_para       => NULL()
PROCEDURE(calc_intr_spec)   , POINTER :: p_calc_intr_spec => NULL()
PROCEDURE(validate_var_spec), POINTER :: p_validate_var_spec => NULL()
!
END MODULE set_sub_generic_mod
