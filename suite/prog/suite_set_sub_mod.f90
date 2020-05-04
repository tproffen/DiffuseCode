MODULE suite_set_sub_mod
!
CONTAINS
!
SUBROUTINE suite_set_sub
!
! Sets the specific SUITE  interfaces for routines that are refecenced in
! LIB_F90 by their generic names
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE suite_mache_kdo (line, lend, length)
!                                                                       
   CHARACTER (LEN= *  ), INTENT(INOUT) :: line
   LOGICAL             , INTENT(  OUT) :: lend
   INTEGER             , INTENT(INOUT) :: length
!
   END SUBROUTINE suite_mache_kdo
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_errlist_appl
   END SUBROUTINE suite_errlist_appl
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!
   USE precision_mod
   CHARACTER (LEN= * )  , INTENT(INOUT) :: string
   INTEGER              , INTENT(IN   ) :: ikl
   INTEGER              , INTENT(IN   ) :: iklz
   INTEGER              , INTENT(INOUT) :: ll
   INTEGER              , INTENT(IN   ) :: maxw
   INTEGER              , INTENT(IN   ) :: ianz
   REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN   ) :: ww
!
   END SUBROUTINE suite_ersetz_para
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr)
!
   USE precision_mod
   CHARACTER (LEN=* ), INTENT(IN   )    :: ctype
   INTEGER           , INTENT(IN   )    :: maxw
   INTEGER           , INTENT(IN   )    :: ianz
   INTEGER           , INTENT(IN   )    :: ww (maxw)
   REAL(KIND=PREC_DP), INTENT(IN   )    :: wert
   CHARACTER (LEN=* ), INTENT(IN   )    :: cstring
   INTEGER, DIMENSION(2), INTENT(IN)    :: substr ! Indices of substring
!
   END SUBROUTINE suite_upd_para
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_get_var_type(line,length, var_is_type)
!
   CHARACTER(LEN=*)     , INTENT(IN)  :: line
   INTEGER              , INTENT(IN)  :: length
   INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
   END SUBROUTINE suite_get_var_type
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_calc_intr_spec (string, line, ikl, iklz, ww, laenge, lp)
!
   USE precision_mod
   CHARACTER (LEN= * ), INTENT(INOUT) :: string
   CHARACTER (LEN= * ), INTENT(INOUT) :: line
   INTEGER            , INTENT(IN   ) :: ikl
   INTEGER            , INTENT(IN   ) :: iklz
   REAL(KIND=PREC_DP) , INTENT(INOUT) :: ww
   INTEGER            , INTENT(INOUT) :: laenge
   INTEGER            , INTENT(INOUT) :: lp
!
   END SUBROUTINE suite_calc_intr_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_calc_intr_log_spec(string, length)
!
   IMPLICIT NONE
   CHARACTER(LEN=*) , INTENT(INOUT) :: string
   INTEGER          , INTENT(INOUT) :: length
!
   END SUBROUTINE suite_calc_intr_log_spec
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_validate_var_spec (string, lp)
!
   CHARACTER (LEN= * ), INTENT(IN   ) :: string
   INTEGER            , INTENT(IN   ) :: lp
!
   END SUBROUTINE suite_validate_var_spec 
END INTERFACE
!
INTERFACE
   SUBROUTINE suite_top(zeile)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
!
   END SUBROUTINE suite_top
END INTERFACE
!
INTERFACE
   SUBROUTINE diffev_loop_mpi(prog_n, prog_l, mac_n, mac_l, out_n, out_l, repeat, nindiv)
!
CHARACTER (LEN=*), INTENT(IN) :: prog_n
CHARACTER (LEN=*), INTENT(IN) :: mac_n
CHARACTER (LEN=*), INTENT(IN) :: out_n
INTEGER          , INTENT(IN) :: prog_l
INTEGER          , INTENT(IN) :: mac_l
INTEGER          , INTENT(IN) :: out_l
LOGICAL          , INTENT(IN) :: repeat
INTEGER          , INTENT(IN) :: nindiv
!
   END SUBROUTINE diffev_loop_mpi
END INTERFACE
!
!
p_mache_kdo         => suite_mache_kdo
p_errlist_appl      => suite_errlist_appl
p_ersetz_para       => suite_ersetz_para
p_upd_para          => suite_upd_para
p_calc_intr_spec    => suite_calc_intr_spec
p_calc_intr_log_spec=> suite_calc_intr_log_spec
p_validate_var_spec => suite_validate_var_spec
p_top               => suite_top
p_get_var_type      => suite_get_var_type
p_loop_mpi          => diffev_loop_mpi
!
END SUBROUTINE suite_set_sub
!
SUBROUTINE  suite_set_sub_cost
!
! Sets the specific SUITE  interfaces for the cost calculation function
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE suite_execute_cost( repeat,    &
                            prog_len,         &
                            prog  ,  prog_l , &
                            mac_len,          &
                            mac   ,  mac_l  , &
                            direc_len,        &
                            direc ,  direc_l, &
                            kid   ,  indiv  , &
                            n_rvalue_i, n_rvalue_o, NRVAL, &
                            rvalue, l_rvalue, &
                            output_len,       &
                            output, output_l, &
                            generation, member, &
                            children, parameters, &
                                    nindiv  , &
                            trial_n,              &
                            trial_v,  NTRIAL,     &
                            l_get_random_state,     &
                            rd_nseeds,rd_seeds,     &
                            l_first_job,            &
                            ierr )
!
   IMPLICIT NONE
   LOGICAL                , INTENT(IN) :: repeat
   INTEGER                , INTENT(IN) :: prog_len
   INTEGER                , INTENT(IN) :: prog_l
   INTEGER                , INTENT(IN) :: mac_len
   INTEGER                , INTENT(IN) :: mac_l
   INTEGER                , INTENT(IN) :: direc_len
   INTEGER                , INTENT(IN) :: direc_l
   INTEGER                , INTENT(IN) :: output_len
   INTEGER                , INTENT(IN) :: output_l
   CHARACTER(LEN=prog_len), INTENT(IN) :: prog
   CHARACTER(LEN=mac_len ), INTENT(IN) :: mac
   CHARACTER(LEN=direc_len ), INTENT(IN) :: direc
   INTEGER                , INTENT(IN) :: kid
   INTEGER                , INTENT(IN) :: indiv
   CHARACTER(LEN=output_len), INTENT(IN) :: output
   INTEGER                , INTENT(IN ):: n_rvalue_i
   INTEGER                , INTENT(OUT):: n_rvalue_o
   INTEGER                , INTENT(IN ):: NRVAL
   REAL, DIMENSION(0:NRVAL), INTENT(OUT):: rvalue
   LOGICAL                , INTENT(OUT):: l_rvalue
   INTEGER                , INTENT(IN) :: generation
   INTEGER                , INTENT(IN) :: member
   INTEGER                , INTENT(IN) :: children
   INTEGER                , INTENT(IN) :: parameters
   INTEGER                , INTENT(IN) :: nindiv
   INTEGER                , INTENT(IN) :: NTRIAL
   CHARACTER(LEN=16),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_n
   REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
   LOGICAL                , INTENT(IN)  :: l_get_random_state
   INTEGER                , INTENT(OUT) :: rd_nseeds
   INTEGER, DIMENSION(64) , INTENT(OUT) :: rd_seeds
   LOGICAL                , INTENT(IN)  :: l_first_job
   INTEGER                , INTENT(OUT):: ierr 
!
   END SUBROUTINE suite_execute_cost
END INTERFACE
!
p_execute_cost      =>  suite_execute_cost
!
END SUBROUTINE suite_set_sub_cost
!
SUBROUTINE suite_set_sub_branch
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE suite_branch(zeile, length, lreset)
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
   END SUBROUTINE suite_branch
END INTERFACE
!
p_branch            =>  suite_branch 
!
END SUBROUTINE suite_set_sub_branch
!
!
END MODULE suite_set_sub_mod
