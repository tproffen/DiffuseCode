MODULE diffev_setup_cost_mod
!
CONTAINS
!
!
!
SUBROUTINE diffev_set_sub_cost
!
! Sets the specific DIFFEV interfaces for the cost calculation function
!
USE set_sub_generic_mod
!
INTERFACE
   SUBROUTINE diffev_execute_cost( repeat,    &
                            prog_len,         &
                            prog  ,  prog_l , &
                            mac_len,          &
                            mac   ,  mac_l  , &
                            direc_len,        &
                            direc ,  direc_l, &
                            kid   ,  indiv  , &
                            n_rvalue_i, n_rvalue_o, NRVAL,  &
                            rvalue, l_rvalue, &
                            output_len,       &
                            output, output_l, &
                            generation, member, &
                            children, parameters, &
                                    nindiv  , &
                            trial_n,         &
                            trial_v, NTRIAL, &
                            l_get_random_state,     &
                            rd_nseeds,rd_seeds,     &
                            l_first_job,            &
                            ierr, ierr_typ,         &
                            ierr_msg_l, ierr_msg_n, ierr_msg )
!
   USE precision_mod
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
   CHARACTER(LEN=prog_len  ), INTENT(IN) :: prog
   CHARACTER(LEN=mac_len   ), INTENT(IN) :: mac
   CHARACTER(LEN=direc_len ), INTENT(IN) :: direc
   INTEGER                , INTENT(IN) :: kid
   INTEGER                , INTENT(IN) :: indiv
   CHARACTER(LEN=output_len), INTENT(IN) :: output
   INTEGER                , INTENT(IN ):: n_rvalue_i
   INTEGER                , INTENT(OUT):: n_rvalue_o
   INTEGER                , INTENT(IN ):: NRVAL
   REAL(kind=PREC_DP), DIMENSION(0:NRVAL)    , INTENT(OUT):: rvalue
   LOGICAL                , INTENT(OUT):: l_rvalue
   INTEGER                , INTENT(IN) :: generation
   INTEGER                , INTENT(IN) :: member
   INTEGER                , INTENT(IN) :: children
   INTEGER                , INTENT(IN) :: parameters
   INTEGER                , INTENT(IN) :: nindiv
   INTEGER                , INTENT(IN) :: NTRIAL
   CHARACTER(LEN=16),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_n
   REAL(kind=PREC_DP),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
   LOGICAL                , INTENT(IN)  :: l_get_random_state
   INTEGER                , INTENT(OUT) :: rd_nseeds
   INTEGER, DIMENSION(64) , INTENT(OUT) :: rd_seeds
   LOGICAL                , INTENT(IN ) :: l_first_job
   INTEGER                , INTENT(OUT):: ierr 
   INTEGER                , INTENT(OUT):: ierr_typ
   integer                , INTENT(in)  :: ierr_msg_l
   integer                , INTENT(in)  :: ierr_msg_n
   character(len=ierr_msg_l),dimension(ierr_msg_n) :: ierr_msg
!
   END SUBROUTINE diffev_execute_cost
END INTERFACE
!
p_execute_cost      => diffev_execute_cost
!
END SUBROUTINE diffev_set_sub_cost
!
end MODULE diffev_setup_cost_mod
