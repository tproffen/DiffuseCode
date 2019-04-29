MODULE refine_reset
!
CONTAINS
!
SUBROUTINE refine_do_reset
!-
!   Reset REFINE to system start
!+
!
USE refine_allocate_appl
USE refine_blk_appl
USE refine_control_mod
USE refine_data_mod
USE refine_params_mod
USE do_variable_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=1024)                   :: zeile
!INTEGER            , PARAMETER         :: MAXW=2
!INTEGER                                :: ianz
!CHARACTER(LEN=1024), DIMENSION(1:MAXW) :: cpara
!INTEGER            , DIMENSION(1:MAXW) :: lpara
!
!LOGICAL, PARAMETER   :: is_refine = .TRUE.
INTEGER              :: lcomm
!INTEGER              :: i
!
! Remove all parameter names from the variable entry
!
zeile ='default'
lcomm = 7
!
CALL refine_do_allocate_appl(zeile,lcomm)
CALL refine_initarrays
!
refine_par_n    = 0 ! number of parameters
refine_fix_n    = 0 
ref_dim(1)      = 1
ref_dim(2)      = 1
refine_cycles   = 1
!
END SUBROUTINE refine_do_reset
!
END MODULE refine_reset
