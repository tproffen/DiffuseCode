MODULE phases_stack_mod
!-
!  Specific tasks for stacking faults 
!+
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE phases_place_stack_form(n_layers, st_nlayer, st_ncunit)
!-
!  Place the form factors for the current layer into the phases section
!
USE discus_allocate_appl_mod
USE crystal_mod
USE crystal_task_mod
USE diffuse_mod
USE phases_mod
USE powder_mod
USE powder_tables_mod
!
use precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: n_layers  ! Number of layers for current layer type
INTEGER, INTENT(IN) :: st_nlayer ! Total number of layers
INTEGER, INTENT(IN) :: st_ncunit ! Total number of layers
INTEGER :: k,iscat, jscat
INTEGER :: n_pha=0    ! Number of phases
INTEGER :: n_pts=0    ! Number of points in powder pattern
INTEGER :: n_scat=0   ! Number of atom types
INTEGER :: npkt
!
REAL(KIND=PREC_DP)            :: signum
!
npkt = NINT((pow_qmax-pow_qmin)/pow_deltaq) + 1
!
IF(cr_nscat+pha_nscat(pha_curr) > PHA_MAXSCAT .OR.  &
   npkt > PHA_MAXPTS                          .OR.  &
   pha_curr > PHA_MAXPHA                            ) THEN
   n_pha  = MAX(n_pha, pha_curr)
   n_pts  = MAX(n_pts, npkt , PHA_MAXPTS)
   n_scat = MAX(n_pha, cr_nscat + pha_nscat(pha_curr))
   CALL alloc_phases(n_pha, (/n_pts,1,1/), n_scat)
ENDIF
!
CALL crystal_calc_mass
!
pha_weight(pha_curr) = pha_weight(pha_curr) + cr_mass   ! accumulate masses of layer types
pha_calc (pha_curr)  = pow_four_type      ! Fourier type Complete / Debye
!
!write(*,*) ' cfact_pure ', allocated(cfact_pure) 
!write(*,*) ' pha_form   ', allocated(pha_form) 
!write(*,*) ' powder_istl', allocated(powder_istl) 
!write(*,*) ' cfact_pure ', allocated(cfact_pure), lbound(cfact_pure), ubound(cfact_pure)
!write(*,*) ' pha_form   ', allocated(pha_form)  , lbound(pha_form)  , ubound(pha_form) 
!write(*,*) ' powder_istl', allocated(powder_istl)  , lbound(powder_istl)  , ubound(powder_istl) 
!write(*,*) ' npkt   ', npkt
DO iscat = 1, cr_nscat
   jscat = pha_nscat(pha_curr) + iscat
   DO k=0, npkt
      signum = 1.0D0
      IF(REAL(cfact_pure(1, iscat))< 0.0D0) signum = -1.0D0
      pha_form(k, jscat, pha_curr) = SQRT( DBLE (       cfact_pure (powder_istl (k), iscat)   *   &
                                                 CONJG (cfact_pure (powder_istl (k), iscat) )  )  &
                                         )                                                        &
                                     * signum
   ENDDO
   pha_adp  (jscat, pha_curr) = cr_dw(iscat)
   pha_occ  (jscat, pha_curr) = cr_occ(iscat)
   pha_niscat(jscat, pha_curr) = cr_niscat(iscat)*n_layers  ! cr_niscat is the number of atoms of type iscat
ENDDO
!
pha_nscat(pha_curr)  = pha_nscat(pha_curr) + cr_nscat           ! Number of atom types for this phase
pha_ncreal(pha_curr) = pha_ncreal(pha_curr) + (cr_ncreal)* &
                       REAL(n_layers, kind=PREC_DP)/REAL(st_nlayer, kind=PREC_DP)*st_ncunit ! Add relative amount of 
pha_nreal(pha_curr)  = pha_nreal(pha_curr) + REAL(cr_nreal, kind=PREC_DP)*REAL(n_layers, kind=PREC_DP)  ! Add absolute amount of 
!write(*,*) 'PHASES_STACK ', pha_n, pha_curr, pha_nscat(pha_curr)
!write(*,*) 'form         ', pha_form(1, 1:pha_nscat(pha_curr), pha_curr)
!write(*,*) 'adp          ', pha_adp (   1:pha_nscat(pha_curr), pha_curr)
!write(*,*) 'NISCAT       ', pha_niscat( 1:pha_nscat(pha_curr), pha_curr)
!write(*,*) 'NCATOMS      ', cr_ncatoms, cr_ncreal, pha_ncreal(pha_curr)
!write(*,*) 'NREAL        ', cr_nreal, pha_nreal(pha_curr), pha_ncreal(pha_curr)
!write(*,*) 'CALC NCREAL  ', cr_ncreal, n_layers, st_nlayer, st_ncunit
!                                                           ! atoms pere unit cell
!
END SUBROUTINE phases_place_stack_form
!
!*******************************************************************************
!
SUBROUTINE phases_place_stack_rese
!-
! Reset entries in case of single powder pattern
!
use phases_mod
!
implicit none
!
pha_curr = 1
pha_weight = 0.0D0
pha_niscat = 0
pha_nscat  = 0
pha_ncreal = 0.0
pha_nreal  = 0.0
pha_adp    = 0.0
pha_occ    = 0.0
!
END SUBROUTINE phases_place_stack_rese
!
!*******************************************************************************
!
END MODULE phases_stack_mod
