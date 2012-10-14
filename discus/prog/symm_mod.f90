MODULE symm_mod
!+
!
!     variables needed for the generalized symmetry operations
!-
!
!
SAVE
!
INTEGER, PARAMETER       :: SYM_RUN_MOLECULE = 0
INTEGER, PARAMETER       :: SYM_RUN_DOMAIN   = 1
!
INTEGER                                      ::  SYM_MAXSCAT = 1
!
LOGICAL,          DIMENSION(:), ALLOCATABLE  ::  sym_latom   ! (0:SYM_MAXSCAT)
!
CHARACTER (LEN = 4)      :: sym_incl       = 'list'
INTEGER                  :: sym_sel_mode   = 0
INTEGER, DIMENSION(0:1)  :: sym_sel_prop   = (/0,0/)
INTEGER                  :: sym_start      = 1
INTEGER                  :: sym_end        = 1
INTEGER                  :: sym_power      = 1
LOGICAL                  :: sym_mode       = .true.
LOGICAL                  :: sym_new        = .false.
LOGICAL                  :: sym_orig_mol   = .true.
LOGICAL                  :: sym_power_mult = .true.
LOGICAL                  :: sym_type       = .true.
LOGICAL                  :: sym_sel_atom   = .true.
LOGICAL                  :: sym_dom_mode_shape = .false.
LOGICAL                  :: sym_dom_mode_atom  = .false.
REAL                     :: sym_angle     = 0.0
REAL, DIMENSION(3)       :: sym_hkl       = (/0.0, 0.0, 1.0/)
REAL, DIMENSION(3)       :: sym_orig      = (/0.0, 0.0, 0.0/)
REAL, DIMENSION(3)       :: sym_or_tr     = (/0.0, 0.0, 0.0/)
REAL, DIMENSION(3)       :: sym_trans     = (/0.0, 0.0, 0.0/)
REAL, DIMENSION(3)       :: sym_uvw       = (/0.0, 0.0, 1.0/)
REAL, DIMENSION(4,4)     :: sym_mat       = 0.0
REAL, DIMENSION(4,4)     :: sym_rmat      = 0.0
!
INTEGER                  :: sym_size_of  ! Bytes allocated for symmetry
!
!     COMMON /symm_par/ sym_sel_prop,                                   &
!    &                  sym_start,sym_end,                              &
!    &                  sym_power,sym_sel_mode,                         &
!    &                  sym_latom,sym_new,sym_orig_mol,                 &
!    &                  sym_mode,sym_power_mult,sym_type,sym_sel_atom,  &
!    &                  sym_angle,sym_hkl,sym_orig,sym_or_tr,           &
!    &                  sym_trans,sym_uvw,                              &
!    &                  sym_mat,sym_rmat,sym_incl,                      &
!    &                  sym_dom_mode_atom,sym_dom_mode_shape
END MODULE symm_mod
