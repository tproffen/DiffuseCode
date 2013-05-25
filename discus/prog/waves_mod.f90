MODULE waves_mod
!+
!     This file contains COMMON block for wave input
!-
!
SAVE
!
INTEGER, PARAMETER  ::  WV_RAND      =  1
INTEGER, PARAMETER  ::  WV_FIX       =  2
!
INTEGER, PARAMETER  ::  WV_LONG      =  1
INTEGER, PARAMETER  ::  WV_TRANS     =  2
INTEGER, PARAMETER  ::  WV_DENS      =  3
INTEGER, PARAMETER  ::  WV_ROT       =  4
!
INTEGER, PARAMETER  ::  WV_SINUS     =  1
INTEGER, PARAMETER  ::  WV_BOX       =  2
INTEGER, PARAMETER  ::  WV_TRIANGLE  =  3
!
CHARACTER(LEN=4)    ::  wv_func      = 'sinu'
!
INTEGER             ::  WV_MAXSCAT = 1
!
INTEGER,          DIMENSION(:), ALLOCATABLE  ::  wv_repl      ! (0:WV_MAXSCAT)
LOGICAL,          DIMENSION(:), ALLOCATABLE  ::  wv_latom     ! (0:WV_MAXSCAT)
LOGICAL,          DIMENSION(:), ALLOCATABLE  ::  wv_latom_rot ! (0:WV_MAXSCAT)
!
INTEGER             ::  wv_iwave         = WV_LONG
INTEGER             ::  wv_ifunc         = WV_SINUS
INTEGER             ::  wv_sel_prop(0:1) = 0
INTEGER             ::  wv_phase_typ     = WV_FIX
REAL                ::  wv_wave(3)       = (/1.,0.,0./)
REAL                ::  wv_swing(3)      = (/0.,1.,0./)
REAL                ::  wv_rot_uvw(3)    = (/0.,0.,1./)
REAL                ::  wv_rot_orig(3)   = 0.0
REAL                ::  wv_amp           = 0.5
REAL                ::  wv_rlam          =50.0
REAL                ::  wv_phase         = 0.0
REAL                ::  wv_amp0          = 0.0
REAL                ::  wv_plow          = 0.0
REAL                ::  wv_phigh         = 0.0
REAL                ::  wv_asym          = 0.5
LOGICAL             ::  wv_sel_atom      = .true.
LOGICAL             ::  wv_lacoust       = .true.
LOGICAL             ::  wv_viceversa     = .false.
!
INTEGER             ::  wv_size_of  = 0 ! Bytes allocated for waves      
!
!     COMMON /waves_cmm/ wv_wave,wv_swing,wv_rlam,wv_amp,wv_amp0,       &
!    &                   wv_iwave,wv_lacoust,wv_phase,wv_plow,          &
!    &                   wv_phigh,wv_phase_typ,wv_latom,wv_repl,        &
!    &                   wv_sel_prop,                                   &
!    &                   wv_md_s,wv_md_nr,wv_func,wv_sel_atom,          &
!    &                   wv_ifunc,wv_rot_uvw,wv_rot_orig,wv_asym
!
END MODULE waves_mod
