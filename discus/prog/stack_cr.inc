!-
!     Definition of all variables for the microdomain - crystal
!+
      CHARACTER*80 st_name
      CHARACTER*16 st_spcgr
      CHARACTER*4  st_at_lis(0:ST_MAX_SCAT)
      CHARACTER*4  sa_at_lis(0:ST_MAX_SCAT)
!
      REAL         st_scat(9,0:ST_MAX_SCAT)
      REAL         st_dw(0:ST_MAX_SCAT)
      REAL         sa_dw(0:ST_MAX_SCAT)
      REAL         ST_delfr(0:ST_MAX_SCAT)
      REAL         st_delfi(0:ST_MAX_SCAT)
      REAL         st_pos(3,ST_MMAX)
      REAL         sa_pos(3,ST_MAX_SCAT)
      REAL         st_a0(3)
      REAL         st_win(3)
      REAL         st_v
      REAL         st_gten(3,3)
      REAL         st_eps(3,3,3)
      REAL         st_ar(3)
      REAL         st_wrez(3)
      REAL         st_vr
      REAL         st_rten(3,3)
      REAL         st_reps(3,3,3)
      REAL         st_dim(3,2)
!
      INTEGER*4    sa_natoms
      INTEGER*4    st_nscat
      INTEGER*4    sa_nscat
      INTEGER*4    st_iscat(ST_MMAX)
      INTEGER*4    st_prop (ST_MMAX)
      INTEGER*4    sa_iscat(ST_MAX_SCAT)
      INTEGER*4    sa_prop (ST_MAX_SCAT)
      INTEGER      st_ndel
      INTEGER      st_natoms
      INTEGER      st_spcgr_ianz
      INTEGER      st_spcgr_para
!
!       COMMON /stack_cr/ st_name,st_spcgr,st_a0,st_win,st_v,
!    &                st_gten,st_eps,
!    &       st_ar,st_wrez,st_vr,st_rten,st_reps,st_natoms,
!    &       st_nscat,st_dw,st_scat,st_delfr,st_delfi,
!    &       st_at_lis,st_atom,st_pos,st_iscat,st_ndel,st_dim
!    &                st_spcgr_ianz,st_spcgr_para
!
