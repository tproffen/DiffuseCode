MODULE mixscat_mod
!
USE config_mod
!
IMPLICIT NONE
PUBLIC
SAVE
!
!*****7****************************************************************
   CHARACTER(LEN=80), DIMENSION(1:MAXDSET) :: exp_name  ! (MAXDSET)
!
   CHARACTER(LEN=4 ), DIMENSION(1:MAXELEM) :: elem_name ! (0:MAXELEM)
!
   INTEGER, DIMENSION(1:2)                 :: elem_diff ! (2)
   INTEGER                                 :: elem_n
!
   INTEGER, DIMENSION(1:MAXDSET)           :: exp_np ! (MAXDSET)
   INTEGER                                 :: exp_nd
!
   INTEGER                                 :: pdf_diff_np
!
   INTEGER                                 :: fit_np
   INTEGER                                 :: npara
   INTEGER                                 :: ncycle
!
   REAL   , DIMENSION(1:MAXDSET,1:MAXDAT)     :: exp_data           ! (MAXDSET,MAXDAT)
   REAL   , DIMENSION(1:MAXDSET,1:MAXDAT)     :: exp_data_err       ! (MAXDSET,MAXDAT)
   REAL   , DIMENSION(1:MAXDSET,1:MAXDAT)     :: exp_wic            ! (MAXDSET,MAXDAT)
   REAL   , DIMENSION(1:MAXDSET)              :: exp_scal           ! (MAXDSET)
   REAL   , DIMENSION(1:MAXDSET)              :: exp_qmax           ! (MAXDSET)
   REAL   , DIMENSION(1:MAXDSET)              :: exp_rmin           ! (MAXDSET)
   REAL   , DIMENSION(1:MAXDSET)              :: exp_rmax           ! (MAXDSET)
   REAL   , DIMENSION(1:MAXDSET)              :: exp_deltar         ! (MAXDSET)
   REAL   , DIMENSION(1:MAXDSET)              :: exp_w              ! (MAXDSET)
!
   REAL   , DIMENSION(0:MAXELEM)              :: elem_conc          ! (0:MAXELEM)
   REAL   , DIMENSION(0:MAXELEM,1:MAXDSET)    :: elem_w             ! (0:MAXELEM,MAXDSET)
!
   REAL                                       :: pdf_xq
   REAL   , DIMENSION(9,0:MAXELEM,1:MAXDSET)  :: pdf_scat           ! (9,0:MAXELEM,MAXDSET)
   REAL   , DIMENSION(1:MAXDSET)              :: pdf_bave2          ! (MAXDSET)
   REAL   , DIMENSION(1:MAXELEM,0:MAXELEM,1:MAXDSET) :: pdf_weights ! (MAXELEM,0:MAXELEM,MAXDSET)
   REAL   , DIMENSION(1:MAXDAT )              :: pdf_diff           ! (MAXDAT)
   REAL   , DIMENSION(1:MAXDAT )              :: pdf_diff_err       ! (MAXDAT)
   REAL                                       :: pdf_diff_rmin
   REAL                                       :: pdf_diff_rmax
   REAL                                       :: pdf_diff_deltar
!
   REAL   , DIMENSION(1:MAXDAT )              :: fit_dat_x          ! (MAXDAT)
   REAL   , DIMENSION(1:MAXDAT )              :: fit_dat_y          ! (MAXDAT)
   REAL   , DIMENSION(1:MAXDAT )              :: fit_dat_w          ! (MAXDAT)
   REAL   , DIMENSION(1:MAXPARA)              :: p                  ! (MAXPARA)
   REAL   , DIMENSION(1:MAXPARA)              :: dp                 ! (MAXPARA)
   REAL   , DIMENSION(1:MAXPARA)              :: pinc               ! (MAXPARA)
   REAL   , DIMENSION(1:MAXPARA,1:MAXPARA)    :: cl                 ! (MAXPARA,MAXPARA)
   REAL                                       :: r4
   REAL                                       :: rexp
   REAL                                       :: urf
!
   LOGICAL, DIMENSION(          1:MAXDSET)    :: exp_lxray          ! (MAXDSET)
   LOGICAL, DIMENSION(0:MAXELEM,1:MAXDSET)    :: pdf_scat_int       ! (0:MAXELEM,MAXDSET)
!
!     COMMON /expt/ exp_name,exp_lxray,exp_nd,exp_np,exp_data,           &
!    &       exp_wic,exp_scal,exp_rmin,exp_rmax,exp_deltar,             &
!    &       exp_qmax,exp_w,exp_data_err
!     COMMON /elem/ elem_name,elem_diff,elem_conc,elem_n,elem_w
!     COMMON /pdf/  pdf_scat,pdf_xq,pdf_diff,pdf_diff_np,               &
!    &       pdf_diff_rmin,pdf_diff_rmax,pdf_diff_deltar,               &
!    &       pdf_bave2,pdf_weights,pdf_diff_err,                        &
!    &       pdf_scat_int
!     COMMON /fit/ fit_dat_x,fit_dat_y,fit_dat_w,                       &
!    &       fit_np,p,dp,pinc,cl,npara,r4,rexp,urf,ncycle
!
END MODULE mixscat_mod
