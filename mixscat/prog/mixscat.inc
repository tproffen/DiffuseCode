!*****7****************************************************************
      CHARACTER*80      exp_name(MAXDSET)
!
      CHARACTER*4      elem_name(0:MAXELEM)
!
      INTEGER elem_diff(2)
      INTEGER elem_n
!
      INTEGER exp_np(MAXDSET)
      INTEGER exp_nd
!
      INTEGER pdf_diff_np
!
      INTEGER fit_np
      INTEGER npara,ncycle
!
      REAL exp_data(MAXDSET,MAXDAT)
      REAL exp_data_err(MAXDSET,MAXDAT)
      REAL exp_wic(MAXDSET,MAXDAT)
      REAL exp_scal(MAXDSET)
      REAL exp_qmax(MAXDSET)
      REAL exp_rmin(MAXDSET)
      REAL exp_rmax(MAXDSET)
      REAL exp_deltar(MAXDSET)
      REAL exp_w(MAXDSET)
!
      REAL    elem_conc(0:MAXELEM)
      REAL    elem_w(0:MAXELEM,MAXDSET)
!
      REAL pdf_xq
      REAL pdf_scat(9,0:MAXELEM,MAXDSET)
      REAL pdf_bave2(MAXDSET)
      REAL pdf_weights(MAXELEM,0:MAXELEM,MAXDSET)
      REAL pdf_diff(MAXDAT)
      REAL pdf_diff_err(MAXDAT)
      REAL pdf_diff_rmin,pdf_diff_rmax,pdf_diff_deltar
!
      REAL fit_dat_x(MAXDAT)
      REAL fit_dat_y(MAXDAT)
      REAL fit_dat_w(MAXDAT)
      REAL p(MAXPARA),dp(MAXPARA),pinc(MAXPARA)
      REAL cl(MAXPARA,MAXPARA)
      REAL r4,rexp,urf
!
      LOGICAL exp_lxray(MAXDSET)
      LOGICAL pdf_scat_int(0:MAXELEM,MAXDSET)
!
      COMMON /expt/ exp_name,exp_lxray,exp_nd,exp_np,exp_data,           &
     &       exp_wic,exp_scal,exp_rmin,exp_rmax,exp_deltar,             &
     &       exp_qmax,exp_w,exp_data_err
      COMMON /elem/ elem_name,elem_diff,elem_conc,elem_n,elem_w
      COMMON /pdf/  pdf_scat,pdf_xq,pdf_diff,pdf_diff_np,               &
     &       pdf_diff_rmin,pdf_diff_rmax,pdf_diff_deltar,               &
     &       pdf_bave2,pdf_weights,pdf_diff_err,                        &
     &       pdf_scat_int
      COMMON /fit/ fit_dat_x,fit_dat_y,fit_dat_w,                       &
     &       fit_np,p,dp,pinc,cl,npara,r4,rexp,urf,ncycle
