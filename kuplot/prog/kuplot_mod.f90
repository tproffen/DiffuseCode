MODULE kuplot_mod
!
USE kuplot_config
USE precision_mod
!
IMPLICIT NONE
PUBLIC
SAVE
!
!****7*****************************************************************
!     This include file contains all COMMON variables for the
!     plotting routines of KUPLOT.
!*****7*****************************************************************
!
      INTEGER, PARAMETER :: ndev = 7
      INTEGER, PARAMETER :: x11  = 1
      INTEGER, PARAMETER :: pic  = 2
      INTEGER, PARAMETER :: ps   = 3
      INTEGER, PARAMETER :: vpic = 4
      INTEGER, PARAMETER :: vps  = 5
      INTEGER, PARAMETER :: png  = 6
      INTEGER, PARAMETER :: lat  = 7
!
!
      INTEGER, PARAMETER :: if_left   = 1
      INTEGER, PARAMETER :: if_centre = 2
      INTEGER, PARAMETER :: if_right  = 3
!
!
      INTEGER, PARAMETER :: nfon = 6
!
!
!------ Device handling
!
      CHARACTER(LEN=80), DIMENSION(NDEV) :: dev_name ! (ndev)
      CHARACTER(LEN=80), DIMENSION(NDEV) :: dev_prn  ! (ndev)
      INTEGER, DIMENSION(MAXWIN,NDEV)    :: dev_id     ! (maxwin,ndev)
      REAL,    DIMENSION(MAXWIN,NDEV)    :: dev_sf     ! (maxwin,ndev)
      REAL,    DIMENSION(MAXWIN,2   )    :: dev_draw   ! (maxwin,2)
      REAL,    DIMENSION(MAXWIN     )    :: dev_width  ! (maxwin)
      REAL,    DIMENSION(MAXWIN     )    :: dev_height ! (maxwin)
!
!     Color map 
      INTEGER, PARAMETER                 :: COL_MAP_NONE=0
      INTEGER, PARAMETER                 :: COL_MAP_GRAY=0
      INTEGER, PARAMETER                 :: COL_MAP_FIRE=1
      INTEGER, PARAMETER                 :: COL_MAP_ICE =2
      INTEGER, PARAMETER                 :: COL_MAP_THER=3
      INTEGER, PARAMETER                 :: COL_MAP_KUPL=4
      REAL,    DIMENSION(MAXWIN,MAXCOL,3):: col_map    ! (maxwin,maxcol,3)
      INTEGER                            :: col_map_type
!
!     COMMON /dev/  dev_name,dev_prn,dev_sf,dev_height,  &
!    &       dev_id,dev_width,dev_draw,col_map
!
!------ Menu button definition variables
!
      CHARACTER(LEN=80), DIMENSION(MAXB) :: btext ! (maxb)
      CHARACTER(LEN= 1)                  :: butt_l,butt_m,butt_r
      REAL, DIMENSION(MAXB)              :: bx    ! (maxb)
      REAL, DIMENSION(MAXB)              :: by    ! (maxb)
      REAL, DIMENSION(MAXB)              :: bw    ! (maxb)
      REAL, DIMENSION(MAXB)              :: bh    ! (maxb)
!
!     COMMON      /butt/      btext,bx,by,bw,bh,butt_l,butt_r,butt_m
!
!------ Data set information
!
      CHARACTER(LEN=200), DIMENSION(MAXKURVTOT) :: fname ! (maxkurvtot)
      CHARACTER(LEN=  2), DIMENSION(MAXKURVTOT) :: fform ! (maxkurvtot)
      REAL   , DIMENSION(MAXKURVTOT) :: xmax ! (maxkurvtot)
      REAL   , DIMENSION(MAXKURVTOT) :: xmin ! (maxkurvtot)
      REAL   , DIMENSION(MAXKURVTOT) :: ymax ! (maxkurvtot)
      REAL   , DIMENSION(MAXKURVTOT) :: ymin ! (maxkurvtot)
      INTEGER, DIMENSION(MAXKURVTOT) :: len  ! (maxkurvtot)
      INTEGER                        :: iz
!
!     COMMON      /info/  fname,fform,xmax,xmin,ymax,ymin,len,iz
!
!------ Plotting PARAMETERs
!
      CHARACTER(LEN=80)      titel(maxwin,maxframe,2)
      CHARACTER(LEN=60)      achse(maxwin,maxframe,3)
      CHARACTER(LEN=40)      ftext(maxwin,maxframe)
      CHARACTER(LEN=40)      clegend(maxwin,maxframe,maxkurvtot)
      CHARACTER(LEN=40)      antext(maxwin,maxframe,maxan)
      REAL colour(maxwin,0:20,3)
      REAL frame(maxwin,maxframe,4)
      REAL ex(maxwin,maxframe,2),ey(maxwin,maxframe,2)
      REAL pex(maxwin,maxframe,2),pey(maxwin,maxframe,2)
      REAL t (maxwin,maxframe,2)
      REAL pt(maxwin,maxframe,2)
      REAL yskal_u(maxwin,maxframe)
      REAL yskal(maxwin,maxframe)
      REAL anx(maxwin,maxframe,maxan)
      REAL any(maxwin,maxframe,maxan)
      REAL antx(maxwin,maxframe,maxan)
      REAL anty(maxwin,maxframe,maxan)
      REAL ibuf(maxwin,maxframe,4)
      REAL shear(maxwin,maxframe)
      REAL info_orig(maxwin,maxframe,2)
      REAL anangle(maxwin,maxframe,maxan)
      REAL lab_angle(maxwin,maxframe,2)
      REAL lab_d(maxwin,maxframe,2)
      REAL ax_d(maxwin,maxframe,2)
      REAL tick_ma_h(maxwin,maxframe,2)
      REAL tick_mi_h(maxwin,maxframe,2)
      INTEGER tick_nsub(maxwin,maxframe,2)
      INTEGER anjust(maxwin,maxframe,maxan)
      INTEGER ilegend(maxwin,maxframe,maxkurvtot)
      INTEGER infra(maxwin,maxframe,maxkurvtot)
      INTEGER ibox(maxwin,maxframe)
      INTEGER icol_old(maxwin)
      INTEGER iaf(maxwin)
      INTEGER iframe,iwin
      LOGICAL lachse(maxwin,maxframe,3)
      LOGICAL lyskal(maxwin,maxframe)
      LOGICAL igrid(maxwin,maxframe)
      LOGICAL orient(maxwin)
      LOGICAL lgui(maxwin)
      LOGICAL sfl(maxwin,maxframe)
      LOGICAL iden(maxwin)
!
!     COMMON      /plpar/      ibox,igrid,ilegend,infra,ex,ey,t,        &
!    &       lgui,lachse,                                               &
!    &       yskal,ibuf,orient,anx,any,anjust,achse,titel,              &
!    &       clegend,ftext,antext,icol_old,yskal_u,lyskal,              &
!    &       shear,sfl,frame,iframe,iaf,info_orig,anangle,              &
!    &       antx,anty,iwin,lab_angle,tick_ma_h,tick_mi_h,              &
!    &       tick_nsub,lab_d,ax_d,iden,pex,pey,pt
!
!------ Plotting style PARAMETERs
!
      REAL sizemark(maxwin,maxframe,maxkurvtot)
      INTEGER rel_mark(maxwin,maxframe,maxkurvtot)
      REAL linewid(maxwin,maxframe,0:maxkurvtot)
      REAL fonsize (maxwin,maxframe,nfon)
      REAL fonscal (maxwin,maxframe)
      REAL frback(maxwin,maxframe,3)
      REAL fillrange(maxwin,maxframe,maxkurvtot,4)
       INTEGER ilinecol(maxwin,maxframe,0:maxkurvtot)
      INTEGER imarkcol(maxwin,maxframe,maxkurvtot)
      INTEGER ierrcol (maxwin,maxframe,maxkurvtot)
       INTEGER ilinetyp(maxwin,maxframe,maxkurvtot)
      INTEGER imarktyp(maxwin,maxframe,maxkurvtot)
      INTEGER ilineart(maxwin,maxframe,maxkurvtot)
      INTEGER imarkmax(maxwin,maxframe,maxkurvtot)
      INTEGER ifillcol(maxwin,maxframe,maxkurvtot)
      INTEGER ifilltyp(maxwin,maxframe,maxkurvtot)
      INTEGER ierr(maxwin,maxframe,maxkurvtot)
      INTEGER frjust(maxwin,maxframe)
      INTEGER fon_id (maxwin,maxframe,nfon)
      INTEGER foncol (maxwin,maxframe,nfon)
      LOGICAL ifname(maxwin,maxframe) 
      LOGICAL tot_frame(maxwin)
!
      INTEGER ifen
!
!     COMMON  /plwer/ ilinecol,imarkcol,ierr,ilinetyp,imarktyp,         &
!    &        ifen,ilineart,imarkmax,ifname,sizemark,                   &
!    &        linewid,fonsize,fon_id,tot_frame,foncol,                  &
!    &        frjust,frback,fonscal,colour,ierrcol,                     &
!    &       ifillcol,ifilltyp,fillrange
!
!------ 3D data set information
!
      REAL zmax(maxkurvtot),zmin(maxkurvtot) 
      REAL z_min(maxwin,maxframe,maxhl)
      REAL z_inc(maxwin,maxframe,maxhl)
      REAL pgmlow,pgmhigh
      INTEGER nx(maxkurvtot),ny(maxkurvtot)
      INTEGER nz(maxwin,maxframe,maxhl)
      INTEGER iho(maxwin,maxframe)
      INTEGER hlinecol(maxwin,maxframe,maxkurvtot,maxhl)
      INTEGER hlinetyp(maxwin,maxframe,maxkurvtot,maxhl)
      INTEGER hlineart(maxwin,maxframe,maxkurvtot)
      INTEGER hlabel(maxwin,maxframe,maxkurvtot)
      INTEGER         exclude9999
      LOGICAL lni(maxkurvtot)
!
!     COMMON /nipl/ lni,iho,nz,z_min,z_inc,nx,ny,zmin,zmax,             &
!    &              hlinecol,hlinetyp,hlineart,pgmlow,pgmhigh,          &
!    &       hlabel,exclude9999
!
!------ Parameter frame 
!
      CHARACTER(LEN=40)      fpara(maxwin)
      INTEGER ipara(maxwin)
      INTEGER ilength(maxwin)
      INTEGER iwidth(maxwin)
!
!     COMMON /para/ ipara,fpara,ilength,iwidth
!
!------ Fit PARAMETERs
!
      CHARACTER(LEN=PREC_STRING) ::     fit_func
      CHARACTER(LEN=4)      ftyp
      CHARACTER(LEN=3)      wtyp
      CHARACTER(LEN=1)      dummy_spacer
      REAL w(maxarray)
      REAL p(maxpara),pinc(maxpara),dp(maxpara)
      REAL p_merk(maxpara),pinc_merk(maxpara)
      REAL cl(maxpara,maxpara)
      REAL urf
      REAL zalt,zwert,zdif
      REAL r4,re,fend
      REAL            wval
      REAL tmp(maxpara)
      REAL p_origin
      INTEGER ncycle,npara,ikfit,ikfit2,ikcal,ikdif
      INTEGER fit_ikcal(maxkurvtot)
      INTEGER fit_ikdif(maxkurvtot)
      INTEGER fit_lfunc,fit_ifen
      INTEGER np1,np2,np3
      INTEGER n_backgrd
      LOGICAL fstart,frall,ikfirst(maxkurvtot)
!
!     COMMON  /fit/ w,p,dp,pinc,cl,urf,ncycle,npara,np1,np2,np3,        &
!    &              zalt,zwert,zdif,fend,r4,re,fstart,p_merk,           &
!    &              pinc_merk,ikfit,ikfit2,ikcal,ikdif,ikfirst,         &
!    &       fit_ifen,fit_ikcal,fit_ikdif,fit_lfunc,                    &
!    &       ftyp,fit_func,tmp,wval,wtyp,dummy_spacer,frall,n_backgrd,  &
!    &       p_origin
!
!------ Bonds
!
      REAL bond_rad(maxwin,maxframe,maxbond)
      REAL bond_sig(maxwin,maxframe,maxbond)
      REAL bond_lwid(maxwin,maxframe,maxbond)
      INTEGER bond_lcol(maxwin,maxframe,maxbond)
      INTEGER bond_ltyp(maxwin,maxframe,maxbond)
!
!     COMMON /bond/ bond_rad,bond_sig,bond_lwid,bond_lcol,              &
!    &              bond_ltyp
!
!------ Data
!
      REAL x(maxarray),dx(maxarray)
      REAL y(maxarray),dy(maxarray)
      REAL z(maxarray)
      INTEGER offxy(0:maxkurvtot)
      INTEGER offz (0:maxkurvtot)
!
!     COMMON  /data/ x,y,z,dx,dy,offxy,offz
!
!------ GSAS related stuff to remember
!
      REAL tof_offset
!
!     COMMON  /gsas/ tof_offset
!
      LOGICAL :: l_two_col = .false.  ! save two column xy-files
END MODULE kuplot_mod
