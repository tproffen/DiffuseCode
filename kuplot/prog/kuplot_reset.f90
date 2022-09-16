module kuplot_reset_mod
!
! Reset KUPLOT to start conditions
!+
!
contains
!
!****7******************************************************************
SUBROUTINE kuplot_do_reset(zeile, lp) 
!+                                                                      
!     Reset KUPLOT                                                      
!-                                                                      
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
      USE param_mod
USE str_comp_mod
      USE take_param_mod
use kuplot_blk_mod
use lib_data_struc_h5
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz 
      INTEGER iw, i, j 
!                                                                       
!------ get parameters                                                  
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ No parameters                                                   
!                                                                       
      IF (ianz.eq.0) then 
         WRITE (output_io, 1000) 
         iz = 1 
         iframe = 1 
!                                                                       
         DO i = 1, maxkurvtot 
         lni (i) = .false. 
         ikfirst (i) = .true. 
         ENDDO 
!                                                                       
         DO iw = 1, maxwin 
         iaf (iw) = 1 
         frame (iw, 1, 1) = 0.0 
         frame (iw, 1, 2) = 0.0 
         frame (iw, 1, 3) = 1.0 
         frame (iw, 1, 4) = 1.0 
!                                                                       
         DO i = 1, maxframe 
         ex (iw, i, 1) = - 9999. 
         ey (iw, i, 1) = - 9999. 
         lyskal (iw, i) = .false. 
         t (iw, i, 1) = - 9999.0 
         t (iw, i, 2) = - 9999.0 
         frback (iw, i, 1) = 1.0 
         frback (iw, i, 2) = 1.0 
         frback (iw, i, 3) = 1.0 
         fonscal (iw, i) = 1.0 
         shear (iw, i) = 90.0 
!                                                                       
         DO j = 1, maxkurvtot 
         infra (iw, i, j) = j 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Subcommand                                                      
!                                                                       
      ELSE 
         IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) ) then 
            WRITE (output_io, 1100) 
            CALL kuplot_initarrays 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
!
dx    = 0.0
dy    = 0.0
x     = 0.0
y     = 0.0
xmin  = 0.0
xmax  = 0.0
ymin  = 0.0
ymax  = 0.0
z     = 0.0
lenc  = 0
offxy = 0
offz  = 0
fname = ' '
fform = ' '
lh5   = .false.
lni   = .false.
ku_ndims = 1
!
CALL dgl5_reset
!                                                                       
 1000 FORMAT     (1x,'Resetting ..') 
 1100 FORMAT     (1x,'Resetting - all parameters ..') 
      END SUBROUTINE kuplot_do_reset                       
!
!****7******************************************************************
!
end module kuplot_reset_mod
