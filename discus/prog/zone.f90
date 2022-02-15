MODULE zone
!
use precision_mod
IMPLICIT NONE
!
PRIVATE
PUBLIC zone_setup, zone_project
!
   REAL(kind=PREC_DP), DIMENSION(3)  :: absc
   REAL(kind=PREC_DP), DIMENSION(3)  :: ordi
   REAL(kind=PREC_DP), DIMENSION(3)  :: top
!
CONTAINS
   SUBROUTINE zone_setup
!
   USE crystal_mod
   USE diffuse_mod 
   USE metric_mod
   USE output_mod
!
   USE param_mod
   USE precision_mod
!
   IMPLICIT NONE     
!
   REAL(kind=PREC_DP), PARAMETER :: EPS = 0.0001
!
   CHARACTER(LEN=PREC_STRING) :: line
   INTEGER             :: length
   REAL(kind=PREC_DP)  :: dstar
   REAL(kind=PREC_DP)  :: tthmax
   REAL(kind=PREC_DP)  :: lmin, lmax
!
!     Transform zone axis into reciprocal space
!
!  CALL trans(zone_uvw, cr_gten, top, 3)
   top = matmul(cr_gten, zone_uvw)
!
!     Define directions in reciprocal space
!
   IF ( zone_uvw(2) < EPS .and. zone_uvw(3) < EPS ) THEN   ! [U 0 0]
!
!        abs = [0 1 0]; ord = zone x abs ~ c*; top = zone
      absc(1)  = 0.0
      absc(2)  = 1.0
      absc(3)  = 0.0
      extr_abs = 2 
      extr_ord = 3 
      extr_top = 1 
   ELSEIF ( zone_uvw(1) < EPS .and. zone_uvw(3) < EPS ) THEN   ! [0 V 0]
!
!        abs = [0 0 1]; ord = zone x abs ~ a*; top = zone
      absc(1)  = 0.0
      absc(2)  = 0.0
      absc(3)  = 1.0
      extr_abs = 3 
      extr_ord = 1 
      extr_top = 2 
   ELSEIF ( zone_uvw(1) < EPS .and. zone_uvw(2) < EPS ) THEN   ! [0 0 W]
!
!        abs = [1 0 0]; ord = zone x abs ~ b*; top = zone
      absc(1)  = 1.0
      absc(2)  = 0.0
      absc(3)  = 0.0
      extr_abs = 1 
      extr_ord = 2 
      extr_top = 3 
   ELSEIF (                         zone_uvw(3) < EPS ) THEN   ! [U V 0]
!
!        abs = [0 0 1]; ord = zone x abs ~  a* - b*; top = zone
      absc(1)  = 0.0
      absc(2)  = 0.0
      absc(3)  = 1.0
      extr_abs = 3 
      extr_ord = 1 
      extr_top = 2 
   ELSEIF (                         zone_uvw(2) < EPS ) THEN   ! [U 0 W]
!
!        abs = [0 1 0]; ord = zone x abs ~ -a* + c*; top = zone
      absc(1)  = 0.0
      absc(2)  = 1.0
      absc(3)  = 0.0
      extr_abs = 2 
      extr_ord = 3 
      extr_top = 1 
   ELSEIF (                         zone_uvw(3) < EPS ) THEN   ! [0 V W]
!
!        abs = [1 0 0]; ord = zone x abs ~ +b* - c*; top = zone
      absc(1)  = 1.0
      absc(2)  = 0.0
      absc(3)  = 0.0
      extr_abs = 1 
      extr_ord = 2 
      extr_top = 3 
   ELSEIF(ABS(ABS(zone_uvw(1))-ABS(zone_uvw(2))) < EPS  .and. &
          ABS(ABS(zone_uvw(1))-ABS(zone_uvw(3))) < EPS       ) THEN ! [UUU]
      IF( ABS(zone_uvw(1)-zone_uvw(2)) < EPS ) THEN  ! [+U+U +-U] or [-U-U +-U]
         absc(1)  = 1.0
         absc(2)  =-1.0
         absc(3)  = 0.0
         extr_abs = 1 
         extr_ord = 3 
         extr_top = 2 
      ELSE                                           ! [+U-U +-U] or [+U-U +-U]
         absc(1)  = 1.0
         absc(2)  = 1.0
         absc(3)  = 0.0
         extr_abs = 1 
         extr_ord = 3 
         extr_top = 2 
      ENDIF
   ELSE
      absc(1)  =  1.0
      absc(2)  = -1.0*zone_uvw(1)/zone_uvw(2)
      absc(3)  =  0.0
      extr_abs = 1 
      extr_ord = 2 
      extr_top = 3 
   ENDIF
   line = ' '
   WRITE(line,2000) zone_uvw(:), absc, 'drr'
   length = 81
   CALL vprod(line,length)
   ordi(1) = res_para(1)
   ordi(2) = res_para(2)
   ordi(3) = res_para(3)
!
   dstar   = sqrt(skalpro(top,top,cr_rten))     ! Normalize top to one A^-1
   top(1)  = top(1)/dstar
   top(2)  = top(2)/dstar
   top(3)  = top(3)/dstar
   dstar   = sqrt(skalpro(absc,absc,cr_rten))   ! normalize absc to resolution
   absc(1) = absc(1)*zone_res/dstar
   absc(2) = absc(2)*zone_res/dstar
   absc(3) = absc(3)*zone_res/dstar
   dstar   = sqrt(skalpro(ordi,ordi,cr_rten))   ! normalize ordi to resolution
   ordi(1) = ordi(1)*zone_res/dstar
   ordi(2) = ordi(2)*zone_res/dstar
   ordi(3) = ordi(3)*zone_res/dstar
   zone_ewald(1) = top(1)/rlambda
   zone_ewald(2) = top(2)/rlambda
   zone_ewald(3) = top(3)/rlambda
!
   tthmax  =             asin(zone_res  /(1./rlambda - zone_delta_d))
   lmax    = 1./rlambda - cos(tthmax)*   (1./rlambda - zone_delta_d)
   lmin    = 1.*zone_delta_d
!
   eck(1,1) = - absc(1) - ordi(1) - top(1)*lmin
   eck(2,1) = - absc(2) - ordi(2) - top(2)*lmin
   eck(3,1) = - absc(3) - ordi(3) - top(3)*lmin
!
   eck(1,2) = + absc(1) - ordi(1) - top(1)*lmin
   eck(2,2) = + absc(2) - ordi(2) - top(2)*lmin
   eck(3,2) = + absc(3) - ordi(3) - top(3)*lmin
!
   eck(1,3) = - absc(1) + ordi(1) - top(1)*lmin
   eck(2,3) = - absc(2) + ordi(2) - top(2)*lmin
   eck(3,3) = - absc(3) + ordi(3) - top(3)*lmin
!
   eck(1,4) = - absc(1) - ordi(1) + top(1)*lmax
   eck(2,4) = - absc(2) - ordi(2) + top(2)*lmax
   eck(3,4) = - absc(3) - ordi(3) + top(3)*lmax
!
   inc(3)   = NINT(inc(1) * (lmax+lmin)/zone_res/2.)
!
   vi(1,1)  = (eck(1,2)-eck(1,1))/(inc(1)-1)
   vi(2,1)  = (eck(2,2)-eck(2,1))/(inc(1)-1)
   vi(3,1)  = (eck(3,2)-eck(3,1))/(inc(1)-1)
!
   vi(1,2)  = (eck(1,3)-eck(1,1))/(inc(2)-1)
   vi(2,2)  = (eck(2,3)-eck(2,1))/(inc(2)-1)
   vi(3,2)  = (eck(3,3)-eck(3,1))/(inc(2)-1)
!
   vi(1,3)  = (eck(1,4)-eck(1,1))/(inc(3)-1)
   vi(2,3)  = (eck(2,4)-eck(2,1))/(inc(3)-1)
   vi(3,3)  = (eck(3,4)-eck(3,1))/(inc(3)-1)
!
   2000  FORMAT(6(F12.6,','),a3)
   END SUBROUTINE zone_setup
!
   SUBROUTINE zone_project
!
   USE crystal_mod
   USE diffuse_mod
   USE discus_allocate_appl_mod
   USE fourier_sup
   USE metric_mod
   USE output_mod
!
   USE param_mod
use matrix_mod
   USE precision_mod
!
   IMPLICIT NONE
!
   CHARACTER(LEN=PREC_STRING)   :: line
   INTEGER               :: i,j,l, ii, k,m
   INTEGER               :: length
   INTEGER               :: hmax, kmax, lmax
   INTEGER              :: n_qxy    ! required size in reciprocal space this run
   INTEGER              :: n_nscat  ! required no of atom types right now
   INTEGER              :: n_natoms ! required no of atoms
!
   REAL(kind=PREC_DP), DIMENSION(1:3)  :: rvec, rproj, rres
real(kind=PREC_DP), dimension(1:3) :: rdif
   REAL(kind=PREC_DP), DIMENSION(3,3)  :: matrix, inverse
   REAL(kind=PREC_DP)    :: dstar
   REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: layer
!
   INTEGER, DIMENSION(1:3)     :: z_inc
   REAL(kind=PREC_DP), dimension(1:3,1:4) :: z_eck
   REAL(KIND=PREC_DP), DIMENSION(1:3,1:3) :: z_vi 
!
   ALLOCATE(layer(1:inc(1),1:inc(2)))
   layer(:,:) = 0.0D0
   ii = 0
   DO i = 1, inc(1)
      DO j = 1, inc(2)
         DO l = 1, inc(3)
            ii = ii+1
            rvec(1) = eck(1,1)+(i-1)*vi(1,1) + (j-1)*vi(1,2) + (l-1)*vi(1,3)
            rvec(2) = eck(2,1)+(i-1)*vi(2,1) + (j-1)*vi(2,2) + (l-1)*vi(2,3)
            rvec(3) = eck(3,1)+(i-1)*vi(3,1) + (j-1)*vi(3,2) + (l-1)*vi(3,3)
            rdif(1) = rvec(1) - zone_ewald(1)
            rdif(2) = rvec(2) - zone_ewald(2)
            rdif(3) = rvec(3) - zone_ewald(3)
            dstar   = sqrt(skalpro(rdif,rdif,cr_rten))   ! 
            IF(abs(dstar-1./rlambda) < zone_delta_d) THEN
               layer(i,j) = layer(i,j) + dsi(ii)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
!
!  Add Bragg reflections in zone
!
   z_inc(:)   = inc(:)
   z_eck(:,:) = eck(:,:)
   z_vi (:,:) = vi (:,:)
   matrix(:,1) = vi(:,1)
   matrix(:,2) = vi(:,2)
   matrix(:,3) = vi(:,3) 
!  CALL invmat(inverse,matrix)
   CALL matinv(matrix, inverse)
!
   eck(:,:) = 0.0
   hmax     = 1*INT(zone_res*cr_a0(1)+0.)
   kmax     = 1*INT(zone_res*cr_a0(2)+0.)
   lmax     = 1*INT(zone_res*cr_a0(3)+0.)
!
   eck(1,1) = REAL(-hmax)
   eck(2,1) = REAL(-kmax)
   eck(3,1) = REAL(-lmax)
!
   eck(1,2) = REAL( hmax)
   eck(2,2) = REAL(-kmax)
   eck(3,2) = REAL(-lmax)
!
   eck(1,3) = REAL(-hmax)
   eck(2,3) = REAL( kmax)
   eck(3,3) = REAL(-lmax)
!
   eck(1,4) = REAL(-hmax)
   eck(2,4) = REAL(-kmax)
   eck(3,4) = REAL( lmax)
!
   inc(1)   = 2*hmax + 1
   inc(2)   = 2*kmax + 1
   inc(3)   = 2*lmax + 1
!
   vi (:,:) = 0.0
   vi(1,1)  = 1.0
   vi(2,2)  = 1.0
   vi(3,3)  = 1.0
!
   fave     = 0.0
   ilots    = LOT_OFF
   nlots    = 1
!
   IF (inc(1) * inc(2) *inc(3) .gt. MAXQXY  .OR.    &
       cr_natoms > DIF_MAXAT                .OR.    &
       cr_nscat>DIF_MAXSCAT              ) THEN
      n_qxy    = MAX(n_qxy,inc(1)*inc(2)*inc(3),MAXQXY)
      n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
      n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
      call alloc_diffuse (n_qxy, n_nscat, n_natoms)
      IF (ier_num.ne.0) THEN
         RETURN
      ENDIF
   ENDIF

   CALL four_run
!
   length = 95
   ii = 0
   DO l = 1, inc(3)
      DO j = 1, inc(2)
         DO i = 1, inc(1)
            ii = ii+1
            rvec(1) = eck(1,1)+(i-1)*vi(1,1) + (j-1)*vi(1,2) + (l-1)*vi(1,3)
            rvec(2) = eck(2,1)+(i-1)*vi(2,1) + (j-1)*vi(2,2) + (l-1)*vi(2,3)
            rvec(3) = eck(3,1)+(i-1)*vi(3,1) + (j-1)*vi(3,2) + (l-1)*vi(3,3)
            rdif(1) = rvec(1) - zone_ewald(1)
            rdif(2) = rvec(2) - zone_ewald(2)
            rdif(3) = rvec(3) - zone_ewald(3)
            dstar   = sqrt(skalpro(rdif,rdif,cr_rten))   ! 
            IF(abs(dstar-1./rlambda) < zone_delta_d) THEN
               ii = (i-1)*inc(3)*inc(2) + (j-1)*inc(3) + l
               WRITE(line,3000) rvec(:), zone_ewald(:), 'rrrr'
               CALL do_proj(line, length)
               rproj(:) = res_para(4:6)
               rres(:)  = 0.0
               DO k=1,3
                  DO m=1,3
                     rres(k) = rres(k) + inverse(k,m)*rproj(m)
                  ENDDO
               ENDDO
               k = NINT(rres(1)+z_inc(1)/2.)
               m = NINT(rres(2)+z_inc(2)/2.)
               IF(k>0 .and. k<z_inc(1) .and. &
                  m>0 .and. m<z_inc(2)       ) THEN
                  layer(k,m) = layer(k,m) + dsi(ii)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   3000 FORMAT(6(F12.6,','),a4)
!
   ii = 0
   l = 1
   DO i = 1, z_inc(1)
      DO j = 1, z_inc(2)
         ii = ii+1
         dsi(ii) = layer(i,j)
      ENDDO
   ENDDO
   DEALLOCATE(layer)
!
!
   eck(1,1) = - zone_res
   eck(2,1) = - zone_res
   eck(3,1) =   0.0
!
   eck(1,2) = + zone_res
   eck(2,2) = - zone_res
   eck(3,2) =   0.0
!
   eck(1,3) = - zone_res
   eck(2,3) = + zone_res
   eck(3,3) =   0.0
!
   eck(1,4) = - zone_res
   eck(2,4) = - zone_res
   eck(3,4) =   0.0
!
   inc(1)   = z_inc(1)
   inc(2)   = z_inc(2)
   inc(3)   = 1
!
   vi(1,1)  = (eck(1,2)-eck(1,1))/(inc(1)-1)
   vi(2,1)  = (eck(2,2)-eck(2,1))/(inc(1)-1)
   vi(3,1)  = (eck(3,2)-eck(3,1))/(inc(1)-1)
!
   vi(1,2)  = (eck(1,3)-eck(1,1))/(inc(2)-1)
   vi(2,2)  = (eck(2,3)-eck(2,1))/(inc(2)-1)
   vi(3,2)  = (eck(3,3)-eck(3,1))/(inc(2)-1)
!
   vi(1,3)  = 0.0
   vi(2,3)  = 0.0
   vi(3,3)  = 0.0
!
   extr_abs = 1
   extr_ord = 2
   extr_top = 3
!
   l_zone = .false.
!
   END SUBROUTINE zone_project
!
END MODULE zone
